suppressPackageStartupMessages({
  library(devtools)
  library(rio)
  library(tidyverse)
  library(reshape2)
  library(exact2x2)
})

# SETTINGS ----
project.root <- file.path("~", "projects", "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

suppressPackageStartupMessages(source(file.path(project.root, "00_utils", "cleanUtils.R")))
dir.create(file.path(project.root, "01_tables", "data"), showWarnings = FALSE)

# FUNCTIONS ----

#' Extract counts (ya, yb, na, nb) for a specific study
#'
#' @param analysis.unique.id Unique ID from KeyTable
#' @param varfun_name Name of the varfun function to use
#' @param Nmin.raw Min participants in a study
#' @param Nmin.cond Min participants in each condition
#' @param subset.type "all" # TODO: now its all using 'all'
#' @param aggregate If TRUE, returns a 1-row summary of totals. If FALSE, returns site-level counts.
extract_counts <- function(analysis.unique.id, varfun_name, Nmin.raw = 0, Nmin.cond = 0, subset.type = "all", aggregate = FALSE) {

  # GET LOOKUP TABLES
  ML2.key.all <- rio::import(file.path(project.root, "00_data", "ML2_KeyTable.csv"))
  ML2.key <- ML2.key.all[!is.na(ML2.key.all$unique.id) & ML2.key.all$unique.id == analysis.unique.id, ]

  if (nrow(ML2.key) == 0) return(NULL)

  # Get data frame
  if (ML2.key$study.slate == 1) {
    df <- rio::import(file.path(OSFdata.root, "!!RawData", "ML2_S1.csv"))
  } else {
    df <- rio::import(file.path(OSFdata.root, "!!RawData", "ML2_S2.csv"))
  }

  df$uID <- seq(1, nrow(df))
  ML2.in <- get.info(ML2.key, colnames(df), subset.type)
  ML2.id <- get.chain(ML2.in)
  df_filtered <- eval(parse(text = paste("df", ML2.id$df)))

  toRun <- decide.analysis(ML2.key, analysis.unique.id, 3, doAll = TRUE)
  if (nrow(df_filtered) <= 0) return(NULL)

  runGroups <- sort(na.exclude(unique(df_filtered[[toRun$ugroup]])))

  results_list <- list()
  for (g in seq_along(runGroups)) {
    ML2.sr_g <- check.nMin(df_filtered, runGroups[g], Nmin.raw, Nmin.cond, ML2.id, ML2.in)
    if (is.null(ML2.sr_g)) next

    varfun <- match.fun(varfun_name)
    ML2.var_g <- varfun(ML2.sr_g)
    ct <- table(Condition = ML2.var_g$Condition, Response = ML2.var_g$Response)

    results_list[[g]] <- data.frame(
      site = runGroups[g],
      n1_yes = as.numeric(ct[1, 1]),
      n2_yes = as.numeric(ct[2, 1]),
      n1_total = sum(ct[1, ]),
      n2_total = sum(ct[2, ]),
      stringsAsFactors = FALSE
    )
  }

  results_df <- do.call(rbind, results_list)
  if (is.null(results_df)) return(NULL)

  if (aggregate) {
    return(data.frame(
      site = "AGGREGATED",
      n1_yes = sum(results_df$n1_yes),
      n2_yes = sum(results_df$n2_yes),
      n1_total = sum(results_df$n1_total),
      n2_total = sum(results_df$n2_total),
      stringsAsFactors = FALSE
    ))
  }

  return(results_df)
}

#' Helper to format results for CSV export
format_results <- function(df) {
  num_cols <- intersect(colnames(df), c("p.value", "OR", "delta", "d"))
  for (col in num_cols) {
    df[[col]] <- format(df[[col]], scientific = TRUE, digits = 10)
  }
  return(df)
}

# MAIN EXECUTION ----

analysis_mapping <- data.frame(
  name = c("Hauser.1", "Hauser.2", "Hauser.3", "Hauser.4", "Hauser.5", "Hauser.6", "Tversky.1", "Rottenstreich.1"),
  id = c(39, 40, 41, 47, 48, 49, 46, 16),
  varfun = c("varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Tversky.1", "varfun.Rottenstreich.1"),
  stringsAsFactors = FALSE
)

all_global_results <- data.frame()
all_site_results <- data.frame()

for (i in 1:nrow(analysis_mapping)) {
  a_name <- analysis_mapping$name[i]
  a_id <- analysis_mapping$id[i]
  a_varfun <- analysis_mapping$varfun[i]

  cat("--------------------------------------------------\n")
  cat(paste0("\tProcessing: ", a_name, " (ID: ", a_id, ")\n"))

  # 1. Extract counts
  site_counts <- extract_counts(a_id, a_varfun, subset.type = "all")
  if (is.null(site_counts)) {
    cat(paste("Skipping", a_name, "- no data found.\n"))
    next
  }

  # Export individual count file (backward compatibility)
  counts_export <- site_counts %>% rename(sites = site, ya = n1_yes, yb = n2_yes, na = n1_total, nb = n2_total)
  export_path <- file.path(project.root, "01_tables", "data", paste0(gsub("\\.", "_", a_name), "_clean_tables.csv"))
  rio::export(counts_export, export_path)
  cat(paste0("Counts saved to ", export_path, "\n"))

  # 2. Run Site Fisher Tests
  # We need the alternative from the keytable
  ML2.key.all <- rio::import(file.path(project.root, "00_data", "ML2_KeyTable.csv"))
  ML2.key <- ML2.key.all[!is.na(ML2.key.all$unique.id) & ML2.key.all$unique.id == a_id, ]
  stat_params <- eval(parse(text = ML2.key$stat.params))
  base_alternative <- if (!is.null(stat_params$alternative)) stat_params$alternative else "two.sided"

  for (j in 1:nrow(site_counts)) {
    row <- site_counts[j, ]
    tab <- matrix(c(row$n1_yes, row$n1_total - row$n1_yes, row$n2_yes, row$n2_total - row$n2_yes), nrow = 2)

    alt_s <- base_alternative
    if (a_name == "Rottenstreich.1" && alt_s == "greater") alt_s <- "less"

    test_s <- tryCatch({ exact2x2::fisher.exact(tab, alternative = alt_s) }, error = function(e) NULL)

    if (!is.null(test_s)) {
      delta_s <- log(as.numeric(test_s$estimate))
      all_site_results <- rbind(all_site_results, data.frame(
        study.name = a_name, unique.id = a_id, site = row$site, alternative = alt_s,
        N = row$n1_total + row$n2_total, p.value = test_s$p.value, OR = as.numeric(test_s$estimate),
        delta = delta_s, d = delta_s / (pi / sqrt(3)), stringsAsFactors = FALSE
      ))
    }
  }

  # 3. Global Aggregation & Test
  global_counts <- extract_counts(a_id, a_varfun, subset.type = "all", aggregate = TRUE)
  tab_g <- matrix(c(global_counts$n1_yes, global_counts$n1_total - global_counts$n1_yes,
                    global_counts$n2_yes, global_counts$n2_total - global_counts$n2_yes), nrow = 2)

  alt_g <- base_alternative
  if (a_name == "Rottenstreich.1") {
    cat("\n[WARNING] Rottenstreich.1: Effect direction is REVERSED to 'less' \n")
    alt_g <- "less"
  }

  test_g <- exact2x2::fisher.exact(tab_g, alternative = alt_g)
  delta_g <- log(as.numeric(test_g$estimate))
  cohen_d_g <- delta_g / (pi / sqrt(3))

  cat(paste0("\nGlobal Results for ", a_name, "\n"))
  cat(paste0("Sample Size (N): ", global_counts$n1_total + global_counts$n2_total, "\n"))
  cat(paste0("Alternative : ", alt_g, "\n"))
  cat(paste0("p-value: ", test_g$p.value, "\n"))
  cat(paste0("Log-odds ratio (delta): ", round(delta_g, 4), "\n"))
  cat(paste0("Cohen's d: ", round(cohen_d_g, 4), "\n"))

  all_global_results <- rbind(all_global_results, data.frame(
    study.name = a_name, unique.id = a_id, alternative = alt_g,
    N = global_counts$n1_total + global_counts$n2_total, p.value = test_g$p.value,
    OR = as.numeric(test_g$estimate), delta = delta_g, d = cohen_d_g, stringsAsFactors = FALSE
  ))
}

# EXPORT MASTER CSVs ----
all_global_fmt <- format_results(all_global_results)
all_site_fmt <- format_results(all_site_results)

rio::export(all_global_fmt, file.path(project.root, "01_tables", "data", "global_fisher.csv"))
rio::export(all_site_fmt, file.path(project.root, "01_tables", "data", "site_fisher.csv"))

cat("\nMaster Fisher results exported to global_fisher.csv and site_fisher.csv\n")
cat("--------------------------------------------------\n")
