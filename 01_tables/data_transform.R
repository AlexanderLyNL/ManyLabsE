library(devtools)
library(rio)
library(tidyverse)
library(reshape2)
library(exact2x2)

project.root <- file.path("~", "projects", "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "cleanUtils.R"))
# source manylabRs_SOURCE to get all varfun functions
# Create output directory if it doesn't exist
dir.create(file.path(project.root, "01_tables", "data"), showWarnings = FALSE)

analysis_mapping <- data.frame(
  analysis.name = c("Hauser.1", "Hauser.2", "Hauser.3", "Hauser.4", "Hauser.5", "Hauser.6", "Tversky.1", "Rottenstreich.1"),
  analysis.unique.id = c(39, 40, 41, 47, 48, 49, 46, 16),
  varfun = c("varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Hauser.1", "varfun.Tversky.1", "varfun.Rottenstreich.1"),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(analysis_mapping)) {
  analysis.name <- analysis_mapping$analysis.name[i]
  analysis.unique.id <- analysis_mapping$analysis.unique.id[i]
  varfun_name <- analysis_mapping$varfun[i]

  cat(paste0("\n\n\tProcessing: ", analysis.name, "\n"))

  # ANALYSIS INFO ----
  analysis.type <- 3 # sites
  Nmin.raw <- 30 # minimal # of participants in a study
  Nmin.cond <- 15 # minimal # of participants in each condition
  subset.type <- "sites"

  # GET LOOKUP TABLES ----
  ML2.key.all <- rio::import(file.path(project.root, "00_data", "ML2_KeyTable.csv"))
  ML2.key <- ML2.key.all[!is.na(ML2.key.all$unique.id) & ML2.key.all$unique.id == analysis.unique.id, ]

  if (nrow(ML2.key) == 0) {
    print(paste("Skipping", analysis.name, "- no key found for unique.id", analysis.unique.id))
    next
  }

  SourceInfoTable <- rio::import(file.path(OSFdata.root, "!!KeyTables", "ML2_SourceInfo - ML2_SourceInfo.csv"))

  # Get the correct slate according to info in ML2.key['study.slate']
  if (ML2.key$study.slate == 1) {
    ML2.df <- rio::import(file.path(OSFdata.root, "!!RawData", "ML2_S1.csv"))
  } else {
    ML2.df <- rio::import(file.path(OSFdata.root, "!!RawData", "ML2_S2.csv"))
  }

  # PREPARE DATA & OUTPUT ----

  # Add a unique ID
  ML2.df$uID <- seq(1, nrow(ML2.df))

  # Get info to create a dataset for the current study
  ML2.in <- get.info(ML2.key, colnames(ML2.df), subset.type)

  # print(ML2.in$stat.params)

  # Get filters for selecting columns
  ML2.id <- get.chain(ML2.in)

  # Apply the df chain to select relevant subset of variables
  ML2.df <- eval(parse(text = paste("ML2.df", ML2.id$df)))

  # Decide which analyses to run on which groups
  toRun <- decide.analysis(
    ML2.key,
    analysis.unique.id,
    analysis.type,
    doAll = TRUE
  )


  if (nrow(ML2.df) <= 0 || (length(toRun$studiess) <= 0 && is.na(toRun$studiess))) {
    print("No tests to run, nothing selected!")
    next
  }

  ML2.sr <- list()

  # Loop over sites in runGroups within a study
  runGroups <- sort(na.exclude(unique(ML2.df[[toRun$ugroup]])))
  print(paste0(length(runGroups), " sites in total for ", analysis.name))

  # Initialize an empty dataframe to store the results
  results_df <- data.frame(
    sites = character(),
    ya = integer(),
    yb = integer(),
    na = integer(), # This is total_a, as per user's fixed logic
    nb = integer(), # This is total_b, as per user's fixed logic
    stringsAsFactors = FALSE
  )

  # START GROUPS ----
  for (g in seq_along(runGroups)) {
    ML2.sr_g <- check.nMin(ML2.df, runGroups[g], Nmin.raw, Nmin.cond, ML2.id, ML2.in)

    if (is.null(ML2.sr_g)) {
      print(paste0("Discarding invalid data in group ", runGroups[[g]]))
      next
    }

    # START ANALYSIS ----------------------------------------
    varfun <- match.fun(varfun_name)
    ML2.var_g <- varfun(ML2.sr_g)

    group_a_name <- levels(ML2.var_g$Condition)[1]
    group_b_name <- levels(ML2.var_g$Condition)[2]

    df_temp <- data.frame(
      condition = ML2.var_g$Condition,
      response = ML2.var_g$Response
    )

    df_temp$response_numeric <- ifelse(df_temp$response == "Yes", 1, 0)

    ya <- sum(df_temp$response_numeric[df_temp$condition == group_a_name])
    na <- sum(df_temp$condition == group_a_name) # This is total_a, as per user's fixed logic

    yb <- sum(df_temp$response_numeric[df_temp$condition == group_b_name])
    nb <- sum(df_temp$condition == group_b_name) # This is total_b, as per user's fixed logic

    site_results <- data.frame(
      sites = runGroups[g],
      ya = ya,
      yb = yb,
      na = na,
      nb = nb,
      stringsAsFactors = FALSE
    )

    results_df <- rbind(results_df, site_results)
  }

    # Add original group names to column names
    if (nrow(results_df) > 0) {
      group_a_name <- levels(ML2.var_g$Condition)[1]
      group_b_name <- levels(ML2.var_g$Condition)[2]
      colnames(results_df) <- c(
        "sites",
        paste0("ya (", group_a_name, ")"),
        paste0("yb (", group_b_name, ")"),
        paste0("na (", group_a_name, ")"),
        paste0("nb (", group_b_name, ")")
      )
    }

  clean_table_filename <- file.path(project.root, "01_tables", "data", paste0(gsub("\\.", "_", analysis.name), "_clean_tables.csv"))
  rio::export(
    results_df,
    clean_table_filename
  )
  print(paste0("csv saved to ", clean_table_filename))
}
