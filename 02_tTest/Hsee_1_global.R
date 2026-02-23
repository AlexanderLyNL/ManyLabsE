library(devtools)
library(plyr)
library(rio)
library(tidyverse)

# TODO: Set up your directory
project.root <- file.path("~", "projects", "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))

# ANALYSIS INFO ----
study.description      <- 'Less is Better (Hsee, 1998)'
analysis.unique.id     <- 61
analysis.name          <- 'Hsee.1'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
# subset -> subset.type to avoid conflicts
subset.type <- "all" # "sites"

# GET LOOKUP TABLES ----
ML2.key <- rio::import(file.path(project.root, "00_data", "ML2_KeyTable.csv"))
ML2.key <- ML2.key[!is.na(ML2.key$unique.id) & ML2.key$unique.id == analysis.unique.id, ]
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

# Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
# Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
ML2.id <- get.chain(ML2.in)

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2, 7, 190, 195, 805, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 938, 939, 940) %>%
  dplyr::filter(is.character(source))

# Decide which analyses to run on which groups
toRun <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

if (nrow(ML2.df) <= 0 || length(toRun$studiess) <= 0) {
  print("No tests to run, nothing selected!")
  stop()
}

# Create a variable indicating the study order for each case
ML2.df$study.order <- NA
stmp <- strsplit(ML2.df$StudyOrderN, "[|]")

# Correct differences in study names
# TODO: a function to fix this
Stud <- ML2.key$study.name

ML2.df$study.order <- plyr::laply(seq_along(stmp), function(o) {
  which(grepl(Stud, stmp[[o]])) %00% NA
})

ML2.sr <- list()
ML2.var <- list()
testVarEqual <- ML2.in$stat.params$var.equal

# analysis.type == 1 no for loop
gID <- rep(TRUE, nrow(ML2.df))
g <- 1

# Check nMin
if (sum(gID, na.rm = TRUE) >= Nmin.raw) {
  nMin1 <- TRUE
  # Get a list containing the data frames to be used in the analysis
  ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
}

# Double-check nMin
if (nMin1) {
  compN <- ML2.sr[[g]]$N
  compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
  compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
  if (any(compN >= Nmin.raw) & (all(compN1 >= Nmin.cond, compN2 >= Nmin.cond))) {
    nMin2 <- TRUE
  }
}

ML2.var[[g]] <- varfun.Hsee.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


t.test(variable ~ factor, data = ML2.var[[g]]$cleanDataFilter, var.equal = FALSE)

#write.csv(ML2.df, file = file.path(project.root,"02_tTest","hseeData"), row.names = FALSE)

# NOTE: we get t-statistic: -34.2, while the original study says 34.2 - the order of the groups is reversed 

extendedCleanDataFilter <- ML2.var[[g]]$cleanDataFilter
extendedCleanDataFilter$study.order <- merge(ML2.df, extendedCleanDataFilter, by = "uID")$study.order

n_studies <- max(extendedCleanDataFilter$study.order)

results <- data.frame(
  study.order = character(n_studies),
  mean_group1 = numeric(n_studies),
  mean_group2 = numeric(n_studies),
  p.value = numeric(n_studies),
  stringsAsFactors = FALSE
)
# the global analysis at index 0
res <- t.test(variable ~ factor, data = extendedCleanDataFilter, var.equal = FALSE)
results$study.order[0] <- 0
results$mean_group1[0] <- res$estimate[1]
results$mean_group2[0] <- res$estimate[2]
results$p.value[0] <- res$p.value 

for (i in seq_len(n_studies)){
  res <- t.test(variable ~ factor, data = extendedCleanDataFilter[extendedCleanDataFilter$study.order == i,], var.equal = FALSE)
  results$study.order[i] <- i
  results$mean_group1[i] <- res$estimate[1]
  results$mean_group2[i] <- res$estimate[2]
  results$p.value[i] <- res$p.value 
}
view(results)
