library(devtools)
library(plyr)
library(rio)
library(tidyverse)
library(safestats)

source(file.path("~", "projects", "manyLabsE","02_tTest","t_test_functions.R"))

# TODO: Set up your directory
project.root <- file.path("~", "projects", "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))

# ANALYSIS INFO ----
study.description      <- 'Intentional Side-Effects (Knobe, 2003)'
analysis.unique.id     <- 79
analysis.name          <- 'Knobe.1'
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
  dplyr::select(2, 7, 180, 188, 805, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 938, 939, 940) %>%
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

ML2.var[[g]] <- varfun.Knobe.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


t.test(variable ~ factor, data = ML2.var[[g]]$cleanDataFilter, var.equal = FALSE)

#write.csv(ML2.df, file = file.path(project.root,"02_tTest","knobeData"), row.names = FALSE)

#############################################################################################################
## extended clean data filter and frequentist analysis

cleanDataFilter <- ML2.var[[g]]$cleanDataFilter
variable <- colnames(cleanDataFilter)[2]
factor <- colnames(cleanDataFilter)[3]
extendedCleanDataFilter <- merge(cleanDataFilter, ML2.df, by = "uID")[,c("uID",variable,factor,"source")] #inner join
extendedCleanDataFilter <- remove_degenerate_sources(extendedCleanDataFilter)
#view(extendedCleanDataFilter)

# THIS IS STUDY-SPECIFFIC
varEqual <- FALSE

frequentist_results <- full_freq_t_test_analysis(extendedCleanDataFilter, var_equal = varEqual)
#view(frequentist_results)

##################################################################################################################
## sequential analysis

# THIS IS STUDY-SPECIFFIC
original_study_estimated_effect_size <- 1.45

esMinFutility <- original_study_estimated_effect_size
deltaMin <- original_study_estimated_effect_size 
alpha <- 0.05
betaFutility <- 0.05

# permute the rows of ECDF to avoid only sampling one group
set.seed(1)
PECDF <- extendedCleanDataFilter[sample(nrow(extendedCleanDataFilter)), ]
sequential_results_list <- full_seq_t_test_analysis(PECDF, alpha, betaFutility, deltaMin, esMinFutility, varEqual = varEqual)
#view(sequential_results_list$sequential_results)
eValueMat  <- sequential_results_list$eValueMat
fValueMat  <- sequential_results_list$fValueMat
metaEType1 <- sequential_results_list$metaEType1
metaFType1 <- sequential_results_list$metaFType1
metaEType2 <- sequential_results_list$metaEType2
metaFType2 <- sequential_results_list$metaFType2
metaEType3 <- sequential_results_list$metaEType3
metaFType3 <- sequential_results_list$metaFType3

worstCaseMetaEType1 <- sequential_results_list$worstCaseMetaEType1
worstCaseMetaFType1 <- sequential_results_list$worstCaseMetaFType1
stoppedMetaEType2 <- sequential_results_list$stoppedMetaEType2
stoppedMetaFType2 <- sequential_results_list$stoppedMetaFType2

#plot(metaEType1, type = 'l')
#lines(metaFType1, col="red")

#plot(metaEType2, type = 'l')
#lines(log(metaFType2), col="red")

#plot(metaEType3, type = 'l')
#lines(metaFType3, col="red")

#plot(stoppedMetaEType2, type = 'l')
#lines(stoppedMetaFType2, col="red")

#plot(worstCaseMetaEType1, type = 'l')
#lines(worstCaseMetaFType1, col="red")

# ANALYSIS OVER MANY PERMUTATIONS OF THE DATA ORDER
avgMetaResultsList <- avg_meta_results(PECDF, alpha, betaFutility, deltaMin, esMinFutility, varEqual, n_permutations = 5, meta_alpha = alpha/n_studies, meta_betaFutility = betaFutility/n_studies)
metaType1StoppingTimes <- avgMetaResultsList$metaType1StoppingTimes
metaType2StoppingTimes <- avgMetaResultsList$metaType2StoppingTimes
metaType3StoppingTimes <- avgMetaResultsList$metaType3StoppingTimes
stoppedMetaType2StoppingTimes <- avgMetaResultsList$stoppedMetaType2StoppingTimes

#print(metaType3StoppingTimes)
#cat(mean(metaType1StoppingTimes), sd(metaType1StoppingTimes))
#hist(metaType1StoppingTimes)

#print(stoppedMetaType2StoppingTimes)
#cat(mean(stoppedMetaType2StoppingTimes), sd(stoppedMetaType2StoppingTimes))

metaEType1Realizations <- avgMetaResultsList$metaEType1Realizations
metaFType1Realizations <- avgMetaResultsList$metaFType1Realizations
metaEType2Realizations <- avgMetaResultsList$metaEType2Realizations
metaFType2Realizations <- avgMetaResultsList$metaFType2Realizations
metaEType3Realizations <- avgMetaResultsList$metaEType3Realizations
metaFType3Realizations <- avgMetaResultsList$metaFType3Realizations

#plot(metaEType1Realizations[[3]], type = 'l')
#lines(metaFType1Realizations[[3]], col="red")

stoppedMetaEType2Realizations <- avgMetaResultsList$stoppedMetaEType2Realizations
stoppedMetaFType2Realizations <- avgMetaResultsList$stoppedMetaFType2Realizations

#plot(stoppedMetaEType2Realizations[[1]], type = 'l')
#lines(stoppedMetaFType2Realizations[[1]], col="red")

all_results <- list(
  extendedCleanDataFilter = extendedCleanDataFilter,
  frequentist_results = frequentist_results,
  sequential_results_list = sequential_results_list,
  avgMetaResultsList = avgMetaResultsList
)

# THIS IS STUDY-SPECIFFIC
#saveRDS(all_results, file = file.path("~", "projects", "manyLabsE","02_tTest","KnobeData.rds"))


