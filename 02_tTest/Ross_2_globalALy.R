library(devtools)
library(plyr)
library(rio)
library(tidyverse)

library(safestats)

#remotes::install_github("AlexanderLyNL/safestats", ref ="futility88")


sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))


# ANALYSIS INFO ----
study.description      <- 'False Consensus 2 (Ross et al., 1977)'
analysis.unique.id     <- 44
analysis.name          <- 'Ross.2'
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
  dplyr::select(2, 7, 19, 25, 805, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 938, 939, 940) %>%
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

# Freq test ------
ML2.var[[g]] <- varfun.Ross.2(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


dat <- ML2.var[[g]]$cleanDataFilter

freqRes <- t.test(variable1 ~ variable2, data = ML2.var[[g]]$cleanDataFilter, var.equal = FALSE)

dat$variable <- dat$variable1
dat$factor <- dat$variable2

studySummary <- dat %>%
  group_by(factor) %>%
  summarise(
    n = n(),
    mean = mean(variable, na.rm = TRUE),
    sd = sd(variable, na.rm=TRUE)
  )

sum(studySummary$n)
studySummary$mean
studySummary$sd

freqRes$statistic
freqRes$p.value
freqRes$statistic*sqrt(sum(studySummary$n)/prod(studySummary$n))

# Alexander ----
dat <- addSources(ML2.var, ML2.df)
dat$variable <- dat$variable1
dat$factor <- dat$variable2
# save(dat, stat.params, file="ross2.RData")

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

dat <- dat[dat$source %in% allSources, ]

if (stat.params$alternative=="two.sided")
  stat.params$alternative <- "twoSided"

# Here -------
alpha <- 0.05
betaFutility <- alpha
deltaMin <- 0.8
varEqual <- stat.params$var.equal
power <- 0.8
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative

designObj <- designSaviT(alpha=alpha, power=power,
                         deltaMin=deltaMin, futility=TRUE,
                         betaFutility=betaFutility,
                         varEqual=varEqual, testType="twoSample",
                         alternative=alternative)

# Scenario 1 ----
res1 <- scenario1T(dat=dat, allSources=allSources, designObj=designObj,
                   nuMin=3, alpha=alpha, betaFutility=betaFutility,
                   nSim=1e3, alternative=alternative)

mean(res1$eValues >= 1/alpha)
mean(res1$eValuesFut <= betaFutility)

res1$nStudiesAlternativeWorstCase
res1$nStudiesFutilityWorstCase

res1$nSamplesAlternativeWorstCase
res1$nSamplesFutilityWorstCase

mean(res1$stopDecision==1)
mean(res1$stopDecision==-1)

mean(res1$nStudies)

mean(res1$logMetaE)
sd(res1$logMetaE)

mean(res1$logMetaEFut)
sd(res1$logMetaEFut)

mean(res1$totalStoppingTimes)
sd(res1$totalStoppingTimes)

# Scenario 2-----
res2 <- scenario2T(dat, allSources, designObj=designObj, seed=1, nSim=1e3)

logMetaE<- rowSums(log(res2$eValues))
mean(logMetaE)
sd(logMetaE)

logMetaEFut <- rowSums(log(res2$eValuesFut))
mean(logMetaEFut)
sd(logMetaEFut)

mean(res2$alternativeProportion)
sd(res2$alternativeProportion)

mean(res2$futilityProportion)
sd(res2$futilityProportion)


mean(res2$totalStoppingTimes)
sd(res2$totalStoppingTimes)

#Scenario 3 ------

res3 <- scenario3T(dat=dat, allSources=allSources, designObj=designObj,
                   alpha=alpha, betaFutility=betaFutility,
                   nuMin=nuMin, nSim=1e3L)

mean(res3$logMetaE)
sd(res3$logMetaE)

mean(res3$logMetaEFut)
sd(res3$logMetaEFut)

mean(res3$alternativeProportion)
sd(res3$alternativeProportion)

mean(res3$futilityProportion)
sd(res3$futilityProportion)

mean(res3$totalStoppingTimes)
sd(res3$totalStoppingTimes)

# save(res1, res2, res3, file="ross2Result.RData")
