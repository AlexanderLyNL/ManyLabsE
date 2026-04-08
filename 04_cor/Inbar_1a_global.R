library(devtools)
library(plyr)
library(rio)
library(tidyverse)


sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# ANALYSIS INFO ----

study.description      <- 'Disgust & Homophobia (Inbar et al., 2009)'
analysis.unique.id     <- 26
analysis.name          <- 'Inbar.1a'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset.type <- "all"
saveAll <- FALSE


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

ML2.id$df

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2,7,165,174,382,383,389,390,391,521,522,523,524,525,526,527,528,529,530,531,532,535,536,537) %>%
  dplyr::filter(is.character(source))


# Decide which analyses to run on which groups
toRun  <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

ML2.df$study.order <- NA
stmp <- strsplit(ML2.df$StudyOrderN,"[|]")

# Correct differences in study names
Stud <- ML2.key$study.name
if(Stud%in%"Tversky"){Stud <- "Tversky.Gati"}
if(Stud%in%"Rottenstreich"){Stud <- "Rottenstrich"}
if(Stud%in%"Ross"&(ML2.key['study.slate'] == 1)){Stud <- "Ross.Slate1"}
if(Stud%in%"Ross"&(ML2.key['study.slate'] == 2)){Stud <- "Ross.Slate2"}
if(Stud%in%"vanLange"){Stud <- "VanLange"}
if(Stud%in%"Giessner"){Stud <- "Geissner"}

ML2.df$study.order <- plyr::laply(seq_along(stmp), function(o){which(grepl(Stud,stmp[[o]]))%00%NA})

ML2.sr       <- list()
ML2.var      <- list()
outputSource <- list()
dataSource   <- list()
raw.df       <- list()
clean.df     <- list()
cleanData    <- list()
testVarEqual <- ML2.in$stat.params$var.equal
g <- 1

# Loop over sites in runGroups within a study
if(analysis.type==1){
  runGroups <- "all"
} else {
  runGroups <- sort(na.exclude(unique(ML2.df[[toRun$ugroup]])))
}

disp(paste(analysis.unique.id, ML2.key$study.analysis,"- START"), header = toupper(ML2.key$study.analysis), footer = FALSE)
cat("\n")


# START GROUPS ----
# HERE -------


# Include only datasets that have N >= Nmin.raw & n.group >= Nmin.cond
listIT     <- FALSE
nMin1      <- FALSE
nMin2      <- FALSE
compN <- compN1 <- compN2 <- 0

gID <- rep(TRUE, nrow(ML2.df))

# Check nMin
if(sum(gID, na.rm=TRUE) >= Nmin.raw){
  nMin1 <- TRUE
  # Get a list containing the data frames to be used in the analysis
  ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
}

# Double-check nMin
compN  <- ML2.sr[[g]]$N
compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
if(any(compN >= Nmin.raw)&(all(compN1>=Nmin.cond, compN2>=Nmin.cond))){nMin2 <- TRUE}


# Freq test ------
ML2.var[[g]] <- varfun.Inbar.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


freqTest <- try.CATCH(with(ML2.var[[g]],
                           cor_test_fisherZ(r1=r1,r2=r2,n1=N[1],n2=N[2], conf.level=stat.params$conf.level, alternative = stat.params$alternative)))

dat <- ML2.var[[g]]$cleanDataFilter



studySummary <- dat %>%
  group_by(factor) %>%
  summarise(
    n = n(),
    mean1 = mean(variable1, na.rm = TRUE),
    sd1 = sd(variable1, na.rm=TRUE),
    mean2 = mean(variable2, na.rm = TRUE),
    sd2 = sd(variable2, na.rm=TRUE),
    cor = cor(variable1, variable2),
    pValue = cor.test(variable1, variable2)$p.value
  )


print(paste0("Different n: ", sum(studySummary$n)))
sum(studySummary$n)
freqTest$value$p.value
freqTest$value$conf.int

diff(-atanh(freqTest$value$estimate))

saviZTestStat(z=sqrt(studySummary$n[1]*studySummary$n[2]/sum(studySummary$n))*diff(-atanh(freqTest$value$estimate)),
              n1=studySummary$n[1], n2=studySummary$n[2],
              parameter=designObj$parameter)

saviFutilityZStat(z=sqrt(studySummary$n[1]*studySummary$n[2]/sum(studySummary$n))*diff(-atanh(freqTest$value$estimate)),
                  n1=studySummary$n[1], n2=studySummary$n[2],
                  parameter=designObj$parameter)




# Alexander ----
dat <- addSources(ML2.var, ML2.df)
# save(dat, stat.params, file="inbar.RData")

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

dat <- dat[dat$source %in% allSources, ]

# Here -------
alpha <- 0.05
betaFutility <- alpha

deltaMin <- 0.7

varEqual <- stat.params$var.equal
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative

set.seed(1234)
designObj <- designSaviZ(meanDiffMin=deltaMin, power=0.8,
                         alpha=alpha, futility=TRUE,
                         testType="twoSample",
                         alternative=alternative)

designObj$testName <- "Correlation"


alphaMeta <- alpha^4
betaFutilityMeta <- alphaMeta

# Scenario 1 ----
res1 <- metaScenario1(dat=dat, allSources=allSources, designObj=designObj,
                      alphaMeta=alphaMeta,
                      betaFutilityMeta=betaFutilityMeta,
                      nSim=1e3)

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
res2 <- metaScenario2(dat=dat, allSources=allSources,
                      designObj=designObj, seed=1, nSim=1e3L,
                      nEffMin=2)

sum(is.infinite(res2$eValues))
sum(is.infinite(res2$eValuesFut))

logMetaE<- rowSums(log(res2$eValues))
mean(logMetaE)
sd(logMetaE)



which(is.infinite(logMetaE))


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
res3 <- metaScenario3(dat=dat, allSources=allSources, designObj=designObj,
                      alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                      nEffMin=2, nuMin=nuMin, nSim=10, wantCi=wantCi, seed=1)

mean(res3$logMetaE >= log(1/alphaMeta))
sd(res3$logMetaE >= log(1/alphaMeta))

mean(res3$logMetaEFut <= log(betaFutilityMeta))
sd(res3$logMetaEFut <= log(betaFutilityMeta))

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

# save(res1, res2, res3, file="inbar1Result.RData")

# Plus -------
alpha <- 0.05
betaFutility <- alpha

varEqual <- stat.params$var.equal
alternative <- "greater"

set.seed(1234)
designObj <- designSaviZ(meanDiffMin=deltaMin, power=0.8,
                         alpha=alpha, futility=TRUE,
                         testType="twoSample",
                         alternative=alternative)

designObj$testName <- "Correlation"

alphaMeta <- alpha/4
betaFutilityMeta <- alphaMeta

# Scenario 1 ----
res1Plus <- metaScenario1(dat=dat, allSources=allSources, designObj=designObj,
                          alphaMeta=alphaMeta,
                          betaFutilityMeta=betaFutilityMeta,
                          nSim=1e3)

mean(res1Plus$eValues >= 1/alphaMeta)
mean(res1Plus$eValuesFut <= betaFutilityMeta)

res1Plus$nStudiesAlternativeWorstCase
res1Plus$nStudiesFutilityWorstCase

res1Plus$nSamplesAlternativeWorstCase
res1Plus$nSamplesFutilityWorstCase

mean(res1Plus$stopDecision==1)
mean(res1Plus$stopDecision==-1)

mean(res1Plus$nStudies)

mean(res1Plus$logMetaE)
sd(res1Plus$logMetaE)

mean(res1Plus$logMetaEFut)
sd(res1Plus$logMetaEFut)

mean(res1Plus$totalStoppingTimes)
sd(res1Plus$totalStoppingTimes)

# Scenario 2-----
res2Plus <- metaScenario2(dat=dat, allSources=allSources,
                          designObj=designObj, seed=1, nSim=1e3L,
                          nEffMin=2)

sum(is.infinite(res2Plus$eValues))
sum(is.infinite(res2Plus$eValuesFut))

logMetaE<- rowSums(log(res2Plus$eValues))
mean(logMetaE)
sd(logMetaE)



logMetaEFut <- rowSums(log(res2Plus$eValuesFut))
mean(logMetaEFut)
sd(logMetaEFut)

mean(res2Plus$alternativeProportion)
sd(res2Plus$alternativeProportion)

mean(res2Plus$futilityProportion)
sd(res2Plus$futilityProportion)


mean(res2Plus$totalStoppingTimes)
sd(res2Plus$totalStoppingTimes)

#Scenario 3 ------
res3Plus <- metaScenario3(dat=dat, allSources=allSources, designObj=designObj,
                          alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                          nEffMin=2, nuMin=nuMin, nSim=1e3L, wantCi=wantCi, seed=1)

mean(res3Plus$logMetaE >= log(1/alphaMeta))
sd(res3Plus$logMetaE >= log(1/alphaMeta))

mean(res3Plus$logMetaEFut <= log(betaFutilityMeta))
sd(res3Plus$logMetaEFut <= log(betaFutilityMeta))

mean(res3Plus$logMetaE)
sd(res3Plus$logMetaE)



mean(res3Plus$logMetaEFut)
sd(res3Plus$logMetaEFut)



mean(res3Plus$alternativeProportion)
sd(res3Plus$alternativeProportion)

mean(res3Plus$futilityProportion)
sd(res3Plus$futilityProportion)

mean(res3Plus$totalStoppingTimes)
sd(res3Plus$totalStoppingTimes)

# save(res1Plus, res2Plus, res3Plus, file="inbar1PlusResult.RData")
