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

study.description      <- 'Direction & Similarity (Tversky & Gati, 1978)'
analysis.unique.id     <- 83
analysis.name          <- 'Gati.2'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset                 <- 'all'
subset.type <- "all"
saveAll <- FALSE

# GET LOOKUP TABLES ----
ML2.key <- rio::import(file.path(project.root, "00_data", "ML2_KeyTable.csv"))
ML2.key <- ML2.key[!is.na(ML2.key$unique.id) & ML2.key$unique.id == analysis.unique.id, ]
SourceInfoTable <- rio::import(file.path(OSFdata.root, "!!KeyTables", "ML2_SourceInfo - ML2_SourceInfo.csv"))

# Get the correct slate according to info in ML2.key['study.slate']
ML2.df <- rio::import(file.path(OSFdata.root, "!!RawData", "ML2_S2.csv"))


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
  dplyr::select(2,7,709,710,711,712,713,714,715,716,717,718,719,720,721,722,723,724,725,726,727,728,729,735,736,737,738,739,740,741,742,743,744,745,746,747,748,749,750,751,752,753,754,755,805,904,905,906,907,908,909,910,911,912,913,914,915,938,939,940) %>%
  dplyr::filter(is.character(source))

# Decide which analyses to run on which groups
toRun  <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)



# Create a variable indicating the study order for each case
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

# Loop over sites in runGroups within a study
runGroups <- "all"

disp(paste(analysis.unique.id, ML2.key$study.analysis,"- START"), header = toupper(ML2.key$study.analysis), footer = FALSE)
cat("\n")


# START GROUPS ----


# Include only datasets that have N >= Nmin.raw & n.group >= Nmin.cond
listIT     <- FALSE
nMin1      <- FALSE
nMin2      <- FALSE
compN <- compN1 <- compN2 <- 0

gID <- rep(TRUE, nrow(ML2.df))

# Check nMin

nMin1 <- TRUE
# Get a list containing the data frames to be used in the analysis
ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)


# Double-check nMin
compN  <- ML2.sr[[g]]$N
compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
if(any(compN >= Nmin.raw)&(all(compN1>=Nmin.cond, compN2>=Nmin.cond))){nMin2 <- TRUE}


# Freq test ------

ML2.var[[g]] <- varfun.Gati.2(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params
stat.test   <- try.CATCH(with(ML2.var[[g]],
                              t.test(x = Asymmetry, mu=0, conf.level=stat.params$conf.level, var.equal = stat.params$var.equal, alternative = stat.params$alternative)))

# t.test(x = Asymmetry, mu=0, conf.level=stat.params$conf.level, var.equal = stat.params$var.equal, alternative = stat.params$alternative))

dat <- ML2.var[[g]]$cleanDataFilter
dim(dat)


freqRes <- stat.test$value
freqRes$statistic
freqRes$p.value
freqRes$conf.int


# Alexander ----
dat <- addSources(ML2.var, ML2.df)

# save(dat, stat.params, file="gati.RData")


dat <- checkUniqueIds(dat)
# tempRes <- removeOneConditionSources(dat)

allSources <- unique(dat$source)

dat <- dat[dat$source %in% allSources, ]

# if (stat.params$alternative=="two.sided")
#   stat.params$alternative <- "twoSided"


# Here -------
alpha <- 0.05
betaMinEffi <- alpha
deltaMin <- 0.48
varEqual <- NULL
power <- 0.8
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative
wantCi <- FALSE

alphaMeta <- alpha^4
betaMinEffiMeta <- alphaMeta

set.seed(1234)
designObj <- designSaviT(alpha=alpha, power=power,
                         deltaMin=deltaMin, minEffiTest=TRUE,
                         betaMinEffi=betaMinEffi,
                         varEqual=varEqual, testType="oneSample",
                         alternative=alternative)
dir()
load("~/Desktop/git/manyLabsE/02_tTest/gati.RData")


# Scenario 1 ----
res1 <- metaScenario1(dat=dat, allSources=allSources, designObj=designObj,
                      nuMin=3, alphaMeta=alphaMeta,
                      betaMinEffiMeta=betaMinEffiMeta,
                      nSim=1e3)

mean(res1$eValues >= 1/alpha)
mean(res1$eValuesMinEffi <= betaMinEffi)

res1$nStudiesAlternativeWorstCase
res1$nStudiesMinEffiWorstCase

res1$nSamplesAlternativeWorstCase
res1$nSamplesMinEffiWorstCase

mean(res1$stopDecision==1)
mean(res1$stopDecision==-1)

mean(res1$nStudies)

mean(res1$logMetaE)
sd(res1$logMetaE)

mean(res1$logMetaEMinEffi)
sd(res1$logMetaEMinEffi)

mean(res1$totalStoppingTimes)
sd(res1$totalStoppingTimes)

# Scenario 2-----
res2 <- metaScenario2(dat=dat, allSources=allSources,
                      designObj=designObj, nSim=10)



logMetaE<- rowSums(log(res2$eValues))
mean(logMetaE)
sd(logMetaE)

logMetaEMinEffi <- rowSums(log(res2$eValuesMinEffi))
mean(logMetaEMinEffi)
sd(logMetaEMinEffi)

mean(res2$alternativeProportion)
sd(res2$alternativeProportion)

mean(res2$minEffiProportion)
sd(res2$minEffiProportion)


mean(res2$totalStoppingTimes)
sd(res2$totalStoppingTimes)

#Scenario 3 ------
nuMin <- 3
res3 <- metaScenario3(dat=dat, allSources=allSources, designObj=designObj,
                      alphaMeta=alphaMeta, betaMinEffiMeta=betaMinEffiMeta,
                      nuMin=nuMin, nSim=10)



mean(res3$logMetaE)
sd(res3$logMetaE)

mean(res3$logMetaEMinEffi)
sd(res3$logMetaEMinEffi)

mean(res3$alternativeProportion)
sd(res3$alternativeProportion)

mean(res3$minEffiProportion)
sd(res3$minEffiProportion)

mean(res3$totalStoppingTimes)
sd(res3$totalStoppingTimes)

# save(res1, res2, res3, file="gati2Result.RData")

