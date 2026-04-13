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

# debugonce(saviTTest)

# bep1$designObjList$knobe$bootObjN1Plan$data


bep1 <- manyLabsMetaScenarios(1, nSim=1e3L)
bep2 <- manyLabsMetaScenarios(2,
                               designObjList=bep1$designObjList,
                               nSim=1e3L)

bep3 <- manyLabsMetaScenarios(3,
                              designObjList=bep1$designObjList,
                              nSim=1e3L)

bep1$resultTable
bep2$resultTable
bep3$resultTable
debugonce(manyLabsMetaScenarios1)

manyLabsMetaScenarios1 <- function(
    metaScenario=1, deltaMinFactor=0.7,
    alternative="greater", nSim=100,
    alpha=0.05, power=0.8, betaFutility=alpha,
    alphaMeta=alpha^4, betaFutilityMeta=alphaMeta,
    wantCi=FALSE, seed=1234, nuMin=3,
    designObjList=NULL, ...)  {

  studyNames <- c("knobe", "ross1", "gray", "ross2",
                  "norenzayan", "hsee", "huang", "kay",
                  "risen", "bauer", "critcher", "giessner",
                  "gati", "zhong", "alter",
                  "zaval", "anderson")

  deltaMinList <- list("knobe"=1.45, "ross1"=0.99, "gray"=0.8,
                       "ross2"=0.8, "norenzayan"=0.35, "hsee"=0.69,
                       "huang"=0.68, "kay"=0.49, "risen"=0.39,
                       "bauer"=0.87, "critcher"=0.3, "giessner"=0.48,
                       "gati"=0.48, "zhong"=1.02, "alter"=0.63, "zaval"=0.31,
                       "anderson"=0.57)

  nStudies <- length(studyNames)



  individualResultList <- vector(mode="list", nStudies)
  names(individualResultList) <- studyNames

  if (is.null(designObjList))
    designObjList <- individualResultList

  if (metaScenario==1) {
    nCol <- 7
  } else if (metaScenario==2) {
    nCol <- 7
  } else if (metaScenario==3) {
    nCol <- 5
  }

  resultTable <- matrix(nrow=length(studyNames), ncol=nCol)

  for (i in seq_along(studyNames)) {
    studyNeem <- studyNames[i]

    if (studyNeem=="gati")
      next()

    testType <- if (studyNeem=="gati") "oneSample" else "twoSample"

    ### Data -------
    # TODO(Alexander): ------
    #     Add the data to the data folder of the package
    #
    load(paste0(myWd, studyNeem, ".RData"))

    dat <- checkUniqueIds(dat)

    if (studyNeem!="gati")
      tempRes <- removeOneConditionSources(dat)

    allSources <- tempRes$allSources
    sampleSize <- tempRes$sampleSize

    dat <- dat[dat$source %in% allSources, ]

    ### Study param setting ----
    varEqual <- stat.params$var.equal

    deltaMin <- deltaMinList[[studyNeem]]
    deltaMin <- deltaMin*deltaMinFactor

    ### designObj ------
    designObj <- designObjList[[studyNeem]]

    if (is.null(designObj)) {
      designObj <- designSaviT(alpha=alpha, power=power,
                               deltaMin=deltaMin, futility=TRUE,
                               betaFutility=betaFutility,
                               varEqual=varEqual, testType=testType,
                               alternative=alternative, seed=seed)

      designObjList[[studyNeem]] <- designObj
    }

    ### analysis ------
    #
    if (metaScenario==1) {
      res <- metaScenario1(
        dat=dat, allSources=allSources,
        designObj=designObj, seed=seed,
        nuMin=nuMin, alphaMeta=alphaMeta,
        betaFutilityMeta=betaFutilityMeta, nSim=nSim)

      resultTable[i, 1] <- mean(res[["logMetaE"]])
      resultTable[i, 2] <- mean(res[["logMetaEFut"]])
      resultTable[i, 3] <- mean(res[["eValues"]] >= 1/alpha)
      resultTable[i, 4] <- mean(res[["eValuesFut"]] <= beta)
      resultTable[i, 5] <- mean(res[["totalStoppingTimes"]])

      resultTable[i, 7] <- dim(dat)[1]
      resultTable[i, 6] <- resultTable[i, 5]/resultTable[i, 7]*100

    } else if (metaScenario==2) {
      res <- metaScenario2(
        dat=dat, allSources=allSources,
        designObj=designObj, seed=seed,
        nuMin=nuMin, nSim=nSim)

      resultTable[i, 1] <- mean(res[["logMetaE"]])
      resultTable[i, 2] <- mean(res[["logMetaEFut"]])
      resultTable[i, 3] <- mean(res[["alternativeProportion"]])
      resultTable[i, 4] <- mean(res[["futilityProportion"]])
      resultTable[i, 5] <- mean(res[["totalStoppingTimes"]])

      resultTable[i, 7] <- dim(dat)[1]
      resultTable[i, 6] <- resultTable[i, 5]/resultTable[i, 7]*100
    } else if (metaScenario==3) {
      res <- metaScenario3(
        dat=dat, allSources=allSources,
        designObj=designObj, alphaMeta=alphaMeta,
        betaFutilityMeta=betaFutilityMeta, nuMin=nuMin,
        nSim=nSim, seed=seed)

      resultTable[i, 1] <- mean(res[["logMetaE"]])
      resultTable[i, 2] <- mean(res[["logMetaEFut"]])
      resultTable[i, 3] <- mean(res[["totalStoppingTimes"]])
      resultTable[i, 5] <- dim(dat)[1]
      resultTable[i, 4] <- resultTable[i, 3]/resultTable[i, 5]*100
    } else {
      stop("Only metaScenario %in% c(1, 2, 3) available")
    }

    individualResultList[[studyNeem]] <- res
  }

  resultTable <- as.data.frame(resultTable)
  rownames(resultTable) <- studyNames

  res <- list(resultTable=resultTable, designObjList=designObjList, individualResultList=individualResultList)

  class(res) <- "saviManyLabs2"
  return(res)
}

# Setting all ----
alpha <- 0.05
betaFutility <- alpha
power <- 0.8
wantCi <- FALSE
alphaMeta <- alpha^4
betaFutilityMeta <- alphaMeta
nSim <- 10



studyNames <- c("knobe", "ross1", "gray", "ross2",
                "norenzayan", "hsee", "huang", "kay",
                "risen", "bauer", "critcher", "giessner",
                "gati", "zhong", "alter",
                "zaval", "anderson")

deltaMinList <- list("knobe"=1.45, "ross1"=0.99, "gray"=0.8,
                     "ross2"=0.8, "norenzayan"=0.35, "hsee"=0.69,
                     "huang"=0.68, "kay"=0.49, "risen"=0.39,
                     "bauer"=0.87, "critcher"=0.3, "giessner"=0.48,
                     "gati"=0.48, "zhong"=1.02, "alter"=0.63, "zaval"=0.31,
                     "anderson"=0.57)

doScenario1 <- function(studyNames, deltaMinFactor=0.7,
                        alternative="greater", deltaMinList,
                        nSim=100, alpha=0.05, betaFutility=alpha,
                        alphaMeta=alpha^4, betaFutilityMeta=alphaMeta)  {

  for (i in seq_along(studyNames)) {
    studyNeem <- studyNames[i]

    ### Data -------
    load(paste0(myWd, studyNeem, "RData"))

    dat <- checkUniqueIds(dat)
    tempRes <- removeOneConditionSources(dat)

    allSources <- tempRes$allSources
    sampleSize <- tempRes$sampleSize

    dat <- dat[dat$source %in% allSources, ]

    ### deltaMin Study setting ----
    deltaMin <- deltaMinList[[dataSetName]]
    deltaMin <- deltaMin*deltaMinFactor

    varEqual <- stat.params$var.equal
  }



  res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                         alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                         nuMin=nuMin, nSim=100, seed=15)




}

# Knobe ------
fileNeem <- paste0(myWd, "knobe.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

varEqual <- stat.params$var.equal



deltaMin <- 1.45*0.7

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative=alternative)


res1 <- metaScenario1(dat=dat, allSources=allSources, designObj=designObjP,
                      nuMin=3, alphaMeta=alphaMeta,
                      betaFutilityMeta=betaFutilityMeta,
                      nSim=1e3)

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=100, seed=15)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion)*100,2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion)*100, 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)



# Ross1 ----
fileNeem <- paste0(myWd, "ross2.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

deltaMin <- 0.99


## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)

round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion)*100,2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion)*100, 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)


# gray -------
fileNeem <- paste0(myWd, "gray.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

## Here -------
deltaMin <- 0.8
varEqual <- stat.params$var.equal
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative

set.seed(1234)
designObj <- designSaviT(alpha=alpha, power=power,
                         deltaMin=deltaMin, futility=TRUE,
                         betaFutility=betaFutility,
                         varEqual=varEqual, testType="twoSample",
                         alternative=alternative)

## Scenario 1 ----
res1 <- metaScenario1(dat=dat, allSources=allSources, designObj=designObj,
                      nuMin=3, alphaMeta=alphaMeta,
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

## Scenario 2-----
res2 <- metaScenario2(dat=dat, allSources=allSources,
                      designObj=designObj, seed=1, nSim=1e3L)

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

##Scenario 3 ------

res3 <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                      alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
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

# save(res1, res2, res3, file="gray1Result.RData")
#
# load(file=paste0(myWd, "gray1Result.RData"))
# res2

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)

# Ross2 ------
fileNeem <- paste0(myWd, "ross2.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

deltaMin <- 0.8

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Norenzayan ------
fileNeem <- paste0(myWd, "norenzayan.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

deltaMin <- 0.35

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Hsee -------
fileNeem <- paste0(myWd, "hsee.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

deltaMin <- 0.69
## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Huang -----------
fileNeem <- paste0(myWd, "huang.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

deltaMin <- 0.68

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Kay --------
fileNeem <- paste0(myWd, "kay.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.49

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Risen ---------
fileNeem <- paste0(myWd, "risen.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.39

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)



# Bauer --------
fileNeem <- paste0(myWd, "bauer.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.87

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Critcher -------
fileNeem <- paste0(myWd, "critcher.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.3

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Giessner -------
deltaMin <- 0.48

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Gati ------

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Zhong -------
fileNeem <- paste0(myWd, "zhong.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 1.02

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Alter ------
fileNeem <- paste0(myWd, "alter.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
# varEqual <- FALSE
deltaMin <- 0.63*0.7


## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Zaval --------
fileNeem <- paste0(myWd, "zaval.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.31

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


round(mean(res3P$logMetaE),2)
# sd(res3P$logMetaE)

round(mean(res3P$logMetaE > log(1/alphaMeta)),2)

round(mean(res3P$logMetaEFut),2)
# sd(res3P$logMetaEFut)

round(mean(res3P$logMetaEFut <= log(betaFutilityMeta)),2)

round(mean(res3P$alternativeProportion),2)
# sd(res3P$alternativeProportion)

round(mean(res3P$futilityProportion), 2)
# sd(res3P$futilityProportion)

round(mean(res3P$totalStoppingTimes), 2)
# sd(res3P$totalStoppingTimes)


# Anderson -------
fileNeem <- paste0(myWd, "anderson.RData")
load(fileNeem)

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize


varEqual <- stat.params$var.equal
varEqual <- FALSE
deltaMin <- 0.57

## Plus -------
alternative <- "greater"
wantCi <- FALSE

set.seed(1234)
designObjP <- designSaviT(alpha=alpha, power=power,
                          deltaMin=deltaMin, futility=TRUE,
                          betaFutility=betaFutility,
                          varEqual=varEqual, testType="twoSample",
                          alternative="greater")

res3P <- metaScenario3(dat=dat, allSources=allSources, designObj=designObjP,
                       alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                       nuMin=nuMin, nSim=nSim, seed=1)


