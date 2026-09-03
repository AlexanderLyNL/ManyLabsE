# Utility ----
addSources <- function(ML2.var, ML2.df) {
  dat <- ML2.var[[1]]$cleanDataFilter
  sourceColumn <- character(length = dim(dat)[1])

  for (i in seq_along(sourceColumn)) {
    iets <- ML2.df[ML2.df$uID == dat$uID[i], ]
    sourceColumn[i] <- iets$source
  }

  dat$source <- sourceColumn
  return(dat)
}

addUniqueIds <- function(ML2.var, ML2.df) {
  return(addSources(ML2.var, ML2.df))
}

checkUniqueIds <- function(dat) {
  allIds <- unique(dat$uID)
  nTotal <- length(allIds)

  if (length(dat$uID) == nTotal) {
    return(dat)
  } else {
    sadf <- table(dat$uID)
    doubleIds <- as.integer(names(which(sadf>1)))
    iets <- setdiff(allIds, doubleIds)

    dat <- dat[dat$uID %in% iets, ]

    for (id in doubleIds)
      cat("Removed double ID: ", id, "\n")

    return(dat)
  }
}

checkEnoughDataInTable <- function(dat, Nmin.raw=30, Nmin.cond=15) {
  allSources <- unique(dat$sites)
  nSources <- length(allSources)

  allBadIndices <- integer(0)

  for (i in seq_along(allSources)) {

    someDat <- dat[dat$sites==allSources[i], ]

    na <- someDat[["na"]]
    nb <- someDat[["nb"]]
    ya <- someDat[["ya"]]
    yb <- someDat[["yb"]]
    faila <- na-ya
    failb <- nb-yb
    nTotal <- na+nb

    if (nTotal < Nmin.raw) {
      allBadIndices <- c(allBadIndices, i)
      next()
    }

    if (any(c(ya, yb, faila, failb) < Nmin.cond)) {
      allBadIndices <- c(allBadIndices, i)
      next()
    }
  }

  if (length(allBadIndices) > 0) {
    for (i in seq_along(allBadIndices))
      cat("Not enough total or per condition samples: ", allSources[allBadIndices[i]], "\n")

    allSources <- allSources[-allBadIndices]
  }

  res <- list("allSources"=allSources)
  return(res)
}

removeOneConditionSources <- function(dat) {
  allSources <- unique(dat$source)

  nSources <- length(allSources)
  factorLevels <- unique(dat$factor)
  nFactors <- length(factorLevels)

  dataLeveled <- vector("list", length=nFactors)
  badIndices <- vector("list", length=nFactors)
  sampleSize <- matrix(nrow=nSources, ncol=nFactors)

  for (i in 1:nFactors)
    dataLeveled[[i]] <- dat %>% dplyr::filter(factor==factorLevels[i])

  for (i in 1:nSources) {
    someSource <- allSources[i]

    for (j in seq_along(dataLeveled)) {
      someDat <- dataLeveled[[j]]
      someDat <- someDat[someDat$source==someSource, ]
      sampleSize[i, j] <- dim(someDat)[1]
    }
  }

  allBadIndices <- integer(0)

  for (i in 1:nFactors) {
    badIndices[[i]] <- which(sampleSize[, i]==0)
    allBadIndices <- union(allBadIndices, badIndices[[i]])
    # allBadIndices <-
  }

  if (length(allBadIndices) > 0) {
    for (i in seq_along(allBadIndices))
      cat("Only one condition for source: ", allSources[allBadIndices[i]], "\n")

    allSources <- allSources[-allBadIndices]
    sampleSize <- sampleSize[-allBadIndices, ]
  }

  res <- list("allSources"=allSources, "sampleSize"=sampleSize)
  return(res)
}


checkTwoSample <- function(x, y, nMin=1) {
  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  if (length(x) >= nMin && length(y) >= nMin)
    return(TRUE)

  if (is.null(x))
    return(FALSE)

  if (is.null(y))
    return(FALSE)

  if (length(x) < nMin)
    return(FALSE)

  if (length(y) < nMin)
    return(FALSE)

  if (is.na(x))
    return(FALSE)

  if (is.na(y))
    return(FALSE)

  return(TRUE)
}

checkOneSample <- function(x, nMin=1) {
  x <- x[!is.na(x)]

  if (length(x) >= nMin)
    return(TRUE)

  if (is.null(x))
    return(FALSE)

  if (length(x) < nMin)
    return(FALSE)

  if (is.na(x))
    return(FALSE)

  return(TRUE)
}


computeWorstCaseScenario1 <- function(
    res, alphaMeta=0.05, betaMinEffiMeta=alphaMeta,
    seed=NULL, nSim=1e3L, ...) {

  eValues <- sort(res$eValues)
  eValuesMinEffi <- sort(res$eValuesMinEffi, decreasing=TRUE)

  nSources <- length(eValuesMinEffi)

  nStudiesAlternativeWorstCase <- min(which(cumsum(log(eValues)) >= log(1/alphaMeta)))
  nStudiesMinEffiWorstCase <- min(which(cumsum(log(eValuesMinEffi)) <= log(betaMinEffiMeta)))

  if (is.infinite(nStudiesAlternativeWorstCase)) {
    nSamplesAlternativeWorstCase <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValues)
    indexStudiesNeeded <- someOrder[1:nStudiesAlternativeWorstCase]

    nSamplesAlternativeWorstCase <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  if (is.infinite(nStudiesMinEffiWorstCase)) {
    nSamplesMinEffiWorstCase <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValuesMinEffi)
    indexStudiesNeeded <- someOrder[1:nStudiesMinEffiWorstCase]

    nSamplesMinEffiWorstCase <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  stopDecision <- nStudies <- totalStoppingTimes <- integer(nSim)
  logMetaE <- logMetaEMinEffi <- integer(nSim)

  set.seed(seed)
  for (i in 1:nSim) {
    someOrder <- sample(nSources, nSources)

    tempEValues <- eValues[someOrder]
    tempEValuesMinEffi <- eValuesMinEffi[someOrder]

    logMetaETemp <- cumsum(log(tempEValues))
    logMetaEMinEffiTemp <- cumsum(log(tempEValuesMinEffi))

    tauForAlt <- min(which(logMetaETemp >= log(1/alphaMeta)))
    tauForMinEffi <- min(which(logMetaEMinEffiTemp <= log(betaMinEffiMeta)))

    if (tauForMinEffi < tauForAlt)
      stopDecision[i] <- -1

    if (tauForAlt < tauForMinEffi)
      stopDecision[i] <- 1

    tauRace <- min(tauForAlt, tauForMinEffi)

    stopIndex <- nStudies[i] <- min(tauRace, nSources)

    logMetaE[i] <- logMetaETemp[stopIndex]
    logMetaEMinEffi[i] <- logMetaEMinEffiTemp[stopIndex]

    indexNeededStudies <- someOrder[1:stopIndex]

    totalStoppingTimes[i] <- sum(res$n1Vec[indexNeededStudies])+sum(res$n2Vec[indexNeededStudies])
  }

  res <- list("nStudiesAlternativeWorstCase"=nStudiesAlternativeWorstCase,
              "nStudiesMinEffiWorstCase"=nStudiesMinEffiWorstCase,
              "nSamplesAlternativeWorstCase"=nSamplesAlternativeWorstCase,
              "nSamplesMinEffiWorstCase"=nSamplesMinEffiWorstCase,
              "stopDecision"=stopDecision,
              "logMetaE"=logMetaE,
              "logMetaEMinEffi"=logMetaEMinEffi,
              "nStudies"=nStudies,
              "totalStoppingTimes"=totalStoppingTimes)

  return(res)
}

# Meta scenarios -------

manyLabsMetaScenarios <- function(
    scenarioNumber=1, deltaMinFactor=1,
    alternative="greater", nSim=100,
    alpha=0.05, power=0.8, betaMinEffi=alpha,
    alphaMeta=alpha^4, betaMinEffiMeta=alphaMeta,
    wantCi=FALSE, seed=1234, nuMin=3,
    designObjList=NULL, analysisType=c("2x2", "tTest", "zTest"), ...)  {

  analysisType <- match.arg(analysisType)

  if (analysisType=="2x2" && scenarioNumber != 1)
    stop("For 2x2 tables on scenarioNumber 1 is available")

  studyNames <- switch(
    analysisType,
    "2x2"=c("hauser1", "tversky", "hauser2", "rottenstreich"),
    "tTest"=c("knobe", "ross1", "gray", "ross2",
              "norenzayan", "hsee", "huang", "kay",
              "risen", "bauer", "critcher", "giessner",
              "gati", "zhong", "alter",
              "zaval", "anderson"),
    "zTest"=c("graham", "inbar", "vanLange", "schwarz", "shafir")
  )

  deltaMinList <- switch(
    analysisType,
    "2x2"=list("hauser1"=2.5*pi/sqrt(3), "tversky"=log(4.96),
               "hauser2"=0.34*pi/sqrt(3), "rottenstreich"=0.74*pi/sqrt(3)),
    "tTest"=list("knobe"=1.45, "ross1"=0.99, "gray"=0.8,
               "ross2"=0.8, "norenzayan"=0.35, "hsee"=0.69,
               "huang"=0.68, "kay"=0.49, "risen"=0.39,
               "bauer"=0.87, "critcher"=0.3, "giessner"=0.48,
               "gati"=0.48, "zhong"=1.02, "alter"=0.63, "zaval"=0.31,
               "anderson"=0.57),
    "zTest"=list("graham"=2*0.25/sqrt(1-0.25^2),
                 "inbar"=0.7, "vanLange"=2*0.25/sqrt(1-0.25^2),
                 "schwarz"=0.48, "shafir"=sqrt(4*.095^2/(1-.095^2)))
  )

  nStudies <- length(studyNames)

  individualResultList <- vector(mode="list", nStudies)
  names(individualResultList) <- studyNames

  if (is.null(designObjList))
    designObjList <- individualResultList

  if (scenarioNumber==1) {
    nCol <- 7
    nColFull <- 11
  } else if (scenarioNumber==2) {
    nCol <- 7
    nColFull <- 13
  } else if (scenarioNumber==3) {
    nCol <- 5
    nColFull <- 13
  }

  resultTable <- matrix(nrow=length(studyNames), ncol=nCol)

  resultTableFull <- matrix(nrow=length(studyNames), ncol=nColFull)

  oneSampleNames <- c("hauser1", "tversky", "hauser2", "rottenstreich",
                      "gati", "graham", "vanLange", "shafir")

  correlationList <- c("graham", "inbar", "vanLange", "schwarz")

  for (i in seq_along(studyNames)) {
    studyNeem <- studyNames[i]

    testType <- if (studyNeem %in% oneSampleNames) "oneSample" else "twoSample"

    ### Data -------
    # TODO(Alexander): ------
    #     Add the data to the data folder of the package
    #
    load(paste0(myWd, studyNeem, ".RData"))

    dat <- checkUniqueIds(dat)

    # TODO(Alexander): Check for other analyses as well
    #
    if (studyNeem %in% oneSampleNames) {
      tempRes <- list(allSources=unique(dat[["source"]]))
    } else {
      tempRes <- removeOneConditionSources(dat)
    }

    allSources <- tempRes[["allSources"]]
    sampleSize <- tempRes[["sampleSize"]]

    dat <- dat[dat[["source"]] %in% allSources, ]

    ### Study param setting ----
    if (analysisType=="tTest")
      varEqual <- stat.params[["var.equal"]]

    deltaMin <- deltaMinList[[studyNeem]]
    deltaMin <- deltaMin*deltaMinFactor

    ### designObj ------
    designObj <- designObjList[[studyNeem]]

    if (is.null(designObj)) {
      if (analysisType=="2x2") {
        designObj <- list(
          "esMin"=deltaMin, "minEffiTestResult"=list("parameter"=deltaMin),
          "alternative"=alternative,
          "testName"="2x2")
      } else if (analysisType=="tTest") {
        designObj <- designSaviT(
          alpha=alpha, power=power,
          deltaMin=deltaMin, minEffiTest=TRUE,
          betaMinEffi=betaMinEffi,
          varEqual=varEqual, testType=testType,
          alternative=alternative, seed=seed)
      } else if (analysisType=="zTest") {
        designObj <- designSaviZ(
          alpha=alpha, power=power,
          meanDiffMin=deltaMin, minEffiTest=TRUE,
          betaMinEffi=betaMinEffi,
          testType=testType,
          alternative=alternative, seed=seed)

        if (studyNeem %in% correlationList) {
          designObj[["testName"]] <- "Correlation"
        } else if (studyNeem=="shafir") {
          designObj[["testName"]] <- "Binomial"
        }
      }

      designObjList[[studyNeem]] <- designObj
    }

    ### analysis ------
    #
    if (scenarioNumber==1) {
      res <- metaScenario1(
        dat=dat, allSources=allSources,
        designObj=designObj, seed=seed,
        nuMin=nuMin, alphaMeta=alphaMeta,
        betaMinEffiMeta=betaMinEffiMeta, nSim=nSim)

      # Table --
      resultTable[i, 7] <- sum(res[["n1Vec"]])+sum(res[["n2Vec"]])

      resultTable[i, 1] <- mean(res[["logMetaE"]])
      resultTable[i, 2] <- mean(res[["logMetaEMinEffi"]])
      resultTable[i, 3] <- mean(res[["eValues"]] >= 1/alpha)*100
      resultTable[i, 4] <- mean(res[["eValuesMinEffi"]] <= betaMinEffi)*100
      resultTable[i, 5] <- mean(res[["totalStoppingTimes"]])
      resultTable[i, 6] <- (resultTable[i, 5]/resultTable[i, 7])*100

      # Full table --
      resultTableFull[i, 9] <- sum(res[["n1Vec"]])+sum(res[["n2Vec"]])

      resultTableFull[i, 1] <- mean(res[["logMetaE"]])
      resultTableFull[i, 2] <- sd(res[["logMetaE"]])

      resultTableFull[i, 3] <- mean(res[["logMetaEMinEffi"]])
      resultTableFull[i, 4] <- sd(res[["logMetaEMinEffi"]])

      resultTableFull[i, 5] <- mean(res[["eValues"]] >= 1/alpha)*100

      resultTableFull[i, 6] <- mean(res[["eValuesMinEffi"]] <= betaMinEffi)*100

      resultTableFull[i, 7] <- mean(res[["totalStoppingTimes"]])

      resultTableFull[i, 8] <- (resultTableFull[i, 7]/resultTableFull[i, 9])*100

      resultTableFull[i, 10] <- sd(res[["totalStoppingTimes"]])


      resultTableFull[i, 11] <- resultTableFull[i, 10]/resultTableFull[i, 9]*100

    } else if (scenarioNumber==2) {
      res <- metaScenario2(
        dat=dat, allSources=allSources,
        designObj=designObj, seed=seed,
        nuMin=nuMin, nSim=nSim)

      # Table --
      logMetaE <- rowSums(log(res[["eValues"]]))
      logMetaEMinEffi <- rowSums(log(res[["eValuesMinEffi"]]))

      resultTable[i, 7] <- dim(dat)[1]

      resultTable[i, 1] <- mean(logMetaE)
      resultTable[i, 2] <- mean(logMetaEMinEffi)
      resultTable[i, 3] <- mean(res[["alternativeProportion"]])*100
      resultTable[i, 4] <- mean(res[["minEffiProportion"]])*100
      resultTable[i, 5] <- mean(res[["totalStoppingTimes"]])
      resultTable[i, 6] <- (resultTable[i, 5]/resultTable[i, 7])*100

      # Table full --
      resultTableFull[i, 11] <- dim(dat)[1]

      resultTableFull[i, 1] <- mean(logMetaE)
      resultTableFull[i, 2] <- sd(logMetaE)
      resultTableFull[i, 3] <- mean(logMetaEMinEffi)
      resultTableFull[i, 4] <- sd(logMetaEMinEffi)
      resultTableFull[i, 5] <- mean(res[["alternativeProportion"]])*100
      resultTableFull[i, 6] <- sd(res[["alternativeProportion"]])*100
      resultTableFull[i, 7] <- mean(res[["minEffiProportion"]])*100
      resultTableFull[i, 8] <- sd(res[["minEffiProportion"]])*100
      resultTableFull[i, 9] <- mean(res[["totalStoppingTimes"]])

      resultTableFull[i, 10] <- (resultTableFull[i, 9]/resultTableFull[i, 11])*100

      resultTableFull[i, 12] <- sd(res[["totalStoppingTimes"]])
      resultTableFull[i, 13] <- resultTableFull[i, 12]/resultTableFull[i, 11]*100
    } else if (scenarioNumber==3) {
      res <- metaScenario3(
        dat=dat, allSources=allSources,
        designObj=designObj, alphaMeta=alphaMeta,
        betaMinEffiMeta=betaMinEffiMeta, nuMin=nuMin,
        nSim=nSim, seed=seed)

      # Table
      #
      resultTable[i, 5] <- dim(dat)[1]

      resultTable[i, 1] <- mean(res[["logMetaE"]])
      resultTable[i, 2] <- mean(res[["logMetaEMinEffi"]])
      resultTable[i, 3] <- mean(res[["totalStoppingTimes"]])
      resultTable[i, 4] <- (resultTable[i, 3]/resultTable[i, 5])*100

      # Table full
      #
      resultTableFull[i, 11] <- dim(dat)[1]

      resultTableFull[i, 1] <- mean(res[["logMetaE"]])
      resultTableFull[i, 2] <- sd(res[["logMetaE"]])
      resultTableFull[i, 3] <- mean(res[["logMetaEMinEffi"]])
      resultTableFull[i, 4] <- sd(res[["logMetaEMinEffi"]])
      resultTableFull[i, 5] <- mean(res[["alternativeProportion"]])*100
      resultTableFull[i, 6] <- sd(res[["alternativeProportion"]])*100
      resultTableFull[i, 7] <- mean(res[["minEffiProportion"]])*100
      resultTableFull[i, 8] <- sd(res[["minEffiProportion"]])*100
      resultTableFull[i, 9] <- mean(res[["totalStoppingTimes"]])

      resultTableFull[i, 10] <- (resultTableFull[i, 9]/resultTableFull[i, 11])*100

      resultTableFull[i, 12] <- sd(res[["totalStoppingTimes"]])
      resultTableFull[i, 13] <- resultTableFull[i, 12]/resultTableFull[i, 11]*100
    } else {
      stop("Only scenarioNumber %in% c(1, 2, 3) available")
    }

    individualResultList[[studyNeem]] <- res
  }

  resultTable <- as.data.frame(resultTable)
  rownames(resultTable) <- studyNames

  if (scenarioNumber %in% 1:2) {
    colnames(resultTable) <- c("logMetaE", "logMetaEMinEffi",
                               "Reject H0", "Reject H1",
                               "nStop", "% of", "nTotal")
  } else if (scenarioNumber==3) {
    colnames(resultTable) <- c("logMetaE", "logMetaEMinEffi",
                               "nStop", "% of", "nTotal")
  }

  resultTableFull <- as.data.frame(resultTableFull)
  rownames(resultTableFull) <- studyNames

  if (scenarioNumber==1) {
    colnames(resultTableFull) <- c("logMetaE", "sd(logMetaE)", "logMetaEMinEffi",
                                   "sd(logMetaEMinEffi)", "Reject H0", "Reject H1",
                                   "nStop", "% of", "nTotal", "sd(nStop)", "sd %")
  } else if (scenarioNumber %in% 2:3) {
    colnames(resultTableFull) <- c("logMetaE", "sd(logMetaE)", "logMetaEMinEffi",
                                   "sd(logMetaEMinEffi)", "Reject H0", "sd(Reject H0)",
                                   "Reject H1", "sd(Reject H1)",
                                   "nStop", "% of", "nTotal",
                                   "sd(nStop)", "sd %")
  }

  res <- list(resultTable=resultTable, resultTableFull=resultTableFull, designObjList=designObjList, individualResultList=individualResultList)

  class(res) <- "saviManyLabs2"
  return(res)
}

metaScenario1 <- function(dat, allSources, designObj,
                          nuMin=3, wantCi=FALSE,
                          alphaMeta=0.05, betaMinEffiMeta=alphaMeta,
                          seed=NULL, nSim=1e3L,
                          alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesMinEffi <- pValues <- numeric(nSources)
  n1Vec <- n2Vec <- integer(nSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat[["source"]]==allSources[i], ]

    if (designObj[["testName"]]=="T-Test") {
      tempRes <- scenario1TTestHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "wantCi"=wantCi, "nuMin"=nuMin)
    } else if (designObj[["testName"]]=="Binomial") {
      tempRes <- scenario1BinomialHelp(
        "someDat"=someDat, "designObj"=designObj)
    } else if (designObj[["testName"]]=="2x2") {
      tempRes <- scenario12x2Help(
        "someDat"=someDat, "designObj"=designObj)
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- scenario1CorHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels)
    }

    n1Vec[i] <- tempRes[["n1"]]
    n2Vec[i] <- tempRes[["n2"]]
    somePValue <- tempRes[["pValue"]]
    pValues[i] <- if (is.null(somePValue)) 1 else somePValue
    eValues[i] <- tempRes[["eValue"]]
    eValuesMinEffi[i] <- tempRes[["eValueMinEffi"]]
  }

  tempRes <- list("eValues"=eValues, "eValuesMinEffi"=eValuesMinEffi,
                  "pValues"=pValues,
                  "n1Vec"=n1Vec, "n2Vec"=n2Vec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaMinEffiMeta"=betaMinEffiMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}

metaScenario2 <- function(dat, allSources, designObj, alphaMeta=0.05,
                          betaMinEffiMeta=alphaMeta, nuMin=3, nSim=1e2L,
                          nMax=NULL, seed=NULL, wantCi=FALSE, nEffMin=2) {

  alternative <- designObj[["alternative"]]
  nSources <- length(allSources)

  nSamples <- eValues <- eValuesMinEffi <- matrix(nrow=nSim, ncol=nSources)

  if (designObj[["testType"]]=="twoSample") {
    factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)
  } else if (designObj[["testType"]]=="oneSample") {
    factorLevels <- NULL
  }

  seedNext <- NULL

  for (i in 1:length(allSources)) {
    someDat <- dat[dat[["source"]]==allSources[i], ]

    if (!is.null(seed)) seedNext <- seed+i

    if (designObj[["testName"]]=="T-Test") {
      tempRes <- scenario2TTestHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "wantCi"=wantCi,
        "nuMin"=nuMin, "nSim"=nSim, "seed"=seedNext)
    } else if (designObj[["testName"]]=="Binomial") {
      tempRes <- scenario2BinomHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "wantCi"=wantCi,
        "nuMin"=nuMin, "nSim"=nSim, "seed"=seedNext)
    } else if (designObj[["testName"]]=="2x2") {
      stop("2x2 not yet done")
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- scenario2CorHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "nSim"=nSim,
        "nEffMin"=nEffMin, "seed"=seedNext)
    }

    nSamples[, i] <- tempRes[["nSamples"]]
    eValues[, i] <- tempRes[["eValues"]]
    eValuesMinEffi[, i] <- tempRes[["eValuesMinEffi"]]
  }

  alternativeProportion <- minEffiProportion <- numeric(length=nSim)

  for (i in 1:nSim) {
    alternativeProportion[i] <- mean(eValues[i, ] >= 1/designObj[["alpha"]])
    minEffiProportion[i] <- mean(eValuesMinEffi[i, ] <= designObj$minEffiTestResult$beta)
  }

  totalStoppingTimes <- rowSums(nSamples)

  res <- list("nSamples"=nSamples, "eValues"=eValues, "eValuesMinEffi"=eValuesMinEffi,
              "alternativeProportion"=alternativeProportion,
              "minEffiProportion"=minEffiProportion,
              "totalStoppingTimes"=totalStoppingTimes)
  return(res)
}


metaScenario3 <- function(dat, allSources, designObj, alphaMeta=0.05,
                          betaMinEffiMeta=alphaMeta, nuMin=3, nSim=1e3L,
                          nMax=NULL, seed=NULL, wantCi=FALSE,
                          nPlanLimit=FALSE, nEffMin=2) {

  nTotal <- length(unique(dat[["uID"]]))
  nSources <- length(allSources)

  logMetaE <- logMetaEMinEffi <- numeric(nSim)

  alternativeProportion <- minEffiProportion <- totalStoppingTimes <-
    integer(nSim)

  nSamples <- nStopDecision <- matrix(nrow=nSim, ncol=nSources)
  logEValues <- logEValuesMinEffi <- matrix(nrow=nSim, ncol=nSources)

  seedNext <- NULL

  for (i in 1:nSim) {

    if (!is.null(seed)) seedNext <- seed+i

    if (designObj[["testName"]]=="T-Test") {
      tempRes <- computeScenario3TOneSim(
        dat=dat, allSources=allSources, designObj=designObj,
        alphaMeta=alphaMeta, betaMinEffiMeta=betaMinEffiMeta,
        nuMin=nuMin, nSim=nSim, seed=seedNext,
        wantCi=wantCi, nPlanLimit=nPlanLimit)
    } else if (designObj[["testName"]]=="Binomial") {
      tempRes <- computeScenario3BinomialOneSim(
        dat=dat, allSources=allSources, designObj=designObj,
        alphaMeta=alphaMeta, betaMinEffiMeta=betaMinEffiMeta,
        nuMin=nuMin, nSim=nSim, wantCi=wantCi, nPlanLimit=nPlanLimit,
        seed=seedNext)
    } else if (designObj[["testName"]]=="2x2") {
      stop("2x2 not yet done")
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- computeScenario3CorOneSim(
        dat=dat, allSources=allSources, designObj=designObj,
        alphaMeta=alphaMeta, betaMinEffiMeta=betaMinEffiMeta,
        nuMin=nuMin, nSim=nSim,
        wantCi=wantCi, nPlanLimit=nPlanLimit, nEffMin=nEffMin,
        seed=seedNext)
    }

    logMetaE[i] <- tempRes[["logMetaE"]]
    logMetaEMinEffi[i] <- tempRes[["logMetaEMinEffi"]]
    logEValues[i, ] <- tempRes[["logEValues"]]
    logEValuesMinEffi[i, ] <- tempRes[["logEValuesMinEffi"]]
    nSamples[i, ] <- tempRes[["nSamples"]]
    nStopDecision[i, ] <- tempRes[["stopDecision"]]

    totalStoppingTimes[i] <- sum(tempRes[["nSamples"]])
    alternativeProportion[i] <- mean(tempRes[["stopDecision"]]==1)
    minEffiProportion[i] <- mean(tempRes[["stopDecision"]]==-1)
  }

  res <- list("logMetaE"=logMetaE, "logMetaEMinEffi"=logMetaEMinEffi,
              "logEValues"=logEValues, "logEValuesMinEffi"=logEValuesMinEffi,
              "nSamples"=nSamples, "nStopDecision"=nStopDecision,
              "totalStoppingTimes"=totalStoppingTimes,
              "alternativeProportion"=alternativeProportion,
              "minEffiProportion"=minEffiProportion)
  return(res)
}

# T-test ----
computeEValuesT <- function(x, y, designObj, nuMin=3, h0=0) {
  res <- list(eValue=NULL, eValueMinEffi=NULL, n1=NULL, n2=NULL)

  sumStats <- computeZTSumStats(
    "x"=x, "y"=y, "sequential"=FALSE,
    "varEqual"=designObj[["varEqual"]], "paired"=FALSE,
    "testType"=designObj[["testType"]])

  list2env(sumStats, envir=environment())

  if (nu <= nuMin || is.na(sdObs)) {
    tStat <- 0
  } else {
    tStat <- try(sqrt(nEff)*(meanObs - h0)/sdObs)
  }

  if (is.na(tStat) && sdObs==0 && meanObs-h0==0)
    tStat <- 0

  if (is.na(tStat))
    stop("Data error: Could not compute the t-statistic")

  names(tStat) <- "t"

  ### Compute: eValue ----
  #
  testResult <- suppressWarnings(
    saviTTestStatNEffNu("t"=tStat, "nEff"=nEff, "nu"=nu,
                        "parameter"=designObj[["parameter"]],
                        "alternative"=designObj[["alternative"]],
                        "paired"=FALSE,
                        "tDensity"=FALSE,
                        "nuMin"=nuMin, "eType"=designObj[["eType"]])
  )



  res[["eValue"]] <- unname(testResult[["eValue"]])


  if (designObj[["minEffiTest"]]) {
    testResultMinEffi <- suppressWarnings(
      saviMinEffiTStatNEffNu("t"=tStat, "nEff"=nEff, "nu"=nu,
                              "parameter"=designObj[["minEffiTestResult"]][["parameter"]],
                              "alternative"=designObj[["alternative"]], "paired"=FALSE,
                              "nuMin"=nuMin)
    )
  }

  res[["eValueMinEffi"]] <- unname(testResultMinEffi[["eValue"]])

  res[["n1"]] <- sumStats[["n1"]]

  res[["n2"]] <- if (is.null(n2)) 0 else sumStats[["n2"]]

  return(res)
}


scenario1TTestHelp <- function(someDat, designObj, factorLevels=NULL,
                               wantCi=FALSE, nuMin=3) {

  res <- list(eValue=NULL, eValueMinEffi=NULL, n1=NULL, n2=NULL, pValue=NULL)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
    y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable

    # Remove non-available entries
    x <- x[!is.na(x)]
    n1 <- length(x)

    y <- y[!is.na(y)]
    n2 <- length(y)
  } else if (designObj[["testType"]]=="oneSample") {
    x <- someDat[["outcome"]]

    if (is.null(x))
      stop("No 'outcome' column in the data set.")

    # Remove non-available entries
    x <- x[!is.na(x)]
    n1 <- length(x)

    y <- NULL
    n2 <- 0
  }

  alternativeOld <- switch(designObj[["alternative"]],
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  tempResult <- t.test(x[1:n1], y[1:n2],
                       alternative=alternativeOld,
                       var.equal=designObj[["varEqual"]])

  res[["pValue"]] <- tempResult[["p.value"]]

  tempRes <- computeEValuesT(
    "x"=x, "y"=y, "designObj"=designObj,
    "nuMin"=nuMin)

  res[["eValue"]] <- tempRes[["eValue"]]
  res[["eValueMinEffi"]] <- tempRes[["eValueMinEffi"]]

  res[["n1"]] <- tempRes[["n1"]]
  res[["n2"]] <- tempRes[["n2"]]

  return(res)
}

scenario2TTestHelp <- function(
    someDat, designObj, factorLevels,
    wantCi=FALSE, nuMin=3, seed=NULL,
    nSim=1e3L, ...) {

  res <- list(nSamples=NULL, eValue=NULL, eValueMinEffi=NULL)

  nSamples <- integer(nSim)
  eValues <- eValuesMinEffi <- numeric(nSim)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
    y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable
  } else if (designObj[["testType"]]=="oneSample") {
    x <- someDat[["outcome"]]
    y <- NULL
  }


  # Remove non-available entries
  x <- x[!is.na(x)]
  n1 <- length(x)

  y <- y[!is.na(y)]
  n2 <- length(y)

  if (!is.null(designObj[["nPlan"]])) {
    n1 <- min(n1, designObj[["nPlan"]][1])
    n2 <- min(n2, designObj[["nPlan"]][2], na.rm=TRUE)
  }

  nParticipants <- n1+n2

  seedNext <- NULL

  for (k in 1:nSim) {
    if (!is.null(seed)) seedNext <- seed + k

    tempRes <- tTestRandomOrder(
      "x"=x, "y"=y, "n1"=n1, "n2"=n2,
      "designObj"=designObj, "nuMin"=nuMin,
      "wantCi"=wantCi, "seed"=seedNext
    )

    nSamples[k] <- tempRes[["nSamples"]]
    eValues[k] <- tempRes[["eValue"]]
    eValuesMinEffi[k] <- tempRes[["eValueMinEffi"]]
  }

  res[["nSamples"]] <- nSamples
  res[["eValues"]] <- eValues
  res[["eValuesMinEffi"]] <- eValuesMinEffi

  return(res)
}

tTestRandomOrder <- function(
    x, y, n1, n2, designObj,
    nuMin=3, wantCi=FALSE, nMax=NULL,
    seed=NULL, ...) {

  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nParticipants <- n1+n2

  xRun <- numeric(0)
  yRun <- numeric(0)

  set.seed(seed)
  someOrder <- sample(nParticipants, nParticipants)
  xTemp <- sample(x, length(x))
  yTemp <- sample(y, length(y))

  xTemp <- xTemp[1:n1]
  yTemp <- yTemp[1:n2]

  totalVar <- c(xTemp, yTemp)

  nMax <- if (is.null(nMax)) nParticipants else min(nParticipants, nMax)

  for (j in seq_along(someOrder)) {
    partId <- someOrder[j]

    if (partId <= n1) {
      xRun <- c(xRun, totalVar[partId])
    } else {
      yRun <- c(yRun, totalVar[partId])
    }

    if (designObj[["testType"]]=="twoSample"){
      someCheck <- checkTwoSample(xRun, yRun)
    } else {
      someCheck <- if (length(xRun) > 1) TRUE else FALSE
    }

    if (someCheck) {

      if (designObj[["testType"]]=="oneSample")
        yRun <- NULL

      tempRes <- computeEValuesT(
        "x"=xRun, "y"=yRun,
        "designObj"=designObj, "nuMin"=nuMin)

      eNow <- tempRes[["eValue"]]
      eMinEffiNow <- tempRes[["eValueMinEffi"]]

      if (eNow >= 1/alpha || eMinEffiNow <= betaMinEffi ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueMinEffi"=eMinEffiNow)
        return(res)
      }
    }
  }
}


computeScenario3TOneSim <- function(
    dat, allSources, designObj, alphaMeta=0.05,
    betaMinEffiMeta=alphaMeta, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nSources <- length(allSources)

  if (designObj[["testType"]]=="twoSample") {
    factorLevels <- if (is.ordered(dat[["factor"]])) levels(dat[["factor"]]) else unique(dat[["factor"]])
  } else if (designObj[["testType"]]=="oneSample") {
    factorLevels <- NULL
  }

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources

  for (neem in allSources)
    sourceDataTracker[[neem]] <- list(x=NULL, y=NULL)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEMinEffiTracker <- logETracker

  nTotal <- length(dat[["uID"]])

  # set.seed(seed)
  someOrder <- sample(unique(dat[["uID"]]), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEMinEffiNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat$uID==someId), ]
    someSource <- someRow$source

    nSamples[[someSource]] <- nSamples[[someSource]] + 1

    sourceDataTemp <- sourceDataTracker[[someSource]]

    # Retrieve old values from state
    #
    x <- sourceDataTemp$x
    y <- sourceDataTemp$y

    # Skip if sample size limit is reached within trial
    #
    if (designObj[["testType"]]=="twoSample" &&
        nPlanLimit && length(x) >= designObj$nPlan[1] &&
        length(y) >= designObj$nPlan[2])
      next()

    if (designObj[["testType"]]=="oneSample" &&
        nPlanLimit && length(x) >= designObj$nPlan[1])
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()

    if (designObj[["testType"]]=="twoSample") {
      if (someRow$factor==factorLevels[1]) {
        sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
      } else if (someRow$factor==factorLevels[2]) {
        sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
      }
    } else if (designObj[["testType"]]=="oneSample") {
      sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$outcome)
    }

    # TODO(Alexander): Hier tTestRandom....

    if (designObj[["testType"]]=="twoSample") {
      someCheck <- checkTwoSample(x, y)
    } else {
      someCheck <- if (length(x) > 1) TRUE else FALSE
    }

    if (someCheck) {
      logEValueOld <- logETracker[[someSource]]
      logEValueMinEffiOld <- logEMinEffiTracker[[someSource]]

      tempRes <- computeEValuesT(
        "x"=x, "y"=y,
        "designObj"=designObj, "nuMin"=nuMin)

      # tempRes <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
      #                      "sequential"=FALSE, "wantCi"=wantCi)

      logEValueNow <- logETracker[[someSource]] <-
        log(tempRes$eValue)
      logEValueMinEffiNow <- logEMinEffiTracker[[someSource]] <-
        log(tempRes$eValueMinEffi)

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueMinEffiNow <= log(betaMinEffi))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEMinEffiAdd <- logEValueMinEffiNow - logEValueMinEffiOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEMinEffiNow <- logMetaEMinEffiNow+logMetaEMinEffiAdd

      if (logMetaENow >= log(1/alphaMeta) || logMetaEMinEffiNow <= log(betaMinEffiMeta)) {
        break
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEMinEffi=logMetaEMinEffiNow,
              logEValuesMinEffi=logEMinEffiTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}



# Tables --------
scenario1Table <- function(dat, allSources, designObj,
                           alpha=0.05, betaMinEffi=alpha,
                           seed=NULL, nSim=1e3L) {

  alternative <- designObj$alternative

  nSources <- length(allSources)

  eValues <- eValuesMinEffi <- pValues <- numeric(nSources)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$sites==allSources[i], ]

    ## Data ---
    na <- someDat[["na"]]
    nb <- someDat[["nb"]]
    ya <- someDat[["ya"]]
    yb <- someDat[["yb"]]

    n1 <- ya+yb

    # TODO(Alexander): Frequen,
    #
    # alternativeOld <- switch(alternative,
    #                          "twoSided"="two.sided",
    #                          "greater"="greater",
    #                          "less"="less")
    #
    # tempResult <- t.test(x[1:n1], y[1:n2],
    #                      alternative=alternativeOld,
    #                      var.equal=varEqual)
    #
    # pValues[i] <- tempResult$p.value

    tempRes <- saviTwoPropConditionalStat(
      ya=ya, na=na, nb=nb, n1=n1,
      logOddsRatio=designObj$esMin,
      alternative=alternative,
      eType="eGauss")

    eValues[i] <- tempRes$eValue

    tempRes <- saviMinEffiTwoPropConditionalStat(
      ya=ya, na=na, nb=nb, n1=n1,
      logOddsRatio=designObj[["minEffiTestResult"]][["parameter"]],
      alternative=alternative)

    eValuesMinEffi[i] <- tempRes$eValue
  }

  tempRes <- list("eValues"=eValues, "eValuesMinEffi"=eValuesMinEffi,
                  "pValues"=pValues)

  print("Broken off function")
  return(tempRes)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaMinEffi"=betaMinEffiMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}

scenario12x2Help <- function(someDat, designObj) {
  res <- list(eValue=NULL, eValueMinEffi=NULL, n1=NULL, n2=0, pValue=NULL)

  ## Data ---
  na <- someDat[["na"]]
  nb <- someDat[["nb"]]
  ya <- someDat[["ya"]]
  yb <- someDat[["yb"]]

  n1 <- ya+yb
  nTotal <- na+nb

  if (length(n1)==0 || nTotal==0)
    return(list(eValue=1, eValueMinEffi=1, n1=0, n2=0, pValue=1))

  alternative <- designObj[["alternative"]]

  # Frequentist analysis
  #
  alternativeOld <- switch(alternative,
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  someMatrix <- matrix(ncol=2, nrow=2)

  someMatrix[1, 1] <- ya
  someMatrix[1, 2] <- yb
  someMatrix[2, 1] <- na - ya
  someMatrix[2, 2] <- nb - yb

  freqTest <- try(fisher.test(x=someMatrix,
                          alternative=alternative))

  if (isTryError(freqTest))
    freqTest <- list("p.value"=NULL)

  res[["pValue"]] <- freqTest[["p.value"]]


  tempRes <- saviTwoPropConditionalStat(
    ya=ya, na=na, nb=nb, n1=n1,
    logOddsRatio=designObj[["esMin"]],
    alternative=designObj[["alternative"]],
    eType=designObj[["eGauss"]])

  res[["eValue"]] <- tempRes[["eValue"]]

  tempRes <- saviMinEffiTwoPropConditionalStat(
    ya=ya, na=na, nb=nb, n1=n1,
    logOddsRatio=designObj[["minEffiTestResult"]][["parameter"]],
    alternative=designObj[["alternative"]])

  res[["eValueMinEffi"]] <- tempRes[["eValue"]]

  res[["n1"]] <- nTotal

  return(res)
}

saviTwoPropConditionalStat <- function(ya, na, nb, n1, logOddsRatio,
                                       eType=c("eGauss", "grow"),
                                       alternative=c("twoSided", "greater", "less")) {
  res <- list("eValue"=1, "eValueApproxError"=NULL)

  alternative <- match.arg(alternative)
  eType <- match.arg(eType)

  marginalNull <- BiasedUrn::dFNCHypergeo(x=ya, m1=na, m2=nb, n=n1,
                                          odds=1)

  if (eType=="grow") {
    if (alternative %in% c("twoSided", "greater")) {
      marginalPlus <- BiasedUrn::dFNCHypergeo(x=ya, m1=na, m2=nb, n=n1,
                                        odds=exp(abs(logOddsRatio)))
    }

    if (alternative %in% c("twoSided", "less")) {
      marginalMin <- BiasedUrn::dFNCHypergeo(x=ya, m1=na, m2=nb, n=n1,
                                       odds=exp(-abs(logOddsRatio)))
    }

    eValue <- switch(alternative,
                     "twoSided"=1/2*marginalPlus+1/2*marginalMin,
                     "greater"=marginalPlus,
                     "less"=marginalMin)/marginalNull

    res[["eValue"]] <- eValue

    return(res)
  }


  # Numerical versions
  #
  integrandBound <- switch(alternative,
                           "twoSided"=c(-Inf, Inf),
                           "greater"=c(0, Inf),
                           "less"=c(-Inf, 0))

  if (eType=="eGauss") {
    integrand <- function(z) {
      priorValue <- gaussPrior(z, sd=abs(logOddsRatio), alternative=alternative)

      if (priorValue==0)
        return(0)

      res <- BiasedUrn::dFNCHypergeo(x=ya, m1=na, m2=nb, n=n1, odds=exp(z))*priorValue

      return(res)
    }

    integrand2 <- Vectorize(integrand, vectorize.args = "z")

    tempRes <- integrate(integrand2, "lower"=integrandBound[1],
                         "upper"=integrandBound[2])

    res[["eValue"]] <- tempRes[["value"]]/marginalNull
    res[["eValueApproxError"]] <- tempRes[["abs.error"]]

    return(res)
  }
}

saviMinEffiTwoPropConditionalStat <- function(
    ya, na, nb, n1, logOddsRatio,
    alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  if (alternative %in% c("twoSided", "greater")) {
    sPlus0 <- saviTwoPropConditionalStat(
      ya=ya, na=na, nb=nb, n1=n1, logOddsRatio=logOddsRatio,
      alternative="greater", eType="grow")
  }

  if (alternative %in% c("twoSided", "less")) {
    sMin0 <- saviTwoPropConditionalStat(
      ya=ya, na=na, nb=nb, n1=n1, logOddsRatio=logOddsRatio,
      alternative="less", eType="grow")
  }

  eValue <- switch(alternative,
                   "greater"=sPlus0$eValue,
                   "less"=sMin0$eValue,
                   "twoSided"=max(sPlus0$eValue, sMin0$eValue))

  res <- list(eValue=eValue)

  return(res)
}

gaussPrior <- function(z, sd=1, alternative) {
  normalisationConstant <- switch(alternative,
                                  "twoSided"=1,
                                  "greater"=1/2,
                                  "less"=1/2)

  if (alternative=="greater" && z < 0)
    return(0)

  if (alternative=="less" && z > 0)
    return(0)

  return(dnorm(z, mean=0, sd=sd)/normalisationConstant)
}

# Binomial Z --------
scenario1BinomialHelp <- function(someDat, designObj, factorLevels=NULL) {

  res <- list(eValue=NULL, eValueMinEffi=NULL, n1=NULL, n2=NULL, pValue=NULL)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    stop("Not yet implemented")
  } else if (designObj[["testType"]]=="oneSample") {
    count <- as.integer(someDat[["variable1"]]=="Parent B")

    datAward <- someDat[someDat[["variable2"]]=="Award", ]
    datDeny <- someDat[someDat[["variable2"]]=="Deny", ]

    meansAward <- datAward %>%
      summarise(mean=mean(count, na.rm=TRUE))
    meansDeny <- datDeny %>%
      summarise(mean=mean(count, na.rm=TRUE))

    n1 <- length(count)

    allMeans <- (meansAward[["mean"]]+meansDeny[["mean"]])/2
    se <- sqrt(0.5*(1-0.5) / length(count))
    zScore <- (allMeans-0.5) / se
  }

  alternativeOld <- switch(designObj[["alternative"]],
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  res[["pValue"]] <- 1-pnorm(abs(zScore))

  tempRes <-  saviZTestStat(z=zScore, n1=n1,
                            parameter=designObj[["parameter"]],
                            alternative=designObj[["alternative"]],
                            sigma=designObj[["sigma"]],
                            eType=designObj[["eType"]])

  tempResMinEffi <- saviMinEffiZStat(z=zScore, n1=n1,
                               parameter=designObj[["minEffiTestResult"]][["parameter"]],
                               alternative=designObj[["alternative"]],
                               sigma=designObj[["sigma"]])

  res[["eValue"]] <- tempRes[["eValue"]]
  res[["eValueMinEffi"]] <- tempResMinEffi[["eValue"]]

  res[["n1"]] <- n1
  res[["n2"]] <- 0

  return(res)
}


scenario2BinomHelp <- function(
    someDat, designObj, factorLevels,
    seed=NULL, nSim=1e3L, ...) {

  res <- list(nSamples=NULL, eValue=NULL, eValueMinEffi=NULL)

  nSamples <- integer(nSim)
  eValues <- eValuesMinEffi <- numeric(nSim)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    stop("Not yet implemented")
  } else if (designObj[["testType"]]=="oneSample") {
    count <- as.integer(someDat[["variable1"]]=="Parent B")
  }


  # Remove non-available entries
  someDat <- someDat[!is.na(someDat[["variable1"]]), ]
  someDat <- someDat[!is.na(someDat[["variable2"]]), ]

  n1 <- length(someDat[["variable1"]])
  n2 <- 0


  if (!is.null(designObj[["nPlan"]]))
    n1 <- min(n1, designObj[["nPlan"]][1])

  nParticipants <- n1+n2

  seedNext <- NULL

  for (k in 1:nSim) {
    if (!is.null(seed)) seedNext <- seed + k
    tempRes <- try(binomialTestRandomOrder(
      "x"=someDat, "n1"=n1,
      "designObj"=designObj,
      "alpha"=alpha, "betaMinEffi"=betaMinEffi
    ))

    nSamples[k] <- tempRes[["nSamples"]]
    eValues[k] <- tempRes[["eValue"]]
    eValuesMinEffi[k] <- tempRes[["eValueMinEffi"]]
  }

  res[["nSamples"]] <- nSamples
  res[["eValues"]] <- eValues
  res[["eValuesMinEffi"]] <- eValuesMinEffi

  return(res)
}

binomialTestRandomOrder <- function(
    x, n1, n2=0, designObj, # nuMin=3,
    alpha=0.05,
    betaMinEffi=alpha, # wantCi=FALSE
    seed=NULL, nMax=NULL) {


  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nParticipants <- n1+n2

  set.seed(seed)
  someOrder <- sample(nParticipants, nParticipants)

  xTemp <- x[sample(nrow(x)), ]

  nMax <- if (is.null(nMax)) nParticipants else min(nParticipants, nMax)

  for (j in seq_along(someOrder)) {
    partId <- someOrder[1:j]

    xRun <- xTemp[partId, ]
    count <- as.integer(xRun[["variable1"]]=="Parent B")

    datAward <- xRun[xRun[["variable2"]]=="Award", ]
    datDeny <- xRun[xRun[["variable2"]]=="Deny", ]

    if (length(datAward[["uID"]]) > 2 && length(datDeny[["uID"]]) > 2) {
      meansAward <- datAward %>%
        summarise(mean=mean(count, na.rm=TRUE))
      meansDeny <- datDeny %>%
        summarise(mean=mean(count, na.rm=TRUE))

      allMeans <- (meansAward$mean+meansDeny$mean)/2
      se <- sqrt(0.5*(1-0.5) / length(count))
      zScore <- (allMeans-0.5) / se

      tempRes <-  saviZTestStat(z=zScore, n1=length(count),
                                parameter=designObj[["parameter"]],
                                alternative=designObj[["alternative"]],
                                sigma=designObj[["sigma"]],
                                eType=designObj[["eType"]])
      eNow <- tempRes[["eValue"]]

      tempRes <- saviMinEffiZStat(z=zScore, n1=length(count),
                                   parameter=designObj[["minEffiTestResult"]][["parameter"]],
                                   alternative=designObj[["alternative"]],
                                   sigma=designObj[["sigma"]])
      eMinEffiNow <- tempRes[["eValue"]]

      if (eNow >= 1/alpha || eMinEffiNow <= betaMinEffi ||
          j==nMax) {

        res <- list("nSamples"=j, "eValue"=eNow, "eValueMinEffi"=eMinEffiNow)
        return(res)
      }
    }
  }
}

computeScenario3BinomialOneSim <- function(
    dat, allSources, designObj, alphaMeta=0.05,
    betaMinEffiMeta=alphaMeta, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nSources <- length(allSources)

  ## Data ---
  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources

  for (neem in allSources)
    sourceDataTracker[[neem]] <- list(x=NULL)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEMinEffiTracker <- logETracker

  # Alexander: Perhaps remove na here
  # Remove non-available entries
  dat <- dat[!is.na(dat[["variable1"]]), ]
  dat <- dat[!is.na(dat[["variable2"]]), ]

  nTotal <- length(dat[["uID"]])

  # set.seed(seed)
  someOrder <- sample(unique(dat[["uID"]]), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEMinEffiNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat$uID==someId), ]
    someSource <- someRow[["source"]]

    nSamples[[someSource]] <- nSamples[[someSource]] + 1

    sourceDataTemp <- sourceDataTracker[[someSource]]

    # Retrieve old values from state
    #
    x <- sourceDataTemp[["x"]]

    # Skip if sample size limit is reached within trial
    #
    if (designObj[["testType"]]=="oneSample" &&
        nPlanLimit && nrow(x) >= designObj[["nPlan"]][1])
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()

    if (designObj[["testType"]]=="twoSample") {
      stop("Not yet implemented")
    } else if (designObj[["testType"]]=="oneSample") {

      sourceDataTracker[[someSource]]$x <- x <- rbind(x, someRow)
    }

    # TODO(Alexander): Hier tTestRandom....

    if (designObj[["testType"]]=="twoSample") {
      someCheck <- checkTwoSample(x, y)
    } else {
      someCheck <- if (nrow(x) > 1) TRUE else FALSE
    }

    if (someCheck) {
      logEValueOld <- logETracker[[someSource]]
      logEValueMinEffiOld <- logEMinEffiTracker[[someSource]]

      count <- as.integer(x[["variable1"]]=="Parent B")

      datAward <- x[x[["variable2"]]=="Award", ]
      datDeny <- x[x[["variable2"]]=="Deny", ]

      if (length(datAward[["uID"]]) > 2 && length(datDeny[["uID"]]) > 2) {
        meansAward <- datAward %>%
          summarise(mean=mean(count, na.rm=TRUE))
        meansDeny <- datDeny %>%
          summarise(mean=mean(count, na.rm=TRUE))

        allMeans <- (meansAward[["mean"]]+meansDeny[["mean"]])/2
        se <- sqrt(0.5*(1-0.5) / length(count))
        zScore <- (allMeans-0.5) / se

        tempRes <-  saviZTestStat(z=zScore, n1=length(count),
                                  parameter=designObj[["parameter"]],
                                  alternative=designObj[["alternative"]],
                                  sigma=designObj[["sigma"]],
                                  eType=designObj[["eType"]])

        tempResMinEffi <- saviMinEffiZStat(z=zScore, n1=length(count),
                                     parameter=designObj[["minEffiTestResult"]][["parameter"]],
                                     alternative=designObj[["alternative"]],
                                     sigma=designObj[["sigma"]])
      } else {
        tempRes <- tempResMinEffi <- list(eValue=1)
      }


      logEValueNow <- logETracker[[someSource]] <-
        log(tempRes[["eValue"]])
      logEValueMinEffiNow <- logEMinEffiTracker[[someSource]] <-
        log(tempResMinEffi[["eValue"]])

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueMinEffiNow <= log(betaMinEffi))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEMinEffiAdd <- logEValueMinEffiNow - logEValueMinEffiOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEMinEffiNow <- logMetaEMinEffiNow+logMetaEMinEffiAdd

      if (logMetaENow >= log(1/alphaMeta) || logMetaEMinEffiNow <= log(betaMinEffiMeta)) {
        break
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEMinEffi=logMetaEMinEffiNow,
              logEValuesMinEffi=logEMinEffiTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}



# Cor-to-Z --------
scenario1CorHelp <- function(someDat, designObj, factorLevels=NULL) {

  res <- list(eValue=NULL, eValueMinEffi=NULL, n1=NULL, n2=NULL, pValue=NULL)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    dat1 <- someDat[which(someDat$factor==factorLevels[1]), ]
    dat2 <- someDat[which(someDat$factor==factorLevels[2]), ]

    n1 <- dim(dat1)[1]-3
    n2 <- dim(dat2)[1]-3

    nEff <- n1*n2/(n1+n2)

    corrie1 <- cor(dat1[["variable1"]], dat1[["variable2"]])
    corrie2 <- cor(dat2[["variable1"]], dat2[["variable2"]])

    zScore <- sqrt(nEff)*(atanh(corrie1)- atanh(corrie2))
  } else if (designObj[["testType"]]=="oneSample") {
    n1 <- dim(someDat)[1]-3
    n2 <- NULL

    nEff <- n1

    corrie <- cor(someDat[["variable1"]], someDat[["variable2"]])
    zScore <- sqrt(nEff)*(atanh(corrie))
  }

  alternativeOld <- switch(designObj[["alternative"]],
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  res[["pValue"]] <- 1-pnorm(abs(zScore))

  tempRes <- saviZTestStat(z=zScore, n1=n1, n2=n2,
                           parameter=designObj[["parameter"]],
                           eType=designObj[["eType"]], alternative=designObj[["alternative"]])

  res[["eValue"]] <- tempRes[["eValue"]]

  tempRes <- saviMinEffiZStat(z=zScore, n1=n1, n2=n2,
                               parameter=designObj[["minEffiTestResult"]][["parameter"]],
                               alternative=designObj[["alternative"]])

  res[["eValueMinEffi"]] <- tempRes[["eValue"]]
  res[["n1"]] <- n1
  res[["n2"]] <- if (!is.null(n2)) n2 else 0

  return(res)
}


scenario2CorHelp <- function(someDat, designObj, factorLevels,
                             nSim=1e3L, seed=NULL, nEffMin=2) {

  res <- list(nSamples=NULL, eValue=NULL, eValueMinEffi=NULL)

  nSamples <- integer(nSim)
  eValues <- eValuesMinEffi <- numeric(nSim)

  ## Data ---
  if (designObj[["testType"]]=="twoSample") {
    dat1 <- someDat[which(someDat$factor==factorLevels[1]), ]
    dat2 <- someDat[which(someDat$factor==factorLevels[2]), ]

    n1 <- dim(dat1)[1]
    n2 <- dim(dat2)[1]
  } else if (designObj[["testType"]]=="oneSample") {
    dat1 <- someDat
    dat2 <- NULL

    n1 <- dim(dat1)[1]
    n2 <- 0
  }

  if (!is.null(designObj[["nPlan"]])) {
    n1 <- min(n1, designObj[["nPlan"]][1])
    n2 <- min(n2, designObj[["nPlan"]][2], na.rm=TRUE)
  }

  nParticipants <- n1+n2

  seedNext <- NULL

  for (k in 1:nSim) {
    if (!is.null(seed)) seedNext <- seed + k

    tempRes <- corTestRandomOrder(dat1, dat2, n1, n2, designObj,
                                  seed=seedNext, nEffMin=nEffMin)

    nSamples[k] <- tempRes[["nSamples"]]
    eValues[k] <- tempRes[["eValue"]]
    eValuesMinEffi[k] <- tempRes[["eValueMinEffi"]]
  }

  res[["nSamples"]] <- nSamples
  res[["eValues"]] <- eValues
  res[["eValuesMinEffi"]] <- eValuesMinEffi

  return(res)
}

corTestRandomOrder <- function(
    dat1, dat2, n1, n2, designObj,
    seed=NULL, nMax=NULL,
    nEffMin=2) {

  res <- list(nSamples=NULL, eValue=NULL, eValueMinEffi=NULL)

  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nParticipants <- n1+n2

  dat1XRun <- numeric(0)
  dat1YRun <- numeric(0)

  if (designObj[["testType"]]=="oneSample") {
    dat2XRun <- NULL
    dat2YRun <- NULL
  } else if (designObj[["testType"]]=="twoSample") {
    dat2XRun <- numeric(0)
    dat2YRun <- numeric(0)
  }


  set.seed(seed)
  someOrder <- sample(nParticipants, nParticipants)

  dat1Order <- sample(n1)
  dat2Order <- sample(n2)

  dat1Temp <- dat1[dat1Order, ]
  dat2Temp <- dat2[dat2Order, ]

  # totalVar <- c(xTemp, yTemp)
  totalXVar <- c(dat1Temp$variable1, dat2Temp$variable1)
  totalYVar <- c(dat1Temp$variable2, dat2Temp$variable2)

  nMax <- if (is.null(nMax)) nParticipants else min(nParticipants, nMax)

  for (j in seq_along(someOrder)) {
    partId <- someOrder[j]

    if (partId <= n1) {
      dat1XRun <- c(dat1XRun, totalXVar[partId])
      dat1YRun <- c(dat1YRun, totalYVar[partId])
    } else {
      dat2XRun <- c(dat2XRun, totalXVar[partId])
      dat2YRun <- c(dat2YRun, totalYVar[partId])
    }

    if (designObj[["testType"]]=="twoSample") {
      someCheck <- checkTwoSample(dat1XRun, dat2XRun, nMin=4)
    } else if (designObj[["testType"]]=="oneSample") {
      someCheck <- checkOneSample(dat1XRun, nMin=4)
    }

    if (someCheck) {
      tempCor <- computeCor(x1=dat1XRun, y1=dat1YRun, x2=dat2XRun, y2=dat2YRun)

      r1 <- tempCor[["r1"]]
      r2 <- tempCor[["r2"]]

      n1Now <- length(dat1XRun)-3
      n2Now <- length(dat2XRun)-3

      if (designObj[["testType"]]=="oneSample") {
        nEffNow <- n1Now
        n2Now <- NULL

        zScore <-  sqrt(nEffNow)*(atanh(r1))
      } else if (designObj[["testType"]]=="twoSample") {
        nEffNow <- n1Now*n2Now/(n1Now+n2Now)

        zScore <-  sqrt(nEffNow)*(atanh(r1)-atanh(r2))
      }

      if (is.na(zScore) && nEffNow < nEffMin) {
        eNow <- 1
        eMinEffiNow <- 1
      } else {
        tempRes <- saviZTestStat(z=zScore, n1=n1Now, n2=n2Now,
                                 parameter=designObj[["parameter"]],
                                 alternative=designObj[["alternative"]],
                                 eType=designObj[["eType"]])

        eNow <- tempRes[["eValue"]]

        tempRes <- try(saviMinEffiZStat(z=zScore, n1=n1Now, n2=n2Now,
                                     parameter=designObj[["minEffiTestResult"]][["parameter"]],
                                     alternative=designObj[["alternative"]]))

        eMinEffiNow <- tempRes[["eValue"]]
      }

      if (eNow >= 1/alpha || eMinEffiNow <= betaMinEffi ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueMinEffi"=eMinEffiNow)
        return(res)
      }
    }
  }
}


computeScenario3CorOneSim <- function(
    dat, allSources, designObj, alphaMeta=0.05,
    betaMinEffiMeta=alphaMeta, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE, nEffMin=2,
    nPlanLimit=TRUE) {

  alpha <- designObj[["alpha"]]
  betaMinEffi <- designObj[["minEffiTestResult"]][["beta"]]

  nSources <- length(allSources)

  if (designObj[["testType"]]=="twoSample") {
    factorLevels <- if (is.ordered(dat[["factor"]])) levels(dat[["factor"]]) else unique(dat[["factor"]])
  } else if (designObj[["testType"]]=="oneSample") {
    factorLevels <- NULL
  }

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources

  for (neem in allSources)
    sourceDataTracker[[neem]] <- list(dat1X=NULL, dat1Y=NULL, dat2X=NULL, dat2Y=NULL)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEMinEffiTracker <- logETracker

  nTotal <- length(dat[["uID"]])

  set.seed(seed)
  someOrder <- sample(unique(dat[["uID"]]), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEMinEffiNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat[["uID"]]==someId), ]
    someSource <- someRow[["source"]]

    nSamples[[someSource]] <- nSamples[[someSource]] + 1

    sourceDataTemp <- sourceDataTracker[[someSource]]

    # Retrieve old values from state
    #
    dat1X <- sourceDataTemp[["dat1X"]]
    dat1Y <- sourceDataTemp[["dat1Y"]]

    dat2X <- sourceDataTemp[["dat2X"]]
    dat2Y <- sourceDataTemp[["dat2Y"]]


    # Skip if sample size limit is reached within trial
    #
    if (designObj[["testType"]]=="twoSample" &&
        nPlanLimit && length(dat1X) >= designObj[["nPlan"]][1] &&
        length(dat2X) >= designObj[["nPlan"]][2]) {
      next()
    }

    if (designObj[["testType"]]=="oneSample" &&
        nPlanLimit && length(dat1X) >= designObj[["nPlan"]][1])
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()

    if (designObj[["testType"]]=="twoSample") {
      if (someRow[["factor"]]==factorLevels[1]) {
        sourceDataTracker[[someSource]][["dat1X"]] <- dat1X <- c(dat1X, someRow[["variable1"]])
        sourceDataTracker[[someSource]][["dat1Y"]] <- dat1Y <- c(dat1Y, someRow[["variable2"]])
      } else if (someRow[["factor"]]==factorLevels[2]) {
        sourceDataTracker[[someSource]][["dat2X"]] <- dat2X <- c(dat2X, someRow[["variable1"]])
        sourceDataTracker[[someSource]][["dat2Y"]] <- dat2Y <- c(dat2Y, someRow[["variable2"]])
      }

      someCheck <- checkTwoSample(dat1X, dat2X, nMin=4)
    } else if (designObj[["testType"]]=="oneSample") {
      sourceDataTracker[[someSource]][["dat1X"]] <- dat1X <- c(dat1X, someRow[["variable1"]])
      sourceDataTracker[[someSource]][["dat1Y"]] <- dat1Y <- c(dat1Y, someRow[["variable2"]])

      dat2X <- dat2Y <- NULL

      someCheck <- checkOneSample(dat1X, nMin=4)
    }

    if (someCheck) {
      logEValueOld <- logETracker[[someSource]]
      logEValueMinEffiOld <- logEMinEffiTracker[[someSource]]

      tempCor <- computeCor(x1=dat1X, y1=dat1Y, x2=dat2X, y2=dat2Y)

      r1Now <- tempCor[["r1"]]
      r2Now <- tempCor[["r2"]]

      n1Now <- length(dat1X)-3
      n2Now <- length(dat2X)-3

      if (designObj[["testType"]]=="oneSample") {
        nEffNow <- n1Now
        n2Now <- NULL

        zScore <-  sqrt(nEffNow)*(atanh(r1Now))
      } else if (designObj[["testType"]]=="twoSample") {
        nEffNow <- n1Now*n2Now/(n1Now+n2Now)

        zScore <-  sqrt(nEffNow)*(atanh(r1Now)-atanh(r2Now))
      }

      if (nEffNow < nEffMin) {
        tempRes <- list(eValue=1)
        tempResMinEffi <- list(eValue=1)
      } else {
        tempRes <- saviZTestStat("z"=zScore, "n1"=n1Now, "n2"=n2Now,
                                 "parameter"=designObj[["parameter"]],
                                 "eType"=designObj[["eType"]],
                                 "alternative"=designObj[["alternative"]])

        tempResMinEffi <- try(saviMinEffiZStat("z"=zScore, "n1"=n1Now, "n2"=n2Now,
                                         "parameter"=designObj[["minEffiTestResult"]][["parameter"]],
                                         "alternative"=designObj[["alternative"]]))
      }

      logEValueNow <- logETracker[[someSource]] <-
        log(tempRes[["eValue"]])

      logEValueMinEffiNow <- logEMinEffiTracker[[someSource]] <-
        log(tempResMinEffi[["eValue"]])

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueMinEffiNow <= log(betaMinEffi))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEMinEffiAdd <- logEValueMinEffiNow - logEValueMinEffiOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEMinEffiNow <- logMetaEMinEffiNow+logMetaEMinEffiAdd

      if (logMetaENow >= log(1/alphaMeta) || logMetaEMinEffiNow <= log(betaMinEffiMeta)) {
        break()
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEMinEffi=logMetaEMinEffiNow,
              logEValuesMinEffi=logEMinEffiTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}

computeCor <- function(x1, y1, x2, y2) {
  r1 <- cor(x1, y1)

  r2 <- NULL

  if (!is.null(x2) && !is.null(y2))
    r2 <- cor(x2, y2)

  if (is.na(r1)) {
    if (length(unique(x1))==1)
      x1 <- perturbX(x1)

    if (length(unique(y1))==1)
      y1 <- perturbX(y1)

    r1 <- cor(x1, y1)
  }

  if (!is.null(r2) && is.na(r2)) {
    if (length(unique(x2))==1)
      x2 <- perturbX(x2)

    if (length(unique(y2))==1)
      y2 <- perturbX(y2)

    r2 <- cor(x2, y2)
  }

  if (abs(r1)==1) {
    x1 <- perturbX(x1)
    y1 <- perturbX(y1)
    r1 <- cor(x1, y1)
  }

  if (!is.null(r2) && abs(r2)==1) {
    x2 <- perturbX(x2)
    y2 <- perturbX(y2)
    r2 <- cor(x2, y2)
  }

  return(list(r1=r1, r2=r2))
}


perturbX <- function(x, sd=0.1) {
  # if (length(unique(x))==1) {
  #   x <- x + rnorm(length(x))
  # }

  x <- x + rnorm(length(x), sd=sd)

  return(x)
}


