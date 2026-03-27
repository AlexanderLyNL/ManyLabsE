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


checkXY <- function(x, y, nMin=1) {
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

# T-test ----
scenario1T <- function(dat, allSources, designObj,
                       nuMin=3, wantCi=FALSE,
                       alphaMeta=0.05, betaFutilityMeta=alphaMeta,
                       seed=NULL, nSim=1e3L) {

  metaScenario1(dat=dat, allSources=allSources, designObj=designObj,
                nuMin=nuMin, wantCi=wantCi,
                alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
                seed=seed, nSim=nSim)
}

metaScenario1 <- function(dat, allSources, designObj,
                          nuMin=3, wantCi=FALSE,
                          alphaMeta=0.05, betaFutilityMeta=alphaMeta,
                          seed=NULL, nSim=1e3L,
                          alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)
  n1Vec <- n2Vec <- integer(nSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$source==allSources[i], ]

    if (designObj[["testName"]]=="T-Test") {
      tempRes <- scenario1TTestHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "wantCi"=wantCi, "nuMin"=nuMin)
    } else if (designObj[["testName"]]=="Z-Test") {
      stop("Z-test not yet done")
    } else if (designObj[["testName"]]=="2x2") {
      stop("2x2 not yet done")
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- scenario1CorHelp("someDat"=someDat, "designObj"=designObj,
                                  "factorLevels"=factorLevels)
    }

    n1Vec[i] <- tempRes[["n1"]]
    n2Vec[i] <- tempRes[["n2"]]
    pValues[i] <- tempRes[["pValue"]]
    eValues[i] <- tempRes[["eValue"]]
    eValuesFut[i] <- tempRes[["eValueFut"]]
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "n1Vec"=n1Vec, "n2Vec"=n2Vec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaFutilityMeta"=betaFutilityMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}

scenario1TTestHelp <- function(someDat, designObj, factorLevels=NULL,
                               wantCi=FALSE, nuMin=3) {

  res <- list(eValue=NULL, eValueFut=NULL, n1=NULL, n2=NULL, pValue=NULL)

  ## Data ---
  if (!is.null(factorLevels)) {
    x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
    y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable

    # Remove non-available entries
    x <- x[!is.na(x)]
    n1 <- length(x)

    y <- y[!is.na(y)]
    n2 <- length(y)
  }

  alternativeOld <- switch(designObj[["alternative"]],
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  tempResult <- t.test(x[1:n1], y[1:n2],
                       alternative=alternativeOld,
                       var.equal=designObj[["varEqual"]])

  res[["pValue"]] <- tempResult[["p.value"]]

  tempRes <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
                       "sequential"=FALSE, "wantCi"=wantCi,
                       "nuMin"=nuMin)
  res[["eValue"]] <- tempRes[["eValue"]]
  res[["eValueFut"]] <- tempRes[["eValueFut"]]

  res[["n1"]] <- n1
  res[["n2"]] <- n2

  return(res)
}

computeWorstCaseScenario1 <- function(
    res, alphaMeta=0.05, betaFutilityMeta=alphaMeta,
    seed=NULL, nSim=1e3L, ...) {

  eValues <- sort(res$eValues)
  eValuesFut <- sort(res$eValuesFut, decreasing=TRUE)

  nSources <- length(eValuesFut)

  nStudiesAlternativeWorstCase <- min(which(cumsum(log(eValues)) >= log(1/alphaMeta)))
  nStudiesFutilityWorstCase <- min(which(cumsum(log(eValuesFut)) <= log(betaFutilityMeta)))

  if (is.infinite(nStudiesAlternativeWorstCase)) {
    nSamplesAlternativeWorstCase <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValues)
    indexStudiesNeeded <- someOrder[1:nStudiesAlternativeWorstCase]

    nSamplesAlternativeWorstCase <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  if (is.infinite(nStudiesFutilityWorstCase)) {
    nSamplesFutilityWorstCase <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValuesFut)
    indexStudiesNeeded <- someOrder[1:nStudiesFutilityWorstCase]

    nSamplesFutilityWorstCase <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  stopDecision <- nStudies <- totalStoppingTimes <- integer(nSim)
  logMetaE <- logMetaEFut <- integer(nSim)

  set.seed(seed)
  for (i in 1:nSim) {
    someOrder <- sample(nSources, nSources)

    tempEValues <- eValues[someOrder]
    tempEValuesFut <- eValuesFut[someOrder]

    logMetaETemp <- cumsum(log(tempEValues))
    logMetaEFutTemp <- cumsum(log(tempEValuesFut))

    tauForAlt <- min(which(logMetaETemp >= log(1/alphaMeta)))
    tauForFutility <- min(which(logMetaEFutTemp <= log(betaFutilityMeta)))

    if (tauForFutility < tauForAlt)
      stopDecision[i] <- -1

    if (tauForAlt < tauForFutility)
      stopDecision[i] <- 1

    tauRace <- min(tauForAlt, tauForFutility)

    stopIndex <- nStudies[i] <- min(tauRace, nSources)

    logMetaE[i] <- logMetaETemp[stopIndex]
    logMetaEFut[i] <- logMetaEFutTemp[stopIndex]

    indexNeededStudies <- someOrder[1:stopIndex]

    totalStoppingTimes[i] <- sum(res$n1Vec[indexNeededStudies])+sum(res$n2Vec[indexNeededStudies])
  }

  res <- list("nStudiesAlternativeWorstCase"=nStudiesAlternativeWorstCase,
              "nStudiesFutilityWorstCase"=nStudiesFutilityWorstCase,
              "nSamplesAlternativeWorstCase"=nSamplesAlternativeWorstCase,
              "nSamplesFutilityWorstCase"=nSamplesFutilityWorstCase,
              "stopDecision"=stopDecision,
              "logMetaE"=logMetaE,
              "logMetaEFut"=logMetaEFut,
              "nStudies"=nStudies,
              "totalStoppingTimes"=totalStoppingTimes)

  return(res)
}

scenario2T <- function(dat, allSources, designObj, alpha=0.05,
                       betaFutility=alpha, nuMin=3, nSim=1e2L,
                       nMax=NULL, seed=NULL, wantCi=FALSE,
                       alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)
  nSources <- length(allSources)

  nSamples <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$source==allSources[i], ]

    ## Data ---
    x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
    y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable

    # Remove non-available entries
    x <- x[!is.na(x)]
    n1 <- length(x)

    y <- y[!is.na(y)]
    n2 <- length(y)

    if (!is.null(designObj$nPlan)) {
      n1 <- min(n1, designObj$nPlan[1])
      n2 <- min(n2, designObj$nPlan[2])
    }

    nParticipants <- n1+n2

    for (k in 1:nSim) {
      tempRes <- twoSampleTTestRandomOrder(
        "x"=x, "y"=y, "n1"=n1, "n2"=n2,
        "designObj"=designObj, "nuMin"=nuMin,
        "wantCi"=wantCi, "nMax"=nMax, "nSim"=nSim,
        "seed"=seed
      )

      nSamples[k, i] <- tempRes$nSamples
      eValues[k, i] <- tempRes$eValue
      eValuesFut[k, i] <- tempRes$eValueFut
    }
  }

  alternativeProportion <- futilityProportion <- numeric(length=nSim)

  for (i in 1:nSim) {
    alternativeProportion[i] <- mean(eValues[i, ] >= 1/alpha)
    futilityProportion[i] <- mean(eValuesFut[i, ] <= betaFutility)
  }

  totalStoppingTimes <- rowSums(nSamples)

  res <- list("nSamples"=nSamples, "eValues"=eValues, "eValuesFut"=eValuesFut,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion,
              "totalStoppingTimes"=totalStoppingTimes)
  return(res)
}

metaScenario2 <- function(dat, allSources, designObj, alphaMeta=0.05,
                          betaFutilityMeta=alphaMeta, nuMin=3, nSim=1e2L,
                          nMax=NULL, seed=NULL, wantCi=FALSE, nEffMin=2) {

  alternative <- designObj[["alternative"]]
  nSources <- length(allSources)

  nSamples <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$source==allSources[i], ]

    if (designObj[["testName"]]=="T-Test") {
      tempRes <- scenario2TTestHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "wantCi"=wantCi,
        "nuMin"=nuMin, "nSim"=nSim, "seed"=seed)
    } else if (designObj[["testName"]]=="Z-Test") {
      stop("Z-test not yet done")
    } else if (designObj[["testName"]]=="2x2") {
      stop("2x2 not yet done")
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- scenario2CorHelp(
        "someDat"=someDat, "designObj"=designObj,
        "factorLevels"=factorLevels, "nSim"=nSim,
        "nEffMin"=nEffMin)
    }

    nSamples[, i] <- tempRes[["nSamples"]]
    eValues[, i] <- tempRes[["eValues"]]
    eValuesFut[, i] <- tempRes[["eValuesFut"]]
  }

  alternativeProportion <- futilityProportion <- numeric(length=nSim)

  for (i in 1:nSim) {
    alternativeProportion[i] <- mean(eValues[i, ] >= 1/designObj[["alpha"]])
    futilityProportion[i] <- mean(eValuesFut[i, ] <= designObj$futilityResult$beta)
  }

  totalStoppingTimes <- rowSums(nSamples)

  res <- list("nSamples"=nSamples, "eValues"=eValues, "eValuesFut"=eValuesFut,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion,
              "totalStoppingTimes"=totalStoppingTimes)
  return(res)
}


scenario2TTestHelp <- function(
    someDat, designObj, factorLevels,
    wantCi=FALSE, nuMin=3, seed=NULL,
    nSim=1e3L, ...) {

  res <- list(nSamples=NULL, eValue=NULL, eValueFut=NULL)

  nSamples <- integer(nSim)
  eValues <- eValuesFut <- numeric(nSim)

  ## Data ---
  x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
  y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable

  # Remove non-available entries
  x <- x[!is.na(x)]
  n1 <- length(x)

  y <- y[!is.na(y)]
  n2 <- length(y)

  if (!is.null(designObj$nPlan)) {
    n1 <- min(n1, designObj$nPlan[1])
    n2 <- min(n2, designObj$nPlan[2])
  }

  nParticipants <- n1+n2

  for (k in 1:nSim) {
    tempRes <- twoSampleTTestRandomOrder(
      "x"=x, "y"=y, "n1"=n1, "n2"=n2,
      "designObj"=designObj, "nuMin"=nuMin,
      "wantCi"=wantCi, "seed"=seed
    )

    nSamples[k] <- tempRes[["nSamples"]]
    eValues[k] <- tempRes[["eValue"]]
    eValuesFut[k] <- tempRes[["eValueFut"]]
  }

  res[["nSamples"]] <- nSamples
  res[["eValues"]] <- eValues
  res[["eValuesFut"]] <- eValuesFut

  return(res)
}

twoSampleTTestRandomOrder <- function(
    x, y, n1, n2, designObj,
    nuMin=3, wantCi=FALSE, nMax=NULL,
    seed=NULL, ...) {

  alpha <- designObj[["alpha"]]
  betaFutility <- designObj[["futilityResult"]][["beta"]]

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

    someCheck <- checkXY(xRun, yRun)

    if (someCheck) {
      tempRes <- saviTTest(
        "x"=xRun, "y"=yRun, "designObj"=designObj,
        "sequential"=FALSE, "nuMin"=nuMin, "wantCi"=wantCi)

      eNow <- tempRes$eValue
      eFutNow <- tempRes$eValueFut

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueFut"=eFutNow)
        return(res)
      }
    }
  }
}

scenario3T <- function(dat, allSources, designObj, alphaMeta=0.05,
                       betaFutilityMeta=alphaMeta, nuMin=3, nSim=1e3L,
                       nMax=NULL, seed=NULL, wantCi=FALSE,
                       nPlanLimit=FALSE) {
  nTotal <- length(unique(dat$uID))
  nSources <- length(allSources)

  logMetaE <- logMetaEFut <- numeric(nSim)

  alternativeProportion <- futilityProportion <- totalStoppingTimes <-
    integer(nSim)

  nSamples <- nStopDecision <- matrix(nrow=nSim, ncol=nSources)
  logEValues <- logEValuesFut <- matrix(nrow=nSim, ncol=nSources)

  for (i in 1:nSim) {
    tempRes <- computeScenario3TOneSim(
      dat=dat, allSources=allSources, designObj=designObj,
      alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta, nuMin=nuMin,
      nSim=nSim, wantCi=wantCi, nPlanLimit=nPlanLimit, seed=seed)

    logMetaE[i] <- tempRes$logMetaE
    logMetaEFut[i] <- tempRes$logMetaEFut
    logEValues[i, ] <- tempRes$logEValues
    logEValuesFut[i, ] <- tempRes$logEValuesFut
    nSamples[i, ] <- tempRes$nSamples
    nStopDecision[i, ] <- tempRes$stopDecision

    totalStoppingTimes[i] <- sum(tempRes$nSamples)
    alternativeProportion[i] <- mean(tempRes$stopDecision==1)
    futilityProportion[i] <- mean(tempRes$stopDecision==-1)
  }

  res <- list("logMetaE"=logMetaE, "logMetaEFut"=logMetaEFut,
              "logEValues"=logEValues, "logEValuesFut"=logEValuesFut,
              "nSamples"=nSamples, "nStopDecision"=nStopDecision,
              "totalStoppingTimes"=totalStoppingTimes,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion)
  return(res)
}

computeScenario3TOneSim <- function(
    dat, allSources, designObj, alphaMeta=0.05,
    betaFutilityMeta=alphaMeta, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

  alpha <- designObj[["alpha"]]
  betaFutility <- designObj[["futilityResult"]][["beta"]]

  nSources <- length(allSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources

  for (neem in allSources)
    sourceDataTracker[[neem]] <- list(x=NULL, y=NULL)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEFutTracker <- logETracker

  nTotal <- length(dat$uID)

  set.seed(seed)
  someOrder <- sample(unique(dat$uID), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEFutNow <- 0

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
    if (nPlanLimit && length(x) >= designObj$nPlan[1] && length(y) >= designObj$nPlan[2])
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()


    if (someRow$factor==factorLevels[1]) {
      sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
    } else if (someRow$factor==factorLevels[2]) {
      sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
    }

    someCheck <- checkXY(x, y)

    if (someCheck) {
      logEValueOld <- logETracker[[someSource]]
      logEValueFutOld <- logEFutTracker[[someSource]]

      tempRes <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
                           "sequential"=FALSE, "wantCi"=wantCi)

      logEValueNow <- logETracker[[someSource]] <-
        log(tempRes$eValue)
      logEValueFutNow <- logEFutTracker[[someSource]] <-
        log(tempRes$eValueFut)

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueFutNow <= log(betaFutility))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEFutAdd <- logEValueFutNow - logEValueFutOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEFutNow <- logMetaEFutNow+logMetaEFutAdd

      if (logMetaENow >= log(1/alphaMeta) || logMetaEFutNow <= log(betaFutilityMeta)) {
        break
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEFut=logMetaEFutNow,
              logEValuesFut=logEFutTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}

metaScenario3 <- function(dat, allSources, designObj, alphaMeta=0.05,
                          betaFutilityMeta=alphaMeta, nuMin=3, nSim=1e3L,
                          nMax=NULL, seed=NULL, wantCi=FALSE,
                          nPlanLimit=FALSE, nEffMin=2) {

  nTotal <- length(unique(dat$uID))
  nSources <- length(allSources)

  logMetaE <- logMetaEFut <- numeric(nSim)

  alternativeProportion <- futilityProportion <- totalStoppingTimes <-
    integer(nSim)

  nSamples <- nStopDecision <- matrix(nrow=nSim, ncol=nSources)
  logEValues <- logEValuesFut <- matrix(nrow=nSim, ncol=nSources)

  # set.seed(seed)
  for (i in 1:nSim) {
    if (designObj[["testName"]]=="T-Test") {
      tempRes <- computeScenario3TOneSim(
        dat=dat, allSources=allSources, designObj=designObj,
        alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
        nuMin=nuMin, nSim=nSim,
        wantCi=wantCi, nPlanLimit=nPlanLimit, seed=seed+i)
    } else if (designObj[["testName"]]=="Z-Test") {
      stop("Z-test not yet done")
    } else if (designObj[["testName"]]=="2x2") {
      stop("2x2 not yet done")
    } else if (designObj[["testName"]]=="Correlation") {
      tempRes <- computeScenario3CorOneSim(
        dat=dat, allSources=allSources, designObj=designObj,
        alphaMeta=alphaMeta, betaFutilityMeta=betaFutilityMeta,
        nuMin=nuMin, nSim=nSim,
        wantCi=wantCi, nPlanLimit=nPlanLimit, nEffMin=nEffMin,
        seed=seed+i)
    }

    logMetaE[i] <- tempRes[["logMetaE"]]
    logMetaEFut[i] <- tempRes[["logMetaEFut"]]
    logEValues[i, ] <- tempRes[["logEValues"]]
    logEValuesFut[i, ] <- tempRes[["logEValuesFut"]]
    nSamples[i, ] <- tempRes[["nSamples"]]
    nStopDecision[i, ] <- tempRes[["stopDecision"]]

    totalStoppingTimes[i] <- sum(tempRes[["nSamples"]])
    alternativeProportion[i] <- mean(tempRes[["stopDecision"]]==1)
    futilityProportion[i] <- mean(tempRes[["stopDecision"]]==-1)
  }

  res <- list("logMetaE"=logMetaE, "logMetaEFut"=logMetaEFut,
              "logEValues"=logEValues, "logEValuesFut"=logEValuesFut,
              "nSamples"=nSamples, "nStopDecision"=nStopDecision,
              "totalStoppingTimes"=totalStoppingTimes,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion)
  return(res)
}


# Tables --------
scenario1Table <- function(dat, allSources, designObj,
                           alpha=0.05, betaFutility=alpha,
                           seed=NULL, nSim=1e3L) {

  alternative <- designObj$alternative

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)

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

    tempRes <- saviFutilityTwoPropConditionalStat(
      ya=ya, na=na, nb=nb, n1=n1,
      logOddsRatio=designObj$futilityResult$parameter,
      alternative=alternative)

    eValuesFut[i] <- tempRes$eValue
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues)

  print("Broken off function")
  return(tempRes)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaFutility"=betaFutilityMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

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

saviFutilityTwoPropConditionalStat <- function(
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
scenario1Z_binomial <- function(dat, allSources, designObj,
                                nuMin=3, wantCi=FALSE,
                                alpha=0.05, betaFutility=alpha,
                                seed=NULL, nSim=1e3L,
                                alternative=c("twoSided", "greater", "less")) {

  #alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)
  nVec <- integer(nSources)

  # factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$source==allSources[i], ]
    count <- as.integer(someDat$variable1=="Parent B")

    datAward <- someDat[someDat$variable2=="Award", ]
    datDeny <- someDat[someDat$variable2=="Deny", ]

    meansAward <- datAward %>%
      summarise(mean=mean(count, na.rm=TRUE))
    meansDeny <- datDeny %>%
      summarise(mean=mean(count, na.rm=TRUE))

    n <- length(count)
    nVec[i] <- n

    allMeans <- (meansAward$mean+meansDeny$mean)/2
    se <- sqrt(0.5*(1-0.5) / length(count))
    z_score <- (allMeans-0.5) / se

    pValues[i] <- 1-pnorm(abs(z_score))

    tempRes <-  saviZTestStat(z=z_score, n1=length(count),
                              parameter=designObj$parameter,
                              alternative=designObj$alternative,
                              sigma=designObj$sigma,
                              eType=designObj$eType)
    eValues[i] <- tempRes$eValue
    tempRes <- saviFutilityZStat(z=z_score, n1=length(count),
                                 parameter=designObj$futilityResult$parameter,
                                 alternative=designObj$alternative,
                                 sigma=designObj$sigma)
    eValuesFut[i] <- tempRes$eValue
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "nVec"=nVec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaFutility"=betaFutilityMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}




scenario2Z_binomial <- function(dat, allSources, designObj, alpha=0.05,
                                betaFutility=alpha, nuMin=3, nSim=1e2L,
                                seed=NULL, wantCi=FALSE,
                                alternative=c("twoSided", "greater", "less")) {

  #alternative <- match.arg(alternative)
  nSources <- length(allSources)

  nSamples <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

  #factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    print(allSources[i])
    someDat <- dat[dat$source==allSources[i], ]

    someDat <- someDat[!is.na(someDat$variable1), ]
    someDat <- someDat[!is.na(someDat$variable2), ]
    n <- length(someDat$variable1)

    if (!is.null(designObj$nPlan)) {
      n <- min(n, designObj$nPlan)
    }

    set.seed(seed)
    for (k in 1:nSim) {
      print(k)
      tempRes <- Z_binomialTestRandomOrder(
        "x"=someDat, "n"=n,
        "designObj"=designObj,# "nuMin"=nuMin,
        "alpha"=alpha, "betaFutility"=betaFutility,
        #"wantCi"=wantCi,
        "nMax"=n
      )

      nSamples[k, i] <- tempRes$nSamples
      eValues[k, i] <- tempRes$eValue
      eValuesFut[k, i] <- tempRes$eValueFut
    }
  }

  alternativeProportion <- futilityProportion <- numeric(length=nSim)

  for (i in 1:nSim) {
    alternativeProportion[i] <- mean(eValues[i, ] >= 1/alpha)
    futilityProportion[i] <- mean(eValuesFut[i, ] <= betaFutility)
  }

  totalStoppingTimes <- rowSums(nSamples)

  res <- list("nSamples"=nSamples, "eValues"=eValues, "eValuesFut"=eValuesFut,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion,
              "totalStoppingTimes"=totalStoppingTimes)
  return(res)
}

Z_binomialTestRandomOrder <- function(
    x, n, designObj, # nuMin=3,
    alpha=0.05,
    betaFutility=alpha, # wantCi=FALSE
    nMax=NULL) {

  #xTemp <- sample(x, length(x))
  xTemp <- x[sample(nrow(x)), ]

  for (j in 1:n) {
    xRun <- xTemp[1:j, ]
    count <- as.integer(xRun$variable1=="Parent B")

    datAward <- xRun[xRun$variable2=="Award", ]
    datDeny <- xRun[xRun$variable2=="Deny", ]
    if (length(datAward$uID) > 2 && length(datDeny$uID) > 2) {

      meansAward <- datAward %>%
        summarise(mean=mean(count, na.rm=TRUE))
      meansDeny <- datDeny %>%
        summarise(mean=mean(count, na.rm=TRUE))

      allMeans <- (meansAward$mean+meansDeny$mean)/2
      se <- sqrt(0.5*(1-0.5) / length(count))
      z_score <- (allMeans-0.5) / se

      #pValues[i] <- 1-pnorm(abs(z_score))
      tempRes <-  saviZTestStat(z=z_score, n1=length(count),
                                parameter=designObj$parameter,
                                alternative=designObj$alternative,
                                sigma=designObj$sigma,
                                eType=designObj$eType)
      eNow <- tempRes$eValue
      tempRes <- saviFutilityZStat(z=z_score, n1=length(count),
                                   parameter=designObj$futilityResult$parameter,
                                   alternative=designObj$alternative,
                                   sigma=designObj$sigma)
      eFutNow <- tempRes$eValue

      #eNow <- tempRes$eValue
      #eFutNow <- tempRes$eValueFut

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueFut"=eFutNow)
        return(res)
      }
    }
  }
}





scenario3ZBinom <- function(dat, allSources, designObj, alpha=0.05,
                            betaFutility=alpha, nuMin=3, nSim=1e3L,
                            nMax=NULL, seed=NULL, wantCi=FALSE,
                            nPlanLimit=FALSE) {
  nTotal <- length(unique(dat$uID))
  nSources <- length(allSources)

  logMetaE <- logMetaEFut <- numeric(nSim)

  alternativeProportion <- futilityProportion <- totalStoppingTimes <-
    integer(nSim)

  nSamples <- nStopDecision <- matrix(nrow=nSim, ncol=nSources)
  logEValues <- logEValuesFut <- matrix(nrow=nSim, ncol=nSources)

  set.seed(seed)
  for (i in 1:nSim) {
    print(i)

    tempRes <- computeScenario3ZbinomOneSim(
      dat=dat, allSources=allSources, designObj=designObj,
      alpha=alpha, betaFutility=betaFutility, nuMin=nuMin, nSim=nSim,
      wantCi=wantCi, nPlanLimit=nPlanLimit)

    logMetaE[i] <- tempRes$logMetaE
    logMetaEFut[i] <- tempRes$logMetaEFut
    logEValues[i, ] <- tempRes$logEValues
    logEValuesFut[i, ] <- tempRes$logEValuesFut
    nSamples[i, ] <- tempRes$nSamples
    nStopDecision[i, ] <- tempRes$stopDecision

    totalStoppingTimes[i] <- sum(tempRes$nSamples)
    alternativeProportion[i] <- mean(tempRes$stopDecision==1)
    futilityProportion[i] <- mean(tempRes$stopDecision==-1)
  }

  res <- list("logMetaE"=logMetaE, "logMetaEFut"=logMetaEFut,
              "logEValues"=logEValues, "logEValuesFut"=logEValuesFut,
              "nSamples"=nSamples, "nStopDecision"=nStopDecision,
              "totalStoppingTimes"=totalStoppingTimes,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion)
  return(res)
}

computeScenario3ZbinomOneSim <- function(
    dat, allSources, designObj, alpha=0.05,
    betaFutility=alpha, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

  nSources <- length(allSources)

  #factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources
  cols <- colnames(dat)
  for (neem in allSources)
    sourceDataTracker[[neem]] <- setNames(vector("list", length(cols)), cols)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEFutTracker <- logETracker

  nTotal <- length(dat$uID)

  someOrder <- sample(unique(dat$uID), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEFutNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat$uID==someId), ]
    someSource <- someRow$source

    nSamples[[someSource]] <- nSamples[[someSource]] + 1

    sourceDataTemp <- sourceDataTracker[[someSource]]

    # Retrieve old values from state
    #
    oldData <- sourceDataTemp

    # Skip if sample size limit is reached within trial
    #
    if (nPlanLimit && length(sourceDataTemp$uID) >= designObj$nPlan)
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()
    someRowList <- as.list(someRow[1, ])
    sourceDataTracker[[someSource]] <- Map(
      c,
      sourceDataTracker[[someSource]],
      someRowList
    )
    #if (someRow$factor==factorLevels[1]) {
    #  sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
    #} else if (someRow$factor==factorLevels[2]) {
    #  sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
    #}

    #someCheck <- checkXY(x, y)
    xRun <- sourceDataTracker[[someSource]]
    xRun_df <- as.data.frame(xRun)
    count <- as.integer(xRun_df$variable1==2)

    datAward <- xRun_df[xRun_df$variable2==1, ]
    datDeny <- xRun_df[xRun_df$variable2==2, ]
    if (length(datAward$uID) > 2 && length(datDeny$uID) > 2) {
      logEValueOld <- logETracker[[someSource]]
      logEValueFutOld <- logEFutTracker[[someSource]]

      meansAward <- datAward %>%
        summarise(mean=mean(count, na.rm=TRUE))
      meansDeny <- datDeny %>%
        summarise(mean=mean(count, na.rm=TRUE))

      allMeans <- (meansAward$mean+meansDeny$mean)/2
      se <- sqrt(0.5*(1-0.5) / length(count))
      z_score <- (allMeans-0.5) / se

      #pValues[i] <- 1-pnorm(abs(z_score))
      tempRes <-  saviZTestStat(z=z_score, n1=length(count),
                                parameter=designObj$parameter,
                                alternative=designObj$alternative,
                                sigma=designObj$sigma,
                                eType=designObj$eType)
      eNow <- tempRes$eValue
      tempRes <- saviFutilityZStat(z=z_score, n1=length(count),
                                   parameter=designObj$futilityResult$parameter,
                                   alternative=designObj$alternative,
                                   sigma=designObj$sigma)
      eFutNow <- tempRes$eValue


      logEValueNow <- logETracker[[someSource]] <-
        log(eNow)
      logEValueFutNow <- logEFutTracker[[someSource]] <-
        log(eFutNow)

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueFutNow <= log(betaFutility))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEFutAdd <- logEValueFutNow - logEValueFutOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEFutNow <- logMetaEFutNow+logMetaEFutAdd

      if (logMetaENow >= log(1/alpha) || logMetaEFutNow <= log(betaFutility)) {
        break
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEFut=logMetaEFutNow,
              logEValuesFut=logEFutTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}



# Corrie Z --------
scenario1ZCorr <- function(dat, allSources, designObj,
                           nuMin=3, wantCi=FALSE,
                           alpha=0.05, betaFutility=alpha,
                           seed=NULL, nSim=1e3L,
                           alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)
  nVec <- integer(nSources)

  # Alexander(TODO): check one or two-sample
  # factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    someDat <- dat[dat$source==allSources[i], ]
    someN <- length(someDat$uID)
    nVec[i] <- someN
    someR  <- cor(someDat$variable1, someDat$variable2)
    someZ  <- atanh(someR)/sqrt(1/(length(someDat$uID)-3))
    pValues[i] <- 1-pnorm(abs(someZ))

    tempRes <- saviZTestStat(
      z=someZ, n1=someN,
      parameter=designObj$parameter,
      eType=designObj$eType, sigma=designObj$sigma)

    eValues[i] <- tempRes$eValue

    tempRes <- saviFutilityZStat(
      z=someZ, n1=someN,
      parameter=designObj$futilityResult$parameter,
      sigma=designObj$sigma)

    eValuesFut[i] <- tempRes$eValue

  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "nVec"=nVec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alphaMeta"=alphaMeta, "betaFutilityMeta"=betaFutilityMeta,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}



scenario2ZCorr  <- function(dat, allSources, designObj, alpha=0.05,
                            betaFutility=alpha, nuMin=3, nSim=1e2L,
                            seed=NULL, wantCi=FALSE,
                            alternative=c("twoSided", "greater", "less")) {

  #alternative <- match.arg(alternative)
  nSources <- length(allSources)

  nSamples <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

  #factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  for (i in 1:length(allSources)) {
    # print(allSources[i])
    someDat <- dat[dat$source==allSources[i], ]

    someDat <- someDat[!is.na(someDat$variable1), ]
    someDat <- someDat[!is.na(someDat$variable2), ]
    n <- length(someDat$variable1)

    if (!is.null(designObj$nPlan))
      n <- min(n, designObj$nPlan)

    set.seed(seed)
    for (k in 1:nSim) {
      # print(k)
      tempRes <- ZCorrTestRandomOrder(
        "x"=someDat, "n"=n,
        "designObj"=designObj,# "nuMin"=nuMin,
        "alpha"=alpha, "betaFutility"=betaFutility,
        #"wantCi"=wantCi,
        "nMax"=n
      )

      nSamples[k, i] <- tempRes$nSamples
      eValues[k, i] <- tempRes$eValue
      eValuesFut[k, i] <- tempRes$eValueFut
    }
  }

  alternativeProportion <- futilityProportion <- numeric(length=nSim)

  for (i in 1:nSim) {
    alternativeProportion[i] <- mean(eValues[i, ] >= 1/alpha)
    futilityProportion[i] <- mean(eValuesFut[i, ] <= betaFutility)
  }

  totalStoppingTimes <- rowSums(nSamples)

  res <- list("nSamples"=nSamples, "eValues"=eValues, "eValuesFut"=eValuesFut,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion,
              "totalStoppingTimes"=totalStoppingTimes)
  return(res)
}

ZCorrTestRandomOrder <- function(
    x, n, designObj, # nuMin=3,
    alpha=0.05,
    betaFutility=alpha, # wantCi=FALSE
    nMax=NULL) {

  #xTemp <- sample(x, length(x))
  xTemp <- x[sample(nrow(x)), ]

  for (j in 1:n) {
    xRun <- xTemp[1:j, ]

    sd1 <- sqrt(var(xRun$variable1))
    sd2 <- sqrt(var(xRun$variable2))
    someN <- nrow(xRun)

    if (someN > 3 && sd1 > 0 && sd2 > 0){
      someR  <- cor(xRun$variable1, xRun$variable2)
      someZ  <- atanh(someR)/sqrt(1/(someN-3))

      tempRes <- saviZTestStat(
        z=someZ, n1=someN,
        parameter=designObj$parameter,
        eType=designObj$eType, sigma=designObj$sigma)

      eNow <- tempRes$eValue

      tempRes <- saviFutilityZStat(
        z=someZ, n1=someN,
        parameter=designObj$futilityResult$parameter,
        sigma=designObj$sigma)

      eFutNow <- tempRes$eValue

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueFut"=eFutNow)
        return(res)
      }
    }
  }
}


scenario3ZCorr <- function(dat, allSources, designObj, alpha=0.05,
                           betaFutility=alpha, nuMin=3, nSim=1e3L,
                           nMax=NULL, seed=NULL, wantCi=FALSE,
                           nPlanLimit=FALSE) {
  nTotal <- length(unique(dat$uID))
  nSources <- length(allSources)

  logMetaE <- logMetaEFut <- numeric(nSim)

  alternativeProportion <- futilityProportion <- totalStoppingTimes <-
    integer(nSim)

  nSamples <- nStopDecision <- matrix(nrow=nSim, ncol=nSources)
  logEValues <- logEValuesFut <- matrix(nrow=nSim, ncol=nSources)

  set.seed(seed)
  for (i in 1:nSim) {
    #print(i)

    tempRes <- computeScenario3ZCorrOneSim(
      dat=dat, allSources=allSources, designObj=designObj,
      alpha=alpha, betaFutility=betaFutility, nuMin=nuMin, nSim=nSim,
      wantCi=wantCi, nPlanLimit=nPlanLimit)

    logMetaE[i] <- tempRes$logMetaE
    logMetaEFut[i] <- tempRes$logMetaEFut
    logEValues[i, ] <- tempRes$logEValues
    logEValuesFut[i, ] <- tempRes$logEValuesFut
    nSamples[i, ] <- tempRes$nSamples
    nStopDecision[i, ] <- tempRes$stopDecision

    totalStoppingTimes[i] <- sum(tempRes$nSamples)
    alternativeProportion[i] <- mean(tempRes$stopDecision==1)
    futilityProportion[i] <- mean(tempRes$stopDecision==-1)
  }

  res <- list("logMetaE"=logMetaE, "logMetaEFut"=logMetaEFut,
              "logEValues"=logEValues, "logEValuesFut"=logEValuesFut,
              "nSamples"=nSamples, "nStopDecision"=nStopDecision,
              "totalStoppingTimes"=totalStoppingTimes,
              "alternativeProportion"=alternativeProportion,
              "futilityProportion"=futilityProportion)
  return(res)
}

computeScenario3ZCorrOneSim <- function(
    dat, allSources, designObj, alpha=0.05,
    betaFutility=alpha, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

  nSources <- length(allSources)

  #factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources
  cols <- colnames(dat)
  for (neem in allSources)
    sourceDataTracker[[neem]] <- setNames(vector("list", length(cols)), cols)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEFutTracker <- logETracker

  nTotal <- length(dat$uID)

  someOrder <- sample(unique(dat$uID), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEFutNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat$uID==someId), ]
    someSource <- someRow$source

    nSamples[[someSource]] <- nSamples[[someSource]] + 1

    sourceDataTemp <- sourceDataTracker[[someSource]]

    # Retrieve old values from state
    #
    oldData <- sourceDataTemp

    # Skip if sample size limit is reached within trial
    #
    if (nPlanLimit && length(sourceDataTemp$uID) >= designObj$nPlan)
      next()

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()
    someRowList <- as.list(someRow[1, ])
    sourceDataTracker[[someSource]] <- Map(
      c,
      sourceDataTracker[[someSource]],
      someRowList
    )
    #if (someRow$factor==factorLevels[1]) {
    #  sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
    #} else if (someRow$factor==factorLevels[2]) {
    #  sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
    #}

    #someCheck <- checkXY(x, y)
    xRun <- sourceDataTracker[[someSource]]
    xRun_df <- as.data.frame(xRun)

    sd1 <- sqrt(var(xRun_df$variable1))
    sd2 <- sqrt(var(xRun_df$variable2))
    someN <- nrow(xRun_df)


    if (someN > 3 && sd1 > 0 && sd2 > 0){
      logEValueOld <- logETracker[[someSource]]
      logEValueFutOld <- logEFutTracker[[someSource]]

      someR  <- cor(xRun$variable1, xRun$variable2)
      someZ  <- atanh(someR)/sqrt(1/(someN-3))

      tempRes <- saviZTestStat(
        z=someZ, n1=someN,
        parameter=designObj$parameter,
        eType=designObj$eType, sigma=designObj$sigma)

      eNow <- tempRes$eValue

      tempRes <- saviFutilityZStat(
        z=someZ, n1=someN,
        parameter=designObj$futilityResult$parameter,
        sigma=designObj$sigma)

      eFutNow <- tempRes$eValue


      logEValueNow <- logETracker[[someSource]] <-
        log(eNow)
      logEValueFutNow <- logEFutTracker[[someSource]] <-
        log(eFutNow)

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueFutNow <= log(betaFutility))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEFutAdd <- logEValueFutNow - logEValueFutOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEFutNow <- logMetaEFutNow+logMetaEFutAdd

      if (logMetaENow >= log(1/alpha) || logMetaEFutNow <= log(betaFutility)) {
        break
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEFut=logMetaEFutNow,
              logEValuesFut=logEFutTracker,
              logEValues=logETracker,
              nSamples=nSamples,
              stopDecision=stopDecision)
  return(res)
}

# 2-Corrie Z --------
scenario1CorHelp <- function(someDat, designObj, factorLevels=NULL) {

  res <- list(eValue=NULL, eValueFut=NULL, n1=NULL, n2=NULL, pValue=NULL)

  ## Data ---
  if (!is.null(factorLevels)) {
    dat1 <- someDat[which(someDat$factor==factorLevels[1]), ]
    dat2 <- someDat[which(someDat$factor==factorLevels[2]), ]

    n1 <- dim(dat1)[1]-3
    n2 <- dim(dat2)[1]-3
  }

  alternativeOld <- switch(designObj[["alternative"]],
                           "twoSided"="two.sided",
                           "greater"="greater",
                           "less"="less")

  corrie1 <- cor(dat1[["variable1"]], dat1[["variable2"]])
  corrie2 <- cor(dat2[["variable1"]], dat2[["variable2"]])

  #TODO(Alexander): Perhaps later

  res[["pValue"]] <- NA

  nEff <- n1*n2/(n1+n2)

  zScore <- sqrt(nEff)*(atanh(corrie1)- atanh(corrie2))

  tempRes <- saviZTestStat(z=zScore, n1=n1, n2=n2, parameter=designObj[["parameter"]],
                           eType=designObj[["eType"]], alternative=designObj[["alternative"]])

  res[["eValue"]] <- tempRes[["eValue"]]

  tempRes <- saviFutilityZStat(z=zScore, n1=n1, n2=n2,
                               parameter=designObj[["futilityResult"]][["parameter"]],
                               alternative=designObj[["alternative"]])

  res[["eValueFut"]] <- tempRes[["eValue"]]
  res[["n1"]] <- n1
  res[["n2"]] <- n2

  return(res)
}


scenario2CorHelp <- function(someDat, designObj, factorLevels,
                             nSim=1e3L, seed=NULL, nEffMin=2) {

  res <- list(nSamples=NULL, eValue=NULL, eValueFut=NULL)

  nSamples <- integer(nSim)
  eValues <- eValuesFut <- numeric(nSim)

  ## Data ---
  if (!is.null(factorLevels)) {
    dat1 <- someDat[which(someDat$factor==factorLevels[1]), ]
    dat2 <- someDat[which(someDat$factor==factorLevels[2]), ]

    n1 <- dim(dat1)[1]
    n2 <- dim(dat2)[1]
  }

  if (!is.null(designObj$nPlan)) {
    n1 <- min(n1, designObj$nPlan[1])
    n2 <- min(n2, designObj$nPlan[2])
  }

  nParticipants <- n1+n2

  for (k in 1:nSim) {

    tempRes <- twoSampleCorTestRandomOrder(dat1, dat2, n1, n2, designObj,
                                           seed=seed, nEffMin=nEffMin)

    nSamples[k] <- tempRes[["nSamples"]]
    eValues[k] <- tempRes[["eValue"]]
    eValuesFut[k] <- tempRes[["eValueFut"]]
  }

  res[["nSamples"]] <- nSamples
  res[["eValues"]] <- eValues
  res[["eValuesFut"]] <- eValuesFut

  return(res)
}

twoSampleCorTestRandomOrder <- function(
    dat1, dat2, n1, n2, designObj,
    seed=NULL, nMax=NULL,
    nEffMin=2) {

  res <- list(nSamples=NULL, eValue=NULL, eValueFut=NULL)

  alpha <- designObj[["alpha"]]
  betaFutility <- designObj[["futilityResult"]][["beta"]]

  nParticipants <- n1+n2

  dat1XRun <- numeric(0)
  dat1YRun <- numeric(0)
  dat2XRun <- numeric(0)
  dat2YRun <- numeric(0)

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

    someCheck <- checkXY(dat1XRun, dat2XRun, nMin=4)

    if (someCheck) {
      tempCor <- computeCor(x1=dat1XRun, y1=dat1YRun, x2=dat2XRun, y2=dat2YRun)

      r1 <- tempCor[["r1"]]
      r2 <- tempCor[["r2"]]

      n1Now <- length(dat1XRun)-3
      n2Now <- length(dat2XRun)-3

      nEffNow <- n1Now*n2Now/(n1Now+n2Now)

      zScore <-  sqrt(nEffNow)*(atanh(r1)-atanh(r2))

      if (is.na(zScore) && nEffNow < nEffMin) {
        eNow <- 1
        eFutNow <- 1
      } else {
        tempRes <- saviZTestStat(z=zScore, n1=n1Now, n2=n2Now,
                                 parameter=designObj[["parameter"]],
                                 alternative=designObj[["alternative"]],
                                 eType=designObj[["eType"]])

        eNow <- tempRes[["eValue"]]

        if (is.infinite(eNow))
          browser()

        tempRes <- try(saviFutilityZStat(z=zScore, n1=n1Now, n2=n2Now,
                                     parameter=designObj[["futilityResult"]][["parameter"]],
                                     alternative=designObj[["alternative"]]))

        if (isTryError(tempRes))
          browser()

        eFutNow <- tempRes[["eValue"]]
      }

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list("nSamples"=j, "eValue"=eNow, "eValueFut"=eFutNow)
        return(res)
      }
    }
  }
}


computeScenario3CorOneSim <- function(
    dat, allSources, designObj, alphaMeta=0.05,
    betaFutilityMeta=alphaMeta, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE, nEffMin=2,
    nPlanLimit=TRUE) {

  alpha <- designObj[["alpha"]]
  betaFutility <- designObj[["futilityResult"]][["beta"]]

  nSources <- length(allSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  sourceDataTracker <- vector(mode="list", length=nSources)
  names(sourceDataTracker) <- allSources

  for (neem in allSources)
    sourceDataTracker[[neem]] <- list(dat1X=NULL, dat1Y=NULL, dat2X=NULL, dat2Y=NULL)

  nSamples <- integer(length=nSources)
  names(nSamples) <- allSources
  stopDecision <- nSamples

  logETracker <- numeric(length=nSources)
  names(logETracker) <- allSources
  logEFutTracker <- logETracker

  nTotal <- length(dat$uID)

  set.seed(seed)
  someOrder <- sample(unique(dat$uID), nTotal)

  # meta eValues are all 1 at the start
  #
  logMetaENow <- logMetaEFutNow <- 0

  for (j in seq_along(someOrder)) {
    someId <- someOrder[j]

    someRow <- dat[which(dat$uID==someId), ]
    someSource <- someRow$source

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
    if (nPlanLimit && length(dat1X) >= designObj$nPlan[1] &&
        length(dat2X) >= designObj$nPlan[2]) {
      next()
    }

    # Skip if already stopped within trial
    #
    if (stopDecision[[someSource]]!=0)
      next()

    if (someRow$factor==factorLevels[1]) {
      sourceDataTracker[[someSource]]$dat1X <- dat1X <- c(dat1X, someRow$variable1)
      sourceDataTracker[[someSource]]$dat1Y <- dat1Y <- c(dat1Y, someRow$variable2)
    } else if (someRow$factor==factorLevels[2]) {
      sourceDataTracker[[someSource]]$dat2X <- dat2X <- c(dat2X, someRow$variable1)
      sourceDataTracker[[someSource]]$dat2Y <- dat2Y <- c(dat2Y, someRow$variable2)
    }

    someCheck <- checkXY(dat1X, dat2X, nMin=4)

    if (someCheck) {
      logEValueOld <- logETracker[[someSource]]
      logEValueFutOld <- logEFutTracker[[someSource]]

      tempCor <- computeCor(x1=dat1X, y1=dat1Y, x2=dat2X, y2=dat2Y)

      r1Now <- tempCor[["r1"]]
      r2Now <- tempCor[["r2"]]

      n1Now <- length(dat1X)-3
      n2Now <- length(dat2X)-3

      nEffNow <- n1Now*n2Now/(n1Now+n2Now)

      if (nEffNow < nEffMin) {
        tempRes <- list(eValue=1)
        tempResFut <- list(eValue=1)
      } else {
        zScore <- sqrt(nEffNow)*(atanh(r1Now)-atanh(r2Now))

        tempRes <- saviZTestStat("z"=zScore, "n1"=n1Now, "n2"=n2Now,
                                 "parameter"=designObj[["parameter"]],
                                 "eType"=designObj[["eType"]],
                                 "alternative"=designObj[["alternative"]])

        tempResFut <- try(saviFutilityZStat("z"=zScore, "n1"=n1Now, "n2"=n2Now,
                                         "parameter"=designObj[["futilityResult"]][["parameter"]],
                                         "alternative"=designObj[["alternative"]]))
      }

      logEValueNow <- logETracker[[someSource]] <-
        log(tempRes[["eValue"]])

      logEValueFutNow <- logEFutTracker[[someSource]] <-
        log(tempResFut[["eValue"]])

      if (logEValueNow >= log(1/alpha))
        stopDecision[[someSource]] <- 1

      if (logEValueFutNow <= log(betaFutility))
        stopDecision[[someSource]] <- -1

      logMetaEAdd <- logEValueNow - logEValueOld
      logMetaEFutAdd <- logEValueFutNow - logEValueFutOld

      logMetaENow <- logMetaENow+logMetaEAdd
      logMetaEFutNow <- logMetaEFutNow+logMetaEFutAdd

      if (logMetaENow >= log(1/alphaMeta) || logMetaEFutNow <= log(betaFutilityMeta)) {
        break()
      }
    }
  }

  res <- list(logMetaE=logMetaENow, logMetaEFut=logMetaEFutNow,
              logEValuesFut=logEFutTracker,
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
