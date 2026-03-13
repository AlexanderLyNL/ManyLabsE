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


checkXY <- function(x, y) {
  if (length(x) >= 1 && length(y) >= 1)
    return(TRUE)

  if (is.null(x))
    return(FALSE)

  if (is.null(y))
    return(FALSE)

  if (length(x)==0)
    return(FALSE)

  if (length(y)==0)
    return(FALSE)

  if (is.na(x))
    return(FALSE)

  if (is.na(y))
    return(FALSE)

  return(TRUE)
}


scenario1T <- function(dat, allSources, designObj,
                       nuMin=3, wantCi=FALSE,
                       alpha=0.05, betaFutility=alpha,
                       seed=NULL, nSim=1e3L,
                       alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)
  n1Vec <- n2Vec <- integer(nSources)

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

    n1Vec[i] <- n1
    n2Vec[i] <- n2

    alternativeOld <- switch(alternative,
                             "twoSided"="two.sided",
                             "greater"="greater",
                             "less"="less")

    tempResult <- t.test(x[1:n1], y[1:n2],
                         alternative=alternativeOld,
                         var.equal=varEqual)

    pValues[i] <- tempResult$p.value

    tempRes <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
                         "sequential"=FALSE, "wantCi"=wantCi,
                         "nuMin"=nuMin)

    eValues[i] <- tempRes$eValue
    eValuesFut[i] <- tempRes$eValueFut
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "n1Vec"=n1Vec, "n2Vec"=n2Vec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alpha"=alpha, "betaFutility"=betaFutility,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}

computeWorstCaseScenario1 <- function(
    res, alpha=0.05, betaFutility=alpha,
    seed=NULL, nSim=1e3L) {

  eValues <- sort(res$eValues)
  eValuesFut <- sort(res$eValuesFut, decreasing=TRUE)

  nSources <- length(eValuesFut)

  nStudiesAlternativeWorstCase <- min(which(cumsum(log(eValues)) >= log(1/alpha)))
  nStudiesFutilityWorstCase <- min(which(cumsum(log(eValuesFut)) <= log(betaFutility)))

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

    tauForAlt <- min(which(logMetaETemp >= log(1/alpha)))
    tauForFutility <- min(which(logMetaEFutTemp <= log(betaFutility)))

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

    set.seed(seed)
    for (k in 1:nSim) {
      tempRes <- twoSampleTTestRandomOrder(
        "x"=x, "y"=y, "n1"=n1, "n2"=n2,
        "designObj"=designObj, "nuMin"=nuMin,
        "alpha"=alpha, "betaFutility"=betaFutility,
        "wantCi"=wantCi, "nMax"=nMax
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



twoSampleTTestRandomOrder <- function(
    x, y, n1, n2, designObj,
    nuMin=3, alpha=0.05,
    betaFutility=alpha, wantCi=FALSE, nMax=NULL) {

  nParticipants <- n1+n2

  xRun <- numeric(0)
  yRun <- numeric(0)

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

scenario3T <- function(dat, allSources, designObj, alpha=0.05,
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

    tempRes <- computeScenario3TOneSim(
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

computeScenario3TOneSim <- function(
    dat, allSources, designObj, alpha=0.05,
    betaFutility=alpha, nuMin=3, nSim=1e3L,
    seed=NULL, wantCi=FALSE,
    nPlanLimit=TRUE) {

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


scenario1Cor <- function(dat, allSources, designObj,
                         nuMin=3, wantCi=FALSE,
                         alpha=0.05, betaFutility=alpha,
                         seed=NULL, nSim=1e3L,
                         alternative=c("twoSided", "greater", "less")) {

  alternative <- match.arg(alternative)

  nSources <- length(allSources)

  eValues <- eValuesFut <- pValues <- numeric(nSources)
  n1Vec <- n2Vec <- integer(nSources)

  factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)

  browser()

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

    n1Vec[i] <- n1
    n2Vec[i] <- n2

    alternativeOld <- switch(alternative,
                             "twoSided"="two.sided",
                             "greater"="greater",
                             "less"="less")

    tempResult <- t.test(x[1:n1], y[1:n2],
                         alternative=alternativeOld,
                         var.equal=varEqual)

    pValues[i] <- tempResult$p.value

    tempRes <- saviTTest("x"=x, "y"=y, "designObj"=designObj,
                         "sequential"=FALSE, "wantCi"=wantCi,
                         "nuMin"=nuMin)

    eValues[i] <- tempRes$eValue
    eValuesFut[i] <- tempRes$eValueFut
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "n1Vec"=n1Vec, "n2Vec"=n2Vec)

  tempRes2 <- computeWorstCaseScenario1(
    tempRes, "alpha"=alpha, "betaFutility"=betaFutility,
    "seed"=seed, "nSim"=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}
