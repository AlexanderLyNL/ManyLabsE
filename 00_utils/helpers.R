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


twoSampleTTestRandomOrder <- function(
    x, y, n1, n2, designObj,
    nuMin=3, alpha=0.05,
    betaFutility=alpha, nMax=n1+n2, wantCi=FALSE) {

  nParticipants <- n1+n2
  someOrder <- sample(nParticipants, nParticipants)

  nMax <- min(nMax, nParticipants)

  xRun <- numeric(0)
  yRun <- numeric(0)

  for (j in seq_along(someOrder)) {
    partId <- someOrder[j]

    xTemp <- sample(x, length(x))
    yTemp <- sample(y, length(y))

    xTemp <- xTemp[1:n1]
    yTemp <- yTemp[1:n2]

    totalVar <- c(xTemp, yTemp)

    if (partId <= n1) {
      xRun <- c(xRun, totalVar[someOrder[j]])
    } else {
      yRun <- c(yRun, totalVar[someOrder[j]])
    }


    someCheck <- checkXY(xRun, yRun)

    if (someCheck) {
      tempRes <- saviTTest(
        xRun, yRun, designObj=designObj,
        sequential=FALSE, nuMin=nuMin, wantCi=wantCi)

      eNow <- tempRes$eValue
      eFutNow <- tempRes$eValueFut

      # length(xRun) >= nMax[1] && length(yRun) >= nMax[2]

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list(nStop=j, eValue=eNow, eValueFut=eFutNow)
        return(res)
      }
    }
  }
}


scenario1T <- function(dat, allSources, designObj,
                       nuMin=3, wantCi=FALSE,
                       alpha=0.05, betaFutility=alpha,
                       seed=NULL, nSim=1e3L) {

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

    tempResult <- t.test(x[1:n1], y[1:n2],
                         var.equal=stat.params$var.equal)
    pValues[i] <- tempResult$p.value

    tempRes <- saviTTest(x=x, y=y, designObj=designObj,
                         sequential=FALSE, wantCi=wantCi,
                         nuMin=nuMin)

    eValues[i] <- tempRes$eValue
    eValuesFut[i] <- tempRes$eValueFut
  }

  tempRes <- list("eValues"=eValues, "eValuesFut"=eValuesFut,
                  "pValues"=pValues,
                  "n1Vec"=n1Vec, "n2Vec"=n2Vec)
  tempRes2 <- computeWorstCaseScenario1(
    tempRes, alpha=alpha, betaFutility=betaFutility,
    seed=seed, nSim=nSim)

  res <- utils::modifyList(tempRes, tempRes2)

  return(res)
}

computeWorstCaseScenario1 <- function(
    res, alpha=0.05, betaFutility=alpha,
    seed=NULL, nSim=1e3L) {

  eValues <- sort(res$eValues)
  eValuesFut <- sort(res$eValuesFut, decreasing=TRUE)

  nSources <- length(eValuesFut)

  mostNStudiesForAlternative <- min(which(cumsum(log(eValues)) >= 1/alpha))
  mostNStudiesForFutility <- min(which(cumsum(log(eValuesFut)) <= betaFutility))

  if (is.infinite(mostNStudiesForAlternative)) {
    mostNSamplesForAlternative <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValues)
    indexStudiesNeeded <- someOrder[1:mostNStudiesForAlternative]

    mostNSamplesForAlternative <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  if (is.infinite(mostNStudiesForFutility)) {
    mostNStudiesForFutility <- sum(res$n1Vec)+sum(res$n2Vec)
  } else {
    someOrder <- order(res$eValuesFut)
    indexStudiesNeeded <- someOrder[1:mostNStudiesForFutility]

    mostNSamplesForFutility <- sum(res$n1Vec[indexStudiesNeeded])+sum(res$n2Vec[indexStudiesNeeded])
  }

  #
  stopForAlt <- stopForFutility <- nStopStudies <- nStopVec <- integer(nSim)
  logMetaEvalues <- logMetaEvaluesFut <- integer(nSim)



  set.seed(seed)
  for (i in 1:nSim) {
    someOrder <- sample(nSources, nSources)

    tempEValues <- eValues[someOrder]
    tempEValuesFut <- eValuesFut[someOrder]

    logMetaE <- cumsum(log(tempEValues))
    logMetaEFut <- cumsum(log(tempEValuesFut))

    tauForAlt <- min(which(logMetaE >= log(1/alpha)))
    tauForFutility <- min(which(logMetaEFut <= log(betaFutility)))

    if (tauForFutility < tauForAlt)
      stopForFutility[i] <- 1

    if (tauForAlt < tauForFutility)
      stopForAlt[i] <- 1

    tauRace <- min(tauForAlt, tauForFutility)

    stopIndex <- nStopStudies[i] <- min(tauRace, nSources)

    logMetaEvalues[i] <- logMetaE[stopIndex]
    logMetaEvaluesFut[i] <- logMetaEFut[stopIndex]

    indexNeededStudies <- someOrder[1:stopIndex]

    nStopVec[i] <- sum(res$n1Vec[indexNeededStudies])+sum(res$n2Vec[indexNeededStudies])
  }

  res <- list("mostNStudiesForAlternative"=mostNStudiesForAlternative,
              "mostNStudiesForFutility"=mostNStudiesForFutility,
              "mostNSamplesForAlternative"=mostNSamplesForAlternative,
              "mostNSamplesForFutility"=mostNSamplesForFutility,
              "stopForAlt"=stopForAlt, "stopForFutility"=stopForFutility,
              "stopTime"=stopTime, "logMetaEvalues"=logMetaEvalues,
              "logMetaEvaluesFut"=logMetaEvaluesFut,
              "nStopVec"=nStopVec)

  return(res)
}


scenario2T <- function(dat, allSources, designObj, alpha=0.05,
                       betaFutility=alpha, nuMin=3, nSim=1e2L,
                       nMax=NULL, seed=NULL, wantCi=FALSE) {

  nSources <- length(allSources)
  nSamples <- integer(nSources)

  nStop <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

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

    nParticipants <- nSamples[i] <- n1+n2

    nMax <- n1+n2

    set.seed(seed)
    for (k in 1:nSim) {
      tempRes <- twoSampleTTestRandomOrder(
        "x"=x, "y"=y, "n1"=n1, "n2"=n2,
        "designObj"=designObj, "nuMin"=nuMin,
        "alpha"=alpha, "betaFutility"=betaFutility,
        "wantCi"=wantCi, "nMax"=nMax
      )

      nStop[k, i] <- tempRes$nStop
      eValues[k, i] <- tempRes$eValue
      eValuesFut[k, i] <- tempRes$eValueFut
    }
  }
  res <- list(nStop=nStop, eValues=eValues, eValuesFut=eValuesFut, pValues=pValues)
  return(res)
}



scenario3T <- function(dat, allSources, designObj, alpha=0.05,
                       betaFutility=alpha, nuMin=3, nSim=1e2L,
                       nMax=NULL, seed=NULL, wantCi=FALSE) {

  nSamples <- nSources <- pValues <- length(allSources)
  nStop <- eValues <- eValuesFut <- matrix(nrow=nSim, ncol=nSources)

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

    tempResult <- t.test(x[1:n1], y[1:n2],
                         var.equal=stat.params$var.equal)
    pValues[i] <- tempResult$p.value

    nParticipants <- nSamples[i] <- n1+n2

    nMax <- if (is.null(nMax)) n1+n2 else sum(designObj$nPlan)

    set.seed(seed)
    for (k in 1:nSim) {
      tempRes <- twoSampleTTestRandomOrder(
        "x"=x, "y"=y, "n1"=n1, "n2"=n2,
        "designObj"=designObj, "nuMin"=nuMin,
        "alpha"=alpha, "betaFutility"=betaFutility,
        "wantCi"=wantCi, "nMax"=nMax
      )

      nStop[k, i] <- tempRes$nStop
      eValues[k, i] <- tempRes$eValue
      eValuesFut[k, i] <- tempRes$eValueFut
    }
  }
  res <- list(nStop=nStop, eValues=eValues, eValuesFut=eValuesFut, pValues=pValues)
  return(res)
}


