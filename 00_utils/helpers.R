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
    betaFutility=alpha, nMax=n1+n2) {

  nParticipants <- n1+n2
  someOrder <- sample(nParticipants, nParticipants)

  nMax <- min(nMax, nParticipants)

  xRun <- numeric(0)
  yRun <- numeric(0)

  totalVar <- c(x, y)

  for (j in seq_along(someOrder)) {
    partId <- someOrder[j]

    if (partId <= n1) {
      xRun <- c(xRun, totalVar[someOrder[j]])
    } else {
      yRun <- c(yRun, totalVar[someOrder[j]])
    }

    someCheck <- checkXY(xRun, yRun)

    if (someCheck) {
      tempRes <- saviTTest(
        xRun, yRun, designObj=designObj,
        sequential=FALSE, nuMin=nuMin)

      eNow <- tempRes$eValue
      eFutNow <- tempRes$eValueFut

      if (eNow >= 1/alpha || eFutNow <= betaFutility ||
          j==nMax) {
        res <- list(nStop=j, eValue=eNow, eValueFut=eFutNow)
        return(res)
      }
    }
  }
}

