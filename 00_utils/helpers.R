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

  dataLeveled <- vector("list", length(levels(dat$factor)))
  badIndices <- vector("list", length(levels(dat$factor)))
  sampleSize <- matrix(nrow=length(allSources), ncol=length(levels(dat$factor)))

  for (i in seq_along(dataLeveled))
    dataLeveled[[i]] <- dat %>% dplyr::filter(factor==levels(dat$factor)[i])

  for (i in seq_along(allSources)) {
    someSource <- allSources[i]

    for (j in seq_along(dataLeveled)) {
      someDat <- dataLeveled[[j]]
      someDat <- someDat[someDat$source==someSource, ]
      sampleSize[i, j] <- dim(someDat)[1]
    }
  }

  allBadIndices <- integer(0)

  for (i in seq_along(badIndices)) {
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

  res <- list("allSources"=allSources, sampleSize=sampleSize)
  return(res)
}


checkXY <- function(x, y) {
  if (length(x) >= 1 && length(y) >= 1)
    return(TRUE)

  if (is.null(x))
    return(FALSE)

  if (is.null(y))
    return(FALSE)

  if (is.na(x))
    return(FALSE)

  if (is.na(y))
    return(FALSE)

  return(TRUE)
}
