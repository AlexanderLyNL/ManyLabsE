addUniqueIds <- function(ML2.var, ML2.df) {
  dat <- ML2.var[[1]]$cleanDataFilter
  sourceColumn <- character(length = dim(dat)[1])

  for (i in seq_along(sourceColumn)) {
    iets <- ML2.df[ML2.df$uID == dat$uID[i], ]
    sourceColumn[i] <- iets$source
  }

  dat$source <- sourceColumn
  return(dat)
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
