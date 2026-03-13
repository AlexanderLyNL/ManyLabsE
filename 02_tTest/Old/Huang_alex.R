library(devtools)
library(plyr)
library(rio)
library(tidyverse)

library(safestats)

# helper fns
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
  
  if (is.null(x))
    return(FALSE)
  
  if (is.null(y))
    return(FALSE)
  
  if (length(unique(x)) >= 2 && length(unique(y)) >= 2)
    return(TRUE)
  
  return(FALSE)
}





# source(file.path("~", "projects", "manyLabsE","02_tTest","t_test_functions.R"))

project.root <- file.path("~", "projects", "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))

# ANALYSIS INFO ----
study.description      <- 'Direction & SES (Huang et al., 2014)'
analysis.unique.id     <- 1
analysis.name          <- 'Huang.1'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
# subset -> subset.type to avoid conflicts
subset.type <- "all" # "sites"

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

#view(ML2.df)

# PREPARE DATA & OUTPUT ----
# Add a unique ID
ML2.df$uID <- seq(1, nrow(ML2.df))

# Get info to create a dataset for the current study
ML2.in <- get.info(ML2.key, colnames(ML2.df), subset.type)

# Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
# Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
ML2.id <- get.chain(ML2.in)

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2, 7, 150, 151, 157, 158, 521, 522, 523, 524, 525, 526, 527, 528, 529, 530, 531, 532, 535, 536, 537) %>%
  dplyr::filter(is.character(source))

# Decide which analyses to run on which groups
toRun <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

if (nrow(ML2.df) <= 0 || length(toRun$studiess) <= 0) {
  print("No tests to run, nothing selected!")
  stop()
}

#view(ML2.df)

# Create a variable indicating the study order for each case
ML2.df$study.order <- NA
stmp <- strsplit(ML2.df$StudyOrderN, "[|]")

# Correct differences in study names
# TODO: a function to fix this
Stud <- ML2.key$study.name

ML2.df$study.order <- plyr::laply(seq_along(stmp), function(o) {
  which(grepl(Stud, stmp[[o]])) %00% NA
})

ML2.sr <- list()
ML2.var <- list()
testVarEqual <- ML2.in$stat.params$var.equal

# analysis.type == 1 no for loop
gID <- rep(TRUE, nrow(ML2.df))
g <- 1

# Check nMin
if (sum(gID, na.rm = TRUE) >= Nmin.raw) {
  nMin1 <- TRUE
  # Get a list containing the data frames to be used in the analysis
  ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
}

# Double-check nMin
if (nMin1) {
  compN <- ML2.sr[[g]]$N
  compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
  compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
  if (any(compN >= Nmin.raw) & (all(compN1 >= Nmin.cond, compN2 >= Nmin.cond))) {
    nMin2 <- TRUE
  }
}

ML2.var[[g]] <- varfun.Huang.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params

t.test(variable ~ factor, data = ML2.var[[g]]$cleanDataFilter, var.equal = FALSE)

dat <- ML2.var[[g]]$cleanDataFilter

x <- dat$variable[dat$factor=="High"]
y <- dat$variable[dat$factor=="Low"]

t.test(x, y)

#write.csv(ML2.df, file = file.path(project.root,"02_tTest","huangData"), row.names = FALSE)

# Alexander -------
dat <- addUniqueIds(ML2.var, ML2.df)
dat <- checkUniqueIds(dat)

# Here -------
xName <- "High"
yName <- "Low"
alpha <- 0.05
betaFutility <- 0.2
deltaMin <- 0.68


allSources <- sort(unique(dat$source)) ## S: added this

nXs <- nYs <- integer(length(allSources))

for (i in seq_along(allSources)) {
  thisSource <- allSources[i]
  
  kip <- dat[dat$source==thisSource, ]
  nXs[i] <- sum(kip$factor==xName)
  nYs[i] <- sum(kip$factor==yName)
}

maxNX <- max(nXs)
maxNY <- max(nYs)


designObj <- designSaviT(deltaMin=deltaMin, power=0.8,
                         futility=TRUE, testType="twoSample",
                         varEqual=FALSE)

# Result containers
#   General data set attributes
n1Vec <- n2Vec <- ratios  <- numeric(length(allSources))


#   sample sizes for p-value based inference
n1VecFreq <- n2VecFreq <- pValues <- numeric(length(allSources))

#   sample sizes for e-value based inference
#
n1VecE <- n2VecE <- firstTimes <- eValues <- numeric(length(allSources))
n1VecEFut <- n2VecEFut <- firstTimesFut <- eValuesFut <- numeric(length(allSources))

allEValueVecs <- matrix(nrow=maxNX,
                        ncol=length(allSources))

allEValueVecsFut <- allEValueVecs


# Analyse data for each source
# loop start -----
for (i in 1:length(allSources)) {
  someDat <- dat[dat$source==allSources[i], ]
  
  ## Data -----
  x <- someDat[which(someDat$factor==xName), ]$variable
  y <- someDat[which(someDat$factor==yName), ]$variable
  
  # Remove non-available entries
  x <- x[!is.na(x)]
  n1 <- length(x)
  
  y <- y[!is.na(y)]
  n2 <- length(y)
  
  # Store valid sample size characteristics
  n1Vec[i] <- n1
  n2Vec[i] <- n2
  
  # ## Freq -----
  # n1Freq <- min(ceiling(freqDesign$n), n1)
  # n1VecFreq[i] <- n1Freq
  #
  # n2Freq <- min(ceiling(freqDesign$n), length(y))
  # n2VecFreq[i] <- n2Freq
  #
  # tempResult <- t.test(x[1:n1Freq], y[1:n2Freq],
  #                      var.equal=TRUE)
  # pValues[i] <- tempResult$p.value
  
  ## e-value ----
  n1EValue <- min(maxNX, n1)
  n2EValue <- min(maxNY, n2)
  
  n1VecE[i] <- n1EValue
  n2VecE[i] <- n2EValue
  
  ratios[i] <- n2EValue/n1EValue
  
  tempResult <- saviTTest(
    x[1:n1EValue],
    y[1:n2EValue],
    designObj=designObj, sequential=TRUE)
  
  eValues[i] <- max(tempResult$eValueVec, na.rm=TRUE)
  
  # Used to fill up an e-value sequence if there is too little data
  nLast <- length(tempResult$eValueVec)
  nRemaining <- maxNX - nLast
  
  if (nRemaining > 0) {
    tempResult$eValueVec <- c(tempResult$eValueVec, rep(tempResult$eValueVec[nLast], nRemaining))
    
    
    tempResult$eValueFutVec <- c(unlist(tempResult$eValueFutVec), rep(unlist(tempResult$eValueFutVec)[nLast], nRemaining))
  }
  
  allEValueVecs[, i] <- tempResult$eValueVec
  allEValueVecsFut[, i] <- tempResult$eValueFutVec
  
  firstTimes[i] <- min(which(tempResult$eValueVec >= 1/alpha))
  firstTimesFut[i] <- min(which(tempResult$eValueFutVec <= betaFutility))
}
# loop end ----

# added by S: after the first passage make all e-values equal to value at first passage time 
allStoppedEValueVecs <- allEValueVecs
allStoppedEValueVecsFut <- allEValueVecsFut
stopped_at <- pmin(firstTimes,firstTimesFut, nXs+nYs)

for (i in 1:length(allSources)) {
    allStoppedEValueVecs[stopped_at[i]:maxNX,i] <- allStoppedEValueVecs[stopped_at[i],i]
    allStoppedEValueVecsFut[stopped_at[i]:maxNX,i] <- allStoppedEValueVecsFut[stopped_at[i],i]
    }

stopped_at
firstTimes
firstTimesFut
pmin(firstTimes,firstTimesFut, nXs+nYs)
nXs
nYs
# Scenario 1 ----
# In the order of how the sources are mentioned, but can use randomisation
#   Also only up to n1=n1End and n2 = ratio*n1,
#   where ratio n2End/n1End, for instance when n1End > 47.
#   When n1End < 47, then the eValue is copied until n1=47
#   For instance:
#
#   print(allEValueVecs[, 3])
#
n1End <- dim(allEValueVecs)[1]

eMeta <- exp(cumsum(log(allStoppedEValueVecs[n1End, ])))
eFutMeta <- exp(cumsum(log(allStoppedEValueVecsFut[n1End, ])))

plot(log(eMeta), type="l")
lines(log(eFutMeta), col="red")

#eMetaAverage <- cumsum(allEValueVecs[n1End, ])/(1:length(allSources))
#eFutMetaAverage <- cumsum(allEValueVecsFut[n1End, ])/(1:length(allSources))

#plot(eMetaAverage, type="l", log="y")
#lines(eFutMetaAverage, col="red")

which(eMeta > 1/alpha)
#which(eMetaAverage > 1/alpha)

which(eFutMeta < betaFutility)
#which(eFutMetaAverage < betaFutility)

# Scenario 2 ----
eMeta <- exp(rowSums(log(allStoppedEValueVecs)))
eFutMeta <- exp(rowSums(log(allStoppedEValueVecsFut)))

plot(log(eMeta), type="l")
lines(eFutMeta, col="red")

#eMetaAverage <- rowMeans(allEValueVecs)
#eFutMetaAverage <- rowMeans(allEValueVecsFut)

#plot(eMetaAverage, type="l", log="y")
#lines(eFutMetaAverage, col="red")

which(eMeta >= 1/alpha)
#which(eMetaAverage > 1/alpha)

which(eFutMeta < betaFutility)
#which(eFutMetaAverage < betaFutility)

# Scenario 3 ---------
nTotal <- length(dat$uID)

set.seed(1)
dat <- dat[sample(nrow(dat)), ]
someOrder <- dat$uID

sourceDataTracker <- vector(mode="list", length=length(allSources))
names(sourceDataTracker) <- allSources

for (neem in allSources)
  sourceDataTracker[[neem]] <- list(x=NULL, y=NULL)

eFutCollection <- eCollection <- as.data.frame(matrix(ncol=length(allSources), nrow=nTotal))
names(eFutCollection) <- names(eCollection) <- allSources
eFutCollection[1, ] <- eCollection[1, ] <- 1


for (i in seq_along(someOrder)) {
  # for (i in 1:1000) {
  someId <- someOrder[i]
  
  someRow <- dat[which(dat$uID==someId), ]
  
  someSource <- someRow$source
  
  sourceDataTemp <- sourceDataTracker[[someSource]]
  
  x <- sourceDataTemp$x
  y <- sourceDataTemp$y
  
  if (someRow$factor==xName) {
    sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
  } else if (someRow$factor==yName) {
    sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
  }
  
  if (i > 1) {
    someCheck <- checkXY(x, y)
    
    if (someCheck) {
      tempRes <- saviTTest(x, y, designObj=designObj, sequential=FALSE)
      
      eCollection[[someSource]][i] <- tempRes$eValue
      eFutCollection[[someSource]][i] <- tempRes$eValueFut
    } else {
      eCollection[[someSource]][i] <- 1
      eFutCollection[[someSource]][i] <- 1
    }
    
    for (source in allSources) {
      if (source!=someSource) {
        eCollection[[source]][i] <- eCollection[[source]][i-1]
        eFutCollection[[source]][i] <- eFutCollection[[source]][i-1]
      }
    }
  }
}

eMatrix <- as.matrix(eCollection)
eFutMatrix <- as.matrix(eFutCollection)

eMeta <- exp(rowSums(log(eMatrix)))
eFutMeta <- exp(rowSums(log(eFutMatrix)))

plot(1:nTotal, log(eMeta), type="l")
lines(1:nTotal, log(eFutMeta), col="red")


#eMetaAverage <- rowMeans(eMatrix)
#eFutMetaAverage <- rowMeans(eFutMatrix)

#plot(1:nTotal, eMetaAverage, type="l")
#lines(1:nTotal, eFutMetaAverage, col="red")

which(eMeta >= 1/alpha)
which(eMetaAverage >= 1/alpha)

which(eFutMeta <= betaFutility)
which(eFutMetaAverage <= betaFutility)