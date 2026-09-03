# Code written by Alexander Ly and Udo Boehm
#   to reproduce
#
#   Figure 1
#
#     of
#
#  Ly, A., Boehm, U., Grunwald, P.D., Ramdas, A., van Ravenzwaaij, D. (2025). A Tutorial on Safe Anytime-Valid Inference: Practical Maximally Flexible Sampling Designs for Experiments Based on e-Values. psyArxiv preprint.

# Gray -------
# Data retrieved from
#
#   https://github.com/ManyLabsOpenScience/ManyLabs2/tree/master/OSFdata/Moral%20Typecasting%20(Gray%20%26%20Wegner%2C%202009)
#
# Data preparation using
#
#   Gray_1_study_global_include_all.R
#
# to generate relevant data frames
#
#   See also: https://github.com/ManyLabsOpenScience/ManyLabs2
#     Folder
#       OSFData/
#
#   and
#       OSFData/Moral Typecasting (Gray & Wegner, 2009)/Gray.1/Global/Gray_1_study_global_include_all.R
#
#
# For convenience, the data frame is saved as grayData.RData
#
#


# To install the package 0.88 when it's not yet on cran
#
# remotes::install_github("AlexanderLyNL/safestats", ref = "088")
#
# library(safestats)

overColour <- "#A6CEE380"
overColourBorder <- "#1F78B4E6"
underColour <- "#FFB90F86"
underColourBorder <- "#FFB90FCC"
continueColour <- "#556B2F4D"
continueColourBorder <- "#556B2FCC"

# eValueColour <- "#A6CEE380"
# eValueColourBorder <- "#1F78B4E6"
#
# underColour <- adjustcolor("darkolivegreen", alpha.f=0.3)
# underColourBorder <- adjustcolor("darkolivegreen", alpha.f=0.8)
# overColour <- adjustcolor("#DAA52066", alpha.f=0.8)
# overColourBorder <- "#DAA52066"

# library("safestats")
pdfWidth <- 14
pdfHeight <- 7

cexFactor <- 1.3
myCexAxis <- 2.25

alpha <- 0.05
alphaRelevance <- alpha

# Gray data --------

# Original study: Lower bound of effect size found in the original study
deltaMin <- (5.29-3.86)/1.86



# Prospective frequentist analysis
freqDesign <- power.t.test(delta=deltaMin, power=0.8,
                           alternative="one.sided")

# Prospective e-value analysis
# designObj1 <- designSaviT(deltaMin=deltaMin, beta=0.2, seed=1,
#                           testType="twoSample", alternative="greater")
designObj1 <- designSaviT(deltaMin=deltaMin,
                          power=0.8, seed=2,
                          testType="twoSample",
                          varEqual=FALSE, relevanceTest=TRUE)

plot(designObj1)

#   Data retrieved from the ManyLabs2 project, see Appendix below
#
# load("~/dropbox/projects/savitutorial/code/grayData.RData")

sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"


load(paste0(myWd, "grayData.RData"))


which(!is.na(grayData$gray1.2))
which(!is.na(grayData$gray2.2))




allSources <- unique(grayData$source)

# Result containers
#   General data set attributes
n1Vec <- n2Vec <- ratios  <- numeric(length(allSources))


#   sample sizes for p-value based inference
n1VecFreq <- n2VecFreq <- pValues <- numeric(length(allSources))

#   sample sizes for e-value based inference
#
n1VecE <- n2VecE <- firstTimes <- eValues <- numeric(length(allSources))
n1VecERelevance <- n2VecERelevance <- firstTimesFut <- eRelevance <- numeric(length(allSources))

allEValueVecs <- matrix(nrow=designObj1$nPlan[1],
                        ncol=length(allSources))
allEValueVecsFut <- allEValueVecs


# i <- 53
# Loop start -----
for (i in seq_along(eValues)) {
  someDat <- grayData[grayData$source==allSources[i], ]

  ## Data -----
  x <- someDat$gray1.2
  y <- someDat$gray2.2

  # Remove non-available entries
  x <- x[!is.na(x)]
  n1 <- length(x)

  y <- y[!is.na(y)]
  n2 <- length(y)

  # Store valid sample size characteristics
  n1Vec[i] <- n1
  n2Vec[i] <- n2

  ## Freq -----
  n1Freq <- min(ceiling(freqDesign$n), n1)
  n1VecFreq[i] <- n1Freq

  n2Freq <- min(ceiling(freqDesign$n), length(y))
  n2VecFreq[i] <- n2Freq

  tempResult <- t.test(x[1:n1Freq], y[1:n2Freq], var.equal=TRUE)
  pValues[i] <- tempResult$p.value

  ## e-Value ----
  n1EValue <- min(designObj1$nPlan[1], n1)
  n2EValue <- min(designObj1$nPlan[2], n2)

  n1VecE[i] <- n1EValue
  n2VecE[i] <- n2EValue

  ratios[i] <- n2EValue/n1EValue


  tempResult <- saviTTest(
    x[1:n1EValue], y[1:n2EValue], designObj=designObj1,
    sequential=TRUE, nuMin=2)

  eValues[i] <- max(tempResult$eValueVec, na.rm=TRUE)

  eRelevance[i] <- min(tempResult$eValuERelevanceVec, na.rm=TRUE)

  # Used to fill up an e-value sequence if there is too little data
  nLast <- length(tempResult$eValueVec)
  nRemaining <- designObj1$nPlan[1] - nLast

  if (nRemaining > 0) {
    tempResult$eValueVec <- c(tempResult$eValueVec, rep(tempResult$eValueVec[nLast], nRemaining))
    tempResult$eValuERelevanceVec <- c(unlist(tempResult$eValuERelevanceVec), rep(unlist(tempResult$eValuERelevanceVec)[nLast], nRemaining))
  }

  allEValueVecs[, i] <- tempResult$eValueVec
  allEValueVecsFut[, i] <- tempResult$eValuERelevanceVec

  firstTimes[i] <- min(which(tempResult$eValueVec >= 20))
  firstTimesFut[i] <- min(which(tempResult$eValuERelevanceVec <= alphaRelevance))
}


firstTimes
firstTimesFut


# Scenario 1 ----
# In the order of how the sources are mentioned, but can use randomisation
#   Also only up to n1=47 and n2 = ratio*n1,
#   where ratio n2End/n1End, for instance when n1End > 47.
#   When n1End < 47, then the eValue is copied until n1=47
#   For instance:
#
#   print(allEValueVecs[, 3])
#

eMeta <- exp(cumsum(log(allEValueVecs[47, ])))
ERelevanceMeta <- exp(cumsum(log(allEValueVecsFut[47, ])))

plot(eMeta, type="l", log="y")
lines(ERelevanceMeta, col="red")

eMetaAverage <- cumsum(allEValueVecs[47, ])/(1:length(allSources))
ERelevanceMetaAverage <- cumsum(allEValueVecsFut[47, ])/(1:length(allSources))

plot(eMetaAverage, type="l", log="y")
lines(ERelevanceMetaAverage, col="red")

which(ERelevanceMetaAverage<0.2)

# Scenario 2 ----
eMeta <- exp(rowSums(log(allEValueVecs)))
ERelevanceMeta <- exp(rowSums(log(allEValueVecsFut)))

plot(eMeta, type="l", log="y")
lines(ERelevanceMeta, col="red")

which(ERelevanceMeta < 0.2)

eMetaAverage <- rowMeans(allEValueVecs)
ERelevanceMetaAverage <- rowMeans(allEValueVecsFut)

plot(eMetaAverage, type="l", log="y")
lines(ERelevanceMetaAverage, col="red")

which(ERelevanceMetaAverage < 0.2)

# Scenario 3 ---------
xNa <- which(is.na(grayData$gray1.2))
yNa <- which(is.na(grayData$gray2.2))

grayDataCleaned <- grayData[-intersect(xNa, yNa), ]

nTotal <- length(unique(grayDataCleaned$uID))

set.seed(1)
someOrder <- sample(unique(grayDataCleaned$uID), nTotal)

sourceDataTracker <- vector(mode="list", length=length(allSources))
names(sourceDataTracker) <- allSources

for (neem in allSources)
  sourceDataTracker[[neem]] <- list(x=NULL, y=NULL)

ERelevanceCollection <- eCollection <- as.data.frame(matrix(ncol=length(allSources), nrow=nTotal))
names(ERelevanceCollection) <- names(eCollection) <- allSources
ERelevanceCollection[1, ] <- eCollection[1, ] <- 1


for (i in seq_along(someOrder)) {
  # for (i in 1:1000) {
  someId <- someOrder[i]

  someRow <- grayDataCleaned[which(grayDataCleaned$uID==someId), ]

  someSource <- someRow$source

  sourceDataTemp <- sourceDataTracker[[someSource]]

  x <- sourceDataTemp$x
  y <- sourceDataTemp$y

  xTemp <- someRow$gray1.2
  yTemp <- someRow$gray2.2

  if (!is.na(xTemp)) {
    sourceDataTracker[[someSource]]$x <- x <- c(x, xTemp)
  } else {
    yTemp <- someRow$gray2.2
    sourceDataTracker[[someSource]]$y <- y <- c(y, yTemp)
  }

  if (i >1) {
    someCheck <- checkXY(x, y)

    if (someCheck) {
      tempRes <- saviTTest(x, y, designObj=designObj1, sequential=FALSE)
      eCollection[[someSource]][i] <- tempRes$eValue
      ERelevanceCollection[[someSource]][i] <- tempRes$eValuERelevance
    } else {
      eCollection[[someSource]][i] <- 1
      ERelevanceCollection[[someSource]][i] <- 1
    }

    for (source in allSources) {
      if (source!=someSource) {
        eCollection[[source]][i] <- eCollection[[source]][i-1]
        ERelevanceCollection[[source]][i] <- ERelevanceCollection[[source]][i-1]
      }
    }
  }
}

eMatrix <- as.matrix(eCollection)
ERelevanceMatrix <- as.matrix(ERelevanceCollection)

eMeta <- exp(rowSums(log(eMatrix)))
ERelevanceMeta <- exp(rowSums(log(ERelevanceMatrix)))

plot(1:nTotal, eMeta, type="l", log="y")
lines(1:nTotal, ERelevanceMeta, col="red")

which(ERelevanceMeta < 0.2)



eMetaAverage <- rowMeans(eMatrix)
ERelevanceMetaAverage <- rowMeans(ERelevanceMatrix)

plot(1:nTotal, eMetaAverage, type="l", log="y")
lines(1:nTotal, ERelevanceMetaAverage, col="red")

which(ERelevanceMetaAverage < 0.2)
