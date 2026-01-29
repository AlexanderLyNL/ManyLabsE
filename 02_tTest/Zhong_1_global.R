# New function to be added to the package in adapted form
#
# remotes::install_github("AlexanderLyNL/safestats", ref = “futility88")
# library("safestats")
#

saviTTestFut <- function(
    x, y=NULL, designObj=NULL, paired=FALSE,
    varEqual=TRUE, ciValue=NULL,
    maxRoot=10, sequential=NULL,
    futility=FALSE, esMinFutility=NULL, ...) {

  result <- constructSaviTestObj("T-Test")

  # Vars for sequential analysis
  eValueVec <- NULL
  confSeqMatrix <- NULL
  n1Vec <- NULL
  n2Vec <- NULL

  ## Def: test type -------
  if (is.null(y)) {
    testType <- "oneSample"
  } else {
    if (paired) {
      testType <- "paired"
    } else {
      testType <- "twoSample"
    }
  }

  ## Check: designObj ----
  if (is.null(designObj)) {
    designObj <- designSaviT(0.5, "eType"="mom",
                             "testType"=testType)
    designObj[["pilot"]] <- TRUE

    warningMessage <- paste("No designObj given. Default test computed based",
                            "on a non-local prior at +1/2 and -1/2.")
    warning(warningMessage)
  }

  if (designObj[["testName"]] != "T-Test")
    warning("The provided design is not constructed for the t-test,",
            "please use designSaviT() instead. The test results might be invalid.")

  if (designObj[["testType"]] != testType)
    warning('The test type of designObj is "', designObj[["testType"]],
            '", whereas the data correspond to a testType "', testType, '"')

  ## Check: Data -----
  #
  if (is.null(y)) {
    ### One-sample -----
    #
    if (isTRUE(paired))
      stop("Data error: Paired analysis requested without specifying the second variable")

    dataName <- deparse1(substitute(x))
    x <- x[!is.na(x)]

    n <- nEff <- n1 <- length(x)
    n2 <- NULL
    nu <- n-1

    meanObs <- estimate <- mean(x)
    sdObs <- stats::sd(x)

    names(estimate) <- "mean of x"
    names(n) <- "n1"

    if (is.null(sequential))
      sequential <- if (n1 <= 200) TRUE else FALSE

    if (sequential) {
      tempN <- defineTTestN("lowN"=1, "highN"=n1,
                            "testType"="oneSample")

      nEffVec <- tempN[["nEff"]]
      n1Vec <- tempN[["n1"]]
      n2Vec <- tempN[["n2"]]
      nuVec <- tempN[["nu"]]

      meanObsVec <- 1/nEffVec*cumsum(x)
      sdObsVec <- sqrt(1/nuVec*(cumsum(x^2)-nEffVec*meanObsVec^2))
    }
  } else {
    dataName <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))

    if (isTRUE(paired))
      xGoodIndeces <- yGoodIndeces  <-
        stats::complete.cases(x, y)
    else {
      yGoodIndeces <- !is.na(y)
      xGoodIndeces <- !is.na(x)
    }

    x <- x[xGoodIndeces]
    y <- y[yGoodIndeces]

    n1 <- length(x)
    n2 <- length(y)

    ### Paired ----
    #
    if (isTRUE(paired)) {
      if (n1 != n2)
        stop("Data error: Error in complete.cases(x, y): Paired analysis requested, ",
             "but the two samples are not of the same size.")

      nEff <- n1
      nu <- n1-1
      meanObs <- estimate <- mean(x-y)
      sdObs <- stats::sd(x-y)
      names(estimate) <- "mean of the differences"

      if (is.null(sequential))
        sequential <- if (n1 <= 200) TRUE else FALSE

      if (sequential) {
        tempN <- defineTTestN("lowN"=1, "highN"=n1, testType="paired")

        nEffVec <- tempN[["nEff"]]
        n1Vec <- tempN[["n1"]]
        n2Vec <- tempN[["n2"]]
        nuVec <- tempN[["nu"]]

        meanObsVec <- 1/nEffVec*cumsum(x-y)
        sdObsVec <- sqrt(1/nuVec*(cumsum((x-y)^2)-nEffVec*meanObsVec^2))
      }
    } else {
      ## Two-sample ----
      nEff <- (1/n1+1/n2)^(-1)
      nu <- n1+n2-2

      sPooledSquared <- ((n1-1)*stats::var(x)+(n2-1)*stats::var(y))/nu

      sdObs <- sqrt(sPooledSquared)

      estimate <- c(mean(x), mean(y))
      names(estimate) <- c("mean of x", "mean of y")
      meanObs <- estimate[1]-estimate[2]

      if (is.null(sequential))
        sequential <- if (n1 <= 200) TRUE else FALSE

      if (sequential) {
        tempN <- defineTTestN(1, n1, n2/n1, testType="twoSample")

        nEffVec <- tempN[["nEff"]]
        nuVec <- tempN[["nu"]]

        # These now serve as an order
        n1Vec <- tempN[["n1"]]
        n2Vec <- tempN[["n2"]]

        xMeanObsRaw <- 1/(1:n1)*cumsum(x)
        yMeanObsRaw <- 1/(1:n2)*cumsum(y)

        xSumsOfSquaresRaw <- (cumsum(x^2)-(1:n1)*xMeanObsRaw^2)
        ySumsOfSquaresRaw <- (cumsum(y^2)-(1:n2)*yMeanObsRaw^2)

        if (n2/n1==1) {
          xMeanObsVec <- xMeanObsRaw
          yMeanObsVec <- yMeanObsRaw
          xSumsOfSquaresVec <- xSumsOfSquaresRaw
          ySumsOfSquaresVec <- ySumsOfSquaresRaw
        } else {
          vecLength <- length(n1Vec)

          xMeanObsVec <- yMeanObsVec <-
            xSumsOfSquaresVec <- ySumsOfSquaresVec <- numeric(vecLength)

          for (j in 1:vecLength) {
            nowN1 <- n1Vec[j]
            nowN2 <- n2Vec[j]

            xMeanObsVec[j] <- xMeanObsRaw[nowN1]
            yMeanObsVec[j] <- yMeanObsRaw[nowN2]
            xSumsOfSquaresVec[j] <- xSumsOfSquaresRaw[nowN1]
            ySumsOfSquaresVec[j] <- ySumsOfSquaresRaw[nowN2]
          }
        }

        sPooledSquaredVec <- (xSumsOfSquaresVec+ySumsOfSquaresVec)/nuVec

        meanObsVec <- xMeanObsVec-yMeanObsVec
        sdObsVec <- sqrt(sPooledSquaredVec)
      }
    }

    n <- c(n1, n2)
    names(n) <- c("n1", "n2")
  }

  alpha <- designObj[["alpha"]]
  alternative <- designObj[["alternative"]]
  h0 <- designObj[["h0"]]

  if (is.null(ciValue))
    ciValue <- 1-alpha

  if (ciValue < 0 || ciValue > 1)
    stop("Can't make a confidence sequence with ciValue < 0 or ciValue > 1, or alpha < 0 or alpha > 1")

  tStat <- tryOrFailWithNA(sqrt(nEff)*(meanObs - h0)/sdObs)

  if (is.na(tStat))
    stop("Data error: Could not compute the t-statistic")

  names(tStat) <- "t"

  ### Compute: eValue ----
  #
  testResult <- suppressWarnings(
    saviTTestStat("t"=tStat, "parameter"=designObj[["parameter"]], "n1"=n1,
                  "n2"=n2, "alternative"=alternative, "paired"=paired,
                  "eType"=designObj[["eType"]])
  )


  ### Compute: confSeq ----
  #
  result[["confSeq"]] <- computeConfidenceIntervalT(
    "meanObs"=meanObs, "sdObs"=sdObs,
    "nEff"=nEff, "nu"=nu,
    "parameter"=designObj[["parameter"]],
    "eType"=designObj[["eType"]], "ciValue"=ciValue, "maxRoot"=maxRoot)

  eValueVecFut <- NULL

  ## Compute: Sequential ----
  if (sequential) {
    tStatVec <- sqrt(nEffVec)*(meanObsVec-h0)/sdObsVec

    mIter <- length(n1Vec)

    eValueVec <- numeric(mIter)
    confSeqMatrix <- matrix(nrow=mIter, ncol=2)

    if (futility)
      eValueVecFut <- eValueVec

    for (i in seq_along(n1Vec)) {
      brie <- suppressWarnings(
        saviTTestStat("t"=tStatVec[i], "parameter"=designObj[["parameter"]],
                      "n1"=n1Vec[i], "n2"=n2Vec[i], "alternative"=alternative,
                      "paired"=paired, "eType"=designObj[["eType"]])
      )

      if (futility) {
        eValueVecFut[i] <- suppressWarnings(unlist(
          saviFutilityTStatNEffNu(
            "t"=tStatVec[i], "nEff"=nEffVec[i],
            "nu"=nuVec[i], "esMinFutility"=esMinFutility,
            "alternative"="twoSided"))
        )
      }

      eValueVec[i] <- unname(brie[["eValue"]])

      kaas <- computeConfidenceIntervalT("meanObs"=meanObsVec[i], "sdObs"=sdObsVec[i],
                                         "nEff"=nEffVec[i], "nu"=nuVec[i],
                                         "parameter"=designObj[["parameter"]],
                                         "eType"=designObj[["eType"]], "ciValue"=ciValue,
                                         "maxRoot"=maxRoot)

      confSeqMatrix[i, ] <- kaas
    }
  }

  ### Fill: Result -----
  #
  result[["statistic"]] <- tStat
  result[["estimate"]] <- estimate
  result[["stderr"]] <- sdObs/sqrt(nEff)
  result[["dataName"]] <- dataName
  result[["designObj"]] <- designObj
  result[["testType"]] <- testType
  result[["n"]] <- n
  result[["ciValue"]] <- ciValue

  result[["eValueVec"]] <- eValueVec
  result[["confSeqMatrix"]] <- confSeqMatrix
  result[["n1Vec"]] <- n1Vec
  result[["n2Vec"]] <- n2Vec

  result[["eValue"]] <- testResult[["eValue"]]
  result[["eValueApproxError"]] <- testResult[["eValueApproxError"]]

  result[["eValueVecFut"]] <- eValueVecFut

  names(result[["statistic"]]) <- "t"

  return(result)
}

# SETUP ENVIRONMENT ----
library(devtools)
library(plyr)
library(rio)
library(tidyverse)
library(reshape2)
library(stats)

# TODO: Set up your directory
# If


sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))

# ANALYSIS INFO ----
study.description <- "Moral Cleansing (Zhong & Liljenquist, 2006)"
analysis.unique.id <- 65
analysis.name <- "Zhong.1"
analysis.type <- 1
analysis.type.name <- "study_global_include"
analysis.type.groups <- "Source.Global"
Nmin.raw <- 30
Nmin.cond <- 15
# subset -> subset.type to avoid conflicts
subset.type <- "all"

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

# Add a unique ID
ML2.df$uID <- seq(1, nrow(ML2.df))

# Get info to create a dataset for the current study
# keytable <- ML2.key
ML2.in <- get.info(ML2.key, colnames(ML2.df), subset.type)

# Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
# Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
ML2.id <- get.chain(ML2.in)

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2, 7, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 805, 904, 905, 906, 907, 908, 909, 910, 911, 912, 913, 914, 915, 934, 935, 938, 939, 940) %>%
  dplyr::filter(is.character(source))

# Decide which analyses to run on which groups
toRun <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

if (nrow(ML2.df) <= 0 || length(toRun$studiess) <= 0) {
  print("No tests to run, nothing selected!")
  stop()
}



# Create a variable indicating the study order for each case
ML2.df$study.order <- NA
stmp <- strsplit(ML2.df$StudyOrderN, "[|]")

# Correct differences in study names
Stud <- ML2.key$study.name

ML2.df$study.order <- plyr::laply(seq_along(stmp), function(o) {
  which(grepl(Stud, stmp[[o]])) %00% NA
})

ML2.sr <- list()
ML2.var <- list()
outputSource <- list()
dataSource <- list()
raw.df <- list()
clean.df <- list()
cleanData <- list()
testVarEqual <- ML2.in$stat.params$var.equal

# Loop over sites in runGroups within a study
if (analysis.type == 1) {
  runGroups <- "all"
} else {
  runGroups <- sort(na.exclude(unique(ML2.df[[toRun$ugroup]])))
}

disp(paste(analysis.unique.id, ML2.key$study.analysis, "- START"), header = toupper(ML2.key$study.analysis), footer = FALSE)
cat("\n")


# START GROUPS ----
g <- 1

# Include only datasets that have N >= Nmin.raw & n.group >= Nmin.cond
listIT <- FALSE
nMin1 <- FALSE
nMin2 <- FALSE
compN <- compN1 <- compN2 <- 0

gID <- rep(TRUE, nrow(ML2.df))

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

# START ANALYSIS ----------------------------------------

if (all(nMin1, nMin2)) {
  # To see the function code type:varfun.Zhong.1, or lookup in manylabRs_SOURCE.R
  ML2.var[[g]] <- varfun.Zhong.1(ML2.sr[[g]])


  # Check equal variance assumption
  if (!is.na(testVarEqual)) {
    if (testVarEqual) {
      logtxt <- paste(analysis.unique.id, ML2.key$study.analysis, "-", runGroups[g])
      ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]], ML2.in$study.vars.labels, ML2.key, group = logtxt) # don't pass the cleanData frame
    }
  }

  # Run the analysis according to ML2.key: 'stat.test'
  stat.params <<- ML2.in$stat.params


  stat.test <- try.CATCH(with(ML2.var[[g]], t.test(x = Ethical, y = Unethical, conf.level = stat.params$conf.level, var.equal = stat.params$var.equal, alternative = stat.params$alternative)))
}

# Alexander -----

sourceColumn <- character(length = dim(ML2.var[[1]]$cleanDataFilter)[1])

for (i in seq_along(sourceColumn)) {
  iets <- ML2.df[ML2.df$uID == ML2.var[[1]]$cleanDataFilter$uID[i], ]
  sourceColumn[i] <- iets$source
}

ML2.var[[1]]$cleanDataFilter$source <- sourceColumn

zhongData <- ML2.var[[1]]


# ASDF---------

# remotes::install_github("AlexanderLyNL/safestats", ref = “futility88")

library(safestats)

eValueColour <- "#A6CEE380"
eValueColourBorder <- "#1F78B4E6"

underColour <- adjustcolor("darkolivegreen", alpha.f=0.3)
underColourBorder <- adjustcolor("darkolivegreen", alpha.f=0.8)
overColour <- adjustcolor("#DAA52066", alpha.f=0.8)
overColourBorder <- "#DAA52066"


# library("safestats")
pdfWidth <- 14
pdfHeight <- 7

cexFactor <- 1.3
myCexAxis <- 2.25
myCex <- 1.5

betaFutility <- 0.2


# Original study: Lower bound of effect size found in the original study
deltaMin <- (4.95-3.75)/1.32

# Prospective frequentist analysis
freqDesign <- power.t.test(delta=deltaMin, alternative="one.sided",
                           power=0.8)

# Prospective e-value analysis
designObj2 <- designSaviT(deltaMin=deltaMin, beta=0.2,
                          testType="twoSample",
                          alternative="greater", seed=1)

designObj2Two <- designSaviT(deltaMin=deltaMin, beta=0.2,
                             testType="twoSample",
                             alternative="twoSided", seed=5)
# designObj2 <- designObj2Two

# Zhong data --------
#
#   Data taken from: https://github.com/ManyLabsOpenScience/ManyLabs2
#     Folder
#       OSFData/
#
#   Extract based on
#     OSFData/Moral Cleansing (Zhong & Liljenquist, 2006)/Zhong.1/Global/Zhong_1_study_global_include_all.R
#
#
#
# load(paste0(myWd, "zhongData.RData"))

dat <- zhongData$cleanDataFilter

head(dat)

allSources <- unique(dat$source)



# Remove "bogota" data, which only has participants i "Ethical" condition
allSources <- allSources[-which(allSources=="bogota")]

dat <- dat[dat$source %in% allSources, ]



# Result containers
#   General data set attributes
n1Vec <- n2Vec <- ratios  <- numeric(length(allSources))


#   sample sizes for p-value based inference
n1VecFreq <- n2VecFreq <- pValues <- numeric(length(allSources))

#   sample sizes for e-value based inference
#
n1VecE <- n2VecE <- firstTimes <- eValues <- numeric(length(allSources))
n1VecEFut <- n2VecEFut <- firstTimesFut <- eValuesFut <- numeric(length(allSources))

allEValueVecs <- matrix(nrow=designObj2$nPlan[1],
                        ncol=length(allSources))

allEValueVecsFut <- allEValueVecs


# Analyse data for each source
# loop start -----
for (i in 1:length(allSources)) {
  someDat <- dat[dat$source==allSources[i], ]

  ## Data -----
  x <- someDat[which(someDat$factor=="Ethical"), ]$variable
  y <- someDat[which(someDat$factor=="Unethical"), ]$variable

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

  ## e-value ----
  n1EValue <- min(designObj2$nPlan[1], n1)
  n2EValue <- min(designObj2$nPlan[2], n2)

  n1VecE[i] <- n1EValue
  n2VecE[i] <- n2EValue

  ratios[i] <- n2EValue/n1EValue

  tempResult <- saviTTestFut(
    x[1:n1EValue],
    y[1:n2EValue],
    designObj=designObj2, sequential=TRUE,
    futility=TRUE, esMinFutility=deltaMin)

  eValues[i] <- max(tempResult$eValueVec, na.rm=TRUE)

  # Used to fill up an e-value sequence if there is too little data
  nLast <- length(tempResult$eValueVec)
  nRemaining <- designObj2$nPlan[1] - nLast

  if (nRemaining > 0) {
    tempResult$eValueVec <- c(tempResult$eValueVec, rep(tempResult$eValueVec[nLast], nRemaining))


    tempResult$eValueVecFut <- c(unlist(tempResult$eValueVecFut), rep(unlist(tempResult$eValueVecFut)[nLast], nRemaining))
  }

  allEValueVecs[, i] <- tempResult$eValueVec
  allEValueVecsFut[, i] <- tempResult$eValueVecFut

  firstTimes[i] <- min(which(tempResult$eValueVec >= 20))
  firstTimesFut[i] <- min(which(tempResult$eValueVecFut <= betaFutility))
}
# loop end ----

# Futilty HERE new -------

firstTimes
firstTimesFut

plot(1:57, exp(cumsum(log(allEValueVecsFut[27, ]))), type="l", log="y")



# Freq result ----
sum(pValues < 0.05)
mean(pValues < 0.05)

# Number of samples used in the frequentist analysis
sum(n1VecFreq)
sum(n2VecFreq)

# Percentage saved by frequentist analysis
(sum(n1Vec)-sum(n1VecFreq))/sum(n1Vec)
(sum(n2Vec)-sum(n2VecFreq))/sum(n2Vec)


# e-value result ----
sum(eValues > 20)
mean(eValues > 20)


# Analysis of stopping time -----
n1Fpt <- firstTimes
n2Fpt <- ceil(firstTimes*ratios)

notStoppedIndex <- which(is.infinite(firstTimes))
n1Fpt[notStoppedIndex] <- n1VecE[notStoppedIndex]
n2Fpt[notStoppedIndex] <- n2VecE[notStoppedIndex]

# Average sample size used by e-value analysis
mean(n1Fpt)
mean(n2Fpt)

# Percentage saved by e-value analysis
(sum(n1Vec)-sum(n1Fpt))/sum(n1Vec)
(sum(n2Vec)-sum(n2Fpt))/sum(n2Vec)


# Plot-----
stoppedTimes <- firstTimes
stoppedTimes[is.infinite(firstTimes)] <- n1VecE[is.infinite(firstTimes)]

fptHist <- hist(stoppedTimes, plot=FALSE,
                breaks=1:designObj2$nPlan[1])

y <- fptHist[["density"]]
nB <- length(fptHist$breaks)
yRange <- range(y, 0)

alpha <- 0.05
ylim <- c(-1*log(20/(2*alpha)), 2.75*log(1/alpha))

someConstant <- (ylim[2]+log(alpha))/yRange[2]
textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

xlim <- c(0, designObj2$nPlan[1])


notStoppedTable <- table(stoppedTimes[which(is.infinite(firstTimes))])


notStoppedN <- as.integer(names(notStoppedTable))
notStoppedBottom <- rep(0, length(notStoppedN))
notStoppedTop <- notStoppedTable/57

overIndexes <- which(is.finite(firstTimes))
underIndexes <- which(is.infinite(firstTimes))


myName <- "zhongExample"
pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)

graphics::par(cex.main=1.5, mar=c(6, 6, 4, 0)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
              font.lab=2, cex.axis=1.3, bty="n", las=1)

plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "",
     cex.lab = 1.3, cex.axis = 1.3, las = 1, main=NULL,
     xaxt = "n", yaxt = "n", bty = "n", type = "p", pch = 15,
     bg = "grey")


abline(h = log(1), col = "darkgrey", lwd = 2, lty = 2)
abline(h = log(1/alpha))

criticalP <- log(c(alpha/10, alpha, 1, 1/alpha))

axis(side = 2, at = c(criticalP), tick = TRUE, las = 2, cex.axis = 1.3,
     labels = c(alpha/10, alpha, "1", 1/alpha), cex.axis=myCexAxis)
# axis(side = 1)
axis(side = 1, at=c(0, 10*(1:3)), cex.axis=myCexAxis)

ylab <- "Evidence"

mtext(ylab, side = 2, line = 2.5, las = 0, cex = cexFactor*myCex,
      adj=0.5, padj=-0.5)

xlab <- "Sample size"

mtext(xlab, side = 1, line = 2.5, las = 1,
      cex = cexFactor*myCex, padj=0.5)


rect(fptHist$breaks[-nB]+0.5, log(1/alpha),
     fptHist$breaks[-1L]+0.5, someConstant*y+log(1/alpha),
     col = eValueColour, border = eValueColourBorder, lwd=2,
     angle = 45, density = NULL, lty = NULL)

for (i in seq_along(notStoppedN)) {
  tempN <- notStoppedN[i]

  rect(xleft=tempN-0.5, ybottom=someConstant*notStoppedBottom[i]+log(1/alpha),
       xright=tempN+0.5, ytop=someConstant*notStoppedTop[i]+log(1/alpha),
       col = underColour, lwd=2, border=underColourBorder,
       angle = 45, density = NULL, lty = NULL)

}

for (j in underIndexes) {
  n1Temp <- stoppedTimes[j]

  lines(1:n1Temp, log(allEValueVecs[1:n1Temp, j]),
        lwd=2, col=underColour)
}
# j <- 26

for (j in underIndexes) {
  n1Temp <- stoppedTimes[j]

  points(n1Temp,
         log(allEValueVecs[n1Temp, j]),
         pch=15, col=underColourBorder)
}

for (i in overIndexes) {
  n1Temp <- stoppedTimes[i]

  lines(1:n1Temp, c(log(allEValueVecs[1:(n1Temp-1), i]), log(1/alpha)),
        lwd=5, col=overColourBorder)

  points(n1Temp,
         log(1/alpha),
         pch=15, col=overColourBorder)
}

mtext("e-value analyses of ManyLabs2 replications of Zhong & Liljenquist, 2006, study 2",
      side = 3, line = 2.5, las = 1, cex = 2, adj=-0.3)

dev.off()


