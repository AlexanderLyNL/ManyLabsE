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

# Gray data --------

# Original study: Lower bound of effect size found in the original study
deltaMin <- (5.29-3.86)/1.86



# Prospective frequentist analysis
freqDesign <- power.t.test(delta=deltaMin, power=0.8,
                           alternative="one.sided")

# Prospective e-value analysis
# designObj1 <- designSaviT(deltaMin=deltaMin, beta=0.2, seed=1,
#                           testType="twoSample", alternative="greater")
designObj1 <- designSaviT(deltaMin=deltaMin, beta=0.2, seed=2,
                                  testType="twoSample")


#   Data retrieved from the ManyLabs2 project, see Appendix below
#
# load("~/dropbox/projects/savitutorial/code/grayData.RData")

sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"


load(paste0(myWd, "grayData.RData"))


allSources <- unique(grayData$source)

# Result containers
#   General data set attributes
n1Vec <- n2Vec <- ratios  <- numeric(length(allSources))


#   sample sizes for p-value based inference
n1VecFreq <- n2VecFreq <- pValues <- numeric(length(allSources))

#   sample sizes for e-value based inference
#
n1VecE <- n2VecE <- firstTimes <- eValues <- numeric(length(allSources))
n1VecEFut <- n2VecEFut <- firstTimesFut <- eValuesFut <- numeric(length(allSources))

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

  tempResult <- saviTTestFut(x[1:n1EValue],
                             y[1:n2EValue],
                             designObj=designObj1, sequential=TRUE,
                             futility=TRUE,
                             esMinFutility=deltaMin,
                             nuMin=2)

  eValues[i] <- max(tempResult$eValueVec, na.rm=TRUE)

  # Used to fill up an e-value sequence if there is too little data
  nLast <- length(tempResult$eValueVec)
  nRemaining <- designObj1$nPlan[1] - nLast

  if (nRemaining > 0) {
    tempResult$eValueVec <- c(tempResult$eValueVec, rep(tempResult$eValueVec[nLast], nRemaining))


    tempResult$eValueVecFut <- c(unlist(tempResult$eValueVecFut), rep(unlist(tempResult$eValueVecFut)[nLast], nRemaining))
  }

  allEValueVecs[, i] <- tempResult$eValueVec
  allEValueVecsFut[, i] <- tempResult$eValueVecFut

  firstTimes[i] <- min(which(tempResult$eValueVec >= 20))
  firstTimesFut[i] <- min(which(tempResult$eValueVecFut <= betaFutility))
}

firstTimes
firstTimesFut



min(allEValueVecsFut)

firstTimes

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
                breaks=1:designObj1$nPlanBatch[1])

plot(fptHist)

y <- fptHist[["density"]]
nB <- length(fptHist$breaks)
yRange <- range(y, 0)

alpha <- 0.05
ylim <- c(-1*log(20/(2*alpha)), 2.75*log(1/alpha))

someConstant <- (ylim[2]+log(alpha))/yRange[2]
textHeightQuant <- (ylim[2]+log(alpha))+log(1/alpha)

xlim <- c(0, designObj1$nPlan[1])

notStoppedTable <- table(stoppedTimes[which(is.infinite(firstTimes))])


notStoppedN <- as.integer(names(notStoppedTable))
notStoppedBottom <- c(0, 2/61, 0, 1/61, 0)
notStoppedTop <- c(1/61, 3/61, 1/61, 2/61, 6/61)

myName <- "grayExample"
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



axis(side = 1, at=c(0, 10, 20, 30, 40), cex.axis=myCexAxis)

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

overIndexes <- which(is.finite(firstTimes))
underIndexes <- which(is.infinite(firstTimes))



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
        lwd=5, col=overColour)

  points(n1Temp,
         log(1/alpha),
         pch=15, col=overColourBorder)
}

mtext("e-value analyses of ManyLabs2 replications of Gray and Wegner, 2009, study 1a", side = 3,
      line = 2.5, las = 1, cex = 2, adj=-0.3)
dev.off()
