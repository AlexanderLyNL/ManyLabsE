# Setup paths ----
library(devtools)
library(plyr)
library(rio)
library(tidyverse)

sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# Run all t-test

res1 <- manyLabsMetaScenarios(1, analysisType="tTest", nSim=1e3L)
# res3 <- manyLabsMetaScenarios(3, analysisType="tTest", designObjList=res1$designObjList, nSim=1e3L)
# res2 <- manyLabsMetaScenarios(2, analysisType="tTest", designObjList=res1$designObjList, nSim=1e3L)

# res3 <- manyLabsMetaScenarios(3, analysisType="tTest", designObjList=res1$designObjList, nSim=1e3L)

# res3 <- manyLabsMetaScenarios(3, analysisType="tTest",
#                               designObjList=res1$designObjList,
#                               nSim=1e3L, wantPaths=25)


res3a <- manyLabsMetaScenarios(3,
                               analysisType="tTest",
                               designObjList=res1$designObjList, nSim=10)
# res3PathsOld <- res3Paths
# res3Paths <- res3



res3Paths$individualResultList$knobe$paths[[1]]$logSamplePaths

object.size(res3)

# k <- 14
# max(which(res3$individualResultList$knobe$paths[[1]]$logSamplePaths[k, ]!=0))

# load(file="~/dropbox/projects/manyLabsE/results/tTestAll.RData")



# load(file="tTestAll.RData")

load(file="~/dropbox/projects/manylabse/results/tTestAll.RData")




hist(res1DeltaMinFactor$individualResultList$risen$totalStoppingTimes)
hist(res2DeltaMinFactor$individualResultList$risen$totalStoppingTimes)
hist(res3DeltaMinFactor$individualResultList$risen$totalStoppingTimes)
quantile(res3DeltaMinFactor$individualResultList$risen$totalStoppingTimes, 0.975)

mean(res3DeltaMinFactor$individualResultList$risen$totalStoppingTimes)+sd(res3DeltaMinFactor$individualResultList$risen$totalStoppingTimes)


res1DeltaMinFactor$resultTable
res2DeltaMinFactor$resultTable
res3DeltaMinFactor$resultTable

res1$resultTableFull

nBand <- quantile(res1DeltaMinFactor$individualResultList$knobe$totalStoppingTimes, c(0.05, 0.95))
nBand

names(res1DeltaMinFactor$individualResultList$knobe)



studyNames

manyLabsMetaScenarios

# # myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"
#
# res1DeltaMinFactor <- manyLabsMetaScenarios(1, analysisType="tTest", nSim=1e3L, deltaMinFactor=1)
# res3DeltaMinFactor <- manyLabsMetaScenarios(3, analysisType="tTest",
#                                             designObjList=res1DeltaMinFactor$designObjList, nSim=1e3L, deltaMinFactor=1)
# res2DeltaMinFactor <- manyLabsMetaScenarios(2, analysisType="tTest",
#                                             designObjList=res1DeltaMinFactor$designObjList, nSim=1e3L, deltaMinFactor=1)
#
# res1TwoSided <- manyLabsMetaScenarios(1, analysisType="tTest", nSim=1e3L)
# res3TwoSided <- manyLabsMetaScenarios(3, analysisType="tTest",
#                                       designObjList=res1TwoSided$designObjList, nSim=1e3L)
# res2TwoSided <- manyLabsMetaScenarios(2, analysisType="tTest",
#                                       designObjList=res1TwoSided$designObjList, nSim=1e3L)

# save(res1, res2, res3,
#      res1DeltaMinFactor, res2DeltaMinFactor, res3DeltaMinFactor,
#      res1TwoSided, res2TwoSided, res3TwoSided,
#      file="tTestAll.RData")


# save(res3Paths, file="tTestRes3Paths.RData")

myHeightEPS <- 6.5
myWidthEPS <- 10

ltyGauss <- 4
someLty <- 6

myHeightEPS <- 6.5
myWidthEPS <- 10

pdfWidth <- 14
pdfHeight <- 7

myWidth <- 1070
myHeight <- 485

myCex <- 2
myLwd <- 7
myCexLab <- 1.5

underColour <- "#FFB90F86"
underColourBorder <- "#FFB90FCC"
overColourBorder <- "#1F78B4E6"
overColour <- "#A6CEE380"

bord <- col <- res3$resultTable$logMetaE

col[res3$resultTable$logMetaE > 0] <- overColour
bord[res3$resultTable$logMetaE > 0] <- overColourBorder

col[res3$resultTable$logMetaE < 0] <- underColour
bord[res3$resultTable$logMetaE < 0] <- underColourBorder


myName <- paste0("aap")
pdf(paste0(myName, ".pdf"), width=pdfWidth, height=pdfHeight)
graphics::par(cex.main=1.5, mar=c(5, 7, 4, 4)+0.1, mgp=c(3.5, 1, 0), cex.lab=1.5,
              font.lab=2, cex.axis=1.3, bty="n", las=1)
barplot(height=res3$resultTable$`% of`,
        names=rownames(res3$resultTable),
        las=1, horiz=TRUE, xlim=c(0, 50), col=col,
        border=bord, xlab="% of realised sample sizes")
dev.off()

bord

# barplot(height=res3$resultTable$`% of`,
#         names=rownames(res3$resultTable),
#         las=1, horiz=FALSE, ylim=c(0, 50))


res3$resultTable

res3$individualResultList$kay$logEValues[1, ]
res3$individualResultList$kay$logEValuesFut[1, ]
