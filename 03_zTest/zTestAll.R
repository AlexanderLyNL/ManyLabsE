# Setup paths ----
library(devtools)
library(plyr)
library(rio)
library(tidyverse)

sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/03_zTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# Run all z-test

res1 <- manyLabsMetaScenarios(1, analysisType="zTest", nSim=1e3L)
res3 <- manyLabsMetaScenarios(3, analysisType="zTest", designObjList=res1$designObjList, nSim=1e3L)
res2 <- manyLabsMetaScenarios(2, analysisType="zTest", designObjList=res1$designObjList, nSim=1e3L)

res1DeltaMinFactor <- manyLabsMetaScenarios(1, analysisType="zTest", nSim=1e3L, deltaMinFactor=1)
res3DeltaMinFactor <- manyLabsMetaScenarios(3, analysisType="zTest",
                                            designObjList=res1DeltaMinFactor$designObjList, nSim=1e3L, deltaMinFactor=1)
res2DeltaMinFactor <- manyLabsMetaScenarios(2, analysisType="zTest",
                                            designObjList=res1DeltaMinFactor$designObjList, nSim=1e3L, deltaMinFactor=1)


res1TwoSided <- manyLabsMetaScenarios(1, analysisType="zTest", nSim=1e3L)
res3TwoSided <- manyLabsMetaScenarios(3, analysisType="zTest",
                                      designObjList=res1TwoSided$designObjList, nSim=1e3L)
res2TwoSided <- manyLabsMetaScenarios(2, analysisType="zTest",
                                      designObjList=res1TwoSided$designObjList, nSim=1e3L)

# save(res1, res2, res3,
#      res1DeltaMinFactor, res2DeltaMinFactor, res3DeltaMinFactor,
#      res1TwoSided, res2TwoSided, res3TwoSided,
#      file="zTestAll.RData")

res1z <- res1
res2z <- res2
res3z <- res3
res1DeltaMinFactorz <- res1DeltaMinFactor
res2DeltaMinFactorz <- res2DeltaMinFactor
res3DeltaMinFactorz <- res3DeltaMinFactor
res1TwoSidedz <- res1TwoSided
res2TwoSidedz <- res2TwoSided
res3TwoSidedz <- res3TwoSided


res1DeltaMinFactor$resultTable
