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

# debugonce(saviTTest)

# bep1$designObjList$knobe$bootObjN1Plan$data
bep1$resultTable

source(file.path(project.root, "00_utils", "helpers.R"))

bep1 <- manyLabsMetaScenarios(
  1, designObjList=bep1$designObjList,
  nSim=1e3L)

bep3 <- manyLabsMetaScenarios(
  3, designObjList=bep1$designObjList,
  nSim=1e3L)

bep2 <- manyLabsMetaScenarios(
  2, designObjList=bep1$designObjList,
  nSim=1e3L)

res1 <- bep1
res2 <- bep2
res3 <- bep3

res1$resultTable
res2$resultTable
res3$resultTable
# save(res1, res2, res3, file="tTestAll.RData")
