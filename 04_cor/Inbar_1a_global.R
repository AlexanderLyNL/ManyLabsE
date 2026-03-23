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

# ANALYSIS INFO ----

study.description      <- 'Disgust & Homophobia (Inbar et al., 2009)'
analysis.unique.id     <- 26
analysis.name          <- 'Inbar.1a'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset.type <- "all"
saveAll <- FALSE


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



# PREPARE DATA & OUTPUT ----
# Add a unique ID
ML2.df$uID <- seq(1, nrow(ML2.df))

# Get info to create a dataset for the current study
ML2.in <- get.info(ML2.key, colnames(ML2.df), subset.type)

# Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
# Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
ML2.id <- get.chain(ML2.in)

ML2.id$df

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2,7,165,174,382,383,389,390,391,521,522,523,524,525,526,527,528,529,530,531,532,535,536,537) %>%
  dplyr::filter(is.character(source))


# Decide which analyses to run on which groups
toRun  <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

ML2.df$study.order <- NA
stmp <- strsplit(ML2.df$StudyOrderN,"[|]")

# Correct differences in study names
Stud <- ML2.key$study.name
if(Stud%in%"Tversky"){Stud <- "Tversky.Gati"}
if(Stud%in%"Rottenstreich"){Stud <- "Rottenstrich"}
if(Stud%in%"Ross"&(ML2.key['study.slate'] == 1)){Stud <- "Ross.Slate1"}
if(Stud%in%"Ross"&(ML2.key['study.slate'] == 2)){Stud <- "Ross.Slate2"}
if(Stud%in%"vanLange"){Stud <- "VanLange"}
if(Stud%in%"Giessner"){Stud <- "Geissner"}

ML2.df$study.order <- plyr::laply(seq_along(stmp), function(o){which(grepl(Stud,stmp[[o]]))%00%NA})

ML2.sr       <- list()
ML2.var      <- list()
outputSource <- list()
dataSource   <- list()
raw.df       <- list()
clean.df     <- list()
cleanData    <- list()
testVarEqual <- ML2.in$stat.params$var.equal

# Loop over sites in runGroups within a study
if(analysis.type==1){
  runGroups <- "all"
} else {
  runGroups <- sort(na.exclude(unique(ML2.df[[toRun$ugroup]])))
}

disp(paste(analysis.unique.id, ML2.key$study.analysis,"- START"), header = toupper(ML2.key$study.analysis), footer = FALSE)
cat("\n")


# START GROUPS ----
# HERE -------


# Include only datasets that have N >= Nmin.raw & n.group >= Nmin.cond
listIT     <- FALSE
nMin1      <- FALSE
nMin2      <- FALSE
compN <- compN1 <- compN2 <- 0

gID <- rep(TRUE, nrow(ML2.df))

# Check nMin
if(sum(gID, na.rm=TRUE) >= Nmin.raw){
  nMin1 <- TRUE
  # Get a list containing the data frames to be used in the analysis
  ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
}

# Double-check nMin
compN  <- ML2.sr[[g]]$N
compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
if(any(compN >= Nmin.raw)&(all(compN1>=Nmin.cond, compN2>=Nmin.cond))){nMin2 <- TRUE}


# START ANALYSIS ----------------------------------------



# Freq test ------
ML2.var[[g]] <- varfun.Inbar.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


freqTest <- try.CATCH(with(ML2.var[[g]],
                           cor_test_fisherZ(r1=r1,r2=r2,n1=N[1],n2=N[2], conf.level=stat.params$conf.level, alternative = stat.params$alternative)))

dat <- ML2.var[[g]]$cleanDataFilter



studySummary <- dat %>%
  group_by(factor) %>%
  summarise(
    n = n(),
    mean1 = mean(variable1, na.rm = TRUE),
    sd1 = sd(variable1, na.rm=TRUE),
    mean2 = mean(variable2, na.rm = TRUE),
    sd2 = sd(variable2, na.rm=TRUE),
    cor = cor(variable1, variable2),
    pValue = cor.test(variable1, variable2)$p.value
  )

studySummary
freqTest$value$p.value


diff(-atanh(freqTest$value$estimate))
