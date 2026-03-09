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

study.description      <- 'Incidental Disfluency (Alter et al., 2007)'
analysis.unique.id     <- 8
analysis.name          <- 'Alter.1'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset                 <- 'all'
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
  dplyr::select(2,7,314,315,316,317,318,319,327,328,329,330,331,332,521,522,523,524,525,526,527,528,529,530,531,532,535,536,537) %>%
  dplyr::filter(is.character(source))

# Decide which analyses to run on which groups
toRun  <- decide.analysis(ML2.key, analysis.unique.id, analysis.type, doAll = TRUE)

if(length(toRun$studiess)>0){

  if(NROW(ML2.df)>0){

    # Create a variable indicating the study order for each case
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

    for(g in seq_along(runGroups)){

      # Include only datasets that have N >= Nmin.raw & n.group >= Nmin.cond
      listIT     <- FALSE
      nMin1      <- FALSE
      nMin2      <- FALSE
      compN <- compN1 <- compN2 <- 0

      if(analysis.type<4){
        if(runGroups[g]=="all"){
          gID <- rep(TRUE, nrow(ML2.df))
        } else {
          gID <- ML2.df$source%in%runGroups[g]
        }
      } else {
        gID <-  ML2.df$study.order%in%runGroups[g]
      }

      # Check nMin
      if(sum(gID, na.rm=TRUE) >= Nmin.raw){
        nMin1 <- TRUE
        # Get a list containing the data frames to be used in the analysis
        ML2.sr[[g]] <- get.sourceData(ML2.id, ML2.df[gID, ], ML2.in)
      }

      # Double-check nMin
      if(nMin1){
        compN  <- ML2.sr[[g]]$N
        compN1 <- sum(ML2.sr[[g]]$RawDataFilter[[1]]$Included, na.rm = TRUE)
        compN2 <- sum(ML2.sr[[g]]$RawDataFilter[[2]]$Included, na.rm = TRUE)
        if(any(compN >= Nmin.raw)&(all(compN1>=Nmin.cond, compN2>=Nmin.cond))){nMin2 <- TRUE}
      }

      # START ANALYSIS ----------------------------------------

      if(all(nMin1,nMin2)){

# To see the function code type:varfun.Alter.1, or lookup in manylabRs_SOURCE.R
ML2.var[[g]] <- varfun.Alter.1(ML2.sr[[g]])


        # Check equal variance assumption
        if(!is.na(testVarEqual)){
          if(testVarEqual){
            logtxt <- paste(analysis.unique.id,ML2.key$study.analysis,'-', runGroups[g])
            ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]],ML2.in$study.vars.labels, ML2.key, group = logtxt) # don't pass the cleanData frame
          }}

        # Run the analysis according to ML2.key: 'stat.test'
        stat.params <<- ML2.in$stat.params


stat.test   <- try.CATCH(with(ML2.var[[g]],t.test(x = DisFluent, y = Fluent, conf.level=stat.params$conf.level, var.equal = stat.params$var.equal, alternative = stat.params$alternative)))


        # Check for errors and warnings
        if(all(is.null(stat.test$warning), grepl("simpleWarning",stat.test$warning),
               !grepl("Error", stat.test$value[[1]]),
               !grepl("message", names(unlist(stat.test))[1]))){
          stat.test  <- stat.test$value
          ConsoleOut <- paste(capture.output(print(stat.test)),collapse="\n")
          listIT     <- TRUE
        }

        if(listIT){
          describe <- get.descriptives(stat.test = stat.test,
                                       vars      = ML2.var[[g]],
                                       keytable  = ML2.key)

          if(any(describe$descr.raw$n<Nmin.cond)){
            listIT<- FALSE
            nMin2 <- FALSE
          }
          rm(describe)
        }

        # START RECORD DATA ----

        if(listIT){

          describe <- get.descriptives(stat.test = stat.test,
                                       vars      = ML2.var[[g]],
                                       keytable  = ML2.key)

          var.lor <- ifelse(grepl("OR",describe$test$estype),
                            sum(1/(table(ML2.var[[g]]$Condition,ML2.var[[g]]$Response)), na.rm = TRUE),
                            NA)

          ESCI  <-   generateOutput(describe        = describe,
                                    var.lor         = var.lor,
                                    runningGroup    = runGroups[g],
                                    runningAnalysis = paste(analysis.unique.id,ML2.key$study.analysis),
                                    stat.params = stat.params)

          # Raw and clean datasets
          if(length(ML2.sr[[g]]$RawDataFilter)>1){
            case.include <- ML2.sr[[g]]$RawDataFilter[[1]]$Included|ML2.sr[[g]]$RawDataFilter[[2]]$Included
            df1 <- ML2.sr[[g]]$RawDataFilter[[1]][ ,-which(colnames( ML2.sr[[g]]$RawDataFilter[[1]])=="Included")]
            raw.df[[g]] <-  cbind.data.frame(df1, analysis.type = c("Global","Primary","Secondary","Order")[analysis.type],subset=subset,case.include = case.include)
          } else {
            case.include <- ML2.sr[[g]]$RawDataFilter[[1]]$Included
            df1 <- ML2.sr[[g]]$RawDataFilter[[1]][ ,-which(colnames( ML2.sr[[g]]$RawDataFilter[[1]])=="Included")]
            raw.df[[g]] <-  cbind.data.frame(df1, analysis.type = c("Global","Primary","Secondary","Order")[analysis.type],subset=subset,cases.include = case.include)
          }





SourceInfo <- raw.df[[g]] %>% dplyr::filter(case.include) %>%
dplyr::summarise(
  N.sources.global    = length(unique(Source.Global)),
  N.sources.primary   = length(unique(Source.Primary)),
  N.sources.secondary = length(unique(Source.Secondary)),
  N.countries         = length(unique(Country)),
  N.locations         = length(unique(Location)),
  N.languages         = length(unique(Language)),
  Pct.WEIRD           = mean(Weird, na.rm=TRUE)*100,
  Tbl.Execution       = paste0(capture.output(table(Execution)),collapse='\n'),
  Tbl.subjectpool     = paste0(capture.output(table(SubjectPool)),collapse='\n'),
  Tbl.setting       = paste0(capture.output(table(Setting)),collapse='\n'),
  Tbl.Tablet        = paste0(capture.output(table(Tablet)),collapse='\n'),
  Tbl.Pencil        = paste0(capture.output(table(Pencil)),collapse='\n'),
  N.studyorders1    = length(unique(StudyOrderN)),
  N.IDiffOrderN     = length(unique(IDiffOrderN)),
  N.uIDs            = length(unique(uID)),
  N.studyorders2    = length(unique(study.order)),
  Tbl.analysistype  = paste0(capture.output(table(analysis.type)),collapse='\n'),
  Tbl.subset        = paste0(capture.output(table(subset)),collapse='\n'),
  N.cases.included  = sum(case.include, na.rm=TRUE),
  N.cases.excluded  = sum(case.include==FALSE,na.rm=TRUE))





          rownames(SourceInfo) <- NULL

          test  <- describe$test
          descr <- describe$descr.raw
          outputSource[[g]] <- get.output(key      = ML2.key,
                                          vars     = ML2.var[[g]],
                                          descr    = descr,
                                          group    = runGroups[g],
                                          analysis = c("Global","Primary","Secondary","Order")[analysis.type],
                                          varEqual = stat.params$var.equal,
                                          test     = test,
                                          ESCI     = ESCI,
                                          test.ConsoleOutput = ConsoleOut,
                                          SourceInfo = SourceInfo,
                                          stat.test = stat.test)

          # Data list for output to spreadsheet
          dataSource[[g]] <- list(
            study.id      = ML2.key$study.id,
            study.slate   = ML2.key$study.slate,
            study.name    = ML2.key$study.name,
            study.source  = runGroups[g],
            analysis.type = c("Global","Primary","Secondary","Order")[analysis.type],
            analysis.name = ML2.key$study.analysis,
            subset        = subset,
            stat.info     = ML2.in,
            stat.data.cleanchain = ML2.id,
            stat.data.raw       = raw.df[[g]],
            stat.data.cleaned   = ML2.sr[[g]][1:length(ML2.sr[[g]])-1],
            stat.data.analysed  = ML2.var[[g]][1:length(ML2.var[[g]])-1],
            stat.test = stat.test)

          suppressMessages(clean.df[[g]] <- plyr::ldply(dataSource[[g]]$stat.data.analysed,reshape2::melt))
          colnames(clean.df[[g]])[colnames(clean.df[[g]])==".id"] <- "Condition"


          cleanData[[g]] <- ML2.var[[g]]$cleanDataFilter


          rm(stat.params)

        } else { # LISTIT
          cat("\nListIT = FALSE\n")
          if(!is.null(stat.test$value)){

            if(grepl("observations",as.character(stat.test$value))){
              disp(paste(analysis.unique.id, ML2.key$study.analysis,'-',
                         runGroups[g],'>> Not enough observations'),
                   header = FALSE, footer = FALSE)}
          } else {
            disp(paste(analysis.unique.id,ML2.key$study.analysis,'-', runGroups[g],'>> stat.test failed:'),
                 header = FALSE, footer = FALSE)
            # disp(paste('value: ',stat.test$value),
            #      header = FALSE, footer = FALSE)
            disp(paste('warning:',stat.test$warning),
                 header = FALSE, footer = FALSE)
          }
          ConsoleOut <- paste(gsub("[[:punct:]]", "", stat.test$warning, perl = TRUE), collapse="\n")
          NN <- lengths(ML2.var[[g]])
          NN <- NN[!names(NN)=="N"]
          N  <- rep(ML2.var[[g]]$N,length.out = length(NN))

          opt <- OutputTemplate()

          outputSource[[g]] <- get.output(key      = ML2.key,
                                          vars     = ML2.var[[g]],
                                          descr    = opt,
                                          group    = runGroups[g],
                                          analysis = c("Global","Primary","Secondary","Order")[analysis.type],
                                          varEqual = stat.params$var.equal,
                                          test     = test,
                                          ESCI     = ESCI,
                                          test.ConsoleOutput = ConsoleOut,
                                          SourceInfo = NA,
                                          stat.test = stat.test)

          cleanData[[g]] <- ML2.var[[g]]$cleanDataFilter

          # Data list for output to spreadsheet
          dataSource[[g]] <- NA


          rm(stat.params)

        } #Listit = FALSE
      } # all nMin 1,2


      # HANDLE ERRORS ----
      if(!nMin1){
        disp(paste0(analysis.unique.id,' ',ML2.key$study.analysis,' - ',
                    runGroups[g],' not included in results >> Cases in source file (',
                    sum(gID, na.rm = TRUE),') < Nmin.raw (',Nmin.raw,')'),
             header = FALSE, footer = FALSE)
      } # Check nMin 1}
      if(!nMin2){
        disp(paste0(analysis.unique.id,' ',ML2.key$study.analysis,' - ',
                    runGroups[g],' not included in results >> Valid cases after varfun (n',
                    c(1,2)[compN < Nmin.cond],"=", compN[compN < Nmin.cond],') < Nmin.cond (',Nmin.cond,')'),
             header = FALSE, footer = FALSE)
      } # Double-check nMin

    } # iterate groups

    disp(paste(analysis.unique.id, ML2.key$study.analysis,"- COMPLETED"), header = FALSE)

    ML2.output  <- plyr::ldply(outputSource)
    ML2.rawdata <- plyr::ldply(raw.df)
    ML2.cleandata <- plyr::ldply(cleanData)

    if(saveAll){
      if(!overWrite){
        basename <- paste0(gsub("[.]","_",analysis.name),"_",analysis.type.name,"_",subset, as.Date(now()))
      } else {
        basename <- paste0(gsub("[.]","_",analysis.name),"_",analysis.type.name,"_",subset)
      }
      if(outdir$Results!=""){
        rio::export(ML2.output, file = file.path(normalizePath(outdir$Results),paste0(basename,"_RESULTS.csv")))
        rio::export(ML2.output, file = file.path(normalizePath(outdir$Results),paste0(basename,"_RESULTS.xlsx")))
      }
      if(outdir$Data!=""){
        if(NCOL(ML2.rawdata)>0){
          rio::export(ML2.rawdata, file = file.path(normalizePath(outdir$Data),paste0(basename,"_RAW_CASE.csv")))
          rio::export(ML2.rawdata, file = file.path(normalizePath(outdir$Data),paste0(basename,"_RAW_CASE.xlsx")))
        }
        if(NCOL(ML2.cleandata)>0){
          if(all(!is.null(ML2.cleandata$uID),!is.null(ML2.rawdata$uID))){
            ML2.cleandata <- dplyr::left_join(ML2.cleandata,ML2.rawdata,by="uID")
          }
          rio::export(ML2.cleandata, file = file.path(normalizePath(outdir$Data),paste0(basename,"_CLEAN_CASE.csv")))
          rio::export(ML2.cleandata, file = file.path(normalizePath(outdir$Data),paste0(basename,"_CLEAN_CASE.xlsx")))
        }
        save(dataSource, file = file.path(normalizePath(outdir$Results),paste0(basename,"_R_OBJECTS.RData")))
      }
    }

    #rm(ML2.in, ML2.var, ML2.id, ML2.df, ML2.sr, outputSource, dataSource, raw.df, clean.df, descr, SourceInfo, nMin1, nMin2, listIT)

  } else { # if nrow > 0

    disp(paste(analysis.unique.id, ML2.key$study.analysis,"- SKIPPED"), header = FALSE)

    ML2.output  <- NULL
    ML2.rawdata <- NULL

    #rm(ML2.in, ML2.var, ML2.id, ML2.df, ML2.sr)
  }
}

# Freq test ------

ML2.var[[g]] <- varfun.Alter.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params


freqRes <- t.test(variable ~ factor, data = ML2.var[[g]]$cleanDataFilter, var.equal = stat.params$var.equal)

dat <- ML2.var[[g]]$cleanDataFilter


studySummary <- dat %>%
  group_by(factor) %>%
  summarise(
    n = n(),
    mean = mean(variable, na.rm = TRUE),
    sd = sd(variable, na.rm=TRUE)
  )

sum(studySummary$n)
studySummary$mean
studySummary$sd


freqRes$statistic
freqRes$p.value
freqRes$statistic*sqrt(sum(studySummary$n)/prod(studySummary$n))

# Alexander ----
dat <- addSources(ML2.var, ML2.df)
# save(dat, stat.params, file="alter.RData")

dat <- checkUniqueIds(dat)
tempRes <- removeOneConditionSources(dat)

allSources <- tempRes$allSources
sampleSize <- tempRes$sampleSize

dat <- dat[dat$source %in% allSources, ]

if (stat.params$alternative=="two.sided")
  stat.params$alternative <- "twoSided"


# Here -------
alpha <- 0.05
betaFutility <- alpha
deltaMin <- 0.63

designObj <- designSaviT(alpha=alpha,
  deltaMin=deltaMin, futility=TRUE, power=0.8,
  varEqual=stat.params$var.equal, testType="twoSample",
  alternative=stat.params$alternative, wantSampling=FALSE)

maxNX <- max(sampleSize[, 1])
maxNY <- max(sampleSize[, 2])

nMax <- max(rowSums(sampleSize))

factorLevels <- if (is.ordered(dat$factor)) levels(dat$factor) else unique(dat$factor)


# Result containers ----
#   General data set attributes
n1Vec <- n2Vec <- ratios  <- numeric(length(allSources))


#   sample sizes for p-value based inference
n1VecFreq <- n2VecFreq <- pValues <- numeric(length(allSources))

#   sample sizes for e-value based inference
#
n1VecE <- n2VecE <- firstTimes <- eValues <- numeric(length(allSources))
n1VecEFut <- n2VecEFut <- firstTimesFut <- eValuesFut <- numeric(length(allSources))

allEValueVecs <- matrix(nrow=nMax,
                        ncol=length(allSources))

allEValueVecsFut <- allEValueVecs

# Analyse data for each source
# loop start -----
for (i in 1:length(allSources)) {
  someDat <- dat[dat$source==allSources[i], ]

  # nParticipants <- length(someDat$uID)
  # someOrder <- sample(someDat$uID, nParticipants)

  ## Data -----
  x <- someDat[which(someDat$factor==factorLevels[1]), ]$variable
  y <- someDat[which(someDat$factor==factorLevels[2]), ]$variable

  # Remove non-available entries
  x <- x[!is.na(x)]
  n1 <- length(x)

  y <- y[!is.na(y)]
  n2 <- length(y)

  # Store valid sample size characteristics
  # n1VecFreq[i] <- n1Freq <- n1Vec[i] <- n1
  # n2VecFreq[i] <- n2Freq <- n2Vec[i] <- n2

  tempResult <- t.test(x[1:n1], y[1:n2],
                       var.equal=stat.params$var.equal)
  pValues[i] <- tempResult$p.value

  ## e-value ----
  # n1VecE[i] <- n1EValue <- sampleSize[i, 1]
  # n2VecE[i] <- n2EValue <- sampleSize[i, 2]

  ratios[i] <- n2/n1

  # debugonce(saviTTest)
  nParticipants <- n1+n2

  set.seed(seed)
  for (k in 1:nSim) {
    tempRes <- twoSampleTTestRandomOrder(
      "x"=x, "y"=y, "n1"=n1, "n2"=n2,
      "designObj"=designObj, "nuMin"=nuMin,
      "alpha"=alpha, "betaFutility"=betaFutility
    )

    tempRes$nStop
    tempRes$eValue
    tempRes$eValueFut
  }
}
# loop end ----

sum(is.finite(firstTimes))
sum(is.finite(firstTimesFut))

for (i in 1:61) {
  if (sum(is.na(allEValueVecs[, i]))>0)
    print(i)
}

allEValueVecs[, 8]

firstTimes
firstTimesFut

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

eMeta <- exp(cumsum(log(allEValueVecs[n1End, ])))
eFutMeta <- exp(cumsum(log(allEValueVecsFut[n1End, ])))

plot(eMeta, type="l", log="y")
lines(eFutMeta, col="red")

eMetaAverage <- cumsum(allEValueVecs[n1End, ])/(1:length(allSources))
eFutMetaAverage <- cumsum(allEValueVecsFut[n1End, ])/(1:length(allSources))

plot(eMetaAverage, type="l", log="y")
lines(eFutMetaAverage, col="red")

which(eMeta > 1/alpha)
which(eMetaAverage > 1/alpha)

which(eFutMeta <= betaFutility)
which(eFutMetaAverage <= betaFutility)


# Scenario 2 ----
eMeta <- exp(rowSums(log(allEValueVecs)))
eFutMeta <- exp(rowSums(log(allEValueVecsFut)))

plot(eMeta, type="l", log="y")
lines(eFutMeta, col="red")

eMetaAverage <- rowMeans(allEValueVecs)
eFutMetaAverage <- rowMeans(allEValueVecsFut)

plot(eMetaAverage, type="l", log="y")
lines(eFutMetaAverage, col="red")

which(eMeta >= 1/alpha)
which(eMetaAverage >= 1/alpha)

which(eFutMeta <= betaFutility)
which(eFutMetaAverage <= betaFutility)

# Scenario 3 ---------
nTotal <- length(unique(dat$uID))


set.seed(1)
someOrder <- sample(unique(dat$uID), nTotal)

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

  if (someRow$factor==factorLevels[1]) {
    sourceDataTracker[[someSource]]$x <- x <- c(x, someRow$variable)
  } else if (someRow$factor==factorLevels[2]) {
    sourceDataTracker[[someSource]]$y <- y <- c(y, someRow$variable)
  }

  if (i >1) {
    someCheck <- checkXY(x, y)

    if (someCheck) {
      tempRes <- saviTTest(x, y, designObj=designObj,
                           sequential=FALSE)

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

plot(1:nTotal, eMeta, type="l", log="y")
lines(1:nTotal, eFutMeta, col="red")


eMetaAverage <- rowMeans(eMatrix)
eFutMetaAverage <- rowMeans(eFutMatrix)

plot(1:nTotal, eMetaAverage, type="l", log="y")
lines(1:nTotal, eFutMetaAverage, col="red")

which(eMeta >= 1/alpha)
which(eMetaAverage >= 1/alpha)

which(eFutMeta <= betaFutility)
which(eFutMetaAverage <= betaFutility)

