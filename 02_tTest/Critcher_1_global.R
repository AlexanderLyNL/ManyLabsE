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

study.description      <- 'Incidental Anchors (Critcher & Gilovich, 2008)'
analysis.unique.id     <- 36
analysis.name          <- 'Critcher.1'
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



# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df %>%
  dplyr::select(2,7,297,302,521,522,523,524,525,526,527,528,529,530,531,532,535,536,537) %>%
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

# To see the function code type:varfun.Critcher.1, or lookup in manylabRs_SOURCE.R
ML2.var[[g]] <- varfun.Critcher.1(ML2.sr[[g]])


        # Check equal variance assumption
        if(!is.na(testVarEqual)){
          if(testVarEqual){
            logtxt <- paste(analysis.unique.id,ML2.key$study.analysis,'-', runGroups[g])
            ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]],ML2.in$study.vars.labels, ML2.key, group = logtxt) # don't pass the cleanData frame
          }}

        # Run the analysis according to ML2.key: 'stat.test'
        stat.params <<- ML2.in$stat.params


stat.test   <- try.CATCH(with(ML2.var[[g]],t.test(x = P97, y = P17, conf.level=stat.params$conf.level, var.equal = stat.params$var.equal, alternative = stat.params$alternative)))


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

ML2.var[[g]] <- varfun.Critcher.1(ML2.sr[[g]])

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
print("FACTOR CHANGE SIGN")
print("https://docs.google.com/document/d/1b7MTOAiB7NPWlYBnwkhNTsj9i-upDMODgQ9t_djjSz8/edit?tab=t.0")
print("Critcher stats are duplicated across columns for the pencil-and-paper sites (variables crit2.1 and crit1.1, which should only have a value for one or the other). For example, site ‘bc’:")


dat <- addSources(ML2.var, ML2.df)
dat$factor <- ordered(dat$factor,
                      sort(unique(dat$factor), decreasing=TRUE))

freqRes2 <- t.test(variable ~ factor, data = dat, var.equal = stat.params$var.equal)
freqRes2$statistic
freqRes$statistic

# save(dat, stat.params, file="critcher.RData")

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
deltaMin <- 0.3
varEqual <- stat.params$var.equal
power <- 0.8
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative

designObj <- designSaviT(alpha=alpha, power=power,
                         deltaMin=deltaMin, futility=TRUE,
                         betaFutility=betaFutility,
                         varEqual=varEqual, testType="twoSample",
                         alternative=alternative)

# Scenario 1 ----
res1 <- scenario1T(dat=dat, allSources=allSources, designObj=designObj,
                   nuMin=3, alpha=alpha, betaFutility=betaFutility,
                   nSim=1e3, alternative=alternative)

mean(res1$eValues >= 1/alpha)
mean(res1$eValuesFut <= betaFutility)

res1$nStudiesAlternativeWorstCase
res1$nStudiesFutilityWorstCase

res1$nSamplesAlternativeWorstCase
res1$nSamplesFutilityWorstCase

mean(res1$stopDecision==1)
mean(res1$stopDecision==-1)

mean(res1$nStudies)

mean(res1$logMetaE)
sd(res1$logMetaE)

mean(res1$logMetaEFut)
sd(res1$logMetaEFut)

mean(res1$totalStoppingTimes)
sd(res1$totalStoppingTimes)

# Scenario 2-----
res2 <- scenario2T(dat, allSources, designObj=designObj, seed=1, nSim=1e3)

logMetaE<- rowSums(log(res2$eValues))
mean(logMetaE)
sd(logMetaE)

logMetaEFut <- rowSums(log(res2$eValuesFut))
mean(logMetaEFut)
sd(logMetaEFut)

mean(res2$alternativeProportion)
sd(res2$alternativeProportion)

mean(res2$futilityProportion)
sd(res2$futilityProportion)


mean(res2$totalStoppingTimes)
sd(res2$totalStoppingTimes)

#Scenario 3 ------

res3 <- scenario3T(dat=dat, allSources=allSources, designObj=designObj,
                   alpha=alpha, betaFutility=betaFutility,
                   nuMin=nuMin, nSim=1e3L)

mean(res3$logMetaE)
sd(res3$logMetaE)

mean(res3$logMetaEFut)
sd(res3$logMetaEFut)

mean(res3$alternativeProportion)
sd(res3$alternativeProportion)

mean(res3$futilityProportion)
sd(res3$futilityProportion)

mean(res3$totalStoppingTimes)
sd(res3$totalStoppingTimes)

# save(res1, res2, res3, file="critcher1Result.RData")
