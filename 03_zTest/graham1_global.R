library(devtools)
library(plyr)
library(rio)
library(tidyverse)

library(safestats)

#remotes::install_github("AlexanderLyNL/safestats", ref ="futility88")


sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "manyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))


# ANALYSIS INFO ----


study.description      <- 'Moral Foundations (Graham et al., 2009)'
analysis.unique.id     <- 14
analysis.name          <- 'Graham.1'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset                 <- 'all'
onlineTables           <- TRUE
staticData             <- TRUE
saveAll                <- FALSE
overWrite              <- FALSE
#OSFdata.root           <- paste0(myWd, "ManyLabsE/", "OSFdata") #file.path('~','OSFdata')
analysis.root          <- file.path(OSFdata.root,study.description,analysis.name,'Global')
outdir                 <- list(Data = file.path(analysis.root,'Data'), Results = file.path(analysis.root,'Results'))


# This function will be used to change the raw dataset to a dataset ready for analysis

varfun.Graham.1

varfun.Graham.1 <- function(vars){

  uID <- vars$Binding$uID

  binding_mean <- vars$Binding %>%
    dplyr::select(-uID) %>%
    rowMeans(na.rm = TRUE)

  cleanDataFilter <- data.frame(
    uID = uID,
    variable1 = vars$Politics$politics,
    variable2 = binding_mean
  )

  return(list(
    Politics        = vars$Politics$politics,
    Binding      = binding_mean,
    N               = sum(complete.cases(binding_mean, vars$Politics$politics)),
    cleanDataFilter = cleanDataFilter
  ))
}

if(dplyr::between(analysis.type,2,3)){subset <- "all"}

# GET LOOKUP TABLES ----

if(onlineTables){
  # Get the Keytable with analysis information
  ML2.key <- get.GoogleSheet(data='ML2masteRkey')$df
  ML2.key <- ML2.key[!is.na(ML2.key$unique.id)&ML2.key$unique.id==analysis.unique.id,]

  # Get info about the sites
  SourceInfoTable    <- get.GoogleSheet(url = "https://docs.google.com/spreadsheets/d/1Qn_kVkVGwffBAmhAbpgrTjdxKLP1bb2chHjBMVyGl1s/pub?gid=1435507167&single=true&output=csv")$df
} else {
  # Get the Keytable with analysis information
  ML2.key <- rio::import(file.path(OSFdata.root,"!!KeyTables","ML2_KeyTable.csv"))
  ML2.key <- ML2.key[!is.na(ML2.key$unique.id)&ML2.key$unique.id==analysis.unique.id,]

  # Get info about the sites
  SourceInfoTable    <- rio::import(file.path(OSFdata.root,"!!KeyTables","ML2_SourceInfoTable.csv"))
}

# GET DATA ----

if(!staticData){
  # CANNOT TEST UNTIL OSF DATA ARE PUBLIC
  # Get the correct slate according to info in ML2.key['study.slate']
  if(ML2.key[study,'study.slate'] == 1){
    data <- osfr::download_files(id = 'cwjp3', path =  getwd())
  } else {
    data <- osfr::download_files(id = 'jg9hc', path =  getwd())
  }
  ML2.df <- rio::import(data)
  disp(paste("Downloaded data from OSF"), header = FALSE, footer = FALSE)
} else {
  # Get the correct slate according to info in ML2.key['study.slate']
  if(ML2.key$study.slate == 1){
    ML2.df <- rio::import(file.path(OSFdata.root,"!!RawData","ML2_Slate1.csv"))
  } else {
    ML2.df <- rio::import(file.path(OSFdata.root,"!!RawData","ML2_Slate2.csv"))
  }
}

# PREPARE DATA & OUTPUT ----

# Add a unique ID
ML2.df$uID = seq(1, nrow(ML2.df))

# Get info to create a dataset for the current study
# keytable <- ML2.key
ML2.in <- get.info(ML2.key, colnames(ML2.df), subset)

# Generate chain to select variables for the data frame and create a filter chain for the variables to use for analysis
# Info based on KeyTable information in study.vars, cases.include, site.include, params.NA
ML2.id <- get.chain(ML2.in)

# Apply the df chain to select relevant subset of variables

ML2.df <- ML2.df  %>% dplyr::select(1,6,135,136,137,138,139,140,141,142,143,430,520,521,522,523,524,525,526,527,528,529,530,531,534,535,536) %>% dplyr::filter(is.character(source))



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

        # To see the function code type:varfun.Graham.1, or lookup in manylabRs_SOURCE.R
        ML2.var[[g]] <- varfun.Graham.1(ML2.sr[[g]])


        # Check equal variance assumption
        if(!is.na(testVarEqual)){
          if(testVarEqual){
            logtxt <- paste(analysis.unique.id,ML2.key$study.analysis,'-', runGroups[g])
            ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]],ML2.in$study.vars.labels, ML2.key, group = logtxt) # don't pass the cleanData frame
          }}

        # Run the analysis according to ML2.key: 'stat.test'
        stat.params <<- ML2.in$stat.params


        stat.test   <- try.CATCH(with(ML2.var[[g]],cor_test_fisherZ(r1=matrix(c(Binding,Politics),ncol=2), r2=NULL, n1=N[1], n2=NULL, conf.level=stat.params$conf.level, alternative = stat.params$alternative)))


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
ML2.var[[g]] <- varfun.Graham.1(ML2.sr[[g]])

stat.params <<- ML2.in$stat.params

freqRes <- try.CATCH(with(ML2.var[[g]],cor_test_fisherZ(r1=matrix(c(Binding,Politics),ncol=2), r2=NULL, n1=N[1], n2=NULL, conf.level=stat.params$conf.level, alternative = stat.params$alternative)))

corrie <- freqRes$value$estimate
cohensD <- 2*corrie/sqrt(1-corrie^2)




# Alexander -------
dat <- addSources(ML2.var, ML2.df)

sum(is.na(dat$variable1))
sum(is.na(dat$variable2))

dat <- checkUniqueIds(dat)

# save(dat, stat.params, file="graham.RData")

# Here -------
alpha <- 0.05
betaFutility <- alpha
deltaMin <- 0.52
power <- 0.8
alternative <- if (stat.params$alternative=="two.sided") "twoSided" else stat.params$alternative

designObj <- designSaviZ(deltaMin, alpha=alpha, power=power,
                         futility=TRUE, alternative=alternative,
                         betaFutility=betaFutility)

allSources <- unique(dat$source)

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



zVec <- rVec <- nVec <- eValueVec <- eValueFutVec <- numeric(length(allSources))
mean1Vec <- mean2Vec <- sd1Vec <- sd2Vec <- zVec

for (i in seq_along(rVec)) {
  someDat <- dat[dat$source==allSources[i], ]

  someN <- nVec[i] <- dim(someDat)[1]
  someR <- rVec[i] <- cor(someDat$variable1, someDat$variable2)
  someZ <- zVec[i] <- atanh(someR)/sqrt(1/(someN-3))

  mean1Vec[i] <- mean(someDat$variable1)
  sd1Vec[i] <- sd(someDat$variable1)
  mean2Vec[i] <- mean(someDat$variable2)
  sd2Vec[i] <- sd(someDat$variable2)

  tempRes <- saviZTestStat(
    z=someZ, n1=someN,
    parameter=designObj$parameter,
    eType=designObj$eType, sigma=designObj$sigma)

  eValueVec[i] <- tempRes$eValue

  tempRes <- saviFutilityZStat(
    z=someZ, n1=someN,
    parameter=designObj$futilityResult$parameter,
    sigma=designObj$sigma)

  eValueFutVec[i] <- tempRes$eValue
}

which(eValueVec > 1/alpha)
which(eValueFutVec < betaFutility)


#
eMeta <- exp(cumsum(log(eValueVec)))
eFutMeta <- exp(cumsum(log(eValueFutVec)))


which(eMeta >= 1/alpha)
which(eFutMeta <= betaFutility)

plot(eMeta, log="y", type="l")
lines(eFutMeta, col="red")
abline(h=betaFutility)
abline(h=1/alpha, col="blue")




# Stopping time

get_stop_time <- function(eVec, futVec, alpha, betaFutility) {

  eMeta <- exp(cumsum(log(eVec)))
  eFutMeta <- exp(cumsum(log(futVec)))

  stop_reject <- which(eMeta >= 1/alpha)[1]
  stop_fut <- which(eFutMeta <= betaFutility)[1]

  stop_reject <- ifelse(is.na(stop_reject), Inf, stop_reject)
  stop_fut <- ifelse(is.na(stop_fut), Inf, stop_fut)

  list(
    reject = stop_reject,
    futility = stop_fut,
    stop = min(stop_reject, stop_fut)
  )
}

# Average stopping time (random permutations)

m <- 1000
n <- length(eValueVec)

stops <- numeric(m)

for(i in 1:m){

  perm <- sample(n)

  e_perm <- eValueVec[perm]
  fut_perm <- eValueFutVec[perm]

  res <- get_stop_time(e_perm, fut_perm, alpha, betaFutility)

  stops[i] <- res$stop
}

mean(stops[is.finite(stops)])
var(stops[is.finite(stops)])
# We stop after only one study!
# Times stopped for futility/Null:
type <- character(m)

for(i in 1:m){

  perm <- sample(n)

  res <- get_stop_time(eValueVec[perm], eValueFutVec[perm], alpha, betaFutility)

  if(res$reject < res$futility) type[i] <- "reject"
  else if(res$futility < res$reject) type[i] <- "futility"
  else type[i] <- "none"
}
mean(type == "futility")# Over 90% times we stopped for futility

e_worst <- sort(eValueVec)              # smallest first
fut_worst <- sort(eValueFutVec, decreasing = TRUE)

res_worst <- get_stop_time(e_worst, fut_worst, alpha, betaFutility)

res_worst




# Scenario 2

n_sources <- length(allSources)
n_max <- max(table(dat$source))

zMat <- rMat <- nMat <- eValueMat <- eValueFutMat <- matrix(NA, n_sources, n_max)

compute_evalue_mats <- function(dat, allSources, designObj) {

  n_sources <- length(allSources)
  n_max <- max(table(dat$source))

  zMat <- rMat <- nMat <- eValueMat <- eValueFutMat <- matrix(NA, n_sources, n_max)

  for (s in seq_along(allSources)) {

    someDat <- dat[dat$source == allSources[s], ]
    n_obs <- nrow(someDat)

    for (i in seq_len(n_obs)) {

      tempDat <- someDat[1:i, ]

      sd1 <- sqrt(var(tempDat$variable1))
      sd2 <- sqrt(var(tempDat$variable2))
      someN <- nMat[s,i] <- nrow(tempDat)

      if (someN > 3 && sd1 > 0 && sd2 > 0) {

        someR <- rMat[s,i] <- cor(tempDat$variable1, tempDat$variable2)

        someZ <- zMat[s,i] <- atanh(someR) / sqrt(1/(someN-3))

        tempRes <- saviZTestStat(
          z = someZ,
          n1 = someN,
          parameter = designObj$parameter,
          eType = designObj$eType,
          sigma = designObj$sigma
        )

        eValueMat[s,i] <- tempRes$eValue

        tempRes <- saviFutilityZStat(
          z = someZ,
          n1 = someN,
          parameter = designObj$futilityResult$parameter,
          sigma = designObj$sigma
        )

        eValueFutMat[s,i] <- tempRes$eValue
      }
    }
  }

  list(
    eValueMat = eValueMat,
    eValueFutMat = eValueFutMat,
    zMat = zMat,
    rMat = rMat,
    nMat = nMat
  )
}

res <- compute_evalue_mats(dat, allSources, designObj)

eValueMat <- res$eValueMat
eValueFutMat <- res$eValueFutMat

# I do a function to fill the rest of the vector with the last e-value for the Scenario 2
locf_row <- function(x) {
  for (i in seq_along(x)) {
    if (is.na(x[i]) && i > 1) {
      x[i] <- x[i-1]
    }
  }
  x
}

eValueMat2 <- t(apply(eValueMat, 1, locf_row))
eValueFutMat2 <- t(apply(eValueFutMat, 1, locf_row))


# Aggregation for Scenario 2
aggregate_mult <- function(mat) {

  n_max <- ncol(mat)
  agg <- numeric(n_max)

  for (t in seq_len(n_max)) {
    vals <- mat[, t]
    vals <- vals[!is.na(vals)]
    agg[t] <- prod(vals)
  }

  agg
}

aggE <- aggregate_mult(eValueMat2)
aggFut <- aggregate_mult(eValueFutMat2)


which(aggE >= 1/alpha)
which(aggFut <= betaFutility)

plot(aggE, log="y", type="l")
lines(aggFut, col="red")
abline(h=betaFutility)
abline(h=1/alpha, col="blue")

plot(aggFut, log="y", type="l")


# Stopping time

permute_sources <- function(dat, allSources) {

  dat_perm <- dat

  for (s in allSources) {

    idx <- which(dat$source == s)
    dat_perm[idx, ] <- dat[idx[sample(length(idx))], ]

  }

  dat_perm
}

M <- 100
results <- vector("list", M)
stops <- numeric(M)
type <- character(M)

for (m in 1:M) {
  print(m)
  dat_perm <- permute_sources(dat, allSources)

  results[[m]] <- compute_evalue_mats(dat_perm, allSources, designObj)
  eValueMat <- results[[m]]$eValueMat
  eValueFutMat <- results[[m]]$eValueFutMat

  eValueMat2 <- t(apply(eValueMat, 1, locf_row))
  eValueFutMat2 <- t(apply(eValueFutMat, 1, locf_row))

  aggE <- aggregate_mult(eValueMat2)
  aggFut <- aggregate_mult(eValueFutMat2)

  stopE <- which(aggE >= 1/alpha)[1]
  stopF <- which(aggFut <= betaFutility)[1]
  stopE <- ifelse(is.na(stopE), Inf, stopE)
  stopF <- ifelse(is.na(stopF), Inf, stopF)

  stops[m] <- min(stopE, stopF)

  if(stopE < stopF) type[m] <- "reject"
  else if(stopF < stopE) type[m] <- "futility"
  else type[m] <- "none"
}
mean(stops[is.finite(stops)])
var(stops[is.finite(stops)])


mean(type == "futility")# Over 90% times we stopped for futility




# Scenario 3
last_obs <- apply(eValueMat, 1, function(x) max(which(!is.na(x)))) # Number of samples per source


for (s in unique(dat$source)){
  print(sum(dat$source == s))
}

aggregate_random <- function(mat, last_obs) {

  n_sources <- nrow(mat)

  idx <- rep(1, n_sources)
  agg <- numeric(sum(last_obs))
  t <- 1

  while (any(idx <= last_obs)) {

    available <- which(idx <= last_obs)
    s <- sample(available, 1)

    vals <- sapply(seq_len(n_sources), function(i) {
      mat[i, min(idx[i], last_obs[i])]
    })

    agg[t] <- prod(vals, na.rm = TRUE)
    t <- t + 1

    idx[s] <- idx[s] + 1
  }

  agg
}

aggRandE <- aggregate_random(eValueMat2, last_obs)
aggRandFut <- aggregate_random(eValueFutMat2, last_obs)


which(aggRandE >= 1/alpha)
which(aggRandFut <= betaFutility)

plot(aggRandE, log="y", type="l")
lines(aggRandFut, col="red")
abline(h=betaFutility)
abline(h=1/alpha, col="blue")


M <- 1000
results <- vector("list", M)
stops <- numeric(M)
type <- character(M)

# A bit heavy computationnally
for (m in 1:M) {
  print(m)
  dat_perm <- permute_sources(dat, allSources)

  results[[m]] <- compute_evalue_mats(dat_perm, allSources, designObj)
  eValueMat <- results[[m]]$eValueMat
  eValueFutMat <- results[[m]]$eValueFutMat

  #eValueMat2 <- t(apply(eValueMat, 1, locf_row))
  #eValueFutMat2 <- t(apply(eValueFutMat, 1, locf_row))
  aggRandE <- aggregate_random(eValueMat2, last_obs)
  aggRandFut <- aggregate_random(eValueFutMat2, last_obs)


  stopE <- which(aggRandE >= 1/alpha)[1]
  stopF <- which(aggRandFut <= betaFutility)[1]
  stopE <- ifelse(is.na(stopE), Inf, stopE)
  stopF <- ifelse(is.na(stopF), Inf, stopF)

  stops[m] <- min(stopE, stopF)

  if(stopE < stopF) type[m] <- "reject"
  else if(stopF < stopE) type[m] <- "futility"
  else type[m] <- "none"
}
mean(stops[is.finite(stops)])# 202.401
var(stops[is.finite(stops)])#14612.77


mean(type == "futility")# 0.094
mean(type == "reject")# 0.906

type
stops



















# The rest is deprecated


n1 <- length(dat$variable1)
atanh(stat.test$estimate) / sqrt(1 / (n1 - 3))

allSources <- unique(dat$source)
cor(dat$variable1, dat$variable2)


r1 <- cor(dat$variable1, dat$variable2)
n1 <- length(dat$variable1)

z <- atanh(r1) / sqrt(1 / (n1 - 3))
alpha <- 0.05
sides <- 2
conf.low <- tanh(atanh(r1) - (qnorm(1 - (alpha / sides)) * sqrt((1 / (n1 - 3)))))
conf.high <- tanh(atanh(r1) + (qnorm(1 - (alpha / sides)) * sqrt((1 / (n1 - 3)))))
p <- 2 * (1 - pnorm(abs(z)))




betaFutility <- 0.2


# Original study: Lower bound of effect size found in the original study
deltaMin <- 0.25*0.68#(4.95-3.75)/1.32 # r = -0.21, d= -0.43 r = 0.25, d= 0.52




deltaMin <- atanh(deltaMin) / sqrt(1 / (1548 - 3))

freqDesign <- power.t.test(delta=atanh(0.29*0.7), sd =1, alternative="two.sided",
                           power=0.8)
# Prospective frequentist analysis
#freqDesign <- power.t.test(delta=deltaMin, alternative="one.sided",
#                          power=0.8)

# Prospective e-value analysis


designObj2Two <- designSaviZ(meanDiffMin=deltaMin, beta=0.2,
                             testType="oneSample",
                             alternative="twoSided", seed=5)


dat <- grahamData$cleanDataFilter
head(ML2.var)

head(dat)
cor(dat$variable1, dat$variable2)
allSources <- unique(dat$source)


z <- atanh(cor(dat$variable1, dat$variable2)) * sqrt(length(dat$variable1) - 3)
e_tot <- saviZTest(z, designObj = designObj2Two)
# Scenario sequential studies


eval <- numeric(length(dat$variable1))
eval[1:3] <- 1
z_vect <- numeric(length(dat$variable1)-3)
e_seq <- numeric(length(dat$variable1)-3)
for (i in 4:length(dat$variable1)) {
  r <- cor(dat$variable1[1:i], dat$variable2[1:i])
  z_vect[i-3] <- atanh(r) * sqrt(i - 3)
  eval[i] <- saviZTest(z_vect[i-3], designObj = designObj2Two)$eValue
  e_seq[i-3] <- saviZTest(z_vect[1:(i-3)], designObj = designObj2Two)$eValue
}
eval[6000]
cor(dat$variable1[6000:6900], dat$variable2[6000:6900])
cor(dat$variable1, dat$variable2)
length(dat$variable2)
e_tot <- saviZTest(z_vect, designObj = designObj2Two)
firstTime <- min(which(e_seq>= 20))
min(which(eval>= 20))
e_tot$eValue
eval

r <- cor(dat$variable1, dat$variable2)
z_vect <- atanh(r) * sqrt(length(dat$variable1) - 3)
saviZTest(z_vect, designObj = designObj2Two)$eValue

for (i in 1:length(allSources)) {
  someDat <- dat[dat$source==allSources[i], ]

  ## Data -----
  x <- someDat$variable1
  y <- someDat$variable2
  print(cor(x, y))
  print(length(x))
}




# CORRELATION STUDY


library(ggplot2)
library(dplyr)

plotDat <- dat %>%
  filter(source %in% allSources) %>%
  select(source, variable1, variable2)

cors <- plotDat %>%
  group_by(source) %>%
  summarise(
    r = cor(variable1, variable2, use = "complete.obs"),
    n = sum(complete.cases(variable1, variable2))
  )

ggplot(plotDat, aes(variable1, variable2)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "steelblue") +
  facet_wrap(~ source, scales = "free") +
  geom_text(
    data = cors,
    aes(
      x = -Inf, y = Inf,
      label = paste0("r = ", round(r, 2), "\nn = ", n)
    ),
    hjust = -0.1, vjust = 1.1,
    inherit.aes = FALSE,
    size = 4
  ) +
  theme_minimal() +
  labs(x = "variable1", y = "variable2")


ggplot(cors, aes(x = reorder(source, r), y = r)) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  theme_minimal() +
  labs(
    x = "Source",
    y = "Correlation (r)",
    title = "Correlation between variable1 and variable2 by source"
  )


threshold <- 20
deltaMin

designObj2 <- designSaviZ(meanDiffMin=atanh(0.29*0.7), beta=0.2,
                          testType="oneSample",
                          alternative="twoSided", seed=5, futility = TRUE)
plot(designObj2)
sources <- unique(dat$source)
results <- lapply(sources, function(g) {

  r <- cor(dat$variable1[dat$source == g], dat$variable2[dat$source == g])
  z <- atanh(r) * sqrt(sum(dat$source == g) - 3)
  ev <- saviZTest(z, design = designObj2)$eValue
  ef <- saviZTest(z, design = designObj2)$eValueFut
  status <- ev > threshold


  data.frame(
    source = g,
    status = status,
    n = sum(dat$source == g),
    ev = ev,
    ef = ef,
    z = z
  )
})

results <- do.call(rbind, results)
results
mean(results$status)
prod(results$ev)
prod(results$ef)
plot(results$ev)


r <- cor(dat$variable1, dat$variable2)
z <- atanh(r) * sqrt(length(dat$variable1) - 3)
ev <- saviZTest(z, design = designObj2)$eValue




atanh(r1)/sqrt(1/(n1-3))
