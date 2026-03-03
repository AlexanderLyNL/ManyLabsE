library(devtools)
library(plyr)
library(rio)
library(tidyverse)
library(reshape2)
library(stats)
# rm(list = ls())



# setwd("/home/areyerol/Bureau/futility/ManyLabsE/03_zTest")

# library(safestats)
# repo.path <- "/home/areyerol/Bureau/git/safestats-futility88"
# load_all(repo.path)



sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("are", "Are")) "Bureau/futility"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("are", "Are")) "/home/areyerol/Bureau/futility/ManyLabsE/03_zTest" #"~/Desktop/git/manyLabsE/02_tTest/"

sourcePath <- if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "/Desktop/git/"
myWd <-  if (substr(system("whoami", intern=TRUE), 1, 3) %in% c("ale", "Ale")) "~/Desktop/git/manyLabsE/02_tTest/"

project.root <- file.path("~", sourcePath, "ManyLabsE")
OSFdata.root <- file.path(project.root, "OSFdata")

source(file.path(project.root, "00_utils", "WYQ_manylabRs_SOURCE.R"))
source(file.path(project.root, "00_utils", "helpers.R"))

# ANALYSIS INFO ----


study.description      <- 'Choosing or Rejecting (Shafir, 1993)'
analysis.unique.id     <- 74
analysis.name          <- 'Shafir.1'
analysis.type          <- 1
analysis.type.name     <- 'study_global_include'
analysis.type.groups   <- 'Source.Global'
Nmin.raw               <- 30
Nmin.cond              <- 15
subset                 <- 'all'
onlineTables           <- TRUE
staticData             <- TRUE
saveAll                <- TRUE
overWrite              <- TRUE
#OSFdata.root           <- file.path('~','OSFdata')
analysis.root          <- file.path(OSFdata.root,study.description,analysis.name,'Global')
outdir                 <- list(Data = file.path(analysis.root,'Data'), Results = file.path(analysis.root,'Results'))


# This function will be used to change the raw dataset to a dataset ready for analysis

varfun.Shafir.1


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

ML2.df <- ML2.df  %>% dplyr::select(1,6,285,290,804,903,904,905,906,907,908,909,910,911,912,913,914,937,938,939) %>% dplyr::filter(is.character(source))



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

        # To see the function code type:varfun.Shafir.1, or lookup in manylabRs_SOURCE.R
        ML2.var[[g]] <- varfun.Shafir.1(ML2.sr[[g]])


        # Check equal variance assumption
        if(!is.na(testVarEqual)){
          if(testVarEqual){
            logtxt <- paste(analysis.unique.id,ML2.key$study.analysis,'-', runGroups[g])
            ML2.in$stat.params$var.equal <- decide.EqualVar(ML2.var[[g]],ML2.in$study.vars.labels, ML2.key, group = logtxt) # don't pass the cleanData frame
          }}

        # Run the analysis according to ML2.key: 'stat.test'
        stat.params <<- ML2.in$stat.params


        stat.test   <- try.CATCH(with(ML2.var[[g]],z.test(x = ParentB, pi = .5, N = sum(N, na.rm=TRUE), proportion = TRUE, alternative = stat.params$alternative)))


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


# Alexander -------
dat <- addUniqueIds(ML2.var, ML2.df)
dat <- checkUniqueIds(dat)

# Here -------


xName <- "High"
yName <- "Low"
alpha <- 0.05
betaFutility <- 0.2
deltaMin <- 0.68


dat$variable1_binary <- as.numeric(dat$variable1 == "Parent B")
head(dat)

(1-pnorm(2.48))*2 # Original p-value
(1-pnorm(6.1))*2 # Replicated p-value

(1-pnorm(abs((((0.476+0.455)/2)-0.5)/sqrt(0.25/7901))))*2 # Replicated p-value
stat.test
z.test(x = shafirData$ParentB, pi = .5, N = 7901, proportion = TRUE)# sigma = 0.25




for (g in unique(shafirData$cleanDataFilter$source)){
  print(sum(dat$variable1_binary[dat$source==g])/sum(dat$source==g))
}

library(dplyr)
library(ggplot2)

study_summary <- dat %>%
  group_by(source) %>%
  summarise(
    n = n(),
    mean = mean(variable1_binary, na.rm = TRUE),
    se = sqrt(mean * (1 - mean) / n),   # binomial SE
    lower = mean - 1.96 * se,
    upper = mean + 1.96 * se
  )
ggplot(study_summary,
       aes(x = reorder(source, mean), y = mean)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.2, alpha = 0.6) +
  geom_hline(yintercept = 0.5,
             linetype = "dashed", color = "red") +
  coord_flip() +
  labs(
    x = "Study",
    y = "Proportion (Mean of Binary Variable)",
    title = "Study-Level Proportions with 95% Confidence Intervals"
  ) +
  theme_minimal()


library(meta)

m <- metagen(

  TE = study_summary$mean,

  seTE = study_summary$se,

  studlab = study_summary$source

)

forest(m) # forest plot

funnel(m) # funnel plot





betaFutility <- 0.2


# Original study: Lower bound of effect size found in the original study


deltaMin <- ((.64+.55)/2-0.5)* 0.7 #((.64+.55)/2-0.5)/sqrt(0.25/170)*0.68

# Prospective frequentist analysis
freqDesign <- power.t.test(delta=deltaMin, sd = 0.5, alternative="two.sided",
                           power=0.8)
((.64+.55)/2-0.5)/sqrt(0.25/170)
# Prospective e-value analysis

designObj2 <- designSaviZ(meanDiffMin=deltaMin, beta=0.2,
                          testType="oneSample", sigma = 0.5,
                          alternative="twoSided", seed=5,
                          futility = TRUE)
plot(designObj2)



dat$variable1_binary_norm = dat$variable1_binary-0.5
e_test = saviZTest(dat$variable1_binary_norm ,designObj = designObj2, futility=TRUE)
e_test$eValue
e_test$eValueFut

e_test_out <- saviZTest(dat$variable1_binary_norm[dat$source!='jos'] ,designObj = designObj2)
e_test_out$eValue


library(dplyr)
library(purrr)
library(ggplot2)

study_results <- dat %>%
  group_by(source) %>%
  summarise(
    n = n(),
    mean = mean(variable1_binary_norm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # classical z test against 0
    z = mean / sqrt(0.25 / n),
    p_value = 2 * (1 - pnorm(abs(z))),

    # e-values from saviZTest
    e_obj = map(source, ~ saviZTest(
      dat$variable1_binary_norm[dat$source == .x],
      design = designObj2,
      futility = TRUE
    )),
    e_value = map_dbl(e_obj, "eValue"),
    e_value_fut = map_dbl(e_obj, "eValueFut")
  )


total_e      <- prod(study_results$e_value)
total_e_fut  <- prod(study_results$e_value_fut)

prop_sig     <- mean(study_results$p_value < 0.05)

total_e > 20
total_e_fut < betaFutility
prop_sig


# DISTRIBUTION OF THE EFFECT
ggplot(study_results, aes(x = mean)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribution of Study-Level Effects",
       x = "Mean Effect",
       y = "Count") +
  theme_minimal()

# E-values

ggplot(study_results,
       aes(x = reorder(source, e_value), y = e_value)) +
  geom_point() +
  geom_hline(yintercept = 20, linetype = "dashed") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "E-values per Study",
       x = "Study",
       y = "E-value (log scale)") +
  theme_minimal()

study_results <- study_results %>%
  arrange(desc(e_value)) %>%
  mutate(cum_e = cumprod(e_value))

# CUMULATIVE

ggplot(study_results,
       aes(x = seq_along(cum_e), y = cum_e)) +
  geom_line() +
  geom_hline(yintercept = 20, linetype = "dashed") +
  scale_y_log10() +
  labs(title = "Cumulative Product of E-values",
       x = "Number of Studies Included",
       y = "Cumulative E-value (log scale)") +
  theme_minimal()


# FIRST TIME ANALYSIS


sources <- unique(dat$source)
first_reject <- sapply(sources, function(g) {
  ev <- saviZTest(dat$variable1_binary_norm[dat$source == g], design = designObj2,
                  futility = TRUE)$eValueVec

  idx <- which(ev > 20)

  if (length(idx) == 0) NA else min(idx)
})

results <- data.frame(
  source = sources,
  first_reject = first_reject,
  rejected = !is.na(first_reject)
)
mean(results$rejected)


first_reject_fut <- sapply(sources, function(g) {
  ev <- saviZTest(dat$variable1_binary_norm[dat$source == g], design = designObj2,
                  futility = TRUE)$eValueFutVec

  idx <- which(ev < betaFutility)

  if (length(idx) == 0) NA else min(idx)
})

results_fut <- data.frame(
  source = sources,
  first_reject = first_reject_fut,
  rejected = !is.na(first_reject_fut)
)
mean(results_fut$rejected)


library(ggplot2)

ggplot(results[results$rejected,],
       aes(x = first_reject)) +
  geom_histogram(bins = 30) +
  labs(title = "Distribution of First Rejection Time",
       x = "Time of First e-value > 20",
       y = "Number of Sources") +
  theme_minimal()



threshold <- 20
sources <- unique(dat$source)

results <- lapply(sources, function(g) {

  x <- dat$variable1_binary_norm[dat$source == g]
  ev <- saviZTest(x, futility = TRUE, design = designObj2)$eValueVec

  idx <- which(ev > threshold)

  if (length(idx) == 0) {
    time <- length(ev)
    status <- 0
  } else {
    time <- min(idx)
    status <- 1
  }

  data.frame(
    source = g,
    time = time,
    status = status,
    n = length(x),
    mean = mean(x, na.rm = TRUE)
  )
})

results <- do.call(rbind, results)

# OUTLIER

study_results$e_value[1] <- 1
study_results$e_value_fut[1] <- 1
total_e      <- prod(study_results$e_value)
total_e_fut  <- prod(study_results$e_value_fut)

prop_sig     <- mean(study_results$p_value < 0.05)

total_e > 20
total_e_fut < betaFutility
prop_sig


# DISTRIBUTION OF THE EFFECT
ggplot(study_results, aes(x = mean)) +
  geom_histogram(bins = 30) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Distribution of Study-Level Effects",
       x = "Mean Effect",
       y = "Count") +
  theme_minimal()

# E-values

ggplot(study_results,
       aes(x = reorder(source, e_value), y = e_value)) +
  geom_point() +
  geom_hline(yintercept = 20, linetype = "dashed") +
  scale_y_log10() +
  coord_flip() +
  labs(title = "E-values per Study",
       x = "Study",
       y = "E-value (log scale)") +
  theme_minimal()

study_results <- study_results %>%
  arrange(desc(e_value)) %>%
  mutate(cum_e = cumprod(e_value))

# CUMULATIVE

ggplot(study_results,
       aes(x = seq_along(cum_e), y = cum_e)) +
  geom_line() +
  geom_hline(yintercept = 20, linetype = "dashed") +
  scale_y_log10() +
  labs(title = "Cumulative Product of E-values",
       x = "Number of Studies Included",
       y = "Cumulative E-value (log scale)") +
  theme_minimal()



dat[dat$source=="jos", ]
