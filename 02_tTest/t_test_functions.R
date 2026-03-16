# Functions for performing the full frequentist and sequential analyses of the t-test studies

remove_degenerate_sources <- function(ECDF){
  source_list <- sort(unique(ECDF$source))
  variable <- colnames(ECDF)[2]
  factor <- colnames(ECDF)[3]
  
  for (i in seq_along(source_list)) {
    current_source <- source_list[i]
    current_factors <- unique(ECDF[ECDF$source == current_source,][[factor]])
      if(length(current_factors)<2){
        print(paste("Source",current_source,"contains < 2 factors - excluded from analysis."))
        ECDF <- ECDF[ECDF$source != current_source, ]
      }
  }
  return(ECDF)
}

# Frequentist analysis
frequentist_t_test <- function(ECDF, subStudyName, global = FALSE, var_equal = FALSE) {
  
  variable <- colnames(ECDF)[2]
  factor <- colnames(ECDF)[3]
  
  # Choose subStudy or Global
  if (!global) {
    df <- ECDF[ECDF$source == subStudyName, ]
  } else {
    df <- ECDF
  }
  
  # Two-sample (Welch) t-test
  res <- t.test(df[[variable]] ~ df[[factor]], var.equal = var_equal)
  
  # Group SDs
  sds <- tapply(df[[variable]], df[[factor]], sd, na.rm = TRUE)
  
  # Group sample sizes
  ns <- tapply(df[[variable]], df[[factor]], function(x) sum(!is.na(x)))
  
  # Build results row
  re <- data.frame(
    source = subStudyName,
    n_group1   = ns[[1]],
    mean_group1 = res$estimate[1],
    sd_group1  = sds[[1]],
    n_group2   = ns[[2]],
    mean_group2 = res$estimate[2],
    sd_group2  = sds[[2]],
    t.statistic = unname(res$statistic),
    p.value     = res$p.value,
    row.names = NULL
  )
  
  return(re)
}

full_freq_t_test_analysis <- function(ECDF, var_equal = FALSE){
  
  source_list <- sort(unique(ECDF$source))
  n_substudies <- length(source_list)
  
  frequentist_results <- data.frame(
    source      = character(n_substudies+1),
    n_group_1   = integer(n_substudies+1),
    mean_group1 = numeric(n_substudies+1),
    sd_group_1  = numeric(n_substudies+1),
    n_group_2   = integer(n_substudies+1),
    mean_group2 = numeric(n_substudies+1),
    sd_group_2  = numeric(n_substudies+1),
    t.statistic = numeric(n_substudies+1),
    p.value     = numeric(n_substudies+1)
  )

  
  for (i in seq_along(source_list)) {
      frequentist_results[i,] <- frequentist_t_test(ECDF, subStudyName = source_list[i], global = FALSE, var_equal = var_equal)
      }

  # add GLOBAL row as the last row, no bootstrapping here !!
  global_row <- frequentist_t_test(ECDF, subStudyName = "GLOBAL", global = TRUE, var_equal = var_equal)
  frequentist_results[n_substudies+1,] <- global_row
  
  return(frequentist_results)
}

##############################################################################################################################################################
# Sequential analysis

metaType3Helper <- function(v, M) {
  nz <- v != 0
  vals <- M[cbind(which(nz), v[nz])]
  return(sum(log(vals)))
}


full_seq_t_test_analysis <- function(PECDF, alpha, betaFutility, deltaMin, esMinFutility, varEqual = FALSE){
  
  original_source_list <- unique(PECDF$source) # only for meta type 1
  source_list <- sort(unique(PECDF$source))
  n_substudies <- length(source_list)
  variable <- colnames(PECDF)[2]
  factor <- colnames(PECDF)[3]
  distinct_factors <- unique(PECDF[[factor]])
  factor1 <- distinct_factors[1]
  factor2 <- distinct_factors[2]
  
  sequential_results <- data.frame(
    source_id      = 1:n_substudies,
    source         = character(n_substudies),
    stopped_times  = integer(n_substudies),
    final_times    = integer(n_substudies),
    percent_saved  = numeric(n_substudies),
    stopped_for_R  = numeric(n_substudies),
    stopped_for_N  = numeric(n_substudies),
    stopped_for_F  = numeric(n_substudies)
  )
  
  # savi design object
  DO <- designSaviT(deltaMin=deltaMin, testType="twoSample", alternative="twoSided", futility = TRUE, wantSampling = FALSE, varEqual = varEqual)
  # matrices containing the e/f-value sequences
  n_max_max <- max(table(PECDF$source))
  eValueMat <- matrix(1, nrow = n_substudies, ncol = n_max_max)
  fValueMat <- matrix(1, nrow = n_substudies, ncol = n_max_max)
  
  for (i in seq_along(source_list)) {
    current_source <- source_list[i]
    subStudy <- PECDF[PECDF$source == current_source, ]
    
    sequential_results$source[i] <- current_source
    sequential_results$final_times[i] <- nrow(subStudy)
    x <- subStudy[subStudy[[factor]] == factor1, ][[variable]]
    y <- subStudy[subStudy[[factor]] == factor2, ][[variable]]
    n1 <- 0
    n2 <- 0
    found_stopping_time <- FALSE
    
    for (j in seq_len(nrow(subStudy))){
      
      if(subStudy[j,factor] == factor1){
        n1 <- n1 + 1
      }
      else{
        n2 <- n2 + 1
      }
      
      if (n1 < 2 || n2 < 2 || all(x[1:n1] == x[1]) || all(y[1:n2] == y[1])){next}

      res <- saviTTest(x[1:n1], y[1:n2], designObj = DO, sequential = FALSE, varEqual = varEqual, futility = TRUE, esMinFutility = esMinFutility)
      eValue <- unname(res$eValue)
      fValue <- unname(res$eValueFut)
      eValueMat[i,j] <- eValue
      fValueMat[i,j] <- fValue
      
      if (!found_stopping_time && eValue >= 1/alpha){
        sequential_results$stopped_times[i] <- j
        sequential_results$stopped_for_R[i] <- 1
        found_stopping_time <- TRUE
      }
      if (!found_stopping_time && fValue <= betaFutility){
        sequential_results$stopped_times[i] <- j
        sequential_results$stopped_for_F[i] <- 1
        found_stopping_time <- TRUE
      }
      if(j == nrow(subStudy)){
        if(!found_stopping_time){
        sequential_results$stopped_times[i] <- j
        sequential_results$stopped_for_N[i] <- 1
        }
        if (j < n_max_max){
          eValueMat[i,(j+1):n_max_max] <- eValue
          fValueMat[i,(j+1):n_max_max] <- fValue
        }
      }
    }
  }
  
  sequential_results$percent_saved <- (1-sequential_results$stopped_times/sequential_results$final_times)*100
  
  source_to_index <- setNames(sequential_results$source_id, sequential_results$source)
  desired_indices <- source_to_index[original_source_list]
  eValueMatReordered <- eValueMat[desired_indices, ]
  fValueMatReordered <- fValueMat[desired_indices, ]
  
  # NOTE: all meta e-processes are in log space
  
  metaEType1 <- cumsum(log(eValueMatReordered[,n_max_max])) # no longer uses alphabetic order
  metaFType1 <- cumsum(log(fValueMatReordered[,n_max_max]))
  
  metaEType2 <- colSums(log(eValueMat))
  metaFType2 <- colSums(log(fValueMat))
  
  metaEType3 <- rep(0,nrow(PECDF))
  metaFType3 <- rep(0,nrow(PECDF))
  
  # Stopped version of meta type 2
  stoppedEValueMat <- eValueMat
  stoppedFValueMat <- fValueMat
  for (i in seq_along(source_list)){
    stoppedEValueMat[i,sequential_results$stopped_times[i]:n_max_max] <- stoppedEValueMat[i,sequential_results$stopped_times[i]]
    stoppedFValueMat[i,sequential_results$stopped_times[i]:n_max_max] <- stoppedFValueMat[i,sequential_results$stopped_times[i]]
  }
  
  stoppedMetaEType2 <- colSums(log(stoppedEValueMat))
  stoppedMetaFType2 <- colSums(log(stoppedFValueMat))
  
  # Worst case permutation version of meta type-1
  worstCaseMetaEType1 <- cumsum(log(sort(eValueMat[,n_max_max])))
  worstCaseMetaFType1 <- cumsum(log(sort(fValueMat[,n_max_max], decreasing = TRUE)))
  
  
  
  # IMPORTANT, TYPE3 SHOULD NOT BE USED on an PECDF which is ordered by factor values !
  
  vector_index <- rep(0,n_substudies)
  for (i in seq_along(PECDF$source)){
    current_source <- PECDF$source[i]
    current_source_id <- sequential_results[sequential_results$source==current_source, ]$source_id
    vector_index[current_source_id] <- vector_index[current_source_id] + 1
    metaEType3[i] <- metaType3Helper(vector_index, eValueMat)
    metaFType3[i] <- metaType3Helper(vector_index, fValueMat)
  }
  
  return(list(sequential_results = sequential_results, eValueMat = eValueMat, fValueMat = fValueMat,
              metaEType1 = metaEType1, metaFType1 = metaFType1, metaEType2 = metaEType2, metaFType2 = metaFType2,
              metaEType3 = metaEType3, metaFType3 = metaFType3,
              worstCaseMetaEType1 = worstCaseMetaEType1, worstCaseMetaFType1 = worstCaseMetaFType1,
              stoppedMetaEType2 = stoppedMetaEType2, stoppedMetaFType2 = stoppedMetaFType2))
}


#########

avg_meta_results <- function(PECDF, alpha, betaFutility, deltaMin, esMinFutility, varEqual = FALSE, n_permutations = 1, meta_alpha=NULL, meta_betaFutility=NULL){
  
  original_source_list <- unique(PECDF$source)
  n_substudies <- length(original_source_list)
  
  if(is.null(meta_alpha)){meta_alpha <- alpha/n_substudies}
  if(is.null(meta_betaFutility)){meta_betaFutility <- betaFutility/n_substudies}
  
  metaEType1Realizations <- vector("list", n_permutations)
  metaEType2Realizations <- vector("list", n_permutations)
  metaEType3Realizations <- vector("list", n_permutations)
  stoppedMetaEType2Realizations <- vector("list", n_permutations)
  
  metaFType1Realizations <- vector("list", n_permutations)
  metaFType2Realizations <- vector("list", n_permutations)
  metaFType3Realizations <- vector("list", n_permutations)
  stoppedMetaFType2Realizations <- vector("list", n_permutations)
  
  
  metaType1StoppingTimes <- rep(-1, n_permutations)
  metaType2StoppingTimes <- rep(-1, n_permutations)
  metaType3StoppingTimes <- rep(-1, n_permutations)
  stoppedMetaType2StoppingTimes <- rep(-1, n_permutations)
  
  for (i in seq_len(n_permutations)){
    cat("Started iteration:",i,"\n")
    
    PECDF <- PECDF[sample(nrow(PECDF)),]
    sequential_results_list <- full_seq_t_test_analysis(PECDF, alpha, betaFutility, deltaMin, esMinFutility, varEqual)
    
    metaEType1Realizations[[i]] <- sequential_results_list$metaEType1
    metaFType1Realizations[[i]] <- sequential_results_list$metaFType1
    meta1FirstPassage <- which(sequential_results_list$metaEType1 >= log(1/meta_alpha) | sequential_results_list$metaFType1 <= log(meta_betaFutility))[1]
    if (is.na(meta1FirstPassage)){
      metaType1StoppingTimes[i] <- length(sequential_results_list$metaEType1)
      }
    else{
      metaType1StoppingTimes[i] <- meta1FirstPassage
    }
    
    metaEType2Realizations[[i]] <- sequential_results_list$metaEType2
    metaFType2Realizations[[i]] <- sequential_results_list$metaFType2
    meta2FirstPassage <- which(sequential_results_list$metaEType2 >= log(1/meta_alpha) | sequential_results_list$metaFType2 <= log(meta_betaFutility))[1]
    if (is.na(meta2FirstPassage)){
      metaType2StoppingTimes[i] <- length(sequential_results_list$metaEType2)
    }
    else{
      metaType2StoppingTimes[i] <- meta2FirstPassage
    }
    
    metaEType3Realizations[[i]] <- sequential_results_list$metaEType3
    metaFType3Realizations[[i]] <- sequential_results_list$metaFType3
    meta3FirstPassage <- which(sequential_results_list$metaEType3 >= log(1/meta_alpha) | sequential_results_list$metaFType3 <= log(meta_betaFutility))[1]
    if (is.na(meta3FirstPassage)){
      metaType3StoppingTimes[i] <- length(sequential_results_list$metaEType3)
    }
    else{
      metaType3StoppingTimes[i] <- meta3FirstPassage
    }
    
    stoppedMetaEType2Realizations[[i]] <- sequential_results_list$stoppedMetaEType2
    stoppedMetaFType2Realizations[[i]] <- sequential_results_list$stoppedMetaFType2
    stoppedMeta2FirstPassage <- which(sequential_results_list$stoppedMetaEType2 >= log(1/meta_alpha) | sequential_results_list$stoppedMetaFType2 <= log(meta_betaFutility))[1]
    if (is.na(stoppedMeta2FirstPassage)){
      stoppedMetaType2StoppingTimes[i] <- length(sequential_results_list$stoppedMetaEType2)
    }
    else{
      stoppedMetaType2StoppingTimes[i] <- stoppedMeta2FirstPassage
    }
  }

  return(list(metaEType1Realizations = metaEType1Realizations,
              metaEType2Realizations = metaEType2Realizations,
              metaEType3Realizations = metaEType3Realizations,
              stoppedMetaEType2Realizations = stoppedMetaEType2Realizations,
              metaFType1Realizations = metaFType1Realizations,
              metaFType2Realizations = metaFType2Realizations,
              metaFType3Realizations = metaFType3Realizations,
              stoppedMetaFType2Realizations = stoppedMetaFType2Realizations,
              metaType1StoppingTimes = metaType1StoppingTimes,
              metaType2StoppingTimes = metaType2StoppingTimes, 
              metaType3StoppingTimes = metaType3StoppingTimes,
              stoppedMetaType2StoppingTimes = stoppedMetaType2StoppingTimes))
}
