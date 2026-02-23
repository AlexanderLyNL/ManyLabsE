# Functions for performing the full frequentist and sequential analyses of the t-test studies

# Frequentist analysis
frequentist_t_test <- function(ECDF,
                               variable,
                               factor,
                               subStudyIndex = 1,
                               global = FALSE,
                               var_equal = FALSE, 
                               bootstrap = FALSE,
                               n_plan = NULL) {
  
  # Choose subStudy or Global
  if (!global) {
    df <- ECDF[ECDF$study.order == subStudyIndex, ]
  } else {
    df <- ECDF
    subStudyIndex <- max(ECDF$study.order)+1
  }
  
  
  # Bootstrapping, take only a 2*n_plan-sized subsample from the empirical distribution, each group n_plan samples
  if(bootstrap){
    distinct_factors <- unique(ECDF[[factor]])
    factor1 <- distinct_factors[1]
    factor2 <- distinct_factors[2]
    
    df1 <- df[df[[factor]] == factor1, ]
    s1 <- df1[sample(nrow(df1), n_plan, replace = TRUE), ]
    df2 <- df[df[[factor]] == factor2, ]
    s2 <- df2[sample(nrow(df2), n_plan, replace = TRUE), ]
    df <- rbind(s1, s2)
  }
  
  # Two-sample (Welch) t-test
  res <- t.test(df[[variable]] ~ df[[factor]], var.equal = var_equal)
  
  # Group SDs
  sds <- tapply(df[[variable]], df[[factor]], sd, na.rm = TRUE)
  
  # Group sample sizes
  ns <- tapply(df[[variable]], df[[factor]], function(x) sum(!is.na(x)))
  
  # Build results row
  re <- data.frame(
    study.order = subStudyIndex,
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

full_freq_t_test_analysis <- function(ECDF, variable, factor, var_equal = FALSE, bootstrap = FALSE, n_bootstrap=1, deltaMin = NULL, beta = NULL){
  
  n_substudies <- max(ECDF$study.order)
  
  frequentist_results <- data.frame(
    study.order = integer(n_substudies+1),
    n_group_1   = integer(n_substudies+1),
    mean_group1 = numeric(n_substudies+1),
    sd_group_1  = numeric(n_substudies+1),
    n_group_2   = integer(n_substudies+1),
    mean_group2 = numeric(n_substudies+1),
    sd_group_2  = numeric(n_substudies+1),
    t.statistic = numeric(n_substudies+1),
    p.value     = numeric(n_substudies+1)
  )
  
  if(bootstrap){
   if( is.null(deltaMin) || is.null(beta)){
     print("Must specify deltaMin and beta for bootstrap method.")
     return(NULL) 
   }
    n_plan <- ceiling(power.t.test(delta=deltaMin, sd= 1, alternative="two.sided", power=1-beta)$n)
  }
  else{
    n_plan <- NULL
  }
  
  ### TODO: think about: does it make sense to just compute this table for all k in seq_len(n_bootstrap) and then average over them - possibly without global analysis row
  
  for (j in seq_len(n_bootstrap)){
    for (i in seq_len(n_substudies)) {
      frequentist_results[i,] <- frequentist_results[i,] + frequentist_t_test(ECDF, variable, factor, subStudyIndex = i, global = FALSE, var_equal = var_equal, bootstrap = bootstrap, n_plan = n_plan)
      }
  }
  frequentist_results <- frequentist_results/n_bootstrap
  
  # add GLOBAL row as the last row, no bootstrapping here !!
  global_row <- frequentist_t_test(ECDF, variable, factor, subStudyIndex = n_substudies + 1, global = TRUE, var_equal = var_equal)
  frequentist_results[n_substudies+1,] <- global_row
  
  return(frequentist_results)
}

##############################################################################################################################################################
# Sequential analysis

full_seq_t_test_analysis <- function(ECDF, alpha, betaFutility, deltaMin, esMinFutility, variable, factor, n_permutations = 10){
  
  n_substudies <- max(ECDF$study.order)
  distinct_factors <- unique(ECDF[[factor]])
  factor1 <- distinct_factors[1]
  factor2 <- distinct_factors[2]
  
  sequential_results <- data.frame(
    study.order    = 1:n_substudies * n_permutations, # because we divide later
    stopped_times  = integer(n_substudies),
    final_times    = integer(n_substudies),
    percent_saved  = numeric(n_substudies),
    stopped_for_R  = numeric(n_substudies),
    stopped_for_N  = numeric(n_substudies),
    stopped_for_F  = numeric(n_substudies)
  )
  
  
  # to optimize speed swap the places of k and i !
  
  for (i in seq_len(n_substudies)) {
    subStudy <- ECDF[ECDF$study.order == i, ]
    sequential_results$final_times[i] <- nrow(subStudy)*n_permutations # because we divide later
  
  for (k in seq_len(n_permutations)){
    subStudy <- subStudy[sample(nrow(subStudy)),]
    x <- subStudy[subStudy[[factor]] == factor1, ][[variable]]
    y <- subStudy[subStudy[[factor]] == factor2, ][[variable]]
    n1 <- 0
    n2 <- 0
    
    for (j in seq_len(nrow(subStudy))){
      
      if(subStudy[j,factor] == factor1){
        n1 <- n1 + 1
      }
      else{
        n2 <- n2 + 1
      }
      
      if (n1 < 2 || n2 < 2){next}

      DO <- designSaviT(deltaMin=deltaMin, testType="twoSample", alternative="twoSided", futility = TRUE, wantSampling = FALSE)
      res <- saviTTest(x[1:n1], y[1:n2], designObj = DO, sequential = FALSE, varEqual = TRUE, futility=TRUE, esMinFutility = esMinFutility)
      eValue <- unname(res$eValue)
      fValue <- unname(res$eValueFut)
      if (eValue >= 1/alpha){
        sequential_results$stopped_times[i] <- sequential_results$stopped_times[i] + j
        sequential_results$stopped_for_R[i] <- sequential_results$stopped_for_R[i] + 1
        break
      }
      if (fValue <= betaFutility){
        sequential_results$stopped_times[i] <- sequential_results$stopped_times[i] + j
        sequential_results$stopped_for_F[i] <- sequential_results$stopped_for_F[i] + 1
        break
      }
      if(j == nrow(subStudy)){
        sequential_results$stopped_times[i] <- sequential_results$stopped_times[i] + j
        sequential_results$stopped_for_N[i] <- sequential_results$stopped_for_N[i] + 1
        break
      }
    }
  }
  }
  sequential_results <- sequential_results/n_permutations
  sequential_results$percent_saved <- (1-sequential_results$stopped_times/sequential_results$final_times)*100
  
  return(sequential_results)
}

