#' disp
#'
#' @param message     A message to be displayed in the Console.
#' @param header     Print a header of '~' symbols (=\code{TRUE}), or '~' symbols with few words of text (=\code{character vector})
#' @param footer     Print a footer '~' symbols.
#'
#' @description Displays easy-to-spot text in the Console.
#'
#' @author Fred Hasselman
#'
#' @export
#'
disp <- function(message = "Hello world!", header = FALSE, footer = FALSE) {
  ps <- "# "

  msg <- textConnection(message)
  mWidth <- max(sapply(readLines(msg), nchar))
  if (!grepl(ps, message)) {
    mWidth <- mWidth + 2
  }

  if (is.character(header)) {
    hWidth <- max(sapply(header, nchar))
    mWidth <- max(hWidth, mWidth)
  }

  dmessage <- list()
  for (m in 1:length(message)) {
    # b <- floor((mWidth-nchar(message[m]))/2)
    e <- mWidth - nchar(message[m])
    dmessage[[m]] <- paste0(ps, message[m])
  }

  banner <- paste0(rep("~", mWidth), collapse = "")
  if (is.character(header)) {
    b <- floor((nchar(banner) - nchar(header)) / 2)
    e <- ceiling((nchar(banner) - nchar(header)) / 2)
    leader <- paste0(
      "\t",
      paste0(rep("~", b), collapse = ""),
      header,
      paste0(rep("~", e), collapse = "")
    )
  }
  if (header == TRUE) {
    leader <- banner
  }
  if (header == FALSE) {
    leader <- banner
  }

  if (footer) {
    cat(paste0("\t", leader, "\n\t", dmessage, "\n\t", banner, "\n"))
  } else {
    cat(paste0("\t", leader, "\n\t", dmessage))
  }
  close(msg)
  return(invisible(message))
}

#' get.chain
#'
#' @param inf Internal
#'
#' @export
#'
get.chain <- function(inf) {
  # Build a filter chain

  filt.vars <- unlist(inf$study.cases.include, recursive = F)

  # include sites
  filt.site <- paste0(
    " %>% dplyr::filter(",
    paste0(inf$study.sites.include),
    ")"
  )
  # filt.site <- paste0(" %>% dplyr::filter(is.character(.id))")

  # Data frame filter
  filt.df <- paste0(
    " %>% dplyr::select(",
    paste0(inf$id.vars, collapse = ","),
    ")",
    filt.site
  )

  # Variables filter
  return(list(df = filt.df, vars = filt.vars))
}

#' get.sourceData
#'
#' @param ML2.id  Internal
#' @param ML2.df  Internal
#' @param ML2.in  Inernal
#'
#' @return A list with fields \code{study.vars} (data organised according to the \code{masteRkey} spreadsheet), \code{study.vars/labels}, \code{N}, and \code{RawDataFilter}(raw data, unfiltered).
#' @export
#'

get.sourceData <- function(ML2.id, ML2.df, ML2.in) {
  N <- numeric(length(ML2.in$study.vars))
  if (is.null(ML2.df$uID)) {
    ML2.df$uID <- seq_along(ML2.df$source)
  }
  # study.vars[1]    <- unlist(ML2.in$study.vars)
  vars <- list()
  dfname <- list()
  RawData <- list()
  # id <- factor(ML2.df$.id)
  for (i in seq_along(ML2.in$study.vars)) {
    dfname[i] <- names(ML2.in$study.vars)[[i]]
    eval(parse(
      text = paste0(
        names(ML2.in$study.vars)[[i]],
        " <- tibble::as_tibble(dplyr::select(ML2.df,",
        paste0(c(ML2.in$study.vars[[i]], "uID"), collapse = ","),
        "))"
      )
    ))

    # Change text to numbers
    suppressWarnings(
      if (
        any(eval(parse(
          text = paste0(
            "apply(",
            names(ML2.in$study.vars)[i],
            ",2,is.numeric)==FALSE"
          )
        )))
      ) {
        eval(parse(
          text = paste0(
            names(ML2.in$study.vars)[i],
            " <- data.frame(sapply(which(apply(",
            names(ML2.in$study.vars)[i],
            ", 2, is.numeric)==FALSE), function(si) as.numeric(",
            names(ML2.in$study.vars)[i],
            "[[si]])))"
          )
        ))
      }
    )
  }

  for (i in seq_along(dfname)) {
    eval(parse(text = paste0(dfname[i], " <- ", unlist(ML2.id$vars)[i])))
    eval(parse(
      text = paste0(
        "attr(",
        names(ML2.in$study.vars)[[i]],
        ",'uID')",
        " <- ",
        paste0(names(ML2.in$study.vars)[[i]], "[['uID']]")
      )
    ))
    RawData[[i]] <- dplyr::mutate(
      ML2.df,
      Included = eval(parse(
        text = paste0("ML2.df$uID %in% ", dfname[i], "$uID")
      ))
    )
    #    N[i]         <- eval(parse(text=paste0("nrow(",dfname[i],")"))
  }

  if (ML2.in$stat.params$within) {
    sameID <- ML2.df$uID[RawData[[1]]$Included & RawData[[2]]$Included]
  }

  #  sameID <- eval(parse(text=paste0(dfname[which.max(N)],"$uID %in% ", dfname[which.min(N)],"$uID") ))
  for (i in seq_along(dfname)) {
    # Check if datasets are equal length for within subject analyses
    if (ML2.in$stat.params$within) {
      eval(parse(
        text = paste0(
          dfname[i],
          " <- ",
          dfname[i],
          "[",
          dfname[i],
          "$uID %in% sameID, ]"
        )
      ))
      RawData[[i]] <- dplyr::mutate(
        ML2.df,
        Included = eval(parse(
          text = paste0(
            "ML2.df$uID %in% ",
            dfname[i],
            "$uID",
            collapse = "&"
          )
        ))
      )
    }
    N[i] <- eval(parse(
      text = paste0("sum(!is.na(", dfname[i], "), na.rm = TRUE)")
    ))
    # eval(parse(text=paste0(dfname[i],"<-", dfname[i]," %>% dplyr::select(which(colnames(",dfname[i],")!='uID'))")))
    eval(parse(text = paste0(dfname[i], " -> vars[[i]]")))
  }

  vars[[length(ML2.in$study.vars) + 1]] <- N
  # if(length(ML2.in$study.vars.labels)==0){ML2.in$study.vars.labels <- list(NoLabels="None")}
  vars[[length(ML2.in$study.vars) + 2]] <- ML2.in$study.vars.labels
  vars[[length(ML2.in$study.vars) + 3]] <- RawData
  names(vars) <- c(names(ML2.in$study.vars), "N", "labels", "RawDataFilter")
  return(vars)
}

#' get.info
#'
#' @param keytable Internal
#' @param cols Internal
#'

#' @export
#'
get.info <- function(keytable, cols, subset) {
  # Read Variables and Parameters from:
  # keytable <- ML2.key[s,]
  study.vars <- eval(parse(text = keytable[, "study.vars"]))
  study.vars.labels <- eval(parse(text = keytable[, "study.vars.labels"]))
  cases.include <- eval(parse(text = keytable[, "study.cases.include"]))
  stat.params <- eval(parse(text = keytable[, "stat.params"]))
  # cases <- plyr::llply(seq_along(cases.include),function(i) get.cases(cases.include[i],study.vars,study.vars.labels,stat.params))
  cases <- get.cases(
    cases.include,
    study.vars,
    study.vars.labels,
    stat.params
  )
  sites.include <- eval(parse(text = keytable[, "study.sites.include"]))
  if (sites.include[[1]][1] == "all") {
    sites.include[[1]] <- "is.character(source)"
  }
  if (subset != "all") {
    W <- ifelse(subset == "WEIRD", 1, 0)
    sites.include[[1]] <- paste0(
      sites.include[[1]][1],
      " & (Weird == ",
      W,
      ")"
    )
  }

  # Find correct columns in this dataset according to ML2.key: 'ML2.in$study.vars'
  id.vars <- which(
    cols %in%
      c(
        unlist(study.vars),
        "uID",
        ".id",
        "age",
        "sex",
        "source",
        "Source.Global",
        "Source.Primary",
        "Source.Secondary",
        "Country",
        "Location",
        "Language",
        "Weird",
        "SubjectPool",
        "Setting",
        "Tablet",
        "Pencil",
        "Execution",
        "StudyOrderN",
        "IDiffOrderN"
      )
  )
  return(list(
    study.vars = study.vars,
    study.vars.labels = study.vars.labels,
    stat.params = stat.params,
    study.cases.include = cases,
    study.sites.include = sites.include,
    id.vars = id.vars
  ))
}


#' decide.analysis
#'
#' @param ML2.key Key
#' @param studies studies
#' @param tp analysis type
#' @param doAll use all studies
#'
#' @return list of studies
#' @export
#'
decide.analysis <- function(ML2.key, studies = NA, tp = NA, doAll = FALSE) {
  analysis <- c(
    "study.global.include",
    "study.primary.include",
    "study.secondary.include",
    "study.global.include"
  )
  groups <- c(
    "Source.Global",
    "Source.Primary",
    "Source.Secondary",
    "study.order"
  )

  if (is.null(tp)) {
    tp <- NA
  }
  if (any(is.na(studies))) {
    studies <- na.exclude(ML2.key$unique.id)
  }
  if (tp == 4) {
    if (doAll) {
      studies <- studies[studies %in% ML2.key$unique.id]
    } else {
      studies <- studies[
        studies %in%
          ML2.key$unique.id[ML2.key$study.figure2.include == 1]
      ]
    }
  }

  if (is.na(tp[1])) {
    disp(
      paste0(
        "Analyzing global, primary and secondary studies marked '1' in the corrsponding columns of the masteRkey spreadsheet"
      ),
      header = FALSE,
      footer = TRUE
    )
    studiess <- c(
      ML2.key$unique.id[ML2.key$study.global.include[studies] == 1],
      ML2.key$unique.id[ML2.key$study.primary.include[studies] == 1],
      ML2.key$unique.id[ML2.key$study.secondary.include[studies] == 1]
    )
    tps <- c(
      rep(
        1,
        length(ML2.key$unique.id[
          ML2.key$study.global.include[studies] == 1
        ])
      ),
      rep(
        2,
        length(ML2.key$unique.id[
          ML2.key$study.primary.include[studies] == 1
        ])
      ),
      rep(
        3,
        length(ML2.key$unique.id[
          ML2.key$study.secondary.include[studies] == 1
        ])
      )
    )
  } else {
    disp(
      paste0("Analyzing studies in ", analysis[tp]),
      header = FALSE,
      footer = TRUE
    )
    if (doAll) {
      studiess <- lapply(analysis[tp], function(c) {
        ML2.key$unique.id[ML2.key[, c] >= 0]
      })
    } else {
      studiess <- lapply(analysis[tp], function(c) {
        ML2.key$unique.id[ML2.key[, c] == 1]
      })
    }
    tps <- unlist(sapply(seq_along(studiess), function(s) {
      rep(tp[s], length(studiess[[s]]))
    }))
    studiess <- unlist(studiess)
  }

  tp <- tps[!is.na(studiess) & studiess %in% studies]

  if (any(is.na(studiess[!is.na(studiess) & studiess %in% studies]))) {
    stop(
      "Analysis ID and analysis type do not agree [e.g. analysis type is 'primary', but analysis ID refers to 'secondary']"
    )
  }

  return(list(
    studiess = studiess[!is.na(studiess) & studiess %in% studies],
    ugroup = groups[tp],
    tp = tp
  ))
}

#' get.cases
#'
#' @param rule Internal
#' @param study.vars Internal
#' @param study.vars.labels Internal
#' @param stat.params Internal
#'
#' @export
#'
get.cases <- function(rule, study.vars, study.vars.labels, stat.params) {
  # rule <- cases.include
  type <- names(rule)
  if (!is.matrix(rule[[1]])) {
    rule <- rbind(rule[[1]])
  } else {
    rule <- rule[[1]]
  }
  Nrule <- nrow(rule)

  isna <- ifelse(
    stat.params$censorNA,
    {
      " & !is.na(X)"
    },
    {
      ""
    }
  )
  do <- list()

  # Rule 1 should always be about the variables in 'study.vars'
  # Rule 2..N can concern a subset of the variables in 'study.vars'.
  # Assumption for type="each": Rule 2 number of variables is a multiple of number of conditions, with variables grouped to fit to each condition.
  # Assumption for type="sep": Rule 2 contains a separate rule independent of variables listed in the first rule.

  r <- 1
  filtvars <- study.vars.labels[names(study.vars.labels) %in% rule[r, 1]]

  # Start with 'pre'
  pre <- plyr::llply(unlist(filtvars), function(v) {
    paste0(v, " %>% dplyr::filter(")
  })

  if (type == "each") {
    if (all(filtvars[[1]] %in% names(study.vars))) {
      filtvars <- study.vars
      do <- plyr::laply(
        seq_along(filtvars),
        function(v) {
          plyr::laply(
            seq_along(filtvars[[v]]),
            function(vv) {
              paste0(
                gsub("X", filtvars[[v]][vv], rule[r, 2]),
                gsub("X", filtvars[[v]][vv], isna)
              )
            }
          )
        }
      )
    }

    if (!is.matrix(do)) {
      do <- as.matrix(do)
    }

    if (Nrule > 1) {
      s <- ncol(do)
      for (r in 2:Nrule) {
        filtvars <- unlist(study.vars.labels[
          names(study.vars.labels) %in% rule[r, 1]
        ])
        do <- cbind(
          do,
          plyr::laply(
            seq_along(filtvars),
            function(vv) {
              paste0(gsub("X", filtvars[[vv]], rule[r, 2]))
            }
          )
        )
      }
    }
    case <- plyr::llply(seq_along(pre), function(fr) {
      paste0(pre[[fr]], paste0(do[fr, ], collapse = " & "), ")")
    })
    names(case) <- names(study.vars)
  }

  if (type == "sep") {
    if (all(filtvars[[1]] %in% names(study.vars))) {
      filtvars <- study.vars
      do <- plyr::llply(
        seq_along(filtvars),
        function(v) {
          plyr::laply(
            seq_along(filtvars[[v]]),
            function(vv) {
              paste0(
                gsub("X", filtvars[[v]][vv], rule[r, 2]),
                gsub("X", filtvars[[v]][vv], isna)
              )
            }
          )
        }
      )
    }

    names(do) <- names(filtvars)

    if (Nrule > 1) {
      for (r in 2:Nrule) {
        if (rule[r, 1] %in% names(do)) {
          filtvars <- unlist(study.vars.labels[
            names(study.vars.labels) %in% rule[r, 1]
          ])
          do[[r]] <- c(
            do[[r]],
            plyr::laply(
              seq_along(filtvars),
              function(vv) {
                paste0(gsub("X", filtvars[[vv]], rule[r, 2]))
              }
            )
          )
        }
      }
    }

    case <- plyr::llply(seq_along(pre), function(fr) {
      paste0(pre[[fr]], paste0(do[[fr]], collapse = " & "), ")")
    })
    names(case) <- names(study.vars)
  }

  return(eval(parse(text = paste0("list(", rule[1], "=case)"))))
}

#' check.nMin
#'
#' @param df The dataframe containing the data.
#' @param group The current group to check.
#' @param Nmin.raw Minimum raw sample size.
#' @param Nmin.cond Minimum sample size per condition.
#' @param id The ML2 id list.
#' @param info The ML2 info list.
#'
#' @return A list containing the source data if the checks pass, otherwise NULL.
#' @export
#'
check.nMin <- function(df, group, Nmin.raw, Nmin.cond, id, info) {
  gID <- df$source %in% group

  if (sum(gID, na.rm = TRUE) < Nmin.raw) {
    # print(paste0("Removing site ", group, " (not enough raw data)"))
    return(NULL)
  }

  ML2.sr_g <- get.sourceData(id, df[gID, ], info)

  compN <- ML2.sr_g$N
  compN1 <- sum(ML2.sr_g$RawDataFilter[[1]]$Included, na.rm = TRUE)
  compN2 <- sum(ML2.sr_g$RawDataFilter[[2]]$Included, na.rm = TRUE)

  if (!(any(compN >= Nmin.raw) & all(compN1 >= Nmin.cond, compN2 >= Nmin.cond))) {
    # print(paste0("Removing site ", group, " (not enough data in one condition)"))
    return(NULL)
  }

  return(ML2.sr_g)
}

#' varfun.Hauser.1
#'
#' @param vars    A list object generated by \code{\link{get.sourceData}} containing cleaned data and variable labels.
#'
#' @return Dataset ready for analysis
#'
#' @export
#'
#' @description \url{https://manylabsopenscience.github.io/ML2_PoPS_proposal#11_hauser}
#'
#' @section Variables:
#'   haus1.1t = timing (side effect scenario);
#'   haus2.1t = timing (greater good scenario);
#'   haus1.2=previous experience (drop if 1 (yes));
#'   haus2.2=previous experience (drop if 1 (yes));
#'   haus1.1=morally permissible (side effect scenario; Yes=1);
#'   haus2.1=morally permissible (greater good scenario; Yes=1).
#'
#'
varfun.Hauser.1 <- function(vars) {
  SE <- unlist(vars$SideEffect[names(vars$SideEffect)[[1]]])
  GG <- unlist(vars$GreaterGood[names(vars$GreaterGood)[[1]]])
  N <- c(nrow(vars$SideEffect), nrow(vars$GreaterGood))

  cleanDataFilter <- data.frame(
    uID = c(vars[[1]]$uID, vars[[2]]$uID),
    variable = factor(
      c(SE, GG),
      levels = c(1, 2),
      labels = vars$labels$Response
    ),
    factor = factor(
      c(rep(1, N[1]), rep(2, N[2])),
      levels = c(1, 2),
      labels = vars$labels$Condition
    )
  )

  return(list(
    Response = factor(
      c(SE, GG),
      levels = c(1, 2),
      labels = vars$labels$Response
    ),
    Condition = factor(
      c(rep(1, N[1]), rep(2, N[2])),
      levels = c(1, 2),
      labels = vars$labels$Condition
    ),
    N = N,
    cleanDataFilter = cleanDataFilter
  ))
}



#' varfun.Tversky.1
#'
#' @param vars     A list object generated by \code{\link{get.sourceData}} containing cleaned data and variable labels.
#'
#' @return Dataset ready for analysis
#'
#' @references    Tversky, A., Kahneman, D. (1981). The framing of decisions and the psychology of choice. Science, 211, 453-458.
#'
#' @description \url{https://manylabsopenscience.github.io/ML2_PoPS_proposal#16_tversky}
#'
#' @section Variables:
#' tver1.1=choice ($250 wall hanging condition, yes=1, no=2);
#' tver2.1=choice ($30 wall hanging cond, yes=1, no=2).
#'
#' tver1.1 Imagine that you are about to purchase a ceramic vase for $30, and a wall hanging for $250. The salesman informs you that the wall hanging you wish to buy is on sale for $240 at the other branch of the store, located 20 minutes drive away. Would you make the trip to the other store?
#' m  Yes, I would go to the other branch. (1)
#' m  No, I would not go to the other branch. (2)
#'
#' tver2.1 Imagine that you are about to purchase a ceramic vase for $250, and a wall hanging for $30. The salesman informs you that the wall hanging you wish to buy is on sale for $20 at the other branch of the store, located 20 minutes drive away. Would you make the trip to the other store?
#' m  Yes, I would go to the other branch. (1)
#' m  No, I would not go to the other branch. (2)
#'
#' @export
#'
varfun.Tversky.1 <- function(vars) {
  cleanDataFilter <- data.frame(
    uID = c(vars$Cheap$uID, vars$Expensive$uID),
    variable = factor(c(vars$Cheap[[1]], vars$Expensive[[1]]), levels = c(1, 2), labels = vars$labels$Response),
    factor = factor(c(rep(1, nrow(vars$Cheap)), rep(2, nrow(vars$Expensive))), levels = c(1, 2), labels = vars$labels$Condition)
  )

  return(list(
    Condition = factor(c(rep(1, nrow(vars$Cheap)), rep(2, nrow(vars$Expensive))), levels = c(1, 2), labels = vars$labels$Condition),
    Response = factor(c(vars$Cheap[[1]], vars$Expensive[[1]]), levels = c(1, 2), labels = vars$labels$Response),
    N = c(nrow(vars$Cheap), nrow(vars$Expensive)),
    cleanDataFilter = cleanDataFilter
  ))
}


#' varfun.Rottenstreich.1
#'
#'
#' @description \url{https://manylabsopenscience.github.io/ML2_PoPS_proposal#5_rottenstreich}
#'
#' @param vars     A list object generated by \code{\link{get.sourceData}} containing cleaned data and variable labels.
#'
#'  @export
#'
#' @return Dataset ready for analysis
#'
varfun.Rottenstreich.1 <- function(vars) {
  cleanDataFilter <- data.frame(
    uID = c(vars[[1]]$uID, vars[[2]]$uID),
    variable = factor(c(vars$Low[[1]], vars$Certain[[1]]), levels = c(1, 2), labels = vars$labels$Response),
    factor = factor(c(rep(1, nrow(vars$Low)), rep(2, nrow(vars$Certain))), levels = c(1, 2), labels = vars$labels$Condition)
  )

  return(list(
    Response = factor(c(vars$Low[[1]], vars$Certain[[1]]), levels = c(1, 2), labels = vars$labels$Response),
    Condition = factor(c(rep(1, nrow(vars$Low)), rep(2, nrow(vars$Certain))), levels = c(1, 2), labels = vars$labels$Condition),
    N = c(nrow(vars$Certain), nrow(vars$Low)),
    cleanDataFilter = cleanDataFilter
  ))
}
