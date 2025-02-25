###### -- Select the true-est ortholog ----------------------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu
# hunger games of predicted pairs

CompetePairs <- function(SynExtendObject,
                         AllowCrossContigConflicts = TRUE, # if EITHER assembly isn't 'complete' this should be FALSE
                         By = "PID",
                         PollContext = TRUE,
                         NormalizeCompetition = TRUE,
                         InflationParameter = .975,
                         Verbose = FALSE) {
  # start with timing
  if (Verbose) {
    pBar <- txtProgressBar(style = 1)
    FunctionTimeStart <- Sys.time()
  }
  # unit normalize if asked for...
  NormScore <- function(vec) {
    vec <- vec / sqrt(sum(vec^2))
    return(vec)
  }
  # overhead checking
  if (!is(object = SynExtendObject,
          class2 = "PairSummaries")) {
    stop ("SynExtendObject must be an object of class 'PairSummaries'.")
  }
  # i need to build a class file for this class, until then we need to do this:
  attr_vals <- attributes(SynExtendObject)
  attr_names <- names(attr_vals)
  
  # if i allow polling of context, i need to make sure that 
  if (PollContext) {
    if (!("Block_UID" %in% colnames(SynExtendObject))) {
      stop ("Block unique IDs are missing")
    }
  }
  if (length(By) != 1) {
    stop ("Competition between conflicting pairs can only currently be resolved with a single colname.")
  }
  # check the compete by argument
  if (!(By %in% colnames(SynExtendObject))) {
    stop ("Competition between conflicting pairs can only be resolved with an available colname.")
  }
  if (PollContext & NormalizeCompetition) {
    # bind a context score to easily access later
    ContextScore <- tapply(X = NormScore(SynExtendObject[, By]),
                           INDEX = SynExtendObject$Block_UID,
                           FUN = function(x) {
                             sum(x) / (length(x) / InflationParameter)
                           })
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = ContextScore[match(x = SynExtendObject$Block_UID,
                                                              table = as.integer(names(ContextScore)))])
  } else if (PollContext & !NormalizeCompetition) {
    # bind a context score to easily access later
    ContextScore <- tapply(X = SynExtendObject[, By],
                           INDEX = SynExtendObject$Block_UID,
                           FUN = function(x) {
                             sum(x) / (length(x) / InflationParameter)
                           })
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = ContextScore[match(x = SynExtendObject$Block_UID,
                                                              table = as.integer(names(ContextScore)))])
  } else if (!PollContext & NormalizeCompetition) {
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = NormScore(SynExtendObject[, By]))
  } else if (!PollContext & !NormalizeCompetition) {
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = SynExtendObject[, By])
  }
  res01 <- vector(mode = "character",
                  length = nrow(SynExtendObject))
  res02 <- rep(TRUE,
               nrow(SynExtendObject))
  
  # divy object up into comparison sections to be evaluated
  if (AllowCrossContigConflicts) {
    # split the object based on genome_contig_
    row_reps <- paste(gsub(pattern = "_[0-9]+$",
                           replacement = "",
                           x = SynExtendObject$p1),
                      gsub(pattern = "_[0-9]+$",
                           replacement = "",
                           x = SynExtendObject$p2),
                      sep = "_")
    row_key <- unique(row_reps)
    # return(list("rep" = row_reps,
    #             "key" = row_key,
    #             "origin" = SynExtendObject))
  } else {
    # split the object based on genome_
    row_reps <- paste(gsub(pattern = "_[0-9+]_[0-9]+$",
                           replacement = "",
                           x = SynExtendObject$p1),
                      gsub(pattern = "_[0-9]+_[0-9]+$",
                           replacement = "",
                           x = SynExtendObject$p2),
                      sep = "_")
    row_key <- unique(row_reps)
    # return(list("rep" = row_reps,
    #             "key" = row_key,
    #             "origin" = SynExtendObject))
  }
  eval_sets <- split(x = SynExtendObject,
                     f = row_reps)
  eval_sets <- lapply(X = eval_sets,
                      FUN = function(x) {
                        class(x) <- c("data.frame",
                                      "PairSummaries")
                        return(x)
                      })
  # return(eval_sets)
  if (Verbose) {
    PBAR <- length(eval_sets)
  }
  # loop through our evaluation sets, return some kind of vector or other object
  # that tells me who to knock out, and *who it was knocked out by*
  for (m1 in seq_along(eval_sets)) {
    current_sets <- DisjointSet(Pairs = eval_sets[[m1]],
                                Verbose = FALSE)
    
    communities01 <- rep(as.integer(names(current_sets)),
                         lengths(current_sets))
    names(communities01) <- unlist(current_sets)
    communities02 <- unname(communities01[match(x = eval_sets[[m1]]$p1,
                                                table = names(communities01))])
    # return(list("set" = current_sets,
    #             "communities01" = communities01,
    #             "communities02" = communities02))
    pairs_by_community <- split(x = eval_sets[[m1]],
                                f = communities02)
    # row names here are relative to ... 
    conflict_pairs <- pairs_by_community[vapply(X = pairs_by_community,
                                                FUN = function(x) {
                                                  nrow(x) > 1
                                                },
                                                FUN.VALUE = vector(mode = "logical",
                                                                   length = 1))]
    # return(list("communities1" = pairs_by_community,
    #             "communities2" = conflict_pairs))
    if (length(conflict_pairs) > 0) {
      conflict_rows <- unlist(unname(lapply(X = conflict_pairs,
                                            FUN = function(x) {
                                              rownames(x)
                                            })))
      keep <- vector(mode = "list",
                     length = length(conflict_pairs))
      winner <- vector(mode = "character",
                       length = length(conflict_pairs))
      # return a list of vectors of logicals and a vector of characters that indicates 
      # WHO won the competition
      # we can propogate these back to the original data frame based on the rownames
      for (m2 in seq_along(conflict_pairs)) {
        # return(list("current" = eval_sets[[m1]],
        #             "conflict" = conflict_pairs[[m2]]))
        winner[m2] <- rownames(conflict_pairs[[m2]])[which.max(conflict_pairs[[m2]]$CompeteBy)]
        keep[[m2]] <- rownames(conflict_pairs[[m2]]) %in% winner[m2]
        # if (m1 == 2 &
        #     m2 == 1) {
        #   return(list("origin" = SynExtendObject,
        #               "current" = eval_sets[[m1]],
        #               "conflicts" = conflict_pairs[[m2]],
        #               "winner" = winner[m2],
        #               "keep" = keep[[m2]]))
        # }
      }
      # return(list("origin" = SynExtendObject,
      #             "current" = eval_sets[[m1]],
      #             "conflicts" = conflict_pairs,
      #             "row_key" = conflict_rows,
      #             "winner" = winner,
      #             "keep" = keep))
      winner <- rep(winner,
                    lengths(keep))
      # fill in which row knocked you out if you were knocked out
      res01[match(x = conflict_rows,
                  table = rownames(SynExtendObject))] <- winner
      # fill in whether you were knocked out
      res02[match(x = conflict_rows,
                  table = rownames(SynExtendObject))] <- unlist(keep)
      # post conflict, see if any of the removed pairs have BOTH PARTNERS
      # missing in the winners set, if so add them back
      # this needs to be rewritten as some kind of competition itself
      # check_conflicts <- do.call(rbind,
      #                            conflict_pairs)
      # retained_p1 <- check_conflicts$p1[unlist(keep)]
      # retained_p2 <- check_conflicts$p2[unlist(keep)]
      # 
      # rescue_p1 <- !(check_conflicts$p1 %in% retained_p1)
      # rescue_p2 <- !(check_conflicts$p2 %in% retained_p2)
      # 
      # if (any(rescue_p1 & rescue_p2)) {
      #   # return(list("a" = SynExtendObject,
      #   #             "b" = eval_sets[[m1]],
      #   #             "c" = check_conflicts[, c(1:2,11:12,16)],
      #   #             "d" = unlist(keep),
      #   #             "e" = list(retained_p1,
      #   #                        retained_p2),
      #   #             "f" = list(rescue_p1,
      #   #                        rescue_p2)))
      #   res02[match(x = conflict_rows,
      #               table = rownames(SynExtendObject))][which(rescue_p1 & rescue_p2)] <- TRUE
      # }
      
    } # end check for checkable conflicts
    
    if (Verbose) {
      setTxtProgressBar(pb = pBar,
                        value = m1 / PBAR)
    }
  }
  if (Verbose) {
    close(pBar)
  }
  
  # until the pair summaries class has a class file i need to do this...
  attrs_to_add <- attr_names[!(attr_names) %in% names(attributes(SynExtendObject))]
  # return(list("a" = attr_names,
  #             "b" = attr_vals,
  #             "c" = attrs_to_add,
  #             "d" = SynExtendObject))
  for (m1 in seq_along(attrs_to_add)) {
    attr(x = SynExtendObject,
         which = attrs_to_add[m1]) <- attr_vals[[attrs_to_add[m1]]]
  }
  attr(x = SynExtendObject,
       which = "RetainByCompetition") <- res02
  attr(x = SynExtendObject,
       which = "Knockout") <- res01
  class(SynExtendObject) <- c("data.frame", "PairSummaries")
  
  if (Verbose) {
    FunctionTimeEnd <- Sys.time()
    print(FunctionTimeEnd - FunctionTimeStart)
  }
  
  return(SynExtendObject)
}


