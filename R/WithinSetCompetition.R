###### -- Select the true-est equivalog ---------------------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu
# hunger games of predicted pairs

WithinSetCompetition <- function(SynExtendObject,
                                 AllowCrossContigConflicts = TRUE,
                                 CompeteBy  = "Delta_Background",
                                 PollContext = TRUE,
                                 ContextInflation = 0.975,
                                 Verbose = FALSE) {
  # start with timing
  if (Verbose) {
    pBar <- txtProgressBar(style = 1)
    FunctionTimeStart <- Sys.time()
  }
  # overhead checking
  if (!is(object = SynExtendObject,
          class2 = "PairSummaries")) {
    stop ("SynExtendObject must be an object of class 'PairSummaries'.")
  }
  # i need to build a class file for this class, until then we need to do this:
  attr_vals <- attributes(SynExtendObject)
  attr_names <- names(attr_vals)
  pair_key <- paste(SynExtendObject[, "p1"],
                    SynExtendObject[, "p2"],
                    sep = "_")
  
  # right now the only context to poll when competing is block uid
  # eventually this should expand to candidate inference source, but for now
  # this is our only context methodology
  if (PollContext) {
    if (!("Block_UID" %in% colnames(SynExtendObject))) {
      stop ("Block unique IDs are missing")
    }
  }
  if (length(CompeteBy) != 1) {
    stop ("Competition between conflicting pairs can only currently be resolved with a single colname.")
  }
  # check the compete by argument
  if (!(CompeteBy %in% colnames(SynExtendObject))) {
    stop ("Competition between conflicting pairs can only be resolved with an available colname.")
  }
  
  if (PollContext) {
    adjustby <- tapply(X = SynExtendObject,
                       INDEX = SynExtendObject[, "Block_UID"],
                       FUN = function(x) {
                         # keep the rownames around for reference later
                         y1 <- x[, CompeteBy, drop = FALSE]
                         y2 <- rownames(y1)
                         y3 <- y1[, 1L] + mean(y1[, 1L])
                         names(y3) <- y2
                         return(y3)
                       },
                       simplify = FALSE)
    adjustby <- unlist(unname(adjustby))
    # return(adjustby)
    adjustby <- adjustby[match(table = names(adjustby),
                               x = rownames(SynExtendObject))]
    # return(adjustby)
    competevector <- unname(adjustby)
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = competevector)
    # SynExtendObject$CompeteBy <- competevector
    # return(competevector)
  } else {
    competevector <- SynExtendObject[, CompeteBy]
    # SynExtendObject$CompeteBy <- competevector
    SynExtendObject <- cbind(SynExtendObject,
                             "CompeteBy" = competevector)
  }
  
  res01 <- vector(mode = "character",
                  length = nrow(SynExtendObject))
  res02 <- rep(TRUE,
               nrow(SynExtendObject))
  
  indexmat <- do.call(rbind,
                      strsplit(x = paste(SynExtendObject$p1,
                                         SynExtendObject$p2,
                                         sep = "_"),
                               split = "_",
                               fixed = TRUE))
  # return(list(indexmat,
  #             pair_key))
  # indexmat <- matrix(data = as.integer(indexmat),
  #                    nrow = nrow(indexmat))
  indexdf <- data.frame("g1" = indexmat[, 1],
                        "i1" = indexmat[, 2],
                        "f1" = indexmat[, 3],
                        "g2" = indexmat[, 4],
                        "i2" = indexmat[, 5],
                        "f2" = indexmat[, 6])
  
  key1 <- apply(X = indexdf[, c(1,2,4,5)],
                MARGIN = 1,
                FUN = function(x) {
                  paste(x,
                        collapse = "_")
                })
  key2 <- apply(X = indexdf[, c(1,4)],
                MARGIN = 1,
                FUN = function(x) {
                  paste(x,
                        collapse = "_")
                })
  
  if (AllowCrossContigConflicts) {
    df1 <- split(x = SynExtendObject,
                 f = key1)
  } else {
    df1 <- split(x = SynExtendObject,
                 f = key2)
  }
  
  df1 <- lapply(X = df1,
                FUN = function(x) {
                  class(x) <- c("data.frame",
                                "PairSummaries")
                  return(x)
                })
  for (m1 in seq_along(df1)) {
    current_sets <- DisjointSet(Pairs = df1[[m1]],
                                Verbose = FALSE)
    communities01 <- rep(as.integer(names(current_sets)),
                         lengths(current_sets))
    names(communities01) <- unlist(current_sets)
    communities02 <- unname(communities01[match(x = df1[[m1]]$p1,
                                                table = names(communities01))])
    # return(list("set" = current_sets,
    #             "communities01" = communities01,
    #             "communities02" = communities02,
    #             "df" = df1[[m1]]))
    # print(length(communities01))
    # print(length(communities02))
    pairs_by_community <- split(x = df1[[m1]],
                                f = communities02)
    # print(length(pairs_by_community))
    # row names here are relative to ... 
    conflict_pairs <- pairs_by_community[vapply(X = pairs_by_community,
                                                FUN = function(x) {
                                                  nrow(x) > 1
                                                },
                                                FUN.VALUE = vector(mode = "logical",
                                                                   length = 1))]
    # print(length(conflict_pairs))
    if (length(conflict_pairs) > 0) {
      conflict_rows <- unlist(unname(lapply(X = conflict_pairs,
                                            FUN = function(x) {
                                              rownames(x)
                                            })))
      keep <- vector(mode = "list",
                     length = length(conflict_pairs))
      winner <- vector(mode = "character",
                       length = length(conflict_pairs))
      # print(length(keep))
      # print(length(winner))
      # return a list of vectors of logicals and a vector of characters that indicates 
      # WHO won the competition
      # we can propogate these back to the original data frame based on the rownames
      for (m2 in seq_along(conflict_pairs)) {
        
        winner[m2] <- rownames(conflict_pairs[[m2]])[which.max(conflict_pairs[[m2]]$CompeteBy)]
        keep[[m2]] <- rownames(conflict_pairs[[m2]]) %in% winner[m2]
        
      }
      winner <- rep(winner,
                    lengths(keep))
      # print(length(keep))
      # print(length(winner))
      # print(sum(lengths(winner)))
      # fill in which row knocked you out if you were knocked out
      res01[match(x = conflict_rows,
                  table = rownames(df1[[m1]]))] <- winner
      # fill in whether you were knocked out
      res02[match(x = conflict_rows,
                  table = rownames(df1[[m1]]))] <- unlist(keep)
      # print("a")
      # print(length(res01))
      # print(length(res02))
    } else {
      # end check for checkable conflicts
      res02 <- rep(TRUE,
                   nrow(df1[[m1]]))
      # print(length(res02))
    }
    # return(list("df" = df1[[m1]],
    #             "r1" = res01,
    #             "r2" = res02))
    df1[[m1]] <- df1[[m1]][res02[match(x = rownames(df1[[m1]]),
                                       table = rownames(SynExtendObject))], ]
    
  } # end of m1 loop
  
  # return(df1)
  df1 <- do.call(rbind,
                 df1)
  rownames(df1) <- NULL
  # return(df1)
  indexmat <- do.call(rbind,
                      strsplit(x = paste(df1$p1,
                                         df1$p2,
                                         sep = "_"),
                               split = "_",
                               fixed = TRUE))
  # return(indexmat)
  # indexmat <- matrix(data = as.integer(indexmat),
  #                    nrow = nrow(indexmat))
  indexdf <- data.frame("g1" = as.integer(indexmat[, 1]),
                        "i1" = as.integer(indexmat[, 2]),
                        "f1" = as.integer(indexmat[, 3]),
                        "g2" = as.integer(indexmat[, 4]),
                        "i2" = as.integer(indexmat[, 5]),
                        "f2" = as.integer(indexmat[, 6]))
  key1 <- apply(X = indexdf[, c(1,2,4,5)],
                MARGIN = 1,
                FUN = function(x) {
                  paste(x,
                        collapse = "_")
                })
  
  df1 <- split(x = df1,
               f = key1)
  df2 <- split(x = indexdf,
               f = key1)
  df1 <- lapply(X = df1,
                FUN = function(x) {
                  class(x) <- c("data.frame",
                                "PairSummaries")
                  return(x)
                })
  
  ph <- vector(mode = "list",
               length = length(df2))
  block_offset <- 0L
  # print("a")
  for (m1 in seq_along(df1)) {
    if (nrow(df2[[m1]]) > 1) {
      # return(df2[[m1]])
      blockres <- BlockByRank(index1 = df2[[m1]]$i1,
                              partner1 = df2[[m1]]$f1,
                              index2 = df2[[m1]]$i2,
                              partner2 = df2[[m1]]$f2)
      block_ph <- blockres$blockidmap
      w1 <- block_ph > 0
      if (any(w1)) {
        block_ph[w1] <- block_ph[w1] + block_offset
        block_offset <- block_offset + max(block_ph)
      }
      # print(block_offset)
    } else {
      # no blocks to add to the offset!
      blockres <- list("absblocksize" = 1L,
                       "blockidmap" = -1L)
      block_ph <- blockres$blockidmap
      # ph[[m1]] <- blockres
    }
    # return(list(df1[[m1]],
    #             df2[[m1]],
    #             blockres))
    
    df1[[m1]]$Block_UID <- block_ph
    df1[[m1]]$blocksize <- blockres$absblocksize
    # return(df1[[m1]])
  }
  
  df1 <- do.call(rbind,
                 df1)
  rownames(df1) <- NULL
  # return(df1)
  w1 <- df1$Block_UID == -1L
  if (any(w1) &
      !all(w1)) {
    df1$Block_UID[w1] <- seq(from = max(df1$Block_UID) + 1L,
                             by = 1,
                             length.out = sum(w1))
  } else if (any(w1) &
             all(w1)) {
    df1$Block_UID <- seq(from = 1,
                         by = 1,
                         length.out = sum(w1))
  }
  
  # return(df1)
  # until the pair summaries class has a class file i need to do this...
  # attrs_to_add <- attr_names[!(attr_names) %in% names(attributes(df1))]
  # # return(list("a" = attr_names,
  # #             "b" = attr_vals,
  # #             "c" = attrs_to_add,
  # #             "d" = SynExtendObject))
  # for (m1 in seq_along(attrs_to_add)) {
  #   attr(x = df1,
  #        which = attrs_to_add[m1]) <- attr_vals[[attrs_to_add[m1]]]
  # }
  
  class(df1) <- c("data.frame",
                  "PairSummaries")
  return(df1)
}

