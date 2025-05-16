###### -- given a matrix-like feature representation return block ids ---------

BlockByRank <- function(index1,
                        partner1,
                        index2,
                        partner2) {
  # the feature indexing doesn't reset on new contigs, so we need to evaluate
  # each contig to contig comparison set individually and not bridge where inappropriate
  
  test_this <- mget(ls())
  if (length(unique(lengths(test_this))) != 1) {
    stop ("all vectors must be the same length")
  }
  
  block_uid <- 1L
  
  # only run block size checks if enough rows are present
  FeaturesMat <- data.frame("i1" = index1,
                            "f1" = partner1,
                            "i2" = index2,
                            "f2" = partner2)
  dr1 <- FeaturesMat[, 2L] + FeaturesMat[, 4L]
  dr2 <- FeaturesMat[, 2L] - FeaturesMat[, 4L]
  InitialBlocks1 <- unname(split(x = FeaturesMat,
                                 f = list(as.integer(FeaturesMat[, 1L]),
                                          as.integer(FeaturesMat[, 3L]),
                                          dr1),
                                 drop = TRUE))
  InitialBlocks2 <- unname(split(x = FeaturesMat,
                                 f = list(as.integer(FeaturesMat[, 1L]),
                                          as.integer(FeaturesMat[, 3L]),
                                          dr2),
                                 drop = TRUE))
  Blocks <- c(InitialBlocks1[vapply(X = InitialBlocks1,
                                    FUN = function(x) {
                                      nrow(x)
                                    },
                                    FUN.VALUE = vector(mode = "integer",
                                                       length = 1)) > 1L],
              InitialBlocks2[vapply(X = InitialBlocks2,
                                    FUN = function(x) {
                                      nrow(x)
                                    },
                                    FUN.VALUE = vector(mode = "integer",
                                                       length = 1)) > 1L])
  L01 <- length(Blocks)
  if (L01 > 0) {
    for (m3 in seq_along(Blocks)) {
      # blocks are guaranteed to contain more than 1 row
      
      sp1 <- vector(mode = "integer",
                    length = nrow(Blocks[[m3]]))
      # we need to check both columns here, this currently is not correct
      # in all cases
      sp2 <- Blocks[[m3]][, 4L]
      sp3 <- Blocks[[m3]][, 2L]
      
      it1 <- 1L
      it2 <- sp2[1L]
      it4 <- sp3[1L]
      # create a map vector on which to split the groups, if necessary
      for (m4 in seq_along(sp1)) {
        it3 <- sp2[m4]
        it5 <- sp3[m4]
        if (abs(it3 - it2 > 1L) |
            abs(it5 - it4 > 1L)) {
          # if predicted pairs are not contiguous, update the iterator
          it1 <- it1 + 1L
        }
        sp1[m4] <- it1
        it2 <- it3
        it4 <- it5
      }
      
      # if the splitting iterator was updated at all, a gap was detected
      if (it1 > 1L) {
        Blocks[[m3]] <- unname(split(x = Blocks[[m3]],
                                     f = sp1))
      } else {
        Blocks[[m3]] <- Blocks[m3]
      }
      
    } # end m3 loop
    # Blocks is now a list where each position is a set of blocked pairs
    Blocks <- unlist(Blocks,
                     recursive = FALSE)
    # drop blocks of size 1, they do not need to be evaluated
    Blocks <- Blocks[vapply(X = Blocks,
                            FUN = function(x) {
                              nrow(x)
                            },
                            FUN.VALUE = vector(mode = "integer",
                                               length = 1)) > 1]
    L01 <- length(Blocks)
    AbsBlockSize <- rep(1L,
                        nrow(FeaturesMat))
    BlockID_Map <- rep(-1L,
                       nrow(FeaturesMat))
    # only bother with this if there are blocks remaining
    # otherwise AbsBlockSize, which is initialized as a vector of 1s
    # will be left as a vector of 1s, all pairs are singleton pairs in that scenario
    if (L01 > 0L) {
      for (m3 in seq_along(Blocks)) {
        # rownames of the Blocks dfs relate to row positions in the original
        # matrix
        pos <- as.integer(rownames(Blocks[[m3]]))
        val <- rep(nrow(Blocks[[m3]]),
                   nrow(Blocks[[m3]]))
        # do not overwrite positions that are in larger blocks
        keep <- AbsBlockSize[pos] < val
        if (any(keep)) {
          AbsBlockSize[pos[keep]] <- val[keep]
          BlockID_Map[pos[keep]] <- rep(block_uid,
                                        sum(keep))
          block_uid <- block_uid + 1L
        }
      } # end m3 loop
    } # end logical check for block size
  } else {
    # no blocks observed, all pairs present are singleton pairs
    AbsBlockSize <- rep(1L,
                        nrow(FeaturesMat))
    BlockID_Map <- rep(-1L,
                       nrow(FeaturesMat))
  }
  res <- list("absblocksize" = AbsBlockSize,
              "blockidmap" = BlockID_Map)
  return(res)
}
