###### -- hit rearrangement by bounds -----------------------------------------

AAHitScoping <- function(hitlist,
                         fstrand1,
                         fstart1,
                         fstop1,
                         fstrand2,
                         fstart2,
                         fstop2) {
  # true is the negative strand
  w1 <- as.logical(fstrand1)
  # true is the negative strand
  w2 <- as.logical(fstrand2)
  
  # we'll need this a few times, just set it
  w3 <- length(hitlist)
  
  # create a big matrix and a column map
  hitmat <- do.call(cbind,
                    hitlist)
  colmap <- rep(seq(w3),
                times = vapply(X = hitlist,
                               FUN = function(x) {
                                 ncol(x)
                               },
                               FUN.VALUE = vector(mode = "integer",
                                                  length = 1L)))
  # re-frame feature bounds
  fstart1 <- fstart1[colmap]
  fstart2 <- fstart2[colmap]
  fstop1 <- fstop1[colmap]
  fstop2 <- fstop2[colmap]
  
  # negative is now FALSE
  # now framed against the hit matrix, as opposed to the hit list
  w1_adj <- !w1[colmap]
  # negative is now FALSE
  # now framed against the hit matrix, as opposed to the hit list
  w2_adj <- !w2[colmap]
  
  # return(list("a" = hitmat,
  #             "b" = w1_adj,
  #             "c" = w2_adj,
  #             "d" = fstart1,
  #             "e" = fstart2,
  #             "f" = fstop1,
  #             "g" = fstop2,
  #             "h" = colmap,
  #             "i" = hitlist))
  
  res1 <- res2 <- res3 <- res4 <- vector(mode = "integer",
                                         length = ncol(hitmat))
  # forward strand first feature
  # get the starts offset them, and then transform
  res1[w1_adj] <- fstart1[w1_adj] + (hitmat[1, w1_adj, drop = TRUE] - 1L) * 3L
  res2[w1_adj] <- fstart1[w1_adj] + (hitmat[2, w1_adj, drop = TRUE] - 1L) * 3L
  # reverse strand first feature
  # get the offsets from the stop, and then transform
  res2[!w1_adj] <- fstop1[!w1_adj] - (hitmat[1, !w1_adj, drop = TRUE] - 1L) * 3L
  res1[!w1_adj] <- fstop1[!w1_adj] - (hitmat[2, !w1_adj, drop = TRUE] - 1L) * 3L
  
  # forward strand second feature
  res3[w2_adj] <- fstart2[w2_adj] + (hitmat[3, w2_adj, drop = TRUE] - 1L) * 3L
  res4[w2_adj] <- fstart2[w2_adj] + (hitmat[4, w2_adj, drop = TRUE] - 1L) * 3L
  
  # reverse strand second feature
  res4[!w2_adj] <- fstop2[!w2_adj] - (hitmat[3, !w2_adj, drop = TRUE] - 1L) * 3L
  res3[!w2_adj] <- fstop2[!w2_adj] - (hitmat[4, !w2_adj, drop = TRUE] - 1L) * 3L
  
  
  
  # remake the matrix
  result <- cbind(res1,
                  res2,
                  res3,
                  res4)
  
  # split recycles down the VECTOR, so we cbind above
  result <- split(x = result,
                  f = colmap)
  
  # then we transform into the appropriate shape for alignpairs
  result <- unname(lapply(X = result,
                          FUN = function(x) {
                            matrix(data = x,
                                   nrow = 4,
                                   byrow = TRUE)
                          }))
  
  return(result)
  
}
