###### -- consensus fun -------------------------------------------------------
# given the bounds of hits and the bounds of the features they live in
# calculate the similarity of positioning of the hits
# hits are *not* required to be the same size in the features they're linking
# and i don't know if the math works out correctly if they're differing sizes

HitConsensus <- function(gene1left,
                         gene2left,
                         gene1right,
                         gene2right,
                         strand1,
                         strand2,
                         hit1left,
                         hit1right,
                         hit2left,
                         hit2right) {
  # minimal overhead checking as the function is built specifically for use internal
  # to SummarizePairs()
  
  # by default this will look at all the objects passed to the function's environment
  # objects not passed to the function's environment will not be captures by ls() in this instance
  test_this <- mget(ls())
  
  if (length(unique(lengths(test_this))) > 1) {
    stop ("all vectors must be the same length")
  }
  
  feat1w <- gene1right - gene1left + 1L
  feat2w <- gene2right - gene2left + 1L
  feat1lb <- (hit1left - gene1left) / feat1w
  feat1rb <- (gene1right - hit1right) / feat1w
  feat2lb <- (hit2left - gene2left) / feat2w
  feat2rb <- (gene2right - hit2right) / feat2w
  
  # return(list("a" = feat1lb,
  #             "b" = feat2lb,
  #             "c" = feat1rb,
  #             "d" = feat2rb,
  #             ))
  l01 <- length(gene1left)
  res <- vector(mode = "numeric",
                length = l01)
  for (m1 in seq_along(gene1left)) {
    
    # for both positive OR both negative
    if (strand1[m1] == strand2[m1]) {
      # mean of the absolute values of the differences between the distances from the bounds
      # when both are in the same direction it's just direct left to right comparisons
      # i.e. is the linking hit's left bound in gene 1 a similar distance from gene 1's left bound
      # as the linking hit's left bound in gene 2 is from gene 2's left bound
      res[m1] <- mean(c(abs(feat1lb[m1] - feat2lb[m1]),
                        abs(feat1rb[m1] - feat2rb[m1])))
      # when one is negative and the other is positive, the comparison is just flipped
      # i.e. is the linking hit's left bound in gene 1 a similar distance from gene 1's left bound
      # as the linking hit's RIGHT bound from gene 2's RIGHT bound
      # and vice versa
    } else {
      # this flip results in the same comparison regardless of which feature is flipped
      res[m1] <- mean(c(abs(feat1lb[m1] - feat2rb[m1]),
                        abs(feat1rb[m1] - feat2lb[m1])))
    } 
  } # end of m1 loop
  return(res)
}
