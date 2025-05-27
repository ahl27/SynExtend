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
  
  # suggestion from Aidan, avoids unique()
  if(!all((y <- lengths(test_this)) == y[1])){
    stop("all vectors must be the same length")
  }
  
  # previous check
  # if (length(unique(lengths(test_this))) > 1) {
  #   stop ("all vectors must be the same length")
  # }
  
  feat1w <- gene1right - gene1left + 1L
  feat2w <- gene2right - gene2left + 1L
  feat1lb <- (hit1left - gene1left) / feat1w
  feat1rb <- (gene1right - hit1right) / feat1w
  feat2lb <- (hit2left - gene2left) / feat2w
  feat2rb <- (gene2right - hit2right) / feat2w
  
  res <- .Call("HitConsensus",
               feat1lb, feat1rb, feat2lb, feat2rb, as.logical(strand1), as.logical(strand2),
               PACKAGE = "SynExtend")
  
  return(res)
}
