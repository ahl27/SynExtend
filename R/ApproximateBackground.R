###### -- background score generation -----------------------------------------
# calculate the background scores for a given set of paired features

ApproximateBackground <- function(p1,
                                  p2,
                                  code1,
                                  code2,
                                  mod1,
                                  mod2,
                                  aa1,
                                  aa2,
                                  nt1,
                                  nt2,
                                  register1,
                                  register2,
                                  aamat,
                                  ntmat) {
  # minimal overhead checking as this function is intended for internal use in
  # SummarizePairs
  
  test_these <- mget(ls())
  if (length(unique(lengths(test_these[c("p1",
                                         "p2",
                                         "code1",
                                         "code2",
                                         "mod1",
                                         "mod2")]))) > 1L) {
    stop ("some input vectors appear to not be the correct length")
  }
  
  l01 <- length(p1)
  res <- vector(mode = "numeric",
                length = l01)
  # all four qualifiers must be true to ask for AA backgrounds
  nt_aa_switch <- code1 & code2 & mod1 & mod2
  for (m1 in seq_len(l01)) {
    if (nt_aa_switch[m1]) {
      # coding
      res[m1] <- sum(aa1[register1[p1[m1]], ] * t(aa2[register2[p2[m1]], ] * aamat))
    } else {
      # non-coding
      res[m1] <- sum(nt1[p1[m1], ] * t(nt2[p2[m1], ] * ntmat))
    }
  }
  return(res)
}

