generate_random_graph <- function(nverts, nedges){
  require(igraph, quietly=TRUE)
  alph <- AA_STANDARD
  num_required <- ceiling(log(nverts, length(alph)))
  num_required <- max(num_required, 3)
  sample_names <- mkAllStrings(alph, num_required)
  labs <- sample(sample_names, nverts)
  g <- sample_gnm(nverts, nedges)
  df <- as_data_frame(g, what="edges")
  data.frame(v1=labs[df[,1]], v2=labs[df[,2]], w=runif(nedges))
}

run_status_tests <- function(){
  if(!require(igraph)){
    cat("Skipping tests, igraph is not available.\n")
    invisible(TRUE)
  }
  require(SynExtend)
  tf1 <- tempfile()
  tf2 <- tempfile()
  WEIGHT_TOLERANCE <- 0.001
  testExo <- SynExtend:::.testExoLabel
  cat("Small graphs:...")
  for(loop in c(0, 0.25, 0.5)){
    df <- generate_random_graph(10, 25)
    if(any(abs(df$w - loop) < WEIGHT_TOLERANCE))
      df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] <- df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] + 3*WEIGHT_TOLERANCE
    write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    testExo(tf1, add_self_loops=loop)
  }
  cat("passed.\n")

  cat("Larger graphs:...")
  for(loop in c(0, 0.5)){
    df <- generate_random_graph(10000, 25000)
    if(any(abs(df$w - loop) < WEIGHT_TOLERANCE))
      df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] <- df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] + 3*WEIGHT_TOLERANCE
    write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    testExo(tf1, add_self_loops=loop)
  }
  cat("passed.\n")


  cat("Directed Edges...")
  testExo(tf1, mode="directed")
  cat("passed.\n")

  cat("No fast sort...")
  testExo(tf1, use_fast_sort=FALSE)
  cat("passed.\n")

  ## I'll just use the same graph here
  cat("Different separator...")
  write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
  testExo(tf1, sep=',')
  cat("passed.\n")

  cat("Multi-file input...")
  tf2 <- tempfile()
  df <- generate_random_graph(50000, 100000)
  write.table(df[1:50000,], tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  write.table(df[50000:100000,], tf2, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  testExo(c(tf1, tf2))
  cat("passed.\n")

  cat("Headers...")
  testExo(c(tf1, tf2), header=TRUE)
  testExo(c(tf1, tf2), header=10L)
  cat("passed.\n")


  cat("Larger weights...")
  df[,3] <- df[,3] * 1000
  write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  testExo(tf1)
  cat("passed.\n")

  file.remove(tf1)
  file.remove(tf2)



  cat("\nAll checks passed!\n")
  invisible(TRUE)
}

run_status_tests()
