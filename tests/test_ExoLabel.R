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

  file_fxns <- list(
    none=\(x) {},
    gz=\(x){ system(paste("gzip", "-f", x))}
  )
  file_endings <- c("", ".gz")

  for(i in seq_along(file_fxns)){
    cat("File Compression:", names(file_fxns)[i], '\n')

    to_process <- paste0(c(tf1, tf2), file_endings[i])
    elf1 <- to_process[1]
    elf2 <- to_process[2]

    cat("\tSmall graphs:...")
    for(loop in c(0, 0.25, 0.5)){
      df <- generate_random_graph(10, 25)
      if(any(abs(df$w - loop) < WEIGHT_TOLERANCE))
        df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] <- df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] + 3*WEIGHT_TOLERANCE
      write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
      file_fxns[[i]](tf1)
      testExo(elf1, add_self_loops=loop)
    }
    cat("passed.\n")

    cat("\tLarger graphs:...")
    for(loop in c(0, 0.5)){
      df <- generate_random_graph(10000, 25000)
      if(any(abs(df$w - loop) < WEIGHT_TOLERANCE))
        df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] <- df$w[abs(df$w - loop) < WEIGHT_TOLERANCE] + 3*WEIGHT_TOLERANCE
      write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
      file_fxns[[i]](tf1)
      testExo(elf1, add_self_loops=loop)
    }
    cat("passed.\n")


    cat("\tDirected Edges...")
    testExo(elf1, mode="directed")
    cat("passed.\n")

    cat("\tNo fast sort...")
    testExo(elf1, use_fast_sort=FALSE)
    cat("passed.\n")

    ## I'll just use the same graph here
    cat("\tDifferent separator...")
    write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    file_fxns[[i]](tf1)
    testExo(elf1, sep=',')
    cat("passed.\n")

    cat("\tMulti-file input...")
    df <- generate_random_graph(50000, 100000)
    write.table(df[1:50000,], tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    write.table(df[50000:100000,], tf2, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    file_fxns[[i]](tf1)
    file_fxns[[i]](tf2)
    testExo(c(elf1, elf2))
    cat("passed.\n")

    cat("\tHeaders...")
    testExo(c(elf1, elf2), header=TRUE)
    testExo(c(elf1, elf2), header=10L)
    cat("passed.\n")

    cat("\tLarger weights...")
    df[,3] <- df[,3] * 1000
    write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
    file_fxns[[i]](tf1)
    testExo(elf1)
    cat("passed.\n")

    file.remove(elf1)
    file.remove(elf2)
  }

  cat("\nAll checks passed!\n")
  invisible(TRUE)
}

run_status_tests()
