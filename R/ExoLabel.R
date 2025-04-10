.safecheck_optional_numeric <- function(vname, value, rep_to=1L){
  if(is.logical(value)){
    value <- as.numeric(value)
  }
  if(!is.numeric(value)){
    stop("'", vname, "' must be a valid numeric value")
  } else if (is.integer(value)){
    value <- as.numeric(value)
  }
  if(any(is.na(value) | is.null(value) | is.infinite(value))){
    stop("'", vname, "' must be a valid numeric value")
  }

  if(rep_to > 1 && length(value) == 1L)
    value <- rep(value, rep_to)
  value
}

.safecheck_optional_integer <- function(vname, value, max=-1,
                                        can_be_invalid=FALSE, rep_to=1L){
  if(!is.numeric(value)){
    stop("'", vname, "' must be integer or numeric.")
  } else {
    value <- as.integer(value)
  }
  if(max > 0 && any(value > max)){
    warning("Coercing invalid values for '", vname, "' to their maximum allowed value of ", max, ".")
    value[value > max] <- max
  }
  if(any(is.na(value) | is.null(value) | is.infinite(value) | value < 0)){
    if(can_be_invalid){
      warning("Invalid value of '", vname, "', will determine automatically.")
      value[is.na(value) | is.null(value) | is.infinite(value) | value < 0] <- 0L
    } else {
      stop("'", vname, "' must be a valid integer value")
    }
  }

  if(rep_to > 1 && length(value) == 1L)
    value <- rep(value, rep_to)
  value
}

ExoLabel <- function(edgelistfiles,
                          outfile=tempfile(tmpdir=tempfiledir),
                          mode=c("undirected", "directed"),
                          add_self_loops=FALSE,
                          attenuation=TRUE,
                          dist_scaling=TRUE,
                          ignore_weights=FALSE,
                          iterations=0L,
                          return_table=FALSE,
                          consensus_cluster=FALSE,
                          use_fast_sort=TRUE,
                          verbose=interactive(),
                          sep='\t',
                          tempfiledir=tempdir()){
  if(return_table){
    maxp <- max(length(add_self_loops),
                length(attenuation),
                length(iterations),
                length(dist_scaling))
    if(!missing(outfile)){
      warning("'outfile' will be ignored since return_table=TRUE")
    }
    if(length(outfile) != maxp){
      outfile <- replicate(maxp, tempfile(tmpdir=tempfiledir))
    }
  }
  iterations <- .safecheck_optional_integer("iterations", iterations, 32767L, TRUE, length(outfile))
  attenuation <- .safecheck_optional_numeric("attenuation", attenuation, length(outfile))
  dist_scaling <- .safecheck_optional_numeric("dist_scaling", dist_scaling, length(outfile))
  add_self_loops <- .safecheck_optional_numeric("add_self_loops", add_self_loops, length(outfile))

  if(any(add_self_loops < 0)){
    warning("self loops weight supplied is negative, setting to zero.")
    add_self_loops[add_self_loops < 0] <- 0
  }
  if(!is.logical(ignore_weights)){
    stop("'ignore_weights' must be logical")
  } else if(is.na(ignore_weights) || is.null(ignore_weights)){
    stop("invalid value for 'ignore_weights' (should be TRUE or FALSE)")
  }
  if(ignore_weights && any(add_self_loops != 0 & add_self_loops != 1)){
    warning("Weight specified for 'add_self_loops' will be ignored")
    add_self_loops <- as.numeric(as.logical(add_self_loops))
  }
  if(!is.logical(use_fast_sort)){
    stop("invalid value for 'use_fast_sort' (should be TRUE or FALSE)")
  }

  if(length(add_self_loops) != length(outfile) ||
     length(attenuation) != length(outfile) ||
     length(iterations) != length(outfile)){
      stop("If more than one outfile is provided, 'add_self_loops', 'iterations'",
          ", 'attenuation', and 'dist_scaling' must be either the same length as 'outfile' ",
          "or length 1.")
  }
  if(!is.logical(verbose) || !is.numeric(verbose)){
    if(as.integer(verbose) != verbose){
      warning("Coercing non-integer argument provided to 'verbose'")
    }
    verbose_int <- 0L
    if(is.numeric(verbose)){
      ## more for internal use
      verbose_int <- as.integer(verbose)
    } else if(verbose){
      verbose_int <- ifelse(interactive(), 2L, 1L)
    }
  }
  # verify that the first few lines of each file are correct
  if(!all(file.exists(edgelistfiles))) stop("edgelist file does not exist")
  edgelistfiles <- normalizePath(edgelistfiles, mustWork=TRUE)
  for(f in edgelistfiles){
    v <- readLines(f, n=10L)
    v <- strsplit(v, sep)
    lv <- lengths(v)
    if(any(lv != lv[1]) || lv[1] < 2) stop("file ", f, " is misformatted")
    lv <- lv[1] # now we know they're all the same
    if(!ignore_weights && lv == 2) stop("file ", f, " is missing weights!")
    if(!ignore_weights && any(vapply(v, \(x) is.na(as.numeric(x[3])), logical(1L))))
      stop("file ", f, " has malformed weights")
  }

  if(is.logical(consensus_cluster)){
    if(consensus_cluster){
      consensus_cluster <- c(0,0.2,0.4,0.6,0.8,1,1.33,1.67,2)
    } else {
      consensus_cluster <- numeric(0L)
    }
  } else {
    if(!is.numeric(consensus_cluster) || any(is.na(consensus_cluster) | is.null(consensus_cluster)))
      stop("'consensus_cluster' must be a logical or numeric vector")
    if(any(consensus_cluster < 0))
      stop("'consensus_cluster' cannot contain negative values")
  }
  tempfiledir <- normalizePath(tempfiledir, mustWork=TRUE)
  tempfiledir <- file.path(tempfiledir, "ExoLabelTemp")
  if(dir.exists(tempfiledir)){
    for(f in list.files(tempfiledir, full.names=TRUE))
      file.remove(f)
  } else {
    dir.create(tempfiledir, recursive = TRUE)
  }

  for(i in seq_along(outfile)){
    d <- dirname(outfile[i])
    if(!dir.exists(d)){
      dir.create(d, recursive=TRUE)
    }
    d <- normalizePath(d, mustWork = TRUE)
    outfile[i] <- file.path(d, basename(outfile[i]))
  }

  ## have to declare this here so we can also clean up the temporary directory
  on.exit(\(){
    .C("cleanup_ondisklp_global_values")
    for(f in list.files(tempfiledir, full.names=TRUE))
      file.remove(f)
    file.remove(tempfiledir)
  })

  mode <- match.arg(mode)
  is_undirected <- mode == "undirected"
  outfile <- file.path(normalizePath(dirname(outfile), mustWork=TRUE), basename(outfile))

  if(verbose_int > 0L) cat("Temporary files stored at ", tempfiledir, "\n")

  seps <- paste(sep, "\n", sep='')

  graph_stats <- .Call("R_LPOOM_cluster",
                       edgelistfiles, length(edgelistfiles),
                       tempfiledir, outfile, seps, iterations,
                       verbose_int, is_undirected,
                       add_self_loops, ignore_weights,
                       consensus_cluster, !use_fast_sort,
                       attenuation, dist_scaling, FALSE)
  names(graph_stats) <- c("num_vertices", "num_edges")
  for(f in list.files(tempfiledir, full.names=TRUE))
    if(file.exists(f)) file.remove(f)
  file.remove(tempfiledir)
  retval <- list()
  for(i in seq_along(outfile)){
    param_vec <- c(add_self_loops=add_self_loops[i],
                   attenuation=attenuation[i],
                   iterations=iterations[i],
                   dist_scaling=dist_scaling[i])
    if(return_table){
      tab <- read.table(outfile[i], sep=sep)
      colnames(tab) <- c("Vertex", "Cluster")
      if(file.exists(outfile[i])) file.remove(outfile[i])
      retval[[i]] <- list(parameters=param_vec,
                          graph_stats=graph_stats,
                          results=tab)
    } else {
      retval[[i]] <- list(parameters=param_vec,
                          graph_stats=graph_stats,
                          results=outfile[i])
    }
  }

  if(length(retval) == 1) return(retval[[1]])
  invisible(retval)
}

EstimateExoLabel <- function(num_v, avg_degree=2, is_undirected=TRUE,
                          num_edges=num_v*avg_degree, node_name_length=10L){
  if(!missing(avg_degree) && !missing(num_edges)){
    warning("Only one of 'avg_degree' and 'num_edges' are needed, ignoring num_edges")
  } else if (missing(avg_degree)){
    avg_degree <- (num_edges*ifelse(is_undirected, 2, 1)) / num_v
  }
  lv <- num_v*node_name_length

  # assuming file is v1 v2 %.3f, which is 2*node_name_len + 3 + 5
  FRCS <- 8192*4
  exp_size_file <- (2*node_name_length+8)*(num_edges)
  exp_size_internal_inplace <- 16*num_edges*ifelse(is_undirected, 2, 1)
  exp_size_internal <- 32*num_edges*ifelse(is_undirected, 2, 1)
  # rough guess at trie size + edge reading buffer
  exp_size_ram <- 24*node_name_length + 18 * num_v * 2 + FRCS*16 + 40960
  if(num_edges > FRCS){
    # internal buffers for mergesorting, not subtracting 1 because we need an
    # additional buffer for copying data around in the in-place merge
    exp_size_ram <- exp_size_ram + FRCS*(64+16)*16
  }
  exp_size_final <- (2+node_name_length+log10(num_v))*num_v
  exp_ratio <- exp_size_internal_inplace / exp_size_file
  v <- c(exp_size_ram, exp_size_file, exp_size_internal_inplace, exp_size_internal, exp_size_final, exp_ratio)
  names(v) <- c("RAM Upper Bound Estimate", "Expected Input File Size", "Expected Internal File Size (SlowSort)",
                "Expected Internal File Size (FastSort)", "Expected Final File Size", "Disk Usage Ratio")

  max_nchar <- max(nchar(names(v)[-length(v)]))
  unitsizes <- c("B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB")
  for(i in seq_along(v)){
    if(i == length(v)){
      if(exp_ratio < 0.001){
        cat("\nExoLabel total disk consumption is <0.001x that of the original files\n")
      } else {
        cat("\nExoLabel total disk consumption is about ", round(exp_ratio, 2), "x that of the initial files.\n", sep='')
      }
      next # skip other stats for disk usage ratio
    }
    unit <- ""
    for(j in seq_along(unitsizes)){
      unit <- unitsizes[j]
      p <- 1024^(j-1)
      if((v[i] / p) < 1024)
        break
    }
    n <- names(v)[i]
    padding_required <- max_nchar - nchar(n)
    pad <- paste(rep(' ', padding_required), collapse='')
    cat(pad, names(v)[i], ": ", sprintf("%5.1f ", v[i]/p), unit, '\n', sep='')
  }
  invisible(v)
}

.testExoLabel <- function(edgelistfiles,
                     outfile=tempfile(tmpdir=tempfiledir),
                     mode=c("undirected", "directed"),
                     add_self_loops=FALSE,
                     attenuation=TRUE,
                     dist_scaling=TRUE,
                     ignore_weights=FALSE,
                     iterations=0L,
                     consensus_cluster=FALSE,
                     use_fast_sort=TRUE,
                     verbose=FALSE,
                     sep='\t',
                     tempfiledir=tempdir(),
                     skip_checks=FALSE){
  verbose_int <- as.integer(verbose)
  maxp <- max(length(add_self_loops),
              length(attenuation),
              length(iterations),
              length(dist_scaling))
  if(length(outfile) != maxp)
    outfile <- replicate(maxp, tempfile(tmpdir=tempfiledir))
  iterations <- .safecheck_optional_integer("iterations", iterations, 32767L, TRUE, length(outfile))
  attenuation <- .safecheck_optional_numeric("attenuation", attenuation, length(outfile))
  dist_scaling <- .safecheck_optional_numeric("dist_scaling", dist_scaling, length(outfile))
  add_self_loops <- .safecheck_optional_numeric("add_self_loops", add_self_loops, length(outfile))

  # verify that the first few lines of each file are correct
  edgelistfiles <- normalizePath(edgelistfiles, mustWork=TRUE)
  for(f in edgelistfiles){
    v <- readLines(f, n=10L)
    v <- strsplit(v, sep)
    lv <- lengths(v)
    if(any(lv != lv[1]) || lv[1] < 2) stop("file ", f, " is misformatted")
    lv <- lv[1] # now we know they're all the same
    if(!ignore_weights && lv == 2) stop("file ", f, " is missing weights!")
    if(!ignore_weights && any(vapply(v, \(x) is.na(as.numeric(x[3])), logical(1L))))
      stop("file ", f, " has malformed weights")
  }

  if(is.logical(consensus_cluster)){
    if(consensus_cluster){
      consensus_cluster <- c(0,0.2,0.4,0.6,0.8,1,1.33,1.67,2)
    } else {
      consensus_cluster <- numeric(0L)
    }
  }
  tempfiledir <- normalizePath(tempfiledir, mustWork=TRUE)
  tempfiledir <- file.path(tempfiledir, "ExoLabelTemp")
  if(dir.exists(tempfiledir)){
    for(f in list.files(tempfiledir, full.names=TRUE))
      file.remove(f)
  } else {
    dir.create(tempfiledir, recursive = TRUE)
  }

  for(i in seq_along(outfile)){
    d <- dirname(outfile[i])
    if(!dir.exists(d)){
      dir.create(d, recursive=TRUE)
    }
    d <- normalizePath(d, mustWork = TRUE)
    outfile[i] <- file.path(d, basename(outfile[i]))
  }

  ## have to declare this here so we can also clean up the temporary directory
  on.exit(\(){
    .C("cleanup_ondisklp_global_values")
    for(f in list.files(tempfiledir, full.names=TRUE))
      file.remove(f)
    file.remove(tempfiledir)
  })

  mode <- match.arg(mode)
  is_undirected <- mode == "undirected"
  outfile <- file.path(normalizePath(dirname(outfile), mustWork=TRUE), basename(outfile))

  if(verbose_int > 0L) cat("Temporary files stored at ", tempfiledir, "\n")

  seps <- paste(sep, "\n", sep='')
  if(verbose){
    cat('\n')
    cat("********************\n")
    cat("* Internal Checks: *\n")
    cat("********************\n")
  }

  result <- .Call("R_LPOOM_cluster",
                       edgelistfiles, length(edgelistfiles),
                       tempfiledir, outfile, seps, iterations,
                       verbose_int, is_undirected,
                       add_self_loops, ignore_weights,
                       consensus_cluster, !use_fast_sort,
                       attenuation, dist_scaling, TRUE)

  if(verbose){
    cat('\n')
    cat("*****************\n")
    cat("* Final Checks: *\n")
    cat("*****************\n")
  }
  graph_stats <- result[[1]]
  all_weights <- result[[2]]
  all_degrees <- result[[3]]
  disjoint_sizes <- result[[4]]

  benchmark_graph <- data.frame(v1=character(0L), v2=character(0L), w=numeric(0L))
  for(f in edgelistfiles){
    tmp <- read.delim(f, header=FALSE, sep=sep)
    benchmark_graph <- rbind(benchmark_graph, tmp)
  }
  colnames(benchmark_graph) <- c("v1", "v2", "weight")
  if(skip_checks){
    cat("SKIPPED!\n")
    names(result) <- c("graph_stats", "all_weights", "all_degrees", "disjoint_sizes")
    result$graph <- benchmark_graph
    return(result)
  }

  names(graph_stats) <- c("num_vertices", "num_edges")
  for(f in list.files(tempfiledir, full.names=TRUE))
    if(file.exists(f)) file.remove(f)
  file.remove(tempfiledir)
  retval <- list()
  cluster_res <- read.table(outfile, sep=sep)
  colnames(cluster_res) <- c("Vertex", "Cluster")

  ## check that results line up
  num_verts_actual <- length(unique(c(benchmark_graph$v1, benchmark_graph$v2)))
  if(num_verts_actual != graph_stats[1])
    stop("Found ", graph_stats[1], " vertices, ",
          "actual count was ", num_verts_actual)
  if(verbose) cat("Correct number of nodes read in.\n")

  ## Checking that node degrees in-degrees are correct
  all_node_labels <- unique(c(benchmark_graph$v1, benchmark_graph$v2))
  actual_degrees <- numeric(length(all_node_labels))
  names(actual_degrees) <- all_node_labels
  if(mode=="undirected"){
    tmp <- table(c(benchmark_graph$v1, benchmark_graph$v2))
  } else {
    tmp <- table(benchmark_graph$v2)
  }
  actual_degrees[names(tmp)] <- tmp
  actual_degrees <- sort(actual_degrees)
  all_degrees <- sort(all_degrees)
  if(!all(actual_degrees == all_degrees))
    stop("Not all node degrees were correct!")
  if(verbose) cat("All node degrees are correct.\n")


  ## Checking that weights are correct
  close_w <- benchmark_graph$w
  close_w <- sort(close_w[close_w > 0])
  all_weights <- sort(all_weights)
  if(mode=='undirected')
    all_weights <- all_weights[seq(1,length(all_weights),2)]
  if(length(close_w) != length(all_weights))
    stop("Found ", length(all_weights), " edges, ",
         "actual count was ", length(close_w))
  TOLERANCE <- 0.001 * close_w # both 0.01% and 0.0001 absolute tolerance
  TOLERANCE[TOLERANCE == 0] <- Inf
  TOLERANCE[TOLERANCE < 0.0001] <- 0.0001
  all_weights <- abs(all_weights - close_w)
  if(any(all_weights > TOLERANCE))
    stop("Misread some edges: maximum difference in read weight was ", max(all_weights))
  if(verbose) cat("All weights were read in correctly.\n")

  ## checking that disjoint sets are correct
  sub_disjoint_validate <- benchmark_graph[benchmark_graph$weight > add_self_loops[1],]
  v1 <- sub_disjoint_validate[,1]
  v2 <- sub_disjoint_validate[,2]
  all_v <- unique(c(v1, v2))
  set_sizes <- FindSets(match(v1, all_v), match(v2, all_v))[,2]
  set_sizes <- match(set_sizes, unique(set_sizes))
  set_sizes <- sort(tabulate(set_sizes))
  disjoint_sizes <- disjoint_sizes[disjoint_sizes > 1]
  if(!all(set_sizes == sort(disjoint_sizes)))
    stop("ExoLabel recorded a different graph! (Disjoint set sizes differ)")
  if(verbose) cat("All disjoint sets are the same size at cutoff ", add_self_loops[1], ".\n", sep='')


  if(nrow(cluster_res) != num_verts_actual)
    stop("ExoLabel reported ", nrow(cluster_res), " clusters, but there were ",
          num_verts_actual, " nodes!")
  if(verbose) cat("Correct number of clusters reported.\n")
  file.remove(outfile)
  if(!verbose) cat("All checks passed.\n")
  return(TRUE)
}
