ExoLabel <- function(edgelistfiles, outfile=tempfile(),
                          mode=c("undirected", "directed"),
                          add_self_loops=FALSE,
                          attenuation_power=1,
                          ignore_weights=FALSE,
                          iterations=0L,
                          return_table=FALSE,
                          consensus_cluster=FALSE,
                          use_fast_sort=FALSE,
                          verbose=interactive(),
                          sep='\t',
                          tempfiledir=tempdir()){
  if(!is.numeric(iterations)){
    stop("'iterations' must be an integer or numeric.")
  } else {
    iterations <- as.integer(iterations)
  }
  if(iterations > 2^15){
    warning("'iterations' currently only supports signed 16-bit numbers, defaulting to max possible value of 32767.")
    iterations <- 32767L
  }
  if(is.na(iterations) || is.null(iterations) || is.infinite(iterations) || iterations < 0){
    warning("Invalid value of 'iterations', will determine automatically from node degree.")
    iterations <- 0L
  }
  if(!is.numeric(add_self_loops) && !is.logical(add_self_loops)){
    stop("value of 'add_self_loops' should be numeric or logical")
  }
  if(any(add_self_loops < 0)){
    warning("self loops weight supplied is negative, setting to zero.")
    add_self_loops[add_self_loops < 0] <- 0
  } else if(is.logical(add_self_loops)){
    add_self_loops <- rep(ifelse(add_self_loops, 1, 0), length(outfile))
  }
  if(!is.logical(ignore_weights)){
    stop("'ignore_weights' must be logical")
  } else if(is.na(ignore_weights) || is.null(ignore_weights)){
    stop("invalid value for 'ignore_weights' (should be TRUE or FALSE)")
  }
  if(ignore_weights && (add_self_loops != 0 || add_self_loops != 1)){
    warning("Weight specified for 'add_self_loops' will be ignored")
    add_self_loops <- 1
  }
  if(!is.logical(use_fast_sort)){
    stop("invalid value for 'use_fast_sort' (should be TRUE or FALSE)")
  }
  if(!is.numeric(attenuation_power)){
    stop("'attenuation_power' must be a valid numeric value")
  } else if (is.integer(attenuation_power)){
    attenuation_power <- as.numeric(attenuation_power)
  }
  if(is.na(attenuation_power) || is.null(attenuation_power) || is.infinite(attenuation_power)){
    stop("'attenuation_power' must be a valid numeric value")
  }
  if(length(add_self_loops) == 1){
    add_self_loops <- rep(add_self_loops, length(outfile))
  }
  if(length(attenuation_power) == 1){
    attenuation_power <- rep(attenuation_power, length(outfile))
  }
  if(length(outfile) != length(add_self_loops) ||
     length(attenuation_power) != length(add_self_loops)){
      stop("If more than one outfile is provided, 'add_self_loops' and ",
          "'attenuation_power' must be either the same length as 'outfile' ",
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
                       attenuation_power)
  names(graph_stats) <- c("num_vertices", "num_edges")
  for(f in list.files(tempfiledir, full.names=TRUE))
    if(file.exists(f)) file.remove(f)
  file.remove(tempfiledir)
  retval <- list()
  for(i in seq_along(outfile)){
    param_vec <- c(add_self_loops=add_self_loops[i],
                   attenuation_power=attenuation_power[i])
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
  return(retval)
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
