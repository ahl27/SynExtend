SuperTree <- function(myDendList, NAMEFUN=NULL, Verbose=TRUE, ...){
  # error checking
  if (!is(myDendList, 'list') ||
      !all(vapply(myDendList, \(x) is(x, 'dendrogram') || is.null(x), FUN.VALUE=TRUE))){
    stop("SuperTree requires input to be a list of dendrograms")
  }

  useNFUN <- is(NAMEFUN, 'function')

  if (useNFUN){
    l1 <- labels(myDendList[[1L]])
    l2 <- NAMEFUN(l1)
    if(length(l1)!=length(l2) || !all(is.atomic(l2))){
      stop("NAMEFUN should operate on a character vector and return a character
           vector of the same size as input.")
    }
  }

  # Get list of species
  allspecies <- character(0)
  if (Verbose){
    cat("  Compiling species...\n")
    pb <- txtProgressBar(max=length(myDendList), style=3)
    start <- Sys.time()
  }
  ctr <- 0
  for ( i in seq_along(myDendList) ){
    if (Verbose) setTxtProgressBar(pb, i)
    lst <- myDendList[[i]]
    if (is.null(lst)) next
    specs <- labels(lst)

    if (useNFUN)
      specs <- NAMEFUN(specs)

    allspecies <- unique(c(allspecies, specs))
    ctr <- ctr + 1
  }

  # Initialize distance matrix and counts matrix
  #dmat <- countmat <- matrix(0, nrow=length(allspecies), ncol=length(allspecies))
  #rownames(dmat) <- colnames(dmat) <- allspecies
  #rownames(countmat) <- colnames(countmat) <- allspecies
  asl <- length(allspecies)
  dholder <- numeric(asl*(asl-1)/2)
  class(dholder) <- "dist"
  attr(dholder, "Size") <- asl
  attr(dholder, "Diag") <- TRUE
  attr(dholder, "Upper") <- TRUE
  attr(dholder, "Labels") <- allspecies
  dmat <- countmat <- dholder

  if (Verbose){
    cat("\n  Done.\n\n  Constructing species-level distance matrix...\n")
    pb <- txtProgressBar(max = ctr, style=3)
  }
  ctr <- 0
  for ( dend in myDendList ){
    if (is.null(dend)) next
    # Calculate Cophenetic
    cp <- Cophenetic(dend)
    if(useNFUN){
      attr(cp, 'Labels') <- NAMEFUN(attr(cp, 'Labels'))
    }
    dmat <- combineDist(dmat, cp)
    cp[] <- 1
    countmat <- combineDist(countmat, cp)

    ctr <- ctr + 1
    if (Verbose) setTxtProgressBar(pb, ctr)
  }

  # Impute missing entries
  if(Verbose) cat("\n  Done.\n")
  posmissing <- which(countmat==0)
  if(length(posmissing) > 0){
    countmat[posmissing] <- 1
    dmat <- dmat / countmat
    dmat[posmissing] <- NA_real_
    # At some point we'll need to refactor dineof to work with dist obj
    dmat <- as.dist(dineof(as.matrix(dmat), verbose=Verbose)$X)
  } else {
    dmat <- dmat / countmat
  }

  # Build species tree with NJ
  if(Verbose){
    cat("\n\n  Building species tree...\n")
  }
  newTree <- Treeline(myDistMatrix=dmat, verbose=Verbose, ...)

  if (Verbose){
    dt <- difftime(start, Sys.time())
    cat("  Done.\n  Time difference of", round(abs(dt), 2), attr(dt, "units"), '\n', sep=' ')
  }
  return(newTree)
}
