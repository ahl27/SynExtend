###### Distance Matrix Methods for EvoWeaver #####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu

#### Implemented Methods: ####
#  - RPMirrorTree
#  - RPContextTree
#  - TreeDistance
#     -> RF, RF w/ pvalue, Jaccard RF, Cluster Info RF
#     -> Nye Similarity
#     -> Kuhner-Felsenstein
##########################

#### S3 Generic Definitions ####
RPMirrorTree <- function(ew, ...) UseMethod('RPMirrorTree')
RPContextTree <- function(ew, ...) UseMethod('RPContextTree')
TreeDistance <- function(ew, ...) UseMethod('TreeDistance')
################################

RPMirrorTree.EvoWeaver <- function(ew, MTCorrection=c(),
                                  Subset=NULL, Verbose=TRUE,
                                  MySpeciesTree=NULL,
                                  precalcProfs=NULL, precalcSubset=NULL,
                                  CombinePVal=TRUE, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  pl <- length(uvals)

  DIM_LENGTH <- min(80L, length(attr(ew, "allOrgs")))
  alllabs <- lapply(uvals, \(x) labels(ew[[x]]))
  MTCorrection <- tolower(MTCorrection)
  useSpecCorr <- FALSE
  if ('speciestree' %in% MTCorrection){
    if ( is.null(MySpeciesTree) || !is(MySpeciesTree, 'dendrogram')){
      MySpeciesTree <- findSpeciesTree(ew, Verbose)
    }
    stopifnot('Missing MySpeciesTree'=!is.null(MySpeciesTree))
    stopifnot('MySpeciesTree must be a dendrogram'=is(MySpeciesTree, 'dendrogram'))
    useSpecCorr <- TRUE
  }

  if (is.null(precalcProfs)){
    if (Verbose) cat('Pre-processing distance matrices...\n')
    spl <- NULL
    if (!is.null(MySpeciesTree)){
      spl <- labels(MySpeciesTree)
      DIM_LENGTH <- min(80L, length(spl))
    }
    #CPs <- CophProfiles(ew, uvals, Verbose=Verbose, speciesList=spl)
    CPs <- RandCophProfiles(ew, uvals, Verbose=Verbose,
                              speciesList=spl, outdim=as.integer(DIM_LENGTH),
                              speciesCorrect=useSpecCorr,
                              mySpeciesTree=MySpeciesTree, ...)

  } else {
    CPs <- precalcProfs
  }

  l <- ncol(CPs)
  if ( l == 1 ){
    v <- ifelse(CombinePVal, 1, 1+1i)
    mat <- simMat(v, nelem=1)
    rownames(mat) <- colnames(mat) <- uvals
    return(mat)
  }

  if (Verbose) cat('Normalizing profiles...\n')
  means <- colMeans(CPs, na.rm=TRUE)
  vars <- apply(CPs, MARGIN=2, var, na.rm=TRUE)
  if (Verbose) pb <- txtProgressBar(max=ncol(CPs), style=3)
  for ( i in seq_len(ncol(CPs)) ){
    CPs[,i] <- (CPs[,i] - means[i]) / (ifelse(vars[i]!=0, sqrt(vars), 1))
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\n')

  if ('satoaverage' %in% MTCorrection){
    means <- rowMeans(CPs, na.rm = TRUE)
    if (Verbose) cat('Calculating Sato projection vectors...\n')

    # Big profiles lead to space issues that crash R
    if (Verbose) pb <- txtProgressBar(max=ncol(CPs), style=3)
    proj_op <- diag(nrow=nrow(CPs)) - (means %*% t(means))
    for ( i in seq_len(ncol(CPs)) ){
      CPs[,i] <- c(CPs[,i] %*% proj_op)
      if ( Verbose ) setTxtProgressBar(pb, i)
    }
    if (Verbose) cat('\n')
  }

  useOverlap <- ('paoverlap' %in% MTCorrection)
  pairscores <- rep(ifelse(CombinePVal, NA_real_, NA_complex_), pl*(pl-1) / 2)
  ctr <- 0
  endOfRow <- 0
  if (Verbose) pb <- txtProgressBar(max=(pl*(pl-1) / 2), style=3)
  for ( i in seq_len(pl-1) ){
    endOfRow <- endOfRow + (pl-i)
    uval1 <- uvals[i]
    v1 <- CPs[,i]
    l1 <- alllabs[[i]]
    sd1 <- sd(v1, na.rm=TRUE)
    # Should only be NA if there's only one entry
    if (is.na(sd1) || sd1 == 0){
      pairscores[(ctr+1):endOfRow] <- ifelse(CombinePVal, 0, 0+0i)
      ctr <- endOfRow
      if (Verbose) setTxtProgressBar(pb, ctr)
    } else {
      for ( j in (i+1):pl ){
        uval2 <- uvals[j]
        accessor <- as.character(min(uval1, uval2))
        entry <- max(uval1, uval2)
        if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
          v2 <- CPs[,j]
          l2 <- alllabs[[j]]
          sd2 <- sd(v2, na.rm=TRUE)
          if (is.na(sd2) || sd2 == 0)
            val <- ifelse(CombinePVal, 0, 0+0i)
          else{
            val <- suppressWarnings(cor(v1, v2,
                                        use='pairwise.complete.obs',
                                        method='spearman'))
            num_branch <- length(v1)
            pval <- 1 - exp(pt(val, num_branch-2, lower.tail=FALSE, log.p=TRUE))
            overlap <- 1L
            if(useOverlap)
              overlap <- length(intersect(l1,l2)) / length(unique(c(l1,l2)))
            val <- ifelse(CombinePVal, val*pval*overlap, complex(real=val*overlap, imaginary=pval))
          }
          if(is.na(val)) val <- ifelse(CombinePVal, 0, 0+0i)
          pairscores[ctr+1] <- val
        }
        ctr <- ctr + 1
        if (Verbose) setTxtProgressBar(pb, ctr)
      }
    }
  }
  if (Verbose) cat('\n')
  pairscores <- as.simMat(pairscores, NAMES=names(ew)[uvals], DIAG=FALSE)
  Diag(pairscores) <- ifelse(CombinePVal, 1, 1+1i)

  return(pairscores)
}

RPContextTree.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE, precalcProfs=NULL,
                                   MySpeciesTree=NULL, CombinePVal=TRUE, ...){

  if ( is.null(MySpeciesTree) || !is(MySpeciesTree, 'dendrogram')){
    MySpeciesTree <- findSpeciesTree(ew, Verbose)
  }

  #MTCorrection <- c('speciestree', 'normalize', 'partialcorrelation')
  MTCorrection <- c('speciestree', 'paoverlap')

  return(RPMirrorTree(ew, MTCorrection=MTCorrection,
                    Verbose=Verbose,
                    precalcCProfs=precalcProfs,
                    MySpeciesTree=MySpeciesTree,
                    CombinePVal=CombinePVal,
                    Subset=Subset))
}

TreeDistance.EvoWeaver <- function(ew, Subset=NULL, Verbose=TRUE,
                                      precalcSubset=NULL,
                                      TreeMethods="RF", JRFk=4, ...){
  if (!is.null(precalcSubset))
    subs <- precalcSubset
  else
    subs <- ProcessSubset(ew, Subset)
  uvals <- subs$uvals
  evalmap <- subs$evalmap
  if(!is.numeric(JRFk)) stop("JRFk value must be numeric")
  if(is.integer(JRFk)) JRFk <- as.numeric(JRFk)

  useColoc <- attr(ew, "useColoc")
  l <- length(uvals)
  n <- names(ew)
  if ( l == 1 ){
    mat <- matrix(1, nrow=1, ncol=1)
    rownames(mat) <- colnames(mat) <- n
    return(mat)
  }


  bmn <- c("CI", "RF", "JRF", "Nye", "KF", "RFPVal")
  if ('all' %in% TreeMethods){
    bitmask <- rep(TRUE, length(bmn))
  } else {
    bitmask <- rep(FALSE, length(bmn))
    bitmask <- vapply(bmn, \(x) x %in% TreeMethods, logical(1))
  }

  if(bitmask[1]){
    CI_DISTANCE_INTERNAL <- NULL
    data('CIDist_NullDist', package="SynExtend", envir=environment())
  }

  bmn <- bmn[bitmask]
  pairscoresList <- vector('list', length=length(bmn))
  for(i in seq_along(bmn))
    pairscoresList[[i]] <- rep(NA_real_, l*(l-1)/2)
  names(pairscoresList) <- bmn

  ctr <- 0
  pArray <- vector('list', length=l)
  labelsArray <- vector('list', length=l)
  for ( i in seq_len(l) ){
    tree <- ew[[uvals[i]]]
    if (useColoc){
      tree <- rapply(tree, \(x){
        if (!is.null(attr(x, 'leaf'))){
          attr(x, 'label') <- gsub("([^_]*)_.*", '\\1', attr(x, 'label'))
        }
        return(x)
      }, how='replace')
    }
    labs <- labels(tree)
    ptr <- .Call("initCDend", tree, PACKAGE="SynExtend")
    pArray[[i]] <- ptr
    labelsArray[[i]] <- labs
  }
  on.exit(vapply(seq_len(l), \(x){
    ptr <- pArray[[x]]
    rm(ptr)
    return(0L)
  }, integer(1)))

  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    p1 <- pArray[[i]]
    uval1 <- uvals[i]
    for ( j in (i+1):l ){
      p2 <- pArray[[j]]
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        interlabs <- intersect(labelsArray[[i]], labelsArray[[j]])
        # GRF/CI
        if (bitmask[1]){
          # GRF/CI
          s <- .Call("GRFInfo", p1, p2, interlabs, FALSE, 0, PACKAGE="SynExtend")
          normval <- 0.5*(s[2] + s[3])

          if (is.na(normval) || normval == 0){
            pairscoresList$CI[ctr+1] <- 0
            #pairscoresList$GRF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          } else {
            s <- (normval - s[1]) / normval
            s <- max(s, 0)
            pv <- 0
            if(!is.null(CI_DISTANCE_INTERNAL)){
              # p-value correction using random trees
              leninter <- length(interlabs) - 3
              leninter <- min(leninter,197)
              rowtocheck <- c(0, CI_DISTANCE_INTERNAL[2:10,leninter], 1)
              pvals <- c(0,1,5,10,25,50,25,10,5,1,0)/100
              # move score to right side
              pvalind <- ifelse(s > 0.5, 1-s, s)
              ploc <- which(pvalind < rowtocheck)[1]
              if(ploc == 1){
                pv <- 0
              } else {
                pv <- (s - rowtocheck[ploc-1]) / (rowtocheck[ploc] - rowtocheck[ploc-1])
                pv <- pv * abs(pvals[ploc] - pvals[ploc-1]) + min(pvals[ploc-1], pvals[ploc])
                pv <- pv * 2
              }
            }
            pairscoresList$CI[ctr+1] <- (1-s)*(1-pv)
          }
        }
        # RF
        if (bitmask[2]){
          s <- .Call("RFDist", p1, p2, interlabs, PACKAGE="SynExtend")
          normval <- s[2] + s[3]
          if (is.na(normval) || normval == 0)
            pairscoresList$RF[ctr+1] <- 0
            #pairscoresList$RF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else{
            s <- s[1] / normval
            pairscoresList$RF[ctr+1] <- 1 - s
            # p-value calculation
            #s <- ppois(length(interlabs) - 3 - 0.5*s[1],
            #                                  lambda=1/8, log.p=TRUE, lower.tail=FALSE)
            #pairscoresList$RF[ctr+1] <- 1 - exp(s)
          }
        }
        # JRF
        if (bitmask[3]){
          s <- .Call("GRFInfo", p1, p2, interlabs, TRUE, JRFk, PACKAGE="SynExtend")
          normval <- (s[2] + s[3])
          if (is.na(normval) || normval == 0)
            pairscoresList$JRF[ctr+1] <- 0
            #pairscoresList$JRF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$JRF[ctr+1] <- 1 - (s[1] / normval)
        }
        # Nye
        if (bitmask[4]){
          s <- .Call("GRFInfo", p1, p2, interlabs, TRUE, 1, PACKAGE="SynExtend")
          normval <- (s[2] + s[3])
          if (is.na(normval) || normval == 0)
            pairscoresList$Nye[ctr+1] <- 0
            #pairscoresList$Nye[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$Nye[ctr+1] <- 1 - (s[1] / normval)
        }
        # KF
        if (bitmask[5]){
          s <- .Call("KFDist", p1, p2, interlabs, PACKAGE="SynExtend")
          normval <- s[2]
          if (is.na(normval) || normval == 0)
            pairscoresList$KF[ctr+1] <- 0
            #pairscoresList$KF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
          else
            pairscoresList$KF[ctr+1] <- s[1] / normval
        }

        # p-value of RF Dist
        if (bitmask[6]){
            s <- .Call("RFDist", p1, p2, interlabs, PACKAGE="SynExtend")
            normval <- s[2] + s[3]
            if (is.na(normval) || normval == 0)
              pairscoresList$RFPVal[ctr+1] <- 0
            #pairscoresList$RF[ctr+1] <- ifelse(s[1] == 0, 0, 1)
            else{
              # p-value calculation
              s <- ppois(length(interlabs) - 3 - 0.5*s[1],
                          lambda=1/8,
                          log.p=TRUE,
                          lower.tail=FALSE)
              # This probability is p(trees are unrelated), so have to take 1-p
              pairscoresList$RFPVal[ctr+1] <- 1 - exp(s)
            }
        }
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')

  if (bitmask[5]){
    pairscoresList$KF <- pairscoresList$KF - min(pairscoresList$KF, na.rm=TRUE)
    pairscoresList$KF <- pairscoresList$KF / max(pairscoresList$KF, na.rm=TRUE)
    pairscoresList$KF <- 1 - pairscoresList$KF
  }
  n <- n[uvals]
  for (i in seq_along(pairscoresList)){
    pairscoresList[[i]] <- as.simMat(pairscoresList[[i]], NAMES=n, DIAG=FALSE)
  }

  return(pairscoresList)
}
