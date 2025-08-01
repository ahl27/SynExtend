#### Helper Functions for EvoWeaver class ####
# author: Aidan Lakshman
# contact: ahl27@pitt.edu
#

#### S3 Generic Definitions ####
PAProfiles <- function(ew, ...) UseMethod('PAProfiles')
CophProfiles <- function(ew, ...) UseMethod('CophProfiles')
RandCophProfiles <- function(ew, ...) UseMethod('RandCophProfiles')
################################

NormArgProcessors <- function(Processors){
  # Normalize argument so it always works
  coresAvailable <- detectCores()
  if(is(Processors, 'numeric')){
    Processors <- as.integer(Processors)
  }
  if(!is(Processors, 'integer')){
    Processors <- 1L
  }
  if(is.null(Processors)){
    Processors <- coresAvailable
  }
  Processors <- max(1L, Processors)
  min(Processors, coresAvailable)
}

BuildSimMatInternal <- function(vecs, uvals, evalmap, l, n, FXN, ARGS, Verbose,
                                 CombinePVal=TRUE, CORRECTION=NULL, InputIsList=FALSE){
  pairscores <- rep(ifelse(CombinePVal, NA_real_, NA_complex_), l*(l-1)/2)
  ctr <- 0
  if (Verbose) pb <- txtProgressBar(max=(l*(l-1) / 2), style=3)
  for ( i in seq_len(l-1) ){
    uval1 <- uvals[i]
    if (InputIsList){
      v1 <- vecs[[i]]
    } else {
      v1 <- vecs[,i]
    }
    for ( j in (i+1):l ){
      uval2 <- uvals[j]
      accessor <- as.character(min(uval1, uval2))
      entry <- max(uval1, uval2)
      if (is.null(evalmap) || entry %in% evalmap[[accessor]]){
        if (InputIsList){
          v2 <- vecs[[j]]
        } else {
          v2 <- vecs[,j]
        }

        pairscores[ctr+1] <- FXN(v1, v2, ARGS, i, j)
      }
      ctr <- ctr + 1
      if (Verbose) setTxtProgressBar(pb, ctr)
    }
  }
  if (Verbose) cat('\n')

  n <- n[uvals]
  if (!is.null(CORRECTION)){
    pairscores <- CORRECTION(pairscores)
  }

  pairscores <- as.simMat(pairscores, NAMES=n, DIAG=FALSE)
  return(pairscores)
}

PAProfiles.EvoWeaver <- function(ew, toEval=NULL, Verbose=TRUE,
                                  speciesList=NULL, ...){
  cols <- names(ew)
  ao <- attr(ew, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(ew, 'useColoc')
  useMT <- attr(ew, 'useMT')
  if (useMT)
    ew <- lapply(ew, labels)
  if (useColoc)
    ew <- lapply(ew, gsub, pattern='([^_]*)_.*', replacement='\\1')

  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  profiles <- matrix(FALSE, nrow=length(allOrgs), ncol=length(ew))
  rownames(profiles) <- allOrgs
  colnames(profiles) <- cols
  if (Verbose) pb <- txtProgressBar(max=length(ew), style=3)
  for ( i in seq_len(length(ew)) ){
    if( !skip || i %in% locs)
      profiles[ew[[i]],i] <- TRUE
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if (Verbose) cat('\n')
  if (!is.null(toEval))
    profiles <- profiles[,locs]
  return(profiles)
}

CophProfiles.EvoWeaver <- function(ew, toEval=NULL, Verbose=TRUE,
                                    speciesList=NULL, ...){
  ## TODO: Some way to handle paralogs
  cols <- names(ew)
  ao <- attr(ew, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }
  useColoc <- attr(ew, 'useColoc')
  useMT <- attr(ew, 'useMT')

  stopifnot('EvoWeaver object must be initialized with dendrograms to run MirrorTree methods'=
              useMT)

  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  l <- length(allOrgs)
  num_entries <- (l * (l-1)) / 2
  outmat <- matrix(0, nrow=num_entries, ncol=length(ew))
  dummycoph <- matrix(NA, nrow=l, ncol=l)
  ut <- upper.tri(dummycoph)
  rownames(dummycoph) <- colnames(dummycoph) <- allOrgs
  if (Verbose) pb <- txtProgressBar(max=length(ew), style=3)
  for ( i in seq_along(ew) ){
    if ( !skip || i %in% locs ){
      dummycoph[] <- NA
      cop <- NA
      # This is occasionally throwing errors that don't affect output for some reason
      cop <- as.matrix(Cophenetic(ew[[i]]))
      copOrgNames <- rownames(cop)
      if (useColoc){
        copOrgNames <- vapply(copOrgNames, gsub, pattern='(.+)_.+_[0-9]+',
                              replacement='\\1', FUN.VALUE=character(1))
        rownames(cop) <- colnames(cop) <- copOrgNames
      }
      dummycoph[copOrgNames, copOrgNames] <- cop
      outmat[,i] <- dummycoph[ut]
    }
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\n')
  colnames(outmat) <- cols
  if (!is.null(toEval)){
    outmat <- outmat[,locs]
    #ltr <- vapply(seq_len(nrow(outmat)), function(x) all(is.na(outmat[x,])),
    #              FUN.VALUE=logical(1))
    #outmat <- outmat[!ltr,]
  }
  return(outmat)
}

RandCophProfiles.EvoWeaver <- function(ew, toEval=NULL, Verbose=TRUE,
                                        speciesList=NULL, outdim=-1,
                                        speciesCorrect=FALSE, mySpeciesTree=NULL,
                                        Processors=1L, ...){
  ## TODO: Some way to handle paralogs
  cols <- names(ew)
  ao <- attr(ew, 'allOrgs')
  if (!is.null(speciesList)){
    stopifnot('Species list is missing species!'=all(ao %in% speciesList))
    allOrgs <- speciesList
  } else {
    allOrgs <- ao
  }

  Processors <- NormArgProcessors(Processors)

  useColoc <- attr(ew, 'useColoc')
  useMT <- attr(ew, 'useMT')

  stopifnot('EvoWeaver object must be initialized with dendrograms to run MirrorTree methods'=
              useMT)

  skip <- FALSE
  if ( !is.null(toEval) ){
    skip <- TRUE
    locs <- unique(c(toEval))
  }
  l <- length(allOrgs)
  outdim <- ifelse(outdim < 1, l, outdim)
  outdim <- as.integer(outdim)
  outmat <- matrix(0, nrow=outdim, ncol=length(ew))
  #dummycoph <- matrix(NA, nrow=l, ncol=l)
  #ut <- upper.tri(dummycoph, diag=FALSE)
  dummycoph <- rep(0, l*(l-1)/2)
  class(dummycoph) <- 'dist'
  attr(dummycoph, "Size") <- l
  attr(dummycoph, "Diag") <- TRUE
  attr(dummycoph, "Upper") <- TRUE
  attr(dummycoph, "Labels") <- allOrgs

  if (speciesCorrect && !is.null(mySpeciesTree)){
    specd <- as.vector(fastCoph(mySpeciesTree))
    #specd[specd==0] <- 1
    spv2 <- specd
    spv2[spv2==0] <- 1
    #nonzeros <- which(specvec != 0)
    #specvec <- .Call("randomProjection", specvec, nonzeros, length(nonzeros), outdim)
  }

  #rownames(dummycoph) <- colnames(dummycoph) <- allOrgs
  if (Verbose) pb <- txtProgressBar(max=length(ew), style=3)
  for ( i in seq_along(ew) ){
    if ( !skip || i %in% locs ){
      dummycoph[] <- 0
      cop <- 0
      # This is occasionally throwing errors that don't affect output for some reason
      #cop <- as.matrix(Cophenetic(ew[[i]]))
      if((!is.null(attr(ew[[i]], 'leaf')) && attr(ew[[i]],'leaf')) || attr(ew[[i]], 'members') <= 1L){
        # edge case where we only have a single leaf in the tree
        copvec <- dummycoph
      } else {
        cop <- fastCoph(ew[[i]])
        #copOrgNames <- rownames(cop)
        copOrgNames <- attr(cop, 'Labels')
        if (useColoc){
          copOrgNames <- vapply(copOrgNames, gsub, pattern='([^_]*)_.*',
                                replacement='\\1', FUN.VALUE=character(1))
          #rownames(cop) <- colnames(cop) <- copOrgNames
          attr(cop, 'Labels') <- copOrgNames
        }
        #dummycoph[copOrgNames, copOrgNames] <- cop
        #copvec <- dummycoph[ut]
        copvec <- combineDist(dummycoph, cop)
        pos <- which(copvec != 0)
        if (speciesCorrect){
          copvec[pos] <- (copvec[pos] - specd[pos]) / spv2[pos]
        }
      }
      copvec <- .Call("randomProjection", copvec,
                      pos, length(pos), outdim,
                      Processors, PACKAGE="SynExtend")
      copvec[copvec==0] <- NA
      outmat[,i] <- copvec
    }
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if(Verbose) cat('\n')
  colnames(outmat) <- cols
  if (!is.null(toEval)){
    outmat <- outmat[,locs]
    #ltr <- vapply(seq_len(nrow(outmat)), function(x) all(is.na(outmat[x,])),
    #              FUN.VALUE=logical(1))
    #outmat <- outmat[!ltr,]
  }
  return(outmat)
}

ProcessSubset <- function(ew, Subset=NULL){
  pl <- length(ew)
  evalmap <- NULL
  uvals <- seq_len(pl)
  if (!is.null(Subset)){
    if (is(Subset, 'data.frame')){
      Subset <- as.matrix(Subset)
    }
    n <- names(ew)
    stopifnot("'Subset' must be either character, numeric, or matrix"=
                (is(Subset, 'character') || is(Subset, 'numeric') || is(Subset, 'matrix')))
    if (is(Subset, 'matrix')){
      if( ncol(Subset) != 2)
        stop('If Subset is a matrix, it must have 2 columns')
      if( is(Subset[1], 'character') ){
        Subset <- matrix(vapply(c(Subset), function(x) {
          val <- which(x==n)
          val <- ifelse(length(val) == 0, -1, val[1])
        }, 0), ncol=2)
        excise <- (Subset[,1] < 0) | (Subset[,2] < 0)
        if (sum(excise) > 0)
          Subset <- Subset[!excise,]
      }
      for ( i in seq_len(nrow(Subset))){
        pos <- Subset[i,]
        i1 <- as.character(min(pos))
        i2 <- max(pos)
        evalmap[[i1]] <- c(evalmap[[i1]], i2)
      }
      uvals <- unique(c(Subset))
    } else {
      if (is(Subset, 'character'))
        Subset <- which(vapply(uvals, function(x) x == n, FUN.VALUE=logical(1)))
      if (length(Subset)==1){
        entry <- Subset[1]
        if (entry > 1){
          evalmap <- lapply(seq_len(entry-1), \(x) entry)
        } else {
          evalmap <- list()
        }

        if (entry < length(n)){
          evalmap[[entry]] <- seq(entry+1, length(n))
        }
        names(evalmap) <- as.character(seq(1, entry))
        uvals <- seq_along(n)
      } else {
        uvals <- unique(Subset)
        evalmap <- lapply(uvals, function(x) uvals)
        names(evalmap) <- as.character(uvals)
      }
    }
  }

  return(list(evalmap=evalmap, uvals=uvals))
}

flatdendrapply <- function(dend, NODEFUN=NULL, LEAFFUN=NODEFUN,
                           INCLUDEROOT=TRUE, ...){
  stopifnot("flatdendrapply only works on objects of class 'dendrogram'"=
              is(dend, 'dendrogram'))
  if (!is(NODEFUN, 'function') && !is(LEAFFUN, 'function'))
    stop("At least one of NODEFUN and LEAFFUN must be a function!")

  val <- lapply(dend,
                \(x){
                  if (is.null(attr(x, 'leaf'))){
                    if (!is(NODEFUN, 'function'))
                      v <- list()
                    else
                      v <- list(NODEFUN(x, ...))
                    for ( child in x ) v <- c(v, Recall(child))
                    return(v)
                  }
                  else if (!is(LEAFFUN, 'function'))
                    return(list())
                  else
                    return(list(LEAFFUN(x, ...)))
                }
  )
  retval <- unlist(val, recursive=FALSE)
  if (!INCLUDEROOT)
    retval[[1]] <- NULL

  lens <- vapply(retval, length, FUN.VALUE=0L)
  atom <- vapply(retval, is.atomic, FUN.VALUE=TRUE)

  if(all(lens == 1L) && all(atom))
    retval <- unlist(retval)

  return(retval)
}

recursive_parentdendrapply <- function(dend, f){
  # f must have two arguments
  stopifnot("function must have two arguments"=length(formals(f)) == 2)
  return(.Call("rpdendrapply", dend, f, new.env(), PACKAGE="SynExtend"))
}

find_dend_distances <- function(dend, useColoc=FALSE){
  if(useColoc){
    dend <- rapply(dend, \(x){
      attr(x, 'label') <- gsub("([^_]*)_.*", '\\1',attr(x, 'label'))
      return(x)
    }, how='replace')
  }
  leafh <- rapply(dend, \(x) c(attr(x, 'label'), attr(x, 'height')))
  leafheights <- as.numeric(leafh[seq(2, length(leafh)+1, 2)])
  names(leafheights) <- leafh[seq(1, length(leafh), 2)]

  roothvec <- rep(attr(dend, 'height'), length(leafheights))
  roothvec <- roothvec - leafheights
  names(roothvec) <- names(leafheights)
  roothvec <- roothvec[order(unlist(dend))]
  attr(dend, 'distances') <- roothvec
  attr(dend, 'branchlen') <- 0

  finddistvec <- function(node, parent){
    curv <- attr(parent, 'distances')
    branchlen <- attr(parent, 'height') - attr(node, 'height')
    pbranchlen <- attr(parent, 'branchlen')
    curv <- curv + 0.5*(branchlen+pbranchlen)
    if (is.null(attr(node, 'leaf'))){
      children <- unlist(node)
    } else {
      children <- attr(node, 'label')
    }
    curv[children] <- curv[children] - branchlen - pbranchlen
    attr(node, 'distances') <- curv
    attr(node, 'branchlen') <- branchlen
    return(node)
  }

  dend <- recursive_parentdendrapply(dend, finddistvec)
  return(dend)
}

AdjMatToDf <- function(preds, Verbose=TRUE, Subset=NULL, CombinePVal=TRUE){
  stopifnot(length(preds) > 0)
  n <- names(preds)
  prednames <- names(preds[[1]])
  lp <- length(prednames)
  v1 <- v2 <- character(lp*(lp+1)/2)
  ctr <- 1
  # can't use expand.grid because of duplicates but this is still fast
  for ( i in seq_len(lp) ){
    for (j in i:lp){
      v1[ctr] <- prednames[i]
      v2[ctr] <- prednames[j]
      ctr <- ctr + 1
    }
  }
  AdjDf <- data.frame(Gene1=v1, Gene2=v2)
  if (Verbose) pb <- txtProgressBar(max=length(n), style=3)
  for ( i in seq_len(length(n))){
    AdjDf[,n[i]] <- unclass(preds[[i]])
    if (Verbose) setTxtProgressBar(pb, i)
  }
  if (Verbose) cat('\n')

  nc <- ncol(AdjDf)
  rtk <- vapply(seq_len(nrow(AdjDf)), function(x) sum(is.na(AdjDf[x,3:nc])) < (nc/2 - 1),
                FUN.VALUE=TRUE)
  AdjDf <- AdjDf[rtk,]
  rownames(AdjDf) <- NULL
  if(!is.null(Subset)){
    ## What if the user wants pair c("B", "A"), but we only have c("A","B")?
    ## This solves that
    AdjDf[,seq_len(2L)] <- t(apply(AdjDf[,seq_len(2L)], 1L, sort))
    Subset[] <- t(apply(Subset, 1L, sort))
    AdjDf <- merge(AdjDf, Subset)
  }

  if(!CombinePVal){
    ## First two columns are the gene names
    pos_comp <- (unlist(lapply(AdjDf, class)) == "complex")[-(c(1,2))]
    AdjDf_names <- AdjDf[,c(1,2)]
    AdjDf_scores <- AdjDf_pvals <- AdjDf[,-(c(1,2))]
    colnames(AdjDf_pvals) <- paste0(colnames(AdjDf_pvals), '.pval')
    colnames(AdjDf_scores) <- paste0(colnames(AdjDf_scores), '.score')

    ## p values only for things that have p-values
    AdjDf_pvals <- AdjDf_pvals[,pos_comp]
    AdjDf_scores[,pos_comp] <- apply(AdjDf_scores[,pos_comp], 2L, Re)
    AdjDf_pvals <- apply(AdjDf_pvals, 2L, Im)

    AdjDf <- cbind(AdjDf_names, AdjDf_scores, AdjDf_pvals)
  }
  return(AdjDf)
}

PAStats <- function(predDf, paps){
  l <- nrow(paps)
  ocv1 <- vapply(predDf[,1], function(x) sum(paps[,x]) / l, 0)
  ocv2 <- vapply(predDf[,2], function(x) sum(paps[,x]) / l, 0)

  d <- abs(ocv1 - ocv2)
  av <- (ocv1 + ocv2) / 2
  return(list(avg=av, diff=d))
}

CheckBifurcating <- function(dend){
  # Checks if a dendrogram is bifurcating
  helperfunc <- function(node){
    if (length(node) == 1) return(TRUE)
    if (length(node) > 2) return(FALSE)
    return(helperfunc(node[[1]]) & helperfunc(node[[2]]))
  }
  return(helperfunc(dend))
}

DCA_minimize_fxn <- function(params, R, spins, i){
  sigi <- spins[,i]
  sigj <- spins[,-i]
  Hi <- params[i]
  Jij <- params[-i]

  sigprods <- sigi * sigj
  firstterm <- colSums(Jij * t(sigprods))
  secondterm <- Hi * sigi
  Si <- mean(exp( -1 * (firstterm + secondterm)))

  regularizer <- R * sum(abs(Jij))
  retval <- log(Si + regularizer)
  if (is.infinite(retval)){
    retval <- -1 * .Machine$double.xmax
  }
  return(retval)
}

DCA_gradient_minimize_fxn <- function(params, R, spins, i){
  grad <- numeric(length(params))

  sigi <- spins[,i]
  sigj <- spins[,-i]
  Hi <- params[i]
  Jij <- params[-i]

  #Calculate gradient
  sigprods <- sigi * sigj
  firstterm <- colSums(Jij * t(sigprods))
  secondterm <- Hi * sigi
  interior <- exp(-1 * firstterm + secondterm)
  entries <- interior * sigi * sigj
  Si <- mean(interior)
  if (is.null(dim(entries))){
    Jpartials <- (mean(entries) / Si) - R
  }
  else
    Jpartials <- (colMeans(entries) / Si) - R
  Hpartial <- mean(interior * sigi) / Si
  grad[-i] <- Jpartials
  grad[i] <- Hpartial
  return(grad)
}

DCA_logrise_run <- function(spins, links, regterm, printProgress=FALSE, Processors=1L){
  Processors <- NormArgProcessors(Processors)
  nnodes <- ncol(spins)
  if (printProgress){
    cat('  Finding topology...\n')
    charsperline <- 73
    charsvec <- diff(floor(seq(0, charsperline, by=charsperline/nnodes))) == 1
  }

  if (printProgress) cat('  |')
  links <- simplify2array(mclapply(seq_len(nnodes),
                                   function(i){
                                     probs <- links[,i]
                                     val <- optim(probs, DCA_minimize_fxn,
                                                  gr=DCA_gradient_minimize_fxn,
                                                  method='BFGS',
                                                  control=list(reltol=1e-6),
                                                  R=regterm, spins=spins, i=i)$par
                                     if (printProgress && charsvec[i])
                                       system2('printf', '=')
                                     return(val)
                                   }, mc.cores=Processors, mc.preschedule = TRUE
  )
  )
  if (printProgress) cat('| ')
  for ( i in seq_len(nrow(links)-1) ){
    for ( j in (i+1):ncol(links)){
      links[i,j] <- links[j,i] <- mean(links[i,j], links[j,i])
    }
  }

  if (printProgress) cat('\n  Done.\n  Eliminating edges close to zero (may take a moment)...\n')
  # Shrink close to zero to zero
  vals <- links[upper.tri(links)]
  h <- hist(vals, breaks=40, plot=FALSE)
  bin0 <- which(h$breaks==0)
  countsNorm <- h$counts == 0
  lb <- length(countsNorm)
  if (length(bin0) != 0 && lb > bin0){
    up <- countsNorm[bin0:lb]
    down <- countsNorm[bin0:1]
    double_up <- which(up[-length(up)] & up[-1])
    double_down <- which(down[-length(down)] & down[-1])
    safeguards <- quantile(vals, c(0.15,0.85))

    if ( length(double_up) == 0 ){
      tmp <- which(up)
      double_up <- ifelse(length(tmp) != 0, which(up)[1], NA)
    } else {
      double_up <- double_up[1]
    }
    uthresh <- ifelse(is.na(double_up), safeguards[2], min(h$breaks[bin0+double_up+1], safeguards[2]))

    if ( length(double_down) == 0 ){
      tmp <- which(down)
      double_down <- ifelse(length(tmp) != 0, which(down)[1], NA)
    } else {
      double_down <- double_down[1]
    }
    lthresh <- ifelse(is.na(double_down), safeguards[1], max(h$breaks[bin0-double_down-1], safeguards[1]))

    d <- diag(links)
    links[(links < uthresh & links > 0) | (links > lthresh & links < 0)] <- 0
    diag(links) <- d
  }

  max_degree <- max(rowSums(links != 0))
  max_val <- max(links) * 2 #conservative estimate on link strength
  bounding <- exp(max_degree * max_val)
  if (printProgress) cat('  Refining edge weights...\n')

  if (printProgress) cat('  |')
  links <- simplify2array(
    mclapply(seq_len(nnodes),
             function(i) {
               probs <- links[,i]
               mask <- seq_len(nnodes) == i
               nonzeros <- (probs != 0) | mask #have to include the self val
               adjustment <- ifelse(i==1, 0, sum(!nonzeros[seq_len(i-1)]))

               if ( sum(nonzeros) > 1 ){
                 pspins <- spins[,nonzeros]
                 pprobs <- probs[nonzeros]
                 val <- optim(pprobs, DCA_minimize_fxn,
                              gr=DCA_gradient_minimize_fxn,
                              method='L-BFGS-B',
                              lower=-bounding, upper=bounding,
                              R=0, spins=pspins, i=(i-adjustment))$par
                 probs[nonzeros] <- val
               }
               if (printProgress && charsvec[i]) system2('printf', '=')
               return(probs)
             },
             mc.cores=Processors, mc.preschedule=TRUE))
  if (printProgress) cat('| ')

  for ( i in seq_len(nrow(links)-1) ){
    for ( j in (i+1):ncol(links)){
      links[i,j] <- links[j,i] <- mean(links[i,j], links[j,i])
    }
  }

  if (max(abs(links)) != 0)
    links <- links / max(abs(links))
  if (printProgress) cat('\n  Done.\n')
  return(links)
}


DCA_logRISE <- function(PAProfiles, niter=1, reg_const=1,
                        Processors=1L, zero_cutoff=0, Verbose=TRUE, ...){

  mult <- 1/niter
  intPA <- PAProfiles + 0
  intPA[intPA==0] <- -1
  nc <- ncol(intPA)
  pp <- FALSE

  if (Verbose & niter==1) pp <- TRUE
  else if (Verbose){
    cat('Running DCA with', niter, 'iterations:\n')
    pb <- txtProgressBar(max=niter, style=3)
  }

  truelinks <- countsmat <- matrix(0, nrow=nc, ncol=nc)
  for ( i in seq_len(niter) ){
    #initlinks <- matrix(0, nrow=nc, ncol=nc)
    initlinks <- matrix(rnorm(nc**2), nrow=nc)
    iterlink <- DCA_logrise_run(intPA, initlinks, reg_const,
                                printProgress=pp, Processors=Processors)
    countsmat <- countsmat + (iterlink != 0)
    truelinks <- truelinks + mult * iterlink
    if(Verbose & !pp) setTxtProgressBar(pb, i)
  }
  if(Verbose & !pp) cat('\n')

  if (zero_cutoff > 0) {
    cutoff <- ifelse(zero_cutoff > 1, zero_cutoff, zero_cutoff * niter)
    countsmat[countsmat < cutoff] <- 0
    countsmat[countsmat > 0] <- 1
    truelinks <- truelinks * countsmat
  }

  return(truelinks)
}

MICalc_C <- function(v1, v2, uv, pseudocount=1L){
  stopifnot("'pseudocount' must be an integer"=is(pseudocount, 'integer'))
  pseudocount <- min(0, pseudocount)
  # Psuedocount=0.01 taken as default by Gerardos et al. (2022) in PLoS
  a <- .Call('calcMIcVec', v1, v2, uv, pseudocount, PACKAGE="SynExtend")
  return(a)
}

CorrComp_C <- function(fm, fsp, ssp, nv, nr){
  on.exit(.C("cleanupFxn"))
  stopifnot(nv == as.integer(nv))
  stopifnot(nr == as.integer(nr))
  stopifnot(all(fsp == as.integer(fsp)))
  stopifnot(all(ssp == as.integer(ssp)))

  a <- .Call('trimCovar', fm, as.integer(fsp),
             as.integer(ssp), as.integer(nv),
             as.integer(nr), PACKAGE="SynExtend")
  return(a)
}

ResidueMIDend <- function(dend1, dend2, cutoff=0.9, comppct=0.25, useColoc, ...){
  if (useColoc){
    l2 <- gsub('([0-9]*)_.*', '\\1', labels(dend1))
    l1 <- gsub('([0-9]*)_.*', '\\1', labels(dend2))
  } else {
    l1 <- labels(dend1)
    l2 <- labels(dend2)
  }
  completeSet <- intersect(l1, l2)
  if (length(completeSet) == 0){
    return(0)
  }

  edges1 <- flatdendrapply(dend1,
                           \(x) list(vals=as.character(unlist(x)),
                                     state=attr(x, 'state')),
                           NULL)
  edges2 <- flatdendrapply(dend2,
                           \(x) list(vals=as.character(unlist(x)),
                                     state=attr(x, 'state')),
                           NULL)


  jsscore <- matrix(Inf, nrow=length(edges1), ncol=length(edges2))
  for ( i in seq_along(edges1) ){
    v1 <- intersect(edges1[[i]]$vals, completeSet)
    for ( j in seq_along(edges2) ){
      v2 <- intersect(edges2[[j]]$vals, completeSet)
      s <- 1 - length(intersect(v1, v2)) / length(union(v1, v2))
      jsscore[i,j] <- ifelse(is.nan(s), 1, s)
    }
  }

  nr <- nrow(jsscore)
  nc <- ncol(jsscore)
  if (nr < nc){
    tm <- edges1
    edges1 <- edges2
    edges2 <- tm
    jsscore <- t(jsscore)
    tm <- nc
    nc <- nr
    nr <- tm
  }
  #now guaranteed to have the larger dimension be nrow

  rownames(jsscore) <- as.character(seq_len(nr))
  colnames(jsscore) <- as.character(seq_len(nc))

  # I'm just using a greedy matching here, couldn't figure out Hungarian
  # and this also scales much better
  pairings <- rep(NA, nc)
  allvals <- rownames(jsscore)
  for ( i in seq_len(nc) ){
    ordered <- allvals[order(jsscore[,i])]
    pos <- which.min(ordered %in% pairings)
    if (jsscore[ordered[pos], i] < cutoff)
      pairings[i] <- ordered[pos]
  }
  # need to catch this here, in case we don't find enough elements just pick
  # random ones
  checksum <- sum(is.na(pairings))
  if (checksum > 0){
    possible <- allvals[!(allvals %in% pairings)]
    pairings[is.na(pairings)] <- sample(possible, checksum)
  }
  names(pairings) <- colnames(jsscore)

  seqset1 <- seqset2 <- NULL
  n <- names(pairings)
  for ( i in seq_along(pairings) ){
    a1 <- as.integer(pairings[i])
    a2 <- as.integer(n[i])
    if (i == 1){
      seqset1 <- BStringSet(edges1[[a1]]$state)
      seqset2 <- BStringSet(edges2[[a2]]$state)
    } else {
      seqset1 <- append(seqset1, edges1[[a1]]$state)
      seqset2 <- append(seqset2, edges2[[a2]]$state)
    }
  }

  names(seqset1) <- names(seqset2) <- seq_len(length(pairings))
  res <- MISeqLevel(seqset1, seqset2, compressionpct=comppct)
  return(res)
}

ResidueMISeqs <- function(seqs1, seqs2, lookup, Processors=1L, CombinePVal, ...){
  baseReturn <- ifelse(CombinePVal, 0, 0+0i)
  if(ncol(seqs1) * ncol(seqs2) == 0){
    return(baseReturn)
  }
  completeSet <- intersect(rownames(seqs1), rownames(seqs2))
  if (length(completeSet) == 0){
    return(baseReturn)
  }
  s1 <- seqs1[completeSet,,drop=FALSE]
  s2 <- seqs2[completeSet,,drop=FALSE]

  nseqs <- length(completeSet)
  if(nseqs == 1){
    return(baseReturn)
  }
  # add gap character
  baseval <- length(lookup)+1L

  AllMIs <- .Call("MIForSequenceSets", s1, s2, nseqs,
                  baseval, baseval, baseval+0.0, Processors)
  AllMIs <- matrix(AllMIs, ncol=ncol(s1))
  # APCm <- matrix(0, nrow=nrow(AllMIs), ncol=ncol(AllMIs))
  # APCm[] <- rowMeans(AllMIs)
  # APCm <- t(t(APCm)*colMeans(AllMIs))
  # AllMIs <- AllMIs - (APCm) / mean(AllMIs)

  #if(ncol(AllMIs) > nrow(AllMIs))
  #  AllMIs <- t(AllMIs)
  meanEnt <- mean(AllMIs)
  sdEnt <- sd(AllMIs)
  #maxVals <- apply(AllMIs, 2L, max)

  ## Greedy pairing
  maxVals <- numeric(min(dim(AllMIs))/2)
  for(i in seq_along(maxVals)){
    p <- arrayInd(which.max(AllMIs), dim(AllMIs))
    maxVals[i] <- AllMIs[p]
    AllMIs <- AllMIs[-p[1],-p[2],drop=FALSE]
  }


  pvals <- pnorm(maxVals, mean=meanEnt, sd=sdEnt, lower.tail=FALSE)
  # Fisher's Method to combine p values
  testStat <- -2 * sum(log(pvals))
  totalP <- pchisq(testStat, df=2*length(pvals), lower.tail=FALSE)

  maxVals <- mean(maxVals)
  totalP <- 1-totalP
  return(ifelse(CombinePVal, maxVals*totalP, complex(real=maxVals, imaginary=totalP)))
}


MISeqLevel <- function(seqSet1, seqSet2, compressionpct=0.25){
  stopifnot('seqSets must be XStringSets'=is(seqSet1, 'XStringSet') && is(seqSet2, 'XStringSet'))
  stopifnot('seqSetq sequences have differing lengths. Ensure you are using an aligned sequence set.'=
              all(width(seqSet1) == width(seqSet1[1])))
  stopifnot('seqSet2 sequences have differing lengths. Ensure you are using an aligned sequence set.'=
              all(width(seqSet2) == width(seqSet2[1])))
  stopifnot('compressionpct must be between 0 and 1'=
              compressionpct < 1 && compressionpct > 0)
  stopifnot('seqSets must be named'=!is.null(names(seqSet1)) && !is.null(names(seqSet2)))
  stopifnot('seqSets must be named'=all(!is.na(c(names(seqSet1), names(seqSet2)))))
  #stopifnot('Both inputs must be DNAStringSets'=
  #            is(seqSet1, 'DNAStringSet') && is(seqSet2, 'DNAStringSet'))
  start2 <- width(seqSet1)[1] + 1
  cali <- ConcatSeqs(seqSet1, seqSet2)
  if (length(cali) == 0){
    #warning('No sequences shared. Check seqSet names!')
    return(0)
  }

  v <- CorrCompressSeqs(cali, start2, mvalpct=compressionpct)
  if (!is.null(v$warn)){
    #warning('Sequences identical.')
    return(1)
  }
  compali <- v$xstrset
  pos <- v$pos
  newstart2 <- which.max(pos >= start2)
  miscore <- CalcMIReduced(compali, newstart2)

  # APC correction
  nr <- nrow(miscore)
  nc <- ncol(miscore)
  if (nr == 0 || nc == 0){
    return(0)
  }
  APC_corr <- matrix(colMeans(miscore), nr, nc, byrow = TRUE) *
    matrix(rowMeans(miscore), nr, nc, byrow = FALSE) / mean(miscore)
  miscore <- miscore - APC_corr

  # scoring
  miscore <- apply(abs(miscore), 2, max)
  retval <- mean(miscore)
  if (is.nan(retval))
    retval <- 0
  return(retval)
}

ConcatSeqs <- function(seqSet1, seqSet2){
  unames <- intersect(names(seqSet1), names(seqSet2))
  concatAli <- xscat(seqSet1[unames], seqSet2[unames])
  names(concatAli) <- unames
  return(concatAli)
}

CorrCompressSeqs <- function(myStringSet, start2, pseudocount=2, mvalpct=0.5,
                             gapLetters=c('-', '.'),
                             uncertainty_cutoff=0.158, MAF_cutoff=0.15){
  freqMat <- consensusMatrix(myStringSet, as.prob=FALSE)
  freqMat <- freqMat[rowSums(freqMat) != 0,]
  freqMat <- freqMat + pseudocount

  freqMat <- t(t(freqMat) / colSums(freqMat))


  to_keep <- rep(FALSE, ncol(freqMat))
  nongaploc <- !(rownames(freqMat) %in% gapLetters)
  for ( i in seq_along(to_keep) ){
    pos <- freqMat[,i]
    pos_no_gap <- pos[nongaploc]

    missing_prob <- sum(pos[gapLetters], na.rm=TRUE)
    vals <- sort(pos_no_gap, decreasing=TRUE)

    MAF <- ifelse(vals[1] == 0, 0, vals[2] / (vals[1] + vals[2]))

    to_keep[i] <- missing_prob < uncertainty_cutoff && MAF > MAF_cutoff
  }
  # Need a guard case here
  if (!any(to_keep)){
    to_keep <- !to_keep
  }
  trimmedFreqMat <- freqMat[,to_keep]
  colnames(trimmedFreqMat) <- which(to_keep)
  fm <- unique(trimmedFreqMat, MARGIN=2)

  nc <- ncol(fm)
  num_vals <- ceiling(mvalpct * nc)

  isInSecondSeq <- as.integer(colnames(fm)) >= start2
  s2 <- which.max(isInSecondSeq)
  firstSeqPos <- seq_len(s2-1)
  if (length(firstSeqPos) == 0){
    return(list(warn=TRUE, 1))
  }
  secondSeqPos <- seq_len(length(isInSecondSeq) - s2 + 1L) + s2 - 1L
  # keeping the entire covariance matrix in memory is really hard
  # this is more code but significantly more memory efficient
  # the cost is slightly more runtime
  corrs <- CorrComp_C(fm, firstSeqPos, secondSeqPos, num_vals, nrow(fm))
  corrs <- unique(corrs)
  if ( length(corrs) < num_vals ){
    fsp <- firstSeqPos[!(firstSeqPos %in% corrs)]
    ssp <- secondSeqPos[!(secondSeqPos %in% corrs)]
    d <- ceiling((num_vals - length(corrs)) / 2)
    l1 <- min(d, length(fsp))
    l2 <- min(d, length(ssp))
    corrs <- c(corrs, sample(fsp, l1), sample(ssp, l2))
  }
  upositions <- sort(as.integer(colnames(fm)[unique(corrs)]))
  subsetXStr <- extractAt(myStringSet, IRanges(upositions, width=1))
  return(list(xstrset=unstrsplit(subsetXStr), pos=upositions))
}

CalcMIReduced <- function(trimmedXStringSet, start2,
                          secondgroupstart=-1){
  matxss <- as.matrix(trimmedXStringSet)
  group1 <- seq_len(start2-1)
  group2 <- seq_len(width(trimmedXStringSet[1]) - start2) + start2

  u <- unique(c(matxss))
  converter <- seq_len(length(u)) - 1L
  names(converter) <- u
  umat <- matrix(converter[matxss], ncol=ncol(matxss))
  scores <- matrix(NA, nrow=length(group1), ncol=length(group2))
  gapnum <- ifelse('-' %in% u, which(u=='-'), -1)
  numunique <- length(u) - (gapnum!=-1)
  ctr <- 0
  for ( i in seq_along(group1) ){
    p1 <- umat[,group1[i]]
    subsetloc <- p1 != gapnum
    for ( j in seq_along(group2) ){
      p2 <- umat[,group2[j]]

      fullsub <- subsetloc & (p2 != gapnum)
      p1p <- p1[fullsub]
      p2p <- p2[fullsub]
      MI <- MICalc_C(p1p, p2p, numunique)
      scores[i,j] <- MI
      ctr <- ctr + 1
    }
  }
  return(scores)
}

predictWithBuiltins <- function(preds, model=c("KEGG", "CORUM")){
  # Key: (val is binary + 1)
  # 000 => PP
  # 001 => PP + PS
  # 010 => PP + GO
  # 011 => PP + PS + GO
  # 100 => PP + SL
  # 101 => PP + PS + SL
  # 110 => PP + GO + SL
  # 111 => PP + PS + GO + SL
  .check_allin <- \(x,y){!(any(match(x,y,nomatch=0L)==0L))}
  model <- match.arg(model)

  pp_algs <- validate_EvoWeaver_methods("PhylogeneticProfiling")$Method
  ps_algs <- validate_EvoWeaver_methods("PhylogeneticStructure")$Method
  go_algs <- validate_EvoWeaver_methods("GeneOrganization")$Method
  sl_algs <- validate_EvoWeaver_methods("SequenceLevel")$Method
  alg_flags <- 0L
  pred_cnames <- colnames(preds)
  if(!all(pp_algs %in% pred_cnames)){
    stop("Phylogenetic profiling algorithms are required for ensemble prediction.")
  }
  if(.check_allin(ps_algs, pred_cnames))
    alg_flags <- alg_flags + 1L
  if(.check_allin(go_algs, pred_cnames))
    alg_flags <- alg_flags + 2L
  if(.check_allin(sl_algs, pred_cnames))
    alg_flags <- alg_flags + 4L

  builtins <- get(data('BuiltInEnsembles', envir=environment()))
  model <- builtins[[model]][[alg_flags+1L]]

  predict(model, preds, type='response')
}

findSpeciesTree <- function(ew, Verbose=TRUE, NameFun=NULL, ...){
  stopifnot("EvoWeaver object must contain dendrograms"=attr(ew, "useMT"))
  if (attr(ew, "useColoc") && is.null(NameFun)){
    NameFun <- function(x) gsub('([^_])_.*', '\\1', x)
  }

  SpecTree <- SuperTree(unclass(ew), NAMEFUN=NameFun,
                        Verbose=Verbose, ...)

  return(SpecTree)
}

## Residue stuff
find_dists_pos <- function(dend, useColoc=FALSE){
  cutoff <- 100L
  dend <- MapCharacters(dend)
  dend <- find_dend_distances(dend, useColoc)
  allPos <- names(sort(table(rapply(dend, \(x){
    as.integer(gsub("[^0-9]([0-9]*)[^0-9]", "\\1", attr(x, 'change')))
  }, how='unlist')), decreasing = TRUE))

  if(length(allPos) > cutoff){
    allPos <- allPos[seq_len(cutoff)]
  }

  # ensure they're always characters
  labs <- as.character(names(attr(dend, 'distances')))
  num_labels <- length(labs)
  posVecs <- lapply(allPos, \(x) {
    y <- rep(Inf, num_labels)
    names(y) <- labs
    return(y)
  })
  names(posVecs) <- allPos

  #names(pm) <- as.character(allPos)
  rapply(dend, \(x) {
    pos <- gsub("[^0-9]([0-9]*)[^0-9]", "\\1", attr(x, 'change'))
    d <- attr(x,'distances')
    n <- names(d)
    for (pv in pos){
      if(!(pv %in% allPos)) next
      curv <- posVecs[[pv]]
      locs <- curv > d
      curv[locs] <- d[locs]
      posVecs[[pv]] <<- curv
    }
    return(NULL)
  }, how='unlist')

  for(i in seq_along(posVecs)){
    p <- posVecs[[i]]
    p[is.infinite(p)] <- NA_real_
    posVecs[[i]] <- p
  }

  return(posVecs)
}

pair_residues <- function(pm1, pm2){
  l1 <- length(pm1)
  l2 <- length(pm2)
  if(l1 >= l2){
    pmL <- pm1
    pmS <- pm2
  } else {
    pmL <- pm2
    pmS <- pm1
    l1 <- l2
    l2 <- length(pm1)
  }
  if(l2 <= 3) return(list(R=0,P=0))
  nameoverlap <- intersect(names(pmL[[1]]), names(pmS[[2]]))
  if(length(nameoverlap) <= 2) return(list(R=0,P=1))
  for(i in seq_len(l1)){
    pmL[[i]] <- pmL[[i]][nameoverlap]
    if(i <= l2)
      pmS[[i]] <- pmS[[i]][nameoverlap]
  }
  ## Greater value should always be columns

  # Max pair because residues can have one to many contacts
  CorrVal <- PVal <- rep(-Inf, l1)
  ns <- names(pmS)
  for(i in seq_len(l1)){
    vCol <- pmL[[i]]
    minCorr <- c(Inf, -1)
    for(j in seq_len(l2)){
      vRow <- pmS[[j]]
      curcorr <- .Call("fastPearsonC", vRow, vCol)

      if (curcorr[1] > CorrVal[i]){
        CorrVal[i] <- curcorr[1]
        # Alternative = 'greater'
        # if it was 'either' it would be
        # p = p < 0.5 ? 2*p : (1-p)*2
        # I don't think 'either' makes sense
        PVal[i] <- 1-pt(curcorr[2], curcorr[3]-2)
      }
      # Negative correlation is meaningless
      # but it shouldn't count against scoring
      if (CorrVal[i] < 0) CorrVal[i] <- 0
    }
  }
  subs <- !is.na(PVal)
  if(!sum(subs)) return(list(R=0,P=1))
  PVal <- PVal[subs]
  CorrVal <- CorrVal[subs]
  # correction for zeros
  PVal[PVal==0] <- 1e-16
  meanP <- 1/sum(1/PVal)
  meanR <- mean(CorrVal)
  return(list(R=meanR, P=meanP))
}

combineDist <- function(dist1, dist2, weightedCombine=FALSE){
  l1 <- attr(dist1, 'Labels')
  l2 <- attr(dist2, 'Labels')
  if(is.null(l1)){
    l1 <- as.character(seq_len(attr(dist1, 'Size')))
    attr(dist1, 'Labels') <- l1
  }
  if(is.null(l2)){
    l2 <- as.character(seq_len(attr(dist2, 'Size')))
    attr(dist2, 'Labels') <- l2
  }
  stopifnot('dist1 must be a superset of dist2'=all(l2 %in% l1))
  pos2 <- vapply(l2, \(x) which(x==l1)[1], integer(1L))
  n1 <- attr(dist1, 'Size')
  n2 <- attr(dist2, 'Size')
  if(weightedCombine){
    mult <- 1/table(l2)
    mult <- mult[l2]
  } else {
    mult <- rep(1.0, length(l2))
  }
  .C('R_combineDistObj', dist1, dist2, pos2, n1, n2, mult)[[1]]
}

fastCoph <- function(dend){
  # DECIPHER:::Cophenetic, but uses new dendrapply for a performance boost
  n <- attr(dend, "members")
  d <- numeric(n*(n - 1)/2)
  u <- unlist(dend)
  o <- order(u)
  u <- u[o]
  labs <- labels(dend)[o]

  # Move unlist() labels into 1:n space
  dend <- rapply(dend,
                 function(y) {
                   y[] <- match(y[1L], u)
                   y
                 },
                 how="replace")

  dendrapply(dend, function(x){
    if(is.leaf(x)) return(x)
    for(k in seq_along(x)){
      h <- attr(x, "height") - attr(x[[k]], "height")
      ## sometimes the rapply loads these as char
      ## not really sure why, but coercing to int is sufficient
      I <- as.integer(unlist(x[[k]]))
      J <- seq_len(n)[-I]
      d <<- .Call("se_cophenetic",
                I,
                J,
                n,
                d,
                h,
                PACKAGE="SynExtend")

    }
    x
  }, how='post.order')

  class(d) <- "dist"
  attr(d, "Size") <- n
  attr(d, "Diag") <- TRUE
  attr(d, "Upper") <- TRUE
  attr(d, "Labels") <- labs

  return(d)
}

rowSumsOfDist <- function(dv){
  n <- attr(dv, 'Size')
  l <- attr(dv, 'Labels')
  outSums <- numeric(n)
  for(i in seq_len(n)){
    lesservals <- seq_len(i-1)
    lesservals <- n*(lesservals-1) - lesservals*(lesservals-1)/2 + i - lesservals
    greatervals <- c()
    if(i+1 <= n){
      greatervals <- seq(i+1, n)
      greatervals <- n*(i-1) - i*(i-1)/2 + greatervals - i
    }
    if(all(is.na(dv[c(lesservals, greatervals)]))){
      outSums[i] <- NA
    } else {
      outSums[i] <- sum(dv[c(lesservals, greatervals)], na.rm=TRUE)
    }
  }
  names(outSums) <- l
  return(outSums)
}

trim_paralogs <- function(xss, verbose=TRUE){
  n <- unique(names(xss))
  lst <- vector('list', length(n))
  names(lst) <- n
  if(verbose) cat('Partitioning sequences...\n')
  if(verbose) pb <- txtProgressBar(style=3, max=length(xss))
  for(i in seq_along(xss)){
    ni <- names(xss)[i]
    if(is.null(lst[[ni]])){
      lst[[ni]] <- xss[i]
    } else {
      lst[[ni]] <- c(lst[[ni]], xss[i])
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\nFinding best sequences...\n')

  outSeqs <- NULL
  if(verbose) pb <- txtProgressBar(style=3, max=length(lst))
  for(i in seq_along(lst)){
    if(length(lst[[i]]) <= 2){
      choice <- 1L
    } else {
      seqs <- lst[[i]]
      st <- seqtype(seqs)
      seqs <- AlignSeqs(seqs, verbose=FALSE, processors=NULL)
      seqs <- MaskAlignment(seqs)
      seqs <- as(seqs, paste0(st, "StringSet"))
      d <- DistanceMatrix(seqs, type='dist', verbose = FALSE)
      d[is.na(d)] <- 2*max(d, na.rm=TRUE)
      sums <- rowSumsOfDist(d)
      choice <- which.min(sums)
    }

    if(is.null(outSeqs)){
      outSeqs <- lst[[i]][choice]
    } else {
      outSeqs <- c(outSeqs, lst[[i]][choice])
    }
    if(verbose) setTxtProgressBar(pb, i)
  }
  if(verbose) cat('\n')

  outSeqs
}
########
