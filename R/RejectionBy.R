###### -- comprehensive and extensible rejection / retention function ---------

RejectionBy <- function(input,
                        criteria = list("fdr" = 1e-5,
                                        "centroidthreshold" = list("globalpid" = 0.3),
                                        "glmforbiddencolumns" = c("alitype"),
                                        "lmforbiddencolumns" = c("response",
                                                                 "alitype"),
                                        "kargs" = list("max" = 15,
                                                       "scalar" = 4,
                                                       "unitnorm" = TRUE)),
                        rankby = "rawscore",
                        method = "direct",
                        supportedcolumns = c("consensus",
                                             "kmerdist",
                                             "featurediff",
                                             "localpid",
                                             "globalpid",
                                             "matchcoverage",
                                             "localscore",
                                             "deltabackground",
                                             "rawscore",
                                             "response"),
                        dropinappropriate = FALSE) {
  # print(criteria$fdr)
  # current overhead checking
  # function built for internal use, so no printing or progress bars
  # return an integer vector identifying rows in the data to retain
  if (!is(object = input,
          class2 = "data.frame")) {
    stop ("only data.frame inputs are supported")
  }
  # if (!all(colnames(input) %in% supportedcolumns)) {
  #   stop ("a column that is present is not supported")
  # }
  if (!(method %in% c("glm",
                      "kmeans",
                      "lm",
                      "direct",
                      "all"))) {
    stop ("method is not supported")
  }
  # print(criteria$fdr)
  check_this <- names(criteria)[!(names(criteria) == "fdr")]
  # criteria has changed a bit and this isn't the right check to use anymore
  # print(check_this)
  # if (!all(check_this %in% colnames(input))) {
  #   stop ("all supplied criteria must be present in the input data")
  # }
  aliindex <- ifelse(test = input$alitype == "AA",
                     yes = "AA",
                     no = "other")
  input <- input[, supportedcolumns]
  
  if (method == "glm") {
    # perform a knockout based rejection scheme with the spiked in negatives
    # response column must exist, and must contain true/false values
    # run a glm, rank predictions on the fitted value (probability of being true)
    # and return all values above the fdr indicated in the arguments while dropping
    # any false positives that are ranked above the fdr
    w1 <- which(!input$response)
    if (length(w1) < 1) {
      stop ("spiked in negatives must be present for ranked rejection")
    }
    if (!is.null(criteria$glmforbiddencolumns)) {
      # do nothing
    } else {
      # check whether the indicated columns are present
      check_this <- which(colnames(input) %in% criteria$forbiddenglmcolumns)
      # if *any* are, remove them
      if (length(check_this) > 0) {
        input <- input[, -check_this]
      }
    }
    g <- glm(formula = response ~ .,
             family = "quasibinomial",
             data = input)
    FDR <- criteria$fdr
    # total number of positives
    N <- sum(input$response)
    # ordered values
    ranking <- order(g$fitted.values,
                     decreasing = TRUE)
    # return(ranking)
    # take the ordered values up to the point where the FDR climbs past the tolerance
    # then drop known falses from the set
    rate <- cumsum(!input$response[ranking]) / seq(nrow(input))
    # return(rate)
    w2 <- min(which(rate >= FDR))
    retained <- ranking[seq(nrow(input)) < w2]
    retained <- retained[!(retained %in% w1)]
    result <- sort(retained,
                   decreasing = FALSE)
    
  } else if (method == "kmeans") {
    # for now we're just going to do aa vs nt alignments
    
    input$rawscore <- NormVec(input$rawscore)
    # return(list("res" = input,
    #             "index" = aliindex))
    input <- tapply(X = input,
                    INDEX = aliindex,
                    FUN = function(x) {
                      x
                    },
                    simplify = FALSE)
    keylist <- vector(mode = "list",
                      length = length(input))
    attrlist <- vector(mode = "list",
                       length = length(input))
    # at this point we have a LIST
    # of data.frames that whose rownames reference the original row order
    # we have some overhead things to do and then we can just loop through these
    # and return hits that pass
    # response column must be removed if it exists,
    # spiked in negatives are not necessary here and this must be tested for performance
    # with both their inlcusion and exclusion
    # a named criteria value must exist and be present
    # run a kmeans routine analogous to the previous setups, with the addition of splitting
    # nt and aa alignments apart,
    # retain clusters with centroids above non-fdr criterias included
    thresholdcolumn <- names(criteria$centroidthreshold)[1L]
    thresholdval <- criteria$centroidthreshold[[1]]
    # inside the running loop, we check whether the current input data.frame fits 
    # check the size of the objects and use the drop logical to choose whether to
    # retain the entire set, or reject the entire set in cases where there are too few candidates
    # for the unlabeled clustering
    
    for (m1 in seq_along(input)) {
      w1_drop <- which(colnames(input[[m1]]) %in% c("response"))
      if (length(w1_drop) > 0) {
        w1_logical <- !input[[m1]]$response
        w1_integer <- which(w1_logical)
        input[[m1]] <- input[[m1]][, -w1_drop]
      } else {
        # assume there are no negative spiked observations
        # create empty values for later comparisons
        w1_logical <- rep(FALSE, nrow(input[[m1]]))
        w1_integer <- vector(mode = "integer",
                             length = 0L)
      }
      
      rowcheck <- nrow(input[[m1]]) <= criteria$kargs$max
      
      if (dropinappropriate & rowcheck) {
        # just drop everything
        keylist[[m1]] <- vector(mode = "integer",
                                length = 0L)
      } else if (!dropinappropriate & rowcheck) {
        # just retain everything 
        keylist[[m1]] <- seq(nrow(input[[m1]]))
      } else {
        
        nclust <- seq(from = 2,
                      by = 1,
                      to = criteria$kargs$max)
        kmc <- vector(mode = "list",
                      length = length(nclust))
        
        for (m2 in seq_along(kmc)) {
          kmc[[m2]] <- suppressWarnings(kmeans(x = input[[m1]],
                                               centers = nclust[m2],
                                               iter.max = 25L,
                                               nstart = 25L))
          # print(m2)
          # print(nclust[m2])
          
        }
        wss <- vapply(X = kmc,
                      FUN = function(x) {
                        x$tot.withinss
                      },
                      FUN.VALUE = vector(mode = "numeric",
                                         length = 1L))
        # dat2 doesn't seem to behave correctly as a data.frame, so ensure that
        # it's a matrix
        dat1 <- cbind("n" = nclust,
                      "wss" = wss)
        dat2 <- cbind("n" = nclust - 2,
                      "wss" = abs(wss - wss[1]))
        fita <- nls(dat2[, 2L]~OneSite(X = dat2[, 1L],
                                       Bmax,
                                       Kd),
                    start = list(Bmax = max(dat2[, 2L]),
                                 Kd = unname(quantile(dat2[, 1L], .25))))
        fitasum <- summary(fita)
        # fit is offset by -2L to plot and fit correctly, re-offset by +1 to select the correct
        # list position
        # because we're taking the ceiling here, we technically don't need to re-offset
        evalclust <- ceiling((fitasum$coefficients["Kd", "Estimate"] + 1L) * criteria$kargs$scalar)
        if (evalclust >= criteria$kargs$max) {
          warning("Evaluated clusters may be insufficient for this task.")
          evalclust <- criteria$kargs$max
        }
        if (evalclust < 1L) {
          warning("scalar selection requested a number of clusters less than 2, defaulting to 2 clusterings.")
          evalclust <- 1L
        }
        # return(list("a" = fitasum,
        #             "b" = evalclust,
        #             "c" = kmc,
        #             "d" = thresholdcolumn,
        #             "e" = thresholdval))
        kmcselect <- kmc[[evalclust]]
        kmcidentifiers <- kmcselect$cluster
        kmcretainclust <- which(kmcselect$centers[, thresholdcolumn] > thresholdval)
        kmcretainobs <- which((kmcselect$cluster %in% kmcretainclust) & !w1_logical)
        keylist[[m1]] <- rownames(input[[m1]])[kmcretainobs]
      }
      
    }
    result <- as.integer(unlist(keylist))
    
  } else if (method == "lm") {
    # perform a knockout based rejection scheme with the spiked in negatives
    # response column must be removed it if exists
    # will calculate the fdr as above and nuke everything below the cutoff
    # flip the logical for illogical reasons
    
    FDR <- criteria$fdr
    w1_logical <- !input$response
    w1_integer <- which(w1_logical)
    if (length(w1_integer) < 1) {
      stop ("spiked in negatives must be present for ranked rejection")
    }
    if (is.null(criteria$lmforbiddencolumns)) {
      # do nothing
    } else {
      # check whether the indicated columns are present
      # we're automatically removing response here from the input, but it's been assigned
      # as w1 above, so that's ok for now
      check_this <- which(colnames(input) %in% criteria$lmforbiddencolumns)
      # if *any* are, remove them
      if (length(check_this) > 0) {
        input <- input[, -check_this]
      }
    }
    g <- lm(formula = rawscore ~ .,
            data = input)
    ranking <- order(g$fitted.values,
                     decreasing = TRUE)
    rate <- cumsum(w1_logical[ranking]) / seq(nrow(input))
    w2 <- min(which(rate >= FDR))
    retained <- ranking[seq(nrow(input)) < w2]
    retained <- retained[!(retained %in% w1_integer)]
    result <- sort(retained)
  } else if (method == "all") {
    # run all methods and combine them
    stop("not supported yet")
  } else if (method == "direct") {
    # simplest method, just rank the column the user asks for and drop based on the FDR
    FDR <- criteria$fdr
    ranking <- order(input[, rankby],
                     decreasing = TRUE)
    rate <- cumsum(!input$response[ranking]) / seq(nrow(input))
    w2 <- min(which(rate >= FDR))
    retained <- ranking[seq(nrow(input)) < w2]
    retained <- retained[!(retained %in% which(!input$response))]
    result <- sort(retained,
                   decreasing = FALSE)
  }
  
  return(result)
}
