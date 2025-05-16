###### -- summarize seqs ------------------------------------------------------
# author: nicholas cooley
# maintainer: nicholas cooley
# contact: npc19@pitt.edu
# given a linked pairs object, and a database connection return a pairsummaries
# object
# this function will always align, unlike PairSummaries
# TODO:
# implement:
# ShowPlot
# Processors -- partially implemented, at least for AlignPairs
# ellipses

# this has undergone a pretty significant rewrite, including the implementation
# of some subfunctions to better compartmentalize some complicated tasks
# still need to implement a showplot argument and leverage processors for more
# than just search index and align pairs

SummarizePairs <- function(SynExtendObject,
                           DataBase01,
                           AlignmentFun = "AlignPairs",
                           DefaultTranslationTable = "11",
                           KmerSize = 5,
                           Verbose = FALSE,
                           ShowPlot = FALSE,
                           Processors = 1,
                           Storage = 2,
                           IndexParams = list("K" = 5),
                           SearchParams = list("perPatternLimit" = 0),
                           SearchScheme = "spike",
                           RejectBy = "rank",
                           RetainInternal = FALSE,
                           ...) {
  
  if (Verbose) {
    TimeStart <- Sys.time()
    pBar <- txtProgressBar(style = 1L)
  }
  
  # overhead checking
  # object types
  if (!is(object = SynExtendObject,
          class2 = "LinkedPairs")) {
    stop ("'SynExtendObject' is not an object of class 'LinkedPairs'.")
  }
  Size <- nrow(SynExtendObject)
  # if (!is(object = FeatureSeqs,
  #         class2 = "FeatureSeqs")) {
  #   stop ("'FeatureSeqs' is not an object of class 'FeatureSeqs'.")
  # }
  # we should only need to talk to the DataBase IF FeatureSeqs is not the right length
  
  if (is.character(DataBase01)) {
    if (!requireNamespace(package = "RSQLite",
                          quietly = TRUE)) {
      stop("Package 'RSQLite' must be installed.")
    }
    if (!("package:RSQLite" %in% search())) {
      print("Eventually character vector access to DECIPHER DBs will be deprecated.")
      requireNamespace(package = "RSQLite",
                       quietly = TRUE)
    }
    dbConn <- dbConnect(dbDriver("SQLite"), DataBase01)
    on.exit(dbDisconnect(dbConn))
  } else {
    dbConn <- DataBase01
    if (!dbIsValid(dbConn)) {
      stop("The connection has expired.")
    }
  }
  if (!is.character(AlignmentFun) |
      length(AlignmentFun) > 1) {
    stop("AlignmentFun must be either 'AlignPairs' or 'AlignProfiles'.")
  }
  if (!(AlignmentFun %in% c("AlignPairs",
                            "AlignProfiles"))) {
    stop("AlignmentFun must be either 'AlignPairs' or 'AlignProfiles'.")
  }
  if (!is.character(DefaultTranslationTable) |
      length(DefaultTranslationTable) > 1) {
    stop("DefaultTranslationTable must be a character of length 1.")
  }
  # check storage
  if (Storage < 0) {
    stop("Storage must be greater than zero.")
  } else {
    Storage <- Storage * 1e9 # conversion to gigabytes
  }
  # deal with Processors, this mimics Erik's error checking
  if (!is.null(Processors) && !is.numeric(Processors)) {
    stop("Processors must be a numeric.")
  }
  if (!is.null(Processors) && floor(Processors) != Processors) {
    stop("Processors must be a whole number.")
  }
  if (!is.null(Processors) && Processors < 1) {
    stop("Processors must be at least 1.")
  }
  if (is.null(Processors)) {
    Processors <- DECIPHER:::.detectCores()
  } else {
    Processors <- as.integer(Processors)
  }
  # deal with user arguments
  # ignore 'verbose'
  UserArgs <- list(...)
  # if (length(UserArgs) > 0) {
  #   UserArgNames <- names(UserArgs)
  #   AlignPairsArgs <- formals(AlignPairs)
  #   AlignPairsArgNames <- names(AlignPairsArgs)
  #   AlignProfilesArgs <- formals(AlignProfiles)
  #   AlignProfilesArgNames <- names(AlignProfilesArgs)
  #   DistanceMatrixArgs <- formals(DistanceMatrix)
  #   DistanceMatrixArgNames <- names(DistanceMatrixArgs)
  #   
  #   APaArgNames <- APrArgNames <- DMArgNames <- vector(mode = "character",
  #                                                    length = length(UserArgNames))
  #   for (m1 in seq_along(UserArgNames)) {
  #     APaArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                  choices = AlignPairsArgNames)
  #     APrArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                  choices = AlignProfilesArgNames)
  #     DMArgNames[m1] <- match.arg(arg = UserArgNames[m1],
  #                                 choices = DistanceMatrixArgNames)
  #   }
  #   # set the user args for AlignProfiles
  #   if (any(APaArgNames)) {
  #     
  #   }
  #   # set the user args for AlignPairs
  #   # set the user args for DistanceMatrix
  # }
  
  GeneCalls <- attr(x = SynExtendObject,
                    which = "GeneCalls")
  GeneCallIDs <- names(GeneCalls)
  ObjectIDs <- rownames(SynExtendObject)
  # when we subset a LinkedPairsObject it doesn't smartly handle the genecalls yet...
  if (!all(ObjectIDs %in% GeneCallIDs)) {
    stop("Function expects all IDs in the SynExtendObject to have supplied GeneCalls.")
  }
  
  AA_matrix <- DECIPHER:::.getSubMatrix("PFASUM50")
  NT_matrix <- DECIPHER:::.nucleotideSubstitutionMatrix(2L, -1L, 1L)
  
  feature_match <- match(x = ObjectIDs,
                         table = GeneCallIDs)
  # These need to all be switched over to feature_match
  # we're only going to scroll through the cells that are supplied in the object
  # the datapool only needs to be as long as the gene calls object
  DataPool <- vector(mode = "list",
                     length = length(GeneCallIDs))
  
  MAT1 <- get(data("HEC_MI1",
                   package = "DECIPHER",
                   envir = environment()))
  MAT2 <- get(data("HEC_MI2",
                   package = "DECIPHER",
                   envir = environment()))
  # set initial progress bars and iterators
  # PH is going to be the container that 'res' eventually gets constructed from
  # while Total
  if (AlignmentFun == "AlignProfiles") {
    # progress bar ticks through each alignment
    # this is the slower non-default option
    Total <- (Size - (Size - 1L)) / 2
    PH <- Attr_PH <- vector(mode = "list",
                            length = Total)
    Total <- sum(sapply(SynExtendObject[upper.tri(SynExtendObject)],
                        function(x) nrow(x),
                        USE.NAMES = FALSE,
                        simplify = TRUE))
    
    structureMatrix <- matrix(c(0.187, -0.8, -0.873,
                                -0.8, 0.561, -0.979,
                                -0.873, -0.979, 0.221),
                              3,
                              3,
                              dimnames=list(c("H", "E", "C"),
                                            c("H", "E", "C")))
    substitutionMatrix <- matrix(c(1.5, -2.134, -0.739, -1.298,
                                   -2.134, 1.832, -2.462, 0.2,
                                   -0.739, -2.462, 1.522, -2.062,
                                   -1.298, 0.2, -2.062, 1.275),
                                 nrow = 4,
                                 dimnames = list(DNA_BASES, DNA_BASES))
    
    # read this message after neighbors and k-mer dists ?
    if (Verbose) {
      cat("Aligning pairs.\n")
    }
  } else if (AlignmentFun == "AlignPairs") {
    # progress bar ticks through each cell
    # this is the faster default option
    # technically these alignments have the possibility of not being
    # as good as align profiles, but that's a hit we're willing to take
    Total <- (Size * (Size - 1L)) / 2L
    PH <- Attr_PH <- vector(mode = "list",
                            length = Total)
    Total <- Total * 2L
    if (Verbose) {
      cat("Collecting pairs.\n")
    }
  }
  
  Count <- 1L
  PBCount <- 0L
  # upper key!
  # QueryGene == 1
  # SubjectGene == 2
  # ExactOverlap == 3
  # QueryIndex == 4
  # SubjectIndex == 5
  # QLeft == 6
  # QRight == 7
  # SLeft == 8
  # SRight == 9
  # MaxKmer == 10
  # TotalKmer == 11
  
  # lower key!
  # same minus max and total, but for individual linking kmers
  # not all linking kmers
  block_uid <- 0L
  Prev_m1 <- 0L
  for (m1 in seq_len(Size - 1L)) {
    for (m2 in (m1 + 1L):Size) {
      # print(c(m1, m2))
      # build our data pool first
      # regardless of how we align or if we're including search index or not, we need to prepare and collect
      # the same basic statistics and look aheads
      if (is.null(DataPool[[feature_match[m1]]])) {
        # the pool position is empty, pull from the DB
        # and generate the AAStructures
        DataPool[[feature_match[m1]]]$DNA <- SearchDB(dbFile = dbConn,
                                                      tblName = "NTs",
                                                      identifier = ObjectIDs[m1],
                                                      verbose = FALSE,
                                                      nameBy = "description",
                                                      type = "DNAStringSet")
        DataPool[[feature_match[m1]]]$AA <- SearchDB(dbFile = dbConn,
                                                     tblName = "AAs",
                                                     identifier = ObjectIDs[m1],
                                                     verbose = FALSE,
                                                     nameBy = "description",
                                                     type = "AAStringSet")
        # # return(list("a" = DataPool[[feature_match[m1]]]$DNA,
        # #             "b" = DataPool[[feature_match[m1]]]$AA))
        # print("a")
        # cds values here aren't really CDS values, but an integer value counting
        # the number of feature ranges pulled from the GFFs
        DataPool[[feature_match[m1]]]$len <- width(DataPool[[feature_match[m1]]]$DNA)
        DataPool[[feature_match[m1]]]$mod <- DataPool[[feature_match[m1]]]$len %% 3L == 0
        DataPool[[feature_match[m1]]]$code <- GeneCalls[[feature_match[m1]]]$Coding
        DataPool[[feature_match[m1]]]$cds <- lengths(GeneCalls[[feature_match[m1]]]$Range)
        
        DataPool[[feature_match[m1]]]$struct <- PredictHEC(myAAStringSet = DataPool[[feature_match[m1]]]$AA,
                                                           type = "probabilities",
                                                           HEC_MI1 = MAT1,
                                                           HEC_MI2 = MAT2)
        DataPool[[feature_match[m1]]]$aa_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[m1]]]$AA)
        DataPool[[feature_match[m1]]]$aa_backgrounds <- DataPool[[feature_match[m1]]]$aa_backgrounds[, colnames(AA_matrix)]
        DataPool[[feature_match[m1]]]$aa_backgrounds <- DataPool[[feature_match[m1]]]$aa_backgrounds / rowSums(DataPool[[feature_match[m1]]]$aa_backgrounds)
        DataPool[[feature_match[m1]]]$dna_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[m1]]]$DNA)
        DataPool[[feature_match[m1]]]$dna_backgrounds <- DataPool[[feature_match[m1]]]$dna_backgrounds[, colnames(NT_matrix)]
        DataPool[[feature_match[m1]]]$dna_backgrounds <- DataPool[[feature_match[m1]]]$dna_backgrounds / rowSums(DataPool[[feature_match[m1]]]$dna_backgrounds)
        DataPool[[feature_match[m1]]]$aa_register <- match(table = which(DataPool[[feature_match[m1]]]$mod &
                                                                           DataPool[[feature_match[m1]]]$code),
                                                           x = seq(length(DataPool[[feature_match[m1]]]$DNA)))
        
        DataPool[[feature_match[m1]]]$index <- do.call(what = "IndexSeqs",
                                                       args = c(list("subject" = DataPool[[feature_match[m1]]]$AA,
                                                                     "verbose" = FALSE),
                                                                IndexParams))
        
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
      }
      
      if (is.null(DataPool[[feature_match[m2]]])) {
        # the pool position is empty, pull from DB
        # and generate the AAStructures
        DataPool[[feature_match[m2]]]$DNA <- SearchDB(dbFile = dbConn,
                                                      tblName = "NTs",
                                                      identifier = ObjectIDs[m2],
                                                      verbose = FALSE,
                                                      nameBy = "description",
                                                      type = "DNAStringSet")
        DataPool[[feature_match[m2]]]$AA <- SearchDB(dbFile = dbConn,
                                                     tblName = "AAs",
                                                     identifier = ObjectIDs[m2],
                                                     verbose = FALSE,
                                                     nameBy = "description",
                                                     type = "AAStringSet")
        DataPool[[feature_match[m2]]]$len <- width(DataPool[[feature_match[m2]]]$DNA)
        DataPool[[feature_match[m2]]]$mod <- DataPool[[feature_match[m2]]]$len %% 3L == 0
        DataPool[[feature_match[m2]]]$code <- GeneCalls[[feature_match[m2]]]$Coding
        DataPool[[feature_match[m2]]]$cds <- lengths(GeneCalls[[feature_match[m2]]]$Range)
        
        DataPool[[feature_match[m2]]]$struct <- PredictHEC(myAAStringSet = DataPool[[feature_match[m2]]]$AA,
                                                           type = "probabilities",
                                                           HEC_MI1 = MAT1,
                                                           HEC_MI2 = MAT2)
        DataPool[[feature_match[m2]]]$aa_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[m2]]]$AA)
        DataPool[[feature_match[m2]]]$dna_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[m2]]]$DNA)
        DataPool[[feature_match[m2]]]$aa_backgrounds <- DataPool[[feature_match[m2]]]$aa_backgrounds[, colnames(AA_matrix)]
        DataPool[[feature_match[m2]]]$aa_backgrounds <- DataPool[[feature_match[m2]]]$aa_backgrounds / rowSums(DataPool[[feature_match[m2]]]$aa_backgrounds)
        DataPool[[feature_match[m2]]]$dna_backgrounds <- alphabetFrequency(x = DataPool[[feature_match[m2]]]$DNA)
        DataPool[[feature_match[m2]]]$dna_backgrounds <- DataPool[[feature_match[m2]]]$dna_backgrounds[, colnames(NT_matrix)]
        DataPool[[feature_match[m2]]]$dna_backgrounds <- DataPool[[feature_match[m2]]]$dna_backgrounds / rowSums(DataPool[[feature_match[m2]]]$dna_backgrounds)
        DataPool[[feature_match[m2]]]$aa_register <- match(table = which(DataPool[[feature_match[m2]]]$mod &
                                                                           DataPool[[feature_match[m2]]]$code),
                                                           x = seq(length(DataPool[[feature_match[m2]]]$DNA)))
        
        DataPool[[feature_match[m2]]]$index <- do.call(what = "IndexSeqs",
                                                       args = c(list("subject" = DataPool[[feature_match[m2]]]$AA,
                                                                     "verbose" = FALSE),
                                                                IndexParams))
      } else {
        # the pool position is not empty, assume that it's populated with all the information
        # that it needs
      }
      
      if (Prev_m1 != m1) {
        QueryDNA <- DataPool[[feature_match[m1]]]$DNA
        QueryAA <- DataPool[[feature_match[m1]]]$AA
        QNTCount <- DataPool[[feature_match[m1]]]$len
        QMod <- DataPool[[feature_match[m1]]]$mod
        QCode <- DataPool[[feature_match[m1]]]$code
        QCDSCount <- DataPool[[feature_match[m1]]]$cds
        QueryStruct <- DataPool[[feature_match[m1]]]$struct
        QueryIndex <- DataPool[[feature_match[m1]]]$index
        Query_Background_AA <- DataPool[[feature_match[m1]]]$aa_backgrounds
        Query_Background_NT <- DataPool[[feature_match[m1]]]$dna_backgrounds
        Q_AA_Register <- DataPool[[feature_match[m1]]]$aa_register
        if (KmerSize < 10L) {
          Q_NucF <- oligonucleotideFrequency(x = QueryDNA,
                                             width = KmerSize,
                                             as.prob = TRUE)
        } else {
          stop ("non-overlapping kmers not implemented")
        }
      } else {
        # do something else?
      }
      
      # m2 never stays the same, it will always change in the current search
      # strategy
      SubjectDNA <- DataPool[[feature_match[m2]]]$DNA
      SubjectAA <- DataPool[[feature_match[m2]]]$AA
      SNTCount <- DataPool[[feature_match[m2]]]$len
      SMod <- DataPool[[feature_match[m2]]]$mod
      SCode <- DataPool[[feature_match[m2]]]$code
      SCDSCount <- DataPool[[feature_match[m2]]]$cds
      SubjectStruct <- DataPool[[feature_match[m2]]]$struct
      SubjectIndex <- DataPool[[feature_match[m2]]]$index
      Subject_Background_AA <- DataPool[[feature_match[m2]]]$aa_backgrounds
      Subject_Background_NT <- DataPool[[feature_match[m2]]]$dna_backgrounds
      S_AA_Register <- DataPool[[feature_match[m2]]]$aa_register
      if (KmerSize < 10L) {
        S_NucF <- oligonucleotideFrequency(x = SubjectDNA,
                                           width = KmerSize,
                                           as.prob = TRUE)
      } else {
        stop ("non-overlapping kmers not implemented")
      }
      
      # return(list("QueryData" = list("dna" = QueryDNA,
      #                                "aa" = QueryAA,
      #                                "len" = QNTCount,
      #                                "mod" = QMod,
      #                                "code" = QCode,
      #                                "CDS" = QCDSCount,
      #                                "struct" = QueryStruct,
      #                                "index" = QueryIndex,
      #                                "aa_background" = Query_Background_AA,
      #                                "nt_background" = Query_Background_NT,
      #                                "register" = Q_AA_Register),
      #             "SubjectData" = list("dna" = SubjectDNA,
      #                                  "aa" = SubjectAA,
      #                                  "len" = SNTCount,
      #                                  "mod" = SMod,
      #                                  "code" = SCode,
      #                                  "CDS" = SCDSCount,
      #                                  "struct" = SubjectStruct,
      #                                  "index" = SubjectIndex,
      #                                  "aa_background" = Subject_Background_AA,
      #                                  "nt_background" = Subject_Background_NT,
      #                                  "register" = S_AA_Register)))
      
      # step 1: build indexes into the data pool if they don't exist already
      # # this is already accomplished
      
      # step 2:
      # changing the strategy here for erik:
      # if RejectBy == "none" run as was previously the case
      # if RejectBy == "Rank" run Erik's ranking scheme, i.e. search the reverse and then rank the hits
      # if RejectBy == "kmeans" run what amounts to ClusterByK, include the reverse search as the 
      
      # a forward search will always be performed
      # we can just search in the forward frame,
      # perform a reciprocal search
      # or perform a search for a negative spike
      
      search_df1 <- do.call(what = "SearchIndex",
                            args = c(list("pattern" = QueryAA,
                                          "invertedIndex" = SubjectIndex,
                                          "subject" = SubjectAA,
                                          "verbose" = FALSE,
                                          "processors" = Processors),
                                     SearchParams))
      if (SearchScheme == "standard") {
        # just create the search_pairs object from 
        search_pairs <- data.frame("p1" = names(QueryAA)[search_df1$Pattern],
                                   "p2" = names(SubjectAA)[search_df1$Subject],
                                   "Pattern" = search_df1$Pattern,
                                   "Subject" = search_df1$Subject,
                                   "group" = rep(1L,
                                                 nrow(search_df1)))
        search_pairs$f_hits <- search_df1$Position
        rownames(search_pairs) <- NULL
      } else if (SearchScheme == "spike") {
        # the ranking strategy and the Kmeans strategy will both use the same template of data
        # start with the search for the reverse
        search_df2 <- do.call(what = "SearchIndex",
                              args = c(list("pattern" = reverse(QueryAA),
                                            "invertedIndex" = SubjectIndex,
                                            "subject" = SubjectAA,
                                            "verbose" = FALSE,
                                            "processors" = Processors),
                                       SearchParams))
        search_pairs <- data.frame("p1" = c(names(QueryAA)[search_df1$Pattern],
                                            names(QueryAA)[search_df2$Pattern]),
                                   "p2" = c(names(SubjectAA)[search_df1$Subject],
                                            names(SubjectAA)[search_df2$Subject]),
                                   "Pattern" = c(search_df1$Pattern,
                                                 search_df2$Pattern),
                                   "Subject" = c(search_df1$Subject,
                                                 search_df2$Subject),
                                   "group" = c(rep(1L,
                                                   nrow(search_df1)),
                                               rep(2L,
                                                   nrow(search_df2))))
        search_pairs$f_hits <- c(search_df1$Position,
                                 search_df2$Position)
        rownames(search_pairs) <- NULL
      } else if (SearchScheme == "reciprocal") {
        # the prior strategy is retained by running no rejection
        search_df2 <- do.call(what = "SearchIndex",
                              args = c(list("pattern" = SubjectAA,
                                            "invertedIndex" = QueryIndex,
                                            "subject" = QueryAA,
                                            "verbose" = FALSE,
                                            "processors" = Processors),
                                       SearchParams))
        if (feature_match[m1] == feature_match[m2]) {
          search_df1 <- search_df1[search_df1$Pattern != search_df1$Subject, ]
          search_df2 <- search_df2[search_df2$Pattern != search_df2$Subject, ]
        }
        #direction 1 is pattern -> subject
        #direction 2 is subject -> pattern
        search_df1 <- search_df1[order(search_df1$Pattern,
                                       search_df1$Subject), ]
        search_df2 <- search_df2[order(search_df2$Subject,
                                       search_df2$Pattern), ]
        search_i1 <- paste(search_df1$Pattern,
                           search_df1$Subject,
                           sep = "_")
        search_i2 <- paste(search_df2$Subject,
                           search_df2$Pattern,
                           sep = "_")
        
        # from here if search pairs is zero rows, we're done
        # these names are in the frame of the combined search space,
        # when we subset later 
        search_pairs <- data.frame("p1" = names(QueryAA)[search_df1$Pattern[search_i1 %in% search_i2]],
                                   "p2" = names(SubjectAA)[search_df1$Subject[search_i1 %in% search_i2]],
                                   "Pattern" = search_df1$Pattern[search_i1 %in% search_i2],
                                   "Subject" = search_df1$Subject[search_i1 %in% search_i2],
                                   "group" = rep(1L,
                                                 sum(search_i1 %in% search_i2)))
        search_pairs$f_hits <- search_df1$Position[search_i1 %in% search_i2]
        rownames(search_pairs) <- NULL
      } else {
        stop ("search scheme not recognized.")
      }
      # print(nrow(search_pairs))
      # drop all duplicated pairings every time they appear that isn't the first
      check_this <- paste(search_pairs$Pattern,
                          search_pairs$Subject,
                          sep = "_")
      check_this <- duplicated(check_this)
      if (any(check_this)) {
        search_pairs <- search_pairs[!check_this, ]
      }
      if (nrow(search_pairs) > 0L) {
        
        place_holder1 <- do.call(rbind,
                                 strsplit(x = search_pairs$p1,
                                          split = "_",
                                          fixed = TRUE))
        search_pairs$i1 <- as.integer(place_holder1[, 2L])
        search_pairs$f1 <- as.integer(place_holder1[, 3L])
        place_holder2 <- do.call(rbind,
                                 strsplit(x = search_pairs$p2,
                                          split = "_",
                                          fixed = TRUE))
        search_pairs$i2 <- as.integer(place_holder2[, 2L])
        search_pairs$f2 <- as.integer(place_holder2[, 3L])
        # search_pairs$f_hits <- search_df1$Position[search_i1 %in% search_i2]
        search_pairs$s1 <- GeneCalls[[feature_match[m1]]]$Strand[search_pairs$f1]
        search_pairs$s2 <- GeneCalls[[feature_match[m2]]]$Strand[search_pairs$f2]
        search_pairs$start1 <- GeneCalls[[feature_match[m1]]]$Start[search_pairs$f1]
        search_pairs$start2 <- GeneCalls[[feature_match[m2]]]$Start[search_pairs$f2]
        search_pairs$stop1 <- GeneCalls[[feature_match[m1]]]$Stop[search_pairs$f1]
        search_pairs$stop2 <- GeneCalls[[feature_match[m2]]]$Stop[search_pairs$f2]
        
        # flip strands on search pairs from the reverse search
        # because we only reverse the query when we create the background search
        # we only need to do this for the query strandedness
        search_pairs$s1[search_pairs$group == 2L] <- as.integer(!(as.logical(search_pairs$s1[search_pairs$group == 2L])))
        
        # slam everything together -- i.e. build out rows that need to be
        # added to the linked pairs object
        hit_adjust_start <- do.call(cbind,
                                    search_pairs$f_hits)
        hit_key <- vapply(X = search_pairs$f_hits,
                          FUN = function(x) {
                            ncol(x)
                          },
                          FUN.VALUE = vector(mode = "integer",
                                             length = 1L))
        hit_q_partner <- rep(search_pairs$f1,
                             times = hit_key)
        hit_s_partner <- rep(search_pairs$f2,
                             times = hit_key)
        
        # should be a list in the same shape as the SearchIndex Positions
        # but with the hits mapped to (mostly) the right spots
        # !!! IMPORTANT !!! This is limited to sequences without introns
        # for this to work correctly with introns, I need to be able to divy
        # these offsets up across CDSs, which though possible now will need
        # some significant infrastructure changes to make work cleanly
        
        # return(list("f_hits" = search_pairs$f_hits,
        #             "s1" = search_pairs$s1,
        #             "s2" = search_pairs$s2,
        #             "start1" = search_pairs$start1,
        #             "start2" = search_pairs$start2,
        #             "stop1" = search_pairs$stop1,
        #             "stop2" = search_pairs$stop2))
        
        # this scoping function will make it easier to implement
        # intron/exon re-framing, but it will need access to cds bounds
        # as opposed to just feature bounds
        hit_arrangement <- AAHitScoping(hitlist = search_pairs$f_hits,
                                        fstrand1 = search_pairs$s1,
                                        fstrand2 = search_pairs$s2,
                                        fstart1 = search_pairs$start1,
                                        fstart2 = search_pairs$start2,
                                        fstop1 = search_pairs$stop1,
                                        fstop2 = search_pairs$stop2)
        # print(length(hit_arrangement))
        # return(list(hitlist = search_pairs$f_hits,
        #             fstrand1 = search_pairs$s1,
        #             fstrand2 = search_pairs$s2,
        #             fstart1 = search_pairs$start1,
        #             fstart2 = search_pairs$start2,
        #             fstop1 = search_pairs$stop1,
        #             fstop2 = search_pairs$stop2,
        #             output = hit_arrangement))
        hit_rearrangement <- t(do.call(cbind,
                                       hit_arrangement))
        block_bounds <- lapply(X = hit_arrangement,
                               FUN = function(x) {
                                 c(min(x[1, ]),
                                   max(x[2, ]),
                                   min(x[3, ]),
                                   max(x[4, ]))
                               })
        block_bounds <- do.call(rbind,
                                block_bounds)
        
        hit_widths <- lapply(X = search_pairs$f_hits,
                             FUN = function(x) {
                               x[2, ] - x[1, ] + 1L
                             })
        hit_totals <- vapply(X = hit_widths,
                             FUN = function(x) {
                               sum(x)
                             },
                             FUN.VALUE = vector(mode = "integer",
                                                length = 1))
        hit_max <- vapply(X = hit_widths,
                          FUN = function(x) {
                            max(x)
                          },
                          FUN.VALUE = vector(mode = "integer",
                                             length = 1))
        # print("a")
        # return(list("QueryGene" = hit_q_partner,
        #             "SubjectGene" = hit_s_partner,
        #             "ExactOverlap" = unlist(hit_widths),
        #             "QueryIndex" = rep(search_pairs$i1,
        #                                times = hit_key),
        #             "SubjectIndex" = rep(search_pairs$i2,
        #                                  times = hit_key),
        #             "QLeftPos" = hit_rearrangement[, 1],
        #             "QRightPos" = hit_rearrangement[, 2],
        #             "SLeftPos" = hit_rearrangement[, 3],
        #             "SRightPos" = hit_rearrangement[, 4]))
        add_by_hit <- data.frame("QueryGene" = hit_q_partner,
                                 "SubjectGene" = hit_s_partner,
                                 "ExactOverlap" = unlist(hit_widths),
                                 "QueryIndex" = rep(search_pairs$i1,
                                                    times = hit_key),
                                 "SubjectIndex" = rep(search_pairs$i2,
                                                      times = hit_key),
                                 "QLeftPos" = hit_rearrangement[, 1],
                                 "QRightPos" = hit_rearrangement[, 2],
                                 "SLeftPos" = hit_rearrangement[, 3],
                                 "SRightPos" = hit_rearrangement[, 4])
        add_by_block <- data.frame("QueryGene" = search_pairs$f1,
                                   "SubjectGene" = search_pairs$f2,
                                   "ExactOverlap" = hit_totals,
                                   "QueryIndex" = search_pairs$i1,
                                   "SubjectIndex" = search_pairs$i2,
                                   "QLeftPos" = block_bounds[, 1],
                                   "QRightPos" = block_bounds[, 2],
                                   "SLeftPos" = block_bounds[, 3],
                                   "SRightPos" = block_bounds[, 4],
                                   "MaxKmerSize" = hit_max,
                                   "TotalKmerHits" = hit_key,
                                   "group" = search_pairs$group)
        
        # using forward hits from here:
        # append pairs that don't appear in the current linked pairs object
        # onto the current linked pairs object
        
        
        # hits are now transposed into the context of the whole sequence
        # as opposed to the feature
        # build out rows to add, and then add them
        if (nrow(SynExtendObject[[m1, m2]]) > 0) {
          select_row1 <- paste(SynExtendObject[[m1, m2]][, 1L],
                               SynExtendObject[[m1, m2]][, 4L],
                               SynExtendObject[[m1, m2]][, 2L],
                               SynExtendObject[[m1, m2]][, 5L],
                               sep = "_")
          select_row2 <- paste(SynExtendObject[[m2, m1]][, 1L],
                               SynExtendObject[[m2, m1]][, 4L],
                               SynExtendObject[[m2, m1]][, 2L],
                               SynExtendObject[[m2, m1]][, 5L],
                               sep = "_")
        } else {
          select_row1 <- select_row2 <- vector(mode = "character",
                                               length = 0)
        }
        # upper diagonal is blocks,
        # lower diagonal is hits
        select_row3 <- paste(add_by_block$QueryGene,
                             add_by_block$QueryIndex,
                             add_by_block$SubjectGene,
                             add_by_block$SubjectIndex,
                             sep = "_")
        select_row4 <- paste(add_by_hit$QueryGene,
                             add_by_hit$QueryIndex,
                             add_by_hit$SubjectGene,
                             add_by_hit$SubjectIndex,
                             sep = "_")
        sr_upper <- !(select_row3 %in% select_row1)
        sr_lower <- !(select_row4 %in% select_row2)
        
        # return(list("upper1" = SynExtendObject[[m1, m2]],
        #             "lower1" = SynExtendObject[[m2, m1]],
        #             "upper2" = as.matrix(add_by_block[sr_upper, ]),
        #             "lower2" = as.matrix(add_by_hit[sr_lower, ])))
        
        # return(list(select_row1,
        #             select_row2,
        #             select_row3,
        #             select_row4,
        #             SynExtendObject[[m1, m2]],
        #             SynExtendObject[[m2, m1]],
        #             add_by_hit,
        #             add_by_block))
        if (any(sr_upper)) {
          SynExtendObject[[m1, m2]] <- cbind(SynExtendObject[[m1, m2]],
                                             "group" = rep(1L,
                                                           nrow(SynExtendObject[[m1, m2]])))
          SynExtendObject[[m1, m2]] <- rbind(SynExtendObject[[m1, m2]],
                                             as.matrix(add_by_block[sr_upper, ]))
          SynExtendObject[[m1, m2]] <- SynExtendObject[[m1, m2]][order(SynExtendObject[[m1, m2]][, 1L],
                                                                       SynExtendObject[[m1, m2]][, 2L]), ]
          rownames(SynExtendObject[[m1, m2]]) <- NULL
          SynExtendObject[[m2, m1]] <- rbind(SynExtendObject[[m2, m1]],
                                             as.matrix(add_by_hit[sr_lower, ]))
          SynExtendObject[[m2, m1]] <- SynExtendObject[[m2, m1]][order(SynExtendObject[[m2, m1]][, 1L],
                                                                       SynExtendObject[[m2, m1]][, 2L]), ]
          rownames(SynExtendObject[[m2, m1]]) <- NULL
        } else {
          SynExtendObject[[m1, m2]] <- cbind(SynExtendObject[[m1, m2]],
                                             "group" = rep(1L,
                                                           nrow(SynExtendObject[[m1, m2]])))
        }
      } else {
        stop ("an unexpected condition occured, please contact the maintainer -- search_pairs is empty when it is not expected to be")
        # it is technically possible to not find anything with search but still have synteny hits to parse,
        # i need to figure out how to handle this and under what circumstances it occurs
      }
      # return(list("upper" = SynExtendObject[[m1, m2]],
      #             "lower" = SynExtendObject[[m2, m1]]))
      # print("b")
      # next step if RejectBy is not none, run the rejection scheme
      
      
      # step 3: morph the searches into the SynExtendObject so other things
      # go smoothly
      
      # alignments happen later and profile vs pairs is chosen by the user
      
      # if (m1 == 1 & m2 == 3) {
      #   return(list(SynExtendObject[[m1, m2]],
      #               SynExtendObject[[m2, m1]]))
      # }
      
      
      ###### -- only evaluate valid positions ---------------------------------
      if (nrow(SynExtendObject[[m1, m2]]) > 0L) {
        # links table is populated, do whatever
        PMatrix <- cbind(SynExtendObject[[m1, m2]][, 1L],
                         SynExtendObject[[m1, m2]][, 2L])
        
        IMatrix <- cbind(SynExtendObject[[m1, m2]][, 4L],
                         SynExtendObject[[m1, m2]][, 5L])
        # find the Consensus of each linking kmer hit
        
        strand1_adj <- rep(SynExtendObject[[m1, m2]][, 12],
                           times = SynExtendObject[[m1, m2]][, 11])
        strand1_ph <- GeneCalls[[feature_match[m1]]]$Strand[SynExtendObject[[m2, m1]][, 1L]]
        strand1_ph[strand1_adj == 2L] <- as.integer(!(as.logical(strand1_ph[strand1_adj == 2L])))
        
        # return(list(gene1left = GeneCalls[[feature_match[m1]]]$Start[SynExtendObject[[m2, m1]][, 1L]],
        #             gene2left = GeneCalls[[feature_match[m2]]]$Start[SynExtendObject[[m2, m1]][, 2L]],
        #             gene1right = GeneCalls[[feature_match[m1]]]$Stop[SynExtendObject[[m2, m1]][, 1L]],
        #             gene2right = GeneCalls[[feature_match[m2]]]$Stop[SynExtendObject[[m2, m1]][, 2L]],
        #             hit1left = SynExtendObject[[m2, m1]][, 6L],
        #             hit1right = SynExtendObject[[m2, m1]][, 7L],
        #             hit2left = SynExtendObject[[m2, m1]][, 8L],
        #             hit2right = SynExtendObject[[m2, m1]][, 9L],
        #             strand1 = strand1_ph,
        #             strand2 = GeneCalls[[feature_match[m2]]]$Strand[SynExtendObject[[m2, m1]][, 2L]],
        #             index1 = SynExtendObject[[m2, m1]][, 1L],
        #             index2 = SynExtendObject[[m2, m1]][, 2L],
        #             group = SynExtendObject[[m1, m2]][, "group"],
        #             upper = SynExtendObject[[m1, m2]],
        #             lower = SynExtendObject[[m2, m1]]))
        
        diff1 <- HitConsensus(gene1left = GeneCalls[[feature_match[m1]]]$Start[SynExtendObject[[m2, m1]][, 1L]],
                              gene2left = GeneCalls[[feature_match[m2]]]$Start[SynExtendObject[[m2, m1]][, 2L]],
                              gene1right = GeneCalls[[feature_match[m1]]]$Stop[SynExtendObject[[m2, m1]][, 1L]],
                              gene2right = GeneCalls[[feature_match[m2]]]$Stop[SynExtendObject[[m2, m1]][, 2L]],
                              hit1left = SynExtendObject[[m2, m1]][, 6L],
                              hit1right = SynExtendObject[[m2, m1]][, 7L],
                              hit2left = SynExtendObject[[m2, m1]][, 8L],
                              hit2right = SynExtendObject[[m2, m1]][, 9L],
                              strand1 = strand1_ph,
                              strand2 = GeneCalls[[feature_match[m2]]]$Strand[SynExtendObject[[m2, m1]][, 2L]])
        # print("c")
        
        
        # get the mean consensus
        diff2 <- vector(mode = "numeric",
                        length = nrow(SynExtendObject[[m1, m2]]))
        hit_relations <- vector(mode = "integer",
                                length = nrow(SynExtendObject[[m2, m1]]))
        hit_iterator <- 0L
        loop_iterator <- 1L
        h1 <- 0L
        h2 <- 0L
        continue <- TRUE
        while (continue) {
          if (SynExtendObject[[m2, m1]][loop_iterator, 1L] == h1 &
              SynExtendObject[[m2, m1]][loop_iterator, 2L] == h2) {
            hit_relations[loop_iterator] <- hit_iterator
            loop_iterator <- loop_iterator + 1L
          } else {
            hit_iterator <- hit_iterator + 1L
            hit_relations[loop_iterator] <- hit_iterator
            h1 <- SynExtendObject[[m2, m1]][loop_iterator, 1L]
            h2 <- SynExtendObject[[m2, m1]][loop_iterator, 2L]
            loop_iterator <- loop_iterator + 1L
          }
          # setTxtProgressBar(pb = pBar,
          #                   value = loop_iterator / nrow(SynExtendObject[[m2, m1]]))
          if (loop_iterator > nrow(SynExtendObject[[m2, m1]])) {
            continue <- FALSE
          }
        }
        
        diff2 <- 1 - unname(tapply(X = diff1,
                                   INDEX = hit_relations,
                                   FUN = function(x) {
                                     mean(x)
                                   }))
        # max match size
        MatchMax <- SynExtendObject[[m1, m2]][, "MaxKmerSize"]
        # total unique matches
        UniqueMatches <- SynExtendObject[[m1, m2]][, "TotalKmerHits"]
        # total matches at all
        TotalMatch <- SynExtendObject[[m1, m2]][, "ExactOverlap"]
        
        # print("d")
        # return(list(p1 = SynExtendObject[[m1, m2]][, 1L],
        #             p2 = SynExtendObject[[m1, m2]][, 2L],
        #             code1 = QCode[SynExtendObject[[m1, m2]][, 1L]],
        #             code2 = SCode[SynExtendObject[[m1, m2]][, 2L]],
        #             mod1 = QMod[SynExtendObject[[m1, m2]][, 1L]],
        #             mod2 = SMod[SynExtendObject[[m1, m2]][, 2L]],
        #             aa1 = Query_Background_AA,
        #             aa2 = Subject_Background_AA,
        #             nt1 = Query_Background_NT,
        #             nt2 = Subject_Background_NT,
        #             register1 = Q_AA_Register,
        #             register2 = S_AA_Register,
        #             aamat = AA_matrix,
        #             ntmat = NT_matrix))
        
        diff3 <- ApproximateBackground(p1 = SynExtendObject[[m1, m2]][, 1L],
                                       p2 = SynExtendObject[[m1, m2]][, 2L],
                                       code1 = QCode[SynExtendObject[[m1, m2]][, 1L]],
                                       code2 = SCode[SynExtendObject[[m1, m2]][, 2L]],
                                       mod1 = QMod[SynExtendObject[[m1, m2]][, 1L]],
                                       mod2 = SMod[SynExtendObject[[m1, m2]][, 2L]],
                                       aa1 = Query_Background_AA,
                                       aa2 = Subject_Background_AA,
                                       nt1 = Query_Background_NT,
                                       nt2 = Subject_Background_NT,
                                       register1 = Q_AA_Register,
                                       register2 = S_AA_Register,
                                       aamat = AA_matrix,
                                       ntmat = NT_matrix)
        # print("e")
        # from here we need to get the kmer differences
        # the PIDs
        # the SCOREs
        
        # if both positions are present in the FeatureSeqs object, do nothing
        # if either or both is missing, they need to be pulled from the DB
        # this functionality needs to be expanded eventually to take in cases where the
        # we're overlaying something with gene calls on something that doesn't have gene calls
        # if (!all(ObjectIDs[c(m1, m2)] %in% FeatureSeqs$IDs)) {
        #   # an object ID does not have a seqs present, pull them
        #   # this is not a priority so we're leaving this blank for a second
        # } else {
        #   TMPSeqs01 <- FALSE
        #   TMPSeqs02 <- FALSE
        # }
        
        # align everyone as AAs who can be, i.e. modulo of 3, is coding, etc
        # then align everyone else as nucs
        # translate the hit locations to the anchor positions
        # every hit is an anchor
        # all hits are stored in the LinkedPairs object in nucleotide space
        # as left/right bounds, orientations will be dependant upon strandedness of the genes
        
        # prepare the kmer distance stuff:
        NucDist <- vector(mode = "numeric",
                          length = nrow(PMatrix))
        
        QueryFeatureLength <- QNTCount[PMatrix[, 1L]]
        SubjectFeatureLength <- SNTCount[PMatrix[, 2L]]
        
        # Perform the alignments
        if (AlignmentFun == "AlignPairs") {
          
          # grab kmer distances ahead of time because AlignPairs doesn't need to loop
          # through anything
          for (m3 in seq_along(NucDist)) {
            it1 <- PMatrix[m3, 1L]
            it2 <- PMatrix[m3, 2L]
            NucDist[m3] <- sqrt(sum((Q_NucF[it1, ] - S_NucF[it2, ])^2)) / ((sum(Q_NucF[it1, ]) + sum(S_NucF[it2, ])) / 2)
            # NucDist[m3] <- sqrt(sum((nuc1[m3, ] - nuc2[m3, ])^2)) / ((sum(nuc1[m3, ]) + sum(nuc2[m3, ])) / 2)
          }
          
          # print("f")
          # spit out the subset vectors and logicals to correctly call both AlignPairs calls
          # and both dfs
          AASelect <- PMatrix[, 1L] %in% which(QCode & QMod) & PMatrix[, 2L] %in% which(SCode & SMod)
          NTSelect <- !AASelect
          
          df_aa <- data.frame("Pattern" = Q_AA_Register[PMatrix[AASelect, 1L]],
                              "Subject" = S_AA_Register[PMatrix[AASelect, 2L]])
          df_nt <- data.frame("Pattern" = PMatrix[NTSelect, 1L],
                              "Subject" = PMatrix[NTSelect, 2L])
          
          check1 <- paste(search_pairs$Pattern,
                          search_pairs$Subject,
                          sep = "_")
          check2 <- paste(df_aa$Pattern,
                          df_aa$Subject,
                          sep = "_")
          aa_pos <- vector(mode = "list",
                           length = nrow(df_aa))
          aa_match1 <- match(x = check2,
                             table = check1)
          aa_match2 <- which(!is.na(aa_match1))
          aa_match3 <- which(is.na(aa_match1))
          aa_pos[aa_match2] <- search_pairs$f_hits[aa_match1[!is.na(aa_match1)]]
          aa_pos[aa_match3] <- mapply(SIMPLIFY = FALSE,
                                      USE.NAMES = FALSE,
                                      FUN = function(y, z) {
                                        cbind(matrix(data = 0L,
                                                     nrow = 4),
                                              matrix(data = c(y,y,z,z),
                                                     nrow = 4))
                                      },
                                      y = width(QueryAA)[df_aa$Pattern[aa_match3]] + 1L,
                                      z = width(SubjectAA)[df_aa$Subject[aa_match3]] + 1L)
          nt_pos <- mapply(SIMPLIFY = FALSE,
                           USE.NAMES = FALSE,
                           FUN = function(y, z) {
                             cbind(matrix(data = 0L,
                                          nrow = 4),
                                   matrix(data = c(y,y,z,z),
                                          nrow = 4))
                           },
                           y = QNTCount[df_nt$Pattern] + 1L,
                           z = SNTCount[df_nt$Subject] + 1L)
          
          # return(list("a" = search_pairs,
          #             "b" = df_aa,
          #             "c" = df_nt,
          #             "e" = Q_AA_Register,
          #             "f" = S_AA_Register,
          #             "g" = check1,
          #             "h" = check2))
          # when patterns are set to zero, almost everything will have a local match
          # just set those from their object above,
          # and if a candidate set doesn't search index hits, just set the anchors to global
          # df_aa$Position <- WithinQueryAAs
          # df_aa$Pattern <- Q_AA_Register[df_aa$Pattern]
          # df_aa$Subject <- S_AA_Register[df_aa$Subject]
          # erik's request is to drop direct global alignment
          # this means that global alignments are now approximations and not true global alignments
          # except in the case where an inferred pair came from find synteny alone and it's bounds couldn't
          # be reconciled for align pairs
          df_aa$Position <- aa_pos
          df_nt$Position <- nt_pos
          # return(list(WithinQueryAAs,
          #             WithinQueryNucs,
          #             QueryAA,
          #             SubjectAA))
          
          # done with anchors at this point
          
          if (sum(AASelect) > 0) {
            # return(list("q" = QueryAA,
            #             "s" = SubjectAA,
            #             "df" = df_aa))
            aapairs <- AlignPairs(pattern = QueryAA,
                                  subject = SubjectAA,
                                  pairs = df_aa,
                                  verbose = FALSE,
                                  processors = Processors)
            current_local_aa_pids <- aapairs$Matches / aapairs$AlignmentLength
            current_global_aa_pids <- aapairs$Matches / pmax(width(QueryAA)[aapairs$Pattern],
                                                             width(SubjectAA)[aapairs$Subject])
            current_local_aa_scores <- aapairs$Score / aapairs$AlignmentLength
            current_global_aa_scores <- aapairs$Score / pmax(width(QueryAA)[aapairs$Pattern],
                                                             width(SubjectAA)[aapairs$Subject])
            current_global_aa_rawscore <- aapairs$Score
          } else {
            current_local_aa_pids <- numeric()
            current_global_aa_pids <- numeric()
            current_local_aa_scores <- numeric()
            current_global_aa_scores <- numeric()
            current_global_aa_rawscore <- numeric()
          }
          if (Verbose) {
            PBCount <- PBCount + 1L
            setTxtProgressBar(pb = pBar,
                              value = PBCount / Total)
          }
          
          if (sum(NTSelect) > 0) {
            ntpairs <- AlignPairs(pattern = QueryDNA,
                                  subject = SubjectDNA,
                                  pairs = df_nt,
                                  verbose = FALSE,
                                  processors = Processors)
            current_local_nt_pids <- ntpairs$Matches / ntpairs$AlignmentLength
            current_global_nt_pids <- ntpairs$Matches / pmax(width(QueryDNA)[ntpairs$Pattern],
                                                             width(SubjectDNA)[ntpairs$Subject])
            current_local_nt_scores <- ntpairs$Score / ntpairs$AlignmentLength
            current_global_nt_scores <- ntpairs$Score / pmax(width(QueryDNA)[ntpairs$Pattern],
                                                             width(SubjectDNA)[ntpairs$Subject])
            current_global_nt_rawscore <- ntpairs$Score
          } else {
            current_local_nt_pids <- numeric()
            current_global_nt_pids <- numeric()
            current_local_nt_scores <- numeric()
            current_global_nt_scores <- numeric()
            current_global_nt_rawscore <- numeric()
          }
          if (Verbose) {
            PBCount <- PBCount + 1L
            setTxtProgressBar(pb = pBar,
                              value = PBCount / Total)
          }
          
          vec1 <- vec2 <- vec3 <- vec4 <- vec5 <- vector(mode = "numeric",
                                                         length = nrow(PMatrix))
          vec1[AASelect] <- current_local_aa_pids
          vec1[NTSelect] <- current_local_nt_pids
          vec2[AASelect] <- current_local_aa_scores
          vec2[NTSelect] <- current_local_nt_scores
          vec3[AASelect] <- current_global_aa_pids
          vec3[NTSelect] <- current_global_nt_pids
          vec4[AASelect] <- current_global_aa_scores
          vec4[NTSelect] <- current_global_nt_scores
          vec5[AASelect] <- current_global_aa_rawscore
          vec5[NTSelect] <- current_global_nt_rawscore
          # For Testing
          # AA_Anchors[[Count]] <- WithinQueryAAs
          # NT_Anchors[[Count]] <- WithinQueryNucs[NTSelect]
          # return(list(aapairs,
          #             ntpairs))
          
        } else if (AlignmentFun == "AlignProfiles") {
          
          stop ("currently not supported due to other priorities...")
          # build out a map of who is being called where
          # if we're aligning nucleotides, just call the position in the DNA
          # stringsets,
          # if you not, use the matched positions
          # ws1 <- match(table = names(QueryAA),
          #              x = names(QueryDNA)[PMatrix[, 1L]])
          # ws2 <- match(table = names(SubjectAA),
          #              x = names(SubjectDNA)[PMatrix[, 2L]])
          # vec1 <- vec2 <- vector(mode = "numeric",
          #                        length = length(NucDist))
          # 
          # for (m3 in seq_along(NucDist)) {
          #   
          #   NucDist[m3] <- sqrt(sum((nuc1[m3, ] - nuc2[m3, ])^2)) / ((sum(nuc1[m3, ]) + sum(nuc2[m3, ])) / 2)
          #   
          #   if (AASelect[m3]) {
          #     # align as amino acids
          #     ph1 <- AlignProfiles(pattern = QueryAA[ws1[m3]],
          #                          subject = SubjectAA[ws2[m3]],
          #                          p.struct = QueryStruct[ws1[m3]],
          #                          s.struct = SubjectStruct[ws2[m3]])
          #     ph2 <- DistanceMatrix(myXStringSet = ph1,
          #                           includeTerminalGaps = TRUE,
          #                           type = "matrix",
          #                           verbose = FALSE)
          #     ph3 <- ScoreAlignment(myXStringSet = ph1,
          #                           structures = PredictHEC(myAAStringSet = ph1,
          #                                                   type = "probabilities",
          #                                                   HEC_MI1 = MAT1,
          #                                                   HEC_MI2 = MAT2),
          #                           structureMatrix = structureMatrix)
          #   } else {
          #     # align as nucleotides
          #     ph1 <- AlignProfiles(pattern = QueryDNA[PMatrix[m3, 1L]],
          #                          subject = SubjectDNA[PMatrix[m3, 2L]])
          #     ph2 <- DistanceMatrix(myXStringSet = ph1,
          #                           includeTerminalGaps = TRUE,
          #                           type = "matrix",
          #                           verbose = FALSE)
          #     ph3 <- ScoreAlignment(myXStringSet = ph1,
          #                           substitutionMatrix = substitutionMatrix)
          #   }
          #   vec1[m3] <- 1 - ph2[1, 2]
          #   vec2[m3] <- ph3
          #   
          #   if (Verbose) {
          #     PBCount <- PBCount + 1L
          #     setTxtProgressBar(pb = pBar,
          #                       value = PBCount / Total)
          #   }
          # } # end m3 loop
          
        } # end if else on alignment function
        
        internal_vals <- data.frame("consensus" = diff2,
                                    "featurediff" = abs(QueryFeatureLength - SubjectFeatureLength) / pmax(QueryFeatureLength,
                                                                                                          SubjectFeatureLength),
                                    "kmerdist" = NucDist,
                                    "localpid" = vec1,
                                    "globalpid" = vec3,
                                    "matchcoverage" = (TotalMatch * 2L) / (QueryFeatureLength + SubjectFeatureLength),
                                    "localscore" = vec2,
                                    "deltabackground" = vec4,
                                    "rawscore" = vec5,
                                    "response" = ifelse(test = SynExtendObject[[m1, m2]][, "group"] == 1L,
                                                        yes = TRUE,
                                                        no = FALSE),
                                    "alitype" = ifelse(test = AASelect,
                                                       yes = "AA",
                                                       no = "NT"))
        if (RejectBy == "glm") {
          w_retain <- RejectionBy(input = internal_vals,
                                  method = "glm")
        } else if (RejectBy == "kmeans") {
          w_retain <- RejectionBy(input = internal_vals,
                                  method = "kmeans")
        } else if (RejectBy == "lm") {
          w_retain <- RejectionBy(input = internal_vals,
                                  method = "lm")
        } else if (RejectBy == "direct") {
          w_retain <- RejectionBy(input = internal_vals,
                                  method = "direct")
        } else {
          w_retain <- seq(nrow(internal_vals))
        }
        Attr_PH[[Count]] <- internal_vals
        
        # if (RejectBy == "Rank") {
        #   return()
        #   z1 <- data.frame("local_PID" = vec1,
        #                    "local_Score" = vec2,
        #                    "approx_global_pid" = vec3,
        #                    # "approx_global_score" = vec4,
        #                    # "Delta_Background" = vec4 - diff3,
        #                    "response" = ifelse(test = SynExtendObject[[m1, m2]][, "group"] == 1L,
        #                                        yes = TRUE,
        #                                        no = FALSE))
        #   g <- glm(response ~ .,
        #            family = "quasibinomial",
        #            data = z1)
        #   
        #   FDR <- RejectionCriteria$FDR
        #   # FDR <- 0.005 # maximum allowed false discovery rate
        #   N <- sum(z1$response)
        #   
        #   ranking <- order(g$fitted.values,
        #                    decreasing = TRUE)
        #   ranking <- ranking[cumprod(cumsum(ranking > N) <= seq_along(ranking)*FDR) == 1L]
        #   w_retain <- sort(ranking[ranking <= N])
        #   # return(list("a" = z1,
        #   #             "b" = ranking,
        #   #             "c" = w_retain))
        #   # return(list("init" = z1,
        #   #             "rank" = sort(ranking[ranking <= N])))
        #   # # z2 <- z1[sort(ranking[ranking <= N]), ]
        #   # # return(list("init" = z1,
        #   # #             "retained" = z2))
        #   
        # } else if (RejectBy == "KMeans") {
        #   FDR <- RejectionCriteria$FDR
        #   
        #   z1 <- data.frame("Consensus" = diff2,
        #                    "KDist" = NucDist,
        #                    "local_PID" = vec1,
        #                    "local_Score" = vec2,
        #                    "approx_global_pid" = vec3,
        #                    "approx_global_score" = vec4,
        #                    "Delta_Background" = vec4 - diff3,
        #                    "response" = ifelse(test = SynExtendObject[[m1, m2]][, "group"] == 1L,
        #                                        yes = TRUE,
        #                                        no = FALSE))
        #   
        #   simplek <- suppressMessages(kmeans(x = z1[, -8L],
        #                                      centers = 10,
        #                                      nstart = 25))
        #   kres <- do.call(rbind,
        #                   tapply(X = as.factor(z1$response),
        #                          INDEX = simplek$cluster,
        #                          FUN = function(x) {
        #                            table(x)
        #                          }))
        #   kvals <- kres[, 1L] / rowSums(kres)
        #   kvals[is.infinite(kvals)] <- 1L
        #   w_kvals <- which(kvals < FDR)
        #   if (length(w_kvals) < 1) {
        #     w_retain <- seq(nrow(z1))
        #     return(list("k" = simplek,
        #                 "z" = z1,
        #                 "kres" = kres,
        #                 "kvals" = kvals,
        #                 "w_kvals" = w_kvals,
        #                 "fdr" = FDR,
        #                 "retain" = w_retain))
        #     warning("kmeans does not appear to return appropriate grouping for rejection, all candidates are being returned")
        #   } else {
        #     w_retain <- which(simplek$cluster %in% which(kvals < FDR))
        #   }
        #   
        #   
        #   
        #   # mimic-ing the old routine here doesn't seem to work as well, and i should
        #   # figure out why
        #   # continue <- TRUE
        #   # kres <- vector(mode = "list",
        #   #                length = 20L)
        #   # wss <- vector(mode = "numeric",
        #   #               length = length(kres))
        #   # count <- 2L
        #   # offsetindex <- 1L
        #   # offsetorigin <- 2L
        #   # calcstart <- 10L
        #   # select_clustering <- 0L
        #   # while (continue) {
        #   #   kres[[count - offsetindex]] <- suppressMessages(kmeans(x = z1[, -8L],
        #   #                                                          centers = count,
        #   #                                                          iter.max = 25L,
        #   #                                                          nstart = 25L))
        #   #   wss[count - offsetindex] <- kres[[count - offsetindex]]$tot.withinss
        #   #   startingBMax <- max(wss)
        #   #   startingHalfMax <- unname(quantile(wss[wss > 0], 0.5))
        #   #   if (count >= calcstart) {
        #   #     dat <- data.frame("n" = seq(count - offsetindex) - 1L,
        #   #                       "wss" = abs(wss[wss > 0] - wss[1]))
        #   #     
        #   #     
        #   #     fitval <- nls(formula = wss ~ OneSite(X = n,
        #   #                                           Bmax,
        #   #                                           Kd),
        #   #                   data = dat,
        #   #                   start = list(Bmax = max(dat$wss),
        #   #                                Kd = unname(quantile(dat$n, 0.25))))
        #   #     fitsum <- summary(fitval)
        #   #     
        #   #     return(list("dat" = dat,
        #   #                 "table" = z1,
        #   #                 "kmeans" = kres,
        #   #                 "fitval" = fitval,
        #   #                 "fitsum" = fitsum))
        #   #     # FitA <- nls(dat[, 2L]~OneSite(X = dat[, 1L],
        #   #     #                               Bmax,
        #   #     #                               Kd),
        #   #     #             start = list(Bmax = startingBMax,
        #   #     #                          Kd = startingHalfMax))
        #   #     # fitasum <- summary(FitA)
        #   #     return(fitasum)
        #   #   }
        #   #   count <- count + 1L
        #   # }
        # } else if (RejectBy == "None") {
        #   w_retain <- seq(nrow(SynExtendObject[[m1, m2]]))
        # }
        
        # block size determination
        if (length(w_retain) > 1) {
          blockres <- BlockByRank(index1 = IMatrix[w_retain, 1L],
                                  partner1 = PMatrix[w_retain, 1L],
                                  index2 = IMatrix[w_retain, 2L],
                                  partner2 = PMatrix[w_retain, 2L])
        } else if (length(w_retain) == 1) {
          blockres <- list("absblocksize" = 1L,
                           "blockidmap" = -1L)
        } else {
          blockres <- list("absblocksize" = vector(mode = "integer",
                                                   length = 0L),
                           "blockidmap" = vector(mode = "integer",
                                                 length = 0L))
        }
        block_ph <- blockres$blockidmap
        w1 <- block_ph > 0
        # if (any(w1)) {
        #   block_ph[w1] <- block_ph[w1] + block_offset
        #   block_offset <- block_offset + max(block_ph)
        # }
        if (any(w1)) {
          block_offset <- block_uid
          blockres$blockidmap[blockres$blockidmap > 0] <- blockres$blockidmap[blockres$blockidmap > 0] + block_offset
          block_uid <- max(blockres$blockidmap[blockres$blockidmap > 0])
        } else {
          # nothing to update or add
        }
        
        # BlockSize evaluation is complete
        # we're eventually moving blocksize stuff down to happen post everything else
        PH[[Count]] <- data.frame("p1" = names(QueryDNA)[PMatrix[w_retain, 1]],
                                  "p2" = names(SubjectDNA)[PMatrix[w_retain, 2]],
                                  "Consensus" = diff2[w_retain],
                                  "p1featurelength" = QueryFeatureLength[w_retain],
                                  "p2featurelength" = SubjectFeatureLength[w_retain],
                                  "blocksize" = blockres$absblocksize, # calculated after retention determination
                                  "KDist" = NucDist[w_retain],
                                  "TotalMatch" = TotalMatch[w_retain],
                                  "MaxMatch" = MatchMax[w_retain],
                                  "UniqueMatches" = UniqueMatches[w_retain],
                                  "Local_PID" = vec1[w_retain],
                                  "Local_Score" = vec2[w_retain],
                                  "Approx_Global_PID" = vec3[w_retain],
                                  "Approx_Global_Score" = vec4[w_retain],
                                  "Alignment" = ifelse(test = AASelect,
                                                       yes = "AA",
                                                       no = "NT")[w_retain],
                                  "Block_UID" = blockres$blockidmap,
                                  "Delta_Background" = (vec4 - diff3)[w_retain])
      } else {
        # link table is not populated
        PH[[Count]] <- data.frame("p1" = character(),
                                  "p2" = character(),
                                  "Consensus" = numeric(),
                                  "p1featurelength" = integer(),
                                  "p2featurelength" = integer(),
                                  "blocksize" = integer(),
                                  "KDist" = numeric(),
                                  "TotalMatch" = integer(),
                                  "MaxMatch" = integer(),
                                  "UniqueMatches" = integer(),
                                  "Local_PID" = numeric(),
                                  "Local_Score" = numeric(),
                                  "Approx_Global_PID" = numeric(),
                                  "Approx_Global_Score" = numeric(),
                                  "Alignment" = character(),
                                  "Block_UID" = integer(),
                                  "Delta_Background" = numeric())
      }
      # Count and PBCount are unlinked,
      # iterate through both separately and correctly
      Count <- Count + 1L
      
      if (object.size(DataPool) > Storage) {
        # ok ... 
        # i need to nuke positions in the pool based on a few things:
        # i can't nuke the current m1,
        # i can't nuke the next m1
        # i can't nuke the next m2
        # return(list("m1" = m1,
        #             "m2" = m2,
        #             "pool" = DataPool))
        sw1 <- sapply(X = DataPool,
                      FUN = function(x) {
                        !is.null(x)
                      })
        sw2 <- which(sw1)
        sw2 <- sw2[!(sw2 %in% c(m1, (m1 + 1L), (m2 + 1L)))]
        if (length(sw2) > 0) {
          # bonk the first one ... this might realistically need to happen in a while loop, but for now we can live with this
          # !!! assigning NULL by default deletes the position, shortening the vector,
          # we can replace the list (the container), with another container (containing NULL) without
          # shortening the list, R inferno reference and relevant examples here:
          # https://stackoverflow.com/questions/7944809/assigning-null-to-a-list-element-in-r
          DataPool[sw2[1L]] <- list(NULL)
        } else {
          # i don't know if this case can happen, but we're putting a print statement here just in case
          print("Please allocate more storage.")
        }
      }
      Prev_m1 <- m1
    } # end m2
  } # end m1
  res <- do.call(rbind,
                 PH)
  Attr_PH <- do.call(rbind,
                     Attr_PH)
  rownames(Attr_PH) <- NULL
  if (nrow(res) > 0) {
    # return(res)
    All_UIDs <- unique(res$Block_UID)
    res$Block_UID[res$Block_UID == -1L] <- seq(from = max(All_UIDs) + 1L,
                                               by = 1L,
                                               length.out = sum(res$Block_UID == -1L))
  }
  attr(x = res,
       which = "GeneCalls") <- GeneCalls
  class(res) <- c("data.frame",
                  "PairSummaries")
  attr(x = res,
       which = "AlignmentFunction") <- AlignmentFun
  attr(x = res,
       which = "KVal") <- KmerSize
  attr(x = res,
       which = "AA_matrix") <- AA_matrix
  attr(x = res,
       which = "NT_matrix") <- NT_matrix
  if (RetainInternal) {
    attr(x = res,
         which = "internal_vals") <- Attr_PH
  }
  # attr(x = res,
  #      which = "NT_Anchors") <- NT_Anchors
  
  # close pBar and return res
  if (Verbose) {
    TimeEnd <- Sys.time()
    close(pBar)
    print(TimeEnd - TimeStart)
  }
  return(res)
}
