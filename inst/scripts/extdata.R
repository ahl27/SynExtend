###### -- External data for SynExtend  ----------------------------------------
# Author: Nicholas Cooley
# email: npc19@pitt.edu
# this pseudocode relies on the ncbi command line utilities, they can be found
# here: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# they must be installed, and R must have access to the executables
# this data was generated on 2022 04 21
# Self note:
# UPDATE DATALIST MANUALLY
# data was regenerated on 2022 09 22

# data was regenerated again on 2024 04 25 for a set of new functions and
# the deprecation of some old ones.
# data was regenerated again on 2025 03 10 for some function adjustments
# source(file = "~/Packages/SynExtend/R/SummarizePairs.R", echo = FALSE)
# source(file = "~/Packages/SynExtend/R/ClusterByK.R", echo = FALSE)
# source(file = "~/Packages/SynExtend/R/ExpandDiagonal.R", echo = FALSE)

suppressMessages(library(SynExtend))
suppressMessages(library(RSQLite))

# data was regenerated again on 2025 02 08 for a new function
source(file = "~/Packages/SynExtend/R/PrepareSeqs.R", echo = FALSE)
source(file = "~/Packages/SynExtend/R/SummarizePairs.R", echo = FALSE)

TODAYSDATE <- paste0(unlist(strsplit(x = as.character(Sys.time()),
                                     split = "-| ")[[1]][1:3]),
                     collapse = "")

###### -- Entrez --------------------------------------------------------------

EntrezQuery <- paste("esearch -db assembly ",
                     "-query '",
                     "endosymbiont[All Fields] ",
                     'AND "complete genome"[filter] ',
                     'AND "RefSeq has annotation"[properties] ',
                     "NOT anomalous[filter]' ",
                     '| ',
                     'esummary ',
                     '| ',
                     'xtract -pattern DocumentSummary -element FtpPath_RefSeq',
                     sep = "")

FtPPaths <- system(command = EntrezQuery,
                   intern = TRUE,
                   timeout = 300L)

# keep the example data small...
if (length(FtPPaths) > 4L) {
  # setting the seed here is a little superfluous, as the data size will almost certainly
  # change before this is run again, or before anyone comes back to run it themselves
  set.seed(1986)
  FtPPaths <- sample(FtPPaths,
                     size = 4,
                     replace = FALSE)
}

adds <- mapply(SIMPLIFY = TRUE,
               USE.NAMES = FALSE,
               FUN = function(x, y) {
                 paste0(x,
                        "/",
                        y[10],
                        c("_genomic.fna.gz",
                          "_genomic.gff.gz",
                          "_protein.faa.gz"))
               },
               x = FtPPaths,
               y = strsplit(x = FtPPaths,
                            split = "/",
                            fixed = TRUE))
fnas <- adds[1, , drop = TRUE]
gffs <- adds[2, , drop = TRUE]
amns <- adds[3, , drop = TRUE]


###### -- Data import ---------------------------------------------------------
# save off gffs as external non-R data
# save off `GeneCalls` as an object for examples

# save off one GFF for `gffToDataFrame's example`
CURLCOMMAND <- paste0("curl --output ",
                      paste0("~/Packages/SynExtend/inst/extdata/",
                             unlist(regmatches(x = gffs[1],
                                               m = gregexpr(pattern = "[^/]+\\.gff\\.gz",
                                                            text = gffs[1])))),
                      " ",
                      gffs[1])

system(command = CURLCOMMAND,
       intern = FALSE)

Endosymbionts_GeneCalls <- vector(mode = "list",
                                    length = length(gffs))

VignetteDB01 <- "~/Packages/SynExtend/inst/extdata/Endosymbionts_v05a.sqlite"
VignetteDB02 <- tempfile()

for (m1 in seq_along(gffs)) {
  Endosymbionts_GeneCalls[[m1]] <- gffToDataFrame(GFF = gffs[m1],
                                                    Verbose = TRUE)
  Seqs2DB(seqs = fnas[m1],
          type = "FASTA",
          dbFile = VignetteDB01,
          identifier = as.character(m1),
          verbose = TRUE)
}

names(Endosymbionts_GeneCalls) <- seq(length(Endosymbionts_GeneCalls))

Endosymbionts_Synteny <- FindSynteny(dbFile = VignetteDB01,
                                       verbose = TRUE)

save(Endosymbionts_Synteny,
     file = "~/Packages/SynExtend/data/Endosymbionts_Synteny.RData",
     compress = "xz")

save(Endosymbionts_GeneCalls,
     file = "~/Packages/SynExtend/data/Endosymbionts_GeneCalls.RData",
     compress = "xz")

###### -- NucleotideOverlap ---------------------------------------------------

Endosymbionts_LinkedFeatures <- NucleotideOverlap(SyntenyObject = Endosymbionts_Synteny,
                                                  GeneCalls = Endosymbionts_GeneCalls,
                                                  Verbose = TRUE)

save(Endosymbionts_LinkedFeatures,
     file = "~/Packages/SynExtend/data/Endosymbionts_LinkedFeatures.RData",
     compress = "xz")

###### -- PrepareSeqs ---------------------------------------------------------

system(command = paste("cp",
                       VignetteDB01,
                       VignetteDB02))

PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
            DataBase = VignetteDB02,
            Verbose = TRUE)

###### -- PairSummaries -------------------------------------------------------

CONN01 <- dbConnect(SQLite(), VignetteDB02)

Endosymbionts_Pairs01 <- SummarizePairs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                        DataBase = CONN01,
                                        Verbose = TRUE)

save(Endosymbionts_Pairs01,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs01.RData",
     compress = "xz")

###### -- Clustering ----------------------------------------------------------

Endosymbionts_Pairs02 <- ClusterByK(SynExtendObject = Endosymbionts_Pairs01,
                                    ClusterScalar = 5,
                                    ShowPlot = TRUE,
                                    Verbose = TRUE)

save(Endosymbionts_Pairs02,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs02.RData",
     compress = "xz")

###### -- BlockReconciliation -------------------------------------------------

Endosymbionts_Pairs03 <- ExpandDiagonal(SynExtendObject = Endosymbionts_Pairs02[Endosymbionts_Pairs02$ClusterID %in% as.integer(names(which(attr(x = Endosymbionts_Pairs02,
                                                                                                                                                 which = "Retain")))), ],
                                        DataBase = CONN01,
                                        Verbose = TRUE)
save(Endosymbionts_Pairs03,
     file = "~/Packages/SynExtend/data/Endosymbionts_Pairs03.RData",
     compress = "xz")

###### -- DisjointSet ---------------------------------------------------------

Endosymbionts_Sets <- DisjointSet(Pairs = Endosymbionts_Pairs03,
                                  Verbose = TRUE)

save(Endosymbionts_Sets,
     file = "~/Packages/SynExtend/data/Endosymbionts_Sets.RData",
     compress = "xz")

###### -- ExtractBy ----------------------------------------------------------- 
# no functions in the pipeline beyond this function
# no need to save off this for examples
# Endosymbionts_Gene_Communities <- ExtractBy(x = Endosymbionts_Pairs03,
#                                             y = VignetteDB02,
#                                             z = Endosymbionts_Sets,
#                                             Verbose = TRUE)
# 
# save(Endosymbionts_Gene_Communities,
#      file = "~/Packages/SynExtend/data/Endosymbionts_Gene_Communities.RData",
#      compress = "xz")

###### -- CompetePairs --------------------------------------------------------

# no functions in the pipeline beyond this function
# so we don't need to save this for now...
# Endosymbionts_Pairs04 <- CompetePairs(SynExtendObject = Endosymbionts_Pairs01,
#                                       Verbose = TRUE)
# save(Endosymbionts_Pairs04,
#      file = "~/Packages/SynExtend/data/Endosymbionts_Pairs04.RData",
#      compress = "xz")



