\name{SummarizePairs}
\alias{SummarizePairs}
\title{
Provide summaries of hypothetical orthologs.
}
\description{
Given a \code{LinkedPairs} object and a \code{DECIPHER} database, return a data.frame of summarized genomic feature pairs. \code{SummarizePairs} will collect all the linked genomic features in the supplied \code{\link{LinkedPairs-class}} object and return descriptions of the alignments of those features.
}
\usage{
SummarizePairs(SynExtendObject,
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
               ...)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{LinkedPairs-class}.
}
  \item{DataBase01}{
A character string pointing to a SQLite database, or a connection to a \code{DECIPHER} database.
}
  \item{AlignmentFun}{
  Character of length 1; a character string of length one specifying a \code{link{DECIPHER}} alignment function. Currently only supports \code{\link{AlignPairs}}.
}
  \item{DefaultTranslationTable}{
  Character of length 1; an identifier that can be recognized by \code{\link{getGeneticCode}} to use as the translation table for translating coding sequences in the case that one is missing from supplied genecalls.
}
  \item{KmerSize}{
  Integer of length 1; Specify the kmer size for assessing kmer distance in nucleotide space between two candidate pairs.
}
  \item{Verbose}{
Logical of length 1; if \code{TRUE} progress bar and function timing will be displayed.
}
  \item{ShowPlot}{
Logical of length 1; if \code{TRUE} provide some plots describing candidate pairs. Currently not implemented.
}
  \item{Processors}{
  Integer of length 1; specify the number of processors available to \code{SummarizePairs} for multithreaded applications. If \code{NULL} all available detectable cores will be requested.
}
  \item{Storage}{
  Numeric of length 1; a soft limit on the memory alloted to \code{SummarizePairs} for the storage of sequence data from the supplied database. In Gb.
}
  \item{IndexParams}{
  A named list of arguments to be passed to \code{\link{IndexSeqs}}. Must be compliant with \code{\link{do.call}}'s expectation for its \code{args} argument.
}
  \item{SearchParams}{
  A named list of arguments to be passed to \code{\link{SearchIndex}}. Must be compliant with \code{\link{do.call}}'s expectation for its \code{args} argument.
}
  \item{SearchScheme}{
  Character of length 1; currently supported arguments include; "spike" indicating to 'spike' in a population of background candidates by searching one set of codings sequences against the reverse of another, "standard" which will only search coding sequences from one genome against the other in the forward direction, and "reciprocal" which will perform a search strategy similar to Reciprocal Best Hits.
}
  \item{RejectBy}{
  Character of length 1; currently supported arguments include; "glm" and "lm" which use the eponymous functions to model the data within a set of candidate pairs and reject candidate pairs below a particular False Discovery Rate as determined from a set of known negatives generated when a "spike" search scheme is used. "kmeans" is a supported method that will run a naive kmeans based routine to cluster candidates within the set and reject candidates below a user supplied threshold. Lastly, "direct" will simply rank all candidate pairs by the user supplied attribute and drop all candidates below a user supplied FDR threshold.
}
  \item{RetainInternal}{
  Logical of length 1; if \code{TRUE} internal values used for candidate pair rejection will be attached to the returned object.
}
  \item{...}{
Additional arguments to pass to interior functions. Currently not implemented.
}
}
\details{
\code{SummarizePairs} collects features describing each linked feature pair. These features include:
  \itemize{
    \item{
    p1: a character identifier for the candidate pair partner in the supplied query.
    }
    \item{
    p2: a character identifier for the candidate pair partner in the supplied subject.
    }
    \item{
    Consensus: a numeric value calculated by \code{\link{HitConsensus}} describing whether relative locations of linking hits are in linearly similar positions in both candidate pair partners.
    }
    \item{
    p1featurelength: length of candidate query feature in nucleotides.
    }
    \item{
    p2featurelength: length of candidate subject feature in nucleotides.
    }
    \item{
    blocksize: integer value indicating the number of shared features blocked together.
    }
    \item{
    KDist: numeric value of the euclidean distance between candidate pairs in kmer space.
    }
    \item{
    TotalMatch: integer value indicating total nucleotides shared between candidates pairs in the original searches.
    }
    \item{
    MaxMatch: integer value indicating the largest kmer shared between candidate pairs in the original searches.
    }
    \item{
    UniqueMatches: integer value indicating the number of kmers shared between candidate pair partners.
    }
    \item{
    Local_PID: numeric value of the local PID for the alignment of the candidate pair.
    }
    \item{
    Local_Score: numeric value of the local alignment score for the candidate pair.
    }
    \item{
    Approx_Global_PID: approximate global PID for the alignment of the candidate pair.
    }
    \item{
    Approx_Global_Score: approximate global score for the alignment of the candidate pair.
    }
    \item{
    Block_UID: integer value giving an identifier number to the feature block that that candidate pair is part of.
    }
    \item{
    Delta_Background: the approximate global score of the alignment modified by the background of the sequences.
    }
  }
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{PrepareSeqs}}, \code{\link{NucleotideOverlap}}, \code{\link{FindSynteny}}, \code{\link{LinkedPairs-class}}
}
\examples{
library(RSQLite)
DBPATH <- system.file("extdata",
                      "Endosymbionts_v05a.sqlite",
                      package = "SynExtend")
tmp01 <- tempfile()
file.copy(from = DBPATH,
          to = tmp01)
data("Endosymbionts_LinkedFeatures", package = "SynExtend")
PrepareSeqs(SynExtendObject = Endosymbionts_LinkedFeatures,
            DataBase = tmp01,
            Verbose = TRUE)
DBCONN <- dbConnect(SQLite(), tmp01)
data("Endosymbionts_LinkedFeatures", package = "SynExtend")
SummarizedPairs <- SummarizePairs(SynExtendObject = Endosymbionts_LinkedFeatures,
                                  DataBase01 = DBCONN,
                                  Verbose = TRUE)
dbDisconnect(DBCONN)
}
