\name{RejectionBy}
\alias{RejectionBy}
\title{
Given an object of candidate pairs, reject false positives.
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
RejectionBy(input,
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
            dropinappropriate = FALSE)
}
\arguments{
  \item{input}{
  A data.frame; currently set up to take in an internal representation of the columns eventually printed by \code{\link{SummarizePairs}}. Supported column names are enumerated in the \code{supportedcolumns} argument.
}
  \item{criteria}{
  A list of named objects that control the nobs of various rejection routines.
}
  \item{rankby}{
  A colname in the the \code{input} data.frame which will be used for the identification of false positives.
}
  \item{method}{
  Character; identify a method by which to reject false positives, currently supported methods include: \code{glm}, \code{lm}, \code{kmeans}, and \code{direct}.
}
  \item{supportedcolumns}{
  Character; a vector of column names to select internal variables for candidate pair evaluation.
}
  \item{dropinappropriate}{
  Logical; \code{FALSE} by default. If \code{TRUE}, should a method imply that it is incapable of logically separating true positives from false positives, no candidates are returned. If \code{FALSE}, all candidates are returned. Only really applicable in cases where very few initial candidates are identified.
}
}
\details{
  \code{RejectionBy} is yet another attempt at building a logical but simple routine for rejecting false positive candidate pairs in a scenario where it is difficult to proxy a distribution of appropriate false-positives.
}
\value{
  A value.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}, \code{\link{WithinSetCompetition}}, \code{\link{RejectionBy}}
}
\examples{
#
}
