\name{WithinSetCompetition}
\alias{WithinSetCompetition}
\title{
Pare down candidate pairs to one-to-one sets.
}
\description{
This function is a work in progress, please be patient.
}
\usage{
WithinSetCompetition(SynExtendObject,
                     AllowCrossContigConflicts = TRUE,
                     CompeteBy  = "Delta_Background",
                     PollContext = TRUE,
                     ContextInflation = 0.975,
                     Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
  A \code{PairSummaries} object created by \code{\link{SummarizePairs}}.
}
  \item{AllowCrossContigConflicts}{
  Logical; return only one candidate per per disjoint set for each contig to contig pair, or for each genome to genome pair.
}
  \item{CompeteBy}{
  Character; a column name from the \code{PairSummaries} object.
}
  \item{PollContext}{
  Logical; when competing candidate pairs, consider block membership.
}
  \item{ContextInflation}{
  Numeric; a value to adjust block membership strength. Lower values increase the strength of block membership.
}
  \item{Verbose}{
  Logical; print a progress bar and timings.
}
}
\details{
  For each assembly to assembly comparison, or each contig to contig comparison, all disjoint sets are collected for candidate pairs. In cases where there are more than 2 nodes, i.e. features, in a set, the connecting edges are competed against each other, and only the strongest edges -- and their resulting nodes -- are retained.
}
\value{
  A \code{PairSummaries} object.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{SummarizePairs}}
}
\examples{
#
}
