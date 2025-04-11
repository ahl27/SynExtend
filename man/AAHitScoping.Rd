\name{AAHitScoping}
\alias{AAHitScoping}
\title{
Adjust the scope of kmer hits between feature and genome space.
}
\description{
This function is designed to work internally to functions within \code{SynExtend} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
AAHitScoping(hitlist,
             fstrand1,
             fstart1,
             fstop1,
             fstrand2,
             fstart2,
             fstop2)
}
\arguments{
  \item{hitlist}{
  A list containing matrices produced by \code{\link{SearchIndex}}.
}
  \item{fstrand1}{
  An integer vector of 0s and 1s describing the strand of features.
}
  \item{fstart1}{
  Integer; a vector of left bounds of features.
}
  \item{fstop1}{
  Integer; a vector of right bounds of features.
}
  \item{fstrand2}{
  An integer vector of 0s and 1s describing the strand of features.
}
  \item{fstart2}{
  Integer; a vector of left bounds of features.
}
  \item{fstop2}{
  Integer; a vector of right bounds of features.
}
}
\details{
  \code{AAHitScoping} converts the hits returned by \code{\link{SearchIndex}} from feature-to-feature context genome-to-genome context.
}
\value{
  A list of matrices.
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
