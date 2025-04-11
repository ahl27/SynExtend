\name{ApproximateBackground}
\alias{ApproximateBackground}
\title{
Return the approximate background alignment score for a series of paired sequences.
}
\description{
This function is designed to work internally to \code{\link{SummarizePairs}} so it works on relatively simple atomic vectors and has little overhead checking.
}
\usage{
ApproximateBackground(p1,
                      p2,
                      code1,
                      code2,
                      mod1,
                      mod2,
                      aa1,
                      aa2,
                      nt1,
                      nt2,
                      register1,
                      register2,
                      aamat,
                      ntmat)
}
\arguments{
  \item{p1}{
  Integer; references positions within nt1 or aa1.
}
  \item{p2}{
  Integer; references positions within nt2 or aa2.
}
  \item{code1}{
  Logical; specifies whether the position referenced by p1 is reported as a coding sequence.
}
  \item{code2}{
  Logical; specifies whether the position referenced by p2 is reported as a coding sequence.
}
  \item{mod1}{
  Logical; specifies whether the position referenced by p1 can be translated without complaint by \code{\link{translate}}.
}
  \item{mod2}{
  Logical; specifies whether the position referenced by p2 can be translated without complaint by \code{\link{translate}}.
}
  \item{aa1}{
  AAStringSet.
}
  \item{aa2}{
  AAStringSet.
}
  \item{nt1}{
  DNAStringSet.
}
  \item{nt2}{
  DNAStringSet.
}
  \item{register1}{
  Integer; a vector that maps which positions in aa1 are the translations of that particular index in nt1. NAs identify positions that are not translated.
}
  \item{register2}{
  Integer; a vector that maps which positions in aa2 are the translations of that particular index in nt2. NAs identify positions that are not translated.
}
  \item{aamat}{
  A substitution matrix for amino acids.
}
  \item{ntmat}{
  A substitution matrix for nucleotides.
}
}
\details{
  \code{ApproximateBackground} generates approximate background alignment scores for sets of sequences.
}
\value{
  A vector of numerics.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{NucleotideOverlap}}, \code{\link{SummarizePairs}}, \code{\link{FindSynteny}}
}
\examples{
#
}
