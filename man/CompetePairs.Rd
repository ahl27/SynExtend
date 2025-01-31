\name{CompetePairs}
\alias{CompetePairs}
\title{
Find the best match pair in cases where ambiguity exists.
}
\description{
A relatively simple routine for identifying a ``best'' pair in cases where many homologous are identified in a single genome-to-genome comparison. Selection is performed with a single collected measure, and can be performed with or without leveraging context of syntenic blocks.
}
\usage{
CompetePairs(SynExtendObject,
             AllowCrossContigConflicts = TRUE,
             By = "PID",
             PollContext = TRUE,
             NormalizeCompetition = TRUE,
             InflationParameter = .975,
             Verbose = FALSE)
}
\arguments{
  \item{SynExtendObject}{
An object of class \code{PairSummaries}.
}
  \item{AllowCrossContigConflicts}{
A logical indicating where pair competition should take place between genomes or contigs.
}
  \item{By}{
A character vector of length 1 indicating which column in the \code{PairSummaries} object to compete pairs with.
}
  \item{PollContext}{
A logical indicating whether to include the context of block membership in the competition.
}

  \item{NormalizeCompetition}{
A logical indicating whether or not to unit normalize the measure being used for competition.
}
  \item{InflationParameter}{
A numeric of length 1 specifying an adjustment for how block context should be penalized, greater than 1 or rewarded, less than 1.
}
  \item{Verbose}{
Logical indicating whether or not to display a progress bar and print the time difference upon completion.
}
}
\details{
\code{CompetePairs} uses a naive competition based approach to select a ``true-est'' ortholog in cases where many competing potential orthologs are present in a set of predicted pairs. The returned value is the previous object with two new attributes, ``RetainByCompetition'' is a vector of logicals specifying which pairs are retained post competition. A second new attribute named ``Knockout'' is a character vector that identifies --by rowname-- the row which knocked out a removed pair.
}
\value{
An object of class \code{PairSummaries}.
}
\author{
Nicholas Cooley \email{npc19@pitt.edu}
}

\seealso{
\code{\link{SummarizePairs}}, \code{\link{NucleotideOverlap}}, \code{\link{FindSynteny}}, \code{\link{ExpandDiagonal}}, \code{\link{ClusterByK}}
}
\examples{
data("Endosymbionts_Pairs01", package = "SynExtend")

Pairs02 <- CompetePairs(SynExtendObject = Endosymbionts_Pairs01)
}
