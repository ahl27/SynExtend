\name{EvoWeaver}
\alias{EvoWeaver}
\alias{EvoWeaver-class}
\alias{EvoWeaver-utils}
\alias{SpeciesTree}
\alias{SpeciesTree.EvoWeaver}
\title{
EvoWeaver: Identifying Gene Functional Associations from Coevolutionary Signals
}
\description{
EvoWeaver is an S3 class with methods for predicting functional association using
protein or gene data. EvoWeaver implements multiple algorithms for analyzing coevolutionary signal between genes, which are combined into overall predictions on functional association. For details on predictions, see \code{\link{predict.EvoWeaver}}.
}
\usage{
EvoWeaver(ListOfData, MySpeciesTree=NULL, NoWarn=FALSE)

\method{SpeciesTree}{EvoWeaver}(ew, Verbose=TRUE, ...)
}
\arguments{
  \item{ListOfData}{
A list of gene data, where each entry corresponds to information on a particular
gene. List must contain either dendrograms or vectors, and cannot contain a mixture.
If list is composed of dendrograms, each dendrogram is a gene tree for the corresponding entry. If list is composed of vectors, vectors should be numeric or character vectors denoting the genomes containing that gene.
}
\item{MySpeciesTree}{
An object of class \code{'dendrogram'} representing the overall species tree for
the list provided in \code{ListOfData}.
}
\item{NoWarn}{
Logical; If \code{FALSE}, displays warnings corresponding to which algorithms are unavailable for given input data format (see Details for more information).
}
\item{ew}{
An object of class \code{EvoWeaver}.
}
\item{Verbose}{
Logical; If \code{TRUE}, displays output when calculating reference tree.
}
\item{...}{
Further arguments passed to \code{\link{SuperTree}} for inferring a reference tree.
}
}
\details{
EvoWeaver expects input data to be a list. All entries must be one of the following cases:
\enumerate{
  \item{
    \code{ListOfData[[i]] = c('ID#1', 'ID#2', ..., 'ID#k')}
  }
  \item{
    \code{ListOfData[[i]] = c('g1_d1_s1_p1', 'g2_d2_s2_p2', ..., 'gk_dk_sk_pk') }}

  \item{
  \code{ListOfData[[i]] = dendrogram(...)}
  }
}
In (1), each ID#i corresponds to the unique identifier for genome #i. For entry #j in the list, the presence of 'ID#i' means genome #i has an ortholog for gene/protein #j.

Case (2) is the same as (1), just with the formatting of names slightly different. Each entry is of the form \code{g_d_p}, where \code{g} is the unique identifier for the genome, \code{d} is which chromosome the ortholog is located, \code{s} indicates whether the gene is on the forward or reverse strand, and \code{p} is what position the ortholog appears in on that chromosome. \code{p} must be a \code{numeric}. \code{s} must be \code{0} or \code{1}, corresponding to whether the gene is on the forward or reverse strand. Whether \code{0} denotes forward or reverse is inconsequential as long as the scheme is consistent. \code{g,d} can be any value as long as they don't contain an underscore (\code{'_'}).

Case (3) expects gene trees for each gene, with labeled leaves corresponding to each source genome. If \code{ListOfData} is in this format, taking \code{labels(ListOfData[[i]])} should produce a character vector that matches the format of one of the previous cases.

\emph{See the Examples section for illustrative examples.}

\emph{Whenever possible, provide a full set of \code{dendrogram} objects with leaf labels in form (2). This will allow the most algorithms to run. What follows is a more detailed description of which inputs allow which algorithms.}

EvoWeaver requires input of scenario (3) to use distance matrix methods, and requires input of scenario (2) (or (3) with leaves labeled according to (2)) for gene organization analyses. Sequence Level methods require dendrograms with sequence information included as the \code{state} attribute in each leaf node.

Note that ALL entries must belong to the same category--a combination of character vectors and dendrograms is not allowed.

Prediction of a functional association network is done using \code{predict(EvoWeaverObject)}. See \code{\link{predict.EvoWeaver}} for more information.

The \code{SpeciesTree} function takes in an object of class \code{EvoWeaver} and returns a species tree. If the object was not initialized with a species tree, it calculates one using \code{\link{SuperTree}}. The species tree for a EvoWeaver object can be set with \code{attr(ew, 'speciesTree') <- ...}.
}
\value{
Returns a EvoWeaver object.
}
\author{
Aidan Lakshman \email{ahl27@pitt.edu}
}
\seealso{
\code{\link{predict.EvoWeaver}},
\code{\link{ExampleStreptomycesData}},
\code{\link{BuiltInEnsembles}},
\code{\link{SuperTree}}
}
\examples{
# I'm using gene to mean either a gene or protein

## Imagine we have the following 4 genomes:
## (each letter denotes a distinct gene)
##    Genome 1: a b c d
##    Genome 2: d c e
##    Genome 3: b a e
##    Genome 4: a e

## We have 5 total genes: (a,b,c,d,e)
##    a is present in genomes 1, 3, 4
##    b is present in genomes 1, 3
##    c is present in genomes 1, 2
##    d is present in genomes 1, 2
##    e is present in genomes 2, 3, 4

## Constructing a EvoWeaver object according to (1):
l <- list()
l[['a']] <- c('1', '3', '4')
l[['b']] <- c('1', '3')
l[['c']] <- c('1', '2')
l[['d']] <- c('1', '2')
l[['e']] <- c('2', '3', '4')

## Each value of the list corresponds to a gene
## The associated vector shows which genomes have that gene
pwCase1 <- EvoWeaver(l)

## Constructing a EvoWeaver object according to (2):
##  Here we need to add in the genome, chromosome, direction, and position
##  As we only have one chromosome,
##  we can just set that to 1 for all.
##  Position can be identified with knowledge, or with
##  FindGenes(...) from DECIPHER.

## In this toy case, genomes are small so it's simple.
l <- list()
l[['a']] <- c('a_1_0_1', 'c_1_1_2', 'd_1_0_1')
l[['b']] <- c('a_1_1_2', 'c_1_1_1')
l[['c']] <- c('a_1_1_3', 'b_1_0_2')
l[['d']] <- c('a_1_0_4', 'b_1_0_1')
l[['e']] <- c('b_1_0_3', 'c_1_0_3', 'd_1_0_2')

pwCase2 <- EvoWeaver(l)

## For Case 3, we just need dendrogram objects for each
# l[['a']] <- dendrogram(...)
# l[['b']] <- dendrogram(...)
# l[['c']] <- dendrogram(...)
# l[['d']] <- dendrogram(...)
# l[['e']] <- dendrogram(...)

## Leaf labels for these will be the same as the
##  entries in Case 1.
}
