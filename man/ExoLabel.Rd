\name{ExoLabel}
\alias{ExoLabel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ExoLabel: Out of Memory Fast Label Propagation
}
\description{
Detects communities in networks with Fast Label Propagation using disk space to drastically reduce memory complexity.
}
\usage{
ExoLabel(edgelistfiles,
              outfile=tempfile(),
              mode=c("undirected", "directed"),
              add_self_loops=FALSE,
              attenuation=TRUE,
              ignore_weights=FALSE,
              iterations=0L,
              return_table=FALSE,
              consensus_cluster=FALSE,
              use_fast_sort=FALSE,
              verbose=interactive(),
              sep='\t',
              tempfiledir=tempdir())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{edgelistfiles}{
Character; vector of files to be processed. Each entry should be a machine-interpretable path to an edgelist file. See Details for expected format.
}

\item{outfile}{
Character; File to write final clusters to. Can be set to a vector of filepaths to run multiple clusterings (see Details).
}

\item{mode}{
Character; Specifies whether edges should be interpreted as undirected (default) or directed. Can be "undirected", "directed", or an unambiguous abbreviation.
}

\item{add_self_loops}{
Logical or Numeric; Determines if self loops should be added to the network. If \code{TRUE}, adds self loops of weight 1.0 to all vertices. If set to numeric value \code{w}, adds self loops of weight \code{w} to all nodes. Can also be set to a vector when running multiple clusterings (see Details).
}

\item{attenuation}{
Logical or Numeric; Determines if label-hop attenuation should be used. If \code{TRUE}, uses attenuation to prevent single clusters from dominating results. Can also be set to a numeric to influence the strength of attenuation. See "Algorithm Convergence" for more information on this parameter. Can also be set to a vector when running multiple clusterings (see Details).
}

\item{ignore_weights}{
Logical; Determines if weights should be ignored. If \code{TRUE}, all edges will be treated as an edge of weight 1. Must be set to \code{TRUE} if any of \code{edgelistfiles} are two-column tables (start->end only, lacking a weights column).
}

\item{iterations}{
Integer; Maximum number of times to process each node. If set to zero or \code{NULL}, automatically uses the square root of the max node degree. See "Algorithm Convergence" for more information.
}

\item{return_table}{
Logical; Determines how the result of clustering is returned. If \code{FALSE}, returns a character vector corresponding to the path of \code{outfile}. If \code{TRUE}, parses \code{outfile} using \code{\link{read.table}} and returns the result (not recommended for very large graphs).
}

\item{consensus_cluster}{
Logical or Numeric; Determines if consensus clustering should be used. If \code{TRUE}, runs the clustering algorithm multiple times and forms a consensus clustering based on the agreement of each run. Can be set to a numeric vector to control the number of iterations. See "Consensus Clustering" below for more information.
}

\item{use_fast_sort}{
Logical; Determines how files should be sorted. If \code{FALSE}, ExoLabel will perform file sorting functions in-place. If \code{TRUE}, ExoLabel will perform its file sorting functions using a second temporary file. This is much faster than the in-place sort, but consumes twice the amount of disk space. The relative disk consumption is about the same size as the input graph for \code{use_fast_sort=FALSE}, and about double the size of the input graph for \code{use_fast_sort=TRUE} (see "Memory Consumption" and the last paragraph of "Warning" below). Set to \code{TRUE} if you're not worried about disk utilization.
}

\item{verbose}{
Logical; Determines if status messages (output, progress, etc.) should be displayed while running. Output messages are reduced if running in non-interactive mode.
}

\item{sep}{
Character; expected character that separates entries on a line in each file in \code{edgelistfiles}. Defaults to tab, as would be expected in a \code{.tsv} formatted file. Set to \code{','} for a \code{.csv} file.
}
\item{tempfiledir}{
Character; Vector corresponding to the location where temporary files used during execution should be stored.
}
}
\details{
ExoLabel identifies communities (clusters) in graph/network structures using a variant of Fast Label Propagation, as proposed by Traag and Subelj (2023).

However, very large graphs require too much RAM for processing on some machines. In a graph containing billions of nodes and edges, loading the entire structure into RAM is rarely feasible. ExoLabel uses disk space for storing representations of graphs. While this is slower than computing on RAM, it allows ExoLabel to scale to graphs of enormous size while only using a comparatively small amount of memory. See "Memory Consumption" for details on the total disk/memory consumption of ExoLabel.

ExoLabel expects a set of edgelist files, provided as a vector of filepaths. Each entry in the file is expected to be in the following:

\code{VERTEX1<sep>VERTEX2<sep>WEIGHT<linesep>}

This line defines a single edge between vertices \code{VERTEX1} and \code{VERTEX2} with weight \code{WEIGHT}. \code{VERTEX1} and \code{VERTEX2} are strings corresponding to vertex names, \code{WEIGHT} is a numeric value that can be interpreted as a \code{double}. The separator \code{<sep>} corresponds to the argument \code{sep} (defaulting to tab for \code{.tsv} format), and \code{linesep} is the newline value \code{'\n'}.

If \code{ignore_weight=TRUE}, the file can be formatted as:

\code{VERTEX1<sep>VERTEX2<linesep>}

Note that the \code{v1 v2 w} format is still accepted for \code{ignore_weight=FALSE}, but the weights will be ignored.
}

\section{Algorithm Convergence}{
One of the main issues of Label Propagation algorithms is that they can fail to converge. Consider an unweighted directed graph with four nodes connected in a loop. That is, \code{A->B, B->C, C->D, D->A}. If \code{A,C} are in cluster 1 and \code{B,D} are in cluster 2, this algorithm could keep processing all the nodes in a loop and never converge. To solve this issue, we introduce an additional measure for convergence controlled by \code{iterations}. If \code{iterations=x}, then we only allow the algorithm to process each node \code{x} times. Once a given node has been seen \code{x} times, it is no longer updated. This can be manually specified, but defaults to the square root of the largest node indegree.

Additionally, ExoLabel incorporates label hop attenuation to reduce the chance of a single massive cluster dominating  results. In short, as a particular label propagates to other nodes, its subsequent contribution diminishes. The farther a particular label is from its original source, the less its contribution. The degree to which its contribution diminishes scales dynamically based on the proportion of nodes that update on each cycle. Each node's attenuated weight is calculated as \eqn{w' = (1-p^ad)w}, where \eqn{w} the node weight, \eqn{p} the proportion of nodes that changed label in the previous iteration, \eqn{a} the attenuation power, and \eqn{d} the distance from the initial label. Larger values of \code{attenuation} create larger clusters, whereas smaller values create smaller clusters. The default value of \code{TRUE} (equivalent to 1.0) recovers the original implementation provided in Leung et al. (2009) linked in the References section. Passing a value of 0 (or \code{FALSE}) disables attenuation entirely rather than returning all singleton clusters.
}

\section{Multiple Clusterings}{
A large portion of the processing time is reading in the graph object. This leads to a lot of duplicated effort when trying to cluster the same network with multiple parameters.

Multiple clusterings on the same network are supported by passing vectors of input to \code{outfile} and \code{add_self_loops} or \code{attenuation}. If the length of \code{outfile} is greater than 1, \code{add_self_loops} and \code{attenuation} can each be set to either a single value or a vector of the same length as \code{outfile}. For a single value, the same parameter value will be used across all clusterings. For multiple values, the corresponding value will be used in each clustering. See "Examples" for example usage.

Note that the order to process each node is randomly initialized, so multiple runs on the same parameters may produce different results if a random seed is not set.
}

\section{Consensus Clustering}{
Consensus clustering can be enabled by setting \code{consensus_cluster=TRUE}. Consensus clustering runs ExoLabel on the input graph multiple times, transforming weight values according to a sigmoid function. By default, this runs nine times for sigmoids with scale 0.5 and shape \code{c(0,0.2,0.4,0.6,0.8,1.0,1.33,1.67,2.0)}, collapsing weights below 0.1 to zero. The resulting clusters form a network such that the edge weight between any two nodes connected in the initial graph is the proportion of clusters they shared over clustering runs. This network is used for a final label propagation run, which identifies the consensus clusters. Users can specify a numeric vector as input to \code{consensus_cluster}, which will override the default shape parameters and number of iterations.
}

\value{
Returns a list object with the parameters and result of the clustering. If using multiple clusterings, the return value is a list of lists, with each entry corresponding to the single-clustering case. This list has two entries, \code{parameters} and \code{results}.

\code{parameters} is a named vector with the values of \code{add_self_loops} and \code{attenuation} used for the clustering. If more parameters are added in the future, they'll be included in this vector.

\code{results} differs depending on the value of \code{return_table}.

If \code{return_table=TRUE}, \code{results} is a \code{\link{data.frame}} object with two columns. The first column contains the name of each vertex, and the second column contains the cluster it was assigned to.

If \code{return_table=FALSE}, \code{results} is a character vector of length 1. This vector contains the path to the file where clusters were written to. The file is formatted as a \code{.tsv}, with each line containing two tab separated columns (vertex name, assigned cluster).
}

\references{
Traag, V.A., and L. Subelj. \emph{Large network community detection by fast label propagation.} Sci. Rep., 2023. \bold{13}(2701). https://doi.org/10.1038/s41598-023-29610-z

Leung, X.Y.I., et al., \emph{Towards real-time community detection in large networks}. Phys. Rev. E, 2009. \bold{79}(066107). https://doi.org/10.1103/PhysRevE.79.066107
}
\author{
Aidan Lakshman <AHL27@pitt.edu>
}
\section{Warning}{
While this algorithm can scale to very large graphs, it does have some internal limitations. First, nodes must be comprised of no more than 254 characters. You can increase the size yourself by changing the definition of \code{MAX_NODE_NAME_SIZE} in \code{src/OnDiskLP.c}. This limitation is provided to decrease memory overhead and improve runtime, but arbitrary values are possible.

Second, nodes are indexed using 54-bit unsigned integers. This means that the maximum possible number of nodes available is 2^54-1, which is about 1.8 quadrillion. You can decrease the size of \code{BITS_FOR_WEIGHT} in \code{src/OnDiskLP.c}, but note that this value determines how many bits to use to represent weights internally, so lower values will lead to less exact weights.

Third, this algorithm uses disk space to store large objects. As such, please ensure you have sufficient disk space for the graph you intend to process. While there are safeguards in the code itself, unhandleable errors can occur when the OS runs out of space. Use \code{\link{EstimateExoLabel}} to estimate the disk consumption of your graph, and see "Memory Consumption" for more details on how the total disk/memory consumption is calculated. Note that using \code{use_fast_sort=TRUE} will double the maximal disk consumption of the algorithm.
}

\section{Memory Consumption}{
  Let \eqn{v} be the number of unique nodes, \eqn{d} the average indegree of nodes, and \eqn{l} the average length of node labels. Note that the number of edges \eqn{e} is equivalent to \eqn{dv}.

  Specific calculations for memory/disk consumption are detailed below. In summary, the absolute worst case memory consumption is roughly \eqn{(24l+18)v} bytes, and the maximum disk consumption during computation is \eqn{16dv} bytes (or \eqn{32dv} bytes if \code{use_fast_sort=TRUE}). The final table returned consumes \eqn{(2+l+\log_{10}{v})v} bytes.

  ExoLabel builds a trie to keep track of vertex names. Each internal node of the trie consumes 24 bytes, and each leaf node consumes 18 bytes. The lowest possible RAM consumption of the trie (if every label is length \eqn{l} and shares the same prefix of length \eqn{l-1}) is roughly \eqn{40v} bytes, and the maximum RAM consumption (if no two node labels have any prefix in common) is \eqn{(24l + 18)v} bytes. We can generalize this to estimate the total memory consumption as roughly \eqn{(24(l-p)+18)v}, where \eqn{p} is the average length of common prefix between any two node labels.

  ExoLabel also uses a number of internal caches to speed up read/writes from files. These caches take around 200MB of RAM in total. Note that this calculation does not include the RAM required for R itself, which is on the order of 300MB on my machine. It also uses an internal queue for processing nodes, which consumes roughly \eqn{10v} bytes, and an internal index of size \eqn{8v} bytes.

  As for disk space, ExoLabel transforms the graph into a CSR-compressed network, which is split across three files: a header, a neighbors list, and a weights list. The header file contains \eqn{v+1} entries of 8 bytes, and the other two files consume a total of 12 bytes per outgoing edge. The number of edges to record is \eqn{vd}. Thus, the total disk consumption in bytes is \eqn{8(v+1) + 12vd \approx 8v+12dv}. However, the initial reading of the edges requires 16 bytes per edge, resulting in a maximum disk consumption of \eqn{16dv} (since \eqn{d > 2} for most graphs). If \code{use_fast_sort=TRUE}, this edge reading maximally consumes 32 bytes per edge (a maximum disk consumption of \eqn{32dv}).

  The final table returned is a tab-separated table containing vertex names and cluster numbers in human-readable format. Each line consumes at most \eqn{l + 2 + \log_{10}{v}} bytes. In the worst case, the number of clusters is equal to the number of vertices, which have \eqn{\log_{10}{v}} digits. The average number of digits is close to the number of digits of the largest number due to how number of digits scale with numbers. The extra two bytes are for the tab and newline characters. Thus, the total size of the file is at most \eqn{(2+l+\log_{10}{v})v} bytes. We remove all intermediate files prior to outputting clusters, so in practical cases this should be smaller than intermediate disk consumption.
}

\examples{
num_verts <- 20L
num_edges <- 20L
all_verts <- sample(letters, num_verts)
all_edges <- vapply(seq_len(num_edges),
      \(i) paste(c(sample(all_verts, 2L),
                   as.character(round(runif(1),3))),
                 collapse='\t'),
                    character(1L))
edgefile <- tempfile()
if(file.exists(edgefile)) file.remove(edgefile)
writeLines(all_edges, edgefile)
res <- ExoLabel(edgefile, return_table=TRUE)
print(res)

###########################
### Multiple Clustering ###
###########################
## Run with multiple add_self_loops values
tfs <- replicate(3, tempfile())
p2 <- ExoLabel(edgefile, tfs,
                add_self_loops=c(0,0.5,1),
                return_table = TRUE)
}

\seealso{
\code{\link{EstimateExoLabel}}
}
