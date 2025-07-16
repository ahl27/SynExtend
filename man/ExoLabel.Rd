\name{ExoLabel}
\alias{ExoLabel}

\title{
ExoLabel: Out-of-Memory Fast Label Propagation
}
\description{
Detects communities in networks with Fast Label Propagation using disk space to drastically reduce memory overhead.
}
\usage{
ExoLabel(edgelistfiles,
              outfile=tempfile(tmpdir=tempfiledir),
              mode=c("undirected", "directed"),
              add_self_loops=FALSE,
              attenuation=TRUE,
              ignore_weights=FALSE,
              iterations=0L,
              return_table=FALSE,
              use_fast_sort=TRUE,
              verbose=interactive(),
              sep='\t',
              header=FALSE,
              tempfiledir=tempdir())
}

\arguments{
\item{edgelistfiles}{
Character; vector of files to be processed. Each entry should be a machine-interpretable path to an edgelist file. See Details for expected format.
}

\item{outfile}{
Character; file to write final clusters to. Can be set to a vector of filepaths to run multiple clusterings (see "Multiple Clusterings").
}

\item{mode}{
Character; specifies whether edges should be interpreted as undirected (default) or directed. If interpreted as directed, each edge \code{V1 V2} is interpreted as \eqn{V_1 \rightarrow V_2}. Can be "undirected", "directed", or an unambiguous abbreviation.
}

\item{add_self_loops}{
Logical or Numeric; determines if a self-loop cutoff should be added to the network. A self-loop cutoff of value \code{w} requires that at least one incoming edge has weight \code{w} in order to assign the node to that cluster (See "Self-Loops" for more information). If \code{TRUE}, adds self-loop cutoffs of weight 1.0 to all vertices. If set to numeric value \code{w}, adds self-loop cutoffs of weight \code{w} to all nodes. Can also be set to a vector when running multiple clusterings (see "Multiple Clusterings").
}

\item{attenuation}{
Logical or Numeric; determines if label-hop attenuation should be used. If \code{TRUE}, uses attenuation to prevent single clusters from dominating results. Can also be set to a numeric to influence the strength of attenuation (larger values produce larger clusters). See "Attenuation" for more information on this parameter. Can also be set to a vector when running multiple clusterings (see "Multiple Clustering").
}

\item{ignore_weights}{
Logical; determines if weights should be ignored. If \code{TRUE}, all edges will be treated as an edge of weight 1. Must be set to \code{TRUE} if any of \code{edgelistfiles} are two-column tables (start->end only, lacking a weights column).
}

\item{iterations}{
Integer; maximum number of times to process each node. If set to zero or \code{NULL}, automatically uses the square root of the max node degree. See "Algorithm Convergence" for more information.
}

\item{return_table}{
Logical; determines how the result of clustering is returned. If \code{FALSE} (default), returns a character vector corresponding to the path of \code{outfile}. If \code{TRUE}, parses \code{outfile} using \code{\link{read.table}} and returns the result (not recommended for very large graphs).
}

\item{use_fast_sort}{
Logical; determines how files should be sorted. If \code{FALSE}, ExoLabel will perform file sorting functions in-place. If \code{TRUE}, ExoLabel will perform its file sorting functions using a second temporary file. This is much faster than the in-place sort, but consumes twice the amount of disk space. The relative disk consumption is about the same size as the input graph for \code{use_fast_sort=FALSE}, and about double the size of the input graph for \code{use_fast_sort=TRUE} (see "Memory Consumption" and the last paragraph of "Warning" below). Set to \code{FALSE} if you're worried about disk utilization.
}

\item{verbose}{
Logical; determines if status messages (output, progress, etc.) should be displayed while running. Output messages are reduced if running in non-interactive mode.
}

\item{sep}{
Character; expected character that separates entries on a line in each file in \code{edgelistfiles}. Defaults to tab, as would be expected in a \code{.tsv} formatted file. Set to \code{','} for a \code{.csv} file. Also determines the separator used in the output table.
}

\item{header}{
Logical or Integer; determines if the first line of edgelist files should be skipped. If logical, \code{TRUE} skips the first line of each file.. If set to an integer \code{n}, skips the first \code{n} lines. Negative values are treated as 0, and decimals are coerced to integer.
}

\item{tempfiledir}{
Character; vector corresponding to the location where temporary files used during execution should be stored. These temporary files are deleted after ExoLabel finishes running.
}
}
\details{
ExoLabel identifies communities (clusters) in graph/network structures using a variant of Fast Label Propagation, as proposed by Traag and Subelj (2023).

However, very large graphs require too much RAM for processing on some machines. In a graph containing billions of nodes and edges, loading the entire structure into RAM is rarely feasible. ExoLabel uses disk space for storing representations of graphs. While this is slower than computing on RAM, it allows ExoLabel to scale to graphs of enormous size while only using a comparatively small amount of memory. See "Memory Consumption" for details on the total disk/memory consumption of ExoLabel.

ExoLabel expects a set of edgelist files, provided as a vector of filepaths. Each entry in the file is expected to be in the following format:

\code{VERTEX1<sep>VERTEX2<sep>WEIGHT<linesep>}

This line defines a single edge between vertices \code{VERTEX1} and \code{VERTEX2} with weight \code{WEIGHT}. \code{VERTEX1} and \code{VERTEX2} are strings corresponding to vertex names, \code{WEIGHT} is a numeric value that can be interpreted as a \code{double}. The separator \code{<sep>} corresponds to the argument \code{sep} (defaulting to tab for \code{.tsv} format), and \code{linesep} is the newline value \code{'\n'}.

If \code{ignore_weight=TRUE}, the file can be formatted as:

\code{VERTEX1<sep>VERTEX2<linesep>}

Note that the \code{VERTEX1<sep>VERTEX2<sep>WEIGHT} format is still accepted for \code{ignore_weight=FALSE}, but the weights will be ignored. Also note that only positive weights are recorded; negative and zero-weighted edges are ignored.
}

\section{Self-Loops}{
  Label Propagation algorithms are susceptible to a large number of small weights outcompeting small numbers of strong edges. While self-loops can be added to mitigate this problem, they fail to scale to larger networks because noise can scale quadratically, whereas self-loops are constants. The standard intepretation of self-loops adds a self-loop edge with fixed weight \eqn{w} to each node, essentially requiring any node's neighboring communities to have at least weight \eqn{w} to propagate. In a setting like orthology detection, spurious similarity scores will eventually outweigh both true similarities and the self-loop edges with increasing graph size.

  To combat this, we treat self-loop values as a "self-loop cutoff" rather than a fixed value. Self-loop cutoffs are a value \eqn{w'} such that all neighboring communities must have at least one edge of weight \eqn{w'} in order to propagate. With this usage, even if a node has many neighbors in the same community with spurious similarities, it must have at least one neighbor in that community with a strong similarity in order for the node to join that community. This approach scales better with the size of graphs compared to the traditional usage of self-loops.

  As an example, consider a node \eqn{N} not yet assigned to a community with 10 neighbors. Neighbors 1-9 are in community 1 with weight 0.1, and neighbor 10 is in community 2 with weight 0.8. Community 1 thus has total weight 0.9, and community 2 has weight 0.8. In the context of orthology detection, values below 0.2 are likely to be spurious. With a standard self-loop of 0.4, \eqn{N} would still be assigned to community 1, despite these being likely spurious. However, with a \emph{self-loop cutoff} of 0.4, \eqn{N} would be assigned to community 2 because no edge in community 1 is at least 0.4.
}

\section{Iterations}{
One of the main issues of Label Propagation algorithms is that they can fail to converge. Consider an unweighted directed graph with four nodes connected in a loop. That is, \code{A->B, B->C, C->D, D->A}. If \code{A,C} are in cluster 1 and \code{B,D} are in cluster 2, this algorithm could keep processing all the nodes in a loop and never converge. To solve this issue, we introduce an additional measure for convergence controlled by \code{iterations}. If \code{iterations=x}, then we only allow the algorithm to process each node \code{x} times. Once a given node has been seen \code{x} times, it is no longer updated. This can be manually specified, but defaults to the square root of the largest node indegree.
}

\section{Attenuation}{
ExoLabel also incorporates label-hop attenuation to reduce the chance of a single massive cluster dominating results, as inspired by Leung et al. (2009). In short, as a particular label propagates to other nodes, its subsequent contribution diminishes. The farther a particular label is from its original source, the less its contribution. The degree to which its contribution diminishes scales dynamically based on the proportion of nodes that update on each cycle. Each node's attenuated weight is calculated as \eqn{w' = w(1-(pd)^a)}, with \eqn{w} the node weight, \eqn{p} the proportion of nodes that changed label in the previous iteration, \eqn{d} the distance from the initial label, and \eqn{a} the attenuation power (as controlled by \code{attenuation}).

Passing a value of \code{FALSE} (equivalent to 0.0) disables attenuation entirely rather than returning all singleton clusters.

The default values of \code{TRUE} for \code{attenuation} (equivalent to 1.0) recovers the original implementation provided in Leung et al. (2009).
}

\section{Multiple Clusterings}{
Reading in the graph object takes a large portion of the processing time. This leads to a lot of duplicated effort when trying to cluster the same network under alternative parameter settings.

Multiple clusterings on the same network are supported by passing vectors of input to \code{outfile} and \code{add_self_loops} or \code{attenuation}. If the length of \code{outfile} is greater than 1, \code{add_self_loops} and \code{attenuation} can each be set to either a single value or a vector of the same length as \code{outfile}. For a single value, the same parameter value will be used across all clusterings. For multiple values, the corresponding value will be used in each clustering. See "Examples" for example usage.

Note that the order to process each node is randomly initialized, so multiple runs on the same parameters may produce different results if a random seed is not set.
}

\value{
Returns a list object with the parameters and result of the clustering. If using multiple clusterings, the return value is a list of lists, with each entry corresponding to the single-clustering case. This list has three entries, \code{parameters}, \code{graph_stats}, and \code{results}.

\code{parameters} is a named vector with the values of \code{add_self_loops}, \code{attenuation}, and \code{iterations} used for the clustering.

\code{graph_stats} is a named numeric vector containing the number of nodes and edges in the input graph.

\code{results} differs depending on the value of \code{return_table}.

If \code{return_table=TRUE}, \code{results} is a \code{\link{data.frame}} object with two columns. The first column contains the name of each vertex, and the second column contains the cluster it was assigned to.

If \code{return_table=FALSE}, \code{results} is a character vector of length 1. This vector contains the path to the file to which the clusters were written. The file is formatted as a \code{.tsv}, with each line containing two tab separated columns (vertex name, assigned cluster). Clusters are numbered from one to the total number of clusters.
}

\references{
Traag, V.A., and L. Subelj. \emph{Large network community detection by fast label propagation.} Sci. Rep., 2023. \bold{13}(2701). https://doi.org/10.1038/s41598-023-29610-z

Leung, X.Y.I., et al., \emph{Towards real-time community detection in large networks}. Phys. Rev. E, 2009. \bold{79}(066107). https://doi.org/10.1103/PhysRevE.79.066107
}
\author{
Aidan Lakshman <AHL27@pitt.edu>
}
\section{Warning}{
While this algorithm can scale to very large graphs, it does have some internal limitations. First, nodes must be comprised of no more than 255 characters. This limitation is provided to decrease memory overhead and improve runtime. This behavior is controlled by the definition of \code{MAX_NODE_NAME_SIZE} in \code{src/OnDiskLP.c}.

Second, nodes are indexed using 44-bit unsigned integers. This means that the maximum possible number of nodes available is \eqn{2^{40}-1}, which is about 17.5 trillion. This is because ExoLabel compresses weights and node labels into a single 64-bit integer to decrease disk consumption during sorting. Weights are rescaled with \eqn{w'=log_2(w+1)}, and the resulting value is transformed into a floating point number with a 16-bit mantissa and 4-bit exponent. This representation maintains a maximum error in precision of less than 0.05\%, but does result in absolute errors getting larger as weights increase in size. For a point of reference, the error in representation is less than 0.00004 for weights in [0,1] and less than 10.5 for weights in [65,000, 70,000]. This error should be undetectable outside of extremely niche scenarios.

Third, this algorithm uses disk space to store large objects. As such, please ensure you have sufficient disk space for the graph you intend to process. While there are safeguards in the code itself, unhandleable errors can occur when the OS runs out of space. Use \code{\link{EstimateExoLabel}} to estimate the disk consumption of your graph, and see "Memory Consumption" for more details on how the total disk/memory consumption is calculated. Note that using \code{use_fast_sort=TRUE} will double the maximal disk consumption of the algorithm.
}

\section{Memory Consumption}{
  Let \eqn{v} be the number of unique nodes, \eqn{d} the average indegree of nodes, and \eqn{l} the average length of node labels. Note that the number of edges \eqn{e} is equivalent to \eqn{dv}.

  Specific calculations for memory/disk consumption are detailed below. In summary, the absolute worst case memory consumption is roughly \eqn{(24l+46)v} bytes, and the maximum disk consumption during computation is \eqn{16dv} bytes (or \eqn{32dv} bytes if \code{use_fast_sort=TRUE}). In practice, the RAM consumption is closer to \eqn{46v} bytes. The final table consumes \eqn{(2+l+\log_{10}{v})v} bytes on disk.

  ExoLabel builds a trie to keep track of vertex names. Each internal node of the trie consumes 24 bytes, and each leaf node consumes 28 bytes. The lowest possible RAM consumption of the trie (if every label is length \eqn{l} and shares the same prefix of length \eqn{l-1}) is roughly \eqn{28v} bytes, and the maximum RAM consumption (if no two node labels have any prefix in common) is \eqn{(24l + 28)v} bytes. We can generalize this to estimate the total memory consumption as roughly \eqn{(24(l-p)+28)v}, where \eqn{p} is the average length of common prefix between any two node labels.

  ExoLabel also uses a number of internal caches to speed up read/writes from files. These caches take around 200MB of RAM in total irrespective of graph size. Note that this calculation does not include the RAM required for R itself. It also uses an internal queue for processing nodes, which consumes roughly \eqn{10v} bytes, and an internal index of size \eqn{8v} bytes.

  As for disk space, ExoLabel transforms the graph into a CSR-compressed network, which is split across two files: a neighbors list, and a weights list. CSR compressions also require an index, which is stored directly in the trie structure. The two files consume a total of 12 bytes per outgoing edge, for a total disk consumption of \eqn{12vd} bytes. However, the initial reading of the edges requires 16 bytes per edge, resulting in a maximum disk consumption of \eqn{16dv}. If \code{use_fast_sort=TRUE}, this edge reading maximally consumes 32 bytes per edge (a maximum disk consumption of \eqn{32dv}). Note that undirected edges are stored as two directed edges, which doubles the disk consumption.

  The final table returned contains vertex names and cluster numbers in human-readable format. Each line is of the format \code{VERTEX<sep>CLUSTER}, where \code{<sep>} is the argument passed to \code{sep}. Each line consumes at most \eqn{l + 2 + \log_{10}{v}} bytes. In the worst case, the number of clusters is equal to the number of vertices, which have at most \eqn{\log_{10}{v}} digits. The average number of digits is close to the number of digits of the largest number due to how the number of digits scales with numbers. The extra two bytes are for the sep and newline characters. Thus, the total size of the file is at most \eqn{(2+l+\log_{10}{v})v} bytes. We remove all intermediate files prior to outputting clusters, so in practical cases this should be smaller than intermediate disk consumption.
}

\examples{
## Build an example edgelist file
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

## Run ExoLabel
res_file <- ExoLabel(edgefile)
clustering <- read.delim(res_file$result, header=FALSE)
colnames(clustering) <- c("Vertex", "Cluster")
clustering


## Can also return the result directly if the network is small enough
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
