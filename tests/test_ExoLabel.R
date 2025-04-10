library(SynExtend)
tf1 <- tempfile()
tf2 <- tempfile()

generate_random_graph <- function(nverts, nedges){
  require(igraph)
  alph <- AA_STANDARD
  num_required <- ceiling(log(nverts, length(alph)))
  num_required <- max(num_required, 3)
  sample_names <- mkAllStrings(alph, num_required)
  labs <- sample(sample_names, nverts)
  g <- sample_gnm(nverts, nedges)
  df <- as_data_frame(g, what="edges")
  data.frame(v1=labs[df[,1]], v2=labs[df[,2]], w=runif(nedges))
}

cat("Small graphs:\n")
for(loop in c(0, 0.25, 0.5)){
  df <- generate_random_graph(10, 25)
  write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  SynExtend:::.testExoLabel(tf1, add_self_loops=loop)
}

cat("Larger graphs:\n")
for(loop in c(0, 0.5)){
  df <- generate_random_graph(10000, 25000)
  write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
  SynExtend:::.testExoLabel(tf1, add_self_loops=loop)
}


cat("Directed Edges\n")
SynExtend:::.testExoLabel(tf1, mode="directed")

cat("No fast sort\n")
SynExtend:::.testExoLabel(tf1, use_fast_sort=FALSE)

## I'll just use the same graph here
cat("Different separator\n")
write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
SynExtend:::.testExoLabel(tf1, sep=',')

cat("Multi-file input\n")
tf2 <- tempfile()
df <- generate_random_graph(50000, 100000)
write.table(df[1:50000,], tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
write.table(df[50000:100000,], tf2, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
SynExtend:::.testExoLabel(c(tf1, tf2))

cat("Larger weights\n")
df[,3] <- df[,3] * 1000
write.table(df, tf1, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')
SynExtend:::.testExoLabel(tf1)

file.remove(tf1)
file.remove(tf2)
