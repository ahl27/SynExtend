###### -- unit normalize a vector ---------------------------------------------

NormVec <- function(vec) {
  vec <- vec / sqrt(sum(vec^2))
  return(vec)
}
