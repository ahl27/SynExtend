###### -- fitting a right hyberbola -------------------------------------------

# calculate a point on a right hyperbola
OneSite <- function(X,
                    Bmax,
                    Kd) {
  Y <- (Bmax * X) / (Kd + X)
  return(Y)
}
