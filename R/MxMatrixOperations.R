
# Define the quadratic product.
# This can probably be optimized in some clever way.
#
"%&%" <- function(A, B) {
  return(A %*% B %*% t(A))
}