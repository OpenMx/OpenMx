#########################################################################/**
# @RdocFunction %&%
#
# @title "Quadratic Product (Matrices)"
# 
# \description{
#   Define the quadratic product.
#   This can probably be optimized in some clever way.
# }
#
# \arguments{
#  \item{A}{The first matrix operand.}
#  \item{B}{The second matrix operand.}
# }
#
#*/######################################################################### 
"%&%" <- function(A, B) {
  return(A %*% B %*% t(A))
}