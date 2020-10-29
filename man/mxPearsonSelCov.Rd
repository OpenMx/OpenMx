\name{mxPearsonSelCov}
\alias{mxPearsonSelCov}
\alias{mxPearsonSelMean}

\title{Perform Pearson Aitken selection}

\description{
\lifecycle{maturing}
   These functions implement the Pearson Aitken selection formulae.
}

\usage{
mxPearsonSelCov(origCov, newCov)
mxPearsonSelMean(origCov, newCov, origMean)
}

\arguments{
\item{origCov}{covariance matrix. The covariance prior to selection.}
\item{newCov}{covariance matrix. A subset of \code{origCov} to replace.}
\item{origMean}{column vector. A mean vector to adjust.}
}

\details{
  Which dimensions to condition on can be communicated in one
  of two ways: (1) \code{newCov} is a submatrix of \code{origCov}.
  The dimnames are matched to determine which
  partition of \code{origCov} to replace with \code{newCov}.
  Or (2) \code{newCov} is the same dimension as \code{origCov}.
  The matrix entries are inspected to determine which entries have
  changed. The changed entries determine which partition of
\code{origCov} to replace with \code{newCov}.

Let the \eqn{n \times n} covariance matrix R (\code{origCov}) be partitioned into non-empty,
disjoint sets p and q.
Let \eqn{R_{ij}} denote the covariance matrix between the p
and q variables where the subscripts denote the variable subsets (e.g. \eqn{R_{pq}}).
Let column vectors \eqn{\mu_p} and \eqn{\mu_q} contain the means of p and q
variables, respectively.
We wish to compute the conditional covariances of the variables in q
for a subset of the population where \eqn{R_{pp}} and \eqn{\mu_p} are known (or partially known)---that is, we wish
to \emph{condition} the covariances and means of q on those of p.
Let \eqn{V_{pp}} (\code{newCov}) be an arbitrary covariance matrix of the same
dimension as \eqn{R_{pp}}.
If we replace \eqn{R_{pp}} by \eqn{V_{pp}} then the mean
of q (\code{origMean}) is transformed as
  \deqn{\mu_q \to \mu_q + R_{qp} R_{pp}^{-1} \mu_p}{See PDF version}
and the covariance of p and q are transformed as
  \deqn{\left[
\begin{array}{c|c}
R_{pp} & R_{pq} \\
\hline
R_{qp} & R_{qq}
\end{array}
\right] \to
\left[
\begin{array}{c|c}
V_{pp} & V_{pp}R_{pp}^{-1}R_{pq} \\
\hline
R_{qp}R_{pp}^{-1}V_{pp} & R_{qq}-R_{qp} (R_{pp}^{-1} - R_{pp}^{-1} V_{pp} R_{pp}^{-1}) R_{pq}
\end{array}
\right]}{See PDF version}

}

\references{
Aitken, A. (1935). Note on selection from a multivariate normal population. \emph{Proceedings of
the Edinburgh Mathematical Society (Series 2), 4}(2), 106-110.
doi:10.1017/S0013091500008063
}

\examples{
library(OpenMx)

m1 <- mxModel(
  'selectionTest',
  mxMatrix('Full', 10, 10, values=rWishart(1, 20, toeplitz((10:1)/10))[,,1],
           dimnames=list(paste0('c',1:10),paste0('c',1:10)), name="m1"),
  mxMatrix('Full', 2, 2, values=diag(2),
           dimnames=list(paste0('c',1:2),paste0('c',1:2)), name="m2"),
  mxMatrix('Full', 10, 1, values=runif(10),
           dimnames=list(paste0('c',1:10),c('v')), name="u1"),
  mxAlgebra(mxPearsonSelCov(m1, m2), name="c1"),
  mxAlgebra(mxPearsonSelMean(m1, m2, u1), name="u2")
)

m1 <- mxRun(m1)
}
