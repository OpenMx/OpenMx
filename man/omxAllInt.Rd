%
%   Copyright 2007-2012 The OpenMx Project
%
%   Licensed under the Apache License, Version 2.0 (the "License");
%   you may not use this file except in compliance with the License.
%   You may obtain a copy of the License at
% 
%        http://www.apache.org/licenses/LICENSE-2.0
% 
%   Unless required by applicable law or agreed to in writing, software
%   distributed under the License is distributed on an "AS IS" BASIS,
%   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%   See the License for the specific language governing permissions and
%   limitations under the License.

\name{omxAllInt}
\alias{omxAllInt}

\title{All Interval Multivariate Normal Integration}

\description{
   \code{omxAllInt} computes the probabilities of a large number of cells of a multivariate normal distribution that has been sliced by a varying number of thresholds in each dimension.  While the same functionality can be achieved by repeated calls to \code{\link{omxMnor}}, \code{omxAllInt} is more efficient for repeated operations on a single covariance matrix.
   \code{omxAllInt} returns an nx1 matrix of probabilities cycling from lowest to highest thresholds in each column with the rightmost variable in \var{covariance} changing most rapidly.
}

\usage{
omxAllInt(covariance, means, ...)
}

\arguments{
   \item{covariance}{the covariance matrix describing the multivariate normal distribution.}
   \item{means}{a row vector containing means of the variables of the underlying distribution.}
   \item{...}{a matrix or set of matrices containing one column of thresholds for each column of \code{covariance}.  Each column must contain a strictly increasing set of thresholds for the corresponding variable of the underlying distribution.  \code{NA} values in these thresholds indicate that the list of thresholds in that column has ended.}
}

\details{
    \var{covariance} and \var{means} contain the covariances and means of the multivariate distribution from which probabilities are to be calculated.  
    
    \var{covariance} must be a square covariance or correlation matrix with one row and column for each variable.  
    
    \var{means} must be a vector of length \code{nrows(covariance)} that contains the mean for each corresponding variable.
    
    All further arguments are considered threshold matrices.  
    
    Threshold matrices contain locations of the hyperplanes delineating the intervals to be calculated.  The first column of the first matrix corresponds to the thresholds for the first variable represented by the covariance matrix.  Subsequent columns of the same matrix correspond to thresholds for subsequent variables in the covariance matrix. If more variables exist in the covariance matrix than in the first threshold matrix, the first column of the second threshold matrix will be used, and so on.  That is, if \var{covariance} is a 4x4 matrix, and the three threshold matrices are specified, one with a single column and the others with two columns each, the first column of the first matrix will contain thresholds for the first variable in \var{covariance}, the two columns of the second matrix will correspond to the second and third variables of \var{covariance}, respectively, and the first column of the third threshold matrix will correspond to the fourth variable.  Any extra columns will be ignored.
    
    Each column in the threshold matrices must contain some number of strictly increasing thresholds, delineating the boundaries of a cell of integration.  That is, if the integral from -1 to 0 and 0 to 1 are required for a given variable, the corresponding threshold column should contain the values -1, 0, and 1, in that order. Thresholds may be set to Inf or -Inf if a boundary at positive or negative infinity is desired.
    
    Within a threshold column, a value of +Inf, if it exists, is assumed to be the largest threshold, and any rows after it are ignored in that column.  A value of NA, if it exists, indicates that there are no further thresholds in that column, and is otherwise ignored. A threshold column consisting of only +Inf or NA values will cause an error.
    
    For all i>1, the value in row i must be strictly larger than the value in row i-1 in the same column.
    
    The return value of \code{omxAllInt} is a matrix consisting of a single column with one row for each combination of threshold levels.
}

\seealso{
\code{\link{omxMnor}}
}

\examples{

data(myFADataRaw)

covariance <- cov(myFADataRaw[,1:5])
means <- colMeans(myFADataRaw[,1:5])
thresholdForColumn1 <- cbind(c(-Inf, 0,   1))               # Integrate from -Infinity to 0 and 0 to 1 on first variable
                                                            # Note: The first variable will never be calculated from 1 to +Infinity.
thresholdsForColumn2 <- cbind(c(-Inf, -1, 0, 1, Inf))       # These columns will be integrated from -Inf to -1, -1 to 0, etc.
thresholdsForColumns3and4 <- cbind(c(-Inf, 1.96, 2.326, Inf), c(-Inf, -1.96, 2.326, Inf))
omxAllInt(covariance, means, thresholdForColumn1, thresholdsForColumn2, thresholdsForColumns3and4, thresholdsForColumn2)
                                                            # Notice that columns 2 and 5 are assigned identical thresholds.

# An alternative specification of the same calculation follows
covariance <- cov(myFADataRaw[,1:5])
means <- colMeans(myFADataRaw[,1:5])
thresholds <- cbind(c(-Inf,     0,     1,  NA,  NA),        # Note NAs to indicate the end of the sequence of thresholds.
                    c(-Inf,    -1,     0,   1, Inf),
                    c(-Inf,  1.96,  2.32, Inf,  NA),
                    c(-Inf, -1.96,  2.32, Inf,  NA),
                    c(-Inf,    -1,     0,   1, Inf))
omxAllInt(covariance, means, thresholds)

}
