\name{omxReadGRMBin}
\alias{omxReadGRMBin}

\title{
Read a GCTA-Format Binary GRM into R.
}
\description{
This simple function is adapted from syntax in the GCTA User Manual.  It loads a binary genomic-relatedness matrix (GRM) from disk into R's workspace.
}
\usage{
omxReadGRMBin(prefix, AllN=FALSE, size=4, returnList=FALSE)
}
\arguments{
  \item{prefix}{Character string: everything in the path, relative to \R's working directory, and filenames of the GRM files preceding '.grm.*'.  See below, under "Details"}
  \item{AllN}{Logical.  If \code{FALSE} (default), then when \code{omxReadGRMBin()} calls \code{\link[base:readBin]{readBin}()}, it passes a value of 1 for argument \code{n}.  if \code{TRUE}, then it instead passes a value equal to the number of nonredundant elements in the GRM.}
  \item{size}{Passed to \code{\link[base:readBin]{readBin}()}.}
  \item{returnList}{Logical.  If \code{FALSE} (default), \code{omxReadGRMBin} returns the GRM.  If \code{TRUE}, then \code{omxReadGRMBin} returns a list as described below, under "Value". }
}
\details{
A GRM calculated in GCTA that is saved to disk in binary format comprises three files, the filenames of which have the same stem but different extensions.  The first, with extension "grm.bin", is the actual binary file containing the GRM elements.  The second, with extension "grm.N.bin", contains information about how many genetic markers were used to calculate the GRM.  The third, with extension "grm.id", is a text file containing two columns of data, respectively, the participant family and individual IDs.  \code{omxReadGRMBin()} is meant to be used with all three files together in the same directory.  Thus, argument \code{prefix} should be everything in the path (relative to \R's working directory) and filenames of those GRM files, up to the first period in their extensions.  In practice, it is simplest to set \R's working directory to whichever directory contains the files, and simply provide the filename stem for argument \code{prefix}.

\code{omxReadGRMBin()} opens three file connections, one for each file.
}
\value{
If \code{returnList=FALSE} (the default), then the GRM itself is returned as a numeric matrix, with each row and column named as the sum of the corresponding participant's family ID and individual ID.  Otherwise, a list of the following four elements is returned:
\enumerate{
	\item \code{"diag"}: Numeric vector containing the GRM's diagonal elements.
	\item \code{"off"}: Numeric vector containing the GRM's off-diagonal elements.
	\item \code{"id"}: Dataframe containing the family and individual IDs corresponding to the rows and columns of the GRM.
	\item \code{"N"}: Numeric; number of markers used to calculate the GRM.
}
}
\references{
Yang J, Lee SH, Goddard ME, Visscher, PM. GCTA: A tool for genome-wide complex trait analysis. American Journal of Human Genetics 2011;88:76-82. doi: 10.1016/j.ajhg.2010.11.011.

Yang J, Lee SH, Goddard ME, Visscher, PM. Genome-wide complex trait analysis (GCTA): Methods, data analyses, and interpretations. In Gondro, C et al. (Eds.), \emph{Genome-Wide Association Studies and Genomic Prediction}. New York: Springer;2013. p. 215-236.

GCTA website:  \url{https://cnsgenomics.com/software/gcta/#Overview}.

Code for \code{omxReadGRMBin()} was adapted from syntax by Jiang Yang in the GCTA User Manual, version 1.24, dated 28 July 2014, retrieved from http://cnsgenomics.com/software/gcta/GCTA_UserManual_v1.24.pdf .}
