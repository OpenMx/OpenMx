\name{imxReportProgress}
\alias{imxReportProgress}
\title{Report backend progress}
\usage{
imxReportProgress(info, eraseLen)
}
\arguments{
\item{info}{the character string to print}

\item{eraseLen}{the number of characters to erase}
}
\description{
Prints a show status string to the console without emitting a
newline.
}
\examples{
library(OpenMx)

previousLen <<- 0

easyReportProcess <- function(msg) {
	imxReportProgress(msg, previousLen)
	previousLen <<- nchar(msg)
}

demo <- function() {
	easyReportProcess("abc123")
	Sys.sleep(1)
	easyReportProcess("this is much longer")
	Sys.sleep(1)
	easyReportProcess("this is short")
	Sys.sleep(1)
	easyReportProcess("almost done")
	Sys.sleep(1)
	easyReportProcess("")
	cat("DONE!", fill=TRUE)
}

\donttest{demo()}
}
