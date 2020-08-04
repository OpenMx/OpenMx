library(devtools)
devtools::install_github("jpritikin/covr")
Sys.setenv(NOT_CRAN="true")
Sys.setenv(OMP_NUM_THREADS=Sys.getenv("NCPUS"))
library(covr)
library(roxygen2)
options(covr.gcov = "gcov")
withr::with_makevars(getOption("covr.flags"), assignment = "+=", {
  utils::install.packages(
    ".", repos = NULL, type = "source",
    INSTALL_opts = c("--example",  "--install-tests", "--with-keep.source",
      "--no-multiarch")  # "--with-keep.parse.data" -- only available w/ R 3.6
  )
})
roxygenize('.', roclets=c('rd'))
options(digits=15)
c1 <- covr::package_coverage(type=c("tests","examples"), quiet=TRUE, pre_clean=FALSE)
pct <- percent_coverage(c1, by="line")
print(pct)
covr::codecov(coverage = c1)
if (pct < 63) {
  print(paste("Coverage dropped to", pct))
  q(status=1)
}
