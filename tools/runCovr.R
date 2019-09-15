library(devtools)
devtools::install_github("jpritikin/covr")
library(covr)
library(roxygen2)
library(withr)

options(covr.gcov = "gcov-9")

withr::with_envvar(c(MAKE="make -j4"), {
  withr::with_makevars(getOption("covr.flags"), assignment = "+=", {
    utils::install.packages(
      ".", repos = NULL, type = "source",
      INSTALL_opts = c("--example",  "--install-tests", "--with-keep.source",
                       "--with-keep.parse.data", "--no-multiarch")
    )
  })
})
roxygenize('.', roclets=c('rd'))
c1 <- covr::package_coverage(type=c("tests","examples"), quiet=TRUE, pre_clean=FALSE)
pl <- covr:::per_line(c1)
s1 <- (t(sapply(pl, function(x) c(total=length(x$coverage),
                                 possible=sum(!is.na(x$coverage)),
                                 covered=sum(!is.na(x$coverage) & x$coverage > 0)))))
s1 <- as.data.frame(s1)
s1$pending = s1$possible - s1$covered
t1 <- colSums(s1[,2:4])
t1['covered'] / t1['possible']

percent_coverage(c1, by="line")
percent_coverage(c1, by="expression")

