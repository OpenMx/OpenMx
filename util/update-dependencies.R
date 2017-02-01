
dir <- "."
#dir <- "~/OpenMx"
rawPkg <- read.dcf(paste(dir,"DESCRIPTION.in", sep="/"),
                fields=c('LinkingTo', 'Depends', 'Suggests'))
pkg <- unlist(strsplit(rawPkg, split=",", fixed=TRUE))
pkg <- sub("\\(.*\\)", "", pkg)
pkg <- gsub("\\s", "", pkg)

rv <- R.Version()
fullversion <- paste(rv$major, rv$minor, sep=".")

ip <- installed.packages()
have <- ip[ip[,'Package'] %in% pkg & ip[,'Built'] == fullversion, 'Version']

ap <- available.packages()
bestVersion <- ap[ap[,'Package'] %in% pkg, 'Version']

for (p1 in pkg) {
  if (p1 == 'R' || p1 == 'methods' || p1 == 'parallel' || p1 == 'Rmpi') next
  bestVersion <- ap[ap[,'Package'] == p1, 'Version']
  haveVersion <- have[p1]
  if (!length(bestVersion)) stop(p1)
  if (!is.na(haveVersion) && haveVersion == bestVersion) {
    cat(paste("package", p1, "version", bestVersion, "is already installed"), fill=T)
    next
  }
  install.packages(p1)
}
