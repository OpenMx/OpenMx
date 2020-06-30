# This is deprecated in favor of the remotes package
# See tools/travis/save_cache

if (length(options("repos")) == 1 && options("repos") == "@CRAN@") {
	# avoid user interaction
	options(repos = c(CRAN = "https://cloud.r-project.org"))
}

options(download.file.method = 'wget')  # easier to route through SOCKS

ip <- installed.packages()
ap <- available.packages()

extractDependencies <- function(desc) {
	rawPkg <- read.dcf(desc,
			   fields=c('LinkingTo', 'Depends', 'Suggests', 'Imports'))
	pkg <- unlist(strsplit(rawPkg, split=",", fixed=TRUE))
	pkg <- sub("\\(.*\\)", "", pkg)
	pkg <- gsub("\\s", "", pkg)
	pkg <- pkg[!is.na(pkg)]
}

updateDependencies <- function(pkg) {
	rv <- R.Version()
	fullversion <- paste(rv$major, rv$minor, sep=".")

	have <- ip[ip[,'Package'] %in% pkg & ip[,'Built'] == fullversion, 'Version']

	for (p1 in pkg) {
		bestVersion <- ap[ap[,'Package'] == p1, 'Version']
		haveVersion <- have[p1]
		if (!length(bestVersion)) cat(paste("Is", p1, "a built-in package?"),fill=T)
		if (!is.na(haveVersion) && haveVersion == bestVersion) {
			cat(paste("package", p1, "version", bestVersion, "is already installed"), fill=T)
			next
		}
#		cat(paste("install", p1),fill=T)
		install.packages(p1)
	}
}

args <- commandArgs(trailingOnly = TRUE)
pkg <- c()
for (arg in args) {
	pkg <- union(pkg, extractDependencies(arg))
}
pkg <- setdiff(pkg, c('R', 'methods', 'parallel', 'Rmpi',
                      'OpenMx', 'stats', 'utils', 'graphics', 'grDevices'))
cat(deparse(pkg), fill=TRUE)
updateDependencies(c('nloptr', 'rjson'))  # mysteriously undeclared dependencies
updateDependencies(pkg)
