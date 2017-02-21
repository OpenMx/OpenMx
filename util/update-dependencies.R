if (length(options("repos")) == 1 && options("repos") == "@CRAN@") {
	# avoid user interaction
	options(repos = c(CRAN = "https://cloud.r-project.org"))
}

options(download.file.method = 'wget')  # easier to route through SOCKS

ip <- installed.packages()
ap <- available.packages()

updateDependencies <- function(desc) {
	rawPkg <- read.dcf(desc,
			   fields=c('LinkingTo', 'Depends', 'Suggests', 'Imports'))
	pkg <- unlist(strsplit(rawPkg, split=",", fixed=TRUE))
	pkg <- sub("\\(.*\\)", "", pkg)
	pkg <- gsub("\\s", "", pkg)
	pkg <- pkg[!is.na(pkg)]

	rv <- R.Version()
	fullversion <- paste(rv$major, rv$minor, sep=".")

	have <- ip[ip[,'Package'] %in% pkg & ip[,'Built'] == fullversion, 'Version']

	bestVersion <- ap[ap[,'Package'] %in% pkg, 'Version']

	for (p1 in pkg) {
		if (p1 == 'R' || p1 == 'methods' || p1 == 'parallel' || p1 == 'Rmpi' || p1 == 'OpenMx') next
		bestVersion <- ap[ap[,'Package'] == p1, 'Version']
		haveVersion <- have[p1]
		if (!length(bestVersion)) stop(p1)
		if (!is.na(haveVersion) && haveVersion == bestVersion) {
			cat(paste("package", p1, "version", bestVersion, "is already installed"), fill=T)
			next
		}
		install.packages(p1)
	}
}

args <- commandArgs(trailingOnly = TRUE)
for (arg in args) {
	updateDependencies(arg)
}
