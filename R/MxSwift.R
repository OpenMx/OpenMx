mxGenSwift <- function(tc, sites, sfile) {
	cmdline <- paste("swift -tc.file ",tc," -sites.file ",sites," ",sfile)
	system(cmdline)
}
