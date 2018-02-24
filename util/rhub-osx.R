library(rhub)

src <- paste0('build/', list.files(path="build", pattern="^OpenMx_.*tar.gz$"))

check(src,
	platform="macos-elcapitan-release",
	env_vars=c("_R_CHECK_FORCE_SUGGESTS_" = "false"))
