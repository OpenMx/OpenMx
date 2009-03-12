
#TODO: replace with .onLoad once namespaces are implemented
.First.lib <- function(libname, pkgname) {
   library.dynam("OpenMx") 
   sfInit(parallel = FALSE)
}
