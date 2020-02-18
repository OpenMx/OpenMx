library(revdepcheck)

args <- commandArgs(trailingOnly = TRUE)

print(getwd())
print(args[1])
print(dir('build'))

e <- c('LD_PRELOAD', paste0('NAMESERVER',1:2), 'FORCE_DNS_OVER_TCP',
  'PATH', 'LD_LIBRARY_PATH', 'PKG_CONFIG_PATH', 'IMX_SAFE',
  'INSTALLMAKEFLAGS', 'MAKEFLAGS', 'R_LIBS_USER')

revdep_check(args[1],
             dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"),
             env=c(revdep_env_vars(force_suggests = FALSE),
                   sapply(e, Sys.getenv)))
