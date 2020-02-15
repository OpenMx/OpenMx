library(revdepcheck)

args <- commandArgs(trailingOnly = TRUE)

revdep_env_vars(force_suggests = FALSE)
revdep_reset(args[1])
revdep_check(args[1],
             dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"))
