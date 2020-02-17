library(revdepcheck)

args <- commandArgs(trailingOnly = TRUE)

print(args[1])

revdep_env_vars(force_suggests = FALSE)
revdep_check(args[1],
             dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"))
