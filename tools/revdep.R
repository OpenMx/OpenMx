library(revdepcheck)

args <- commandArgs(trailingOnly = TRUE)

print(getwd())
print(args[1])
print(dir('build'))

revdep_env_vars(force_suggests = FALSE)
revdep_check(args[1],
             dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"))
