library(revdepcheck)

revdep_reset()
revdep_env_vars(force_suggests = FALSE)
revdep_check(dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"))

