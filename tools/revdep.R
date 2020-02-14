library(revdepcheck)

revdep_reset()
revdep_env_vars(force_suggests = FALSE)
revdep_check(dependencies = c("Depends", "Imports", "LinkingTo"), quiet=FALSE,
             timeout=as.difftime(15, units = "mins"))
revdep_reset()
