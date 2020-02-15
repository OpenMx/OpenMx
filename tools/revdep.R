library(revdepcheck)

revdep_env_vars(force_suggests = FALSE)
revdep_reset('OpenMx')
revdep_check('OpenMx',
             dependencies = c("Depends", "Imports", "LinkingTo"),
             timeout=as.difftime(15, units = "mins"))
