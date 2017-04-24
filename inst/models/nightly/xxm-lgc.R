# http://xxm.times.uh.edu/learn-xxm/latent-growth-curve-model/

library(OpenMx)

options(width=120)
got <- suppressWarnings(try(load("models/nightly/data/reisby.wide.xxm.RData")))
if (is(got, "try-error")) load("data/reisby.wide.xxm.RData")

reisby1 <- mxModel(
  "reisby", type="RAM",
  mxData(reisby.wide, "raw"),
  manifestVars=paste0('depression',0:5), latentVars=c('int', 'slope'),
  mxPath(paste0('depression',0:5), arrows=2, values=1, labels="residual"),
  mxPath('int', paste0('depression',0:5), free=FALSE, values=1),
  mxPath('slope', paste0('depression',0:5), free=FALSE, values=0:5),
  mxPath(c('int', 'slope'), arrows=2, connect="unique.pairs",
         values=c(1,0,1)),
  mxPath('one', c('int', 'slope'))
)

reisby1 <- mxRun(reisby1)
omxCheckCloseEnough(reisby1$output$fit, 2219.038, 1e-2)

# -------------------------------

got <- suppressWarnings(try(load("models/nightly/data/reisby.long.xxm.RData")))
if (is(got, "try-error")) load("data/reisby.long.xxm.RData")

perSubject <- mxModel(
  "perSubject", type="RAM",
  mxData(response[!duplicated(response$subject),'subject',drop=FALSE],
         "raw", primaryKey="subject"),
  latentVars=c('int', 'slope'),
  mxPath('one', c('int', 'slope')),
  mxPath(c('int', 'slope'), arrows=2, connect="unique.pairs",
         values=c(1,0,1)))

reisby <- mxModel(
  "reisby", type="RAM", perSubject,
  mxData(response[,c('subject','depression','week')], "raw"),
  manifestVars="depression",
  mxPath('depression', arrows=2, values=1),
  mxPath('perSubject.int', 'depression', free=FALSE,
         values=1.0, joinKey="subject"),
  mxPath('perSubject.slope', 'depression', free=FALSE,
         labels='data.week', joinKey="subject"),
  mxPath('one', 'depression', free=FALSE))

reisby <- mxRun(reisby)

omxCheckCloseEnough(reisby$output$fit, reisby1$output$fit, 1e-4)
