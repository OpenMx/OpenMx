# http://xxm.times.uh.edu/learn-xxm/bivariate-cross-classified-model/

library(OpenMx)

options(width=120)
got <- suppressWarnings(try(load("models/nightly/data/faces.xxm.RData")))
if (is(got, "try-error")) load("data/faces.xxm.RData")

#head(faces.response)
faces.response$rater <- as.integer(faces.response$rater)
faces.response$target <- as.integer(faces.response$target)

# way too slow with OpenMx
faces.response <- subset(faces.response, response < 500)

any(table(faces.response$rater, faces.response$target) > 1)  # Rampart can't help

raterModel <- mxModel(
  "raterModel", type="RAM",
  latentVars = c("SYM", "PA"),
  mxData(data.frame(rater=unique(faces.response$rater)), "raw",
         primaryKey="rater"),
  mxPath(c("SYM", "PA"), arrows=2, connect="unique.pairs", values=c(1,0,1)))

targetModel <- mxModel(
  "targetModel", type="RAM",
  latentVars = c("SYM", "PA"),
  mxData(data.frame(target=unique(faces.response$target)), "raw",
         primaryKey="target"),
  mxPath(c("SYM", "PA"), arrows=2, connect="unique.pairs", values=c(1,0,1)))

faces <- mxModel(
  "faces", type="RAM", raterModel, targetModel,
  manifestVars=c("SYM", "PA"),
  mxData(faces.response, "raw"),
  mxPath("one", c("SYM", "PA")),
  mxPath(c("SYM", "PA"), arrows=2, connect="unique.pairs", values=c(1,0,1)),
  mxPath(paste0('raterModel.', c("SYM", "PA")), c("SYM", "PA"),
         free=FALSE, values=1, joinKey="rater"),
  mxPath(paste0('targetModel.', c("SYM", "PA")), c("SYM", "PA"),
         free=FALSE, values=1, joinKey="target"))

#faces <- mxRun(mxModel(faces, mxComputeOnce('fitfunction', 'fit')))

faces <- mxRun(faces)

summary(faces)
