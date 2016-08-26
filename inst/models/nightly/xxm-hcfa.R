# http://xxm.times.uh.edu/learn-xxm/three-level-hierarchical-model-with-observed-and-latent-variables-at-multiple-levels/

library(OpenMx)

options(width=120)
got <- suppressWarnings(try(load("models/nightly/data/hcfa.xxm.RData")))
if (is(got, "try-error")) load("data/hcfa.xxm.RData")

hcfa.school$school <- as.integer(hcfa.school$school)
hcfa.teacher <- as.data.frame(lapply(hcfa.teacher, as.integer))
for (col in c('student', 'teacher', 'school')) {
	hcfa.student[[col]] <- as.integer(hcfa.student[[col]])
}

schoolModel <- mxModel(
    "schoolModel", type="RAM",
    manifestVars=paste0('q',1:3), latentVars=paste0('eta',1:2),
    mxData(hcfa.school, 'raw', primaryKey="school"),
    mxPath('one', paste0('q',1:3)),
    mxPath('eta2', paste0('q',1:3), free=c(FALSE,rep(TRUE,2)), values=1),
    mxPath(paste0('q',1:3), arrows=2, values=1),
    mxPath(paste0('eta',1:2), arrows=2, values=1),
    mxPath('eta2', 'eta1'))    

teacherModel <- mxModel(
    "teacherModel", type="RAM", schoolModel,
    latentVars='eta',
    mxData(hcfa.teacher, 'raw', primaryKey="teacher"),
    mxPath('eta', arrows=2, values=1),
    mxPath('schoolModel.eta1', 'eta', free=FALSE, values=1, joinKey="school"))

studentModel <- mxModel(
    "studentModel", type="RAM", teacherModel,
    latentVars='eta', manifestVars=paste0('y',1:3),
    mxData(hcfa.student, 'raw'),
    mxPath('one', paste0('y',1:3)),
    mxPath('eta', paste0('y',1:3), free=c(FALSE,rep(TRUE,2)), values=1),
    mxPath(paste0('y',1:3), arrows=2, values=1),
    mxPath('eta', arrows=2, values=1),
    mxPath('teacherModel.eta', 'eta', values=1, free=FALSE, joinKey='teacher'))

studentModel <- mxRun(studentModel)
omxCheckCloseEnough(studentModel$output$fit, 7107.828, 1e-2)
