# http://xxm.times.uh.edu/learn-xxm/two-level-confirmatory-factor-analysis/

library(OpenMx)

options(width=120)
got <- suppressWarnings(try(load("models/nightly/data/mlcfa.xxm.RData")))
if (is(got, "try-error")) load("data/mlcfa.xxm.RData")

#write.csv(mlcfa.student, file="/tmp/mlcfa.csv", row.names=FALSE)

indicators <- colnames(mlcfa.student)[3:6]

teacherModel <- mxModel(
  "teacherModel", type="RAM",
  latentVars="psi",
  mxData(mlcfa.student[!duplicated(mlcfa.student$teacher),], "raw",
         primaryKey="teacher"),
  mxPath("psi", arrows=2, values=1, lbound=1e-2, ubound=3))

studentModel <- mxModel(
  "studentModel", type="RAM", teacherModel,
  latentVars="psi", manifestVars=indicators,
  mxData(mlcfa.student, "raw"),
  mxPath("psi", arrows=2, values=1, lbound=1e-2, ubound=5),
  mxPath(indicators, arrows=2, values=1, lbound=1e-4, ubound=3),
  mxPath("psi", indicators, free=c(FALSE, rep(TRUE, 3)),
         values=1, lbound=0, ubound=5),
  mxPath("teacherModel.psi", indicators,
         free=c(FALSE, rep(TRUE, 3)), values=1, lbound=0, ubound=5,
         joinKey="teacher"),
  mxPath('one', indicators))

studentModel$expectation$.useSufficientSets <- FALSE
studentModel <- mxRun(studentModel)
summary(studentModel)

omxCheckCloseEnough(studentModel$output$fit, 14463.33, 1e-2)

studentModel$expectation$.useSufficientSets <- TRUE
altFit <- mxRun(mxModel(studentModel, 
                        mxComputeSequence(list(
                          mxComputeOnce('fitfunction', 'fit'),
                          mxComputeReportExpectation()))))
omxCheckCloseEnough(altFit$output$fit, 14463.33, 1e-2)

if(0) {
  layout <- studentModel$expectation$debug$layout
  head(layout[layout$group==2,],n=20)
}

f1 <- studentModel
f1$expectation$.rampartCycleLimit <- 0L
f1 <- mxRun(mxModel(f1, mxComputeOnce('fitfunction', 'fit')))
omxCheckCloseEnough(f1$output$fit, 14463.33, 1e-2)

f2 <- omxSetParameters(studentModel, labels=names(coef(studentModel)),
                 values=c(0.8838, 1.0983, 0.7628,   # student loadings
                          0.7915, 1.2093, 0.8297, 0.5558,  # student indicator variances
                          1.3788, # student psi
                          1.1687, 1.2239, 0.6455,  # teacher loadings
                          46.0726, 47.0572, 46.4774, 48.0057, # indicator means
                          0.5651))  # teacher psi
omxCheckCloseEnough(max(abs(coef(f2) - coef(altFit))), 0, 1e-2)
