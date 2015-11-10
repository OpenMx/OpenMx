library(OpenMx)

set.seed(1)
numTeachers <- 3
teacherData <- data.frame(skill=rnorm(numTeachers), teacherID=1:numTeachers)

numStudentsPerTeacher <- 3
studentData <- NULL
for (sx in 1:numStudentsPerTeacher) {
  studentData <- rbind(studentData, data.frame(
    teacherID=sx, skill=rnorm(numStudentsPerTeacher, mean=teacherData[sx, 'skill'])))
}

createIndicators <- function(indicatorVariance, latentSkill) {
  ind <- matrix(NA, length(latentSkill), length(indicatorVariance))
  for (ix in 1:length(latentSkill)) {
    ind[ix,] <- sapply(indicatorVariance,
                       function(sd) rnorm(1, mean=latentSkill[ix], sd=sd))
  }
  colnames(ind) <- paste0('i', 1:length(indicatorVariance))
  as.data.frame(ind)
}
numIndicators <- 3
teacherData <- cbind(teacherData, createIndicators(
  rlnorm(numIndicators) / 8, teacherData$skill))
studentData <- cbind(studentData, createIndicators(
  rlnorm(numIndicators) / 8, studentData$skill))

flattenBy <- function(df, key, cols, prefix='s') {
  kcol <- df[,key]
  flat <- NULL
  for (kx in unique(kcol)) {
    part <- df[kcol == kx, cols]
    flat1 <- data.frame(t(c(t(part))))
    colnames(flat1) <- apply(expand.grid(cols, 1:nrow(part), prefix)[3:1], 1, paste0, collapse='')
    flat <- rbind(flat, flat1)
  }
  flat
}

flatData <- flattenBy(studentData, 'teacherID', paste0('i',1:numIndicators))
wideData <- cbind(teacherData, flatData)

benchMod <- mxModel("bench", type="RAM",
                    mxData(type="raw", observed=wideData),
        manifestVars = c(paste0('i',1:numIndicators), colnames(flatData)),
        latentVars = c("teacherSkill", paste0('student', 1:numStudentsPerTeacher)),
        mxPath(c(paste0('i',1:numIndicators), colnames(flatData)), arrows=2,
               values=1, labels="err"),
        mxPath("one", paste0('i',1:numIndicators), labels=paste0('tm', 1:numIndicators)),
        mxPath("one", colnames(flatData), labels=paste0('sm', kronecker(rep(1,3), 1:3))),
        mxPath("teacherSkill", arrows=2, values=1, labels="teacherVar"),
        mxPath(paste0('student', 1:numStudentsPerTeacher), arrows=2,
               values=1, labels="studentVar"),
        mxPath('teacherSkill', paste0('i',1:numIndicators), labels=paste0('tl', 1:numIndicators),
               values=1,free=c(FALSE, rep(TRUE,numIndicators-1))),
        mxPath('teacherSkill', paste0('student', 1:numStudentsPerTeacher), labels="regr"),
        mxFitFunctionML(fellner=TRUE))

for (sx in 1:numStudentsPerTeacher) {
  benchMod <- mxModel(benchMod,
                      mxPath(paste0('student', sx), paste0('s',sx,'i',1:numIndicators),
                             values=1, free=c(FALSE, rep(TRUE,numIndicators-1)),
                             labels=paste0('sl', 1:numIndicators)))
}

if (1) {
  #cat(deparse(round(coef(benchFit), 3)))
  p1 <- structure(c(1.003, 2.059, 1.118, 1.009, 0.988, 0.02, 0.175, 0.574,  -0.391, -0.423, -0.442, 0.06, 0.074, 0.081),
            .Names = c("tl2",  "tl3", "regr", "sl2", "sl3", "err", "teacherVar", "studentVar",  "tm1", "tm2", "tm3", "sm1", "sm2", "sm3"))
  benchMod <- omxSetParameters(benchMod, names(p1), values=p1)
}

if (1) {
  benchFit <- mxRun(benchMod)  # 13.90182
} else {
  benchFit <- mxRun(mxModel(benchMod, mxComputeOnce('fitfunction', 'fit')))
}
#summary(benchFit)  # 13.90182

# --------------------------------------

tMod <- mxModel("teacher", type="RAM",
                mxData(type="raw", observed=teacherData,
                       primaryKey="teacherID"),  # the foreign key is matched against his column
                manifestVars = paste0('i', 1:numIndicators),
                latentVars = "teacherSkill",
                mxPath('teacherSkill', arrows=2, labels="teacherVar", values=1),
                mxPath(paste0('i',1:numIndicators), arrows=2,
                       values=1, labels="err"),
                mxPath("one", paste0('i',1:numIndicators),
                       labels=paste0('tm',1:numIndicators)),
                mxPath('teacherSkill', paste0('i',1:numIndicators),
                       labels=paste0('tl',1:numIndicators),
                       values=1,free=c(FALSE, rep(TRUE,numIndicators-1))))

sMod <- mxModel("student", type="RAM",
                # Need to ensure sort is disabled TODO
                mxData(type="raw", observed=studentData, sort=FALSE),
                manifestVars = paste0('i', 1:numIndicators),
                latentVars = "studentSkill",
                mxPath("studentSkill", arrows=2, labels="studentVar", values=1),
                mxPath(paste0('i',1:numIndicators), arrows=2,
                       values=1, labels="err"),
                mxPath("one", paste0('i',1:numIndicators),
                       labels=paste0('sm',1:numIndicators)),
                mxPath('studentSkill', paste0('i',1:numIndicators),
                       labels=paste0('sl',1:numIndicators),
                       values=1,free=c(FALSE, rep(TRUE,numIndicators-1))),

                # this is the between level regression
                mxMatrix(name="Z", nrow=1, ncol=1, free=TRUE, labels="regr",
                         dimnames=list("studentSkill", "teacherSkill")))

sMod$expectation$join <- list(mxJoin(foreignKey='teacherID',
                                     expectation='teacher',
                                     regression='Z'))

container <- mxModel("container", tMod, sMod,
                     mxFitFunctionMultigroup(c('student')))
container <- omxSetParameters(container, names(p1), values=p1)

omxCheckError(mxRun(container), "Join mapping matrix student.Z must have 4 rows: 'i1', 'i2', 'i3', and 'studentSkill'")

map <- mxMatrix(name="Z", nrow=ncol(sMod$A), ncol=1,
                dimnames=list(colnames(sMod$A), "teacherSkill"))
map$free['studentSkill', 'teacherSkill'] <- TRUE
map$labels['studentSkill', 'teacherSkill'] <- 'regr'
container$student$Z <- map
omxCheckError(mxRun(container), "Join mapping matrix student.Z must have 4 columns: 'i1', 'i2', 'i3', and 'teacherSkill'")

map <- mxMatrix(name="Z", nrow=nrow(sMod$A), ncol=ncol(tMod$A),
                dimnames=list(rownames(sMod$A), colnames(tMod$A)))
map$free['studentSkill', 'teacherSkill'] <- TRUE
map$labels['studentSkill', 'teacherSkill'] <- 'regr'
container$student$Z <- map

if (1) {
  container <- mxRun(container)
#  summary(container)
} else {
  cFit <- mxRun(mxModel(container, mxComputeOnce('fitfunction', 'fit')))
  print(cFit$output$fit)
}

omxCheckCloseEnough(logLik(container), logLik(benchFit), 1e-8)
