libraries <- rownames(installed.packages())
if (!all(c("lme4","nlme") %in% libraries)) stop("SKIP")

library(lme4)
data(Orthodont, package="nlme")
Orthodont$nsex <- as.numeric(Orthodont$Sex=="Male")
Orthodont$nsexage <- with(Orthodont, nsex*age)
fm1 <- lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                (0 + nsexage|Subject), data=Orthodont, REML=FALSE)

library(OpenMx)

bySubj <- mxModel(
    model="subj", type="RAM",
    latentVars = c('intercept', paste0(c("age", 'nsex', "nsexage"), "L")),
    mxData(data.frame(Subject=unique(Orthodont$Subject)),
	   type="raw", primaryKey="Subject"),
    mxPath(from=c('intercept', 'ageL'), to=c('intercept', 'ageL'),
	   arrows=2, "unique.pairs", values=c(1,.1,1),
	   labels=c('subjInt', 'subjIntAge', 'subjAge')),
    mxPath(from=c('nsexL', 'nsexageL'), arrows=2, values=1))

ortho <- mxModel(
    model="ortho", bySubj, type="RAM", manifestVars=c("distance"),
    latentVars = c("ageL"),
    mxData(type="raw", observed=Orthodont[,c('distance', 'age',
			   'Subject', 'nsex', "nsexage")]),
    mxPath(from=c("one"), to="distance"),
    mxPath(from=c("one"), to="ageL", free=FALSE, labels="data.age"),
    mxPath(from="ageL", to="distance"),
    mxPath(from="distance", arrows=2, values=1),
    mxPath(from="subj.intercept", to="distance", values=1, free=FALSE,
	   joinKey="Subject"),
    mxPath(from=paste0("subj.", c("ageL", "nsexL", "nsexageL")),
	   to="distance",
           labels=paste0("data.", c("age", "nsex", "nsexage")),
           free=FALSE, joinKey="Subject"))

if (1) {
  # load lme4's parameters
    p1 <- ortho
    p1$subj$S$values[c('intercept', 'ageL'),c('intercept', 'ageL')] <- 
        VarCorr(fm1)$Subject
    p1$subj$S$values[c('nsexL'),c('nsexL')] <- 
        VarCorr(fm1)$Subject.1
    p1$subj$S$values[c('nsexageL'),c('nsexageL')] <- 
        VarCorr(fm1)$Subject.2

    p1$A$values['distance','ageL'] <- fixef(fm1)['age']
    p1$M$values[,'distance'] <- fixef(fm1)['(Intercept)']
    p1$S$values['distance','distance'] <- getME(fm1, "sigma")^2

    pt1 <- mxRun(mxModel(p1, mxComputeSequence(list(
        mxComputeOnce('fitfunction', 'fit'),
        mxComputeReportExpectation()))))

    omxCheckCloseEnough(logLik(pt1), logLik(fm1), 1e-6)
}

orthoFit <- mxRun(ortho)

# OpenMx finds a better solution
omxCheckCloseEnough(orthoFit$output$fit, 436.73, 1e-2)

# ------------------------------

fm2 <- lmer(distance ~ age + (age|Subject) + (0+nsex|Subject) +
                (0 + nsexage|Subject), data=Orthodont, REML=TRUE)

ortho$fitfunction$profileOut <- c("ortho.A[1,2]", "ortho.M[1,1]")

if (1) {
  # load lme4's parameters
    p1 <- ortho
    p1$subj$S$values[c('intercept', 'ageL'),c('intercept', 'ageL')] <-
        VarCorr(fm2)$Subject
    p1$subj$S$values[c('nsexL'),c('nsexL')] <-
        VarCorr(fm2)$Subject.1
    p1$subj$S$values[c('nsexageL'),c('nsexageL')] <-
        VarCorr(fm2)$Subject.2

    p1$A$values['distance','ageL'] <- fixef(fm2)['age']
    p1$M$values[,'distance'] <- fixef(fm2)['(Intercept)']
    p1$S$values['distance','distance'] <- getME(fm2, "sigma")^2

    pt1 <- mxRun(mxModel(p1, mxComputeSequence(list(
        mxComputeOnce('fitfunction', 'fit'),
        mxComputeReportExpectation()))))

    omxCheckCloseEnough(logLik(pt1), logLik(fm2), 1e-6)
}

orthoFit <- mxRun(ortho)

omxCheckCloseEnough(orthoFit$output$fit, 440.43, .01)
