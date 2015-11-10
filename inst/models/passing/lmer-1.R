libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")

library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE)

library(OpenMx)

if (is.factor(sleepstudy$Subject)) {
  sleepstudy$Subject <- as.integer(levels(sleepstudy$Subject)[unclass(sleepstudy$Subject)])
}

m1 <- mxModel(model="sleep", type="RAM", manifestVars=c("Reaction"), latentVars = "DayEffect",
        mxData(type="raw", observed=sleepstudy, sort = FALSE),
        mxPath(c("one"), "Reaction"),
        mxPath(c("one"), "DayEffect", free=FALSE, labels="data.Days"),
        mxPath("DayEffect", "Reaction"),
        mxPath(c("Reaction"), arrows=2, values=1),
        
        # this is the between level mapping
        mxMatrix(name="Z", nrow=1, ncol=2, values=1, labels=c('data.Days', NA),
                 dimnames=list(c("Reaction"), c("slope", "intercept"))),
        mxFitFunctionML(fellner=FALSE))

m1$expectation$join <- list(mxJoin(foreignKey="Subject",
                                   expectation="bySubject",
                                   regression='Z'))

m1 <- mxModel(m1, mxModel(
  model="bySubject", type="RAM",
  latentVars=c("slope", "intercept"),
  mxData(type="raw",
         # no data here, only primary key
         observed=data.frame(Subject=unique(sleepstudy$Subject)),
         sort=FALSE, primaryKey = "Subject"),
  mxPath(c("intercept", "slope"), arrows=2, values=1),
  mxPath("intercept", "slope", arrows=2, values=.25, labels="cov1")))

m1$bySubject$fitfunction <- NULL

omxCheckError(mxRun(m1), "Join mapping matrix sleep.Z must have 2 rows: 'Reaction' and 'DayEffect'")

# fix map matrix
map <- mxMatrix(name="Z", nrow=2, ncol=2, 
                dimnames=list(c("Reaction", 'DayEffect'), c("slope", "intercept")))
map$labels['Reaction','slope'] <- 'data.Days'
map$values['Reaction','intercept'] <- 1
m1 <- mxModel(m1, map)

omxCheckError(mxRun(m1), "sleep.fitfunction: fellner=TRUE is required for sleep.expectation")
m1$fitfunction$fellner <- TRUE

m1 <- mxRun(m1)

omxCheckCloseEnough(logLik(m1), logLik(fm1), 1e-8)
