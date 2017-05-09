libraries <- rownames(installed.packages())
if (!("lme4" %in% libraries)) stop("SKIP")

library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy, REML=FALSE)

library(OpenMx)

# ------------------- hybrid matrix/path spec
if(1) {
m1 <- mxModel(model="sleep", type="RAM", manifestVars=c("Reaction"), latentVars = "DayEffect",
        mxData(type="raw", observed=sleepstudy),
        mxPath(c("one"), "Reaction"),
        mxPath(c("one"), "DayEffect", free=FALSE, labels="data.Days"),
        mxPath("DayEffect", "Reaction"),
        mxPath(c("Reaction"), arrows=2, values=1),
        
        # this is the between level mapping
        mxMatrix(name="Z", nrow=1, ncol=2, values=1, labels=c('data.Days', NA),
                 dimnames=list(c("Reaction"), c("slope", "intercept")),
                 joinKey = "Subject", joinModel = "bySubject"),
        mxFitFunctionML(fellner=FALSE))

m1$expectation$between <- "Z"

m1 <- mxModel(m1, mxModel(
  model="bySubject", type="RAM",
  latentVars=c("slope", "intercept"),
  mxData(type="raw",
         # no data here, only primary key
         observed=data.frame(Subject=unique(sleepstudy$Subject)),
         primaryKey = "Subject"),
  mxPath(c("intercept", "slope"), arrows=2, values=1),
  mxPath("intercept", "slope", arrows=2, values=.25, labels="cov1")))

m1$bySubject.S[1,1]$ubound <- 100
m1$bySubject$fitfunction <- NULL

omxCheckError(mxRun(m1), "Join mapping matrix sleep.Z must have 2 rows: 'Reaction' and 'DayEffect'")

# fix map matrix
map <- mxMatrix(name="Z", nrow=2, ncol=2, 
                dimnames=list(c("Reaction", 'DayEffect'), c("slope", "intercept")),
                joinKey = "Subject", joinModel = "bySubject")
map$labels['Reaction','slope'] <- 'data.Days'
map$values['Reaction','intercept'] <- 1
m1 <- mxModel(m1, map)

omxCheckError(mxRun(m1), "sleep.fitfunction: fellner=TRUE is required for sleep.expectation")
m1$fitfunction$fellner <- TRUE

m1 <- mxRun(m1)

omxCheckCloseEnough(logLik(m1), logLik(fm1), 1e-6)
}
# ------------------- all path spec

m2 <- mxModel(model="sleep", type="RAM", manifestVars=c("Reaction"), latentVars = "DayEffect",
              mxData(type="raw", observed=sleepstudy),
              mxPath(c("one"), "Reaction"),
              mxPath(c("one"), "DayEffect", free=FALSE, labels="data.Days"),
              mxPath("DayEffect", "Reaction", values=1),
              mxPath(c("Reaction"), arrows=2, values=1))

omxCheckError(mxModel(m2, mxPath('by_Subject.intercept', 'Reaction',
                   values=1, free=FALSE, joinKey="Subject")),
              "Nice try. You need to create an upper level RAM model called 'by_Subject' and add it as a submodel of 'sleep' before you can create paths between these models.")

bySub <- mxModel(
  model="by_Subject", type="RAM",
  latentVars=c("slope", "intercept"),
  mxData(type="raw",
         # no data here, only primary key
         observed=data.frame(Subject=unique(sleepstudy$Subject)),
         primaryKey = "Subject"),
  mxPath(c("intercept", "slope"), arrows=2, values=1),
  mxPath("intercept", "slope", arrows=2, values=.25, labels="cov1"),
  mxPath('one', c("intercept", "slope"), free=FALSE))

m2 <- mxModel(m2, bySub)

omxCheckError(mxModel(m2, mxPath('intercept', 'Reaction',
                         values=1, free=FALSE, joinKey="Subject")),
              "Between level paths must be from modelName.variableName, not 'intercept'")

omxCheckError(mxModel(m2, mxPath(c('Bathtub.intercept', 'Sink.intercept'), 'Reaction',
                                 values=1, free=FALSE, joinKey="Subject")),
              "Nice try. You need to create an upper level RAM model called 'Bathtub' and add it as a submodel of 'sleep' before you can create paths between these models.")

omxCheckError(mxModel(m2, mxPath('by_Subject.intersept', 'Reaction',
                                 values=1, free=FALSE, joinKey="Subject")),
              "Nice try, you need to add 'intersept' to either manifestVars or latentVars in model 'by_Subject' before you can use them in a path.")

omxCheckError(mxModel(m2, mxPath('by_Subject.intercept', 'React',
                                 values=1, free=FALSE, joinKey="Subject")),
              "Nice try, you need to add 'React' to either manifestVars or latentVars before you can use them in a path.")

m2 <- mxModel(m2,
              mxPath('by_Subject.intercept', 'Reaction',
                     values=1, free=FALSE, joinKey="Subject"),
              mxPath('by_Subject.slope', 'Reaction', labels='data.Days',
                     free=FALSE, joinKey="Subject"))

m2 <- mxRun(m2)
omxCheckCloseEnough(logLik(m2), logLik(fm1), 1e-6)

oldBetween <- m2$expectation$between
m2$expectation$between <- c(m2$expectation$between, "whatever")
omxCheckError(mxRun(m2), "Level transition matrix 'whatever' listed in 'sleep.expectation' is not found")
m2$expectation$between <- oldBetween
