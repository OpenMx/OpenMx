library(OpenMx)

mxOption(NULL, "Number of Threads", 1L)

data(twinData)

options(scipen=5)

for (bmi in paste0('bmi', 1:2)) {
    # adjust thresholds to better distribute sample by ordinal category
    adj1 <- 1.5
    adj2 <- 1
    twinData[[paste0('o',bmi)]] <- cut(twinData[[bmi]], c(0,18.5+adj1,25-adj2,30),
                                       labels=c("underweight","healthy","overweight"),
                                       ordered_result = TRUE)
}

summary(twinData$obmi1)
summary(twinData$obmi2)

manifests <- apply(expand.grid(c('wt','ht','obmi'),1:2),1,paste0, collapse="")
m1 <- mxModel("example", type="RAM",
        manifestVars=manifests,
        mxData(twinData, "raw"),
        mxPath(manifests, arrows=2, free=TRUE, connect="unique.pairs"),
        mxPath("one", paste0('ht',1:2), labels="ht"),
        mxPath("one", paste0('wt',1:2), labels="wt"),
        mxPath("one", paste0('obmi',1:2), free=FALSE),
        mxThreshold(paste0('obmi',1:2), 2, free=TRUE, values=mxNormalQuantiles(2), labels=paste0("th",1:2)))

m1$S$values[,] <- diag(6)
omxCheckWarning(omxCheckError(mxRun(m1),
              "The job for model 'example' exited abnormally with the error message: fit is not finite (In data 'example.data' row 883 continuous variables are too far from the model implied distribution. Details:
resid = t( matrix(c(    # 4x1
 85.9, 1.4999, 43.9, 1.4999), byrow=TRUE, nrow=1, ncol=4))
inverse covariance =  matrix(c(    # 4x4
 1.03415, -0.0769577, -0.0769577, -0.0769577
, -0.0769577, 1.03415, -0.0769577, -0.0769577
, -0.0769577, -0.0769577, 1.03415, -0.0769577
, -0.0769577, -0.0769577, -0.0769577, 1.03415), byrow=TRUE, nrow=4, ncol=4)
)"),
"In model 'example' Optimizer returned a non-zero status code 10. Starting values are not feasible. Consider mxTryHard()")
