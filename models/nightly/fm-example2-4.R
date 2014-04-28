# flexMIRT 1.88 Example 2-4 "Graded Model Combined Calibration and Scoring Syntax with Recoding"

library(OpenMx)
library(rpf)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/NCSsim.dat"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/NCSsim.dat")

m2.spec <- list()
m2.spec[1:18] <- rpf.grm(outcomes=5)
m2.numItems <- length(m2.spec)

for (c in 1:m2.numItems) {
  m2.data[[c]] <- mxFactor(m2.data[[c]], levels=1:m2.spec[[c]]$outcomes)
}

m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="ItemParam", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, seq(1,-1,length.out=4)), free=TRUE)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
#  ip.mat$values <- m2.fmfit$G1$param

m.mat <- mxMatrix(name="mean", nrow=1, ncol=1, values=0, free=FALSE)
rownames(m.mat) <- 'f1'
cov.mat <- mxMatrix(name="cov", nrow=1, ncol=1, values=1, free=FALSE)
dimnames(cov.mat) <- list('f1','f1')

m2 <- mxModel(model="m2", m.mat, cov.mat, ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(mean="mean", cov="cov",
                                ItemSpec=m2.spec,
                                ItemParam="ItemParam"),
              mxFitFunctionML(),
              mxComputeEM('expectation', 'scores',
                          mxComputeNewtonRaphson(freeSet='ItemParam')))
#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2$fitfunction$result, 140199.13, .01)

#print(m2$matrices$ItemParam$values - fmfit)
print(m2$output$backendTime)

grp <- list(spec=m2.spec,
            param=m2$matrices$ItemParam$values,
            free=apply(ip.mat$free, 2, sum),
            mean=m2$matrices$mean$values, cov=m2$matrices$cov$values)
colnames(grp$param) <- paste("i", 1:dim(grp$param)[2], sep="")
#got <- chen.thissen.1997(grp, m2.data)
#got$std
