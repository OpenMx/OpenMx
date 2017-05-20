# flexMIRT 1.88 Example 2-4 "Graded Model Combined Calibration and Scoring Syntax with Recoding"

library(OpenMx)
library(rpf)

set.seed(1)
m2.data <- suppressWarnings(try(read.table("models/nightly/data/NCSsim.dat"), silent=TRUE))
if (is(m2.data, "try-error")) m2.data <- read.table("data/NCSsim.dat")

m2.spec <- list()
m2.spec[1:18] <- list(rpf.grm(outcomes=5))
m2.numItems <- length(m2.spec)

for (c in 1:m2.numItems) {
  m2.data[[c]] <- mxFactor(m2.data[[c]], levels=1:m2.spec[[c]]$outcomes)
}

m2.maxParam <-max(sapply(m2.spec, rpf.numParam))

ip.mat <- mxMatrix(name="item", nrow=m2.maxParam, ncol=m2.numItems,
                   values=c(1, seq(1,-1,length.out=4)), free=TRUE)
rownames(ip.mat) <- c('f1', paste('b', 1:4, sep=""))
colnames(ip.mat) <- colnames(m2.data)

#  m2.fmfit <- read.flexmirt("~/2012/sy/fm/ms-rasch-prm.txt")
# cat(deparse(round(m2.fmfit$G1$param,6)))
#  ip.mat$values <- m2.fmfit$G1$param

m2 <- mxModel(model="m2", ip.mat,
              mxData(observed=m2.data, type="raw"),
              mxExpectationBA81(ItemSpec=m2.spec),
              mxFitFunctionML(),
              mxComputeEM('expectation', 'scores',
                          mxComputeNewtonRaphson(freeSet='item')))
#  m2 <- mxOption(m2, "Number of Threads", 1)
m2 <- mxRun(m2, silent=TRUE)
omxCheckCloseEnough(m2$fitfunction$result, 140199.13, .01)

emstat <- m2$compute$output
omxCheckCloseEnough(emstat$EMcycles, 19, 2)
omxCheckCloseEnough(emstat$totalMstep, 54, 10)
omxCheckCloseEnough(m2$output$evaluations, 109, 5)

#print(m2$matrices$item$values - fmfit)
print(m2$output$backendTime)

grp <- list(spec=m2.spec,
            param=m2$matrices$item$values,
            free=apply(ip.mat$free, 2, sum))
colnames(grp$param) <- paste("i", 1:dim(grp$param)[2], sep="")
#got <- chen.thissen.1997(grp, m2.data)
#got$std
