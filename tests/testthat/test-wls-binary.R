library(OpenMx)
library(MASS)
library(testthat)

mxOption(key='Number of Threads', value=1) 

context("binary as continuous")

set.seed(1)

AA <- matrix(c(0.000000000, 0.000000000, 0.000000000, -0.093530939, 
	           0.000000000, 0.000000000, 0.000000000, -0.063013490, 
			   0.000000000, 0.000000000, 0.000000000, 0.002058237, 
			   0.000000000, 0.000000000, 0.000000000, 0.000000000), 4, 4)

II <- diag(4)

SS <- matrix(c(0.6425661, 0.0000000, 0.0000000, 0.0000000,
               0.0000000, 0.2484829, 0.0000000, 0.0000000,
               0.0000000, 0.0000000, 0.1328543, 0.0000000,
               0.0000000, 0.0000000, 0.0000000, 0.2774592),4,4)

impCov <- solve(II - AA) %&% SS

mu <-  c(5.670912, 0.4610513, 0.177309, 1.967147)

dat <- as.data.frame(mvrnorm(1000, mu = mu, Sigma = impCov))
colnames(dat) <- c("age", "sex", "snp", "sumsc")
phenoData<- dat

phenoData$sex <- ifelse(phenoData$sex > .46, 1 , 0)  # binary but treat as continuous
covariates<- c("sex")

DepVar <- "sumsc"

snpMu     <- mxPath(from = "one", to = covariates , free = T, labels = paste0(covariates ,"Mean"))
snpBeta   <- mxPath(from = covariates, to = DepVar, labels = paste0(covariates ,"Reg"), values = 0, free = T)
snpres    <- mxPath(from = covariates, arrows=2, values=1, free = T, labels = paste(covariates, "res", sep = "_"))
resid     <- mxPath(from = DepVar, arrows=2, values=1, free = T, labels = paste(DepVar, "res", sep = "_"))
itemMean  <- mxPath(from = 'one', to = DepVar, free= T, values = 1, labels = paste0(DepVar, "Mean"))
dat       <- mxData(observed=phenoData, type="raw")
fun <- mxFitFunctionWLS(type = "WLS", allContinuousMethod= "marginals")
wlsTest4 <- mxModel("GWAS", type="RAM", manifestVars = c(covariates, DepVar), snpMu, snpBeta, snpres, resid, itemMean, dat, fun  )

expect_error(mxRun(wlsTest4),
             "correlated gradients: [3,1]", fixed=TRUE)
