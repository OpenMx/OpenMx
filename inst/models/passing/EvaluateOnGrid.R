library(OpenMx)
library(mvtnorm)

set.seed(1)

sigma <- diag(4) + matrix(rnorm(16,0,.1),4,4)
raw_data <- as.data.frame(rmvnorm(500, rep(0,4), .5*(sigma+t(sigma))))
colnames(raw_data) <- c(paste0('Y',1:2), paste0('X',1:2))
raw_data$Y1 <- unclass(cut(raw_data$Y1, 3, ordered_result=TRUE))-1
raw_data$Y2 <- unclass(cut(raw_data$Y2, 3, ordered_result=TRUE))-1

NOBS = nrow(raw_data)

#scale for Y1 and Y2 are 0 or 1 or 2 
set=matrix(0,2,9);set[1,1:3]=0;set[1,4:6]=1;set[1,7:9]=2;set[2,]=rep(c(0,1,2),3)
OB <- matrix(rep(0, (9*NOBS)),nrow = 9,ncol = NOBS)

#get (Y1,Y2)
y = as.matrix(raw_data[,c(1,2)])


#OB matrix: Organize 2000 observations into a 9*NOBS matrix
#9 means 3*3 possible categories, for example (0,0), (1,0)...
#each col is for one observation, for example if observation number i is (0,0), the ith col will be (1,0,0,0,0,0,0,0,0)
for (ob_num in 1:NOBS) {
  OB[which(colSums(set==y[ob_num,])==2),ob_num]=1
}



#save OB as O1 mxMatrix 
O1         <- mxMatrix( type="Full", nrow=9, ncol=NOBS,
                       free=FALSE, values=OB, name="O1")

A         <- mxMatrix( type="Lower", nrow=4, ncol=4,
                       free=FALSE, values=c(1,0,0,0,
                                            0,1,1,0,
                                            0,0,1,0,
                                            0,0,0,1), name="A")
#two covariates
NBeta_X <- 2

#read covariates
X_value <- as.matrix(raw_data[,c(3,4)])

#initial values for coefficients of covariates
Beta_X0 <- matrix(c(1,0,0,1),nrow=NBeta_X, ncol=2)*0.5

Xcov         <- mxMatrix( type="Full", nrow=NOBS, ncol=NBeta_X,
                          free=FALSE, 
                          values = X_value, name="Xcov")

Beta_X         <- mxMatrix( type="Full", nrow= NBeta_X
                            , ncol=2,
                            free=c(TRUE, FALSE, FALSE,TRUE),
                            values = Beta_X0,
                            name="Beta_X")

#mean for underlining bivariate distribution
MiuVec <- mxMatrix(type="Full", nrow=1, ncol=2,
                   free=c(FALSE,TRUE),
                   values=c(0, 0.5),
                   lbound =c(NA,-2),
                   ubound =c(NA,2),
                   labels = c(NA, "MiuT"), name="MiuVec")

#thresholds to cut underlining bivariate distribution
Y         <- mxMatrix( type="Full", nrow=4, ncol=2,
                       free=c(FALSE,TRUE,TRUE,FALSE,
                              FALSE,TRUE,TRUE,FALSE),
                       labels = c(NA,"TH1", "TH2", NA,
                                  NA,"TH1", "TH2",NA), 
                       lbound = c(NA,-10,0.000001,NA,
                                  NA,-10,0.000001,NA),
                       ubound = c(NA,10,20,NA,
                                  NA,10,20,NA),
                       values=c(-1000,-0.25, 0.25, 1000,
                                -1000,-0.25, 0.25, 1000), name="Y")

#covariance matrix for underlining bivariate distribution
S         <- mxMatrix( type="Symm", nrow=2, ncol=2,
                       free=c(FALSE,TRUE,TRUE,TRUE),
                       values=c(1, 0.2 , 0.2, 1),
                       lbound = c(NA,NA,NA,0.000001),
                       labels = c(NA,"COV21", "COV21", "COV22"),
                       name="S")

OneNobsRowVec <- mxMatrix( type="Unit", nrow=1, ncol=NOBS,
                           free=FALSE, 
                           name="OneNobsRowVec")

#transformed thresholds to cut underlining bivariate distribution
thres <- mxAlgebra(expression =   (A%*%Y), name="thres")

#my questions are here. Thanks for the comments last time
currentAbscissa <- mxMatrix(nrow=2, ncol=1, labels=c("abscissa1","abscissa2"), free=c(TRUE,TRUE), name="currentAbscissa")

#mean
MeanStr <- mxAlgebra(expression = t(MiuVec)%*%OneNobsRowVec+ t(Xcov%*%Beta_X), name="abscissa",dimnames=list(c('abscissa1','abscissa2'), NULL))
stuff <- mxAlgebra(omxAllInt(S,cvectorize(currentAbscissa),thres), name="stuff")
R1 <- mxAlgebra(mxEvaluateOnGrid(stuff, abscissa), name="R1")

#mle function
ds_cov_mle <- mxAlgebra(expression = -sum(log(R1)*O1), name="log_Like")
                        
# User-defined objective
funAl <- mxFitFunctionAlgebra("log_Like")
                        
OneLocusModel <- mxModel("ds_cov_mle", Xcov, Beta_X, MiuVec, OneNobsRowVec, A, Y, S,  MeanStr, currentAbscissa, stuff, R1, O1, thres, ds_cov_mle, funAl)                       
OneLocusFit <- mxRun(OneLocusModel)

omxCheckCloseEnough(OneLocusFit$output$fit, 855.887, .01)
omxCheckCloseEnough(
  OneLocusFit$output$estimate,
  c(-0.058, -0.051, 0.424, -0.986, 2.045, 0.069, 1.207, 0.1, 0.1 ), .01)
omxCheckCloseEnough(OneLocusFit$output$gradient, rep(0,7), 1e-3)
omxCheckCloseEnough(log(det(OneLocusFit$output$hessian)), 38.17, .1)
se <- c(OneLocusFit$output$standardErrors)
omxCheckCloseEnough(se, c(0.065, 0.083, 0.11, 0.095, 0.124, 0.094, 0.216), .01)
