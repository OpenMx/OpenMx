# http://openmx.ssri.psu.edu/thread/2363

set.seed(1)

library(OpenMx)
N <- 2000
u <- rbinom(N,1,.5)
x <- .5*u+rnorm(N)
y <- mxFactor( rbinom(N,1,pnorm(-2+u)) , levels=c(0,1) )
 
model <- mxModel( 'BinCont',      
  mxMatrix('Full',nrow=1,ncol=2,free=c(T,T),name='Betas'),
  mxMatrix('Full',nrow=1,ncol=1,free=F,label='data.u',name='U'),
  mxMatrix('Full',nrow=1,ncol=2,free=c(T,F),name='Means'),
  mxAlgebra( Means + Betas%x%U , name='eMean'),
  mxMatrix('Full',nrow=1,ncol=1,free=T,name='Thresh'),
  mxMatrix('Symm',nrow=2,ncol=2,free=c(T,T,F),values=c(1,0,1),name='Cov'),
  mxData( data.frame(x,y,u) , type='raw'),
  mxExpectationNormal( means='eMean', covariance='Cov',thresholds='Thresh',threshnames='y',dimnames=c('x','y') ),
  mxFitFunctionML(vector=TRUE)
)                                        

maxThreads <- 6
rowLik <- matrix(NA, N, maxThreads)
for (rep in 1:maxThreads) {
  model <- mxOption(model, "Number of Threads", rep)  # failed with 3 or more!

  fit <- mxRun( model )
  rowLik[,rep] <- fit$fitfunction$result

#  print(fit$output$fit)
  
  if(0) {
    omxCheckCloseEnough(fit$output$Minus2LogLikelihood, 6903.199, .1)
    est <- fit$output$estimate
    #print(est)
    omxCheckCloseEnough(est[1], .1, .001)
    omxCheckCloseEnough(est[2], .1, .001)
    omxCheckCloseEnough(est[3], .239, .001)
    omxCheckCloseEnough(est[4], 1.331, .001)
    omxCheckCloseEnough(est[5], 1.1219, .001)
    omxCheckCloseEnough(est[6], .1135, .001)
  }
}

likRange <- t(apply(rowLik, 1, range))
omxCheckCloseEnough(max(abs(likRange[,1] - likRange[,2])), 0, 1e-5)

#------------------------------------------------------------------------------
# Add test for ordinal data generation with definition variables

omxCheckError(mxGenerateData(fit, N-1), 'Definition variable(s) found, but the number of rows in the data do not match the number of rows requested for data generation.')


fakeData <- mxGenerateData(fit, N)

omxCheckEquals(colnames(fakeData), colnames(fit$data$observed))
omxCheckTrue(all(fakeData$u == u))
omxCheckEquals(levels(fakeData$y), levels(fit$data$observed$y))

model2 <- mxModel(model, mxData(fakeData, 'raw'))
fit2 <- mxRun(model2)

omxCheckTrue(cor(coef(fit), coef(fit2)) > .99)
rms <- function(x, y){sqrt(mean((x-y)^2))}
omxCheckTrue(rms(coef(fit), coef(fit2)) < .05)


