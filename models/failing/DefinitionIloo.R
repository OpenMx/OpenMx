# http://openmx.psyc.virginia.edu/thread/2363

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
  mxFIMLObjective( means='eMean', covariance='Cov',thresholds='Thresh',threshnames='y',dimnames=c('x','y') )
)                                        

for (rep in 1:2) {
  if (rep == 2) model <- mxOption(model, "Number of Threads", 1)

  fit <- mxRun( model )
  omxCheckCloseEnough(fit$output$Minus2LogLikelihood, 6578.69, .1)

  est <- fit$output$estimate
  #print(est)
  omxCheckCloseEnough(est[1], 0.9998, .001)
  omxCheckCloseEnough(est[2], 0.7177, .001)
  omxCheckCloseEnough(est[3], -0.13521308, .001)
  omxCheckCloseEnough(est[4], 1.68312081, .001)
  omxCheckCloseEnough(est[5], 0.88751570, .001)
  omxCheckCloseEnough(est[6], -0.08009792, .001)
}
