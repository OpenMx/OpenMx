########################################################################
########################################################################
########################################################################

######################################
###
### Optimization problems in OpenMx
###
######################################

### Libraries used
library( OpenMx )
library( polycor ) # For independent calculation of tetrachoric correlation

### Generate data
datS <- mxFactor( 
  data.frame(
    x1=c(rep(0,1600000),rep(1,800),rep(0,810),rep(1,195),rep(0,810),rep(1,3),rep(0,9),rep(1,2),rep(0,780),rep(1,6),rep(0,6),rep(1,2),rep(0,200),rep(1,2),rep(1,1)) ,
    x2=c(rep(0,1600000),rep(0,800),rep(1,810),rep(1,195),rep(0,810),rep(0,3),rep(1,9),rep(1,2),rep(0,780),rep(0,6),rep(1,6),rep(1,2),rep(0,200),rep(0,2),rep(1,1)) ,
    x3=c(rep(0,1600000),rep(0,800),rep(0,810),rep(0,195),rep(1,810),rep(1,3),rep(1,9),rep(1,2),rep(0,780),rep(0,6),rep(0,6),rep(0,2),rep(1,200),rep(1,2),rep(1,1)) ,
    x4=c(rep(0,1600000),rep(0,800),rep(0,810),rep(0,195),rep(0,810),rep(0,3),rep(0,9),rep(0,2),rep(1,780),rep(1,6),rep(1,6),rep(1,2),rep(1,200),rep(1,2),rep(1,1))
  ) , levels=c(0,1)
)
varNames <- names(datS)

### Observed tetrachoric correlations
obsMat <- matrix(1,4,4)
for(i in 1:3){
  for(j in (i+1):4){
    print( table(datS[,c(i,j)] ) ) # Print each table
    obsMat[i,j] <- obsMat[j,i] <- polychor(table(datS[,c(i,j)]) , ML=T )
  }
}
round( obsMat , 3 )

### Observed prevalences
colMeans( 1*(datS==1) )

### Problem inverting observed correlation matrix?
round( solve( obsMat ) , 4 )

### Very simple model
simpMod <- mxModel('SimpMod',
### Expected covariance matrix (variance==1 due to binary variables)
  mxMatrix(type='Stand',nrow=4,ncol=4,free=T,values=.3,name='eC'),
### Expected Means 
  mxMatrix(type='Full',nrow=1,ncol=4,free=F,values=0,name='eM'),
### Expected thresholds
  mxMatrix(type='Full',nrow=1,ncol=4,free=T,values=3,name='eT'),
### Data and likelihood
  mxModel('Datmod',
    mxData(datS , type='raw'),
    mxExpectationNormal(means='SimpMod.eM',covariance='SimpMod.eC',thresholds='SimpMod.eT',threshnames=varNames,dimnames=varNames),
    mxFitFunctionML()
  ),
### Add the likelihoods (only one)
  mxFitFunctionMultigroup( c('Datmod.fitfunction') )
)
### Fit model
mxOption( NULL , 'Number of Threads' , 6 )
simpModFit <- mxRun( simpMod , intervals=F )
#simpModFit <- mxTryHard( simpMod , intervals=F )
summary(simpModFit)

# Modelled correlations
round( simpModFit$eC$values , 4 )
# Compare with independently calculated tetrachoric correlations
round( obsMat , 3 )

########################################################################
########################################################################
########################################################################
