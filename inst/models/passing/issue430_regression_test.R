#
#   Copyright 2007-2026 by the individuals mentioned in the source code history
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#        http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#   
#   
#   Regression test for issue #430, based on syntax contributed by Kristen Kelly.


library('MASS') # For simulating correlated variables
library('OpenMx')
# mxVersion()

# Polyserial correlation with OpenMx
#    dat: A data.frame containing continvar and binvar
#    continvar: The name of the continuous variable
#    binvar: The name of the binary/ordinal variable. If binary, must be 0/1, not logical.
#    return_details: T/F. If true, return a list contaning "model" and "fit".
polyserial_with_omx <- function(dat,binvar,continvar,fitfunction=mxFitFunctionML(),return_details=T) {
	# Set up the model
	dat[,binvar] <- mxFactor(dat[,binvar],levels=c(0,1))
	covmat <- mxMatrix('Symm',2,2,byrow=T,name='covmat',values=c(1,0.2,0.2,1),free=c(T,T,T,F),
										 labels=c('var_contin','covar','covar','fixed_latent_variance'),
										 dimnames=list(c(continvar,binvar),c(continvar,binvar)))
	meanmat <- mxMatrix('Full',nrow=1,ncol=2,values=c(0,0),free=c(T,F),dimnames=list(NULL,c(continvar,binvar)),name='meanmat')
	z_thresh <- qnorm(mean(dat[,binvar]=='0',na.rm=T)) # Estimate z-score corresponding to threshold. Just a start value, OpenMx will estimate the final threshold.
	threshmat <- mxMatrix('Full',1,1,free=T,values=z_thresh,dimnames=list(binvar,NULL),name='threshmat')
	expectation <- mxExpectationNormal(covariance='covmat',means='meanmat',dimnames=c(continvar,binvar),threshold='threshmat',threshnames=binvar)
	model <- mxModel('polyserial_ML',mxData(dat, type = "raw"),covmat,meanmat,threshmat,expectation,fitfunction)
	fit <- mxRun(model)
	
	# Show results
	cat('Covariance matrix:\n')
	print(mxEval(covmat,fit))
	rho_hat <- mxEval(covar,fit)/sqrt(mxEval(var_contin,fit)) # cor(x,y) = cov(x,y)/(se(x)*se(y)). We force variance of latent binvar liability to 1, so only need variance/SE of the continuous variable.
	cat('\nrho_hat: ',rho_hat,'\n',sep='')
	cat('\n-2logL: ',fit$output$fit,'\n',sep='')
	
	# Return details if requested
	if(return_details) {
		return(list(rho_hat=rho_hat,model=model,fit=fit))
	}
}


########## Make the test data ##########

# What to simulate
n_obs <- 50000
sim_rho <- 0.2
binvar_prevalence <- 0.1

# Simulate it
set.seed(42)
testdat <- as.data.frame(mvrnorm(n_obs,mu=c(orig_continvar=0,binvar_liab=0),Sigma=matrix(c(1,sim_rho,sim_rho,1),nrow=2),empirical=T))
testdat$orig_binvar <- 1*(testdat$binvar_liab > quantile(testdat$binvar_liab,1 - binvar_prevalence)) # The top X% for liability become cases

# Make high-missingness versions of both variables (roughly 50% missingness)
testdat$highmiss_continvar <- ifelse(sample(c(T,F),n_obs,replace=T),NA,testdat$orig_continvar)
testdat$highmiss_binvar <- ifelse(sample(c(T,F),n_obs,replace=T),NA,testdat$orig_binvar)


########## Tests ##########

# The case of no missing data:
r0 <- polyserial_with_omx(testdat,'orig_binvar','orig_continvar')
omxCheckCloseEnough(r0$rho_hat,0.2,0.02) #<--Check closeness to true value.
omxCheckCloseEnough(r0$rho_hat,0.1974,5e-5) #<--Check point estimate's reproducibility.

# The case of high missingness on the binary variable:
r1 <- polyserial_with_omx(testdat,'highmiss_binvar','orig_continvar')
omxCheckCloseEnough(r1$rho_hat,0.2,0.02) #<--Check closeness to true value.
omxCheckCloseEnough(r1$rho_hat,0.184,5e-4) #<--Check point estimate's reproducibility.

# The case of rows that are missing the binary variable getting dropped:
r2 <- polyserial_with_omx(testdat[!is.na(testdat$highmiss_binvar),],'highmiss_binvar','orig_continvar')
omxCheckCloseEnough(r2$rho_hat,0.2,0.02) #<--Check closeness to true value.
omxCheckCloseEnough(r2$rho_hat,0.1860,5e-5) #<--Check point estimate's reproducibility.

# Same result when conditioning on ordinal versus continuous variables?:
r3 <- polyserial_with_omx(
	testdat,'highmiss_binvar','orig_continvar',fitfunction=mxFitFunctionML(jointConditionOn='ordinal',rowDiagnostics=T)) 
r4 <- polyserial_with_omx(
	testdat,'highmiss_binvar','orig_continvar',fitfunction=mxFitFunctionML(jointConditionOn='continuous',rowDiagnostics=T)) 
omxCheckCloseEnough(r3$rho_hat,r4$rho_hat,5e-5)
omxCheckCloseEnough(r3$fit$output$fit,r4$fit$output$fit,0.05)
m3 <- mxRun(r3$model,useOptimizer=F)
m4 <- mxRun(r4$model,useOptimizer=F)
omxCheckCloseEnough(
	cor(m3$fitfunction$info$likelihoods, m4$fitfunction$info$likelihoods), 1, 0.0001
)
omxCheckCloseEnough(
	mean(m3$fitfunction$info$likelihoods - m4$fitfunction$info$likelihoods), 0, 1e-6
)

# Same result when conditioning on ordinal versus continuous variables, with lots of missing data on the continuous variable instead?:
r5 <- polyserial_with_omx(testdat,'orig_binvar','highmiss_continvar',fitfunction=mxFitFunctionML(jointConditionOn='ordinal'))
r6 <- polyserial_with_omx(testdat,'orig_binvar','highmiss_continvar',fitfunction=mxFitFunctionML(jointConditionOn='continuous')) 
omxCheckCloseEnough(r5$rho_hat,r6$rho_hat,1e-4)
omxCheckCloseEnough(r5$fit$output$fit,r6$fit$output$fit,0.1)
