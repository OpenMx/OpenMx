#
#   Copyright 2007-2018 The OpenMx Project
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


require(OpenMx)
require(mvtnorm)

#mxOption(NULL, 'Default optimizer', 'SLSQP')

if (mxOption(NULL, 'Default optimizer') != "SLSQP") stop("SKIP")

# Aug 11, 2016 by Hao Wu

# Updates based on May 12 version: 

# 1. add CIs for A and E parameters.
# 2. calculate regular lower limit of CI only when it is used.
# 3. modified starting values for LRTLB and LRTUB.

# Output of function ACECI

# 1. $CI is the output you want. 
#     errorLB.a through errorUB.e indicate whether the CI can be trusted, OK if 0. 
# 2. $alpha is the alpha level specified by user.
# 3. $ACE.model.fits$ACE gives the output when fitting a regular ACE model:
#     -2LL is -2 log likelihood; a and c are proportions; v is total variance; m is phenotype mean; 
#     errorACE is error code directly from OpenMx.
# 4. $ACE.model.fits$ACE.CI gives the regular CIs from a regular ACE model.
#     Some LBs are NA because they were not calculated. 
#     Error codes are directly from OpenMx.
# 5. $ACE.model.fits$AE, ~$CE and ~$E are outputs when fitting submodels.
# 6. $other.search$C, ~$A and ~$E are additional search employed for upper or lower limits of CI.
#     Some of them are NA because they were not calculated.
#     No search for LB if the parameter estimate is on boundary (LB=0 automatically).
#                   or if the parameter estimate is too far from boundary (LB is regular LB).
#     No search for UB if the parameter estimate is not on boundary (UB is regular UB).
#     When a search was performed, the LB.p.value/UB.p.value slot gives the p.value produced by the LB/UB solution.
#         It should be the alpha level specified by the user, or be larger.
#         If it is smaller than the alpha level, the LB/UB cannot be trusted.

isBadStatus <- function(st) {
	if (st == 0 || st == 6) {
		return(0)
	} else {
		return(6)
	}
}

ACECI<-function(ACEModelTwin0,alpha=0.05,silent=FALSE)
{
  crit90<-qchisq(p=1-2*alpha,df=1);
  sqrtcrit90<-sqrt(crit90);
  crit95<-qchisq(p=1-alpha,df=1);
  sqrtcrit95<-sqrt(crit95);

  # Common components
  
  varnames<-c("twin1","twin2")
  varnamelist<-list(varnames,varnames);

  A<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=.5,	label="a",lbound=0,name="A");
  C<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=.3,	label="c",lbound=0,name="C");
  V<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=1,	label="v",lbound=1e-2,name="V");
  rMZ<-mxAlgebra(expression=A+C,name="rMZ");
  rDZ<-mxAlgebra(expression=.5%x%A+C,name="rDZ");
  pdMZ<-mxConstraint(cbind(1-rMZ,1+rMZ)>cbind(1e-5,1e-5),name="pdMZ");
  pdDZ<-mxConstraint(cbind(1-rDZ,1+rDZ)>cbind(1e-5,1e-5),name="pdDZ");
  M<-mxMatrix(type="Full",nrow=1,ncol=2,free=TRUE,values=c(0,0),label=c("m","m"),name="M");
  ACESigmaMZ<-mxAlgebra(expression=rbind(cbind(V,V*rMZ),cbind(V*rMZ,V)),name="ACESigmaMZ",dimnames=varnamelist);
  ACESigmaDZ<-mxAlgebra(expression=rbind(cbind(V,V*rDZ),cbind(V*rDZ,V)),name="ACESigmaDZ",dimnames=varnamelist);
  ACEModelMZ<- ACEModelTwin0$ACEModelMZ
  ACEModelDZ<-ACEModelTwin0$ACEModelDZ

    ######### different ACE Models ############
  
  AEModelTwin<-omxSetParameters(ACEModelTwin0,"c",values=0,free=FALSE,name="AEModelTwin");
  CEModelTwin<-omxSetParameters(ACEModelTwin0,"a",values=0,free=FALSE,name="CEModelTwin");
   EModelTwin<-omxSetParameters(ACEModelTwin0,c("a","c"),values=0,free=FALSE,name="EModelTwin");
   
	##########	Fitting the unconstrained ACE Model 	#############

  fit0<-mxRun(ACEModelTwin0,intervals=TRUE,silent=silent);
	errACE0<- isBadStatus(fit0$output$status[[1]])
  if (errACE0<0||errACE0>1) {warning("ACE model failed to converge.");return(list(CI=c(NA,NA)));}

	LLACE0<-fit0$output$fit # mxEval(objective, fit) does not work here
	estACE0<-fit0$output$estimate;
  c<-estACE0["c"];
  a<-estACE0["a"];
  r<-a+c;
  parameters<-names(estACE0);

  errACE0.UB<-fit0$output$confidenceIntervalCodes[,"ubound"];
  estACE0.UB<-fit0$compute$steps[['CI']]$output[['detail']];
  estUB95.c<-c(as.matrix(estACE0.UB[estACE0.UB$parameter=="c",parameters])[1,]);
  estUB95.a<-c(as.matrix(estACE0.UB[estACE0.UB$parameter=="a",parameters])[1,]);
  estUB95.r<-c(as.matrix(estACE0.UB[!estACE0.UB$parameter%in%c("c","a"),c("value",parameters)])[1,]);

  ##########   Fitting AE Model #############
    
  fitAE<-mxRun(AEModelTwin,silent=silent);
  estAE<-fitAE$output$estimate;
  errAE<-isBadStatus(fitAE$output$status[[1]])
  LLAE<-mxEval(objective,fitAE);
  
  ##########   Fitting CE Model #############

  fitCE<-mxRun(CEModelTwin,silent=silent);
  estCE<-fitCE$output$estimate;
  errCE<-isBadStatus(fitCE$output$status[[1]])
  LLCE<-mxEval(objective,fitCE);

  ##########   Fitting E Model #############

  fitE<-mxRun(EModelTwin,silent=silent);
  estE<-fitE$output$estimate;
  errE<-isBadStatus(fitE$output$status[[1]])
  LLE<-mxEval(objective,fitE);

  ################
  ### CI of c ####
  ################
  
  ##### CI of c: case of boundary estimate ##########

	if (c<1e-5)
	{
    LB.c<-0; 
    errLB.c<-FALSE;
    LB95.c<-errLB95.c<-NA;
    LBlrt.c<-errLBlrt.c<-pL.c<-NA;
    
    ########## Obtain upper limit of CI from an unconstrained model #############
            
    ACEModelTwinC<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinC");
    ACEModelTwinC<-mxModel(ACEModelTwinC,C,mxCI("C",interval=1-alpha,type="upper",boundAdj = FALSE));
    ACEModelTwinC<-omxSetParameters(ACEModelTwinC,labels=parameters,values=estACE0);
    ACEModelTwinC<-omxSetParameters(ACEModelTwinC,labels="c",values=1e-3,lbound=NA);
    fit<-mxRun(ACEModelTwinC);
    LLACE<-fit$output$fit;
#    print(c(LLACE, coef(fit)))
    if (LLAE-LLACE>crit95) {estLB<-estUB95.c;estLB[names(estAE)]<-estAE;estLB["c"]<-0;} else
    {
      fit<-mxRun(fit,intervals=TRUE);
      estLB<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
    }
    ########## Define and run Model UB ######
    
    L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLAE,name="L")
		L2<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE,name="L2")
		bd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=sqrtcrit95,name="bd")
		alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
    twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
		d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
		d2<-mxAlgebra(expression=sqrt(max(twin-L2,0)),name="d2");
		pUB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,d2,Inf),name="pUB");
		pvalueUB<-mxConstraint(pUB>alphalevel,name="pvalueUB");
		rangeUB<-mxConstraint(cbind(d,-d2)<cbind(bd,-bd),name="rangeUB");
		negC<-mxAlgebra(expression=-C,name="negC");
		LRTUB<-mxModel("LRTUB",L,L2,bd,alphalevel,C,negC,ACEModelMZ,ACEModelDZ,twin,d,d2,pUB,rangeUB,pvalueUB,mxFitFunctionAlgebra("negC"));
		LRTUB<-omxSetParameters(LRTUB,labels="c",lbound=estLB["c"],ubound=estUB95.c["c"]);
    LRTUB<-mxOption(LRTUB,"Line search tolerance",0.1)
		startUB<-(estLB+estUB95.c)/2;
		LRTUB<-omxSetParameters(LRTUB,labels=parameters,values=startUB);

		fitUB<-mxRun(LRTUB,silent=silent);
	  errUBlrt.c<-isBadStatus(fitUB$output$status[[1]])
		UB.c<-UBlrt.c<-mxEval(C,fitUB);
		pU.c<-mxEval(pUB,fitUB);		
    errUB.c<-!((errUBlrt.c==0||errUBlrt.c==1)&&pU.c>0.0499);
	}else
    
  ##### CI of c: case of non-boundary estimate ###################
	{
		UB.c<-estUB95.c["c"];
    errUB.c<-!(errACE0.UB["c"]==0||errACE0.UB["c"]==1);
    UBlrt.c<-errUBlrt.c<-pU.c<-NA;
      
		d0<-sqrt(max(LLAE-LLACE0,0));
		
		###### If regular 90% CI's lower limit is negative, LB=0 ############
		
    if (d0<=sqrtcrit90)  {LB.c<-0;errLB.c<-!(errAE==0||errAE==1);LBlrt.c<-errLBlrt.c<-pL.c<-NA;LB95.c<-errLB95.c<-NA;} else 
    {
      ########## Obtain lower limit of CI #############
    
      ACEModelTwinC<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinC");
      ACEModelTwinC<-mxModel(ACEModelTwinC,C,mxCI("C",interval=1-alpha,type="lower",boundAdj=FALSE));
      ACEModelTwinC<-omxSetParameters(ACEModelTwinC,labels=parameters,values=estACE0+1e-3);
      fit<-mxRun(ACEModelTwinC,silent=silent,intervals=TRUE);
      errLB95.c<-fit$output$confidenceIntervalCodes[,"lbound"];
      estLB95<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
      LB95.c<-estLB95["c"];
    
      ###### If the 95% CI's lower limit is greater than the middle point, LB=LB95 ###############
      
      if (sqrtcrit95<=d0/2) 
      {
        LB.c<-LB95.c;
        errLB.c<-!((errAE==0||errAE==1)&&(errLB95.c==0||errLB95.c==1));
        LBlrt.c<-errLBlrt.c<-pL.c<-NA;
      }else
	    {    
	      ###### Find an appropriate starting value for LRTLB ############

  	    ACEModelTwinC <- mxModel(ACEModelTwinC,mxCI("C",interval=1-2*alpha,type="lower",boundAdj=FALSE));
	      if (sqrtcrit90<d0/2)  ACEModelTwinC$intervals$C$lowerdelta<-as.double(d0^2/4)
    
        fit<-mxRun(ACEModelTwinC,silent=silent,intervals=TRUE);
  	    estLB90<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
    
    	  ### Define LRTLB ################
  
      	lbd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=max(d0/2,sqrtcrit90)-1e-3,name="lbd");
      	ubd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=min(d0,sqrtcrit95),name="ubd");
	      D0<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=d0,name="D0");  
        L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE0,name="L")
       alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
        twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
      	d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
       	pLB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,(D0-d)/2+d^2/2/max(D0-d,.001*d*d),Inf),name="pLB");
       	pvalueLB<-mxConstraint(pLB>alphalevel,name="pvalueLB");
      	rangeLB1<-mxConstraint(d<ubd,name="rangeLB1");
       	rangeLB2<-mxConstraint(d>lbd,name="rangeLB2");
       	LRTLB<-mxModel("LRTLB",C,lbd,ubd,D0,L,alphalevel,ACEModelMZ,ACEModelDZ,twin,d,pLB,rangeLB1,rangeLB2,pvalueLB,mxFitFunctionAlgebra("C"));
      	LRTLB<-omxSetParameters(LRTLB,labels="c",lbound=max(LB95.c,0),ubound=estLB90["c"]);
        LRTLB<-mxOption(LRTLB,"Line search tolerance",0.1)
	
       ####################### Find LB ############
    
  	    startLB<-(estLB90+estLB95)/2;
        LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=startLB);
    	  fitLB<-mxRun(LRTLB,silent=silent);
	errLBlrt.c<-isBadStatus(fitLB$output$status[[1]])
        pL.c<-mxEval(pLB,fitLB)
        if (errLBlrt.c>1||pL.c<0.0499) 
        {
          LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=estLB90);
          fitLB<-mxRun(LRTLB,silent=silent);
          errLBlrt.c<-isBadStatus(fitLB$output$status[[1]])
          pL.c<-mxEval(pLB,fitLB);
        }
        LB.c<-LBlrt.c<-mxEval(C,fitLB);
        errLB.c<-!((errAE==0||errAE==1)&&(errLBlrt.c==0||errLBlrt.c==1)&&pL.c>=0.0499)
      }
  	}
	}

  ################
  ### CI of a ####
  ################

  ##### CI of a: case of boundary estimate ##########

  if (a<1e-5)
  {
    LB.a<-0; 
    errLB.a<-FALSE;
    LB95.a<-errLB95.a<-NA;
    LBlrt.a<-errLBlrt.a<-pL.a<-NA;
  
    ########## Obtain upper limit of CI from an unconstrained model #############
  
    ACEModelTwinA<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinA");
    ACEModelTwinA<-mxModel(ACEModelTwinA,A,mxCI("A",interval=1-alpha,type="upper", boundAdj=FALSE));
    ACEModelTwinA<-omxSetParameters(ACEModelTwinA,labels=parameters,values=estACE0);
    ACEModelTwinA<-omxSetParameters(ACEModelTwinA,labels="a",values=1e-3,lbound=NA);
    fit<-mxRun(ACEModelTwinA);
    LLACE<-fit$output$fit;
    if (LLCE-LLACE>crit95) {estLB<-estUB95.a;estLB[names(estCE)]<-estCE;estLB["a"]<-0;} else
    {
      fit<-mxRun(fit,intervals=TRUE);
      estLB<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
    }
  
    ########## Define and run Model UB ######
  
    L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLCE,name="L")
    L2<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE,name="L2")
    bd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=sqrtcrit95,name="bd")
    alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
    twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
    d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
    d2<-mxAlgebra(expression=sqrt(max(twin-L2,0)),name="d2");
    pUB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,d2,Inf),name="pUB");
    pvalueUB<-mxConstraint(pUB>alphalevel,name="pvalueUB");
    rangeUB<-mxConstraint(cbind(d,-d2)<cbind(bd,-bd),name="rangeUB");
    negA<-mxAlgebra(expression=-A,name="negA");
    LRTUB<-mxModel("LRTUB",L,L2,bd,alphalevel,A,negA,ACEModelMZ,ACEModelDZ,twin,d,d2,pUB,rangeUB,pvalueUB,mxFitFunctionAlgebra("negA"));
    LRTUB<-omxSetParameters(LRTUB,labels="a",lbound=estLB["a"],ubound=estUB95.a["a"]);
    LRTUB<-mxOption(LRTUB,"Line search tolerance",0.1)
    startUB<-(estLB+estUB95.a)/2;
    LRTUB<-omxSetParameters(LRTUB,labels=parameters,values=startUB);
  
    fitUB<-mxRun(LRTUB,silent=silent);
    errUBlrt.a<-isBadStatus(fitUB$output$status[[1]])
    UB.a<-UBlrt.a<-mxEval(A,fitUB);
    pU.a<-mxEval(pUB,fitUB);		
    errUB.a<-!((errUBlrt.a==0||errUBlrt.a==1)&&pU.a>0.0499);
  }else
  
    ##### CI of a: case of non-boundary estimate ###################
  {
    UB.a<-estUB95.a["a"];
    errUB.a<-!(errACE0.UB["a"]==0||errACE0.UB["a"]==1);
    UBlrt.a<-errUBlrt.a<-pU.a<-NA;
  
    d0<-sqrt(max(LLCE-LLACE0,0));
    
    ###### If regular 90% CI's lower limit is negative, LB=0 ############
   
    if (d0<=sqrtcrit90) {LB.a<-0;errLB.a<-!(errCE==0||errCE==1);LBlrt.a<-errLBlrt.a<-pL.a<-NA;LB95.a<-errLB95.a<-NA;} else 
    {        
      ########## Obtain lower limit of CI #############
      
      ACEModelTwinA<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinA");
      ACEModelTwinA<-mxModel(ACEModelTwinA,A,mxCI("A",interval=1-alpha,type="lower", boundAdj=FALSE));
      ACEModelTwinA<-omxSetParameters(ACEModelTwinA,labels=parameters,values=estACE0+1e-3);
      fit<-mxRun(ACEModelTwinA,silent=silent,intervals=TRUE);
      errLB95.a<-fit$output$confidenceIntervalCodes[,"lbound"];
      estLB95<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
      LB95.a<-estLB95["a"];
      
      ###### If the 95% CI's lower limit is greater than the middle point, LB=LB95 ###############
      
      if (sqrtcrit95<=d0/2) 
      {
        LB.a<-LB95.a;
        errLB.a<-!((errCE==0||errCE==1)&&(errLB95.a==0||errLB95.a==1));
        LBlrt.a<-errLBlrt.a<-pL.a<-NA;
      }else
      {    
        ###### Find an appropriate starting value for LRTLB ############
    
        ACEModelTwinA <- mxModel(ACEModelTwinA,mxCI("A",interval=1-2*alpha,type="lower",boundAdj=FALSE));
        if (sqrtcrit90<d0/2)  ACEModelTwinA$intervals$A$lowerdelta<-as.double(d0^2/4)
    
        fit<-mxRun(ACEModelTwinA,silent=silent,intervals=TRUE);
        estLB90<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][parameters])[1,]);
    
        ### Define LRTLB ################
    
        lbd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=max(d0/2,sqrtcrit90)-1e-3,name="lbd");
        ubd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=min(d0,sqrtcrit95),name="ubd");
        D0<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=d0,name="D0");  
        L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE0,name="L")
        alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
        twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
        d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
        pLB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,(D0-d)/2+d^2/2/max(D0-d,.001*d*d),Inf),name="pLB");
        pvalueLB<-mxConstraint(pLB>alphalevel,name="pvalueLB");
        rangeLB1<-mxConstraint(d<ubd,name="rangeLB1");
        rangeLB2<-mxConstraint(d>lbd,name="rangeLB2");
        LRTLB<-mxModel("LRTLB",A,lbd,ubd,D0,L,alphalevel,ACEModelMZ,ACEModelDZ,twin,d,pLB,rangeLB1,rangeLB2,pvalueLB,mxFitFunctionAlgebra("A"));
        LRTLB<-omxSetParameters(LRTLB,labels="a",lbound=c(LB95.a,0),ubound=estLB90["a"]);
        LRTLB<-mxOption(LRTLB,"Line search tolerance",0.1)
      
        ####################### Find LB ############
    
        startLB<-(estLB90+estLB95)/2;
        LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=startLB);
        fitLB<-mxRun(LRTLB,silent=silent);
        errLBlrt.a<-isBadStatus(fitLB$output$status[[1]])
        pL.a<-mxEval(pLB,fitLB)
        if (errLBlrt.a>1||pL.a<0.0499) 
        {
         LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=estLB90);
         fitLB<-mxRun(LRTLB,silent=silent);
         errLBlrt.a<-isBadStatus(fitLB$output$status[[1]])
         pL.a<-mxEval(pLB,fitLB);
         }
        LB.a<-LBlrt.a<-mxEval(A,fitLB);
        errLB.a<-!((errCE==0||errCE==1)&&(errLBlrt.a==0||errLBlrt.a==1)&&pL.a>=0.0499)
      }
    }
  }
  
  ################  
  ### CI of e ####
  ################

  ##### CI of e: case of boundary estimate ##########
  
  if (a+c<1e-5)
  {
    LB.r<-0; 
    errLB.r<-FALSE;
    LB95.r<-errLB95.r<-NA;
    LBlrt.r<-errLBlrt.r<-pL.r<-NA;
  
    ########## Obtain upper limit of CI from an unconstrained model #############
  
    ACEModelTwinE<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinE");
    ACEModelTwinE<-mxModel(ACEModelTwinE,A,C,rMZ,mxCI("rMZ",interval=1-alpha,type="upper",boundAdj=FALSE));
    ACEModelTwinE<-omxSetParameters(ACEModelTwinE,labels=parameters,values=estACE0);
    ACEModelTwinE<-omxSetParameters(ACEModelTwinE,labels=c("a","c"),values=1e-3,lbound=NA);
    LLACE<-fit$output$fit;
    fit<-mxRun(fit,intervals=TRUE);
    estLB<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][c("value",parameters)])[1,]);
    if (any(estLB[c("a","c")]<0)) {estLB[names(estE)]<-estE; estLB[c("value","a","c")]<-0;} else
      
    ########## Define and run Model UB ######
  
    L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLE,name="L")
    L2<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE,name="L2")
    bd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=sqrtcrit95,name="bd")
    alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
    twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
    d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
    d2<-mxAlgebra(expression=sqrt(max(twin-L2,0)),name="d2");
    pUB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,d2,Inf),name="pUB");
    pvalueUB<-mxConstraint(pUB>alphalevel,name="pvalueUB");
    rangeUB<-mxConstraint(cbind(d,-d2)<cbind(bd,-bd),name="rangeUB");
    negR<-mxAlgebra(expression=-rMZ,name="negR");
    LRTUB<-mxModel("LRTUB",L,L2,bd,alphalevel,A,C,rMZ,negR,ACEModelMZ,ACEModelDZ,
                   mxMatrix(type="Full",nrow=1,ncol=2,free=FALSE,values=c(estUB95.r["value"],-estLB['value']),name="rMZbd"),
                   mxConstraint(cbind(rMZ,-rMZ)<rMZbd),
                   twin,d,d2,pUB,rangeUB,pvalueUB,mxFitFunctionAlgebra("negR"));
    LRTUB<-mxOption(LRTUB,"Line search tolerance",0.1)
    startUB<-(estLB+estUB95.r)/2;
    LRTUB<-omxSetParameters(LRTUB,labels=parameters,values=startUB[parameters]);
  
    fitUB<-mxRun(LRTUB,silent=silent);
    errUBlrt.r<-isBadStatus(fitUB$output$status[[1]])
    UB.r<-UBlrt.r<-mxEval(rMZ,fitUB);
    pU.r<-mxEval(pUB,fitUB);  	
    errUB.r<-!((errUBlrt.r==0||errUBlrt.r==1)&&pU.r>0.0499);
  }else
  {
    ##### CI of e: case of non-boundary estimate ##########
    
    UB.r<-estUB95.r["value"];
    errUB.r<-!(errACE0.UB["ACEModelTwin0.rMZ[1,1]"]==0||errACE0.UB["ACEModelTwin0.rMZ[1,1]"]==1);
    UBlrt.r<-errUBlrt.r<-pU.r<-NA;
    
    d0<-sqrt(max(LLE-LLACE0,0));
    
    ###### If regular 90% CI's lower limit is negative, LB=0 ############
    
    if (d0<=sqrtcrit90) {LB.r<-0;errLB.r<-!(errE==0||errE==1);LBlrt.r<-errLBlrt.r<-pL.r<-NA;LB95.r<-errLB95.r<-NA;} else
    {
      ########## Obtain lower limit of CI #############
      
      ACEModelTwinE<-mxModel(model=ACEModelTwin0,ACEModelTwin0$intervals,remove=TRUE,name="ACEModelTwinE");
      ACEModelTwinE<-mxModel(ACEModelTwinE,rMZ,mxCI("rMZ",interval=1-alpha,type="lower", boundAdj=FALSE));
      ACEModelTwinE<-omxSetParameters(ACEModelTwinE,labels=parameters,values=estACE0+1e-3);
      fit<-mxRun(ACEModelTwinE,silent=silent,intervals=TRUE);
      errLB95.r<-fit$output$confidenceIntervalCodes[,"lbound"];
      estLB95<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][c("value",parameters)])[1,]);
      LB95.r<-estLB95["value"];
      
      ###### If the 95% CI's lower limit is greater than the middle point, LB=LB95 ###############
      
      if (sqrtcrit95<=d0/2)
      {
        LB.r<-LB95.r;
        errLB.r<-!((errE==0||errE==1)&&(errLB95.r==0||errLB95.r==1));
        LBlrt.r<-errLBlrt.r<-pL.r<-NA;
      }else
      {    
        ###### Find an appropriate starting value for LRTLB ############
      
        ACEModelTwinE <- mxModel(ACEModelTwinE,mxCI("rMZ",interval=1-2*alpha,type="lower",boundAdj=FALSE));
        if (sqrtcrit90<d0/2)  ACEModelTwinE$intervals$rMZ$lowerdelta<-as.double(d0^2/4)
      
        fit<-mxRun(ACEModelTwinE,silent=silent,intervals=TRUE);
        estLB90<-c(as.matrix(fit$compute$steps[['CI']]$output[['detail']][c("value",parameters)])[1,]);
      
        ### Define LRTLB ################
      
        lbd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=max(d0/2,sqrtcrit90)-1e-3,name="lbd");
        ubd<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=min(d0,sqrtcrit95),name="ubd");
        D0<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=d0,name="D0");  
        L<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=LLACE0,name="L")
        alphalevel<-mxMatrix(type="Full",nrow=1,ncol=1,free=FALSE,values=alpha,name="alphalevel")
        twin<-mxAlgebra(expression=ACEModelMZ.objective+ACEModelDZ.objective,name="twin");
        d<-mxAlgebra(expression=sqrt(max(twin-L,0)),name="d");
        pLB<-mxAlgebra(expression=omxMnor(1,0,d,Inf)+omxMnor(1,0,(D0-d)/2+d^2/2/max(D0-d,.001*d*d),Inf),name="pLB");
        pvalueLB<-mxConstraint(pLB>alphalevel,name="pvalueLB");
        rangeLB1<-mxConstraint(d<ubd,name="rangeLB1");
        rangeLB2<-mxConstraint(d>lbd,name="rangeLB2");
        LRTLB<-mxModel("LRTLB",A,C,rMZ,lbd,ubd,D0,L,alphalevel,ACEModelMZ,ACEModelDZ,
                    mxMatrix(type="Full",nrow=1,ncol=2,free=FALSE,values=cbind(estLB90['value'],-max(0,LB95.r)),name="rMZbd"),
                    mxConstraint(cbind(rMZ,-rMZ)<rMZbd),
                     twin,d,pLB,rangeLB1,rangeLB2,pvalueLB,mxFitFunctionAlgebra("rMZ"));
        LRTLB<-mxOption(LRTLB,"Line search tolerance",0.1)
      
        ####################### Find LB ############
        
        startLB<-(estLB90+estLB95)/2;
        LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=startLB[parameters]);
        fitLB<-mxRun(LRTLB,silent=silent);
        errLBlrt.r<-isBadStatus(fitLB$output$status[[1]])
        pL.r<-mxEval(pLB,fitLB)
        if (errLBlrt.r>1||pL.r<0.0499) 
        {
          LRTLB<-omxSetParameters(LRTLB,labels=parameters,values=estLB90[parameters]);
          fitLB<-mxRun(LRTLB,silent=silent);
          errLBlrt.r<-isBadStatus(fitLB$output$status[[1]])
          pL.r<-mxEval(pLB,fitLB);
        }
        LB.r<-LBlrt.r<-mxEval(rMZ,fitLB);
        errLB.r<-!((errE==0||errE==1)&&(errLBlrt.r==0||errLBlrt.r==1)&&pL.r>=0.0499)
      }
    }
  }

  ## output ## 
  #stop()
  CI<-c(LB.c,UB.c,LB.a,UB.a,1-UB.r,1-LB.r,errLB.c,errUB.c,errLB.a,errUB.a,errUB.r,errLB.r)
  names(CI)<-c("C.LB","C.UB","A.LB","A.UB","E.LB","E.UB","errorLB.c","errorUB.c","errorLB.a","errorUB.a","errorLB.e","errUB.e")  
  ACE<-c(LLACE0,estACE0,errACE0);
  names(ACE)<-c("-2LL",names(estACE0),"errorACE");
  ACE.CI<-c(LB95.c,estUB95.c["c"],LB95.a,estUB95.a["a"],1-estUB95.r["value"],1-LB95.r,
            errLB95.c,errACE0.UB["c"],errLB95.a,errACE0.UB["a"],errACE0.UB["ACEModelTwin0.rMZ[1,1]"],errLB95.r);
  names(ACE.CI)<-c("C.LB","C.UB","A.LB","A.UB","E.LB","E.UB","errorLB.c","errorUB.c","errorLB.a","errorUB.a","errorLB.e","errUB.e")
  AE=c(LLAE,estAE,errAE);
  names(AE)<-c("-2LL",names(estAE),"errorAE");
  CE=c(LLCE,estCE,errCE);
  names(CE)<-c("-2LL",names(estCE),"errorCE");
  E=c(LLE,estE,errE);
  names(E)<-c("-2LL",names(estE),"errorE");
  search.C=c(LBlrt.c,errLBlrt.c,pL.c,UBlrt.c,errUBlrt.c,pU.c);
  names(search.C)<-c("C.LB","error.C.LB","C.LB.p.value","C.UB","error.C.UB","C.UB.p.value");
  search.A=c(LBlrt.a,errLBlrt.a,pL.a,UBlrt.a,errUBlrt.a,pU.a);
  names(search.A)<-c("A.LB","error.A.LB","A.LB.p.value","A.UB","error.A.UB","A.UB.p.value");
  search.E=c(1-UBlrt.r,errUBlrt.r,pU.r,1-LBlrt.r,errLBlrt.r,pL.r);
  names(search.E)<-c("E.LB","error.E.LB","E.LB.p.value","E.UB","error.E.UB","E.UB.p.value");
  return(list(CI=CI,alpha=alpha,
              ACE.model.fits=list(ACE=ACE,ACE.CI=ACE.CI,AE=AE,CE=CE,E=E),
              other.search=list(C=search.C,A=search.A,E=search.E)));
}

mkmodel <- function(SigmaMZ0,SigmaDZ0, alpha=0.05) {
  ### Generating data ###
  DataMZ<-rmvnorm(100,sigma=SigmaMZ0);
  DataDZ<-rmvnorm(100,sigma=SigmaDZ0);
  
  varnames<-c("twin1","twin2")
  varnamelist<-list(varnames,varnames);
  colnames(DataMZ)<-varnames;
  colnames(DataDZ)<-varnames;
  
  ### ACE Model	components ################
  
  A<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=.5,	label="a",lbound=0,name="A");
  C<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=.3,	label="c",lbound=0,name="C");
  V<-mxMatrix(type="Full",nrow=1,ncol=1,free=TRUE,values=1,	label="v",lbound=1e-2,name="V");
  rMZ<-mxAlgebra(expression=A+C,name="rMZ");
  rDZ<-mxAlgebra(expression=.5%x%A+C,name="rDZ");
  pdMZ<-mxConstraint(cbind(1-rMZ,1+rMZ)>cbind(1e-5,1e-5),name="pdMZ");
  pdDZ<-mxConstraint(cbind(1-rDZ,1+rDZ)>cbind(1e-5,1e-5),name="pdDZ");
  M<-mxMatrix(type="Full",nrow=1,ncol=2,free=TRUE,values=c(0,0),label=c("m","m"),name="M");
  ACESigmaMZ<-mxAlgebra(expression=rbind(cbind(V,V*rMZ),cbind(V*rMZ,V)),name="ACESigmaMZ",dimnames=varnamelist);
  ACESigmaDZ<-mxAlgebra(expression=rbind(cbind(V,V*rDZ),cbind(V*rDZ,V)),name="ACESigmaDZ",dimnames=varnamelist);
  ACEModelMZ<-mxModel("ACEModelMZ",A,C,V,rMZ,M,pdMZ,ACESigmaMZ,
                      mxExpectationNormal(covariance="ACESigmaMZ",means="M",dimnames=varnames),
                      mxFitFunctionML(),mxData(observed=DataMZ,type="raw"));
  ACEModelDZ<-mxModel("ACEModelDZ",A,C,V,rDZ,M,pdDZ,ACESigmaDZ,
                      mxExpectationNormal(covariance="ACESigmaDZ",means="M",dimnames=varnames),
                      mxFitFunctionML(),mxData(observed=DataDZ,type="raw"));	
  ACEModelTwin0<-mxModel("ACEModelTwin0",A,C,rMZ,mxCI(c("A","C","rMZ"),interval=1-alpha,type="upper",boundAdj=FALSE),
                         ACEModelMZ,ACEModelDZ,mxFitFunctionMultigroup(c("ACEModelMZ","ACEModelDZ")));
  ACEModelTwin0
}

cmpInterval <- function(model, iname, param, col) {
  m2 <- model
  m2 <- mxModel(m2, remove=TRUE, m2$intervals)
  unadj <- try(mxRun(mxModel(m2, mxCI(param,boundAdj=FALSE)), intervals=TRUE))
  d1 <- unadj$compute$steps[['CI']]$output$detail
  if (is(unadj, "try-error")) return(rep(NA,3))
  adj <- try(mxRun(mxModel(m2, mxCI(param,boundAdj=TRUE)),
               intervals=TRUE,checkpoint=TRUE))
  if (is(adj, "try-error")) return(rep(NA,3))
  d2 <- unadj$compute$steps[['CI']]$output$detail
  ref <- ACECI(model)
  c(unadj=unadj$output$confidenceIntervals[param, col],
    adj=adj$output$confidenceIntervals[param, col],
    ref=ref$CI[[iname]])
}

### a = 0.9, c = 0, e = 0.1 Test C UB ### 
set.seed(3)
SigmaMZ0<-array(c(1,.9,.9,1),dim=c(2,2));
SigmaDZ0<-array(c(1,.45,.45,1),dim=c(2,2));
ace <- mkmodel(SigmaMZ0,SigmaDZ0)
got <- cmpInterval(ace, 'C.UB', 'c', 'ubound')
print(got)
omxCheckCloseEnough(got[1], .17223, 1e-4)
omxCheckCloseEnough(got[2:3], rep(0.16244, 2), 1e-4)

### a = 0.6, c = 0.3, e = 0.1 Test C LB ###
set.seed(5)
SigmaMZ0<-array(c(1,.9,.9,1),dim=c(2,2));
SigmaDZ0<-array(c(1,.6,.6,1),dim=c(2,2));
ace <- mkmodel(SigmaMZ0,SigmaDZ0)
got <- cmpInterval(ace, 'C.LB', 'c', 'lbound')
print(got)
omxCheckCloseEnough(got[1], 0.093102, 1e-4)
omxCheckCloseEnough(got[2:3], rep(0.122418, 2), 1e-4)

### a = 0, c = 0.9, e = 0.1 Test A UB ###
set.seed(19)
SigmaMZ0<-array(c(1,.9,.9,1),dim=c(2,2));
SigmaDZ0<-array(c(1,.9,.9,1),dim=c(2,2));
ace <- mkmodel(SigmaMZ0,SigmaDZ0)
got <- cmpInterval(ace, 'A.UB', 'a', 'ubound')
print(got)
omxCheckCloseEnough(got[1], 0.040689, 1e-4)
omxCheckCloseEnough(got[2:3], rep(0.033402, 2), 1e-5)

### a = 0.2, c = 0.7, e = 0.1 Test A LB ###
set.seed(8)
SigmaMZ0<-array(c(1,.9,.9,1),dim=c(2,2));
SigmaDZ0<-array(c(1,.8,.8,1),dim=c(2,2));
ace <- mkmodel(SigmaMZ0,SigmaDZ0)
got <- cmpInterval(ace, 'A.LB', 'a', 'lbound')
print(got)
omxCheckCloseEnough(got[1], 0.0687491, 1e-4)
omxCheckCloseEnough(got[2:3], rep(0.07511647, 2), 1e-5)

if (FALSE) {
  ### a = 0, c = 0, e = 1 Test E LB ###
  SigmaMZ0<-array(c(1,0,0,1),dim=c(2,2));
  SigmaDZ0<-array(c(1,0,0,1),dim=c(2,2));
  
  ### a = 0.2, c = 0.2, e = 0.6 Test E UB ###
  SigmaMZ0<-array(c(1,0.4,0.4,1),dim=c(2,2));
  SigmaDZ0<-array(c(1,0.3,0.3,1),dim=c(2,2));
}
