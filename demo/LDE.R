# ---------------------------------------------------------------------
# Program: LDE_2ndOrder_OneIndicator_5D.R
#  Author: Steven Boker and Pascal Deboeck
#    Date: Sat Sep  5 12:23:32 EDT 2009
#
# This program runs an example single indicator 2nd order LDE model
#   with 5 time-delay columns in the embedded state space matrix
#
#
# ---------------------------------------------------------------------
# Revision History
#    -- Sat Sep  5 12:23:32 EDT 2009
#      Created untitled.
#
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# Variables 
# ---------------------------------------------------------------------
#
# ---------------------------------------------------------------------

# ----------------------------------
# Read libraries and set options.

require(OpenMx)

# ----------------------------------
# Read demo data.

data(Oscillator)

plot(Oscillator[,1], type='l')


# ----------------------------------
# Set constants.

tau <- 1     # The lag between subsequent columns in the embedded matrix
deltaT <- .3  # The amount of time elapsed between subsequent observations
embedD <- 7  # The number of columns in the time-delay embedded matrix

# ----------------------------------
# Time delay embed the demo data.

Embed <- function(x, E, tau) {  # create a time delay matrix from x with embedding dimension E and lag tau
	len <- length(x)
	out <- x[1:(len-(E*tau)+tau)]
	for(i in 2:E) { out <- cbind(out,x[(1+((i-1)*tau)):(len-(E*tau)+(i*tau))]) }
	return(out)
}
	
embeddedOscillator <- Embed(Oscillator$x, embedD, tau)
dimnames(embeddedOscillator) <- list(NULL, paste("x", 1:embedD, sep=""))

# ----------------------------------
# Create the fixed LDE loading matrix.

L1 <- rep(1,embedD)
L2 <- c(1:embedD)*tau*deltaT-mean(c(1:embedD)*tau*deltaT)
L3 <- (L2^2)/2
LDE.Original <- cbind(L1,L2,L3)

# ----------------------------------
# Create 2nd order LDE model.

manifestVars <- dimnames(embeddedOscillator)[[2]]
latentVars <- paste0('d', 0:(ncol(LDE.Original)-1), 'y')

ldeModel1 <- mxModel("LDE_Model_1",
    mxMatrix("Full",  
        values=LDE.Original, 
        free=FALSE, 
        name="L", 
        byrow=TRUE,
        dimnames=list(manifestVars, latentVars)
    ),
    mxMatrix("Full", 3, 3, 
        values=c(  0,  0, 0,
                   0,  0, 0,
                 -.2,-.2, 0), 
        free=c(FALSE,FALSE,FALSE,
               FALSE,FALSE,FALSE,
	       TRUE, TRUE,FALSE),
	lbound=-2,
	ubound=2,
        name="A",
        byrow=TRUE,
        dimnames=list(latentVars, latentVars)
    ),
    mxMatrix("Symm", 3, 3,
        values=c(.8, -.1, 0,
                  -.1,.8, 0,
                  0, 0,.8), 
        free=c( TRUE, TRUE, FALSE,
               TRUE, TRUE,FALSE,
               FALSE,FALSE, TRUE), 
        name="S", 
        byrow=TRUE,
        lbound=c(0.000001, -10000, 0.000001,
                 -10000, 0.000001, 0.000001,
                 0.000001, 0.000001, 0.000001),
        dimnames=list(latentVars, latentVars)
    ),
    mxMatrix("Diag", embedD, embedD, 
        values=.2, 
        free=TRUE, 
        name="U",
        lbound=0.000001,
        dimnames=list(manifestVars, manifestVars)
    ),
    #mxExpectationNormal(cov="R"),
    mxExpectationLISREL(LY='L', BE='A', PS='S', TE='U'),
    mxFitFunctionML(),
    mxData(cov(embeddedOscillator),
        type="cov", 
        numObs=dim(embeddedOscillator)[1]
    )
)

# ----------------------------------
# Fit the LDE model and examine the summary results.

ldeModel1Fit <- mxRun(ldeModel1)

summary(ldeModel1Fit)


# ----------------------------------
# Specify the same model, but add means
#  and use raw data

ldeModel2 <- mxModel(ldeModel1,
    mxMatrix('Full', name='M', nrow=nrow(LDE.Original), ncol=1, dimnames=list(manifestVars, 'one'), free=TRUE),
    mxMatrix('Full', name='N', nrow=ncol(LDE.Original), ncol=1, dimnames=list(latentVars, 'one')),
    mxExpectationLISREL(LY='L', BE='A', PS='S', TE='U', TY='M', AL='N'),
    mxData(embeddedOscillator, 'raw')
)


# ----------------------------------
# Fit the raw data model with means

ldeModel2Fit <- mxRun(ldeModel2)

summary(ldeModel2Fit)



# ----------------------------------
# Check that parameter estimates are good

a <- c(-0.292118167, -0.108630533,  0.142931593, -0.012848260,  0.044742836,
        0.011126886,  0.007458037,  0.008658088,  0.008505263,  0.008332743,
        0.008255869,  0.008158013,  0.007594966)

b <- c(-0.292118480, -0.108629216,  0.142931822, -0.012848367,  0.044742889,
        0.011126912,  0.007458019,  0.008658099,  0.008505253,  0.008332736,
        0.008255873,  0.008158008,  0.007594939,  0.023681632,  0.011287883,
       -0.002021944, -0.012821669, -0.026326066, -0.035579923, -0.040730626)
omxCheckCloseEnough(coef(ldeModel1Fit), a, 0.001)
omxCheckCloseEnough(coef(ldeModel2Fit), b, 0.001)

