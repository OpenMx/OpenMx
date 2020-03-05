#------------------------------------------------------------------------------
# Author: Michael D. Hunter
# Date: 2020-03-05 16:53:37
# Filename: RampartDimnames.R
# Purpose: Define a behavior genetics ACE model as a 
#  Relational SEM (rSEM).  This is done for one phenotype (height).
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
require(OpenMx)


#------------------------------------------------------------------------------
# Prepare Data
# Take wide format data, put it into tall format for multilevel model.
# Also, add coefficient of relatedness "rel" and its inverse "relu"

# Wide to tall
data("twinData", package="OpenMx")
selVars <- c('ht1', 'ht2', 'zyg')
wideData <- subset(twinData, zyg %in% c(1, 3), selVars)
wideData$rel <- c(1, NA, .5)[wideData$zyg]
tallData <- reshape(wideData, varying=list(c('ht1', 'ht2')), v.names=c('ht'), timevar='twin', times=1:2, idvar='fam', direction='long')

# Add a couple new variables, sort, and drop vars that are not needed.
tallData$personID <- 1:nrow(tallData)
tallData$relu <- 1-tallData$rel
tallData <- tallData[order(tallData$fam, tallData$twin), c('fam', 'personID', 'twin', 'rel', 'relu', 'ht')]

head(tallData)

# Split data into "between" and "within" data sets.
wData <- tallData
bData <- tallData[!duplicated(tallData$fam), c('fam', 'rel', 'relu')]

dim(wData)
dim(bData)


dzData <- wideData[wideData$zyg == 3,]
mzData <- wideData[wideData$zyg == 1,]


#------------------------------------------------------------------------------
# Between Model
#Version 1 that writes everything out

bModel1 <- mxModel('between', type="RAM",
                  mxData(type="raw", observed=bData, primaryKey="fam"),
                  latentVars = c("C1", "AC1"),
                  mxPath(c("C1", "AC1"), arrows=2, values=1,
                      free=FALSE, labels=c(NA, "data.rel")))


dn1 <- list(c("C1", "AC1"), c("C1", "AC1"))
bModel1Matrix <- mxModel('mbetween',
                         mxData(type="raw", observed=bData, primaryKey="fam"),
                         mxMatrix(name='A', type='Zero', nrow=2, ncol=2),
                         mxMatrix(name='S', type='Symm', nrow=2, ncol=2, labels=c(NA, NA, 'data.rel'), values=c(1, 0, 1)),
                         mxMatrix(name='F', type='Zero', nrow=0, ncol=2),
                         mxExpectationRAM(A='A', S='S', F='F', dimnames=c("C1", "AC1")),
                         mxFitFunctionML()
                        )


#------------------------------------------------------------------------------
# Within Model
#Version 1 that writes everything out


wModel1 <- mxModel('within', type="RAM", bModel1,
                  mxData(type="raw", observed=wData, sort=FALSE),
                  manifestVars = c('ht'),
                  latentVars = c("E1", "AU1"),
                  mxPath(from="one", to=c('ht'), arrows=1, free=TRUE, values=1.5, labels=c("meanHt")),
                  mxPath(c("E1", "AU1"), arrows=2, values=1,
                      free=FALSE, labels=c(NA, "data.relu")),
                  mxPath('AU1', c('ht'), values=.8, labels=c('a1Ht'), lbound=1e-6),
                  mxPath('E1', c('ht'), values=.8, labels=c('e1Ht'), lbound=1e-6),
                  mxPath('between.C1', c('ht'), values=.8,
                         labels=c('c1Ht'), lbound=1e-6, joinKey="fam"),
                  mxPath('between.AC1', c('ht'), values=.8, arrows=1,
                         labels=c('a1Ht'), lbound=1e-6, joinKey="fam")
)

dn2 <- list(c("ht", "AU1", "E1"), c("ht", "AU1", "E1"))
wModel1Matrix <- mxModel('mwithin', bModel1Matrix,
                         mxData(type="raw", observed=wData, sort=FALSE),
                         mxMatrix(name="A", type="Full", nrow=3, ncol=3,
                           values=c(0, .8, .8,
                                    0, 0, 0,
                                    0, 0, 0),
                           free=c(FALSE, TRUE, TRUE,
                                  FALSE, FALSE, FALSE,
                                  FALSE, FALSE, FALSE),
                           labels=c(NA, "a1Ht", "e1Ht", NA, NA, NA, NA, NA, NA),
                           byrow=TRUE, dimnames=dn2),
                         mxMatrix(name="S", type="Diag", nrow=3, ncol=3, label=c(NA, "data.relu", NA), values=c(0, 1, 1), , dimnames=dn2),
                         mxMatrix(name="F", type="Full", nrow=1, ncol=3, values=c(1, 0, 0)),#, dimnames=list(c("ht"), c("ht", "AU1", "E1"))),
                         mxMatrix(name="M", type="Full", nrow=1, ncol=3, values=c(1.5, 0, 0), labels=c("meanHt", NA, NA), free=c(TRUE, FALSE, FALSE), dimnames=list(NULL, c("ht", "AU1", "E1"))),
                         mxMatrix(name="cross_model", type="Full", nrow=3, ncol=2,
                           values=c(.8, .8,
                                     0, 0,
                                     0, 0),
                           free=c(TRUE, TRUE,
                                  FALSE, FALSE,
                                  FALSE, FALSE),
                           labels=c("c1Ht", "a1Ht", NA, NA, NA, NA),
                           byrow=TRUE, dimnames=list(c("ht", "AU1", "E1"), c("C1", "AC1")),
                           joinKey="fam", joinModel="mbetween"),
                         mxExpectationRAM(A='A', S='S', F='F', M='M', between="cross_model", dimnames=c("ht", "AU1", "E1")),
                         mxFitFunctionML()
)


#------------------------------------------------------------------------------
# Run 'em
wRun1 <- mxRun(wModel1)
# Works fine

wRun1Matrix <- mxRun(wModel1Matrix)
# Fails
# Error: Join mapping matrix mwithin.cross_model must have 0 columns: NULL
