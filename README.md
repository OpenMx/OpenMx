# OpenMx <img src="https://openmx.ssri.psu.edu/sites/default/files/guineapig.png" align="right" width="120" />

<!-- badges: start -->
[![Build Status](https://api.travis-ci.com/OpenMx/OpenMx.svg?branch=master)](https://app.travis-ci.com/github/OpenMx/OpenMx)
[![Codecov test coverage](https://codecov.io/gh/OpenMx/OpenMx/branch/master/graph/badge.svg)](https://app.codecov.io/gh/OpenMx/OpenMx?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/OpenMx)](https://cran.r-project.org/package=OpenMx)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/OpenMx)](https://cranlogs.r-pkg.org/badges/OpenMx)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/OpenMx)](https://cranlogs.r-pkg.org/badges/grand-total/OpenMx)
[![DOI](https://img.shields.io/badge/doi-110.1007/s11336--014--9435--8-yellow.svg?style=flat)](https://doi.org/10.1007/s11336-014-9435-8)
<!-- badges: end -->

> OpenMx is a [Structural Equation Modeling](https://en.wikipedia.org/wiki/Structural_equation_modeling)
package that encourages users to treat model specifications as something to be generated
and manipulated programmatically.

# TOC

- [Overview](#overview)
- [Installation](#installation)
  * [Development versions](#development-versions)
- [Documentation](#documentation)
- [Quick usage examples](#quick-usage-examples)
  * [Simple regression in path specification](#simple-regression-in-path-specification)
  * [Simple regression using matrix algebra](#simple-regression-using-matrix-algebra)
- [Related work](#related-work)
- [Community and getting help](#community-and-getting-help)
  * [Training](#training)
- [Contributing](#contributing)


# Overview

OpenMx is the next generation of the Mx structural equation modeling tool. It is an R package activelly maintained and supported with the work of developers around the globe. It is designed to allow the user the most freedom possible while specifying structural equation models, therefore providing minimal defaults. This helps the user know that each specification/optimization decision comes with their own assumptions and influences model interpretation.


# Installation

The package is on CRAN and should be installed with:

```r
install.packages("OpenMx")
```

## Development versions

Developers commit to the `master` branch.  Intrepid users are encouraged to install the `master` branch. In order to install locally clone this repo and run:

```r
make cran-install # for the CRAN version
make install      # for the version with the proprietary NPSOL optimizer
```


The `stable` branch can be considered our current alpha release.

The `stable` branch is updated automatically when all `models/passing`
and `models/nightly` tests pass along with `make cran-check`.


On macOS, this can be installed as a binary via travis:

```r

install.packages("https://vipbg.vcu.edu/vipbg/OpenMx2/software/bin/macosx/travis/OpenMx_latest.tgz")

```


# Documentation


OpenMx can fit  everything from confirmatory factor analyses,
through multiple group, mixture distribution, categorical threshold,
modern test theory, differential equations, state space, and many others. Models may be specified as RAM or LISREL paths, or directly in matrix algebra. Fit functions include ML (summary and full information) and WLS.


The package manual can be accessed online in the [link](https://vipbg.vcu.edu/vipbg/OpenMx2/docs//OpenMx/latest/).



# Quick usage examples


Path specifications are matematically complete and is often considered an easier approach to teaching and analysis. The path below represents a simple regression:

![path](https://vipbg.vcu.edu/vipbg/OpenMx2/docs//OpenMx/latest/_images/SimpleRegression.png)



## Simple regression in path specification

One can specify the above model using the following code:

```r
require(OpenMx)

data(myRegDataRaw)  # load data
names(myRegDataRaw)  # get names
SimpleDataRaw <- myRegDataRaw[,c("x","y")]  # take only what is needed

dataRaw      <- mxData( observed=SimpleDataRaw,  type="raw" )
# variance paths
varPaths     <- mxPath( from=c("x","y"), arrows=2,
                        free=TRUE, values = c(1,1), labels=c("varx","residual") )
# regression weights
regPaths     <- mxPath( from="x", to="y", arrows=1,
                        free=TRUE, values=1, labels="beta1" )
# means and intercepts
means        <- mxPath( from="one", to=c("x","y"), arrows=1,
                        free=TRUE, values=c(1,1), labels=c("meanx","beta0") )

uniRegModel  <- mxModel(model="Simple Regression Path Specification", type="RAM",
                        dataRaw, manifestVars=c("x","y"), varPaths, regPaths, means)

uniRegFit <- mxRun(uniRegModel) # run it
summary(uniRegFit)
```

And the following output should appear in your R environment:

```
Summary of Simple Regression Path Specification

free parameters:
      name matrix row col   Estimate  Std.Error A
1    beta1      A   y   x 0.48311962 0.07757687
2     varx      S   x   x 1.10531952 0.15631652
3 residual      S   y   y 0.66520320 0.09407411
4    meanx      M   1   x 0.05415975 0.10513428
5    beta0      M   1   y 2.54776414 0.08166814

Model Statistics:
               |  Parameters  |  Degrees of Freedom  |  Fit (-2lnL units)
       Model:              5                    195              536.8226
   Saturated:              5                    195                    NA
Independence:              4                    196                    NA
Number of observations/statistics: 100/200

Information Criteria:
      |  df Penalty  |  Parameters Penalty  |  Sample-Size Adjusted
AIC:       146.8226               546.8226                 547.4609
BIC:      -361.1856               559.8484                 544.0572
CFI: NA
TLI: 1   (also known as NNFI)
RMSEA:  0  [95% CI (NA, NA)]
Prob(RMSEA <= 0.05): NA
To get additional fit indices, see help(mxRefModels)
timestamp: 2022-05-01 09:53:24
```

## Simple regression using matrix algebra

Since OpenMx is considered the specialist tool, you are probably more interested in the flexibility provided by the fact that you can build your own formulas. So going back to the simple regression, now in the formula (equivalent to the path specified in previous section):

![simple regression](https://vipbg.vcu.edu/vipbg/OpenMx2/docs//OpenMx/latest/_images/math/363ea6ab84c0c97a7d183a4621616c7753040acb.png)

It can be implemented with the following code:


```r
require(OpenMx)

data(myRegDataRaw)  # load data
SimpleDataRaw <- myRegDataRaw[,c("x","y")]  # take only what is needed

# create a data object
dataRaw      <- mxData( observed=SimpleDataRaw, type="raw" )

# A matrix
matrA        <- mxMatrix( type="Full", nrow=2, ncol=2,
                          free=c(F,F,T,F), values=c(0,0,1,0),
                          labels=c(NA,NA,"beta1",NA), byrow=TRUE, name="A" )

# S matrix
matrS        <- mxMatrix( type="Symm", nrow=2, ncol=2,
                          free=c(T,F,F,T), values=c(1,0,0,1),
                          labels=c("varx",NA,NA,"residual"), byrow=TRUE, name="S" )

# Filter matrix
matrF        <- mxMatrix( type="Iden", nrow=2, ncol=2, name="F" )

# M matrix
matrM        <- mxMatrix( type="Full", nrow=1, ncol=2,
                          free=c(T,T), values=c(0,0),
                          labels=c("meanx","beta0"), name="M")

# Which expectation? RAM in this case
expRAM       <- mxExpectationRAM("A","S","F","M", dimnames=c("x","y"))

# Run a maximum likelihood
funML        <- mxFitFunctionML()

# Name it, pass the objects on to a final model object
uniRegModel  <- mxModel("Simple Regression Matrix Specification",
                        dataRaw, matrA, matrS, matrF, matrM, expRAM, funML)

# Run it!
uniRegFit <- mxRun(uniRegModel)

summary(uniRegFit)
```


Now the output looks like:


```
Running Simple Regression Matrix Specification with 5 parameters
Summary of Simple Regression Matrix Specification

free parameters:
      name matrix row col   Estimate  Std.Error A
1    beta1      A   2   1 0.48311963 0.07757699
2     varx      S   1   1 1.10531937 0.15631498
3 residual      S   2   2 0.66520312 0.09407369
4    meanx      M   1   x 0.05416001 0.10513400
5    beta0      M   1   y 2.54776424 0.08166812

Model Statistics:
               |  Parameters  |  Degrees of Freedom  |  Fit (-2
       Model:              5                    195
   Saturated:              5                    195
Independence:              4                    196
Number of observations/statistics: 100/200

Information Criteria:
      |  df Penalty  |  Parameters Penalty  |  Sample-Size Adju
AIC:       146.8226               546.8226                 547.
BIC:      -361.1856               559.8484                 544.
CFI: NA
TLI: 1   (also known as NNFI)
RMSEA:  0  [95% CI (NA, NA)]
Prob(RMSEA <= 0.05): NA
To get additional fit indices, see help(mxRefModels)
timestamp: 2022-05-01 10:04:52
Wall clock time: 0.3507566 secs
optimizer:  SLSQP
OpenMx version number: 2.19.6.6
Need help?  See help(mxSummary)
```

# Related work

[umx()](https://github.com/tbates/umx) is a sister R package that bridges the gap between [lavaan](https://github.com/yrosseel/lavaan) and OpenMx. If you are coming from lavaan it is perhaps useful to check umx() too. [Onyx](https://onyx-sem.com/) is a software that allows you to design nice diagrams, and syncs (exports and imports) the diagrams with OpenMx code.

# Community and getting help

1. The support communication is centered around the OpenMx [forum](https://openmx.ssri.psu.edu/forums)
2. Also, but less often, at the StackOverflow OpenMx [tag](https://stackoverflow.com/questions/tagged/openmx/).


## Training

We gather annually in beautiful [Boulder, CO](https://www.colorado.edu/ibg/workshop) for the international workshop for traning in behavioral genetics applications of OpenMx.


# Contributing

<details>
<summary>How can I contribute to this project?</summary>
<br>
OpenMx is maintained by a small team and all help is appreciated.

First read the team's conduct policy [here](https://github.com/OpenMx/OpenMx/blob/master/CONTRIBUTING). If you agree with it you can choose one of the below paths:

1. Do you have a well documented script (from one of our several workshops) that would make a great vignette? Great, because you don't even need to know how to use git. Simply go to the vignette folder and click in add file. This will automate the forking and uploading. 
2. There are several issues that can be handled by new users. Go over to our oldest issues [here](https://github.com/OpenMx/OpenMx/issues?q=is%3Aissue+is%3Aopen+sort%3Acreated-asc), browse until something you find an issue you feel you can contribute to, and announce that you are planning to tackle it in the issue thread.
3. Have a completely new functionality that you want to discuss? Just create a PR and we will discuss whether it aligns with the package direction. In this case please add proper documentation for the new functionality. If you use RStudio there is a stub at File > New File > R Documentation. Also create a test unit in the tests/testthat folder, we currently use [testthat](https://testthat.r-lib.org/) to manage this. 
<br>

</details>

