# OpenMx

<!-- badges: start -->
[![Build Status](https://travis-ci.org/OpenMx/OpenMx.svg?branch=master)](https://travis-ci.org/OpenMx/OpenMx)
[![Codecov test coverage](https://codecov.io/gh/OpenMx/OpenMx/branch/master/graph/badge.svg)](https://codecov.io/gh/OpenMx/OpenMx?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/OpenMx)](https://cran.r-project.org/package=OpenMx)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/OpenMx)](https://cranlogs.r-pkg.org/badges/OpenMx)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/OpenMx)](https://cranlogs.r-pkg.org/badges/grand-total/OpenMx)
[![Rdoc](https://www.rdocumentation.org/badges/version/OpenMx)](https://www.rdocumentation.org/packages/OpenMx)
[![DOI](https://img.shields.io/badge/doi-110.1007/s11336--014--9435--8-yellow.svg?style=flat)](https://doi.org/10.1007/s11336-014-9435-8)
<!-- badges: end -->

OpenMx is a [Structural Equation Modeling](https://en.wikipedia.org/wiki/Structural_equation_modeling)
package that encourages users to treat model specifications as something to be generated
and manipulated programmatically.

Example models which OpenMx can fit include everything from confirmatory factor,
through multiple group, mixture distribution, categorical threshold,
modern test theory, differential equations, state space, and many others.

Models may be specified as RAM or LISREL paths, or directly in matrix algebra.

Fit functions include ML (summary and full information) and WLS.

The package is on CRAN, so the easiest way to install is just:

`install.packages("OpenMx")`

**Website**: [https://openmx.ssri.psu.edu](https://openmx.ssri.psu.edu) (with forums, example models, documentation)

#### Development versions

Developers commit to the `master` branch.  Intrepid users are encouraged to install the `master` branch.

On Mac OS, this can be installed as a binary via travis:

```R

install.packages("https://vipbg.vcu.edu/vipbg/OpenMx2/software/bin/macosx/travis/OpenMx_latest.tgz")

```

The `stable` branch can be considered our current alpha release.

The `stable` branch is updated automatically when all `models/passing`
and `models/nightly` tests pass along with `make cran-check`.
