# OpenMx ![Build Status](https://travis-ci.org/OpenMx/OpenMx.svg?branch=master)

A Structural Equation Modeling package that encourages users to treat model
specifications as something to be generated and manipulated programmatically.

Example models which OpenMx can fit include confirmatory factor, 
multiple group, mixture distribution, categorical threshold, 
modern test theory, differential equations, state space, and many others.

Models may be specified as RAM or LISREL mxPaths, or 
directly in matrix algebra.

Fit functions include ML (summary and full information) and WLS.

Developers commit to the `master` branch.  Intrepid users are
encouraged to install the `master` branch.

The `stable` branch can be considered our
current alpha release.

The `stable` branch is updated automatically when all `models/passing`
and `models/nightly` tests pass along with `make cran-check`.
