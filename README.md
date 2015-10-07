# OpenMx

A Structural Equation Modeling package encouraging users to treat model
specifications as something able to be generated and manipulated programmatically.

Example models which OpenMx can fit include confirmatory factor, 
multiple group, mixture distribution, categorical threshold, 
modern test theory, differential equations, state space, and many others.

Models may be specified as RAM or LISREL mxPaths, or 
directly in matrix algebra.

Fit functions include FIML, ML and WLS.

We have an active development branch. on github.

The `stable` branch can be considered our
current alpha release.

An easy way to install the `stable` branch is
to use `devtools`:

```R
require(devtools)
install_github("OpenMx/OpenMx", ref="stable")
```

Developers commit to the `master` branch and this is accessible
to more intrepid users.

Commits should be tested using `make test`,
which runs all the tests in `models/passing`, and
also `make cran-check`. buildbot updates the `stable` branch
only when all `models/passing` and `models/nightly` tests pass
along with `make cran-check`.
