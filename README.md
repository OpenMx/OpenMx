# OpenMx

Developers commit to the `master` branch.
These commits should be tested using `make test`,
which runs all the tests in `models/passing`, and
also `make cran-check`.
However, bugs can slip through the cracks.
Therefore, the `master` branch should only
be installed by intrepid developers.

The `stable` branch is updated by our buildbot.
The buildbot will update the `stable` branch
only when all `models/passing` and
`models/nightly` tests pass.
It will also ensure that `make cran-check` passes.
The `stable` branch can be considered our
current alpha release.
An easy way to install the `stable` branch is
to use `devtools`:

```R
require(devtools)
install_github("OpenMx/OpenMx", ref="stable")
```
