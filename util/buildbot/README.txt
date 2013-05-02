At UVa, we have a Linux supercomputer cluster. Job are submitted using
qsub. The OpenMx buildbot takes advantage of this by running 4
buildbot slaves on the supercomputer front-end. These slaves do
nothing except submit jobs to the supercomputer.
