
Mixture Model
====================

This example will demonstrate how to specify a growth mixture model using path specification. Unlike other examples, this application will not be demonstrated with covariance data, as this model can only be fit to raw data. The script for this example can be found in the following file:

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/GrowthMixtureModel_Path.R

A parallel example using matrix specification can be found here:

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/GrowthMixtureModel_Matrix.R

The latent growth curve used in this example is the same one fit in the latent growth curve example. Path and matrix versions of that example for raw data can be found here: 

** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/LatentGrowthCurveModel_PathRaw.R
** http://openmx.psyc.virginia.edu/repoview/1/trunk/demo/LatentGrowthCurveModel_MatrixRaw.R

Mixture Model
-------------

Mixture models

.. math::
   :nowrap:
   
   \begin{eqnarray*}

x_{ij} = p_1 (Intercept_{i1} + \lambda_1 Slope_{i1} + \epsilon_{i1}) + p_2 (Intercept_{i2} + \lambda_2 Slope_{i2} + \epsilon_{i2}) \\ 
\sum_{i=1}^k p_i = 1 
   \end{eqnarray*}

Data
^^^^



Model Specification
^^^^^^^^^^^^^^^^^^^