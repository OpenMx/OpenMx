 /*
 *  Copyright 2007-2015 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
 
 typedef struct {
  omxMatrix *cov, *invcov, *means, *X, *logdetV_om, *cholV_fail_om;
  omxData *y, *data2;
  int alwaysComputeMeans, numcases2drop, cholquadX_fail;
  std::vector< int > dropcase;
  //double logdetV;
  Eigen::VectorXd cholV_vectorD;
  Eigen::VectorXd cholquadX_vectorD;
  Eigen::MatrixXd XtVinv, quadXinv;
  std::vector< const char* > yXcolnames;
  //const char **yXcolnames;
  //std::vector< omxMatrix* > dV;
  //std::vector< const char* > dVnames;
  //int dVlength;
} omxGREMLExpectation;

void omxInitGREMLExpectation(omxExpectation* ox);
void omxComputeGREMLExpectation(omxExpectation* ox, const char *, const char *);
void omxDestroyGREMLExpectation(omxExpectation* ox);
void omxPopulateGREMLAttributes(omxExpectation *ox, SEXP algebra);
void dropCasesAndEigenize(omxMatrix* om, Eigen::MatrixXd &em, int num2drop, std::vector< int > todrop);
omxMatrix* omxGetGREMLExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component);

