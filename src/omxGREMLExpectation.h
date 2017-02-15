 /*
 *  Copyright 2007-2017 The OpenMx Project
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
 
struct omxGREMLExpectation : public omxExpectation {
  omxMatrix *cov, *invcov, *means, *X, *logdetV_om, *cholV_fail_om, *origVdim_om;
  omxData *y, *data2;
  int alwaysComputeMeans, numcases2drop, cholquadX_fail;
  std::vector< int > dropcase;
  Eigen::VectorXd cholV_vectorD;
  Eigen::VectorXd cholquadX_vectorD;
  Eigen::MatrixXd XtVinv, quadXinv;
  std::vector< const char* > yXcolnames;

  virtual ~omxGREMLExpectation();
  virtual void init();
  virtual void compute(FitContext *fc, const char *what, const char *how);
  virtual void populateAttr(SEXP expectation);
  virtual omxMatrix *getComponent(const char*);
};

void dropCasesAndEigenize(omxMatrix* om, Eigen::MatrixXd &em, int num2drop, std::vector< int > todrop, 
	int symmetric, int origDim);
void dropCasesFromAlgdV(omxMatrix* om, int num2drop, std::vector< int > todrop, int symmetric, int origDim);

