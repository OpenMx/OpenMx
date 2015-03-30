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
  omxMatrix *cov, *means, *X, *y;
  int alwaysComputeMeans, numcases2drop, cholV_fail, cholquadX_fail;
  std::vector< int > dropcase;
  Eigen::VectorXd cholV_vectorD;
  Eigen::VectorXd cholquadX_vectorD;
  Eigen::MatrixXd Vinv, XtVinv, quadXinv, b;
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

static double omxAliasedMatrixElement(omxMatrix *om, int row, int col)
{
  int index = 0;
  if(row >= om->originalRows || col >= om->originalCols) {
		char *errstr = (char*) calloc(250, sizeof(char));
		sprintf(errstr, "Requested improper value (%d, %d) from (%d, %d) matrix.", 
			row + 1, col + 1, om->originalRows, om->originalCols);
		Rf_error(errstr);
		free(errstr);  // TODO not reached
        return (NA_REAL);
	}
	if(om->colMajor) {
		index = col * om->originalRows + row;
	} else {
		index = row * om->originalCols + col;
	}
	return om->data[index];
}

