 /*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
	typedef omxExpectation super;
  omxMatrix *cov, *invcov, *means, *X, *logdetV_om, *cholV_fail_om, *origVdim_om;
  omxData *y, *data2;
  int alwaysComputeMeans, numcases2drop, cholquadX_fail;
  std::vector< int > dropcase;
  Eigen::VectorXd cholV_vectorD;
  Eigen::VectorXd cholquadX_vectorD;
  Eigen::MatrixXd XtVinv, quadXinv;
  std::vector< const char* > yXcolnames;

	omxGREMLExpectation(omxState *st, int num) : super(st, num) {}
  virtual ~omxGREMLExpectation();
  virtual void init();
  virtual void connectToData();
  virtual void compute(FitContext *fc, const char *what, const char *how);
  virtual void populateAttr(SEXP expectation);
  virtual omxMatrix *getComponent(const char*);
	virtual bool usesDataColumnNames() const { return false; }
};

double omxAliasedMatrixElement(omxMatrix *om, int row, int col, int origDim);

template <typename T1>
void dropCasesAndEigenize(omxMatrix* om, Eigen::MatrixBase<T1> &em, int num2drop, std::vector< int > &todrop,
                          int symmetric, int origDim){

	if(OMX_DEBUG) { mxLog("Trimming out cases with missing data..."); }

	if(num2drop < 1){ return; }

	omxEnsureColumnMajor(om);

	if(om->algebra == NULL){ //i.e., if omxMatrix is from a frontend MxMatrix

		em.derived().setZero(om->rows - num2drop, om->cols - num2drop);

		int nextCol = 0;
		int nextRow = 0;

		for(int j = 0; j < om->cols; j++) {
			if(todrop[j]) continue;
			nextRow = (symmetric ? nextCol : 0);
			for(int k = (symmetric ? j : 0); k < om->rows; k++) {
				if(todrop[k]) continue;
				em(nextRow,nextCol) = omxAliasedMatrixElement(om, k, j, origDim);
				nextRow++;
			}
			nextCol++;
		}
	}
	else{ /*If the omxMatrix is from an algebra, then copying is not necessary; it can be resized directly
		and Eigen-mapped, since the algebra will be recalculated back to its original dimensions anyhow.*/
 if(origDim==0){mxThrow("Memory not allocated for algebra %s at downsize time",
    om->name());}
 if(om->rows != origDim || om->cols != origDim){
 	//Not sure if there are cases where this should be allowed
 	mxThrow("More than one attempt made to downsize algebra %s", om->name());
 	//return;
 }

 //int oldRows = om->originalRows;
 //int oldCols = om->originalCols;

 int nextCol = 0;
 int nextRow = 0;

 om->rows = origDim - num2drop;
 om->cols = origDim - num2drop;

 for(int j = 0; j < origDim; j++){ //<--j indexes columns
 	if(todrop[j]) continue;
 	nextRow = (symmetric ? nextCol : 0);
 	for(int k = (symmetric ? j : 0); k < origDim; k++){ //<--k indexes rows
 		if(todrop[k]) continue;
 		omxSetMatrixElement(om, nextRow, nextCol, omxAliasedMatrixElement(om, k, j, origDim));
 		nextRow++;
 	}
 	nextCol++;
 }
 em = Eigen::Map< Eigen::MatrixXd >(om->data, om->rows, om->cols);
 omxMarkDirty(om); //<--Need to mark it dirty so that it eventually gets recalculated back to original dimensions.
 //^^^Algebras that do not depend upon free parameters, and upon which V does not depend, will not be
 //recalculated back to full size until optimization is complete (the GREML fitfunction is smart about that).
	}
	if(OMX_DEBUG) { mxLog("Finished trimming out cases with missing data..."); }
}

void dropCasesFromAlgdV(omxMatrix* om, int num2drop, std::vector< int > &todrop, int symmetric, int origDim);
