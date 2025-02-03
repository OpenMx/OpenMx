 /*
 *  Copyright 2007-2021 by the individuals mentioned in the source code history
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
  omxMatrix *cov, *invcov, *means, *X, *logdetV_om;
  omxData *y, *data2;
  int origVdim, numcases2drop; 
  bool alwaysComputeMeans, cholquadX_fail, cholV_fail, doREML, didUserProvideYhat;
  std::vector< bool > dropcase;
  Eigen::VectorXd cholV_vectorD;
  Eigen::VectorXd cholquadX_vectorD;
  Eigen::MatrixXd XtVinv, quadXinv, EigV_filtered;
  std::vector< const char* > yXcolnames;

	omxGREMLExpectation(omxState *st, int num) :
    super(st, num), cov(0), invcov(0), means(0), X(0), logdetV_om(0),
    y(0), data2(0) {}
  virtual ~omxGREMLExpectation();
  virtual void init() override;
  virtual void connectToData() override;
  virtual void compute(FitContext *fc, const char *what, const char *how) override;
  virtual void populateAttr(SEXP expectation) override;
  virtual omxMatrix *getComponent(const char*) override;
	virtual bool usesDataColumnNames() const override { return false; }
	virtual int numObservedStats() override { return 1; }
};

double omxAliasedMatrixElement(omxMatrix *om, int row, int col, int origDim);

template <typename T1>
void dropCasesAndEigenize(
		/*An omxMatrix, from which rows and columns corresponding to missing observations must be dropped,
		  and which must be made usable with the Eigen API:*/
		omxMatrix* om,
		//An Eigen object of appropriate type.  If `om` comes from a frontend MxMatrix, its elements will be copied into `em`:
		Eigen::MatrixBase<T1> &em,
		double* &ptrToMatrix, //<--Pointer to data array, for use with Eigen Map subsequent to this function call.
		int num2drop, //<--How many rows and columns to drop.
		std::vector< bool > &todrop, //<--Should row and column i be dropped?
		bool symmetric, //<--Is the matrix symmetric?
		int origDim, //<--Original pre-filtering dimensions of the matrix.
		bool copyAlg //<--Should the elements of `om` be copied to `em`, even if `om` came from an algebra?
		){

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
		ptrToMatrix = em.derived().data();
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
 if(copyAlg){
 	//Note that this line serves to copy `om`'s data to `em`.  It does not change the type of `em`:
 	em = Eigen::Map< Eigen::MatrixXd >(om->data, om->rows, om->cols);
 }
 ptrToMatrix = omxMatrixDataColumnMajor(om);
 omxMarkDirty(om); //<--Need to mark it dirty so that it eventually gets recalculated back to original dimensions.
 //^^^Algebras that do not depend upon free parameters, and upon which V does not depend, will not be
 //recalculated back to full size until optimization is complete (the GREML fitfunction is smart about that).
	}
	if(OMX_DEBUG) { mxLog("Finished trimming out cases with missing data..."); }
}

void dropCasesFromAlgdV(omxMatrix* om, int num2drop, std::vector< bool > &todrop, int symmetric, int origDim);
