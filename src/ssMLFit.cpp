 /*
 *  Copyright 2016 The OpenMx Project
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

#include "omxFitFunction.h"
#include "Compute.h"

#ifdef SHADOW_DIAG
#pragma GCC diagnostic warning "-Wshadow"
#endif

struct ssMLFitState {
	bool returnRowLikelihoods;
	bool populateRowDiagnostics;
	omxMatrix *cov;
	omxMatrix *smallRow;
	omxMatrix *contRow;
	omxMatrix *rowLikelihoods;
	omxMatrix *RCX;
};

void populate(omxFitFunction *off, SEXP algebra)
{
	ssMLFitState *argStruct = ((ssMLFitState*)off->argStruct);

	if(argStruct->populateRowDiagnostics){
		SEXP rowLikelihoodsExt;
		omxMatrix *rowLikelihoodsInt = argStruct->rowLikelihoods;
		Rf_protect(rowLikelihoodsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		for(int row = 0; row < rowLikelihoodsInt->rows; row++)
			REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);
		Rf_setAttrib(algebra, Rf_install("likelihoods"), rowLikelihoodsExt);
	}
}

static void compute(omxFitFunction *oo, int want, FitContext *fc)
{
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) return;

	ssMLFitState *state = (ssMLFitState *) oo->argStruct;
	omxExpectation *expectation = oo->expectation;
	auto dataColumns	= expectation->getDataColumns();
	omxData *data = expectation->data;
	int rowcount = data->rows;
	omxMatrix *smallRow = state->smallRow;
	omxMatrix *contRow = state->contRow;
	omxMatrix *cov = state->cov;
	omxMatrix *RCX = state->RCX;
	omxMatrix *rowLikelihoods = state->rowLikelihoods;

	omxSetExpectationComponent(expectation, "Reset", NULL);
	Eigen::VectorXi contRemove(cov->cols);

	for (int row=0; row < rowcount; ++row) {
		mxLogSetCurrentRow(row);
		
		omxDataRow(data, row, dataColumns, smallRow);
		
		omxSetExpectationComponent(expectation, "y", smallRow);

		expectation->loadDefVars(row);

		omxExpectationCompute(fc, expectation, NULL);

		int numCont = 0;
		for(int j = 0; j < dataColumns.size(); j++) {
			if (omxDataElementMissing(data, row, dataColumns[j])) {
				contRemove[j] = 1;
			} else {
				contRemove[j] = 0;
				++numCont;
			}
		}
		if (!numCont) {
			omxSetMatrixElement(rowLikelihoods, row, 0, 1.0);
			continue;
		}

		omxMatrix *smallMeans = omxGetExpectationComponent(expectation, "means");
		omxRemoveElements(smallMeans, contRemove.data());
		omxMatrix *smallCov = omxGetExpectationComponent(expectation, "inverse");
		if(OMX_DEBUG_ROWS(row)) { omxPrint(smallCov, "Inverse of Local Covariance Matrix in state space model"); }
		//Get covInfo from state space expectation
		int info = (int) omxGetExpectationComponent(expectation, "covInfo")->data[0];
		if(info!=0) {
			EigenMatrixAdaptor Ecov(smallCov);
			std::string empty = std::string("");
			std::string buf = mxStringifyMatrix("covariance", Ecov, empty);
			if (fc) fc->recordIterationError("%s: expected covariance matrix is "
							 "not positive-definite in data row %d; Details:\n%s",
							 oo->name(), row, buf.c_str());
			omxSetMatrixElement(oo->matrix, 0, 0, NA_REAL);
			return;
		}

		double determinant = *omxGetExpectationComponent(expectation, "determinant")->data;
		if(OMX_DEBUG_ROWS(row)) { mxLog("0.5*log(det(Cov)) is: %3.3f", determinant);}

		omxCopyMatrix(contRow, smallRow);
		omxRemoveElements(contRow, contRemove.data()); 	// Reduce the row to just continuous.
		double minusoned = -1.0;
		int onei = 1;
		F77_CALL(daxpy)(&(contRow->cols), &minusoned, smallMeans->data, &onei, contRow->data, &onei);

		/* Calculate Row Likelihood */
		/* Mathematically: (2*pi)^cols * 1/sqrt(determinant(ExpectedCov)) * (dataRow %*% (solve(ExpectedCov)) %*% t(dataRow))^(1/2) */
		//EigenMatrixAdaptor EsmallCov(smallCov);
		//mxPrintMat("smallcov", EsmallCov);
		//omxPrint(contRow, "contRow");
		double zerod = 0.0;
		char u = 'U';
		double oned = 1.0;
		F77_CALL(dsymv)(&u, &(smallCov->rows), &oned, smallCov->data, &(smallCov->cols), contRow->data, &onei, &zerod, RCX->data, &onei);       // RCX is the continuous-column mahalanobis distance.
		double Q = F77_CALL(ddot)(&(contRow->cols), contRow->data, &onei, RCX->data, &onei);

		double rowLikelihood = pow(2 * M_PI, -.5 * numCont) * (1.0/exp(determinant)) * exp(-.5 * Q);

		omxSetMatrixElement(rowLikelihoods, row, 0, rowLikelihood);
	}

	if (!state->returnRowLikelihoods) {
		double sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			double prob = omxVectorElement(state->rowLikelihoods, i);
			//mxLog("[%d] %g", i, -2.0 * log(prob));
			sum += log(prob);
		}	
		if(OMX_DEBUG) {mxLog("%s: total likelihood is %3.3f", oo->name(), sum);}
		omxSetMatrixElement(oo->matrix, 0, 0, -2.0 * sum);
	} else {
		omxCopyMatrix(oo->matrix, state->rowLikelihoods);
	}
}

static void destroy(omxFitFunction *oo)
{
	ssMLFitState *state = (ssMLFitState *) oo->argStruct;
	omxFreeMatrix(state->smallRow);
	omxFreeMatrix(state->contRow);
	omxFreeMatrix(state->rowLikelihoods);
	delete state;
}

void ssMLFitInit(omxFitFunction* oo)
{
	ssMLFitState *state = new ssMLFitState;
	oo->argStruct = state;
	oo->computeFun = compute;
	oo->destructFun = destroy;
	oo->populateAttrFun = populate;
	oo->openmpUser = false;
	oo->canDuplicate = true;

	state->returnRowLikelihoods = Rf_asInteger(R_do_slot(oo->rObj, Rf_install("vector")));
	state->populateRowDiagnostics = Rf_asInteger(R_do_slot(oo->rObj, Rf_install("rowDiagnostics")));

	omxExpectation *expectation = oo->expectation;

	omxState *currentState = oo->matrix->currentState;
	state->rowLikelihoods = omxInitMatrix(expectation->data->rows, 1, TRUE, currentState);
	state->cov = omxGetExpectationComponent(expectation, "cov");

	int covCols = state->cov->cols;
	state->smallRow = omxInitMatrix(1, covCols, TRUE, currentState);
	state->contRow = omxInitMatrix(covCols, 1, TRUE, currentState);
	state->RCX = omxInitMatrix(1, covCols, TRUE, currentState);
}
