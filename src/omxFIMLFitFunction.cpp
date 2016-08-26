/*
 *  Copyright 2007-2016 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxDefines.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"

#pragma GCC diagnostic warning "-Wshadow"

/* FIML Function body */
void omxDestroyFIMLFitFunction(omxFitFunction *off) {
	if(OMX_DEBUG) { mxLog("Destroying FIML fit function object."); }
	omxFIMLFitFunction *argStruct = (omxFIMLFitFunction*) (off->argStruct);

	omxFreeMatrix(argStruct->smallMeans);
	omxFreeMatrix(argStruct->ordMeans);
	omxFreeMatrix(argStruct->contRow);
	omxFreeMatrix(argStruct->ordCov);
	omxFreeMatrix(argStruct->ordContCov);
	omxFreeMatrix(argStruct->halfCov);
	omxFreeMatrix(argStruct->reduceCov);

	omxFreeMatrix(argStruct->smallRow);
	omxFreeMatrix(argStruct->smallCov);
	omxFreeMatrix(argStruct->RCX);
	omxFreeMatrix(argStruct->rowLikelihoods);
	delete argStruct;
}

void omxPopulateFIMLAttributes(omxFitFunction *off, SEXP algebra) {
	if(OMX_DEBUG) { mxLog("Populating FIML Attributes."); }
	omxFIMLFitFunction *argStruct = ((omxFIMLFitFunction*)off->argStruct);
	SEXP expCovExt, expMeanExt, rowLikelihoodsExt;
	omxMatrix *expCovInt, *expMeanInt;
	expCovInt = argStruct->cov;
	expMeanInt = argStruct->means;
	
	Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
	for(int row = 0; row < expCovInt->rows; row++)
		for(int col = 0; col < expCovInt->cols; col++)
			REAL(expCovExt)[col * expCovInt->rows + row] =
				omxMatrixElement(expCovInt, row, col);
	if (expMeanInt != NULL && expMeanInt->rows > 0  && expMeanInt->cols > 0) {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
		for(int row = 0; row < expMeanInt->rows; row++)
			for(int col = 0; col < expMeanInt->cols; col++)
				REAL(expMeanExt)[col * expMeanInt->rows + row] =
					omxMatrixElement(expMeanInt, row, col);
	} else {
		Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, 0, 0));		
	}

	Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
	Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
	
	if(argStruct->populateRowDiagnostics){
		omxMatrix *rowLikelihoodsInt = argStruct->rowLikelihoods;
		Rf_protect(rowLikelihoodsExt = Rf_allocVector(REALSXP, rowLikelihoodsInt->rows));
		for(int row = 0; row < rowLikelihoodsInt->rows; row++)
			REAL(rowLikelihoodsExt)[row] = omxMatrixElement(rowLikelihoodsInt, row, 0);
		Rf_setAttrib(algebra, Rf_install("likelihoods"), rowLikelihoodsExt);
	}
}

struct FIMLCompare {
	omxData *data;
	omxExpectation *ex;

	bool compareData(int la, int ra, bool &mismatch) const
	{
		mismatch = true;
		auto dc = ex->getDataColumns();
		for (int cx=0; cx < dc.size(); ++cx) {
			int col = dc[cx];
			if (!omxDataColumnIsFactor(data, col)) continue;
			double lv = omxDoubleDataElement(data, la, col);
			double rv = omxDoubleDataElement(data, ra, col);
			if (doubleEQ(lv, rv)) continue;
			return lv < rv;
		}
		for (int cx=0; cx < dc.size(); ++cx) {
			int col = dc[cx];
			if (omxDataColumnIsFactor(data, col)) continue;
			double lv = omxDoubleDataElement(data, la, col);
			double rv = omxDoubleDataElement(data, ra, col);
			if (doubleEQ(lv, rv)) continue;
			return lv < rv;
		}

		mismatch = false;
		return false;
	}

	bool compareMissingness(int la, int ra, bool &mismatch) const
	{
		mismatch = true;
		auto dc = ex->getDataColumns();
		for (int cx=0; cx < dc.size(); ++cx) {
			int col = dc[cx];
			bool lm = omxDataElementMissing(data, la, col);
			bool rm = omxDataElementMissing(data, ra, col);
			if (lm == rm) continue;
			return lm < rm;
		}

		mismatch = false;
		return false;
	}

	bool compareAllDefVars(int la, int ra, bool &mismatch) const
	{
		mismatch = true;

		for (size_t k=0; k < data->defVars.size(); ++k) {
			int col = data->defVars[k].column;
			double lv = omxDoubleDataElement(data, la, col);
			double rv = omxDoubleDataElement(data, ra, col);
			if (doubleEQ(lv, rv)) continue;
			return lv < rv;
		}

		mismatch = false;
		return false;
	}

	bool operator()(int la, int ra) const
	{
		bool mismatch;
		bool got = compareAllDefVars(la, ra, mismatch);
		if (mismatch) return got;
		got = compareMissingness(la, ra, mismatch);
		if (mismatch) return got;
		got = compareData(la, ra, mismatch);
		if (mismatch) return got;
		return false;
	}
};

static void recordGap(int rx, int &prev, std::vector<int> &identical)
{
	int gap = rx - prev;
	for (int gx=0; gx < gap; ++gx) identical[prev + gx] = gap - gx;
	prev = rx;
}

static void sortData(omxFitFunction *off, FitContext *fc)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);
	auto& indexVector = ofiml->indexVector;

	if (fc->isClone() || indexVector.size()) return;

	omxData *data = ofiml->data;
	indexVector.reserve(data->rows);
	for (int rx=0; rx < data->rows; ++rx) indexVector.push_back(rx);

	auto& identicalDefs = ofiml->identicalDefs;
	auto& identicalMissingness = ofiml->identicalMissingness;
	auto& identicalRows = ofiml->identicalRows;
	identicalDefs.assign(data->rows, 1);
	identicalMissingness.assign(data->rows, 1);
	identicalRows.assign(data->rows, 1);
	if (!data->needSort || ofiml->isStateSpace) return;

	if (OMX_DEBUG) mxLog("sort %s for %s", data->name, off->name());
	FIMLCompare cmp;
	cmp.data = data;
	cmp.ex = off->expectation;
	std::sort(indexVector.begin(), indexVector.end(), cmp);

	int prevDefs=0;
	int prevMissingness=0;
	int prevRows=0;
	for (int rx=1; rx < data->rows; ++rx) {
		bool mismatch;
		cmp.compareAllDefVars(indexVector[prevDefs], indexVector[rx], mismatch);
		if (mismatch) recordGap(rx, prevDefs, identicalDefs);
		cmp.compareMissingness(indexVector[prevMissingness], indexVector[rx], mismatch);
		if (mismatch) {
			recordGap(rx, prevMissingness, identicalMissingness);
		}
		cmp.compareData(indexVector[prevRows], indexVector[rx], mismatch);
		if (mismatch) recordGap(rx, prevRows, identicalRows);
	}
	recordGap(data->rows - 1, prevDefs, identicalDefs);
	recordGap(data->rows - 1, prevMissingness, identicalMissingness);
	recordGap(data->rows - 1, prevRows, identicalRows);
}

static void CallFIMLFitFunction(omxFitFunction *off, int want, FitContext *fc)
{
	omxFIMLFitFunction* ofiml = ((omxFIMLFitFunction*)off->argStruct);

	// TODO: Figure out how to give access to other per-iteration structures.
	// TODO: Current implementation is slow: update by filtering correlations and thresholds.
	
	if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) {
		if (fc->isClone()) {
			omxMatrix *pfitMat = fc->getParentState()->getMatrixFromIndex(off->matrix);
			ofiml->parent = (omxFIMLFitFunction*) pfitMat->fitFunction->argStruct;
		}
		off->openmpUser = !ofiml->isStateSpace && ofiml->rowwiseParallel != 0;
		sortData(off, fc);
		return;
	}

	if(OMX_DEBUG) mxLog("%s: joint FIML; openmpUser=%d", off->name(), off->openmpUser);

	omxMatrix* fitMatrix  = off->matrix;
	int numChildren = fc? fc->childList.size() : 0;

	omxMatrix *cov 		= ofiml->cov;
	omxMatrix *means	= ofiml->means;
	omxData* data           = ofiml->data;                            //  read-only

	int returnRowLikelihoods = ofiml->returnRowLikelihoods;   //  read-only
	omxExpectation* expectation = off->expectation;
	if (!means) {
		if (want & FF_COMPUTE_FINAL_FIT) return;
		complainAboutMissingMeans(expectation);
		return;
	}
	auto dataColumns	= expectation->getDataColumns();
	std::vector< omxThresholdColumn > &thresholdCols = expectation->thresholds;

	bool failed = false;

	if (data->defVars.size() == 0 && !ofiml->isStateSpace) {
		if(OMX_DEBUG) {mxLog("Precalculating cov and means for all rows.");}
		omxExpectationRecompute(fc, expectation);
		// MCN Also do the threshold formulae!
		
		omxMatrix* nextMatrix = expectation->thresholdsMat;
		if (nextMatrix) omxRecompute(nextMatrix, fc);
		for(int j=0; j < dataColumns.size(); j++) {
			int var = dataColumns[j];
			if (!omxDataColumnIsFactor(data, var)) continue;
			if (j < int(thresholdCols.size()) && thresholdCols[j].numThresholds > 0) { // j is an ordinal column
				failed |= !thresholdsIncreasing(nextMatrix, thresholdCols[j].column,
								thresholdCols[j].numThresholds, fc);
				for(int index = 0; index < numChildren; index++) {
					FitContext *kid = fc->childList[index];
					omxMatrix *target = kid->lookupDuplicate(nextMatrix);
					omxCopyMatrix(target, nextMatrix);
				}
			} else {
				Rf_error("No threshold given for ordinal column '%s'",
					 omxDataColumnName(data, j));
			}
		}

		for(int index = 0; index < numChildren; index++) {
			FitContext *kid = fc->childList[index];
			omxMatrix *childFit = kid->lookupDuplicate(fitMatrix);
			omxFIMLFitFunction* childOfiml = ((omxFIMLFitFunction*) childFit->fitFunction->argStruct);
			omxCopyMatrix(childOfiml->cov, cov);
			omxCopyMatrix(childOfiml->means, means);
		}
		if(OMX_DEBUG) { omxPrintMatrix(cov, "Cov"); }
		if(OMX_DEBUG) { if (means) omxPrintMatrix(means, "Means"); }
	}

	int parallelism = (numChildren == 0 || !off->openmpUser) ? 1 : numChildren;

	if (parallelism > data->rows) {
		parallelism = data->rows;
	}

	if (parallelism > 1) {
		int stride = (data->rows / parallelism);

		if (OMX_DEBUG) {
			FitContext *kid = fc->childList[0];
			omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			omxFIMLFitFunction* ofo = ((omxFIMLFitFunction*) childFit->argStruct);
			if (!ofo->parent) Rf_error("oops");
		}

#pragma omp parallel for num_threads(parallelism) reduction(||:failed)
		for(int i = 0; i < parallelism; i++) {
			FitContext *kid = fc->childList[i];
			omxMatrix *childMatrix = kid->lookupDuplicate(fitMatrix);
			omxFitFunction *childFit = childMatrix->fitFunction;
			if (i == parallelism - 1) {
				failed |= omxFIMLSingleIterationJoint(kid, childFit, ofiml->rowLikelihoods,
								      stride * i, data->rows - stride * i);
			} else {
				failed |= omxFIMLSingleIterationJoint(kid, childFit, ofiml->rowLikelihoods, stride * i, stride);
			}
		}
	} else {
		failed |= omxFIMLSingleIterationJoint(fc, off, ofiml->rowLikelihoods, 0, data->rows);
	}
	if (failed) {
		omxSetMatrixElement(off->matrix, 0, 0, NA_REAL);
		return;
	}

	if(!returnRowLikelihoods) {
		double sum = 0.0;
		// floating-point addition is not associative,
		// so we serialized the following reduction operation.
		for(int i = 0; i < data->rows; i++) {
			double val = log(omxVectorElement(ofiml->rowLikelihoods, i));
			//mxLog("[%d] %g", i, val);
			sum += val;
		}	
		if(OMX_DEBUG) {mxLog("%s: total likelihood is %3.3f", off->name(), sum);}
		omxSetMatrixElement(off->matrix, 0, 0, -2.0 * sum);
	} else {
		omxCopyMatrix(off->matrix, ofiml->rowLikelihoods);
	}
}

void omxInitFIMLFitFunction(omxFitFunction* off)
{
	if(OMX_DEBUG) {
		mxLog("Initializing FIML fit function function.");
	}
	off->canDuplicate = TRUE;

	omxExpectation* expectation = off->expectation;
	if(expectation == NULL) {
		omxRaiseError("FIML cannot fit without model expectations.");
		return;
	}

	omxFIMLFitFunction *newObj = new omxFIMLFitFunction;
	newObj->parent = 0;
	newObj->isStateSpace = strEQ(expectation->expType, "MxExpectationStateSpace");

	int numOrdinal = 0;

	omxMatrix *cov = omxGetExpectationComponent(expectation, "cov");
	if(cov == NULL) { 
		omxRaiseError("No covariance expectation in FIML evaluation.");
		return;
	}

	omxMatrix *means = omxGetExpectationComponent(expectation, "means");
	
    newObj->cov = cov;
    newObj->means = means;
    newObj->smallMeans = NULL;
    newObj->ordMeans   = NULL;
    newObj->contRow    = NULL;
    newObj->ordCov     = NULL;
    newObj->ordContCov = NULL;
    newObj->halfCov    = NULL;
    newObj->reduceCov  = NULL;
    
    off->computeFun = CallFIMLFitFunction;
	
	off->destructFun = omxDestroyFIMLFitFunction;
	off->populateAttrFun = omxPopulateFIMLAttributes;

	if(OMX_DEBUG) {
		mxLog("Accessing data source.");
	}
	newObj->data = off->expectation->data;

	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood option.");
	}
	SEXP rObj = off->rObj;
	newObj->rowwiseParallel = Rf_asLogical(R_do_slot(rObj, Rf_install("rowwiseParallel")));
	newObj->returnRowLikelihoods = Rf_asInteger(R_do_slot(rObj, Rf_install("vector")));
	newObj->rowLikelihoods = omxInitMatrix(newObj->data->rows, 1, TRUE, off->matrix->currentState);
	
	
	if(OMX_DEBUG) {
		mxLog("Accessing row likelihood population option.");
	}
	newObj->populateRowDiagnostics = Rf_asInteger(R_do_slot(rObj, Rf_install("rowDiagnostics")));


	if(OMX_DEBUG) {
		mxLog("Accessing variable mapping structure.");
	}
	auto dc = off->expectation->getDataColumns();

	if(OMX_DEBUG) {
		mxLog("Accessing Threshold matrix.");
	}
	numOrdinal = off->expectation->numOrdinal;

	omxSetContiguousDataColumns(&(newObj->contiguous), newObj->data, dc);
	
    /* Temporary storage for calculation */
    int covCols = newObj->cov->cols;
	if(OMX_DEBUG){mxLog("Number of columns found is %d", covCols);}
    // int ordCols = omxDataNumFactor(newObj->data);        // Unneeded, since we don't use it.
    // int contCols = omxDataNumNumeric(newObj->data);
    newObj->smallRow = omxInitMatrix(1, covCols, TRUE, off->matrix->currentState);
    newObj->smallCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
    newObj->RCX = omxInitMatrix(1, covCols, TRUE, off->matrix->currentState);
//  newObj->zeros = omxInitMatrix(1, newObj->cov->cols, TRUE, off->matrix->currentState);

    omxCopyMatrix(newObj->smallCov, newObj->cov);          // Will keep its aliased state from here on.
    if (means) {
	    newObj->smallMeans = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
	    omxCopyMatrix(newObj->smallMeans, newObj->means);
	    newObj->ordMeans = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
	    omxCopyMatrix(newObj->ordMeans, newObj->means);
    }
    newObj->contRow = omxInitMatrix(covCols, 1, TRUE, off->matrix->currentState);
    omxCopyMatrix(newObj->contRow, newObj->smallRow );
    newObj->ordCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
    omxCopyMatrix(newObj->ordCov, newObj->cov);

    off->argStruct = (void*)newObj;

    if(numOrdinal > 0) {
        if(OMX_DEBUG) mxLog("Ordinal Data detected");

        newObj->ordContCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->halfCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        newObj->reduceCov = omxInitMatrix(covCols, covCols, TRUE, off->matrix->currentState);
        omxCopyMatrix(newObj->ordContCov, newObj->cov);
    }
}
