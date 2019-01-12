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

#include <sys/stat.h>

#include "omxDefines.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "glue.h"
#include "omxExportBackendState.h"
#include "EnableWarnings.h"

void omxState::omxExportResults(MxRList *out, FitContext *fc)
{
	SEXP matrices;
	SEXP algebras;
	SEXP datums;

	loadDefinitionVariables(false);

	Rf_protect(matrices = Rf_allocVector(VECSXP, matrixList.size()));
	Rf_protect(algebras = Rf_allocVector(VECSXP, algebraList.size()));
	Rf_protect(datums = Rf_allocVector(VECSXP, dataList.size()));

	SEXP nextMat, algebra;
	for(size_t index = 0; index < matrixList.size(); index++) {
		if(OMX_DEBUG) { mxLog("Final Calculation and Copy of Matrix %d.", (int) index); }
		omxMatrix* nextMatrix = matrixList[index];
		nextMat = omxExportMatrix(nextMatrix);
		SET_VECTOR_ELT(matrices, index, nextMat);
	}

	setWantStage(FF_COMPUTE_FIT | FF_COMPUTE_FINAL_FIT);

	for(size_t index = 0; index < algebraList.size(); index++) {
		if(OMX_DEBUG) { mxLog("Final Calculation and Copy of Algebra %d.", (int) index); }
		omxMatrix* nextAlgebra = algebraList[index];
		// If a model has algebra that depend on free parameters
		// but the fitfunction does not depend on those algebra
		// then they need to be recomputed based on the final
		// estimates.
		if (!isErrorRaised()) omxRecompute(nextAlgebra, fc);
		algebra = omxExportMatrix(nextAlgebra);
		/* If an fit function, populate attributes.  Will skip if not fit function. */
		omxFitFunction* currentFit = nextAlgebra->fitFunction;
		if(currentFit != NULL) {
			if(OMX_DEBUG) { mxLog("Algebra %d is a fit function.", (int) index); }
			currentFit->populateAttr(algebra);
		}

		if(OMX_DEBUG) { mxLog("Final Calculation of Algebra %d Complete.", (int) index); }
		SET_VECTOR_ELT(algebras, index, algebra);
	}
	if(OMX_DEBUG) { mxLog("All Algebras complete."); }
	
	for(size_t index = 0; index < dataList.size(); ++index) {
		omxData* dat = dataList[index];
		MxRList tmp;
		dat->reportResults(tmp);
		SET_VECTOR_ELT(datums, index, tmp.asR());
	}

	out->add("matrices", matrices);
	out->add("algebras", algebras);
	out->add("data", datums);
}

void omxPopulateFitFunction(omxMatrix *om, MxRList *result) // deprecated
{
	omxFitFunction* off = om->fitFunction;
	if (!off) return;
	off->addOutput(result);
}
