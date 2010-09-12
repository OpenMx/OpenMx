/*
 *  Copyright 2007-2009 The OpenMx Project
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

#include "omxObjectiveMetadata.h"

/* Specifics for RAM FIML */
extern void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);

void omxInitRAMMetaData(SEXP rObj,  omxObjectiveMetadataContainer* omc, omxState* currentState);

const omxObjectiveMetadata omxObjectiveMetadataTable[numObjectiveMetadatas] = 
							{{"MxRAMMetaData", &omxInitRAMMetaData}};

/* FIML RAM MetaData Pieces */
typedef struct omxRAMMetadata {

	omxMatrix *cov, *means, *I;
	omxMatrix *A, *S, *F, *M;
	omxMatrix *C, *X, *Y, *Z, *Ax;

	int numIters;

} omxRAMMetadata;

void omxCallRAMSubObjective(void* subObjective, omxMatrix* cov, omxMatrix* means) {
	omxRAMMetadata* oro = (omxRAMMetadata*)(subObjective);
	
	omxRecompute(oro->A);
	omxRecompute(oro->S);
	omxRecompute(oro->F);
	if(oro->M != NULL)
	    omxRecompute(oro->M);

	omxCalculateRAMCovarianceAndMeans(oro->A, oro->S, oro->F, oro->M, cov, means, oro->numIters, oro->I, oro->Z, oro->Y, oro->X, oro->Ax);
}

void omxDestroyRAMSubObjective(void* subObjective, omxObjectiveMetadataContainer* omc) {
	if(OMX_DEBUG) { Rprintf("Destroying RAM Metadata.\n"); }
	
	omxRAMMetadata* argStruct = (omxRAMMetadata*)(subObjective);

	/* We allocated 'em, so we destroy 'em. */
	omxFreeMatrixData(omc->cov);

	if(omc->means != NULL)
		omxFreeMatrixData(omc->means);

	omc->cov = NULL;
	omc->means = NULL;
	
	omxFreeMatrixData(argStruct->I);
	omxFreeMatrixData(argStruct->X);
	omxFreeMatrixData(argStruct->Y);
	omxFreeMatrixData(argStruct->Z);
	omxFreeMatrixData(argStruct->Ax);
}

void omxInitRAMMetaData(SEXP rObj,  omxObjectiveMetadataContainer* omc, omxState* currentState) {
	int l, k;

	SEXP slotValue;

	omxRAMMetadata *newObj = (omxRAMMetadata*) R_alloc(1, sizeof(omxRAMMetadata));

	if(OMX_DEBUG) { Rprintf("Initializing RAM Meta Data for objective function.\n"); }

	if(OMX_DEBUG) { Rprintf("Processing M.\n"); }
	newObj->M = omxNewMatrixFromIndexSlot(rObj, currentState, "M");

	if(OMX_DEBUG) { Rprintf("Processing A.\n"); }
	newObj->A = omxNewMatrixFromIndexSlot(rObj, currentState, "A");

	if(OMX_DEBUG) { Rprintf("Processing S.\n"); }
	newObj->S = omxNewMatrixFromIndexSlot(rObj, currentState, "S");

	if(OMX_DEBUG) { Rprintf("Processing F.\n"); }
	newObj->F = omxNewMatrixFromIndexSlot(rObj, currentState, "F");

	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { Rprintf("Generating I.\n"); }
	newObj->I = omxNewIdentityMatrix(newObj->A->rows, currentState);
	omxRecompute(newObj->I);

	if(OMX_DEBUG) { Rprintf("Processing expansion iteration depth.\n"); }
	PROTECT(slotValue = GET_SLOT(rObj, install("depth")));
	newObj->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { Rprintf("Using %d iterations.", newObj->numIters); }
	UNPROTECT(1);

	l = newObj->F->rows;
	k = newObj->A->cols;

	if(OMX_DEBUG) { Rprintf("Generating internals for computation.\n"); }
			
	newObj->Z = 	omxInitMatrix(NULL, k, k, TRUE, currentState);
	newObj->Ax = 	omxInitMatrix(NULL, k, k, TRUE, currentState);
	newObj->Y = 	omxInitMatrix(NULL, l, k, TRUE, currentState);
	newObj->X = 	omxInitMatrix(NULL, l, k, TRUE, currentState);
	
	if(omc->cov == NULL)
		omc->cov = 		omxInitMatrix(NULL, l, l, TRUE, currentState);
	newObj->cov = 	omc->cov;
	if(newObj->M != NULL) {
		if(omc->means ==  NULL)
			omc->means = 	omxInitMatrix(NULL, 1, l, TRUE, currentState);
	} else omc->means = NULL;
	newObj->means = omc->means;

	omc->subObjective = newObj;
	/* Register functions */
	omc->destroySubObjective = omxDestroyRAMSubObjective;
	omc->covarianceMeansFunction = omxCallRAMSubObjective;

}
