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

/***********************************************************
*
*  omxState.cc
*
*  Created: Timothy R. Brick 	Date: 2009-06-05
*
*	omxStates carry the current optimization state
*
**********************************************************/

#include "omxState.h"

/* Initialize and Destroy */
	void omxInitState(omxState* state, omxState *parentState, int numThreads) {
		int i;
		state->numMats = 0;
		state->numAlgs = 0;
		state->numData = 0;
		state->numFreeParams = 0;
		if (numThreads > 1) {
			state->numChildren = numThreads;
			state->childList = (omxState**) Calloc(numThreads, omxState*);
			for(i = 0; i < numThreads; i++) {
				state->childList[i] = (omxState*) R_alloc(1, sizeof(omxState));
				omxInitState(state->childList[i], state, 1);
			}
		} else {
	        state->numChildren = 0;
			state->childList = NULL;
		}
		state->matrixList = NULL;
		state->algebraList = NULL;
        state->parentState = parentState;
        state->parentMatrix = NULL;
        state->parentAlgebra = NULL;
		state->parentConList= NULL;
		state->dataList = NULL;
		state->objectiveMatrix = NULL;
		state->hessian = NULL;
		state->conList = NULL;
		state->freeVarList = NULL;
		state->optimizerState = NULL;
		state->optimalValues = NULL;
		state->optimum = 9999999999;

		state->majorIteration = 0;
		state->minorIteration = 0;
		state->startTime = 0;
		state->endTime = 0;
		state->numCheckpoints = 0;
		state->checkpointList = NULL;
		state->chkptText1 = NULL;
		state->chkptText2 = NULL;

		state->computeCount = -1;
		state->currentRow = -1;

		state->statusCode = 0;
		strncpy(state->statusMsg, "", 1);
	}

	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList,
						omxMatrix** algebraList, omxData** dataList, omxMatrix* objective) {
		error("NYI: Can't fill a state from outside yet. Besides, do you really need a single function to do this?");
	}
	
	omxState* omxGetState(omxState* os, int stateNumber) {
		// TODO: Need to implement a smarter way to enumerate children
		if(stateNumber == 0) return os;
		if((stateNumber-1) < os->numChildren) {
			return(os->childList[stateNumber-1]);
		} else {
			// TODO: Account for unequal numbers of grandchild states
			int subState = (stateNumber - os->numChildren - 1);
			return omxGetState(os->childList[subState % os->numChildren], subState / os->numChildren);
		}
	}

	void omxUpdateState(omxState* tgt, omxState* src) {
		tgt->computeCount 		= src->computeCount;
		tgt->currentRow 		= src->currentRow;
		tgt->optimalValues 		= src->optimalValues;
	}
	
	void omxDuplicateState(omxState* tgt, omxState* src, unsigned short fullCopy) {
		tgt->numMats 			= src->numMats;
		tgt->numAlgs 			= src->numAlgs;
		tgt->numData 			= src->numData;
		tgt->dataList			= src->dataList;
		tgt->numChildren 		= 0;
		
		// Duplicate matrices and algebras and build parentLists.
		tgt->parentState 		= src;
		tgt->parentMatrix 		= src->matrixList;
		tgt->parentAlgebra 		= src->algebraList;
		tgt->matrixList			= (omxMatrix**) R_alloc(tgt->numMats, sizeof(omxMatrix*));
		for(int j = 0; j < tgt->numMats; j++) {
			// TODO: Smarter inference for which matrices to duplicate
			tgt->matrixList[j] = omxDuplicateMatrix(src->matrixList[j], tgt, fullCopy);
		}
				
		tgt->parentConList 		= src->conList;
		tgt->numConstraints     = src->numConstraints;
		tgt->conList			= (omxConstraint*) R_alloc(tgt->numConstraints, sizeof(omxConstraint));
		for(int j = 0; j < tgt->numConstraints; j++) {
			tgt->conList[j].size   = src->conList[j].size;
			tgt->conList[j].opCode = src->conList[j].opCode;
			tgt->conList[j].lbound = src->conList[j].lbound;
			tgt->conList[j].ubound = src->conList[j].ubound;
			tgt->conList[j].result = omxDuplicateMatrix(src->conList[j].result, tgt, fullCopy);
		}

		tgt->algebraList		= (omxMatrix**) R_alloc(tgt->numAlgs, sizeof(omxMatrix*));

		for(int j = 0; j < tgt->numAlgs; j++) {
			// TODO: Smarter inference for which algebras to duplicate
			tgt->algebraList[j] = omxDuplicateMatrix(src->algebraList[j], tgt, fullCopy);
		}

		for(int j = 0; j < tgt->numAlgs; j++) {
			omxDuplicateAlgebra(tgt->algebraList[j], src->algebraList[j], tgt, fullCopy);
		}

		
		tgt->childList 			= NULL;

		tgt->objectiveMatrix	= omxLookupDuplicateElement(tgt, src->objectiveMatrix);
		tgt->hessian 			= src->hessian;

		tgt->numFreeParams			= src->numFreeParams;
		tgt->freeVarList 		= (omxFreeVar*) R_alloc(tgt->numFreeParams, sizeof(omxFreeVar));
		for(int j = 0; j < tgt->numFreeParams; j++) {
			tgt->freeVarList[j].lbound			= src->freeVarList[j].lbound;
			tgt->freeVarList[j].ubound			= src->freeVarList[j].ubound;
			tgt->freeVarList[j].numLocations	= src->freeVarList[j].numLocations;
			
			int nLocs 							= tgt->freeVarList[j].numLocations;
			tgt->freeVarList[j].matrices		= (int*) R_alloc(nLocs, sizeof(int));
			tgt->freeVarList[j].row				= (int*) R_alloc(nLocs, sizeof(int));
			tgt->freeVarList[j].col				= (int*) R_alloc(nLocs, sizeof(int));

			for(int k = 0; k < nLocs; k++) {
				int theMat 						= src->freeVarList[j].matrices[k];
				int theRow 						= src->freeVarList[j].row[k];
				int theCol						= src->freeVarList[j].col[k];

				tgt->freeVarList[j].matrices[k] = theMat;
				tgt->freeVarList[j].row[k]		= theRow;
				tgt->freeVarList[j].col[k]		= theCol;
								
				tgt->freeVarList[j].name		= src->freeVarList[j].name;
			}
		}
		
		if (src->optimizerState) {
			tgt->optimizerState 					= (omxOptimizerState*) R_alloc(1, sizeof(omxOptimizerState));
			tgt->optimizerState->currentParameter	= src->optimizerState->currentParameter;
			tgt->optimizerState->offset				= src->optimizerState->offset;
			tgt->optimizerState->alpha				= src->optimizerState->alpha;
		}
		
		tgt->optimalValues 		= src->optimalValues;
		tgt->optimum 			= 9999999999;
                                  
		tgt->majorIteration 	= 0;
		tgt->minorIteration 	= 0;
		tgt->startTime 			= src->startTime;
		tgt->endTime			= 0;
		
		// TODO: adjust checkpointing based on parallelization method
		tgt->numCheckpoints     = 0;
		tgt->checkpointList 	= NULL;
		tgt->chkptText1 		= NULL;
		tgt->chkptText2 		= NULL;
                                  
		tgt->computeCount 		= src->computeCount;
		tgt->currentRow 		= src->currentRow;

		tgt->statusCode 		= 0;
		strncpy(tgt->statusMsg, "", 1);
	}

    omxMatrix* omxLookupDuplicateElement(omxState* os, omxMatrix* element) {
        if(os == NULL || element == NULL) return NULL;

		if (element->hasMatrixNumber) {
			int matrixNumber = element->matrixNumber;
			if (matrixNumber >= 0) {
				return(os->algebraList[matrixNumber]);
			} else {
				return(os->matrixList[-matrixNumber - 1]);
			}
		}

        for(int i = 0; i < os->numConstraints; i++) {
            if(os->parentConList[i].result == element) {
				if(os->conList[i].result != NULL)   // Not sure of proper failure behavior here.
                    return(os->conList[i].result);
                else
                    omxRaiseError(os, -2, "Initialization Copy Error: Constraint required but not yet processed.");
            }
        }

        return NULL;
    }

	void omxFreeState(omxState *oo) {
		int k;

		if (oo->numChildren > 0) {
			for(k = 0; k < oo->numChildren; k++) {
				omxFreeState(oo->childList[k]);
			}
			Free(oo->childList);
			oo->childList = NULL;
			oo->numChildren = 0;
		}

		if(OMX_DEBUG) { Rprintf("Freeing %d Algebras.\n", oo->numAlgs);}
		for(k = 0; k < oo->numAlgs; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Algebra %d at 0x%x.\n", k, oo->algebraList[k]); }
			omxFreeAllMatrixData(oo->algebraList[k]);
		}

		if(OMX_DEBUG) { Rprintf("Freeing %d Matrices.\n", oo->numMats);}
		for(k = 0; k < oo->numMats; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Matrix %d at 0x%x.\n", k, oo->matrixList[k]); }
			omxFreeAllMatrixData(oo->matrixList[k]);
		}

		if(OMX_DEBUG) { Rprintf("Freeing %d Data Sets.\n", oo->numData);}
		for(k = 0; k < oo->numData; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Data Set %d at 0x%x.\n", k, oo->dataList[k]); }
			omxFreeData(oo->dataList[k]);
		}

        if(OMX_DEBUG) {Rprintf("Freeing %d Children.\n", oo->numChildren);}
        for(k = 0; k < oo->numChildren; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Child State %d at 0x%x.\n", k, oo->childList[k]); }
			omxFreeState(oo->childList[k]);            
        }

		if(OMX_DEBUG) { Rprintf("Freeing %d Checkpoints.\n", oo->numCheckpoints);}
		for(k = 0; k < oo->numCheckpoints; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Data Set %d at 0x%x.\n", k, oo->checkpointList[k]); }
			omxCheckpoint oC = oo->checkpointList[k];
			switch(oC.type) {
				case OMX_FILE_CHECKPOINT:
					fclose(oC.file);
					break;
				case OMX_SOCKET_CHECKPOINT:	// NYI :::DEBUG:::
					// TODO: Close socket
					break;
				case OMX_CONNECTION_CHECKPOINT:	// NYI :::DEBUG:::
					// Do nothing: this should be handled by R upon return.
					break;
			}
			if(oo->chkptText1 != NULL) {
				Free(oo->chkptText1);
			}
			if(oo->chkptText2 != NULL) {
				Free(oo->chkptText2);
			}
			// Checkpoint list itself is freed by R.
		}

		if(OMX_DEBUG) { Rprintf("State Freed.\n");}
	}

	void omxSaveState(omxState *os, double* freeVals, double minimum) {
		if(os->optimalValues == NULL) {
			os->optimalValues = (double*) R_alloc(os->numFreeParams, sizeof(double));
		}

		for(int i = 0; i < os->numFreeParams; i++) {
			os->optimalValues[i] = freeVals[i];
		}
		os->optimum = minimum;
		os->optimumStatus = os->statusCode;
		strncpy(os->optimumMsg, os->statusMsg, 250);
	}

	void omxRaiseError(omxState *oo, int errorCode, char* errorMsg) {
		if(OMX_DEBUG && errorCode) { Rprintf("Error %d raised: %s", errorCode, errorMsg);}
		if(OMX_DEBUG && !errorCode) { Rprintf("Error status cleared."); }
		oo->statusCode = errorCode;
		strncpy(oo->statusMsg, errorMsg, 249);
		oo->statusMsg[249] = '\0';
		if(oo->computeCount <= 0 && errorCode < 0) {
			oo->statusCode--;			// Decrement status for init errors.
		}
	}

	void omxStateNextRow(omxState *oo) {
		oo->currentRow++;
	};
	void omxStateNextEvaluation(omxState *oo) {
		oo->currentRow = 0;
		oo->computeCount++;
	};

	void omxSaveCheckpoint(omxState *os, double* x, double* f) {
		time_t now = time(NULL);
		int soFar = now - os->startTime;		// Translated into minutes
		int n;
		for(int i = 0; i < os->numCheckpoints; i++) {
			n = 0;
			omxCheckpoint* oC = &(os->checkpointList[i]);
			// Check based on time
			if(oC->time > 0 && (soFar - oC->lastCheckpoint) >= oC->time) {
				oC->lastCheckpoint = soFar;
				n = 1;
			}
			// Or iterations
			if(oC->numIterations > 0 && (os->majorIteration - oC->lastCheckpoint) >= oC->numIterations) {
				oC->lastCheckpoint = os->majorIteration;
				n = 1;
			}

			if(n) {		//In either case, save a checkpoint.
				if(os->chkptText1 == NULL) {	// First one: set up output
					// FIXME: Is it faster to allocate this on the stack?
					os->chkptText1 = (char*) Calloc((24+15*os->numFreeParams), char);
					os->chkptText2 = (char*) Calloc(1.0+15.0*os->numFreeParams*
														(os->numFreeParams + 1.0)/2.0, char);
					if (oC->type == OMX_FILE_CHECKPOINT) {
						fprintf(oC->file, "iterations\ttimestamp\tobjective\t");
						for(int j = 0; j < os->numFreeParams; j++) {
							if(strcmp(os->freeVarList[j].name, CHAR(NA_STRING)) == 0) {
								fprintf(oC->file, "%s", os->freeVarList[j].name);
							} else {
								fprintf(oC->file, "\"%s\"", os->freeVarList[j].name);
							}
							if (j != os->numFreeParams - 1) fprintf(oC->file, "\t");
						}
						fprintf(oC->file, "\n");
						fflush(oC->file);
					}
				}
				char tempstring[25];
				sprintf(tempstring, "%d", os->majorIteration);

				if(strncmp(os->chkptText1, tempstring, strlen(tempstring))) {	// Returns zero if they're the same.
					struct tm * nowTime = localtime(&now);						// So this only happens if the text is out of date.
					strftime(tempstring, 25, "%b %d %Y %I:%M:%S %p", nowTime);
					sprintf(os->chkptText1, "%d \"%s\" %9.5f", os->majorIteration, tempstring, f[0]);
					for(int j = 0; j < os->numFreeParams; j++) {
						sprintf(tempstring, " %9.5f", x[j]);
						strncat(os->chkptText1, tempstring, 14);
					}

					double* hessian = os->hessian;
					if(hessian != NULL) {
						for(int j = 0; j < os->numFreeParams; j++) {
							for(int k = 0; k <= j; k++) {
								sprintf(tempstring, " %9.5f", hessian[j]);
								strncat(os->chkptText2, tempstring, 14);
							}
						}
					}
				}

				if(oC->type == OMX_FILE_CHECKPOINT) {
					fprintf(oC->file, "%s", os->chkptText1);
					if(oC->saveHessian)
						fprintf(oC->file, "%s", os->chkptText2);
					fprintf(oC->file, "\n");
					fflush(oC->file);
				} else if(oC->type == OMX_SOCKET_CHECKPOINT) {
					n = write(oC->socket, os->chkptText1, strlen(os->chkptText1));
					if(n != strlen(os->chkptText1)) warning("Error writing checkpoint.");
					if(oC->saveHessian) {
						n = write(oC->socket, os->chkptText2, strlen(os->chkptText2));
						if(n != strlen(os->chkptText1)) warning("Error writing checkpoint.");
					}
					n = write(oC->socket, "\n", 1);
					if(n != 1) warning("Error writing checkpoint.");
				} else if(oC->type == OMX_CONNECTION_CHECKPOINT) {
					warning("NYI: R_connections are not yet implemented.");
					oC->numIterations = 0;
					oC->time = 0;
				}
			}
		}
	}
