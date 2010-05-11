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
	void omxInitState(omxState* state) {
		state->numMats = 0;
		state->numAlgs = 0;
		state->numData = 0;
		state->matrixList = NULL;	
		state->algebraList = NULL;
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
	
	void omxFreeState(omxState *oo) {
		int k;
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
		oo->statusCode = errorCode;
		strncpy(oo->statusMsg, errorMsg, 249);
		oo->statusMsg[249] = '\0';
	}
	
	void omxStateNextRow(omxState *oo) { 
		oo->currentRow++; 
	};
	void omxStateNextEvaluation(omxState *oo) { 
		oo->currentRow = 0;
		oo->computeCount++; 
	};
	
	void omxSaveCheckpoint(omxState *os, double* x) {
		time_t now = time(NULL);
		int soFar = (now - os->startTime) / 60;		// Translated into minutes
		int n;
		for(int i = 0; i < os->numCheckpoints; i++) {
			n = 0;
			omxCheckpoint oC = os->checkpointList[i];
			// Check based on time
			if(oC.time != 0 && soFar - oC.lastCheckpoint >= oC.time) {
				oC.lastCheckpoint = now;
				n = 1;
			}
			// Or iterations
			if(oC.numIterations != 0 && os->majorIteration - oC.lastCheckpoint >= oC.numIterations) {
				oC.lastCheckpoint = os->majorIteration;
				n = 1;
			}

			if(n) {		//In either case, save a checkpoint.
				if(os->chkptText1 == NULL) {	// First one: set up output
					// FIXME: Is it faster to allocate this on the stack?
					os->chkptText1 = (char*) Calloc((24+15*os->numFreeParams), char);
					os->chkptText2 = (char*) Calloc(1.0+15.0*os->numFreeParams*
														(os->numFreeParams + 1.0)/2.0, char);
				}

				char tempstring[25];
				sprintf(tempstring, "%d", os->majorIteration);

				if(strncmp(os->chkptText1, tempstring, strlen(tempstring))) {// Returns zero if they're the same.
					struct tm * nowTime = localtime(&now);
					strftime(tempstring, 25, "%c", nowTime);
					sprintf(os->chkptText1, "%d %s", os->majorIteration, tempstring);
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
				
				if(oC.type == OMX_FILE_CHECKPOINT) {
					fprintf(oC.file, "%s", os->chkptText1);
					if(oC.saveHessian) 
						fprintf(oC.file, "%s", os->chkptText2);
					fprintf(oC.file, "\n");
					fflush(oC.file);
				} else if(oC.type == OMX_SOCKET_CHECKPOINT) {
					n = write(oC.socket, os->chkptText1, strlen(os->chkptText1));
					if(n != strlen(os->chkptText1)) warning("Error writing checkpoint.");
					if(oC.saveHessian) {
						n = write(oC.socket, os->chkptText2, strlen(os->chkptText2));
						if(n != strlen(os->chkptText1)) warning("Error writing checkpoint.");
					}
					n = write(oC.socket, "\n", 1);
					if(n != 1) warning("Error writing checkpoint.");
				} else if(oC.type == OMX_CONNECTION_CHECKPOINT) {
					warning("NYI: R_connections are not yet implemented.");
					oC.numIterations = 9999999;
					oC.time = 99999999999;
				}
			}
		}
	}
