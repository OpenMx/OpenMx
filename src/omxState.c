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
		state->conList = NULL;
		state->freeVarList = NULL;
		state->optimizerState = NULL;

		state->computeCount = -1;
		state->currentRow = -1;

		state->statusCode = 0;
		strncpy(state->statusMsg, "", 1);
	}

	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList, omxMatrix** algebraList, omxData** dataList, omxMatrix* objective) {
		error("NYI: Can't fill a state from outside yet.  Besides, do you really need a single function to do this?");
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
		
		if(OMX_DEBUG) { Rprintf("State Freed.\n");}
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
	
