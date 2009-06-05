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
		state->matrixList = NULL;	
		state->algebraList = NULL;
		state->objective = NULL;
		state->conList = NULL;
		state->freeVarList = NULL;
		state->optimizerState;			

		state->computeCount = 0;
		state->currentRow = 0;

		state->statusCode = 0;
		strncpy(state->statusMsg, "", 1);
	}

	void omxFillState(omxState* state, /*omxOptimizer *oo,*/ omxMatrix** matrixList, omxMatrix** algebraList, omxMatrix* objective) {
		error("NYI: Can't fill a state from outside yet.  Besides, do you really need a single function to do this?");
	}
	
	void omxFreeState(omxState *oo) {
		int k;
		if(OMX_DEBUG) { Rprintf("Freeing Algebras.\n");}
		for(k = 0; k < oo->numAlgs; k++) {
			if(OMX_DEBUG) { Rprintf("Freeing Algebra %d.\n", k); }
			omxFreeAllMatrixData(oo->algebraList[k]);
		}

		if(OMX_DEBUG) { Rprintf("Freeing Matrices.\n");}
		for(k = 0; k < oo->numMats; k++) {
			omxFreeAllMatrixData(oo->matrixList[k]);
		}
		
	}
	
	void omxStateNextRow(omxState *oo) { 
		oo->currentRow++; 
	};
	void omxStateNextEvaluation(omxState *oo) { 
		oo->computeCount++; 
	};
	