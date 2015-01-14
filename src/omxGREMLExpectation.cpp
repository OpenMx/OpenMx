 /*
 *  Copyright 2007-2014 The OpenMx Project
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
 
#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxDefines.h"
#include "omxGREMLExpectation.h"
 
void omxInitGREMLExpectation(omxExpectation* ox){
  
  SEXP rObj = ox->rObj;
  omxState* currentState = ox->currentState;
  
  if(OMX_DEBUG) { mxLog("Initializing GREML expectation."); }
  
  omxGREMLExpectation *oge = (omxGREMLExpectation*) R_alloc(1, sizeof(omxGREMLExpectation));
  
  /* Set Expectation Calls and Structures */
  ox->computeFun = omxComputeGREMLExpectation;
	ox->destructFun = omxDestroyGREMLExpectation;
	ox->componentFun = omxGetGREMLExpectationComponent;
	ox->populateAttrFun = omxPopulateGREMLAttributes;
	ox->argStruct = (void*) oge;
  
    /* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Processing V."); }
	oge->V = omxNewMatrixFromSlot(rObj, currentState, "V");

	if(OMX_DEBUG) { mxLog("Processing X."); }
	oge->X = omxNewMatrixFromSlot(rObj, currentState, "X");

	oge->y = omxNewMatrixFromSlot(rObj, currentState, "y");
  
}


void omxComputeGREMLExpectation(omxExpectation* ox, const char *, const char *) {
  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
	omxRecompute(oge->V, NULL);
}


void omxDestroyGREMLExpectation(omxExpectation* ox) {
	if(OMX_DEBUG) { mxLog("Destroying GREML Expectation."); } //Is this all??
}


void omxPopulateGREMLAttributes(omxExpectation *ox, SEXP algebra) {
  if(OMX_DEBUG) { mxLog("Populating GREML expectation attributes."); }

  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
    
	omxMatrix* V = oge->V;
	omxRecompute(V, NULL);
}


omxMatrix* omxGetGREMLExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component){
/* Return appropriate parts of Expectation to the Fit Function */
  if(OMX_DEBUG) { mxLog("GREML expectation: %s requested--", component); }

	omxGREMLExpectation* oge = (omxGREMLExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(strEQ("V", component)) {
		retval = oge->V;
	} else if(strEQ("X", component)) {
		retval = oge->X;
	} else if(strEQ("y", component)) {
		retval = oge->y;
	}
	if (retval) omxRecompute(retval, NULL);
	
	return retval;
}
