 /*
 *  Copyright 2007-2015 The OpenMx Project
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
  SEXP dV, dVnames;
  int i=0;
  omxState* currentState = ox->currentState;
  
  if(OMX_DEBUG) { mxLog("Initializing GREML expectation."); }
  
  //omxGREMLExpectation *oge = (omxGREMLExpectation*) R_alloc(1, sizeof(omxGREMLExpectation));
  omxGREMLExpectation *oge = new omxGREMLExpectation;
  
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
  
  {ScopedProtect p1(dV, R_do_slot(rObj, Rf_install("dV")));
	ScopedProtect p2(dVnames, R_do_slot(rObj, Rf_install("dVnames")));
  oge->dVlength = Rf_length(dV);  
  oge->dV.resize(oge->dVlength);
  oge->dVnames.resize(oge->dVlength);
	if(oge->dVlength){
    if(OMX_DEBUG) { mxLog("Processing derivatives of V."); }
		int* dVint = INTEGER(dV);
    for(i=0; i < oge->dVlength; i++){
      oge->dV[i] = omxMatrixLookupFromState1(dVint[i], currentState);
      SEXP elem;
      {ScopedProtect p3(elem, STRING_ELT(dVnames, i));
			oge->dVnames[i] = CHAR(elem);}
	}}
  }
}


void omxComputeGREMLExpectation(omxExpectation* ox, const char *, const char *) {
  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
	omxRecompute(oge->V, NULL);
  omxRecompute(oge->y, NULL);
}


void omxDestroyGREMLExpectation(omxExpectation* ox) {
	if(OMX_DEBUG) { mxLog("Destroying GREML Expectation."); } //Is this all??
}


void omxPopulateGREMLAttributes(omxExpectation *ox, SEXP algebra) {
  if(OMX_DEBUG) { mxLog("Populating GREML expectation attributes."); }

  omxGREMLExpectation* oge = (omxGREMLExpectation*) (ox->argStruct);
  
  Rf_setAttrib(algebra, Rf_install("numStats"), Rf_ScalarReal(oge->y->rows));
  Rf_setAttrib(algebra, Rf_install("numFixEff"), Rf_ScalarInteger(oge->X->cols));
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


omxMatrix* omxMatrixLookupFromState1(int matrix, omxState* os) {
	omxMatrix* output = NULL;
	if(matrix == NA_INTEGER){return NULL;}
	if (matrix >= 0) {
		output = os->algebraList[matrix];
	} 
  else {
		output = os->matrixList[~matrix];
	}
	return output;
}