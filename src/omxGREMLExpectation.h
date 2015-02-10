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
 
 typedef struct {
  omxMatrix* y;
  omxMatrix* V;
  omxMatrix* X;
  std::vector< omxMatrix* > dV;
  std::vector< const char* > dVnames;
  int dVlength;
  int do_fixeff;
} omxGREMLExpectation;

void omxInitGREMLExpectation(omxExpectation* ox);
void omxComputeGREMLExpectation(omxExpectation* ox, const char *, const char *);
void omxDestroyGREMLExpectation(omxExpectation* ox);
void omxPopulateGREMLAttributes(omxExpectation *ox, SEXP algebra);
omxMatrix* omxGetGREMLExpectationComponent(omxExpectation* ox, omxFitFunction* off, const char* component);
omxMatrix* omxMatrixLookupFromState1(int matrix, omxState* os);