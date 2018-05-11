/*
 *  Copyright 2007-2018 by the individuals mentioned in the source code history
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

#ifndef _OMX_R_FITFUNCTION_
#define _OMX_R_FITFUNCTION_

struct omxRFitFunction : omxFitFunction {

	SEXP fitfun;
	SEXP model;
	SEXP flatModel;
	SEXP state;

	virtual void init();
	virtual void compute(int ffcompute, FitContext *fc);
};

void omxInitRFitFunction(omxFitFunction* oo);

#endif /* _OMX_R_FITFUNCTION_ */
