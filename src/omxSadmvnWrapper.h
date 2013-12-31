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
 *
 */
 
#ifndef _OMXSADMVNWRAPPER_H
#define _OMXSADMVNWRAPPER_H

#ifdef  __cplusplus
extern "C" {
#endif

void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*,
		     double*, double*, double*, double*, int*, int*);

#ifdef  __cplusplus
}
#endif

void omxSadmvnWrapper(omxFitFunction *oo, omxMatrix *cov, omxMatrix *ordCov, 
	double *corList, double *lThresh, double *uThresh, int *Infin, double *likelihood, int *inform);

#endif 
