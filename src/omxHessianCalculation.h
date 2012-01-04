/*
 *  Copyright 2007-2012 The OpenMx Project
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

/***********************************************************
* 
*  omxMatrix.h
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Contains header information for the omxMatrix class
*   omxDataMatrices hold necessary information to simplify
* 	dealings between the OpenMX back end and BLAS.
*
**********************************************************/

#ifndef _OMXHESSIAN_CALCULATION_H
#define _OMXHESSIAN_CALCULATION_H

unsigned short omxEstimateHessian(int numHessians, double functionPrecision, 
	int r, omxState* currentState);


#endif // _OMXHESSIAN_CALCULATION_H
