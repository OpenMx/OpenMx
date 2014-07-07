/*
 *  Copyright 2007-2014 The OpenMx Project
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

#include "omxDefines.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLFitFunction.h"
#include "omxOpenmpWrap.h"
#include "omxSadmvnWrapper.h"

void omxSadmvnWrapper(omxFitFunction *oo, omxMatrix *cov, omxMatrix *ordCov, 
	double *corList, double *lThresh, double *uThresh, int *Infin, double *likelihood, int *inform) {
    // SADMVN calls Alan Genz's sadmvn.f--see appropriate file for licensing info.
   	// TODO: Check with Genz: should we be using sadmvn or sadmvn?
   	// Parameters are:
   	// 	N 		int			# of vars
   	//	Lower	double*		Array of lower bounds
   	//	Upper	double*		Array of upper bounds
   	//	Infin	int*		Array of flags: 0 = (-Inf, upper] 1 = [lower, Inf), 2 = [lower, upper]
   	//	Correl	double*		Array of correlation coeffs: in row-major lower triangular order
   	//	MaxPts	int			Maximum # of function values (use 1000*N or 1000*N*N)
   	//	Abseps	double		Absolute Rf_error tolerance.  Yick.
   	//	Releps	double		Relative Rf_error tolerance.  Use EPSILON.
   	//	Error	&double		On return: absolute real Rf_error, 99% confidence
   	//	Value	&double		On return: evaluated value
   	//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
   	// TODO: Separate block diagonal covariance matrices into pieces for integration separately
   	double Error;
	double absEps = Global->absEps;
	double relEps = Global->relEps;
	int MaxPts = Global->maxptsa + Global->maxptsb * cov->rows + Global->maxptsc * cov->rows * cov->rows;
	int numVars = ordCov->rows;
	int fortranThreadId = omx_absolute_thread_num() + 1;
   	/* FOR DEBUGGING PURPOSES */
    /*	numVars = 2;
   	lThresh[0] = -2;
   	uThresh[0] = -1.636364;
   	Infin[0] = 2;
   	lThresh[1] = 0;
   	uThresh[1] = 0;
   	Infin[1] = 0;
   	smallCor[0] = 1.0; smallCor[1] = 0; smallCor[2] = 1.0; */
   	F77_CALL(sadmvn)(&numVars, lThresh, uThresh, Infin, corList, &MaxPts, 
		&absEps, &relEps, &Error, likelihood, inform, &fortranThreadId);

   	if(OMX_DEBUG && !oo->matrix->currentState->currentRow) {
   		char infinCodes[3][20];
   		strcpy(infinCodes[0], "(-INF, upper]");
   		strcpy(infinCodes[1], "[lower, INF)");
   		strcpy(infinCodes[2], "[lower, upper]");
   		mxLog("Input to sadmvn is (%d rows):", numVars); //:::DEBUG:::
		omxPrint(ordCov, "Ordinal Covariance Matrix"); //:::DEBUG:::
		for(int i = 0; i < numVars; i++) {
			mxLog("Row %d: %f, %f, %d(%s)", i, lThresh[i], uThresh[i], Infin[i], infinCodes[Infin[i]]);
		}

		mxLog("Cor: (Lower %d x %d):", cov->rows, cov->cols); //:::DEBUG:::
		for(int i = 0; i < cov->rows*(cov->rows-1)/2; i++) {
			// mxLog("Row %d of Cor: ", i);
			// for(int j = 0; j < i; j++)
			mxLog(" %f", corList[i]); // (i*(i-1)/2) + j]);
			// mxLog("");
		}
	}

	if(OMX_DEBUG) {
		mxLog("Output of sadmvn is %f, %f, %d.", Error, *likelihood, *inform); 
	}
} 
