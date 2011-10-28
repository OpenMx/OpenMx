/*
 *  Copyright 2007-2011 The OpenMx Project
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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDefines.h"
#include "omxAlgebraFunctions.h"
#include "omxSymbolTable.h"
#include "omxData.h"
#include "omxFIMLObjective.h"
#include "omxOpenmpWrap.h"

#ifdef _OPENMP

omp_lock_t sadmvn_lock;

#else

void* sadmvn_lock = NULL;

#endif

extern void F77_SUB(sadmvn)(int*, double*, double*, int*, double*, int*, double*, double*, double*, double*, int*);


void omxSadmvnWrapper(omxObjective *oo, omxMatrix *cov, omxMatrix *ordCov, 
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
   	//	Abseps	double		Absolute error tolerance.  Yick.
   	//	Releps	double		Relative error tolerance.  Use EPSILON.
   	//	Error	&double		On return: absolute real error, 99% confidence
   	//	Value	&double		On return: evaluated value
   	//	Inform	&int		On return: 0 = OK; 1 = Rerun, increase MaxPts; 2 = Bad input
   	// TODO: Separate block diagonal covariance matrices into pieces for integration separately
   	double Error;
   	double absEps = 1e-3;
   	double relEps = 0;
   	int MaxPts = 100000*cov->rows;
   	int numVars = ordCov->rows;
   	/* FOR DEBUGGING PURPOSES */
    /*	numVars = 2;
   	lThresh[0] = -2;
   	uThresh[0] = -1.636364;
   	Infin[0] = 2;
   	lThresh[1] = 0;
   	uThresh[1] = 0;
   	Infin[1] = 0;
   	smallCor[0] = 1.0; smallCor[1] = 0; smallCor[2] = 1.0; */
	omx_omp_set_lock(&sadmvn_lock);
   	F77_CALL(sadmvn)(&numVars, lThresh, uThresh, Infin, corList, &MaxPts, &absEps, &relEps, &Error, likelihood, inform);
	omx_omp_unset_lock(&sadmvn_lock);

   	if(OMX_DEBUG && !oo->matrix->currentState->currentRow) {
   		char infinCodes[3][20];
   		strcpy(infinCodes[0], "(-INF, upper]");
   		strcpy(infinCodes[1], "[lower, INF)");
   		strcpy(infinCodes[2], "[lower, upper]");
   		Rprintf("Input to sadmvn is (%d rows):\n", numVars); //:::DEBUG:::
		omxPrint(ordCov, "Ordinal Covariance Matrix"); //:::DEBUG:::
		for(int i = 0; i < numVars; i++) {
			Rprintf("Row %d: %f, %f, %d(%s)\n", i, lThresh[i], uThresh[i], Infin[i], infinCodes[Infin[i]]);
		}

		Rprintf("Cor: (Lower %d x %d):", cov->rows, cov->cols); //:::DEBUG:::
		for(int i = 0; i < cov->rows*(cov->rows-1)/2; i++) {
			// Rprintf("Row %d of Cor: ", i);
			// for(int j = 0; j < i; j++)
			Rprintf(" %f", corList[i]); // (i*(i-1)/2) + j]);
			// Rprintf("\n");
		}
		Rprintf("\n");
	}

	if(OMX_DEBUG) {
		Rprintf("Output of sadmvn is %f, %f, %d.\n", Error, *likelihood, *inform); 
	}
} 
