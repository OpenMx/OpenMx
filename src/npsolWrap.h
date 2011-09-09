#ifndef _NPSOLWRAP_H
#define _NPSOLWRAP_H

#include "omxState.h"

/* Functions for Export */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state);  // Calls NPSOL.  Duh.

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options);
SEXP findIdenticalRowsData(SEXP data, SEXP missing, SEXP defvars,
	SEXP skipMissingness, SEXP skipDefvars);

extern omxState* currentState;

#endif // #define _NPSOLWRAP_H
