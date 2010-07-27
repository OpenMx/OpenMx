#ifndef _NPSOLWRAP_H
#define _NPSOLWRAP_H

/* Functions for Export */
SEXP callNPSOL(SEXP objective, SEXP startVals, SEXP constraints,
	SEXP matList, SEXP varList, SEXP algList,
	SEXP data, SEXP intervalList, SEXP checkpointList, SEXP options, SEXP state);  // Calls NPSOL.  Duh.

SEXP omxCallAlgebra(SEXP matList, SEXP algNum, SEXP options);
SEXP findIdenticalRowsData(SEXP data);

/* Set up R .Call info */
R_CallMethodDef callMethods[] = {
{"callNPSOL", (void*(*)())&callNPSOL, 10},
{"omxCallAlgebra", (void*(*)())&omxCallAlgebra, 3},
{"findIdenticalRowsData", (void*(*)())&findIdenticalRowsData, 1},
{NULL, NULL, 0}
};

void
R_init_mylib(DllInfo *info)
{
/* Register routines, allocate resources. */
R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

void R_unload_mylib(DllInfo *info)
{
/* Release resources. */
}

#endif // #define _NPSOLWRAP_H
