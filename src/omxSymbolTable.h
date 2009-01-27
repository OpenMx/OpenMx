#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxAlgebraFunctions.h"

#ifndef _OMX_SYMBOL_TABLE_
#define _OMX_SYMBOL_TABLE_ TRUE

typedef struct {

	unsigned int number;
	const char opName[250];
	const char rName[250];
	short int numArgs;
	void* funWrapper;

} omxAlgebraTableEntry;

const omxAlgebraTableEntry omxAlgebraSymbolTable[14];


#endif