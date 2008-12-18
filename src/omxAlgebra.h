#ifndef _OMX_ALGEBRA_H_
#define _OMX_ALGEBRA_H_ TRUE

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxMatrix.h"
#include "omxSymbolTable.h"
#include "omxAlgebraFunctions.h"

class omxAlgebra: public omxMatrix 
{

protected:
	void* funWrapper;
	omxMatrix** args;
	int numArgs;
	
	omxMatrix* MxAlgebraParseHelper(int* &spec);
	
public:
	/* Constructor */
	omxAlgebra();
//	omxAlgebra(SEXP alg);
	/* Overloads of dataMatrix functions */
	void compute();
	void recompute();
	bool needsUpdate();
	void fillFromMxAlgebra(SEXP alg);
	void fillFromTableEntry(const omxAlgebraTableEntry* oate);
	
public:
	/* Factory Method */
//	static omxMatrix* fillFromMxAlgebra(SEXP alg);
};

omxMatrix* omxMatrixFromMxMatrixPtr(SEXP s);

#endif