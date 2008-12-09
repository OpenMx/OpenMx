#ifndef _OMX_ALGEBRA_H_
#define _OMX_ALGEBRA_H_ TRUE

#include <R.h> 
#include <Rinternals.h> 
#include <Rdefines.h>
#include <R_ext/Rdynload.h> 
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "omxDataMatrix.h"
#include "omxSymbolTable.h"
#include "omxAlgebraFunctions.h"

class omxAlgebra: public omxDataMatrix 
{

protected:
	void* funWrapper;
	omxDataMatrix* args;
	int numArgs;
	bool isDirty;
	
public:
	/* Constructor */
	omxAlgebra();
	omxAlgebra(SEXP alg);
	/* Overloads of dataMatrix functions */
	void compute();
	void recompute();
	bool needsUpdate();
	void fillFromMxAlgebra(SEXP alg);
};

#endif