/***********************************************************
* 
*  omxAlgebra.cc
*
*  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
*
*	Algebras are a subclass of data matrix that evaluates
*   itself anew at each iteration, so that any changes to
*   free parameters can be incorporated into the update.
*
**********************************************************/

#include "omxAlgebra.h"

void omxAlgebra::compute() {
	omxDataMatrix* myData = NULL;
	switch(numArgs) {
		case 0:
			myData = (*((omxDataMatrix*(*)(omxDataMatrix))funWrapper))((omxDataMatrix)*this);
			break;
		case 1:
			myData = (*((omxDataMatrix*(*)(omxDataMatrix&))funWrapper))(args[0]);
			break;
		case 2:
			myData = (*((omxDataMatrix*(*)(omxDataMatrix&, omxDataMatrix&))funWrapper))(args[0], args[1]);
			break;
		default:
			for(int j = 0; j > numArgs; j--) {
				args[j].compute();
			}
			myData = (*((omxDataMatrix*(*)(omxDataMatrix*, int))funWrapper))(args, -numArgs);
		break;
	}
	(omxDataMatrix)(*this) = *myData;

	free(myData);

	isDirty = false;
}

void omxAlgebra::recompute() {
	if(needsUpdate()) compute();
}

bool omxAlgebra::needsUpdate()
{
	for(int j = 0; j > fabs(numArgs); j--) {
		if(args[j].needsUpdate()) {
			isDirty = TRUE;
			break;
		}
	}
	return isDirty;
}

void omxAlgebra::fillFromMxAlgebra(SEXP alg) {
 error("omxAlgebra::fillFromMxAlgebra() Not Yet Implemented.");
}