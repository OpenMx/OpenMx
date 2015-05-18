#ifndef _OMXBLAS_H
#define _OMXBLAS_H

#include <R_ext/BLAS.h>

#ifdef  __cplusplus
extern "C" {
#endif

void F77_CALL(omxunsafedgemm)(int*, int*, int*, int*,
	int*, double*, double*, int*, double*, int*, double*,
	double*, int*);

void F77_CALL(omxunsafedgemv)(int*, int*, int*, double*, double*,
	int*, double*, int*, double*, double*, int*);

#ifdef  __cplusplus
}
#endif

#endif //_OMXBLAS_H
