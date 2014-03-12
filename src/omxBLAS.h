#ifndef _OMXBLAS_H
#define _OMXBLAS_H

#include <R_ext/BLAS.h>

#ifdef  __cplusplus
extern "C" {
#endif

void F77_CALL(omxunsafedgemm)(const char*, const char*, int*, int*,
	int*, double*, double*, int*, double*, int*, double*,
	double*, int*);

void F77_CALL(omxunsafedgemv)(const char*, int*, int*, double*, double*,
	int*, double*, int*, double*, double*, int*);

#ifdef  __cplusplus
}
#endif

#endif //_OMXBLAS_H
