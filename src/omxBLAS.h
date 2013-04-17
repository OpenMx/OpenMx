#ifndef _OMXBLAS_H
#define _OMXBLAS_H

#ifdef  __cplusplus
extern "C" {
#endif

void F77_SUB(omxunsafedgemm)(const char*, const char*, int*, int*,
	int*, double*, double*, int*, double*, int*, double*,
	double*, int*);

void F77_SUB(omxunsafedgemv)(const char*, int*, int*, double*, double*,
	int*, double*, int*, double*, double*, int*);

#ifdef  __cplusplus
}
#endif

#endif //_OMXBLAS_H
