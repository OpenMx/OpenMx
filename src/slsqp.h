#ifndef SLSQP_H
#define SLSQP_H

#include "nlopt.h"
#include "nlopt-util.h"

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

	typedef struct nlopt_slsqp_wdump{
		double *realwkspc;
		int lengths[8];
		int M;
	} nlopt_slsqp_wdump;
	
	nlopt_result nlopt_slsqp(unsigned n, nlopt_func f, void *f_data,
                          unsigned m, nlopt_constraint *fc,
                          unsigned p, nlopt_constraint *h,
                          const double *lb, const double *ub,
                          double *x, double *minf, nlopt_slsqp_wdump *wkspc,
                          nlopt_stopping *stop);
#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif
