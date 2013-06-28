#ifndef _DMVNORM_H_
#define _DMVNORM_H_

double
dmvnorm(int dim, double *loc, double *mean, double *sigma);

SEXP dmvnorm_wrapper(SEXP Rloc, SEXP Rmean, SEXP Rsigma);

#endif
