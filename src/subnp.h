#ifndef _SUBNP_H_
#define _SUBNP_H_

#include <math.h>
#include <float.h>               
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "matrix.h"

Matrix solnp(Matrix solPars, double (*solFun)(Matrix),  Matrix solEqB, Matrix (*solEqBFun)(Matrix), Matrix (*myineqFun)( Matrix), Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,  Matrix solctrl, bool debugToggle);

int subnp(Matrix pars, double (*solFun)( Matrix), Matrix (*solEqBFun)(Matrix) , Matrix (*myineqFun)( Matrix),  Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,  Matrix ctrl);

#endif