#ifndef _RAMINTERNAL_H_
#define _RAMINTERNAL_H_

struct join {
	int foreignKey;
	struct omxExpectation *ex;
	omxMatrix *regression;
};

struct omxRAMExpectation {

	omxMatrix *cov, *means; // observed covariance and means
	omxMatrix *A, *S, *F, *M, *I;
	omxMatrix *X, *Y, *Z, *Ax;

	int numIters;
	double logDetObserved;
	double n;
	double *work;
	int lwork;

	std::vector<join> joins;

	void ensureTrivialF();
};

#endif
