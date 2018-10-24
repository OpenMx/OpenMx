#ifndef _NLOPTCPP_H_
#define _NLOPTCPP_H_

void omxInvokeNLOPT(GradientOptimizerContext &goc);

struct UnconstrainedObjective {
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;

	virtual double *getParamVec()=0;
	virtual double getFit(const double *)=0;
	virtual void getGrad(double *)=0;
};

class UnconstrainedSLSQPOptimizer {
	const char *name;
	int maxIter;
	double tolerance;
	int verbose;
	int iter;
	UnconstrainedObjective *uo;

	double evaluate(const double *x, double *grad);
	static double obj(unsigned n, const double *x, double *grad, void *mythis);
public:
	UnconstrainedSLSQPOptimizer(const char *_name, int _maxIter, double tol, int _verbose) :
		name(_name), maxIter(_maxIter), tolerance(tol), verbose(_verbose) {};
	void operator()(UnconstrainedObjective &_uo);
	int getIter() { return iter; };
};

#endif
