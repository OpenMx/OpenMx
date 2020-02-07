#ifndef _NLOPTCPP_H_
#define _NLOPTCPP_H_

#include "finiteDifferences.h"

void omxInvokeNLOPT(GradientOptimizerContext &goc);

struct nlopt_opt_dtor {
	void operator()(struct nlopt_opt_s *opt);
};

struct nlopt_opt_ptr : std::unique_ptr< struct nlopt_opt_s, struct nlopt_opt_dtor > {
	typedef std::unique_ptr< struct nlopt_opt_s, struct nlopt_opt_dtor > super;
 	nlopt_opt_ptr() : super() {};
 	nlopt_opt_ptr(struct nlopt_opt_s *ptr) : super(ptr) {};
	operator struct nlopt_opt_s *() { return get(); };
};

struct UnconstrainedObjective {
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;
	// If no analytic gradient then can use numerical approx,
	GradientWithRef *gwrContext;

	UnconstrainedObjective();
	virtual ~UnconstrainedObjective();
	virtual double *getParamVec()=0;
	virtual double getFit(const double *)=0;
	virtual void getGrad(const double *, double *);
	virtual void panic(const char *why) { mxThrow("%s", why); };
};

class UnconstrainedSLSQPOptimizer {
	const char *name;
	int maxIter;
	double tolerance;
	int verbose;
	int iter;
	UnconstrainedObjective *uo;
	nlopt_opt_ptr opt;

	double evaluate(const double *x, double *grad);
	static double obj(unsigned n, const double *x, double *grad, void *mythis);
public:
	UnconstrainedSLSQPOptimizer(const char *_name, int _maxIter, double tol, int _verbose) :
		name(_name), maxIter(_maxIter), tolerance(tol), verbose(_verbose) {};
	void operator()(UnconstrainedObjective &_uo);
	int getIter() { return iter; };
};

#endif
