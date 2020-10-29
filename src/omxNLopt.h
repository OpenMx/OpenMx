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

#endif
