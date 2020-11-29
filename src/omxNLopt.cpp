#include <limits>
#include <cstdlib>
#include <ctype.h>
#include "omxDefines.h"
#include "omxState.h"
#include "omxMatrix.h"
#include "glue.h"
#include "omxNLopt.h"
#include "Compute.h"
#include "ComputeGD.h"
#include "ComputeNM.h"
//#include "finiteDifferences.h"

#include <Eigen/LU>

#include "nlopt.h"
#include "slsqp.h"
#include "nlopt-internal.h"
#include "EnableWarnings.h"

void nlopt_opt_dtor::operator()(struct nlopt_opt_s *opt)
{
	opt->work = 0;
	nlopt_destroy(opt);
}

namespace SLSQP {

	struct nlopt_slsqp_wdump_dtor {
		void operator()(struct nlopt_slsqp_wdump *nsw) {
			free(nsw->realwkspc);
			delete nsw;
		}
	};

	typedef std::unique_ptr< struct nlopt_slsqp_wdump, nlopt_slsqp_wdump_dtor > nlopt_slsqp_wdump_ptr;

static double nloptObjectiveFunction(unsigned n, const double *x, double *grad, void *f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	nlopt_opt opt = (nlopt_opt) goc->extraData;
	int mode = grad != 0;
	double fit = goc->solFun((double*) x, &mode);
	if (grad) {
		goc->iterations += 1;
		if (goc->maxMajorIterations != -1 && goc->iterations >= goc->maxMajorIterations) {
			nlopt_force_stop(opt);
		}
	}
	if (grad && goc->verbose >= 2) {
		mxLog("major iteration fit=%.12f", fit);
	}
	if (mode == -1) {
		if (!goc->feasible) {
			nlopt_force_stop(opt);
		}
		return nan("infeasible");
	}
	if (!grad) return fit;

	Eigen::Map< Eigen::VectorXd > Epoint((double*) x, n);
	Eigen::Map< Eigen::VectorXd > Egrad(grad, n);
	Egrad = goc->grad;
	if (goc->verbose >= 3) {
		mxPrintMat("gradient", Egrad);
	}
	return fit;
}

static void nloptEqualityFunction(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
{
  GradientOptimizerContext &goc = *(GradientOptimizerContext*)f_data;
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobianDest(grad, n, m);
	Eigen::MatrixXd jacobian(m, n);
  goc.copyFromOptimizer(x);

	/*I don't think nloptEqualityFunction is ever called when 'grad' is a null pointer,
	but I don't want to assume that:*/
	if(!grad) {
    goc.evalEq(Eresult.data());
  } else{
    goc.evalEq(Eresult.data(), jacobian.data());
		goc.eqNorm = Eresult.array().abs().sum();
    jacobianDest = jacobian.transpose();
	}
	if (goc.verbose >= 3 && grad) {
		mxPrintMat("eq jacobian", jacobian);
	}
}

static void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobianDest(grad, n, m);
	Eigen::MatrixXd jacobian(m, n);
	/*nloptInequalityFunction() is routinely called when 'grad' is a null pointer:*/
  goc->copyFromOptimizer(x);
	if(!grad) {
		goc->evalIneq(result);
  } else {
		goc->evalIneq(result, jacobian.data());

		if (goc->verbose >= 2) {
			mxPrintMat("major iteration ineq", Eresult);
		}
		goc->ineqNorm = Eresult.array().abs().sum();
		if (!std::isfinite(goc->ineqNorm)) { // infeasible at start of major iteration
			nlopt_opt opt = (nlopt_opt) goc->extraData;
			nlopt_force_stop(opt);
      return;
		}

		jacobianDest = jacobian.transpose();
		if (goc->verbose >= 3) {
			mxPrintMat("inequality jacobian", jacobian);
		}
	}
}

static void omxExtractSLSQPConstraintInfo(nlopt_slsqp_wdump &wkspc, nlopt_opt opt, GradientOptimizerContext &goc){
	int n = opt->n;
	int  n1 = n+1;
	int M = wkspc.M; //<--Total number of constraints (i.e. constraint function elements, not MxConstraint objects)
	if(M <= 0){return;} //<--I can't think of any reason why this function would need to be run if there are no MxConstraints.
	int* lengths = wkspc.lengths;
	double* realwkspc = wkspc.realwkspc;
	int i=0, ro=0, co=0;

	Eigen::MatrixXd Lmat;
	Eigen::MatrixXd Dmat;
	Lmat.setZero(n,n);
	Dmat.setZero(n,n);

	/*SLSQP Jacobian, starting at beginning of realwkspc, is in column-major order; uninitialized "padding" comes after valid elements.
	Jacobian starting at realwkspc[sum(lengths[0:5])] is in row-major order and has no padding(?)*/
	goc.constraintJacobianOut.resize(M,n);
	for(i=0; i<M*n; i++){
		goc.constraintJacobianOut(ro,co) = realwkspc[i];
		ro++;
		if(ro==M){
			ro = 0;
			co++;
		}
	}
	ro=0;
	co=0;

	goc.constraintFunValsOut.resize(lengths[1]);
	for(i=0; i<lengths[1]; i++){
		goc.constraintFunValsOut[i] = realwkspc[lengths[0]+i];
	}

	int Wst = lengths[0]+lengths[1]+lengths[2]+lengths[3]+lengths[4]+lengths[5]+lengths[6];

	goc.LagrMultipliersOut.resize(M);
	for(i=0; i<M; i++){
		goc.LagrMultipliersOut[i] = realwkspc[Wst+i];
	}

	M = M>0 ? M : 1;
	for(i=M; i<M+(n*n1/2); i++){
		if(ro==co){
			Dmat(ro,co) = realwkspc[Wst+i];
			Lmat(ro,co) = 1.0;
		}
		else{Lmat(ro,co) = realwkspc[Wst+i];}
		ro++;
		if(ro==n){
			co++;
			ro = co;
		}
	}
	goc.LagrHessianOut = Lmat * Dmat * Lmat.transpose();
}

//Ideally, this function should recalculate the gradient and Jacobian
//to make sure they're as numerically accurate as possible:
static int constrainedSLSQPOptimalityCheck(GradientOptimizerContext &goc, const double feasTol){
	int code = 0, i=0, arows=0;
	//Nothing to do here if there are no MxConstraints:
	if(!goc.constraintFunValsOut.size() || goc.doingCI()){return(code);}
	//First see if bounds are satisfied:
	for(i=0; i<goc.est.size(); i++){
		if(goc.solLB[i]-goc.est[i] > feasTol || goc.est[i]-goc.solUB[i] > feasTol){
			code = 3;
			break;
		}
	}
	//Now see if constraints are satisfied:
	for(i=0; i<goc.constraintFunValsOut.size(); i++){
		//This works because we set the values of inactive inequality constraints to zero before SLSQP sees them:
		if(fabs(goc.constraintFunValsOut[i]) > 0){
			if(fabs(goc.constraintFunValsOut[i])>feasTol){
				code = 3;
			}
			arows++;
		}
	}
	if(arows){ //&& false){
		int j=0;
		//Per its documentation, this is NPSOL's criterion for first-order optimality conditions:
		double gradThresh = 10 * sqrt(goc.ControlTolerance) * ( 1 + fmax(1 + fabs(goc.getFit()), sqrt(goc.grad.dot(goc.grad)) ) );
		if (goc.verbose >= 2) {
			mxLog("gradient 'threshold': %f", gradThresh);
		}
		Eigen::MatrixXd A(arows, goc.constraintJacobianOut.cols()); //<--Jacobian of active constraints
		A.setZero(arows, goc.constraintJacobianOut.cols());
		for(i=0; i<goc.constraintFunValsOut.size(); i++){
			if(fabs(goc.constraintFunValsOut[i]) > 0){
				A.row(j) = goc.constraintJacobianOut.row(i);
				j++;
			}
		}
		if( !((A.array()==0).all()) ){
			Eigen::FullPivLU< Eigen::MatrixXd > lua(A);
			Eigen::MatrixXd Z = lua.kernel();
			Eigen::VectorXd gz = Z.transpose() * goc.grad;
			double rgnorm = sqrt(gz.dot(gz));
			if (goc.verbose >= 2) {
				mxLog("reduced-gradient norm: %f", rgnorm);
			}
			if(rgnorm > gradThresh){code=6;}
		}
	}

	//Finish:
	return(code);
}

}

void omxInvokeNLOPT(GradientOptimizerContext &goc)
{
	double *est = goc.est.data();
	goc.setEngineName("SLSQP");
	goc.setupSimpleBounds();

	int oldWanted = goc.getWanted();
	goc.setWanted(0);

	nlopt_opt_ptr opt(nlopt_create(NLOPT_LD_SLSQP, goc.numFree));
	goc.extraData = opt.get();
	//local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm

	//nlopt_set_local_optimizer(opt, local_opt);
	nlopt_set_lower_bounds(opt, goc.solLB.data());
	nlopt_set_upper_bounds(opt, goc.solUB.data());

	int eq = goc.EqC.getCount();
	int ieq = goc.IneqC.getCount();

	// The *2 is there to roughly equate accuracy with NPSOL.
	nlopt_set_ftol_rel(opt, goc.ControlTolerance * 2);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());

	nlopt_set_min_objective(opt, SLSQP::nloptObjectiveFunction, &goc);

	double feasibilityTolerance = Global->feasibilityTolerance;
  if (ieq > 0){
    if (goc.verbose >= 2) mxLog("SLSQP: %d inequality constraints", ieq);
    std::vector<double> tol(ieq, feasibilityTolerance);
    nlopt_add_inequality_mconstraint(opt, ieq, SLSQP::nloptInequalityFunction, &goc, tol.data());
  }
  if (eq > 0){
    if (goc.verbose >= 2) mxLog("SLSQP: %d equality constraints", eq);
    std::vector<double> tol(eq, feasibilityTolerance);
    nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &goc, tol.data());
  }

	//The following four lines are only sensible if using SLSQP (noted in case we ever use a different optimizer from the NLOPT collection):
	SLSQP::nlopt_slsqp_wdump_ptr wkspc(new nlopt_slsqp_wdump);
	for(int li=0; li<8; li++){
		wkspc->lengths[li] = 0;
	}
	wkspc->realwkspc = (double*)calloc(1, sizeof(double)); //<--Just to initialize it; it'll be resized later.
	opt->work = wkspc.get();

	double fit = 0;
	int code = nlopt_optimize(opt, est, &fit);

	if (goc.verbose >= 2) mxLog("nlopt_optimize returned %d", code);
	if(code > 0){
		SLSQP::omxExtractSLSQPConstraintInfo(*wkspc, opt, goc);
	}

	goc.setWanted(oldWanted);

	int constrainedCode = SLSQP::constrainedSLSQPOptimalityCheck(goc, feasibilityTolerance);

	if (code == NLOPT_INVALID_ARGS) {
		mxThrow("NLOPT invoked with invalid arguments");
	} else if (code == NLOPT_OUT_OF_MEMORY) {
		mxThrow("NLOPT ran out of memory");
	} else if (code == NLOPT_FORCED_STOP) {
		if (!goc.feasible) {
			goc.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
		} else {
			goc.informOut = INFORM_ITERATION_LIMIT;
		}
	} else if (code == NLOPT_ROUNDOFF_LIMITED) {
		if (goc.eqNorm > feasibilityTolerance || goc.ineqNorm > feasibilityTolerance) {
			goc.informOut = INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE;
		} else if (goc.iterations <= 2) {
			mxThrow("%s: Failed due to singular matrix E or C in LSQ subproblem or "
              "rank-deficient equality constraint subproblem or "
              "positive directional derivative in line search "
				"(eq %.4g ineq %.4g); do you have linearly dependent (i.e., redundant) MxConstraints in your model?",
				goc.getOptName(), goc.eqNorm, goc.ineqNorm);
		} else {
			goc.informOut = INFORM_NOT_AT_OPTIMUM;  // is this correct? TODO
		}
	} else if (code < 0) {
		// Optimizer got stuck.
		if (!goc.feasible) {
			goc.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
		} else {
			goc.informOut = INFORM_NOT_AT_OPTIMUM;  // Can happen when box constraints are too tight.
    }
	} else if (code == NLOPT_MAXEVAL_REACHED) {
		goc.informOut = INFORM_ITERATION_LIMIT;
	} else if(constrainedCode==6){
		goc.informOut = INFORM_NOT_AT_OPTIMUM;
	} else if(constrainedCode==3){
		goc.informOut = INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE;
	}	else {
		goc.informOut = INFORM_CONVERGED_OPTIMUM;
	}
}

void omxInvokeSLSQPfromNelderMead(NelderMeadOptimizerContext* nmoc, Eigen::VectorXd &gdpt)
{
	double *est = gdpt.data();

	nlopt_opt_ptr opt(nlopt_create(NLOPT_LD_SLSQP, nmoc->numFree));
	nmoc->extraData = opt;
	nmoc->subsidiarygoc.extraData = opt;
	nlopt_set_lower_bounds(opt, nmoc->solLB.data());
	nlopt_set_upper_bounds(opt, nmoc->solUB.data());
	nlopt_set_ftol_rel(opt, nmoc->subsidiarygoc.ControlTolerance);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
	nlopt_set_min_objective(opt, nmgdfso, nmoc);

	double feasibilityTolerance = nmoc->NMobj->feasTol;
  int numIneqC = nmoc->IneqC.getCount();
  int numEqC = nmoc->EqC.getCount();
  if (numIneqC > 0){
    std::vector<double> tol(numIneqC, feasibilityTolerance);
    nlopt_add_inequality_mconstraint(opt, numIneqC, SLSQP::nloptInequalityFunction, &(nmoc->subsidiarygoc), tol.data());
  }
  if (numEqC > 0){
    std::vector<double> tol(numEqC, feasibilityTolerance);
    nlopt_add_equality_mconstraint(opt, numEqC, SLSQP::nloptEqualityFunction, &nmoc->subsidiarygoc, tol.data());
  }

	SLSQP::nlopt_slsqp_wdump_ptr wkspc(new nlopt_slsqp_wdump);
	for(int li=0; li<8; li++){
		wkspc->lengths[li] = 0;
	}
	wkspc->realwkspc = (double*)calloc(1, sizeof(double)); //<--Just to initialize it; it'll be resized later.
	opt->work = wkspc.get();

	double fit = 0;
	int code = nlopt_optimize(opt, est, &fit);
	if(nmoc->verbose){
		mxLog("subsidiary SLSQP job returned NLOPT code %d", code);
	}
}
