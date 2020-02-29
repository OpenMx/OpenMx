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

struct context {
	GradientOptimizerContext &goc;
	int origeq;
	int eqredundent;
	std::vector<bool> eqmask;
	context(GradientOptimizerContext &_goc) : goc(_goc) {
		eqredundent = 0;
	};
};

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
	goc->numericalGradientWithRef(Epoint);
	Eigen::Map< Eigen::VectorXd > Egrad(grad, n);
	Egrad = goc->grad;
	if (goc->verbose >= 3) {
		mxPrintMat("gradient", Egrad);
	}
	return fit;
}

struct equality_functional {
	GradientOptimizerContext &goc;
	
	equality_functional(GradientOptimizerContext &_goc) : goc(_goc) {};
	
	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.solEqBFun(false);
		result = goc.equality;
	}
	
	template <typename T1, typename T2, typename T3>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result, Eigen::MatrixBase<T3> &jacobian) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.analyticEqJacTmp.resize(jacobian.rows(), jacobian.cols());
		goc.solEqBFun(true);
		result = goc.equality;
		jacobian = goc.analyticEqJacTmp;
		
	}
};

static void nloptEqualityFunction(unsigned m, double* result, unsigned n, const double* x, double* grad, void* f_data)
{
	context &ctx = *(context *)f_data;
	GradientOptimizerContext &goc = ctx.goc;
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::VectorXd Eresult(ctx.origeq);
	Eigen::MatrixXd jacobian(ctx.origeq, n);
	equality_functional ff(goc);
	/*I don't think nloptEqualityFunction is ever called when 'grad' is a null pointer, 
	but I don't want to assume that:*/
	if(!grad){ff(Epoint, Eresult);}
	else{
		ff(Epoint, Eresult, jacobian);
		goc.eqNorm = Eresult.array().abs().sum();
		fd_jacobian<true>(goc.gradientAlgo, goc.gradientIterations, goc.gradientStepSize,
              ff, Eresult, Epoint, jacobian);
		//int rank;
		//filterJacobianRows(jacobian, &rank);
		if (ctx.eqmask.size() == 0) {
			ctx.eqmask.assign(m, false);
			for (int c1=0; c1 < int(m-1); ++c1) {
				for (int c2=c1+1; c2 < int(m); ++c2) {
					bool match = (Eresult[c1] == Eresult[c2] &&
                   (jacobian.row(c1) == jacobian.row(c2)));
					if (match && !ctx.eqmask[c2]) {
						ctx.eqmask[c2] = true;
						++ctx.eqredundent;
						if (goc.verbose >= 2) {
							mxLog("nlopt: eq constraint %d is redundant with %d",
             c1, c2);
						}
					}
				}
			}
			for (int c1=0; c1 < int(m); ++c1) {
				if (ctx.eqmask[c1]) continue;
				if ((jacobian.row(c1).array() == 0.0).all()) {
					ctx.eqmask[c1] = true;
					++ctx.eqredundent;
					if (goc.verbose >= 2) {
						mxLog("nlopt: eq constraint %d is never active", c1);
					}
				}
			}
			if (ctx.eqredundent) {
				if (goc.verbose >= 1) {
					mxLog("nlopt: detected %d redundant equality constraints; retrying",
           ctx.eqredundent);
				}
				nlopt_opt opt = (nlopt_opt) goc.extraData;
				nlopt_force_stop(opt);
			}
		}
	}
	Eigen::Map< Eigen::VectorXd > Uresult(result, m);
	Eigen::Map< Eigen::MatrixXd > Ujacobian(grad, n, m);
	int dx=0;
	for (int cx=0; cx < int(m); ++cx) {
		if (ctx.eqmask[cx]) continue;
		Uresult[dx] = Eresult[cx];
		if (grad) {
			Ujacobian.col(dx) = jacobian.row(cx);
		}
		++dx;
	}
	if (goc.verbose >= 3 && grad) {
		mxPrintMat("eq jacobian", Ujacobian);
	}
}

struct inequality_functional {
	GradientOptimizerContext &goc;
	
	inequality_functional(GradientOptimizerContext &_goc) : goc(_goc) {};
	
	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.myineqFun(false);
		result = goc.inequality;
	}
	
	template <typename T1, typename T2, typename T3>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result, Eigen::MatrixBase<T3> &jacobian) const {
		goc.copyFromOptimizer(x.derived().data());
		goc.analyticIneqJacTmp.resize(jacobian.rows(), jacobian.cols());
		goc.myineqFun(true);
		result = goc.inequality;
		jacobian = goc.analyticIneqJacTmp;
	}
	
};

static void nloptInequalityFunction(unsigned m, double *result, unsigned n, const double* x, double* grad, void* f_data)
{
	GradientOptimizerContext *goc = (GradientOptimizerContext *) f_data;
	Eigen::Map< Eigen::VectorXd > Epoint((double*)x, n);
	Eigen::Map< Eigen::VectorXd > Eresult(result, m);
	Eigen::Map< Eigen::MatrixXd > jacobianDest(grad, n, m);
	Eigen::MatrixXd jacobian(m, n);
	inequality_functional ff(*goc);
	/*nloptInequalityFunction() is routinely called when 'grad' is a null pointer:*/
	if(!grad){ff(Epoint,Eresult);}
	else{
		ff(Epoint, Eresult, jacobian);
		if (goc->verbose >= 2) {
			mxPrintMat("major iteration ineq", Eresult);
		}
		goc->ineqNorm = Eresult.array().abs().sum();
		fd_jacobian<true>(goc->gradientAlgo, goc->gradientIterations, goc->gradientStepSize,
              ff, Eresult, Epoint, jacobian);
		if (!std::isfinite(Eresult.sum())) {
			// infeasible at start of major iteration
			nlopt_opt opt = (nlopt_opt) goc->extraData;
			nlopt_force_stop(opt);
		}
		jacobianDest = jacobian.transpose();
		if (goc->verbose >= 3) {
			mxPrintMat("inequality jacobian", jacobianDest);
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
	omxState *globalState = goc.getState();
	
	nlopt_opt_ptr opt(nlopt_create(NLOPT_LD_SLSQP, goc.numFree));
	goc.extraData = opt.get();
	//local_opt = nlopt_create(NLOPT_LD_SLSQP, n); // Subsidiary algorithm
	
	//nlopt_set_local_optimizer(opt, local_opt);
	nlopt_set_lower_bounds(opt, goc.solLB.data());
	nlopt_set_upper_bounds(opt, goc.solUB.data());
	
	int eq, ieq;
	globalState->countNonlinearConstraints(globalState->numEqC, globalState->numIneqC, false);
	eq = globalState->numEqC;
	ieq = globalState->numIneqC;
	
	// The *2 is there to roughly equate accuracy with NPSOL.
	nlopt_set_ftol_rel(opt, goc.ControlTolerance * 2);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
	
	nlopt_set_min_objective(opt, SLSQP::nloptObjectiveFunction, &goc);
	
	double feasibilityTolerance = Global->feasibilityTolerance;
	SLSQP::context ctx(goc);
	if (eq + ieq) {
		ctx.origeq = eq;
		if (ieq > 0){
			if (goc.verbose >= 2) mxLog("%d inequality constraints", ieq);
			goc.inequality.resize(ieq);
			std::vector<double> tol(ieq, feasibilityTolerance);
			nlopt_add_inequality_mconstraint(opt, ieq, SLSQP::nloptInequalityFunction, &goc, tol.data());
		}
		
		if (eq > 0){
			if (goc.verbose >= 2) mxLog("%d equality constraints", eq);
			goc.equality.resize(eq);
			std::vector<double> tol(eq, feasibilityTolerance);
			nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());
		}
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
	if (ctx.eqredundent) {
		nlopt_remove_equality_constraints(opt);
		eq -= ctx.eqredundent;
		std::vector<double> tol(eq, feasibilityTolerance);
		nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());
		
		code = nlopt_optimize(opt, est, &fit);
	}
	
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
		// No idea what this code means, but often resolved
		// different starting values.
		goc.informOut = INFORM_STARTING_VALUES_INFEASIBLE;
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

double UnconstrainedSLSQPOptimizer::obj(unsigned n, const double *x, double *grad, void *mythis)
{ return ((UnconstrainedSLSQPOptimizer*)mythis)->evaluate(x, grad); }

void UnconstrainedSLSQPOptimizer::operator()(UnconstrainedObjective &_uo)
{
	uo = &_uo;
	opt = nlopt_opt_ptr(nlopt_create(NLOPT_LD_SLSQP, uo->lbound.size()));
	nlopt_set_lower_bounds(opt, uo->lbound.data());
	nlopt_set_upper_bounds(opt, uo->ubound.data());
	nlopt_set_ftol_rel(opt, tolerance);
	nlopt_set_ftol_abs(opt, std::numeric_limits<double>::epsilon());
	nlopt_set_min_objective(opt, obj, this);

	SLSQP::nlopt_slsqp_wdump_ptr wkspc(new nlopt_slsqp_wdump);
	for(int li=0; li<8; li++){
		wkspc->lengths[li] = 0;
	}
	wkspc->realwkspc = (double*)calloc(1, sizeof(double)); //<--Just to initialize it; it'll be resized later.
	opt->work = (nlopt_slsqp_wdump*)&wkspc;
	
	double fit = 0;
	int code = nlopt_optimize(opt, uo->getParamVec(), &fit);

	if (code == NLOPT_INVALID_ARGS) {
		_uo.panic("NLOPT invoked with invalid arguments");
	} else if (code == NLOPT_OUT_OF_MEMORY) {
		_uo.panic("NLOPT ran out of memory");
	} else if (code == NLOPT_ROUNDOFF_LIMITED) {
		_uo.panic("NLOPT_ROUNDOFF_LIMITED"); // only relevant to constrained optimization
	} else if (code < 0) {
		_uo.panic("STARTING_VALUES_INFEASIBLE");
	} else if (code == NLOPT_MAXEVAL_REACHED) {
		_uo.panic("ITERATION_LIMIT");
	}
	if (iter > maxIter) {
		_uo.panic("ITERATION_LIMIT");
	}
}

double UnconstrainedSLSQPOptimizer::evaluate(const double *x, double *grad)
{
	double fit = uo->getFit(x);
	if (grad) {
		uo->getGrad(x, grad);
		Eigen::Map< Eigen::ArrayXd > Egrad(grad, uo->ubound.size());
		if ((!Egrad.isFinite()).any()) uo->panic("gradient has non-finite entries");
		if (verbose >= 2) {
			mxLog("%f", fit);
			mxPrintMat("grad", Egrad);
		}
		if (++iter > maxIter) nlopt_force_stop(opt);
	} else {
		if (verbose >= 3) mxLog("%f", fit);
	}
	return fit;
}

UnconstrainedObjective::UnconstrainedObjective()
	: gwrContext(0)
{}

UnconstrainedObjective::~UnconstrainedObjective()
{
	if (gwrContext) delete gwrContext;
}

void UnconstrainedObjective::getGrad(const double *x, double *out)
{
	if (!gwrContext) {
		gwrContext = new GradientWithRef(1,lbound.size(),GradientAlgorithm_Central,3,1e-3);
	}
	Eigen::Map<const Eigen::VectorXd> Epoint(x, lbound.size());
	Eigen::Map<Eigen::VectorXd> Egrad(out, lbound.size());
	(*gwrContext)([&](double *myPars, int thrId)->double{ return getFit(myPars); },
		   getFit(x), Epoint, Egrad);
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
	SLSQP::context ctx(nmoc->subsidiarygoc);
	if (nmoc->numEqC + nmoc->numIneqC) {
		ctx.origeq = nmoc->numEqC;
		if (nmoc->numIneqC > 0){
			nmoc->subsidiarygoc.inequality.resize(nmoc->numIneqC);
			std::vector<double> tol(nmoc->numIneqC, feasibilityTolerance);
			nlopt_add_inequality_mconstraint(opt, nmoc->numIneqC, SLSQP::nloptInequalityFunction, &(nmoc->subsidiarygoc), tol.data());
		}
		
		if (nmoc->numEqC > 0){
			nmoc->subsidiarygoc.equality.resize(nmoc->numEqC);
			std::vector<double> tol(nmoc->numEqC, feasibilityTolerance);
			nlopt_add_equality_mconstraint(opt, nmoc->numEqC, SLSQP::nloptEqualityFunction, &ctx, tol.data());
		}
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
	if (ctx.eqredundent) {
		nlopt_remove_equality_constraints(opt);
		int eq = nmoc->numEqC - ctx.eqredundent;
		std::vector<double> tol(eq, feasibilityTolerance);
		nlopt_add_equality_mconstraint(opt, eq, SLSQP::nloptEqualityFunction, &ctx, tol.data());
		
		code = nlopt_optimize(opt, est, &fit);
	}
}
