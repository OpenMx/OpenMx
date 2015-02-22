#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "matrix.h"
#include "omxCsolnp.h"
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

struct CSOLNP {

    int flag, flag_NormgZ, minr_rec;
    Eigen::VectorXd LB_e;
    Eigen::VectorXd UB_e;
    Eigen::MatrixXd resP;
    double resLambda;
    Eigen::MatrixXd resHessv;
    Eigen::MatrixXd resY;
    Eigen::MatrixXd sx_Matrix; // search direction
    Eigen::RowVectorXd resGrad;
    int mode;
    GradientOptimizerContext &fit;

	CSOLNP(GradientOptimizerContext &_fit) : fit(_fit) {};
	~CSOLNP() { freeMatrices(); };
    
	void solnp(double *pars, int verbose);

    template <typename T1, typename T2>
    void subnp(Eigen::MatrixBase<T2>& pars, Eigen::MatrixBase<T1>& yy_e, Eigen::MatrixBase<T1>& ob_e, Eigen::MatrixBase<T1>& hessv_e, double lambda, Eigen::MatrixBase<T2>& vscale_e,
                 const Eigen::Array<double, 4, 1> &ctrl, int verbose);

	enum indParam {
		indNumParam=0,
		indHasGradient,
		indHasHessian,
		indHasIneq,
		indHasJacobianIneq,
		indHasEq,
		indHasJacobianEq,
		indVectorLength  // must be last
	};

	Eigen::Array<double, int(indVectorLength), 1> ind;
};

void solnp(double *solPars, GradientOptimizerContext &fit)
{
	CSOLNP context(fit);
	fit.setupIneqConstraintBounds();
	context.solnp(solPars, fit.verbose);
}

void CSOLNP::solnp(double *solPars, int verbose)
{
	fit.informOut = -1;
    LB_e = fit.solLB;
    UB_e = fit.solUB;
    //verbose = 3;
    int i;
    
    flag = 0;
    flag_NormgZ = 0; minr_rec = 0;
    
    double funv;
    double resultForTT;
    double solnp_nfn = 0;
    
    //time_t sec;
    //sec = time (NULL);
    
    int maxit_trace = 0;
    //free(matrices.front().t);
    
    Eigen::Map< Eigen::RowVectorXd > pars(solPars, LB_e.size());

    const int np = pars.size();
    
    ind.setZero();
    ind[indNumParam] = np;
    
    // does not have a function gradient (currently not supported in Rsolnp)
    ind[indHasGradient] = 0;
    //# do function checks and return starting value
    
    mode = 1;
    funv = fit.solFun(pars.data(), &mode);
    
    // does not have a hessian (currently not supported in Rsolnp)
    ind[indHasHessian] = 0;
    // no jacobian (not implemented)
    ind[indHasJacobianIneq] = 0;
    
    // do inequality checks and return starting values
    const int nineq = fit.inequality.size();
    ind[indHasIneq] = nineq > 0;
    
    const int neq = fit.equality.size();
    
    ind[indHasEq] = neq > 0;
    ind[indHasJacobianEq] = 0;
    
    if (verbose >= 2){
        mxLog("ind is: \n");
        for (i = 0; i < ind.size(); i++) mxLog("%f",ind[i]);
    }
    
    Eigen::RowVectorXd ineqx0_e(nineq); ineqx0_e.setZero();
    Eigen::MatrixXd pb_e;
    
    if(nineq) {
        pb_e.setZero(nineq, 2);
        pb_e.col(1) = Eigen::VectorXd::Constant(pb_e.rows(), INF);
        Eigen::MatrixXd pb_cont_e;
        pb_cont_e.setZero(np, 2);
        pb_cont_e.col(0) = LB_e;
        pb_cont_e.col(1) = UB_e;
        pb_e.transposeInPlace();
        pb_cont_e.transposeInPlace();
        if (verbose == 3){
            cout << "pb_e:\n" << pb_e << endl;
            cout << "pb_cont_e:\n" << pb_cont_e << endl;
        }
        Eigen::MatrixXd pbJoined(2, nineq + np);
        pbJoined << pb_e, pb_cont_e;
        pbJoined.transposeInPlace();
        if (verbose == 3){
            cout << "pbJoined:\n" << pbJoined << endl;
        }
        pb_e.resize(pbJoined.rows(), pbJoined.cols());
        pb_e = pbJoined;
        
    } else {
        pb_e.setZero(np, 2);
        pb_e.col(0) = LB_e;
        pb_e.col(1) = UB_e;
    }
    
    double rho   = fit.ControlRho;
    int maxit = fit.ControlMajorLimit;
    int minit = fit.ControlMinorLimit;
    double delta = fit.ControlFuncPrecision;
    double tol   = fit.ControlTolerance;
    
    int tc = nineq + neq;
    
    double j = funv;
    Eigen::VectorXd tt_e(3); tt_e.setZero();
    
    Eigen::MatrixXd constraint_e;
    
    fit.solEqBFun();
    Eigen::RowVectorXd eqv_e(neq);
    eqv_e = fit.equality;
    
    fit.myineqFun();
    Eigen::RowVectorXd ineqv_e(nineq);
    ineqv_e= fit.inequality;
    
    Eigen::MatrixXd lambda_e;
    
    if (tc > 0){
        lambda_e.setZero(tc, 1);
        
        if (nineq){
            if(neq)
            {
                constraint_e.resize(1, eqv_e.size() + ineqv_e.size());
                constraint_e << eqv_e, ineqv_e;
            }
            else{
                constraint_e = ineqv_e;
            }
        }
        else
            constraint_e = eqv_e;
        
        if( ind[indHasIneq] > 0 ) {
            
            // 	tmpv = cbind(constraint[ (neq[0]):(tc[0]-1) ] - .fit.solIneqLB, .fit.solIneqUB - constraint[ (neq + 1):tc ] )
            Eigen::MatrixXd diff1 = constraint_e.block(0, neq, 1, tc-neq).transpose();
            Eigen::MatrixXd diff2 = constraint_e.block(0, neq, 1, tc-neq).transpose();
            Eigen::VectorXd infVec = Eigen::VectorXd::Constant(diff2.rows(), INF);
            diff2 = infVec - diff2;
            Eigen::MatrixXd tmpv_e(nineq, 2);
            tmpv_e.col(0) = diff1;
            tmpv_e.col(1) = diff2;
            Eigen::MatrixXd testMin_e = tmpv_e.rowwise().minCoeff();
            
            if ((testMin_e.array() > 0).all()) {
                ineqx0_e = constraint_e.block(0, neq, 1, tc-neq);
            }
            
            constraint_e.block(0, neq, 1, tc-neq) = constraint_e.block(0, neq, 1, tc-neq) - ineqx0_e;
        }
        
        tt_e[1] = constraint_e.squaredNorm();
        double zeroCheck = tt_e[1] - (10 * tol);
        if( max(zeroCheck, nineq) <= 0 ) {
            rho = 0;
        }
    } // end if tc > 0
    else {
        lambda_e.setZero(1, 1);
    }
    
    Eigen::RowVectorXd p_e;
    
    if (nineq){
        p_e.resize(1, ineqx0_e.size() + pars.size());
        p_e << ineqx0_e, pars;
    }
    else{
        p_e = pars;
    }
    
    Eigen::MatrixXd hessv_e(np + nineq, np + nineq);
    hessv_e.setIdentity();
    
    double mu = np;
    
    int solnp_iter = 0;
    
    Eigen::MatrixXd ob_e;
    Eigen::MatrixXd funvMatrix_e(1,1);
    funvMatrix_e(0, 0) = funv;
    fit.solEqBFun();
    eqv_e = fit.equality;
    
    if ( nineq){
        fit.myineqFun();
        ineqv_e = fit.inequality;
        if(neq){
            ob_e.resize(1, 1 + eqv_e.size() + ineqv_e.size());
            ob_e << funvMatrix_e, eqv_e, ineqv_e;
        }
        else{
            ob_e.resize(1, 1 + ineqv_e.size());
            ob_e << funvMatrix_e, ineqv_e;
        }
        
    }
    else if (neq) {
        ob_e.resize(1, 1 + eqv_e.size());
        ob_e << funvMatrix_e, eqv_e;
    }
    else {
        ob_e.resize(1, 1);
        ob_e = funvMatrix_e;
    }
    
    if(verbose >= 3){
        cout << "ob_e is:" << ob_e << endl;
    }
    
    Eigen::RowVectorXd vscale_e;
    
    while(solnp_iter < maxit){
        solnp_iter = solnp_iter + 1;
        Eigen::Array<double, 4, 1> subnp_ctrl;
        subnp_ctrl[0] = rho;
        subnp_ctrl[1] = minit;
        subnp_ctrl[2] = delta;
        subnp_ctrl[3] = tol;
        
        if ( ind[indHasEq] > 0){
            
            double max = ob_e.block(0, 1, 1, neq).cwiseAbs().maxCoeff();
            Eigen::MatrixXd temp2_e(1, neq);
            temp2_e = temp2_e.setOnes() * max;
            Eigen::MatrixXd temp1_e(1, 1); temp1_e(0, 0) = ob_e(0, 0);
            vscale_e.resize(1, temp1_e.cols() + temp2_e.cols());
            vscale_e << temp1_e, temp2_e;
        }
        else{
            vscale_e.resize(1, 1); vscale_e.setOnes();
        }
        Eigen::RowVectorXd onesMatrix(1, p_e.size()); onesMatrix.setOnes();
        Eigen::RowVectorXd vscale_t = vscale_e;
        vscale_e.resize(1, vscale_t.cols() + onesMatrix.cols());
        vscale_e << vscale_t, onesMatrix;
        
        minMaxAbs(vscale_e, tol);
        
        if (verbose >= 1){
            mxLog("------------------------CALLING SUBNP------------------------");
            cout<< "p_e is:" << p_e<< endl;
            mxLog("lambda information: ");
            cout<< lambda_e << endl;
            mxLog("ob information: ");
            cout<< ob_e << endl;
            mxLog("hessv information: ");
            cout<< hessv_e << endl;
            mxLog("mu information: ");
            mxLog("%f", mu);
            mxLog("vscale information: ");
            cout<< vscale_e << endl;
            mxLog("subnp_ctrl information: ");
            cout<< subnp_ctrl << endl;
            mxLog("------------------------END CALLING SUBNP------------------------");
        }
        
        if (mode == -1)
        {
            fit.informOut = 0;
            memcpy(pars.data(), p_e.data(), pars.size() * sizeof(double));
            return;
        }
        
        sx_Matrix.setZero(p_e.rows(), p_e.cols());
        
        subnp(p_e, lambda_e, ob_e, hessv_e, mu, vscale_e, subnp_ctrl, verbose);
        
        p_e = resP;
        
        if (flag == 1)
        {
            mode = 0;
            funv = fit.solFun(p_e.data(), &mode);
            funvMatrix_e(0, 0) = funv;
            fit.solEqBFun();
            eqv_e = fit.equality;
            if ( nineq){
                fit.myineqFun();
                ineqv_e = fit.inequality;
                if(neq){
                    ob_e.resize(1, 1 + eqv_e.size() + ineqv_e.size());
                    ob_e << funvMatrix_e, eqv_e, ineqv_e;
                }
                else{
                    ob_e.resize(1, 1 + ineqv_e.size());
                    ob_e << funvMatrix_e, ineqv_e;
                }
                
            }
            else if (neq) {
                ob_e.resize(1, 1 + eqv_e.size());
                ob_e << funvMatrix_e, eqv_e;
            }
            else {
                ob_e.resize(1, 1);
                ob_e = funvMatrix_e;
            }

            if ( ind[indHasEq] > 0){
                
                double max = ob_e.block(0, 1, 1, neq).cwiseAbs().maxCoeff();
                Eigen::MatrixXd temp2_e(1, neq);
                temp2_e = temp2_e.setOnes() * max;
                Eigen::MatrixXd temp1_e(1, 1); temp1_e(0, 0) = ob_e(0, 0);
                vscale_e.resize(1, temp1_e.cols() + temp2_e.cols());
                vscale_e << temp1_e, temp2_e;
            }
            else{
                vscale_e.resize(1, 1); vscale_e.setOnes();
            }
            Eigen::RowVectorXd onesMatrix(1, p_e.size()); onesMatrix.setOnes();
            Eigen::RowVectorXd vscale_t = vscale_e;
            vscale_e.resize(1, vscale_t.cols() + onesMatrix.cols());
            vscale_e << vscale_t, onesMatrix;
            
            lambda_e = resY;
            hessv_e = resHessv;
            mu = resLambda;
            subnp(p_e, lambda_e, ob_e, hessv_e, mu, vscale_e, subnp_ctrl, verbose);
        }
        
        lambda_e = resY;
        hessv_e = resHessv;
        mu = resLambda;
        
        Eigen::MatrixXd temp;
        temp = p_e.block(0, nineq, 1, np);
        
        mode = 1;
        funv = fit.solFun(temp.data(), &mode);
        if (mode == -1)
        {
            fit.informOut = 0;
            memcpy(pars.data(), p_e.data(), pars.size() * sizeof(double));
            return;
        }
        
        solnp_nfn = solnp_nfn + 1;

        fit.solEqBFun();
        eqv_e = fit.equality;
        funvMatrix_e(0, 0) = funv;
        
        if ( nineq){
            fit.myineqFun();
            ineqv_e = fit.inequality;
            if(neq){
                ob_e.resize(1, 1 + eqv_e.size() + ineqv_e.size());
                ob_e << funvMatrix_e, eqv_e, ineqv_e;
            }
            else{
                ob_e.resize(1, 1 + ineqv_e.size());
                ob_e << funvMatrix_e, ineqv_e;
            }
            
        }
        else if (neq) {
            ob_e.resize(1, 1 + eqv_e.size());
            ob_e << funvMatrix_e, eqv_e;
        }
        else {
            ob_e.resize(1, 1);
            ob_e = funvMatrix_e;
        }
        
        if (verbose >= 1){
            mxLog("j2 in while: \n");
            mxLog("%f", j);
            cout<< "ob_e is" << ob_e << endl;
        }
        
        resultForTT = (j - ob_e(0, 0)) / max(ob_e.cwiseAbs().maxCoeff(), 1.0);
        tt_e[0] = resultForTT;
        if (verbose >= 1){
            mxLog("resultForTT \n");
            mxLog("%f", resultForTT);
        }
        j = ob_e(0, 0);
        
        if (tc > 0){
            // constraint = ob[ 2:(tc + 1) ]
            constraint_e = ob_e.block(0, 1, 1, tc);
            
            if ( ind[indHasIneq] > 0.5){
                //tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ], pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
                Eigen::MatrixXd subsetOne = constraint_e.block(0, neq, 1, tc-neq) - pb_e.col(0).transpose().block(0, 0, 1, nineq);
                Eigen::MatrixXd subsetTwo = pb_e.col(1).transpose().block(0, 0, 1, nineq);
                subsetTwo -= subsetOne;
                Eigen::MatrixXd tempv(2, nineq);
                tempv.row(0) = subsetOne;
                tempv.row(1) = subsetTwo;
                
                if (tempv.minCoeff() > 0){
                    p_e.block(0, 0, 1, nineq) = constraint_e.block(0, neq, 1, tc-neq);
                }
                constraint_e.block(0, neq, 1, tc-neq) = constraint_e.block(0, neq, 1, tc-neq) - p_e.block(0, 0, 1, nineq);
                
            } // end if (ind[0][3] > 0.5){
            
            tt_e[2] = constraint_e.squaredNorm();
            
            
            if ( tt_e[2] < (10 *tol)){
                rho =0;
                mu = min(mu, tol);
            }
            
            if ( tt_e[2] < (5 * tt_e[1])){
                rho = rho/5;
            }
            
            if ( tt_e[2] > (10 * tt_e[1])){
                rho = 5 * max(rho, sqrt(tol));
            }
            
            Eigen::MatrixXd llist(1, 2);
            
            llist(0, 0) = tol + tt_e[0];
            llist(0, 1) = tt_e[1] - tt_e[2];
            
            
            if (llist.maxCoeff() <= 0){
                lambda_e.resize(1, 1); lambda_e(0, 0) = 0;
                hessv_e = hessv_e.diagonal().asDiagonal();
            }
            
            tt_e[1] = tt_e[2];
            
        } // end if (tc > 0){
                
        Eigen::VectorXd tempTTVals(2);
        tempTTVals[0] = tt_e[0];
        tempTTVals[1] = tt_e[1];
        double vnormValue = tempTTVals.squaredNorm();
        
        if (vnormValue <= tol){
            maxit_trace = maxit;
            maxit = solnp_iter;
        }
        
        if (verbose >= 3)
        {
            mxLog("vnormValue in while \n");
            mxLog("%f", vnormValue);
        }
    } // end while(solnp_iter < maxit){
    
    Eigen::RowVectorXd p_e_copy = p_e;
    p_e.resize(np);
    p_e = p_e_copy.block(0, nineq, 1, np);

    {
        Eigen::VectorXd tempTTVals(2);
        tempTTVals[0] = tt_e[0];
        tempTTVals[1] = tt_e[1];
        double vnormValue = tempTTVals.squaredNorm();

        if (verbose >= 3){
            mxLog("sx_Matrix (search direction) is: \n");
            cout<< sx_Matrix << endl;
        }
        
        if (verbose >= 1) {
		mxLog("vnormValue %.20f, flag_NormgZ=%d, minr_rec=%d",
		      vnormValue, flag_NormgZ, minr_rec);
	}
        if (vnormValue <= tol && flag_NormgZ == 1 && minr_rec == 1){
		double iterateConverge = delta * pow(sx_Matrix.squaredNorm(),(double)2.0);
		double iterateConvergeCond = sqrt(tol) * ((double)1.0 + pow(p_e.squaredNorm(), (double)2.0));
		if (verbose >= 1) {
			mxLog("vnorm(sx_Matrix) is %.20f, iterateConverge is %.20f, iterateConvergeCond is: %.20f",
			      sx_Matrix.squaredNorm(), iterateConverge, iterateConvergeCond);
		}

            if (iterateConverge <= iterateConvergeCond){
		    if (verbose >= 1) { mxLog("The solution converged in %d iterations", solnp_iter); }
		    fit.informOut = INFORM_CONVERGED_OPTIMUM;
            } else {
                if (verbose >= 1){
                    mxLog("The final iterate x satisfies the optimality conditions to the accuracy requested, but the sequence of iterates has not yet converged. CSOLNP was terminated because no further improvement could be made in the merit function.");}
                fit.informOut = INFORM_UNCONVERGED_OPTIMUM;
            }
        }
        else{
            if (solnp_iter == maxit_trace) {
                if (verbose >= 1){
                    mxLog("Exiting after maximum number of iterations. Tolerance not achieved\n");}
                fit.informOut = INFORM_ITERATION_LIMIT;
            } else {
                if (verbose >= 1) { mxLog("Solution failed to converge."); }
                fit.informOut = INFORM_NOT_AT_OPTIMUM;
            }
        }
    }
    
    memcpy(pars.data(), p_e.data(), pars.size() * sizeof(double));
    fit.gradOut.resize(resGrad.size());
    memcpy(fit.gradOut.data(), resGrad.data(), fit.gradOut.size() * sizeof(double));
    fit.hessOut.resize(hessv_e.rows(), hessv_e.cols());
    memcpy(fit.hessOut.data(), hessv_e.data(), fit.hessOut.size() * sizeof(double));
}

template <typename T1, typename T2>
void CSOLNP::subnp(Eigen::MatrixBase<T2>& pars, Eigen::MatrixBase<T1>& yy_e, Eigen::MatrixBase<T1>& ob_e, Eigen::MatrixBase<T1>& hessv_e,
                     double lambda, Eigen::MatrixBase<T2>& vscale_e, const Eigen::Array<double, 4, 1> &ctrl, int verbose)
{
    int yyRows = yy_e.rows();
    //int yyCols = yy.cols;
    double j;
    
    //mxLog("ctrl is: ");
    //for (int ilog = 0; ilog < ctrl.cols; ilog++) mxLog("%f",ctrl.t[ilog]);
    double rho   = ctrl[0];
    int maxit = ctrl[1];
    double delta = ctrl[2];
    double tol =   ctrl[3];
    
    int neq =  fit.equality.size();
    int nineq = fit.inequality.size();
    int np = (int)ind[indNumParam];
    
    double ch = 1;
    
    if (verbose >= 2){
        mxLog("ind inside subnp is: \n");
        for (int i = 0; i < ind.size(); i++) mxLog("%f",ind[i]);
    }
    
    Eigen::Array<double, 3, 1> alp;
    alp.setZero();
    
    int nc = neq + nineq;
    int npic = np + nineq;
    
    Eigen::RowVectorXd p0_e = pars;
    
    if (verbose >= 3)
    {
        mxLog("p0 p0 is: \n");
        for (int i = 0; i < p0_e.size(); i++) mxLog("%f",p0_e[i]);
    }
    
    Eigen::MatrixXd pb_e;
    /*Eigen::Map< Eigen::VectorXd > LB_e(LB.t, LB.cols);
    Eigen::Map< Eigen::VectorXd > UB_e(UB.t, UB.cols);*/
    
    if(nineq) {
        pb_e.setZero(nineq, 2);
        pb_e.col(1) = Eigen::VectorXd::Constant(pb_e.rows(), INF);
        Eigen::MatrixXd pb_cont_e;
        pb_cont_e.setZero(np, 2);
        pb_cont_e.col(0) = LB_e;
        pb_cont_e.col(1) = UB_e;
        pb_e.transposeInPlace();
        pb_cont_e.transposeInPlace();
        if (verbose == 3){
            cout << "pb_e:\n" << pb_e << endl;
            cout << "pb_cont_e:\n" << pb_cont_e << endl;
        }
        Eigen::MatrixXd pbJoined(2, nineq + np);
        pbJoined << pb_e, pb_cont_e;
        pbJoined.transposeInPlace();
        if (verbose == 3){
            cout << "pbJoined:\n" << pbJoined << endl;
        }
        pb_e.resize(pbJoined.rows(), pbJoined.cols());
        pb_e = pbJoined;
    } else {
        pb_e.setZero(np, 2);
        pb_e.col(0) = LB_e;
        pb_e.col(1) = UB_e;
    }
    
    if (verbose >= 3){
        cout<< "pb_e is: \n" << pb_e<< endl;
    }
    
    Eigen::Array<double, 3, 1> sob;
    sob.setZero();
    
    //Matrix yyMatrix = duplicateIt(yy);
    
    ob_e = ob_e.cwiseQuotient(vscale_e.block(0, 0, 1, nc + 1));
    p0_e = p0_e.cwiseQuotient(vscale_e.block(0, neq + 1, 1, nc + np - neq));
    
    if (verbose >= 3){
        cout<< "p0_e: \n"<< p0_e<< endl;
        cout<< "vscale_e: \n"<< vscale_e<< endl;
    }
    
    int mm = 0;
    {
        mm=npic;
        Eigen::MatrixXd pbCopied;
        pbCopied.setZero(pb_e.rows(), pb_e.cols());
        pbCopied.col(0) = vscale_e.block(0, neq + 1, 1, mm).transpose();
        pbCopied.col(1) = vscale_e.block(0, neq + 1, 1, mm).transpose();
        pb_e = pb_e.cwiseQuotient(pbCopied);
    }
    
    if (verbose >= 3){
        mxLog("pb is: \n");
        cout<< pb_e << endl;
    }
    
    // scale the lagrange multipliers and the Hessian
    if( nc > 0) {
        // yy [total constraints = nineq + neq]
        // scale here is [tc] and dot multiplied by yy
        //yy = vscale[ 2:(nc + 1) ] * yy / vscale[ 1 ]
        
        yy_e = vscale_e.block(0, 1, 1, nc).transpose().array() * yy_e.array();
        yy_e = yy_e / vscale_e[0];
        if (verbose >= 3){
            cout << "yy_e:\n" << yy_e << endl;
        }
    }
    
    // hessv [ (np+nineq) x (np+nineq) ]
    // hessv = hessv * (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1)]) ) / vscale[ 1 ]
    
    Eigen::MatrixXd result_e;
    result_e = vscale_e.block(0, neq + 1, 1, nc + np - neq).transpose() * vscale_e.block(0, neq + 1, 1, nc + np - neq);
    hessv_e = hessv_e.cwiseProduct(result_e);
    hessv_e = hessv_e / vscale_e[0];

    j = ob_e(0, 0);
    if (verbose >= 3){
        mxLog("j j is: \n");
        mxLog("%f", j);
    }
    Eigen::MatrixXd a_e;
    if( ind[indHasIneq] > 0){
        if ( ind[indHasEq] <= 0)
        {
            // arrays, rows, cols
            Eigen::MatrixXd negDiag;
            negDiag.setIdentity(nineq, nineq);
            negDiag.diagonal() *= -1;
            //std::cout << "Here is the matrix negDiag:\n" << negDiag << std::endl;
            Eigen::MatrixXd zeroMatrix(nineq, np);
            zeroMatrix.setZero();
            //std::cout << "Here is the matrix  zeroMatrix:\n" << zeroMatrix << std::endl;
            a_e.resize(nineq, np + nineq);
            a_e << negDiag, zeroMatrix;
            //std::cout << "Here is the matrix a_e:\n" << a_e << std::endl;
        }
        else{
            // [ (neq+nineq) x (nineq+np)]
            //a = rbind( cbind( 0 * .ones(neq, nineq), matrix(0, ncol = np, nrow = neq) ),
            //      cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) )
            
            Eigen::MatrixXd zeroMatrix(nineq, np);
            zeroMatrix.setZero();
            //Matrix zeroMatrix = fill(np, nineq, (double)0.0);
            Eigen::MatrixXd firstHalf_e(neq, nineq + np);
            firstHalf_e.setZero();
            //Matrix firstHalf = copy(fill(nineq, neq, (double)0.0), fill(np, neq, (double)0.0));
            Eigen::MatrixXd negDiag;
            negDiag.setIdentity(nineq, nineq);
            negDiag.diagonal() *= -1;
            //Matrix onesMatrix = fill(nineq, 1, (double)-1.0);
            //Matrix negDiag = diag(onesMatrix);
            
            //Matrix secondHalf = copy(negDiag, zeroMatrix);
            
            firstHalf_e.transpose();
            Eigen::MatrixXd secondHalf_e(nineq, np + nineq);
            secondHalf_e << negDiag, zeroMatrix;
            a_e.resize(nineq + np, neq + nineq);
            a_e << firstHalf_e.transpose(), secondHalf_e.transpose();
            a_e.transposeInPlace();
            //a = transpose(copy(transpose(firstHalf), transpose(secondHalf)));
        }
    }	// end 	if(ind[0][3] > 0){
    
    if ( (ind[indHasEq] > 0) && ind[indHasIneq] <= 0 ){
        a_e.resize(neq, np);
        a_e.setZero();
        //a = fill(np, neq, (double)0.0);
    }
    if (ind[indHasEq]<= 0 && (ind[indHasIneq] <= 0)){
        a_e.resize(1, np);
        a_e.setZero();
        //a = fill(np, 1, (double)0.0);
    }

    Eigen::RowVectorXd g_e; g_e.setZero(npic);
    Eigen::RowVectorXd p_e;
    p_e = p0_e.block(0, 0, 1, npic);
    
    Eigen::MatrixXd b_e;
    double funv;
    
    int solnp_nfn = 0;
    double go, reduce = 1e-300;
    int minit;
    double lambdaValue = lambda;
    
    Eigen::MatrixXd constraint_e(1, nc);
    Eigen::MatrixXd y_e;
    
    if (nc > 0) {
        constraint_e = ob_e.block(0, 1, 1, nc);

        for (int i=0; i<np; i++){
            int index = nineq + i;
            p0_e[index] = p0_e[index] + delta;
            Eigen::MatrixXd tmpv_e;
            tmpv_e = p0_e.block(0, nineq, 1, npic - nineq);
            tmpv_e = tmpv_e.array() * vscale_e.block(0, nc+1, 1, np).array();
            if (verbose >= 2){
                mxLog("7th call is \n");
            }
            funv = fit.solFun(tmpv_e.data(), &mode);
            
            fit.solEqBFun();
            fit.myineqFun();
            
            solnp_nfn = solnp_nfn + 1;
            
            Eigen::MatrixXd firstPart_e;
            
            Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
            Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
            Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
            
            if (nineq){
                if(neq)
                {
                    firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                    firstPart_e << funv_e, eqv_e, ineqv_e;
                }
                else{
                    firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                    firstPart_e << funv_e, ineqv_e;
                }
            }
            else if (neq){
                firstPart_e.resize(1, funv_e.size() + eqv_e.size());
                firstPart_e << funv_e, eqv_e;
            }
            else
            {
                firstPart_e.resize(1, funv_e.size());
                firstPart_e << funv_e;
            }
            Eigen::RowVectorXd secondPart_e;
            secondPart_e = vscale_e.block(0, 0, 1, nc+1);
            firstPart_e = firstPart_e * secondPart_e.asDiagonal().inverse();
            ob_e = firstPart_e;
            
            g_e[index] = (ob_e(0, 0)-j) / delta;
            
            if (verbose >= 3){
                mxLog("g is: \n");
                for (int ilog = 0; ilog < g_e.size(); ilog++) mxLog("%f",g_e[ilog]);
                mxLog("a is: \n");
                cout<< a_e << endl;
            }
            
            a_e.col(index) = (ob_e.block(0, 1, 1, nc) - constraint_e).transpose() / delta;
            p0_e[index] = p0_e[index] - delta;
        } // end for (int i=0; i<np, i++){
        
        if (mode == -1)
        {
            funv = 1e24;
            mode = 0;
        }
        
        if(ind[indHasIneq] > 0){
            //constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
            constraint_e.block(0, neq, 1, nineq) = (constraint_e.block(0, neq, 1, nineq) - p0_e.block(0, 0, 1, nineq)).block(0, 0, 1, nineq);
        }
        
        if (false && solvecond(a_e) > 1/DBL_EPSILON) { // this can't be the cheapest way to check TODO
            Rf_error("Redundant constraints were found. Poor intermediate results may result. "
                     "Remove redundant constraints and re-OPTIMIZE.");
        }
        
        b_e = (a_e * p0_e.transpose()).transpose();
        //  b [nc,1]
        b_e -= constraint_e;
        ch = -1;
        alp[0] = tol - constraint_e.cwiseAbs().maxCoeff();
        if (alp[0] <= 0){
            
            ch = 1;
            
        } // end if (alp[0][0] <= 0){
        
        if (alp[0] <= 0){
            int npic_int = npic;
            Eigen::RowVectorXd onesMatrix_e;
            onesMatrix_e.setOnes(1, 1);
            Eigen::RowVectorXd p0_e_copy = p0_e;
            p0_e.resize(p0_e_copy.rows(), p0_e_copy.cols() + onesMatrix_e.cols());
            p0_e << p0_e_copy, onesMatrix_e;
            constraint_e *= (-1.0);
            Eigen::MatrixXd a_e_copy = a_e;
            a_e.resize(a_e.rows(), a_e.cols() + constraint_e.transpose().cols());
            a_e << a_e_copy, constraint_e.transpose();
            Eigen::MatrixXd firstMatrix_e(1, npic);
            firstMatrix_e.setZero();
            Eigen::MatrixXd cx_e(firstMatrix_e.rows(), firstMatrix_e.cols() + onesMatrix_e.cols());
            cx_e << firstMatrix_e, onesMatrix_e;
            Eigen::MatrixXd dx_e(npic + 1, 1);
            dx_e.setOnes();
            go = 1;
            minit = 0;
            
            while(go >= tol)
            {
                minit = minit + 1;
                Eigen::MatrixXd gap_e(mm, 2);
                gap_e.setZero();
                gap_e.col(0) = p0_e.block(0, 0, 1, mm).transpose() - pb_e.col(0);
                gap_e.col(1) = p0_e.block(0, 0, 1, mm).transpose() - pb_e.col(1);
                rowSort_e(gap_e);
                dx_e.transpose().block(0, 0, 1, mm) = gap_e.col(0).transpose().block(0, 0, 1, mm);
                dx_e(npic_int, 0) = p0_e(0, npic_int);
                Eigen::MatrixXd argum1_e;
                argum1_e = a_e * dx_e.asDiagonal();
                argum1_e.transposeInPlace();
                Eigen::MatrixXd argum2_e;
                argum2_e = cx_e.asDiagonal() * dx_e;
                y_e = argum1_e.colPivHouseholderQr().solve(argum2_e);
                Eigen::MatrixXd cx_e_r;
                cx_e_r = cx_e.transpose() - (a_e.transpose() * y_e);
                dx_e = (cx_e_r.asDiagonal() * dx_e).asDiagonal() * dx_e;
                Eigen::MatrixXd v_e = dx_e.transpose();
                int indexx = npic;
                
                if (v_e(0, indexx) > 0)
                {
                    double z = p0_e(indexx)/v_e(0, indexx);
                    
                    for (int i=0; i<mm; i++)
                    {
                        if(v_e(0, i) < 0)
                        {
                            z = min(z, -(pb_e(i, 1) - p0_e(i))/v_e(0, i));
                            
                        }
                        else if(v_e(0, i) > 0)
                        {
                            
                            z = min(z, (p0_e(i) - pb_e(i, 0))/v_e(0, i));
                        }
                    }
                    
                    if(z < (p0_e(indexx)/v_e(0, indexx))) {
                        z *= 0.9;
                    }
                    
                    p0_e -= v_e * z;
                    go = p0_e(indexx);
                    if(minit >= 10){
                        go = 0;
                    }
                }
                else{
                    go = 0;
                    minit = 10;
                }
            }// end while(go >= tol)
            
            if (minit >= 10){
                mxLog("The linearized problem has no feasible solution. The problem may not be feasible.");
            }
            
            int h;
            Eigen::MatrixXd a_e_c(nc, npic);
            
            for (h = 0; h<a_e.rows(); h++)
            {
                a_e_c.row(h) = a_e.row(h).block(0, 0, 1, npic);
            }
            a_e.resize(a_e_c.rows(), a_e_c.cols());
            a_e = a_e_c;
            b_e = (a_e * p0_e.block(0, 0, 1, npic).transpose()).transpose();
        }// end if(M(alp, 0, 0) <= 0)
    } // end if (nc > 0){
    
    p_e = p0_e.block(0, 0, 1, npic);
    
    if (nc == 0){
        y_e.resize(1,1);
        y_e(0, 0) = 0;
    }

    if (verbose >= 3){
        mxLog("p is: \n");
        cout << p_e << endl;
    }
    
    if (ch > 0){
        
        Eigen::MatrixXd tmpv_e;
        tmpv_e = p_e.block(0, nineq, 1, npic-nineq);
        tmpv_e = tmpv_e.array() * vscale_e.block(0, nc+1, 1, np).array();
        funv = fit.solFun(tmpv_e.data(), &mode);
        if (verbose >= 3){
            mxLog("funv is: \n");
            mxLog("%f", funv);
        }
        
        if (mode == -1)
        {
            funv = 1e24;
            mode = 0;
        }
        
        fit.solEqBFun();
        
        fit.myineqFun();
        
        solnp_nfn = solnp_nfn + 1;
        Eigen::MatrixXd firstPart_e;
        Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
        Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
        Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
        
        if (nineq){
            if (neq){
                firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                firstPart_e << funv_e, eqv_e, ineqv_e;
            }
            else{
                firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                firstPart_e << funv_e, ineqv_e;
            }
        }
        else if (neq){
            firstPart_e.resize(1, funv_e.size() + eqv_e.size());
            firstPart_e << funv_e, eqv_e;
        }
        else
        {
            firstPart_e.resize(1, funv_e.size());
            firstPart_e << funv_e;
        }
        
        Eigen::RowVectorXd secondPart_e;
        secondPart_e = vscale_e.block(0, 0, 1, nc+1);
        firstPart_e = firstPart_e * secondPart_e.asDiagonal().inverse();
        ob_e = firstPart_e;
        
    } // end of if (ch>0)
    
    if (verbose >= 3){
        mxLog("ob is: \n");
        cout<< ob_e << endl;
    }
    
    j = ob_e(0, 0);
    
    if (ind[indHasIneq] > 0){
        ob_e.block(0, neq+1, 1, nc-neq) -= p_e.block(0, 0, 1, nineq);
    }
    
    if (nc > 0){
        Eigen::MatrixXd result_e = ob_e.block(0, 1, 1, nc);
        result_e -= (a_e * p_e.transpose()).transpose();
        result_e += b_e;
        ob_e.block(0, 1, 1, nc) = result_e;
        double vnormTerm = ob_e.block(0, 1, 1, nc).squaredNorm() * ob_e.block(0, 1, 1, nc).squaredNorm();
        double dotProductTerm = yy_e.transpose().row(0).dot(ob_e.block(0, 1, 1, nc).row(0));
        j = ob_e(0, 0) - dotProductTerm + rho * vnormTerm;
    }
    
    minit = 0;
    Eigen::MatrixXd yg_e;
    Eigen::MatrixXd yg_rec(1, 2);
    Eigen::MatrixXd sx_e;
    sx_e.setZero(p_e.rows(), p_e.cols());
    Eigen::MatrixXd obm_e;
    
    while (minit < maxit){
        minit = minit + 1;
        if (ch > 0){
            
            for (int i=0; i<np; i++){
                int index = nineq + i;
                p_e[index] = p_e[index] + delta;
                Eigen::MatrixXd tmpv_e = p_e.block(0, nineq, 1, npic - nineq).array() * vscale_e.block(0, nc+1, 1, np).array();
                if (verbose >= 3){
                    mxLog("9th call is \n");
                }
                mode = 0;
                funv = fit.solFun(tmpv_e.data(), &mode);
                if (verbose >= 3){
                    mxLog("funv is: \n");
                    mxLog("%f", funv);
                }
                
                if (mode == -1)
                {
                    funv = 1e24;
                    mode = 0;
                }
                fit.solEqBFun();
                fit.myineqFun();
                
                solnp_nfn = solnp_nfn + 1;
                
                Eigen::MatrixXd firstPart_e;
                Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
                Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
                Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
                
                if (nineq){
                    if (neq){
                        firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                        firstPart_e << funv_e, eqv_e, ineqv_e;
                    }
                    else{
                        firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                        firstPart_e << funv_e, ineqv_e;
                    }
                }
                else if (neq){
                    firstPart_e.resize(1, funv_e.size() + eqv_e.size());
                    firstPart_e << funv_e, eqv_e;
                }
                else
                {
                    firstPart_e.resize(1, funv_e.size());
                    firstPart_e << funv_e;
                }
                
                Eigen::RowVectorXd secondPart_e;
                secondPart_e = vscale_e.block(0, 0, 1, nc+1);
                
                firstPart_e = firstPart_e * secondPart_e.asDiagonal().inverse();
                obm_e = firstPart_e;
                
                if (verbose >= 3){
                    mxLog("j is: \n");
                    mxLog("%f", j);
                }
                
                if (ind[indHasIneq] > 0.5){
                    obm_e.block(0, neq+1, 1, nc-neq) -= p_e.block(0, 0, 1, nineq);
                }
                
                double obm = obm_e(0, 0);
                if (nc > 0){
                    Eigen::MatrixXd result_e = obm_e.block(0, 1, 1, nc);
                    result_e -= (a_e * p_e.transpose()).transpose();
                    result_e += b_e;
                    obm_e.block(0, 1, 1, nc) = result_e;
                    double vnormTerm = obm_e.block(0, 1, 1, nc).squaredNorm() * obm_e.block(0, 1, 1, nc).squaredNorm();
                    double dotProductTerm = yy_e.transpose().row(0).dot(obm_e.block(0, 1, 1, nc).row(0));
                    obm = obm_e(0, 0) - dotProductTerm + rho * vnormTerm;
                }
                
                if (verbose >= 3)   mxLog("obm is: %.20f", obm);
                
                g_e[index] = (obm - j)/delta;
                p_e[index] = p_e[index] - delta;
                
                if (verbose >= 3){
                    mxLog("g is: \n");
                    for (int ilog = 0; ilog < g_e.size(); ilog++) mxLog("%f",g_e[ilog]);
                    mxLog("p is: \n");
                    for (int ilog = 0; ilog < p_e.size(); ilog++) mxLog("%f",p_e[ilog]);
                }
            } // end for (i=0; i<np; i++){
            
            if (ind[indHasIneq] > 0.5){
                Eigen::RowVectorXd temp;
                temp.setZero(1, nineq);
                g_e.block(0, 0, 1, nineq) = temp;
            }
        } // end if (ch > 0){
        
        if (verbose >= 3){
            mxLog("yg is: \n");
            cout<< yg_e << endl;
        }
        
        if (minit > 1){
            if (verbose >= 3){
                cout<< "hessv_e is: \n"  << hessv_e << endl;
                cout<< "yg_e is: \n"  << yg_e << endl;
                cout<< "sx_e is: \n"  << sx_e << endl;
                cout<< "g_e is: \n"  << g_e << endl;
                cout<< "p_e is: \n"  << p_e << endl;
            }
            
            yg_e = g_e - yg_e;
            sx_e = p_e - sx_e;
            if (verbose >= 3){
                cout<< "sx_e rows is: \n"<< sx_e.rows() << endl;
                cout<< "sx_e cols is: \n"<< sx_e.cols() << endl;
            }
            Eigen::MatrixXd sc_m1 = (sx_e * hessv_e) * sx_e.transpose();
            if (verbose >= 3)
                cout<< "sc_:m1 \n" << sc_m1 << endl;
            Eigen::MatrixXd sc_m2 = sx_e * yg_e.transpose();
            if (verbose >= 3)
                cout<< "sc_m2: \n" << sc_m2 << endl;
            Eigen::RowVectorXd sc_e(2);
            sc_e[0] = sc_m1(0, 0);
            sc_e[1] = sc_m2(0, 0);
            if (verbose >= 3)
                cout<< "sc_e: \n" << sc_e << endl;
            if ((sc_e[0] * sc_e[1]) > 0){
                //hessv  = hessv - ( sx %*% t(sx) ) / sc[ 1 ] + ( yg %*% t(yg) ) / sc[ 2 ]
                Eigen::MatrixXd sx_t = sx_e.transpose();
                sx_e.resize(hessv_e.rows(), sx_t.cols());
                sx_e = hessv_e * sx_t;
                
                if (verbose >= 3)
                    cout<< "sx_e: \n" << sx_e << endl;
                Eigen::MatrixXd sxMatrix = sx_e * sx_e.transpose();
                sxMatrix /= sc_e[0];
                if (verbose >= 3)
                    cout<< "sxMatrix: \n" << sxMatrix << endl;
                Eigen::MatrixXd ygMatrix = yg_e.transpose() * yg_e;
                ygMatrix /= sc_e[1];
                if (verbose >= 3)
                    cout<< "ygMatrix: \n" << ygMatrix << endl;
                hessv_e -= sxMatrix;
                hessv_e += ygMatrix;
            }
        }
        if (verbose >= 3){
            cout<< "yg_e: \n" << yg_e << endl;
            cout<< "sx_e: \n" << sx_e << endl;
            cout<< "hessv_e: \n" << hessv_e << endl;
        }
        
        Eigen::MatrixXd dx_e(1, npic);
        dx_e.setOnes();
        dx_e *= 0.01;
        
        {
            Eigen::MatrixXd gap_e(pb_e.rows(), pb_e.cols());
            gap_e.setZero();
            gap_e.col(0) = p_e.block(0, 0, 1, mm).transpose() - pb_e.col(0);
            gap_e.col(1) = pb_e.col(1) - p_e.block(0, 0, 1, mm).transpose();
            rowSort_e(gap_e);
            Eigen::MatrixXd temp(mm, 1);
            temp.setOnes();
            Eigen::MatrixXd gap_eTemp(mm, 1);
            gap_eTemp = gap_e.col(0) + (temp * sqrt(DBL_EPSILON));
            gap_e.resize(mm, 1);
            gap_e = gap_eTemp;
            dx_e.block(0, 0, 1, mm) = temp.cwiseQuotient(gap_e).transpose();
            if (verbose >= 3){
                mxLog("dx is: \n");
                cout << dx_e << endl;
            }
        }
        
        go = -1;
        lambdaValue = lambdaValue/10.0;
        
        if (verbose >= 3){
            mxLog("lambdaValue is: \n");
            mxLog("%.20f", lambdaValue);
            mxLog("hessv is: \n");
            cout<< hessv_e << endl;
        }
        
        while(go <= 0){
            Eigen::RowVectorXd dxDiagValues(dx_e.cols());
            dxDiagValues = dx_e.cwiseProduct(dx_e);
            Eigen::MatrixXd cz_e;
            cz_e = dxDiagValues.asDiagonal();
            if (verbose >= 3){
                cout<< "cz_e size is:	"<< cz_e.rows()<< endl;
                cout<< "cz_e size is:	"<< cz_e.cols()<< endl;
                cout<< "cz_e diag is: \n"  << cz_e << endl;
                cout<< "dxDiagValues.asDiagonal() is: \n" << dxDiagValues.asDiagonal().rows() << endl;
                cout<< "dxDiagValues.asDiagonal() is: \n" << dxDiagValues.asDiagonal().cols() << endl;
                cout<< "cz_e lambda is: \n"  << cz_e * lambdaValue << endl;
            }
            cz_e = hessv_e + (cz_e * lambdaValue);
            if (verbose >= 3){
                cout<< "cz_e is: \n"  << cz_e << endl;
            }
            Eigen::MatrixXd cz_chol = cz_e.llt().matrixL();
            cz_chol.transposeInPlace();

            if (verbose >= 3){
                cout<< "cz_chol is: \n"  << cz_chol << endl;
            }
            
            if (!R_FINITE((cz_e.maxCoeff())))
            {
                mxLog("here in findMax");
                flag = 1;
                p_e = p_e.cwiseProduct(vscale_e.block(0, neq+1, 1, nc+np-neq));
                if (nc > 0){ y_e.resize(1, 1); y_e(0, 0) = 0;}
                hessv_e = hessv_e.cwiseQuotient(vscale_e.block(0, neq+1, 1, nc+np-neq) * vscale_e.block(0, neq+1, 1, nc+np-neq).transpose()) *vscale_e(0);
                resP = p_e;
                resY = y_e;
                resHessv = hessv_e;
                resLambda = lambda;
                resGrad = g_e;
            }
            
            Eigen::MatrixXd cz_inv;
            cz_inv = cz_chol.inverse();
            
            if (verbose >= 3){
                mxLog("cz.rows: %d", cz_chol.rows());
                mxLog("cz.cols: %d", cz_chol.cols());
                cout<< "cz_chol.inverse is: \n"  << cz_inv << endl;
            }
            
            if (verbose >= 3){
                cout<< "g_e is: \n"  << g_e << endl;
            }
            yg_e.resize(cz_inv.cols(), g_e.rows());
            yg_e = cz_inv.transpose() * g_e.transpose();
            if (minit == 1) yg_rec(0, 0) = yg_e.squaredNorm();
            
            if (verbose >= 3){
                cout<< "yg_e is: \n"  << yg_e << endl;
            }
            Eigen::MatrixXd u_e;
            if (nc <= 0){
                u_e = (cz_inv * (-1.0)) * yg_e;
                u_e.transposeInPlace();
                if (verbose >= 3){
                    cout<< "u_e is: \n"  << u_e << endl;
                }
            }
            else{
                //y = qr.solve(t(cz) %*% t(a), yg)
                Eigen::MatrixXd argum1_e;
                argum1_e = cz_inv.transpose() * a_e.transpose();
                if (verbose >= 3){
                    cout<< "argum1_e is: \n"  << argum1_e << endl;
                }
                Eigen::MatrixXd solution;
                
                solution = QRdsolve(argum1_e, yg_e);
                
                if (verbose >= 3){
                    cout<< "solution is: \n"  << solution << endl;
                }
                y_e.resize(solution.cols(), solution.rows());
                y_e = solution.transpose();
                u_e = (cz_inv * (-1.0)) * (yg_e - (argum1_e * solution));
                u_e.transposeInPlace();
                if (verbose >= 3){
                    cout<< "u_e is: \n" << u_e << endl;
                }
            }
            
            p0_e.resize(npic);
            p0_e = u_e.block(0, 0, 1, npic) + p_e;
            
            if (verbose >= 3){
                cout<< "p0_e is: \n"  << p0_e << endl;
                cout<< "pb_e is: \n"  << pb_e.col(0) << endl;
                cout<< "mm is: \n" << mm << endl;
            }
            
            {
                Eigen::MatrixXd listPartOne = p0_e.block(0, 0, 1, mm).transpose() - pb_e.col(0);
                Eigen::MatrixXd listPartTwo = pb_e.col(1) - p0_e.block(0, 0, 1, mm).transpose();
                Eigen::MatrixXd llist(listPartOne.rows(), listPartOne.cols() + listPartTwo.cols());
                llist << listPartOne, listPartTwo;
                go = llist.minCoeff();
                lambdaValue = 3 * lambdaValue;
                if (verbose >= 3){
                    mxLog("go is: \n");
                    mxLog("%f", go);
                    mxLog("lambdaValue is: \n");
                    mxLog("%f", lambdaValue);
                    
                }
            }
        } // end while(go <= 0){
        
        alp[0] = 0;
        Eigen::MatrixXd ob1_e = ob_e;
        Eigen::MatrixXd ob2_e = ob1_e;
        sob[0] = j;
        sob[1] = j;
        
        if (verbose >= 3){
            mxLog("sob is is: \n");
            for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f", sob[ilog]);
        }
        
        Eigen::MatrixXd ptt_e(p_e.cols(), p_e.rows() + p_e.rows());
        ptt_e << p_e.transpose(), p_e.transpose();
        alp[2] = 1.0;
        
        Eigen::MatrixXd ptt_temp(ptt_e.rows(), ptt_e.cols() + p0_e.rows());
        ptt_temp << ptt_e, p0_e.transpose();
        ptt_e.resize(ptt_temp.rows(), ptt_temp.cols());
        ptt_e = ptt_temp;
        if (verbose >= 3){
            cout<< "ptt_e is: \n"  << ptt_e << endl;
        }
        Eigen::MatrixXd pttCol;
        pttCol = ptt_e.col(2);
        Eigen::MatrixXd tmpv_e = pttCol.transpose().block(0, nineq, 1, npic - nineq).cwiseProduct(vscale_e.block(0, nc+1, 1, np));
        
        if (verbose >= 3){
            cout<< "tmpv_e is: \n"  << tmpv_e << endl;
        }
        
        mode = 1;
        funv = fit.solFun(tmpv_e.data(), &mode);
        if (verbose >= 3){
            mxLog("hessv is: \n");
            cout<< hessv_e << endl;
            mxLog("g is: \n");
            for (int ilog = 0; ilog < g_e.size(); ilog++) mxLog("%f",g_e[ilog]);
            mxLog("funv is: \n");
            mxLog("%f", funv);
        }
        
        if (mode == -1)
        {
            funv = 1e24;
            mode = 0;
        }
        
        fit.solEqBFun();
        fit.myineqFun();
        
        solnp_nfn = solnp_nfn + 1;
        
        Eigen::MatrixXd firstPart_e;
        Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
        Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
        Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
        if (nineq){
            if(neq)
            {
                firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                firstPart_e << funv_e, eqv_e, ineqv_e;
            }
            else{
                firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                firstPart_e << funv_e, ineqv_e;
            }
        }
        else if (neq){
            firstPart_e.resize(1, funv_e.size() + eqv_e.size());
            firstPart_e << funv_e, eqv_e;
        }
        else
        {
            firstPart_e.resize(1, funv_e.size());
            firstPart_e << funv_e;
        }
        Eigen::RowVectorXd secondPart_e;
        secondPart_e = vscale_e.block(0, 0, 1, nc+1);
        firstPart_e = firstPart_e * secondPart_e.asDiagonal().inverse();
        Eigen::MatrixXd ob3_e = firstPart_e;
        if (verbose >= 3){
            cout<< "ob3_e is: \n"  << ob3_e << endl;
        }
        sob[2] = ob3_e(0, 0);
        
        if (ind[indHasIneq] > 0.5){
            // ob3[ (neq + 2):(nc + 1) ] = ob3[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq, 3 ]
            Eigen::MatrixXd partOne = ob3_e.block(0, neq+1, 1, nc-neq);
            Eigen::MatrixXd partTwo = ptt_e.col(2).transpose().block(0, 0, 1, nineq);
            ob3_e.block(0, neq+1, 1, nc-neq) = partOne - partTwo;
        }
        
        if (nc > 0){
            //sob[ 3 ] = ob3[ 1 ] - t(yy) %*% ob3[ 2:(nc + 1) ] + rho * .vnorm(ob3[ 2:(nc + 1) ]) ^ 2
            Eigen::MatrixXd result_e = ob3_e.block(0, 1, 1, nc);
            result_e -= (a_e * ptt_e.col(2)).transpose();
            result_e += b_e;
            ob3_e.block(0, 1, 1, nc) = result_e;
            double vnormTerm = ob3_e.block(0, 1, 1, nc).squaredNorm() * ob3_e.block(0, 1, 1, nc).squaredNorm();
            double dotProductTerm = yy_e.transpose().row(0).dot(ob3_e.block(0, 1, 1, nc).row(0));
            sob[2] = ob3_e(0, 0) - dotProductTerm + (rho * vnormTerm);
        }

        go = 1;
        
        while(go > tol){
            alp[1] = (alp[0] + alp[2]) / 2.0;
            
            ptt_e.col(1) = (p_e * (1 - alp[1])) + p0_e * alp[1];
            Eigen::MatrixXd tmpv_e = ptt_e.col(1).transpose().block(0, nineq, 1, npic - nineq).cwiseProduct(vscale_e.block(0, nc+1, 1, np));

            if (verbose >= 3){
                cout<< "tmpv_e is: \n"  << tmpv_e << endl;
            }
            
            if (verbose >= 3){
                mxLog("11th call is \n");
            }
            
            mode = 0;
            funv = fit.solFun(tmpv_e.data(), &mode);
            if (verbose >= 3){
                mxLog("funv is: \n");
                mxLog("%f", funv);
            }
            
            if (mode == -1)
            {
                funv = 1e24;
                mode = 0;
            }
            
            fit.solEqBFun();
            fit.myineqFun();
            
            solnp_nfn = solnp_nfn + 1;
            Eigen::MatrixXd firstPart_e;
            Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
            Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
            Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
            if (nineq){
                if(neq)
                {
                    firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                    firstPart_e << funv_e, eqv_e, ineqv_e;
                }
                else{
                    firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                    firstPart_e << funv_e, ineqv_e;
                }
            }
            else if (neq){
                firstPart_e.resize(1, funv_e.size() + eqv_e.size());
                firstPart_e << funv_e, eqv_e;
            }
            else
            {
                firstPart_e.resize(1, funv_e.size());
                firstPart_e << funv_e;
            }
            Eigen::RowVectorXd secondPart_e;
            secondPart_e = vscale_e.block(0, 0, 1, nc+1);
            firstPart_e = firstPart_e * secondPart_e.asDiagonal().inverse();
            ob2_e = firstPart_e;
            
            if (verbose >= 3){
                cout<< "ob2_e is: \n"  << ob2_e << endl;
            }
            
            sob[1] = ob2_e(0, 0);
            if (verbose >= 3){
                mxLog("sob is: \n");
                for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
            }
            if (ind[indHasIneq] > 0.5){
                Eigen::MatrixXd partOne = ob2_e.block(0, neq+1, 1, nc-neq);
                Eigen::MatrixXd partTwo = ptt_e.col(1).transpose().block(0, 0, 1, nineq);
                ob2_e.block(0, neq+1, 1, nc-neq) = partOne - partTwo;
            }
            if (verbose >= 3){
                cout<< "ob2_e is: \n"  << ob2_e << endl;
            }
            if (nc > 0){
                Eigen::MatrixXd result_e = ob2_e.block(0, 1, 1, nc);
                result_e -= (a_e * ptt_e.col(1)).transpose();
                result_e += b_e;
                ob2_e.block(0, 1, 1, nc) = result_e;
                double vnormTerm = ob2_e.block(0, 1, 1, nc).squaredNorm() * ob2_e.block(0, 1, 1, nc).squaredNorm();
                Eigen::MatrixXd temp = ob2_e.block(0, 1, 1, nc);
                double dotProductTerm = yy_e.transpose().row(0).dot(temp.row(0));
                sob[1] = ob2_e(0, 0) - dotProductTerm + (rho * vnormTerm);
            }
            
            if (verbose >= 3){
                cout<< "sob is: \n"  << sob << endl;
            }
            const double sobMax = sob.maxCoeff();
            if (verbose >= 3){
                mxLog("sobMax is: %f", sobMax);
            }
            if (sobMax < j){
                go = tol * (sobMax - sob.minCoeff()) / (j - sobMax);
            }
            
            const bool condif1 = (sob[1] >= sob[0]);
            const bool condif2 = (sob[0] <= sob[2]) && (sob[1] < sob[0]);
            const bool condif3 = (sob[1] <  sob[0]) && (sob[0] > sob[2]);
            
            if (condif1){
                sob[2] = sob[1];
                ob3_e = ob2_e;
                alp[2] = alp[1];
                ptt_e.col(2) = ptt_e.col(1);
                
                if (verbose >= 3){
                    cout<< "sob is: \n"  << sob << endl;
                    cout<< "ob3_e is: \n"  << ob3_e << endl;
                    cout<< "alp[2] is: \n"  << alp << endl;
                    cout<< "ptt_e is: \n"  << ptt_e << endl;
                }
            }
            
            if (condif2){
                sob[2] = sob[1];
                ob3_e = ob2_e;
                alp[2] = alp[1];
                ptt_e.col(2) = ptt_e.col(1);
                
                if (verbose >= 3){
                    cout<< "sob is: \n"  << sob << endl;
                    cout<< "ob3_e is: \n"  << ob3_e << endl;
                    cout<< "alp[2] is: \n"  << alp << endl;
                    cout<< "ptt_e is: \n"  << ptt_e << endl;
                }
            }
            
            if (condif3){
                sob[0] = sob[1];
                ob1_e = ob2_e;
                alp[0] = alp[1];
                ptt_e.col(0) = ptt_e.col(1);
                
                if (verbose >= 3){
                    cout<< "sob is: \n"  << sob << endl;
                    cout<< "ob1_e is: \n"  << ob1_e << endl;
                    cout<< "alp[2] is: \n"  << alp << endl;
                    cout<< "ptt_e is: \n"  << ptt_e << endl;
                }
            }
            
            if (go >= tol){
                go = alp[2] - alp[0];
                if (verbose >= 3){
                    mxLog("go is: \n");
                    mxLog("%f", go);
                }
            }
            
        } // 	while(go > tol){
        
        if (verbose >= 3){
            mxLog("go is: \n");
            mxLog("%.16f", go);
        }

        sx_Matrix = sx_e;
        sx_e.resize(p_e.rows(), p_e.cols());
        sx_e = p_e;
        yg_e.resize(g_e.rows(), g_e.cols());
        yg_e = g_e;
        
        if (verbose >= 3){
            mxLog("sx is: \n");
            cout<< sx_e << endl;
            mxLog("yg is: \n");
            cout<< yg_e << endl;
        }
        ch = 1;
        
        double obn = sob.minCoeff();
        if (verbose >= 3){
            mxLog("obn is: \n");
            mxLog("%f", obn);
        }
        if (j <= obn){
            maxit = minit;
        }
        if (verbose >= 3){
            mxLog("j is: \n");
            mxLog("%f", j);
        }
        double reduce = (j - obn) / ((double)1.0 + (double)fabs(j));
        if (verbose >= 3){
            mxLog("reduce is: \n");
            mxLog("%f", reduce);
        }
        if (reduce < tol){
            maxit = minit;
        }
        
        const bool condif1 = (sob[0] <  sob[1]);
        const bool condif2 = (sob[2] <  sob[1]) && (sob[0] >= sob[1]);
        const bool condif3 = (sob[0] >= sob[1]) && (sob[2] >= sob[1]);
        
        if (condif1){
            j = sob[0];
            p_e = ptt_e.col(0).transpose();
            ob_e = ob1_e;
            if (verbose >= 3){
                mxLog("condif1\n");
                mxLog("j is: \n");
                mxLog("%f", j);
                cout<< "p_e is: \n"  << p_e << endl;
                cout<< "ob_e is: \n"  << ob_e << endl;
            }
        }
        
        if (condif2){
            
            j = sob[2];
            p_e = ptt_e.col(2).transpose();
            ob_e = ob3_e;
            if (verbose >= 3){
                mxLog("condif2\n");
                mxLog("j is: \n");
                mxLog("%f", j);
                cout<< "p_e is: \n"  << p_e << endl;
                cout<< "ob_e is: \n"  << ob_e << endl;
            }
            
        }
        
        if (condif3){
            j = sob[1];
            p_e = ptt_e.col(1).transpose();
            ob_e = ob2_e;
            if (verbose >= 3){
                mxLog("condif3\n");
                mxLog("j is: \n");
                mxLog("%f", j);
                cout<< "p_e is: \n"  << p_e << endl;
                cout<< "ob_e is: \n"  << ob_e << endl;
            }
        }
        if (verbose >= 3){
            cout<< "yg_e is: \n"  << yg_e << endl;
        }
    } // end while (minit < maxit){
    
    yg_rec(0, 1) = yg_e.squaredNorm();
    if(yg_rec(0, 0) / yg_rec(0, 1) > 1000)  flag_NormgZ = 1;
    
    minr_rec = minit;

    p_e = p_e.cwiseProduct(vscale_e.block(0, neq+1, 1, nc+np-neq));
    // I need vscale, p, y, hessv
    if (nc > 0){
        y_e *= vscale_e(0);
        y_e = y_e.cwiseQuotient(vscale_e.block(0, 1, 1, nc));
    }
    
    // hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*%
    //                                t(vscale[ (neq + 2):(nc + np + 1) ]) )
    
    Eigen::MatrixXd transposePart;
    transposePart = vscale_e.block(0, neq+1, 1, nc+np-neq).transpose() * vscale_e.block(0, neq+1, 1, nc+np-neq);
    hessv_e = hessv_e.cwiseQuotient(transposePart);
    hessv_e = hessv_e * vscale_e(0);
    
    if (verbose >= 1 && reduce > tol) {
        mxLog("m3 solnp Rf_error message being reported.");
    }

    resP = p_e;
    resY = y_e.block(0, 0, 1, yyRows).transpose();
    resHessv = hessv_e;
    resLambda = lambdaValue;
    resGrad = g_e;
    
    if (verbose >= 3){
        mxLog("------------------------RETURNING FROM SUBNP------------------------");
        mxLog("p information: ");
        cout<< resP << endl;
        mxLog("y information: ");
        cout<< resY << endl;
        mxLog("hessv information: ");
        cout<< resHessv << endl;
        mxLog("lambda information: ");
        mxLog("%.20f", resLambda);
        mxLog("minit information: ");
        mxLog("%d", minit);
        mxLog("------------------------END RETURN FROM SUBNP------------------------");
    }
} // end subnp
