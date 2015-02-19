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

    int flag, flag_NormgZ, flag_step, minr_rec;
    Eigen::VectorXd LB_e;
    Eigen::VectorXd UB_e;
    Matrix resP;
    double resLambda;
    Matrix resHessv;
    Matrix resY;
	Matrix sx_Matrix; // search direction
    int mode;
	GradientOptimizerContext &fit;

	CSOLNP(GradientOptimizerContext &_fit) : fit(_fit) {};
	~CSOLNP() { freeMatrices(); };
    
	void solnp(double *pars, int verbose);

    Matrix subnp(Matrix pars, Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,
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
    flag_NormgZ = 0; flag_step = 0; minr_rec = 0;
    
    double funv;
    double resultForTT;
    double solnp_nfn = 0;
    
    //time_t sec;
    //sec = time (NULL);
    
    int maxit_trace = 0;
    
    Matrix grad = fill(LB_e.size(), 1, (double)0.0);
    //free(matrices.front().t);
    Matrix pb_cont;
    Matrix difference1, difference2, tmpv, testMin, firstCopied, subnp_ctrl, subsetMat, temp2, temp1, temp, funv_mat, tempdf, firstPart, copied, subsetOne, subsetTwo, subsetThree, diff1, diff2, copyValues, diff, llist, tempTTVals, searchD;
    
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
    Matrix pb;
    Eigen::MatrixXd pb_e;
    
    if(nineq) {
        pb = new_matrix(2, nineq + np);
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
        Eigen::Map< Eigen::MatrixXd >(pb.t, pbJoined.rows(), pbJoined.cols()) = pbJoined;
        
    } else {
        pb = fill(2, np, (double)0.0);
        pb_e.setZero(np, 2);
        pb_e.col(0) = LB_e;
        pb_e.col(1) = UB_e;
        Eigen::Map< Eigen::MatrixXd >(pb.t, pb_e.rows(), pb_e.cols()) = pb_e;
    }
    
    double rho   = fit.ControlRho;
    int maxit = fit.ControlMajorLimit;
    int minit = fit.ControlMinorLimit;
    double delta = fit.ControlFuncPrecision;
    double tol   = fit.ControlTolerance;
    
    int tc = nineq + neq;
    
    double j = funv;
    Matrix jh = fill(1, 1, funv);
    Matrix tt = fill(1, 3, (double)0.0);
    Eigen::VectorXd tt_e(3); tt_e.setZero();
    
    Matrix constraint;
    Eigen::MatrixXd constraint_e;
    
    fit.solEqBFun();
    Matrix eqv(fit.equality);
    Eigen::RowVectorXd eqv_e(neq);
    eqv_e = fit.equality;
    
    fit.myineqFun();
    Matrix ineqv(fit.inequality);
    Eigen::RowVectorXd ineqv_e(nineq);
    ineqv_e= fit.inequality;
    
    Matrix lambda;
    
    if (tc > 0){
        lambda = fill(1, tc, (double)0.0);
        
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
        lambda = fill(1, 1, (double)0.0);
    }
    
    constraint = new_matrix(constraint_e.cols(), constraint_e.rows());
    Eigen::Map< Eigen::MatrixXd > (constraint.t, constraint_e.rows(), constraint_e.cols()) = constraint_e;
    
    Eigen::Map< Eigen::VectorXd > (tt.t, tt_e.size()) = tt_e;
    Matrix tempv;
    Matrix p;
    Eigen::RowVectorXd p_e;
    
    if (nineq){
        p_e.resize(1, ineqx0_e.size() + pars.size());
        p_e << ineqx0_e, pars;
    }
    else{
        p_e = pars;
    }
    
    p = new_matrix(p_e.cols(), p_e.rows());
    Eigen::Map< Eigen::RowVectorXd > (p.t, p_e.size()) = p_e;
    
    Matrix hessv = diag(fill((np+nineq), 1, (double)1.0));
    Eigen::Map< Eigen::MatrixXd > hessv_e(hessv.t, hessv.rows, hessv.cols);
    hessv_e.setIdentity();
    double mu = np;
    
    int solnp_iter = 0;
    
    Eigen::MatrixXd ob_e;
    Matrix funvMatrix = fill(1, 1, funv);
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
    
    Matrix ob = new_matrix(ob_e.cols(), ob_e.rows());
    Eigen::Map< Eigen::MatrixXd > (ob.t, ob_e.rows(), ob_e.cols()) = ob_e;
    
    Matrix vscale;
    
    while(solnp_iter < maxit){
        solnp_iter = solnp_iter + 1;
        Eigen::Array<double, 4, 1> subnp_ctrl;
        subnp_ctrl[0] = rho;
        subnp_ctrl[1] = minit;
        subnp_ctrl[2] = delta;
        subnp_ctrl[3] = tol;
        
        if ( ind[indHasEq] > 0){
            Matrix subsetMat = subset(ob, 0, 1, neq);
            double max = matrixMaxAbs(subsetMat);
            
            Matrix temp2 = fill(neq, 1, max);
            Matrix temp1 = fill(1, 1, M(ob, 0, 0));
            vscale = copy(temp1, temp2);
        }
        else{
            vscale = fill(1, 1, (double)1.0);
        }
	vscale = copy(vscale, fill(p.cols, 1, (double)1.0));

        minMaxAbs(vscale, tol);
        
        if (verbose >= 1){
            mxLog("------------------------CALLING SUBNP------------------------");
            mxLog("p information: ");
            for (i = 0; i < p.cols; i++) mxLog("%f",p.t[i]);
            mxLog("lambda information: ");
            for (i = 0; i < lambda.cols; i++) mxLog("%f",lambda.t[i]);
            mxLog("ob information: ");
            for (i = 0; i < ob.cols; i++) mxLog("%f",ob.t[i]);
            mxLog("hessv information: ");
            for (i = 0; i < hessv.cols*hessv.cols; i++) mxLog("%f",hessv.t[i]);
            mxLog("mu information: ");
            mxLog("%2f", mu);
            mxLog("vscale information: ");
            for (i = 0; i < vscale.cols; i++) mxLog("%f",vscale.t[i]);
            mxLog("subnp_ctrl information: ");
            for (i = 0; i < subnp_ctrl.size(); i++) mxLog("%f",subnp_ctrl[i]);
            mxLog("------------------------END CALLING SUBNP------------------------");
        }
        
        if (mode == -1)
        {
		fit.informOut = 0;
		memcpy(pars.data(), p.t, pars.size() * sizeof(double));
		return;
        }
        
        if (sx_Matrix.t == NULL) sx_Matrix = fill(p.cols, p.rows, (double)0.0);
        
        grad = subnp(p, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        
        p = duplicateIt(resP);
        if (flag == 1)
        {
	    mode = 0;
            funv = fit.solFun(p.t, &mode);
            funvMatrix = fill(1, 1, funv);
            fit.solEqBFun();
            if ( nineq) {
		    fit.myineqFun();
		    Matrix ineqv(fit.inequality);
                if(eqv.cols)
                {
                    Matrix firstCopied = copy(funvMatrix, eqv);
                    ob = copy(firstCopied, ineqv);
                }
                else{
                    ob = copy(funvMatrix, ineqv);
                }
            }
            else if (eqv.cols){
                ob = copy(funvMatrix, eqv);
            }
            else ob = duplicateIt(funvMatrix);
            
            if ( ind[indHasEq] > 0){
                Matrix subsetMat = subset(ob, 0, 1, neq);
                double max = matrixMaxAbs(subsetMat);
                
                Matrix temp2 = fill(neq, 1, max);
                Matrix temp1 = fill(1, 1, M(ob, 0, 0));
                vscale = copy(temp1, temp2);
                
            }
            else{
                vscale = fill(1, 1, (double)1.0);
            }
	    vscale = copy(vscale, fill(p.cols, 1, (double)1.0));
            minMaxAbs(vscale, tol);
            lambda = duplicateIt(resY);
            hessv = duplicateIt(resHessv);
            mu = resLambda;
            grad = subnp(p, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        }
        
        lambda = duplicateIt(resY);
        
        hessv = duplicateIt(resHessv);
        
        mu = resLambda;
        
        
        Matrix temp = subset(p, 0, nineq, (nineq+np-1));
        mode = 1;
        funv = fit.solFun(temp.t, &mode);
        if (mode == -1)
        {
		fit.informOut = 0;
		memcpy(pars.data(), p.t, pars.size() * sizeof(double));
		return;
        }
        
        solnp_nfn = solnp_nfn + 1;
        
        //Matrix funv_mat = fill(1, 1, funv);
        //Matrix tempdf = copy(temp, funv_mat);
        fit.solEqBFun();
        
        Matrix firstPart, copied;
        if (nineq){
		fit.myineqFun();
		Matrix ineqv(fit.inequality);
            if(eqv.cols){
                copied = copy(fill(1, 1, funv), eqv);
                ob = copy(copied, ineqv);
                
            }
            else{
                ob = copy(fill(1, 1, funv), ineqv);
            }
        }
        else if (eqv.cols){
            ob = copy(fill(1, 1, funv), eqv);
        }
        else ob = fill(1, 1, funv);
        
        if (verbose >= 1){
            mxLog("j2 in while: \n");
            mxLog("%.20f", j);
            mxLog("M(ob2, 0, 0) \n");
            mxLog("%.20f", M(ob, 0, 0));
        }
        
        resultForTT = (j - M(ob, 0, 0)) / max(fabs(M(ob, 0, 0)), 1.0);
        M(tt, 0, 0) = resultForTT;
        if (verbose >= 1){
            mxLog("resultForTT \n");
            mxLog("%.20f", resultForTT);
        }
        j = M(ob, 0, 0);
        
        if (tc > 0){
            // constraint = ob[ 2:(tc + 1) ]
            constraint = subset(ob, 0, 1, tc);
            
            if ( ind[indHasIneq] > 0.5){
                //tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ], pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
                Matrix subsetOne = subset(constraint, 0, neq, tc-1);
                Matrix subsetTwo = subset(getColumn(pb, 0), 0, 0, nineq-1);
                Matrix subsetThree = subset(getColumn(pb, 1), 0, 0, nineq-1);
                subtractEigen(subsetOne, subsetTwo);
                subtractEigen(subsetThree, subsetOne);
                Matrix tempv = fill(nineq, 2, (double)0.0);
                setRowInplace(tempv, 0, subsetOne);
                setRowInplace(tempv, 1, subsetThree);
                
                if (findMin(tempv) > 0){
                    Matrix copyValues = subset(constraint, 0, neq, tc-1);
                    copyIntoInplace(p, copyValues, 0, 0, nineq-1);
                }
                Matrix diff = subset(constraint, 0, neq, tc-1);
                subtractEigen(diff, subset(p, 0, 0, nineq-1));
                
                copyIntoInplace(constraint, diff, 0, neq, tc-1);
            } // end if (ind[0][3] > 0.5){
            
            M(tt, 0, 2) = vnorm(constraint);
            
            
            if ( M(tt, 0, 2) < (10 *tol)){
                rho =0;
                mu = min(mu, tol);
            }
            
            if ( M(tt, 0, 2) < (5 * M(tt, 0, 1))){
                rho = rho/5;
            }
            
            if ( M(tt, 0, 2) > (10 * M(tt, 0, 1))){
                rho = 5 * max(rho, sqrt(tol));
            }
            
            llist = fill(2, 1, (double)0.0);
            
            M(llist, 0, 0) = tol + M(tt, 0, 0);
            M(llist, 1, 0) = M(tt, 0, 1) - M(tt, 0, 2);
            
            
            if (findMax(llist) <= 0){
                //hessv = diag( diag ( hessv ) )
                lambda = fill(1, 1, (double)0.0);
                hessv = diag(diag2(hessv));
            }
            
            M(tt, 0, 1) = M(tt, 0, 2);
            
        } // end if (tc > 0){
        
        if (verbose >= 3){
            mxLog("tt is \n");
            for (i = 0; i < tt.cols; i++) mxLog("%f",tt.t[i]);
        }
        double vnormValue;
        Matrix tempTTVals = fill(2, 1, (double)0.0);
        M(tempTTVals, 0, 0) = M(tt, 0, 0);
        M(tempTTVals, 1, 0) = M(tt, 0, 1);
        
        vnormValue = vnorm(tempTTVals);
        
        if (vnormValue <= tol){
            maxit_trace = maxit;
            maxit = solnp_iter;
        }
        
        if (verbose >= 3)
        {
            mxLog("vnormValue in while \n");
            mxLog("%.20f", vnormValue);
        }
        jh = copy(jh, fill(1, 1, j));
        
    } // end while(solnp_iter < maxit){
    
    p = subset(p, 0, nineq, (nineq + np -1));

    {
        Matrix tempTTVals = fill(2, 1, (double) 0.0);
        M(tempTTVals, 0, 0) = M(tt, 0, 0);
        M(tempTTVals, 1, 0) = M(tt, 0, 1);
        double vnormValue = vnorm(tempTTVals);
        if (verbose >= 3){
            mxLog("sx_Matrix (search direction) is: \n");
            for (i = 0; i < sx_Matrix.cols; i++) mxLog("%f",sx_Matrix.t[i]);
        }
        
        if (verbose >= 1) {
		mxLog("vnormValue %.20f, flag_NormgZ=%d, minr_rec=%d, flag_step=%d",
		      vnormValue, flag_NormgZ, minr_rec, flag_step);
	}
        if (vnormValue <= tol && flag_NormgZ == 1 && minr_rec == 1 && flag_step == 1){
		double iterateConverge = delta * pow(vnorm(sx_Matrix),(double)2.0);
		double iterateConvergeCond = sqrt(tol) * ((double)1.0 + pow(vnorm(p), (double)2.0));
		if (verbose >= 1) {
			mxLog("vnorm(sx_Matrix) is %.20f, iterateConverge is %.20f, iterateConvergeCond is: %.20f",
			      vnorm(sx_Matrix), iterateConverge, iterateConvergeCond);
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
    
    memcpy(pars.data(), p.t, pars.size() * sizeof(double));
    fit.gradOut.resize(grad.cols);
    memcpy(fit.gradOut.data(), grad.t, fit.gradOut.size() * sizeof(double));
    fit.hessOut.resize(hessv.rows, hessv.cols);
    memcpy(fit.hessOut.data(), hessv.t, fit.hessOut.size() * sizeof(double));
}

Matrix CSOLNP::subnp(Matrix pars, Matrix yy,  Matrix ob,  Matrix hessv,
		     double lambda,  Matrix vscale, const Eigen::Array<double, 4, 1> &ctrl, int verbose)
{
    if (verbose >= 3)
    {
        mxLog("pars in subnp is: \n");
        for (int ilog = 0; ilog < pars.cols; ilog++) mxLog("%f",pars.t[ilog]);
    }
    int yyRows = yy.rows;
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
    Matrix t_sol;
    
    int nc = neq + nineq;
    int npic = np + nineq;
    
    Matrix p0 = duplicateIt(pars);
    
    if (verbose >= 3)
    {
        mxLog("p0 p0 is: \n");
        for (int i = 0; i < p0.cols; i++) mxLog("%f",p0.t[i]);
    }
    
    Matrix pb;
    Eigen::MatrixXd pb_e;
    /*Eigen::Map< Eigen::VectorXd > LB_e(LB.t, LB.cols);
    Eigen::Map< Eigen::VectorXd > UB_e(UB.t, UB.cols);*/
    
    if(nineq) {
        pb = new_matrix(2, nineq + np);
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
        Eigen::Map< Eigen::MatrixXd >(pb.t, pbJoined.rows(), pbJoined.cols()) = pbJoined;
        
    } else {
        pb = fill(2, np, (double)0.0);
        pb_e.setZero(np, 2);
        pb_e.col(0) = LB_e;
        pb_e.col(1) = UB_e;
        Eigen::Map< Eigen::MatrixXd >(pb.t, pb_e.rows(), pb_e.cols()) = pb_e;
    }
    
    if (verbose >= 3){
        cout<< "pb_e is: \n" << pb_e<< endl;
    }
    
    Eigen::Array<double, 3, 1> sob;
    sob.setZero();
    
    //Matrix yyMatrix = duplicateIt(yy);
    
    //Eigen::Map< Eigen::MatrixXd > ob_e(ob.t, ob.rows, ob.cols);
    Eigen::Map< Eigen::RowVectorXd > vscale_e(vscale.t, vscale.cols);
    //ob_e = ob_e.cwiseQuotient(vscale_e.block(0, 0, 1, nc + 1));
    divideEigen(ob, subset(vscale, 0, 0, nc));
    Eigen::Map< Eigen::RowVectorXd > p0_e(p0.t, p0.cols);
    p0_e = p0_e.cwiseQuotient(vscale_e.block(0, neq + 1, 1, nc + np - neq));
    
    if (verbose >= 3){
        cout<< "p0_e: \n"<< p0_e<< endl;
        cout<< "vscale_e: \n"<< vscale_e<< endl;
    }
    
    int mm = 0;
    {
        mm=npic;
        Eigen::MatrixXd pbCopied;
        pbCopied.setZero(pb.rows, pb.cols);
        pbCopied.col(0) = vscale_e.block(0, neq + 1, 1, mm).transpose();
        pbCopied.col(1) = vscale_e.block(0, neq + 1, 1, mm).transpose();
        pb_e = pb_e.cwiseQuotient(pbCopied);
    }
    
    if (verbose >= 3){
        mxLog("pb is: \n");
        for (int i = 0; i < pb.cols; i++) mxLog("%f",pb.t[i]);
    }
    
    Eigen::Map < Eigen::MatrixXd > yy_e(yy.t, yy.rows, yy.cols);
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
    
    Eigen::Map < Eigen::MatrixXd > hessv_e(hessv.t, hessv.rows, hessv.cols);
    Eigen::MatrixXd result_e;
    result_e = vscale_e.block(0, neq + 1, 1, nc + np - neq).transpose() * vscale_e.block(0, neq + 1, 1, nc + np - neq);
    hessv_e = hessv_e.cwiseProduct(result_e);
    hessv_e = hessv_e / vscale_e[0];
    Eigen::Map< Eigen::MatrixXd > ob_e(ob.t, ob.rows, ob.cols);
    j = ob_e(0, 0);
    if (verbose >= 3){
        mxLog("j j is: \n");
        mxLog("%f", j);
    }
    Eigen::MatrixXd a_e;
    Matrix a;
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
            a = new_matrix(np + nineq, nineq);
            Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
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
            a = new_matrix(a_e.cols(), a_e.rows());
            Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
            //a = transpose(copy(transpose(firstHalf), transpose(secondHalf)));
        }
    }	// end 	if(ind[0][3] > 0){
    
    if ( (ind[indHasEq] > 0) && ind[indHasIneq] <= 0 ){
        a_e.resize(neq, np);
        a_e.setZero();
        a = new_matrix(a_e.cols(), a_e.rows());
        Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
        //a = fill(np, neq, (double)0.0);
    }
    if (ind[indHasEq]<= 0 && (ind[indHasIneq] <= 0)){
        a_e.resize(1, np);
        a_e.setZero();
        a = new_matrix(a_e.cols(), a_e.rows());
        Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
        //a = fill(np, 1, (double)0.0);
    }

    Matrix g = fill(npic, 1, (double)0.0);
    Eigen::Map< Eigen::RowVectorXd > g_e(g.t, g.cols);
    Matrix p = subset(p0, 0, 0, (npic-1));
    Eigen::Map< Eigen::RowVectorXd > p_e(p.t, p.cols);
    
    Matrix dx;
    Matrix b;
    Eigen::MatrixXd b_e;
    double funv;
    Matrix eqv(fit.equality);
    Matrix ineqv(fit.inequality);
    Matrix constraint;
    Matrix gap;
    
    int solnp_nfn = 0;
    double go, reduce = 1e-300;
    int minit;
    double lambdaValue = lambda;
    
    Eigen::MatrixXd constraint_e(1, nc);
    Eigen::MatrixXd y_e;
    
    if (nc > 0) {
        constraint_e = ob_e.block(0, 1, 1, nc);
        Matrix constraint = new_matrix(nc, 1);
        Eigen::Map< Eigen::MatrixXd > (constraint.t, constraint_e.rows(), constraint_e.cols()) = constraint_e;
        
        for (int i=0; i<np; i++){
            int index = nineq + i;
            M(p0, index, 0) = M(p0, index, 0) + delta;
            Eigen::MatrixXd tmpv_e;
            tmpv_e = p0_e.block(0, nineq, 1, npic - nineq);
            tmpv_e = tmpv_e.array() * vscale_e.block(0, nc+1, 1, np).array();
            Matrix tmpv = new_matrix(1, npic - nineq);
            Eigen::Map< Eigen::MatrixXd > (tmpv.t, tmpv_e.rows(), tmpv_e.cols()) = tmpv_e;
            if (verbose >= 2){
                mxLog("7th call is \n");
            }
            funv = fit.solFun(tmpv.t, &mode);
            
            fit.solEqBFun();
            fit.myineqFun();
            
            solnp_nfn = solnp_nfn + 1;
            
            Eigen::MatrixXd firstPart_e;
            
            Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
            Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
            Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
            
            if (nineq){
                if(eqv.cols)
                {
                    firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                    firstPart_e << funv_e, eqv_e, ineqv_e;
                }
                else{
                    firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                    firstPart_e << funv_e, ineqv_e;
                }
            }
            else if (eqv.cols){
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
            
            M(g, index, 0) = (ob_e(0, 0)-j) / delta;
            
            if (verbose >= 3){
                mxLog("g is: \n");
                for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
                mxLog("a is: \n");
                for (int ilog = 0; ilog < a.cols * a.rows; ilog++) mxLog("%f",a.t[ilog]);
                
            }
            
            a_e.col(index) = (ob_e.block(0, 1, 1, nc) - constraint_e).transpose() / delta;
            Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
            M(p0, index, 0) = M(p0, index, 0) - delta;
        } // end for (int i=0; i<np, i++){
        
        if (mode == -1)
        {
            funv = 1e24;
            mode = 0;
        }
        
        if(ind[indHasIneq] > 0){
            //constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
            Matrix firstPart, secondPart;
            constraint_e.block(0, neq, 1, nineq) = (constraint_e.block(0, neq, 1, nineq) - p0_e.block(0, 0, 1, nineq)).block(0, 0, 1, nineq);
            Eigen::Map< Eigen::MatrixXd > (constraint.t, constraint_e.rows(), constraint_e.cols()) = constraint_e;
        }
        
        if (false && solvecond(a) > 1/DBL_EPSILON) { // this can't be the cheapest way to check TODO
            Rf_error("Redundant constraints were found. Poor intermediate results may result. "
                     "Remove redundant constraints and re-OPTIMIZE.");
        }
        
        b_e = (a_e * p0_e.transpose()).transpose();
        //  b [nc,1]
        b_e -= constraint_e;
        b = new_matrix(b_e.cols(), b_e.rows());
        Eigen::Map< Eigen::MatrixXd > (b.t, b_e.rows(), b_e.cols()) = b_e;
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
            p0 = copy(p0, fill(1, 1, (double)1.0));
            new (&p0_e) Eigen::Map<Eigen::RowVectorXd>(p0.t, p0.rows, p0.cols);
            p0_e.resize(p0_e_copy.rows(), p0_e_copy.cols() + onesMatrix_e.cols());
            p0_e << p0_e_copy, onesMatrix_e;
            constraint_e *= (-1.0);
            Eigen::Map< Eigen::MatrixXd > (constraint.t, constraint_e.rows(), constraint_e.cols()) = constraint_e;
            Eigen::MatrixXd a_e_copy = a_e;
            a_e.resize(a_e.rows(), a_e.cols() + constraint_e.transpose().cols());
            a_e << a_e_copy, constraint_e.transpose();
            a = new_matrix(a_e.cols(), a_e.rows());
            Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
            Eigen::MatrixXd firstMatrix_e(1, npic);
            firstMatrix_e.setZero();
            Eigen::MatrixXd cx_e(firstMatrix_e.rows(), firstMatrix_e.cols() + onesMatrix_e.cols());
            cx_e << firstMatrix_e, onesMatrix_e;
            Matrix cx = new_matrix(cx_e.cols(), cx_e.rows());
            Eigen::Map< Eigen::MatrixXd > (cx.t, cx_e.rows(), cx_e.cols()) = cx_e;
            Eigen::MatrixXd dx_e(npic + 1, 1);
            dx_e.setOnes();
            dx = new_matrix(dx_e.cols(), dx_e.rows());
            Eigen::Map< Eigen::MatrixXd > (dx.t, dx_e.rows(), dx_e.cols()) = dx_e;
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
                gap = new_matrix(2, mm);
                Eigen::Map< Eigen::MatrixXd > (gap.t, gap_e.rows(), gap_e.cols()) = gap_e;
                dx_e.transpose().block(0, 0, 1, mm) = gap_e.col(0).transpose().block(0, 0, 1, mm);
                dx_e(npic_int, 0) = p0_e(0, npic_int);
                Eigen::MatrixXd argum1_e;
                argum1_e = a_e * dx_e.asDiagonal();
                argum1_e.transposeInPlace();
                Eigen::MatrixXd argum2_e;
                argum2_e = cx_e.asDiagonal() * dx_e;
                Matrix argum1 = new_matrix(argum1_e.cols(), argum1_e.rows());
                Eigen::Map< Eigen::MatrixXd > (argum1.t, argum1_e.rows(), argum1_e.cols()) = argum1_e;
                Matrix argum2 = new_matrix(argum2_e.cols(), argum2_e.rows());
                Eigen::Map< Eigen::MatrixXd > (argum2.t, argum2_e.rows(), argum2_e.cols()) = argum2_e;
                Eigen::Map< Eigen::MatrixXd > (dx.t, dx_e.rows(), dx_e.cols()) = dx_e;
                y_e = argum1_e.colPivHouseholderQr().solve(argum2_e);
                Matrix y = new_matrix(y_e.cols(), y_e.rows());
                Eigen::Map< Eigen::MatrixXd > (y.t, y_e.rows(), y_e.cols()) = y_e;
                Eigen::MatrixXd cx_e_r;
                cx_e_r = cx_e.transpose() - (a_e.transpose() * y_e);
                dx_e = (cx_e_r.asDiagonal() * dx_e).asDiagonal() * dx_e;
                Eigen::Map< Eigen::MatrixXd > (dx.t, dx_e.rows(), dx_e.cols()) = dx_e;
                Eigen::MatrixXd v_e = dx_e.transpose();
                Matrix v = new_matrix(v_e.cols(), v_e.rows());
                Eigen::Map< Eigen::MatrixXd > (v.t, v_e.rows(), v_e.cols()) = v_e;
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
            Matrix aMatrix = fill(npic, nc, (double)0.0);
            Eigen::MatrixXd a_e_c(nc, npic);
            
            for (h = 0; h<a_e.rows(); h++)
            {
                a_e_c.row(h) = a_e.row(h).block(0, 0, 1, npic);
            }
            a_e.resize(a_e_c.rows(), a_e_c.cols());
            a_e = a_e_c;
            a = duplicateIt(aMatrix);
            Eigen::Map< Eigen::MatrixXd > (a.t, a_e.rows(), a_e.cols()) = a_e;
            b_e = (a_e * p0_e.block(0, 0, 1, npic).transpose()).transpose();
            Eigen::Map< Eigen::MatrixXd > (b.t, b_e.rows(), b_e.cols()) = b_e;
        }// end if(M(alp, 0, 0) <= 0)
    } // end if (nc > 0){
    
    p_e = p0_e.block(0, 0, 1, npic);
    
    if (nc == 0){
        y_e.resize(1,1);
        y_e(0, 0) = 0;
    }

    if (verbose >= 3){
        mxLog("p is: \n");
        for (int i = 0; i < p.cols; i++) mxLog("%f",p.t[i]);
    }
    
    if (ch > 0){
        
        Eigen::MatrixXd tmpv_e;
        tmpv_e = p_e.block(0, nineq, 1, npic-nineq);
        tmpv_e = tmpv_e.array() * vscale_e.block(0, nc+1, 1, np).array();
        Matrix tmpv = new_matrix(tmpv_e.cols(), tmpv_e.rows());
        Eigen::Map< Eigen::MatrixXd > (tmpv.t, tmpv_e.rows(), tmpv_e.cols()) = tmpv_e;
        funv = fit.solFun(tmpv.t, &mode);
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
        Matrix firstPart, secondPart, firstPartt;
        Eigen::MatrixXd firstPart_e;
        Eigen::RowVectorXd funv_e(1); funv_e[0] = funv;
        Eigen::RowVectorXd eqv_e(neq); eqv_e = fit.equality;
        Eigen::RowVectorXd ineqv_e(nineq); ineqv_e= fit.inequality;
        
        if (nineq){
            if (eqv.cols){
                firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                firstPart_e << funv_e, eqv_e, ineqv_e;
            }
            else{
                firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                firstPart_e << funv_e, ineqv_e;
            }
        }
        else if (eqv.cols){
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
    Matrix yg = fill(npic, 1, (double)0.0);
    Eigen::MatrixXd yg_e;
    Eigen::MatrixXd yg_rec(1, 2);
    Matrix sx;
    Eigen::MatrixXd sx_e;
    sx_e.setZero(p.rows, p.cols);
    Matrix sc = fill(2, 1, (double)0.0);
    Matrix cz;
    Matrix czSolution;
    Matrix u, t_u;
    Matrix vscale_t;
    Matrix t1;
    Matrix firstPart, secondPart, firstPartt, obm_t;
    Matrix result, result1;
    Matrix first_part;
    Matrix p_t, atime, t_atime;
    Matrix temp_t, t_yy, row_tt_y, row_temp_t, t2, t3, t4, t5, t6;
    Matrix m_sc2, m_sc, t_sx, t_yg, sxMatrix, ygMatrix, sx2, yg2;
    Matrix gap1, gap2, gap_c, col_pb, col_pb2, t7, t8, t9, res;
    Matrix t11, t10;
    Matrix dx_c, dxDiag, t12, t13, t14, t15;
    Matrix aTranspose, t_cz, firstMatrix, solution;
    Matrix toSubtract, y, uu;
    Matrix listPartOne, listPartTwo, llist, t16, t17;
    Matrix p0_1;
    Matrix ptt, t18, p0_1_t, ptt2, t19, pttCol, t20, ob3;
    Matrix partOne, partTwo, tempPttCol;
    Matrix yyTerm, firstp, t21, t22, t23, t24, t25, t26, t27;
    Matrix p_copy, p0_copy, pttColOne, t29, t30, t31, t32, t33, t34, temp, tempCol;
    Matrix input, rhs;
    Matrix tmpv;
    
    
    Eigen::MatrixXd obm_e;
    
    while (minit < maxit){
        minit = minit + 1;
        if (ch > 0){
            
            for (int i=0; i<np; i++){
                int index = nineq + i;
                p_e[index] = p_e[index] + delta;
                Eigen::MatrixXd tmpv_e = p_e.block(0, nineq, 1, npic - nineq).array() * vscale_e.block(0, nc+1, 1, np).array();
                Matrix tmpv = new_matrix(tmpv_e.cols(), tmpv_e.rows());
                Eigen::Map< Eigen::MatrixXd > (tmpv.t, tmpv_e.rows(), tmpv_e.cols()) = tmpv_e;
                if (verbose >= 3){
                    mxLog("9th call is \n");
                }
                mode = 0;
                funv = fit.solFun(tmpv.t, &mode);
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
                    if (eqv.cols){
                        firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                        firstPart_e << funv_e, eqv_e, ineqv_e;
                    }
                    else{
                        firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                        firstPart_e << funv_e, ineqv_e;
                    }
                }
                else if (eqv.cols){
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
                
                M(g, index, 0) = (obm - j)/delta;
                p_e[index] = p_e[index] - delta;
                
                if (verbose >= 3){
                    mxLog("g is: \n");
                    for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
                    mxLog("p is: \n");
                    for (int ilog = 0; ilog < p.cols; ilog++) mxLog("%f",p.t[ilog]);
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
            for (int ilog = 0; ilog < yg.cols; ilog++) mxLog("%f",yg.t[ilog]);
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
            dx = new_matrix(dx_e.cols(), dx_e.rows());
            Eigen::Map< Eigen::MatrixXd > (dx.t, dx_e.rows(), dx_e.cols()) = dx_e;
            if (verbose >= 3){
                mxLog("dx is: \n");
                for (int ilog = 0; ilog < dx.cols * dx.rows; ilog++) mxLog("%.50f",dx.t[ilog]);
            }
        }
        
        go = -1;
        lambdaValue = lambdaValue/10.0;
        
        Matrix yMatrix;
        
        if (verbose >= 3){
            mxLog("lambdaValue is: \n");
            mxLog("%.20f", lambdaValue);
            mxLog("hessv is: \n");
            for (int ilog = 0; ilog < hessv.cols; ilog++) mxLog("%f",hessv.t[ilog]);
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
            cz = new_matrix(cz_e.cols(), cz_e.rows());
            Eigen::Map< Eigen::MatrixXd > (cz.t, cz_chol.rows(), cz_chol.cols()) = cz_chol;
            if (verbose >= 3){
                cout<< "cz_chol is: \n"  << cz_chol << endl;
            }
            
            if (!R_FINITE((cz_e.maxCoeff())))
            {
                mxLog("here in findMax");
                flag = 1;
                p_e = p_e.cwiseProduct(vscale_e.block(0, neq+1, 1, nc+np-neq));
                if (nc > 0){ y = fill(1, 1, (double)0.0);}
                hessv_e = hessv_e.cwiseQuotient(vscale_e.block(0, neq+1, 1, nc+np-neq) * vscale_e.block(0, neq+1, 1, nc+np-neq).transpose()) *vscale_e(0);
                resP = duplicateIt(p);
                resY = duplicateIt(y);
                resHessv = duplicateIt(hessv);
                resLambda = lambda;
                return g;
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
            
            Eigen::RowVectorXd p0_e_copy;
            p0_e_copy = u_e.block(0, 0, 1, npic) + p_e;
            p0_copy = new_matrix(p0_e_copy.cols(), 1);
            Eigen::Map< Eigen::RowVectorXd > (p0_copy.t, p0_e_copy.cols()) = p0_e_copy;
            new (&p0_e) Eigen::Map<Eigen::RowVectorXd>(p0_copy.t, p0_copy.cols);
            
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
        Matrix ob1 = new_matrix(ob1_e.cols(), ob1_e.rows());
        Eigen::Map< Eigen::MatrixXd > (ob1.t, ob1_e.rows(), ob1_e.cols()) = ob1_e;
        Matrix ob2 = new_matrix(ob2_e.cols(), ob2_e.rows());
        Eigen::Map< Eigen::MatrixXd > (ob2.t, ob2_e.rows(), ob2_e.cols()) = ob2_e;
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
        Matrix tmpv = new_matrix(tmpv_e.cols(), tmpv_e.rows());
        Eigen::Map< Eigen::MatrixXd > (tmpv.t, tmpv_e.rows(), tmpv_e.cols()) = tmpv_e;
        Matrix ptt2 = new_matrix(ptt_e.cols(), ptt_e.rows());
        Eigen::Map< Eigen::MatrixXd > (ptt2.t, ptt_e.rows(), ptt_e.cols()) = ptt_e;
        
        if (verbose >= 3){
            cout<< "tmpv_e is: \n"  << tmpv_e << endl;
        }
        
        mode = 1;
        funv = fit.solFun(tmpv.t, &mode);
        if (verbose >= 3){
            mxLog("hessv is: \n");
            for (int ilog = 0; ilog < hessv.cols; ilog++) mxLog("%f",hessv.t[ilog]);
            
            mxLog("g is: \n");
            for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
            
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
            if(eqv.cols)
            {
                firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                firstPart_e << funv_e, eqv_e, ineqv_e;
            }
            else{
                firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                firstPart_e << funv_e, ineqv_e;
            }
        }
        else if (eqv.cols){
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
        ob3 = new_matrix(ob3_e.cols(), ob3_e.rows());
        Eigen::Map< Eigen::MatrixXd > (ob3.t, ob3_e.rows(), ob3_e.cols()) = ob3_e;
        
        go = 1;
        
        while(go > tol){
            alp[1] = (alp[0] + alp[2]) / 2.0;
            
            ptt_e.col(1) = (p_e * (1 - alp[1])) + p0_e * alp[1];
            Eigen::MatrixXd tmpv_e = ptt_e.col(1).transpose().block(0, nineq, 1, npic - nineq).cwiseProduct(vscale_e.block(0, nc+1, 1, np));
            Matrix tmpv = new_matrix(tmpv_e.cols(), tmpv_e.rows());
            Eigen::Map< Eigen::MatrixXd > (tmpv.t, tmpv_e.rows(), tmpv_e.cols()) = tmpv_e;
            
            if (verbose >= 3){
                cout<< "tmpv_e is: \n"  << tmpv_e << endl;
            }
            
            if (verbose >= 3){
                mxLog("11th call is \n");
            }
            
            mode = 0;
            funv = fit.solFun(tmpv.t, &mode);
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
                if(eqv.cols)
                {
                    firstPart_e.resize(1, funv_e.size()+ eqv_e.size()+ ineqv_e.size());
                    firstPart_e << funv_e, eqv_e, ineqv_e;
                }
                else{
                    firstPart_e.resize(1, funv_e.size() + ineqv_e.size());
                    firstPart_e << funv_e, ineqv_e;
                }
            }
            else if (eqv.cols){
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
            Eigen::MatrixXd ob2_e = firstPart_e;
            Matrix ob2 = new_matrix(ob2_e.cols(), ob2_e.rows());
            Eigen::Map< Eigen::MatrixXd > (ob2.t, ob2_e.rows(), ob2_e.cols()) = ob2_e;
            
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
        
        Eigen::Map< Eigen::MatrixXd > (ptt2.t, ptt_e.rows(), ptt_e.cols()) = ptt_e;
        Eigen::Map< Eigen::MatrixXd > (ob1.t, ob1_e.rows(), ob1_e.cols()) = ob1_e;
        Eigen::Map< Eigen::MatrixXd > (ob2.t, ob2_e.rows(), ob2_e.cols()) = ob2_e;
        Eigen::Map< Eigen::MatrixXd > (ob3.t, ob3_e.rows(), ob3_e.cols()) = ob3_e;

        if (verbose >= 3){
            mxLog("go is: \n");
            mxLog("%.16f", go);
        }

        sx = new_matrix(sx_e.cols(), sx_e.rows());
        Eigen::Map< Eigen::MatrixXd > (sx.t, sx_e.rows(), sx_e.cols()) = sx_e;
        sx_Matrix = sx;
        sx_e.resize(p_e.rows(), p_e.cols());
        sx_e = p_e;
        yg_e.resize(g_e.rows(), g_e.cols());
        yg_e = g_e;
        
        if (verbose >= 3){
            mxLog("sx is: \n");
            for (int ilog = 0; ilog < sx.cols; ilog++) mxLog("%f",sx.t[ilog]);
            mxLog("yg is: \n");
            for (int ilog = 0; ilog < yg.cols; ilog++) mxLog("%f",yg.t[ilog]);
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
    
    Matrix result2 = getColumn(ptt2, 1);
    subtractEigen(result2, getColumn(ptt2, 0));
    Matrix result3 = getColumn(ptt2, 1);
    subtractEigen(result3, getColumn(ptt2, 2));
    if (all(result2) || all(result3)) flag_step = 1;
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
    y = new_matrix(y_e.cols(), y_e.rows());
    Eigen::Map< Eigen::MatrixXd > (y.t, y_e.rows(), y_e.cols()) = y_e;
    
    resP = duplicateIt(p);
    resY = transpose(subset(y, 0, 0, (yyRows-1)));
    resHessv = duplicateIt(hessv);
    resLambda = lambdaValue;
    
    if (verbose >= 3){
        mxLog("------------------------RETURNING FROM SUBNP------------------------");
        mxLog("p information: ");
        for (int ilog = 0; ilog < resP.cols; ilog++) mxLog("%f",resP.t[ilog]);
        mxLog("y information: ");
        for (int ilog = 0; ilog < resY.cols; ilog++) mxLog("%f",resY.t[ilog]);
        mxLog("hessv information: ");
        for (int ilog = 0; ilog < resHessv.cols * resHessv.rows; ilog++) mxLog("%f",resHessv.t[ilog]);
        mxLog("lambda information: ");
        mxLog("%f", resLambda);
        mxLog("minit information: ");
        mxLog("%d", minit);
        mxLog("------------------------END RETURN FROM SUBNP------------------------");
    }
    
    return g;
    
} // end subnp
