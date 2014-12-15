/*
 *  Copyright 2007-2012 The OpenMx Project
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *   Unless required by applicable law or agreed to in writing, software
 *   distributed under the License is distributed on an "AS IS" BASIS,
 *   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#ifndef _OMX_CSOLNP_SPECIFIC_H
#define _OMX_CSOLNP_SPECIFIC_H

#include "omxMatrix.h"
#include "matrix.h"
#include <Eigen/Core>

struct CSOLNPFit {
	Eigen::VectorXd solIneqLB;
	Eigen::VectorXd solIneqUB;

	virtual double solFun(double *myPars, int* mode, int verbose) = 0;
	virtual Matrix solEqBFun(int verbose) = 0;
	virtual Matrix myineqFun(int verbose) = 0;
};

struct CSOLNP {

    int flag, flag_L, flag_U, index_flag_L, index_flag_U, flag_NormgZ, flag_step, minr_rec;
    Matrix LB;
    Matrix UB;
    Matrix resP;
    double resLambda;
    Matrix resHessv;
    Matrix resY;
    Matrix sx_Matrix;
    int mode_val;
    int* mode;
	CSOLNPFit &fit;

	CSOLNP(CSOLNPFit &_fit) : fit(_fit) {};
    
	template <typename ParType>
	Param_Obj solnp(Eigen::MatrixBase<ParType> &solPars, Matrix solLB, Matrix solUB, Matrix solctrl, bool debugToggle, int verbose);

    Matrix subnp(Matrix pars, Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,
		 const Eigen::Array<double, 5, 1> &ctrl, int verbose);

	enum indParam {
		indNumParam=0,
		indHasGradient,
		indHasHessian,
		indHasIneq,
		indIneqLength,
		indHasJacobianIneq,
		indHasEq,
		indEqLength,
		indHasJacobianEq,
		indHasBounds,
		indHasBoundsOrIneq,
		indVectorLength  // must be last
	};

	Eigen::Array<double, int(indVectorLength), 1> ind;
};

void omxInvokeCSOLNP(omxMatrix *fitMatrix, FitContext *fc, int *inform_out,
                     FreeVarGroup *freeVarGroup, int verbose, double *hessOut,
                     double tolerance);

void omxCSOLNPConfidenceIntervals(omxMatrix *fitMatrix, FitContext *fc, int verbose,
                                  double tolerance);

void CSOLNPOpt_majIter(const char *optionValue);

void CSOLNPOpt_minIter(const char *optionValue);

void CSOLNPOpt_FuncPrecision(const char *optionValue);

template <typename ParType>
Param_Obj CSOLNP::solnp(Eigen::MatrixBase<ParType> &solPars,
			Matrix solLB, Matrix solUB,
			Matrix solctrl, bool debugToggle, int verbose)
{
    mode_val = 0;
    mode = &mode_val;
    Matrix inform;
    Matrix hessi;
    Matrix p_hess;
    Matrix p_grad;
    struct Param_Obj pfunv;
    int i;
    if(verbose >= 3){
        mxLog("solPars is: \n");
        for (i = 0; i < solPars.size(); i++) mxLog("%f", solPars[i]);
        mxLog("4th call is: \n");
        mxLog("%2f", fit.solFun(solPars.derived().data(), mode, 0));
	Matrix eqv = fit.solEqBFun(0);
        mxLog("solEqBFun is: \n");
        for (i = 0; i < eqv.cols; i++) mxLog("%f",eqv.t[i]);
        mxLog("myineqFun is: \n");
        for (i = 0; i < fit.solIneqLB.size(); i++) mxLog("%f", fit.myineqFun(0).t[i]);
        mxLog("solLB is: \n");
        for (i = 0; i < solLB.cols; i++) mxLog("%f", solLB.t[i]);
        mxLog("solUB is: \n");
        for (i = 0; i < solUB.cols; i++) mxLog("%f", solUB.t[i]);
    }
    
    flag = 0; flag_L = 0; flag_U = 0;
    flag_NormgZ = 0; flag_step = 0; minr_rec = 0;
    
    double funv;
    double resultForTT;
    double solnp_nfn = 0;
    
    //time_t sec;
    //sec = time (NULL);
    
    double EMPTY = -999999.0;
    int maxit_trace = 0;
    
    Matrix grad = fill(solPars.size(), 1, (double)0.0);
    //free(matrices.front().t);
    LB.cols = solLB.cols;
    UB.cols = solUB.cols;
    Matrix pb_cont;
    Matrix difference1, difference2, tmpv, testMin, firstCopied, subnp_ctrl, subsetMat, temp2, temp1, temp, funv_mat, tempdf, firstPart, copied, subsetOne, subsetTwo, subsetThree, diff1, diff2, copyValues, diff, llist, tempTTVals, searchD;
    
    Matrix pars(solPars);
    
    if(verbose >= 2){
        mxLog("control is: \n");
        for (i = 0; i < solctrl.cols; i++) mxLog("%f",solctrl.t[i]);
    }
    
    
    if (LB.cols > 1){
        LB = duplicateIt(solLB);
        if (UB.cols < 1){
            UB = fill(LB.cols, 1, (double) DBL_MAX/2);
        }
    }
    else
    {
        LB = fill(1, 1, (double) 0.0);
        M(LB, 0, 0) = EMPTY;
    }
    
    if (UB.cols > 1){
        UB = duplicateIt(solUB);
        if (LB.cols < 1){
            LB = fill(UB.cols, 1, (double) -DBL_MAX/2);
        }
        
    }
    else{
        UB = fill(1, 1, (double) 0.0);
        M(UB, 0, 0) = EMPTY;
    }
    
    
    if(LB.cols > 1)
    {
        for (int i = 0; i < LB.cols; i++)
        {
            if (M(LB, i, 0) < M(pars, i, 0) && M(UB, i, 0) > M(pars, i, 0))
            { continue;  }
            else if (M(pars, i, 0) <= M(LB, i, 0))
            {   inform = fill(1, 1, 9);
                flag_L = 1;
                index_flag_L = i;
                M(pars, i, 0) = M(LB, i, 0) + M(solctrl, 4, 0);
            }
            else if (M(pars, i, 0) >= M(UB, i, 0))
            {   inform = fill(1, 1, 9);
                flag_U = 1;
                index_flag_U = i;
                M(pars, i, 0) = M(UB, i, 0) - M(solctrl, 4, 0);
            }
        }
    }
    
    if (verbose >= 2){
        mxLog("LB is: \n");
        for (i = 0; i < LB.cols; i++) mxLog("%f",LB.t[i]);
        
        mxLog("UB is: \n");
        for (i = 0; i < UB.cols; i++) mxLog("%f",UB.t[i]);
        
        mxLog("pars is: \n");
        for (i = 0; i < pars.cols; i++) mxLog("%f",pars.t[i]);
    }
    
    int np = pars.cols;
    inform = new_matrix(1, 1);
    hessi = new_matrix(np*np, 1);
    p_hess = new_matrix(np+(np*np), 1);
    p_grad = new_matrix(np+(np*np)+np, 1);
    
    ind.setZero();
    ind[indNumParam] = pars.cols;
    
    if (M(LB, 0, 0) != EMPTY || M(UB, 0, 0) != EMPTY){
        ind[indHasBounds] = 1;
    }
    
    // does not have a function gradient (currently not supported in Rsolnp)
    ind[indHasGradient] = 0;
    //# do function checks and return starting value
    
    funv = fit.solFun(pars.t, mode, verbose);
    
    // does not have a hessian (currently not supported in Rsolnp)
    ind[indHasHessian] = 0;
    // no jacobian (not implemented)
    ind[indHasJacobianIneq] = 0;
    
    // do inequality checks and return starting values
    const int nineq = fit.solIneqLB.size();
    ind[indIneqLength] = nineq;
    ind[indHasIneq] = nineq > 0;    

    Eigen::VectorXd Eineqx0(nineq);
    Eineqx0.setZero();
    Matrix ineqx0(Eineqx0);
    
    int neq;
    Matrix eqv = fit.solEqBFun(verbose);
    
    if(eqv.cols){
        ind[indHasEq] = 1;
        neq = eqv.cols;
        ind[indEqLength] = neq;
        ind[indHasJacobianEq] = 0;
    } else{
        ind[indHasEq] = 0;
        neq = 0;
        ind[indEqLength] = 0;
        ind[indHasJacobianEq] = 0;
    }
    if ( (ind[indHasBounds] > 0) || (ind[indHasIneq] > 0) ){
        ind[indHasBoundsOrIneq] = 1;
    }
    
    if (verbose >= 2){
        mxLog("ind is: \n");
        for (i = 0; i < ind.size(); i++) mxLog("%f",ind[i]);
    }
    
    
    Matrix pb;
    
    
    if(ind[indHasBoundsOrIneq])
    {   if((M(LB, 0, 0) != EMPTY) && nineq)
    {   pb = fill(2, nineq, (double)0.0);
        setColumnInplace(pb, fit.solIneqLB, 0);
        setColumnInplace(pb, fit.solIneqUB, 1);
        pb_cont = fill(2, np, (double)0.0);
        setColumnInplace(pb_cont, LB, 0);
        setColumnInplace(pb_cont, UB, 1);
        pb = transpose(copy(transpose(pb), transpose(pb_cont)));
    }
    else if((M(LB, 0, 0) == EMPTY) && nineq)
    {
        pb = fill(2, nineq, (double)0.0);
        setColumnInplace(pb, fit.solIneqLB, 0);
        setColumnInplace(pb, fit.solIneqUB, 1);
    }
    else if((M(LB, 0, 0) != EMPTY) && nineq==0)
    {
        pb = fill(2, np, (double)0.0);
        setColumnInplace(pb, LB, 0);
        setColumnInplace(pb, UB, 1);
    }
    }
    
    else    {pb = fill(1, 1, EMPTY);}
    
    double rho   = M(solctrl, 0, 0);
    int maxit = M(solctrl, 1, 0);
    int minit = M(solctrl, 2, 0);
    double delta = M(solctrl, 3, 0);
    double tol   = M(solctrl, 4, 0);
    double trace = M(solctrl, 5, 0);
    
    int tc = nineq + neq;
    
    double j = funv;
    Matrix jh = fill(1, 1, funv);
    Matrix tt = fill(1, 3, (double)0.0);
    
    Matrix lambda;
    Matrix constraint;
    
    if (tc > 0){
	    Matrix ineqv = fit.myineqFun(verbose);
        lambda = fill(1, tc, (double)0.0);
        
        if (nineq){
            if(eqv.cols)
            {
                constraint = copy(eqv, ineqv);
            }
            else{
                constraint = duplicateIt(ineqv);
            }
        }
        else    constraint = duplicateIt(eqv);
        
        if( ind[indHasIneq] > 0 ) {
            
            // 	tmpv = cbind(constraint[ (neq[0]):(tc[0]-1) ] - .fit.solIneqLB, .fit.solIneqUB - constraint[ (neq + 1):tc ] )
            Matrix difference1 = subset(constraint, 0, neq, tc-1);
            subtractEigen(difference1, fit.solIneqLB);
            Matrix difference2 = subset(constraint, 0, neq, tc-1);
            negate(difference2);
            addEigen(difference2, fit.solIneqUB);
            tmpv = fill(2, nineq, (double)0.0);
            setColumnInplace(tmpv, difference1, 0);
            setColumnInplace(tmpv, difference2, 1);
            testMin = rowWiseMin(tmpv);
            
            if (allGreaterThan(testMin, 0)) {
                ineqx0 = subset(constraint, 0, neq, tc-1);
            }
            
            Matrix diff = subset(constraint, 0, neq, tc-1);
            subtractEigen(diff, ineqx0);
            copyIntoInplace(constraint, diff, 0, neq, tc-1);
        }
        
        M(tt, 0, 1) = vnorm(constraint);
        double zeroCheck = M(tt, 0, 1) - (10 * tol);
        if( max(zeroCheck, nineq) <= 0 ) {
            rho = 0;
        }
    } // end if tc > 0
    else {
        lambda = fill(1, 1, (double)0.0);
    }
    
    
    Matrix tempv;
    Matrix p;
    
    if ( M(ineqx0, 0, 0) != EMPTY){
        p = copy(ineqx0, pars);
    }
    else{
        p = duplicateIt(pars);
    }
    
    Matrix hessv = diag(fill((np+nineq), 1, (double)1.0));
    
    double mu = np;
    
    int solnp_iter = 0;
    
    Matrix ob;
    Matrix funvMatrix = fill(1, 1, funv);
    
    if ( nineq){
	    Matrix ineqv = fit.myineqFun(verbose);
        if(eqv.cols){
            Matrix firstCopied = copy(funvMatrix, eqv);
            ob = copy(firstCopied, ineqv);
        }
        else{
            ob = copy(funvMatrix, ineqv);
        }
        
    }
    else if (eqv.cols) {
        ob = copy(funvMatrix, eqv);
    }
    else ob = duplicateIt(funvMatrix);
    
    
    if(verbose >= 3){
        mxLog("ob is: \n");
        for (i = 0; i < ob.cols; i++) mxLog("%f",ob.t[i]);
    }
    
    Matrix vscale;
    
    while(solnp_iter < maxit){
        solnp_iter = solnp_iter + 1;
        Eigen::Array<double, 5, 1> subnp_ctrl;
        subnp_ctrl[0] = rho;
        subnp_ctrl[1] = minit;
        subnp_ctrl[2] = delta;
        subnp_ctrl[3] = tol;
        subnp_ctrl[4] = trace;
        
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
        if ( ind[indHasBoundsOrIneq] <= 0){
            vscale = copy(vscale, p);
        }
        else{
            vscale = copy(vscale, fill(p.cols, 1, (double)1.0));
        }
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
        
        if (*mode == -1)
        {
            M(inform, 0, 0) = 0;
            hessi = MatrixToVector(fill(np, np, (double)0.0));
            p_hess = copy(p, hessi);
            p_grad = copy(p_hess, grad);
            pfunv.parameter = copy(p_grad, inform);
            pfunv.objValue = funv;
            return pfunv;
        }
        
        if (sx_Matrix.t == NULL) sx_Matrix = fill(p.cols, p.rows, (double)0.0);
        
        grad = subnp(p, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        
        if (flag == 1)
        {
            p = duplicateIt(resP);
            funv = fit.solFun(p.t, mode, verbose);
            funvMatrix = fill(1, 1, funv);
            eqv = fit.solEqBFun(verbose);
            if ( nineq) {
		    Matrix ineqv = fit.myineqFun(verbose);
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
            if ( ind[indHasBoundsOrIneq] <= 0){
                vscale = copy(vscale, p);
            }
            else{
                vscale = copy(vscale, fill(p.cols, 1, (double)1.0));
            }
            minMaxAbs(vscale, tol);
            lambda = duplicateIt(resY);
            hessv = duplicateIt(resHessv);
            mu = resLambda;
            grad = subnp(p, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        }
        p = duplicateIt(resP);
        
        lambda = duplicateIt(resY);
        
        hessv = duplicateIt(resHessv);
        
        mu = resLambda;
        
        
        Matrix temp = subset(p, 0, nineq, (nineq+np-1));
        funv = fit.solFun(temp.t, mode, verbose);
        if (*mode == -1)
        {
            M(inform, 0, 0) = 0;
            hessi = MatrixToVector(hessv);
            p_hess = copy(p, hessi);
            p_grad = copy(p_hess, grad);
            pfunv.parameter = copy(p_grad, inform);
            pfunv.objValue = funv;
            return pfunv;
        }
        
        solnp_nfn = solnp_nfn + 1;
        
        //Matrix funv_mat = fill(1, 1, funv);
        //Matrix tempdf = copy(temp, funv_mat);
        eqv = fit.solEqBFun(verbose);
        
        Matrix firstPart, copied;
        if (nineq){
		Matrix ineqv = fit.myineqFun(verbose);
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
    
    
    if ( ind[indHasIneq] > 0.5){
        ineqx0 = subset(p, 0, 0, nineq-1);
    }
    p = subset(p, 0, nineq, (nineq + np -1));
    
    if (false){
        /* TODO: LIST ERROR MESSAGES HERE */
    }
    else{
        double vnormValue;
        Matrix searchD;
        double iterateConverge;
        double iterateConvergeCond;
        Matrix tempTTVals = fill(2, 1, (double) 0.0);
        M(tempTTVals, 0, 0) = M(tt, 0, 0);
        M(tempTTVals, 1, 0) = M(tt, 0, 1);
        vnormValue = vnorm(tempTTVals);
        if (verbose >= 1)
        {
            mxLog("vnormValue \n");
            mxLog("%.20f", vnormValue);
        }
        searchD = duplicateIt(sx_Matrix);
        if (verbose >= 3){
            mxLog("searchD is: \n");
            for (i = 0; i < searchD.cols; i++) mxLog("%f",searchD.t[i]);
        }
        iterateConverge = delta * pow(vnorm(searchD),(double)2.0);
        if (verbose >= 1)
        {
            mxLog("vnorm(searchD) is: \n");
            mxLog("%.20f", vnorm(searchD));
            mxLog("iterateConverge is: \n");
            mxLog("%.20f", iterateConverge);
        }
        iterateConvergeCond = sqrt(tol) * ((double)1.0 + pow(vnorm(p), (double)2.0));
        if (verbose >= 1)
        {   mxLog("iterateConvergeCond is: \n");
            mxLog("%.20f", iterateConvergeCond);
        }
        
        if (vnormValue <= tol && flag_NormgZ == 1 && minr_rec == 1 && flag_step == 1){
            if (iterateConverge <= iterateConvergeCond){
                if (verbose >= 1){
                    mxLog("The solution converged in %d iterations. It is:", solnp_iter);}
                inform = fill(1, 1, 0);
            }
            else {
                if (verbose >= 1){
                    mxLog("The final iterate x satisfies the optimality conditions to the accuracy requested, but the sequence of iterates has not yet converged. CSOLNP was terminated because no further improvement could be made in the merit function.");}
                inform = fill(1, 1, 1);
                
            }
        }
        else{
            if (solnp_iter == maxit_trace)
            {
                if (verbose >= 1){
                    mxLog("Exiting after maximum number of iterations. Tolerance not achieved\n");}
                inform = fill(1, 1, 4);
            }
            else
            {
                if (verbose >= 1){
                    mxLog("Solution failed to converge.");
                }
                inform = fill(1, 1, 6);
            }
        }
    }
    
    hessi = fill(hessv.cols*hessv.cols,1 , (double)0.0);
    
    int ii;
    int ind_hess = 0;
    for (i = 0; i < hessv.cols; i++)
    {
        for (ii = 0; ii < hessv.rows; ii++)
        {
            M(hessi, ind_hess, 0) = M(hessv, i, ii);
            ind_hess = ind_hess + 1;
        }
    }
    
    if (verbose >= 1){
        mxLog("inform in subnp is: \n");
        for (i = 0; i < inform.cols; i++) mxLog("%f",inform.t[i]);
    }
    
    //hessi = MatrixToVector(hessv);
    p_hess = copy(p, hessi);
    p_grad = copy(p_hess, grad);
    pfunv.parameter = copy(p_grad, inform);
    pfunv.objValue = funv;
    
    //exit(0);
    return pfunv;
    
}

#endif // #define _OMX_CSOLNP_SPECIFIC_H
