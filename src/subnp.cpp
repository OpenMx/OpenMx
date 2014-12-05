#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include "matrix.h"
#include "omxCsolnp.h"

Param_Obj CSOLNP::solnp(Matrix solPars, double (*solFun)(Matrix, int*, int), Matrix solEqB, Matrix (*solEqBFun)(int), Matrix (*myineqFun)(int), Matrix solLB, Matrix solUB, Matrix solIneqUB, Matrix solIneqLB, Matrix solctrl, bool debugToggle, int verbose){
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
        for (i = 0; i < solPars.cols; i++) mxLog("%f", solPars.t[i]);
        mxLog("4th call is: \n");
        mxLog("%2f", solFun(solPars, mode, 0));
        mxLog("solEqB is: \n");
        for (i = 0; i < solEqB.cols; i++) mxLog("%f", solEqB.t[i]);
        mxLog("solEqBFun is: \n");
        for (i = 0; i < solEqB.cols; i++) mxLog("%f",solEqBFun(0).t[i]);
        mxLog("myineqFun is: \n");
        for (i = 0; i < solIneqLB.cols; i++) mxLog("%f", myineqFun(0).t[i]);
        mxLog("solLB is: \n");
        for (i = 0; i < solLB.cols; i++) mxLog("%f", solLB.t[i]);
        mxLog("solUB is: \n");
        for (i = 0; i < solUB.cols; i++) mxLog("%f", solUB.t[i]);
        mxLog("solIneqUB is: \n");
        for (i = 0; i < solIneqUB.cols; i++) mxLog("%f", solIneqUB.t[i]);
        mxLog("solIneqLB is: \n");
        for (i = 0; i < solIneqLB.cols; i++) mxLog("%f", solIneqLB.t[i]);
    }
    
    flag = 0; flag_L = 0; flag_U = 0;
    flag_NormgZ = 0; flag_step = 0; minr_rec = 0;
    
    Matrix pars;
    double funv;
    double resultForTT;
    double solnp_nfn = 0;
    
    //time_t sec;
    //sec = time (NULL);
    ind = fill(11, 1, (double) 0.0);
    
    double EMPTY = -999999.0;
    int maxit_trace = 0;
    
    Matrix grad = fill(solPars.cols, 1, (double)0.0);
    //free(matrices.front().t);
    ineqLB.cols = solIneqLB.cols;
    ineqUB.cols = solIneqUB.cols;
    LB.cols = solLB.cols;
    UB.cols = solUB.cols;
    Matrix ineqLBx;
    Matrix ineqUBx;
    Matrix pb_cont;
    Matrix difference1, difference2, tmpv, testMin, firstCopied, subnp_ctrl, subsetMat, temp2, temp1, temp, funv_mat, tempdf, firstPart, copied, subsetOne, subsetTwo, subsetThree, diff1, diff2, copyValues, diff, llist, tempTTVals, searchD;
    
    pars = duplicateIt(solPars);
    
    if(verbose >= 2){
        mxLog("control is: \n");
        for (i = 0; i < solctrl.cols; i++) mxLog("%f",solctrl.t[i]);
    }
    
    
    if (ineqLB.cols > 1){
        ineqLB = duplicateIt(solIneqLB);
        if (ineqUB.cols < 1){
            ineqUB = fill(ineqLB.cols, 1, (double) DBL_MAX/2);
            
        }
    }
    else
    {
        ineqLB = fill(1, 1, (double) 0.0);
        M(ineqLB, 0, 0) = EMPTY;
    }
    
    if (ineqUB.cols > 1){
        ineqUB = duplicateIt(solIneqUB);
        if (ineqLB.cols < 1){
            ineqLB = fill(ineqUB.cols, 1, (double) -DBL_MAX/2);
        }
    }
    else
    {
        ineqUB = fill(1, 1, (double) 0.0);
        M(ineqUB, 0, 0) = EMPTY;
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
    // [0] Rf_length of pars
    // [1] has function gradient?
    // [2] has hessian?
    // [3] has ineq?
    // [4] ineq Rf_length
    // [5] has jacobian (inequality)
    // [6] has eq?
    // [7] eq Rf_length
    // [8] has jacobian (equality)
    // [9] has upper / lower bounds
    // [10] has either lower/upper bounds or ineq
    
    M(ind, 0, 0) = pars.cols;
    
    if (M(LB, 0, 0) != EMPTY || M(UB, 0, 0) != EMPTY){
        M(ind, 9, 0) = 1;
    }
    
    // does not have a function gradient (currently not supported in Rsolnp)
    M(ind, 1, 0) = 0;
    //# do function checks and return starting value
    
    funv = solFun(pars, mode, verbose);
    
    // does not have a hessian (currently not supported in Rsolnp)
    M(ind, 2, 0) = 0;
    
    // do inequality checks and return starting values
    int nineq;
    Matrix ineqx0 = fill(ineqLB.cols, 1, (double)0.0);
    
    Matrix ineqv = myineqFun(verbose);
    
    if ( M(ineqv, 0, 0) != EMPTY){
        
        M(ind, 3, 0) = 1;
        nineq = ineqLB.cols;
        
        M(ind, 4, 0) = nineq;
        
        // check for infitnites/nans
        
        ineqLBx = ineqLB;
        ineqUBx = ineqUB;
        
        int i;
        for (i = 0; i<ineqLBx.cols; i++)
        {
            if (M(ineqLBx,i,0) <= -99999999.0){
                M(ineqLBx,i,0) = -1.0 * (1e10);
            }
            if (M(ineqUBx,i,0) >= DBL_MAX){
                M(ineqUBx,i,0) = 1e10;
            }
        }
        
        
        for (i = 0; i < ineqLBx.cols; i++)
        {
            M(ineqx0, i, 0) = (M(ineqLBx, i, 0) + M(ineqUBx, i, 0)) / 2.0;
        }
        
        // no jacobian
        M(ind, 5, 0) = 0;
        
    }
    else{
        
        M(ineqv, 0, 0) = EMPTY;
        M(ind, 3, 0) = 0;
        nineq = 0;
        M(ind, 4, 0) = 0;
        M(ind, 5, 0) = 0;
        M(ineqx0, 0, 0) = EMPTY;
    }
    
    int neq;
    Matrix eqv = solEqBFun(verbose);
    
    if( M(eqv, 0, 0) != EMPTY){
        M(ind, 6, 0) = 1;
        neq = solEqB.cols;
        M(ind, 7, 0) = neq;
        M(ind, 8, 0) = 0;
    } else{
        M(eqv, 0, 0) = EMPTY;
        M(ind, 6, 0) = 0;
        neq = 0;
        M(ind, 7, 0) = 0;
        M(ind, 8, 0) = 0;
    }
    if ( (M(ind, 9, 0) > 0) || (M(ind, 3, 0) > 0) ){
        M(ind, 10, 0) = 1;
    }
    
    if (verbose >= 2){
        mxLog("ind is: \n");
        for (i = 0; i < ind.cols; i++) mxLog("%f",ind.t[i]);
    }
    
    
    Matrix pb;
    
    
    if(M(ind, 10, 0))
    {   if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
    {   pb = fill(2, nineq, (double)0.0);
        setColumnInplace(pb, ineqLB, 0);
        setColumnInplace(pb, ineqUB, 1);
        pb_cont = fill(2, np, (double)0.0);
        setColumnInplace(pb_cont, LB, 0);
        setColumnInplace(pb_cont, UB, 1);
        pb = transpose(copy(transpose(pb), transpose(pb_cont)));
    }
    else if((M(LB, 0, 0) == EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
    {
        pb = fill(2, nineq, (double)0.0);
        setColumnInplace(pb, ineqLB, 0);
        setColumnInplace(pb, ineqUB, 1);
    }
    else if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) == EMPTY))
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
        lambda = fill(1, tc, (double)0.0);
        
        if (M(ineqv, 0, 0) != EMPTY){
            if(M(eqv,0,0) != EMPTY)
            {
                constraint = copy(eqv, ineqv);
            }
            else{
                constraint = duplicateIt(ineqv);
            }
        }
        else    constraint = duplicateIt(eqv);
        
        if( M(ind, 3, 0) > 0 ) {
            
            // 	tmpv = cbind(constraint[ (neq[0]):(tc[0]-1) ] - .ineqLB, .ineqUB - constraint[ (neq + 1):tc ] )
            Matrix difference1 = subset(constraint, 0, neq, tc-1);
            subtractEigen(difference1, ineqLB);
            Matrix difference2 = subset(constraint, 0, neq, tc-1);
            negate(difference2);
            addEigen(difference2, ineqUB);
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
    
    if ( M(ineqv, 0, 0) != EMPTY){
        if(M(eqv,0,0) != EMPTY){
            Matrix firstCopied = copy(funvMatrix, eqv);
            ob = copy(firstCopied, ineqv);
        }
        else{
            ob = copy(funvMatrix, ineqv);
        }
        
    }
    else if (M(eqv,0,0) != EMPTY){
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
        
        if ( M(ind, 6, 0) > 0){
            Matrix subsetMat = subset(ob, 0, 1, neq);
            double max = matrixMaxAbs(subsetMat);
            
            Matrix temp2 = fill(neq, 1, max);
            Matrix temp1 = fill(1, 1, M(ob, 0, 0));
            vscale = copy(temp1, temp2);
        }
        else{
            vscale = fill(1, 1, (double)1.0);
        }
        if ( M(ind, 10, 0) <= 0){
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
        
        grad = subnp(p, solFun, solEqBFun, myineqFun, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        
        if (flag == 1)
        {
            p = duplicateIt(resP);
            funv = solFun(p, mode, verbose);
            funvMatrix = fill(1, 1, funv);
            eqv = solEqBFun(verbose);
            ineqv = myineqFun(verbose);
            if ( M(ineqv, 0, 0) != EMPTY)
            {
                if(M(eqv,0,0) != EMPTY)
                {
                    Matrix firstCopied = copy(funvMatrix, eqv);
                    ob = copy(firstCopied, ineqv);
                }
                else{
                    ob = copy(funvMatrix, ineqv);
                }
            }
            else if (M(eqv,0,0) != EMPTY){
                ob = copy(funvMatrix, eqv);
            }
            else ob = duplicateIt(funvMatrix);
            
            if ( M(ind, 6, 0) > 0){
                Matrix subsetMat = subset(ob, 0, 1, neq);
                double max = matrixMaxAbs(subsetMat);
                
                Matrix temp2 = fill(neq, 1, max);
                Matrix temp1 = fill(1, 1, M(ob, 0, 0));
                vscale = copy(temp1, temp2);
                
            }
            else{
                vscale = fill(1, 1, (double)1.0);
            }
            if ( M(ind, 10, 0) <= 0){
                vscale = copy(vscale, p);
            }
            else{
                vscale = copy(vscale, fill(p.cols, 1, (double)1.0));
            }
            minMaxAbs(vscale, tol);
            lambda = duplicateIt(resY);
            hessv = duplicateIt(resHessv);
            mu = resLambda;
            grad = subnp(p, solFun, solEqBFun, myineqFun, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        }
        p = duplicateIt(resP);
        
        lambda = duplicateIt(resY);
        
        hessv = duplicateIt(resHessv);
        
        mu = resLambda;
        
        
        Matrix temp = subset(p, 0, nineq, (nineq+np-1));
        funv = solFun(temp, mode, verbose);
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
        eqv = solEqBFun(verbose);
        
        ineqv = myineqFun(verbose);
        
        Matrix firstPart, copied;
        if (M(ineqv, 0, 0) != EMPTY){
            if(M(eqv,0,0) != EMPTY){
                copied = copy(fill(1, 1, funv), eqv);
                ob = copy(copied, ineqv);
                
            }
            else{
                ob = copy(fill(1, 1, funv), ineqv);
            }
        }
        else if (M(eqv,0,0) != EMPTY){
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
            
            if ( M(ind, 3, 0) > 0.5){
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
    
    
    if ( M(ind, 3, 0) > 0.5){
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

template <typename T1>
Matrix CSOLNP::subnp(Matrix pars, double (*solFun)(Matrix, int*, int), Matrix (*solEqBFun)(int) ,  Matrix(*myineqFun)(int), Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,  Eigen::ArrayBase<T1> &ctrl, int verbose)
{
    double EMPTY = -999999.0;
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
    
    int neq =  (int)M(ind, 7, 0);
    int nineq = (int)M(ind, 4, 0);
    int np = (int)M(ind, 0, 0);
    
    double ch = 1;
    
    if (verbose >= 2){
        mxLog("ind inside subnp is: \n");
        for (int i = 0; i < ind.cols; i++) mxLog("%f",ind.t[i]);
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
    Matrix col1_pb;
    Matrix pb_cont;
    
    if(M(ind, 10, 0))
    {
        if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
        {   pb = fill(2, nineq, (double)0.0);
            setColumnInplace(pb, ineqLB, 0);
            setColumnInplace(pb, ineqUB, 1);
            
            pb_cont = fill(2, np, (double)0.0);
            setColumnInplace(pb_cont, LB, 0);
            setColumnInplace(pb_cont, UB, 1);
            
            pb = transpose(copy(transpose(pb), transpose(pb_cont)));//MAHSA
            
        }
        else if((M(LB, 0, 0) == EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
        {
            pb = fill(2, nineq, (double)0.0);
            setColumnInplace(pb, ineqLB, 0);
            setColumnInplace(pb, ineqUB, 1);
            
        }
        else if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) == EMPTY))
        {
            pb = fill(2, np, (double)0.0);
            setColumnInplace(pb, LB, 0);
            setColumnInplace(pb, UB, 1);
            
        }
    }
    else{
        pb = fill(1, 1, EMPTY);
    }
    
    if (verbose >= 3){
        mxLog("pb is: \n");
        for (int i = 0; i < pb.cols; i++) mxLog("%f",pb.t[i]);
    }
    
    Eigen::Array<double, 3, 1> sob;
    sob.setZero();
    
    //Matrix yyMatrix = duplicateIt(yy);
    
    divideEigen(ob, subset(vscale, 0, 0, nc));
    
    divideEigen(p0, subset(vscale, 0, (neq+1), (nc + np)));
    
    if (verbose >= 3){
        mxLog("p0_1 is: \n");
        for (int i = 0; i < p0.cols; i++) mxLog("%f",p0.t[i]);
    }
    
    int mm = -1;
    if (M(ind, 10, 0) > 0){
        if (M(ind, 9, 0) <= 0){
            mm = nineq;
        }
        else{
            mm=npic;
        }
        Matrix vscaleSubset = subset(vscale, 0, neq+1, neq+mm);
        //double vscaleSubsetLength = (neq+mm) - (neq+1) + 1;
        Matrix vscaleTwice = fill(pb.cols, pb.rows, (double)0.0);
        setColumnInplace(vscaleTwice, vscaleSubset, 0);
        setColumnInplace(vscaleTwice, vscaleSubset, 1);
        
        if (M(pb, 0, 0) != EMPTY){
            divideEigen(pb, vscaleTwice);
        }
    } // end if (ind [0][10] > 0)
    
    if (verbose >= 3){
        mxLog("pb is: \n");
        for (int i = 0; i < pb.cols; i++) mxLog("%f",pb.t[i]);
    }
    
    // scale the lagrange multipliers and the Hessian
    if( nc > 0) {
        // yy [total constraints = nineq + neq]
        // scale here is [tc] and dot multiplied by yy
        //yy = vscale[ 2:(nc + 1) ] * yy / vscale[ 1 ]
        Matrix result = transpose(subset(vscale, 0, 1, nc));
        multiplyEigen(result, yy);
        yy = duplicateIt(result);
        divideByScalar2D(yy, M(vscale, 0, 0));
    }
    
    // hessv [ (np+nineq) x (np+nineq) ]
    // hessv = hessv * (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1)]) ) / vscale[ 1 ]
    
    Matrix vscaleSubset = subset(vscale, 0, (neq+1), (nc + np));
    transposeDP(vscaleSubset);
    multiplyEigen(hessv, vscaleSubset);
    divideByScalar2D(hessv, M(vscale, 0, 0));
    
    j = M(ob, 0, 0);
    if (verbose >= 3){
        mxLog("j j is: \n");
        mxLog("%2f", j);
    }
    Matrix a;
    
    if( M(ind, 3, 0) > 0){
        if ( M(ind, 6, 0) <= 0)
        {
            // arrays, rows, cols
            Matrix onesMatrix = fill(nineq, 1, (double)-1.0);
            Matrix negDiag = diag(onesMatrix);
            Matrix zeroMatrix = fill(np, nineq, (double)0.0);
            // a = cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) )
            a = copy(negDiag, zeroMatrix);
        }
        else{
            // [ (neq+nineq) x (nineq+np)]
            //a = rbind( cbind( 0 * .ones(neq, nineq), matrix(0, ncol = np, nrow = neq) ),
            //      cbind( -diag(nineq), matrix(0, ncol = np, nrow = nineq) ) )
            
            Matrix zeroMatrix = fill(np, nineq, (double)0.0);
            Matrix firstHalf = copy(fill(nineq, neq, (double)0.0), fill(np, neq, (double)0.0));
            Matrix onesMatrix = fill(nineq, 1, (double)-1.0);
            Matrix negDiag = diag(onesMatrix);
            Matrix secondHalf = copy(negDiag, zeroMatrix);
            a = transpose(copy(transpose(firstHalf), transpose(secondHalf)));
        }
    }	// end 	if(ind[0][3] > 0){
    
    if ( (M(ind, 6, 0) > 0) && M(ind, 3, 0) <= 0 ){
        a = fill(np, neq, (double)0.0);
    }
    if (M(ind, 6, 0)<= 0 && (M(ind, 3, 0) <= 0)){
        a = fill(np, 1, (double)0.0);
    }
    Matrix g = fill(npic, 1, (double)0.0);
    Matrix p = subset(p0, 0, 0, (npic-1));
    
    Matrix dx;
    Matrix b;
    double funv;
    Matrix eqv;
    Matrix ineqv;
    Matrix tmpv;
    Matrix constraint;
    Matrix gap;
    
    int solnp_nfn = 0;
    double go, reduce = 1e-300;
    int minit;
    double lambdaValue = lambda;
    
    if (nc > 0) {
        constraint = subset(ob, 0, 1, nc);
        
        for (int i=0; i<np; i++){
            int index = nineq + i;
            M(p0, index, 0) = M(p0, index, 0) + delta;
            Matrix tmpv = subset(p0, 0, nineq, (npic-1));
            multiplyEigen(tmpv, subset(vscale, 0, (nc+1), (nc+np)));
            
            if (verbose >= 2){
                mxLog("7th call is \n");
            }
            funv = solFun(tmpv, mode, verbose);
            
            eqv = solEqBFun(verbose);
            
            ineqv = myineqFun(verbose);
            
            solnp_nfn = solnp_nfn + 1;
            Matrix firstPart;
            Matrix firstPartt;
            Matrix secondPart;
            
            if (M(ineqv,0,0) != EMPTY){
                if(M(eqv,0,0) != EMPTY)
                {
                    firstPartt = copy(fill(1, 1, funv), eqv);
                    firstPart = copy(firstPartt, ineqv);
                }
                else{
                    firstPart = copy(fill(1, 1, funv), ineqv);
                }
            }
            else if (M(eqv,0,0) != EMPTY){
                firstPart = copy(fill(1, 1, funv), eqv);
            }
            else firstPart = fill(1, 1, funv);
            secondPart = subset(vscale, 0, 0, nc);
            divideEigen(firstPart, secondPart);
            ob = duplicateIt(firstPart);
            
            M(g, index, 0) = (M(ob, 0, 0)-j) / delta;
            
            if (verbose >= 3){
                mxLog("g is: \n");
                for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
                mxLog("a is: \n");
                for (int ilog = 0; ilog < a.cols; ilog++) mxLog("%f",a.t[ilog]);
                
            }
            Matrix colValues = subset(ob, 0, 1, nc);
            subtractEigen(colValues, constraint);
            divideByScalar2D(colValues, delta);
            setColumnInplace(a, colValues, index);
            M(p0, index, 0) = M(p0, index, 0) - delta;
        } // end for (int i=0; i<np, i++){
        
        if (*mode == -1)
        {
            funv = 1e24;
            *mode = 0;
        }
        
        if(M(ind, 3, 0) > 0){
            //constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
            Matrix firstPart, secondPart;
            firstPart  = subset(constraint, 0, neq, (neq+nineq-1));
            secondPart = subset(p0, 0, 0, (nineq-1));
            subtractEigen(firstPart, secondPart);
            copyIntoInplace(constraint, firstPart, 0, neq, (neq+nineq-1));
            
        }
        
        if (false && solvecond(a) > 1/DBL_EPSILON) { // this can't be the cheapest way to check TODO
            Rf_error("Redundant constraints were found. Poor intermediate results may result. "
                     "Remove redundant constraints and re-OPTIMIZE.");
        }
        
        //b = fill(nc, 1, (double)0.0);
        
        b = transpose(timess(a, transpose(p0)));
        //  b [nc,1]
        subtractEigen(b, constraint);
        
        ch = -1;
        alp[0] = tol - matrixMaxAbs(constraint);
        if (alp[0] <= 0){
            
            ch = 1;
            
            if ( M(ind, 10, 0) < 0.5){
                //Matrix dotProd = transposeDotProduct(a); //Mahsa: this is equal to "a %*% t(a)"
                //Matrix solution = solve(dotProd, constraint);
                
                //p0 = subtract(p0, matrixDotProduct(transpose(a), solution));
                
                alp[0] = 1;
            }
            
        } // end if (alp[0][0] <= 0){
        
        if (alp[0] <= 0){
            int npic_int = npic;
            p0 = copy(p0, fill(1, 1, (double)1.0));
            
            multiplyByScalar2D(constraint, -1.0);
            a = copy(a, transpose(constraint));
            Matrix cx = copy(fill(npic, 1, (double)0.0), fill(1, 1, (double)1.0));
            
            dx = fill(1, npic+1, (double)1.0);
            
            go = 1;
            minit = 0;
            
            while(go >= tol)
            {
                minit = minit + 1;
                gap = fill(2, mm, (double)0.0);
                Matrix result = subset(p0, 0, 0, mm-1);
                subtractEigen(result, getColumn(pb, 0));
                setColumnInplace(gap, result, 0);
                Matrix result1 = getColumn(pb, 1);
                subtractEigen(result1, subset(p0, 0, 0, mm-1));
                setColumnInplace(gap, result1, 1);
                rowSort(gap);
                Matrix dx_t = transpose(dx);
                copyInto(dx_t, getColumn(gap,0), 0, 0, mm-1);
                dx = duplicateIt(dx_t);
                
                M(dx, npic_int, 0) = M(p0, npic_int, 0);
                
                if(M(ind, 9, 0) <= 0)
                {
                    Matrix argum = fill(1, npic-mm, (double)1.0);
                    multiplyByScalar2D(argum, max(findMax(subset(dx, 0, 0, mm-1)), 100));
                    
                    copyIntoInplace(dx, argum, 0, mm, npic-1);
                    
                }
                
                dx = transpose(dx);
                
                Matrix argum1 = transpose(timess(a, transpose(diag(dx))));
                
                Matrix argum2 = duplicateIt(dx);
                multiplyEigen(argum2, transpose(cx));
                
                
                Matrix y = QRdsolve(argum1, argum2);
                
                Matrix t_cx = transpose(cx);
                subtractEigen(t_cx, timess(transpose(a),y));
                Matrix dx_copy = duplicateIt(dx);//MAHSA
                multiplyEigen(dx_copy, t_cx);
                multiplyEigen(dx, dx_copy);
                Matrix v = transpose(dx);//MAHSA
                
                int indexx = npic;
                
                if (M(v, indexx, 0) > 0)
                {
                    double z = M(p0, indexx, 0)/M(v, indexx, 0);
                    
                    for (int i=0; i<mm; i++)
                    {
                        if(M(v, i, 0) < 0)
                        {
                            z = min(z, -(M(pb, 1, i) - M(p0, i, 0))/M(v, i, 0));
                            
                        }
                        else if(M(v, i, 0) > 0)
                        {
                            
                            z = min(z, (M(p0, i, 0) - M(pb, 0, i))/M(v, i, 0));
                        }
                    }
                    
                    if(z < (M(p0, indexx, 0)/M(v, indexx, 0))) {
                        z *= 0.9;
                    }
                    
                    Eigen::Map< Eigen::VectorXd > Ep0(p0.t, p0.cols);
                    Eigen::Map< Eigen::VectorXd > Ev(v.t, v.cols);
                    Ep0 -= Ev * z;
                    
                    go = M(p0, indexx, 0);
                    
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
            
            for (h = 0; h<a.rows; h++)
            {
                setRowInplace(aMatrix, h, subset(getRow(a, h), 0, 0, npic-1));
            }
            a = duplicateIt(aMatrix);
            
            b = timess(a, transpose(subset(p0, 0, 0, npic-1)));
            
        }// end if(M(alp, 0, 0) <= 0)
    } // end if (nc > 0){
    
    p = subset(p0, 0, 0, npic-1);
    if (nc == 0)    t_sol = fill(1, 1, (double)0.0);
    
    if (verbose >= 3){
        mxLog("p is: \n");
        for (int i = 0; i < p.cols; i++) mxLog("%f",p.t[i]);
    }
    
    if (ch > 0){
        tmpv = subset(p, 0, nineq, (npic-1));
        multiplyEigen(tmpv, subset(vscale, 0, (nc+1), (nc+np)));
        if (verbose >= 2){
            mxLog("tmpv is: \n");
            for (int i = 0; i < tmpv.cols; i++) mxLog("%f",tmpv.t[i]);
            mxLog("8th call is \n");
        }
        funv = solFun(tmpv, mode, verbose);
        if (verbose >= 3){
            mxLog("funv is: \n");
            mxLog("%2f", funv);
        }
        
        if (*mode == -1)
        {
            funv = 1e24;
            *mode = 0;
        }
        
        eqv = solEqBFun(verbose);
        if (verbose >= 3){
            mxLog("eqv is: \n");
            for (int i = 0; i < eqv.cols; i++) mxLog("%f",eqv.t[i]);
        }
        ineqv = myineqFun(verbose);
        if (verbose >= 3){
            mxLog("ineqv is: \n");
            for (int i = 0; i < ineqv.cols; i++) mxLog("%f",ineqv.t[i]);
        }
        solnp_nfn = solnp_nfn + 1;
        Matrix firstPart, secondPart, firstPartt;
        if ( M(ineqv,0,0) != EMPTY){
            if (M(eqv,0,0) != EMPTY){
                firstPartt = copy(fill(1, 1, funv), eqv);
                firstPart = copy(firstPartt, ineqv);
            }
            else{
                firstPart = copy(fill(1, 1, funv), ineqv);
            }
        }
        else if (M(eqv,0,0) != EMPTY){
            firstPart = copy(fill(1, 1, funv), eqv);
        }
        else firstPart = fill(1, 1, funv);
        secondPart = subset(vscale, 0, 0, nc);
        divideEigen(firstPart, secondPart);
        ob = duplicateIt(firstPart);
        
    } // end of if (ch>0)
    
    if (verbose >= 3){
        mxLog("ob is: \n");
        for (int i = 0; i < ob.cols; i++) mxLog("%f",ob.t[i]);
    }
    
    j = M(ob, 0, 0);
    
    if (M(ind, 3, 0) > 0){
        Matrix result = subset(ob, 0, neq+1, nc);
        subtractEigen(result, subset(p, 0, 0, nineq-1));
        copyIntoInplace(ob, result, 0, neq+1, nc);
    }
    
    if (nc > 0){
        Matrix result = subset(ob, 0, 1, nc);
        subtractEigen(result, transpose(timess(a, transpose(p))));
        addEigen(result, b);
        copyIntoInplace(ob, result, 0, 1, nc);
        
        Matrix temp = subset(ob, 0, 1, nc);
        
        double vnormTerm = vnorm(temp) * vnorm(temp);
        Matrix yyTerm = transpose(yy);
        
        double dotProductTerm = dotProduct(getRow(yyTerm, 0), getRow(temp, 0));
        
        j = M(ob, 0, 0) - dotProductTerm + rho * vnormTerm;
        
    } // end if (nc > 0)
    
    minit = 0;
    Matrix obm;
    Matrix yg = fill(npic, 1, (double)0.0);
    Matrix yg_rec = fill(2, 1, (double)0.0);
    Matrix sx;
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
    Matrix ob1, ob2, ptt, t18, p0_1_t, ptt2, t19, pttCol, t20, ob3;
    Matrix partOne, partTwo, tempPttCol;
    Matrix yyTerm, firstp, t21, t22, t23, t24, t25, t26, t27;
    Matrix p_copy, p0_copy, pttColOne, t29, t30, t31, t32, t33, t34, temp, tempCol;
    Matrix input, rhs;
    
    //mxLog("Maxit is: %d", maxit);
    //mxLog("minit is: %d", minit);
    
    
    while (minit < maxit){
        minit = minit + 1;
        //mxLog("minit is: %d", minit);
        if (ch > 0){
            
            for (int i=0; i<np; i++){
                int index = nineq + i;
                M(p, index, 0) = M(p, index, 0) + delta;
                if (tmpv.t == NULL) {tmpv = new_matrix(npic - nineq, 1);
                }
                subsetEigen(tmpv, p, 0, nineq, (npic-1));
                if (vscale_t.t == NULL) {vscale_t = new_matrix(nc+np-nc, 1);
                }
                subsetEigen(vscale_t, vscale, 0, (nc+1), (nc+np));
                multiplyEigen(tmpv, vscale_t);
                if (verbose >= 3){
                    mxLog("9th call is \n");
                    
                }
                funv = solFun(tmpv, mode, verbose);
                if (verbose >= 3){
                    mxLog("funv is: \n");
                    mxLog("%2f", funv);
                }
                
                if (*mode == -1)
                {
                    funv = 1e24;
                    *mode = 0;
                }
                eqv = solEqBFun(verbose);
                
                ineqv = myineqFun(verbose);
                
                solnp_nfn = solnp_nfn + 1;
                
                if (M(ineqv, 0, 0) != EMPTY){
                    if(M(eqv,0,0) != EMPTY)
                    {
                        if(t1.t == NULL) t1 = new_matrix(1, 1);
                        fill_t(t1, 1, 1, funv);
                        if (firstPartt.t == NULL) firstPartt = new_matrix(t1.cols + eqv.cols, t1.rows);
                        copyEigen(firstPartt, t1, eqv);
                        
                        if (firstPart.t == NULL) firstPart = new_matrix(firstPartt.cols + ineqv.cols, firstPartt.rows);
                        copyEigen(firstPart, firstPartt, ineqv);
                        
                    }
                    else{
                        if(t1.t == NULL) t1 = new_matrix(1, 1);
                        fill_t(t1, 1, 1, funv);
                        if (firstPart.t == NULL) firstPart = new_matrix(t1.cols + ineqv.cols, t1.rows);
                        copyEigen(firstPart, t1, ineqv);
                        
                    }
                }
                else if (M(eqv,0,0) != EMPTY){
                    if(t1.t == NULL) t1 = new_matrix(1, 1);
                    fill_t(t1, 1, 1, funv);
                    if (firstPart.t == NULL) firstPart = new_matrix(t1.cols + eqv.cols, t1.rows);
                    copyEigen(firstPart, t1, eqv);
                }
                else {
                    if(firstPart.t == NULL) firstPart = new_matrix(1, 1);
                    fill_t(firstPart, 1, 1, funv);
                    
                }
                
                if (secondPart.t == NULL) secondPart = new_matrix(nc+1, 1);
                subsetEigen(secondPart, vscale, 0, 0, nc);
                divideEigen(firstPart, secondPart);
                if (obm_t.t == NULL) obm_t = new_matrix(firstPart.cols, firstPart.rows);
                
                if (obm.t == NULL) obm = new_matrix(obm_t.cols, obm_t.rows);
                duplicateIt_t(obm_t, firstPart);
                M(obm, 0, 0) = M(obm_t, 0, 0);
                
                if (verbose >= 3){
                    mxLog("obm is: \n");
                    for (int ilog = 0; ilog < obm.cols; ilog++) mxLog("%f",obm.t[ilog]);
                    mxLog("obm_t is: \n");
                    for (int ilog = 0; ilog < obm_t.cols; ilog++) mxLog("%f",obm_t.t[ilog]);
                    mxLog("j is: \n");
                    mxLog("%.20f", j);
                }
                
                if (M(ind, 3, 0) > 0.5){
                    if (result.t == NULL) result = new_matrix(nc - neq, 1);
                    subsetEigen(result, obm_t, 0, neq+1, nc);
                    if (result1.t == NULL) result1 = new_matrix(nineq, 1);
                    subsetEigen(result1, p, 0, 0, nineq-1);
                    subtractEigen(result, result1);
                    copyIntoInplace(obm_t, result, 0, neq+1, nc);
                }
                
                if (nc > 0){
                    if (first_part.t == NULL) first_part = new_matrix(nc, 1);
                    subsetEigen(first_part, obm_t, 0, 1, nc);
                    
                    if (p_t.t == NULL) p_t = new_matrix(p.rows, p.cols);
                    transpose_t(p_t, p);
                    
                    if (atime.t == NULL) atime = new_matrix(p_t.cols, a.rows);
                    timessEigen(atime, a, p_t);
                    
                    if (t_atime.t == NULL) t_atime = new_matrix(atime.rows, atime.cols);
                    transpose_t(t_atime, atime);
                    
                    subtractEigen(first_part, t_atime);
                    
                    addEigen(first_part, b);
                    
                    copyIntoInplace(obm_t, first_part, 0, 1, nc);
                    if (temp_t.t == NULL) temp_t = new_matrix(nc, 1);
                    subsetEigen(temp_t, obm_t, 0, 1, nc);
                    double vnormTerm = vnorm(temp_t) * vnorm(temp_t);
                    if (t_yy.t == NULL) t_yy = new_matrix(yy.rows, yy.cols);
                    transpose_t(t_yy, yy);
                    if (row_tt_y.t == NULL) row_tt_y = new_matrix(t_yy.cols, 1);
                    getRow_t(row_tt_y, t_yy, 0);
                    if (row_temp_t.t == NULL) row_temp_t = new_matrix(temp_t.cols, 1);
                    getRow_t(row_temp_t, temp_t, 0);
                    double dotProductTerm = dotProduct(row_tt_y, row_temp_t);
                    double newOBMValue = M(obm_t, 0, 0) - dotProductTerm + rho * vnormTerm;
                    if (obm.t == NULL) obm = new_matrix(1, 1);
                    M(obm, 0, 0) = newOBMValue;
                }
                //exit(0);
                
                M(g, index, 0) = (M(obm, 0, 0) - j)/delta;
                M(p, index, 0) = M(p, index, 0) - delta;
                
                if (verbose >= 3){
                    mxLog("g is: \n");
                    for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
                    mxLog("p is: \n");
                    for (int ilog = 0; ilog < p.cols; ilog++) mxLog("%f",p.t[ilog]);
                }
            } // end for (i=0; i<np; i++){
            
            if (M(ind, 3, 0) > 0.5){
                if (t2.t == NULL) t2 = new_matrix(nineq, 1);
                fill_t(t2, nineq, 1, (double)0.0);
                copyIntoInplace(g, t2, 0, 0, (nineq-1));
            }
        } // end if (ch > 0){
        
        if (verbose >= 3){
            mxLog("yg is: \n");
            for (int ilog = 0; ilog < yg.cols; ilog++) mxLog("%f",yg.t[ilog]);
        }
        
        if (minit > 1){
            negate(yg);
            addEigen(yg, g);
            negate(sx);
            addEigen(sx, p);
            negate(yg);
            negate(sx);
            
            if (t3.t == NULL) t3 = new_matrix(sx.rows, sx.cols);
            transpose_t(t3, sx);
            if (t4.t == NULL) t4 = new_matrix(hessv.cols, sx.rows);
            timessEigen(t4, sx, hessv);
            if (m_sc.t == NULL) m_sc = new_matrix(t3.cols, t4.rows);
            timessEigen(m_sc, t4, t3);
            M(sc, 0, 0) = M(m_sc, 0, 0);
            if (t5.t == NULL) t5 = new_matrix(yg.rows, yg.cols);
            transpose_t(t5, yg);
            if (m_sc2.t == NULL) m_sc2 = new_matrix(t5.cols, sx.rows);
            timessEigen(m_sc2, sx, t5);
            M(sc, 1, 0) = M(m_sc2, 0, 0);
            
            if ((M(sc, 0, 0) * M(sc, 1, 0)) > 0){
                //hessv  = hessv - ( sx %*% t(sx) ) / sc[ 1 ] + ( yg %*% t(yg) ) / sc[ 2 ]
                if (sx2.t == NULL) sx2 = new_matrix(t3.cols, hessv.rows);
                timessEigen(sx2, hessv, t3);
                if (t6.t == NULL) t6 = new_matrix(sx2.rows, sx2.cols);
                transpose_t(t6, sx2);
                if(sxMatrix.t == NULL) sxMatrix = new_matrix(t6.cols, sx2.rows);
                timessEigen(sxMatrix, sx2, t6);
                divideByScalar2D(sxMatrix, M(sc, 0, 0));
                if (t_yg.t == NULL) t_yg = new_matrix(yg.rows, yg.cols);
                transpose_t(t_yg, yg);
                if (ygMatrix.t == NULL) ygMatrix = new_matrix(yg.cols, t_yg.rows);
                timessEigen(ygMatrix, t_yg, yg);
                divideByScalar2D(ygMatrix, M(sc, 1, 0));
                subtractEigen(hessv, sxMatrix);
                addEigen(hessv, ygMatrix);
            }
        }
        if (verbose >= 3){
            mxLog("yg is: \n");
            for (int ilog = 0; ilog < yg.cols * yg.rows; ilog++) mxLog("%f",yg.t[ilog]);
            mxLog("sx is: \n");
            for (int ilog = 0; ilog < sx.cols * sx.rows; ilog++) mxLog("%f",sx.t[ilog]);
            mxLog("sc is: \n");
            for (int ilog = 0; ilog < sc.cols * sc.rows; ilog++) mxLog("%f",sc.t[ilog]);
            mxLog("hessv is: \n");
            for (int ilog = 0; ilog < hessv.cols * hessv.rows; ilog++) mxLog("%f",hessv.t[ilog]);
        }
        
        dx = fill(npic, 1, 0.01);
        
        if (M(ind, 10, 0) > 0.5){
            if (gap1.t == NULL) gap1 = new_matrix(pb.cols, pb.rows);
            fill_t(gap1, pb.cols, pb.rows, (double)0.0);
            if (res.t == NULL) res = new_matrix(mm, 1);
            subsetEigen(res, p, 0, 0, mm-1);
            if (col_pb.t == NULL) col_pb = new_matrix(pb.rows, 1);
            getColumn_t(col_pb, pb, 0);
            subtractEigen(res, col_pb);
            setColumnInplace(gap1, res, 0);
            if (col_pb2.t == NULL) col_pb2 = new_matrix(pb.rows, 1);
            getColumn_t(col_pb2, pb, 1);
            if (t7.t == NULL) t7 = new_matrix(mm, 1);
            subsetEigen(t7, p, 0, 0, mm-1);
            subtractEigen(col_pb2, t7);
            setColumnInplace(gap1, col_pb2, 1);
            rowSort(gap1);
            if (t8.t == NULL) t8 = new_matrix(1, mm);
            fill_t(t8, 1, mm, (double)1.0);
            multiplyByScalar2D(t8, sqrt(DBL_EPSILON));
            if (gap_c.t == NULL) gap_c = new_matrix(gap1.rows, 1);
            getColumn_t(gap_c, gap1, 0);
            addEigen(gap_c, t8);
            
            if (gap2.t == NULL) gap2 = new_matrix(gap_c.cols, gap_c.rows);
            duplicateIt_t(gap2, gap_c);
            if (t9.t == NULL) t9 = new_matrix(mm, 1);
            fill_t(t9, mm, 1, (double)1.0);
            divideEigen(t9, gap2);
            copyIntoInplace(dx, t9, 0, 0, mm-1);
            
            if (verbose >= 3){
                mxLog("dx is: \n");
                for (int ilog = 0; ilog < dx.cols; ilog++) mxLog("%f",dx.t[ilog]);
            }
            
            if(M(ind, 9, 0) <= 0)
            {
                if (t11.t == NULL) t11 = new_matrix(1, npic-mm);
                fill_t(t11, 1, npic-mm, (double)1.0);
                if (t10.t == NULL) t10 = new_matrix(mm, 1);
                subsetEigen(t10, dx, 0, 0, mm-1);
                multiplyByScalar2D(t11, min(findMin(t10), 0.01));
                copyIntoInplace(dx, t11, 0, mm, npic-1);
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
            if (dx_c.t == NULL) dx_c = new_matrix(dx.cols, dx.rows);
            duplicateIt_t(dx_c, dx);
            
            multiplyEigen(dx_c, dx_c);
            if (verbose >= 3){
                mxLog("dx_c is: \n");
                for (int ilog = 0; ilog < dx_c.cols * dx_c.rows; ilog++) mxLog("%f",dx_c.t[ilog]);
            }
            
            if (dxDiag.t == NULL)
            {
                if (dx_c.cols > dx_c.rows)
                {   dxDiag = new_matrix(dx_c.cols, dx_c.cols);
                    fill_t(dxDiag, dx_c.cols, dx_c.cols, (double)0.0);
                }
                else
                {   dxDiag = new_matrix(dx_c.rows, dx_c.rows);
                    fill_t(dxDiag, dx_c.rows, dx_c.rows, (double)0.0);
                }
            }
            
            diag_t(dxDiag, dx_c);
            
            if (verbose >= 3){
                mxLog("dxDiag is: \n");
                for (int ilog = 0; ilog < dxDiag.cols * dxDiag.rows; ilog++) mxLog("%f",dxDiag.t[ilog]);
            }
            multiplyByScalar2D(dxDiag, lambdaValue);
            if (verbose >= 3){
                mxLog("dxDiag is: \n");
                for (int ilog = 0; ilog < dxDiag.cols * dxDiag.cols; ilog++) mxLog("%f",dxDiag.t[ilog]);
            }
            
            if (cz.t == NULL) cz = new_matrix(hessv.cols, hessv.rows);
            duplicateIt_t(cz, hessv);
            
            addEigen(cz, dxDiag);
            if (verbose >= 3){
                mxLog("cz is: \n");
                for (int ilog = 0; ilog < cz.cols * cz.cols; ilog++) mxLog("%f",cz.t[ilog]);
            }
            chol_lpk(cz);
            if (verbose >= 3){
                mxLog("cz is: \n");
                for (int ilog = 0; ilog < cz.cols * cz.cols; ilog++) mxLog("%f",cz.t[ilog]);
            }
            if (!R_FINITE(findMax(cz)))
            {
                mxLog("here in findMax");
                double nudge = 1.490116e-08;
                flag = 1;
                if (t12.t == NULL) t12 = new_matrix(nc + np - neq, 1);
                subsetEigen(t12, vscale, 0, (neq+1), (nc + np));
                multiplyEigen(p, t12);
                if (flag_L) {M(p, index_flag_L, 0) = M(p, index_flag_L, 0) + nudge;}
                if (flag_U) {M(p, index_flag_U, 0) = M(p, index_flag_U, 0)- nudge;}
                if (nc > 0)
                {
                    if (y.t == NULL) y = new_matrix(1, 1);
                    fill_t(y, 1, 1, (double)0.0);
                }
                if (t13.t == NULL) t13 = new_matrix(t12.cols, t12.cols);
                transposeDP_t(t13, t12);
                divideEigen(hessv, t13);
                multiplyByScalar2D(hessv, M(vscale, 0, 0));
                /*if (resP.t == NULL) resP = new_matrix(p.cols, p.rows);
                 duplicateIt_t(resP, p);
                 if (resY.t == NULL) resY = new_matrix(y.cols, y.rows);
                 duplicateIt_t(resY, y);
                 if (resHessv.t == NULL) resHessv = new_matrix(hessv.cols, hessv.rows);
                 duplicateIt_t(resHessv, hessv);*/
                resP = duplicateIt(p);
                resY = duplicateIt(y);
                resHessv = duplicateIt(hessv);
                resLambda = lambda;
                return g;
            }
            
            //Matrix identityMatrix = diag(fill(hessv.cols, 1, (double)1.0));
            
            solveinv(cz);
            //cz = luSolve(cz, identityMatrix);
            
            if (verbose >= 3){
                mxLog("cz.rows: %d", cz.rows);
                mxLog("cz.cols: %d", cz.cols);
                
                mxLog("cz is: \n");
                for (int ilog = 0; ilog < cz.cols*cz.cols; ilog++) mxLog("%f",cz.t[ilog]);
            }
            
            //Matrix getRowed = getRow(cz, 0);
            if (verbose >= 3){
                mxLog("g is: \n");
                for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
            }
            //Matrix getRowedtwo = getRow(g, 0);
            //double rr = dotProduct(getRowed, getRowedtwo);
            if (t14.t == NULL) t14 = new_matrix(g.rows, g.cols);
            transpose_t(t14, g);
            if (t15.t == NULL) t15 = new_matrix(cz.rows, cz.cols);
            transpose_t(t15, cz);
            if (verbose >= 3){
                mxLog("yg.rows: %d", yg.rows);
                mxLog("yg.cols: %d", yg.cols);
                mxLog("yg is: \n");
                for (int ilog = 0; ilog < yg.cols * yg.rows; ilog++) mxLog("%f",yg.t[ilog]);
            }
            if (yg2.t == NULL) yg2 = new_matrix(t14.cols, t15.rows);
            timessEigen(yg2, t15, t14);
            
            if (minit == 1) M(yg_rec, 0, 0) = vnorm(yg2);
            
            if (verbose >= 3){
                mxLog("yg is: \n");
                for (int ilog = 0; ilog < yg.cols * yg.rows; ilog++) mxLog("%f", yg.t[ilog]);
            }
            if (nc <= 0){
                divideByScalar2D(cz, -1.0);
                if (u.t == NULL) u = new_matrix(yg2.cols, cz.rows);
                timessEigen(u, cz, yg2);
                divideByScalar2D(cz, -1.0);
                if (t_u.t == NULL) t_u = new_matrix(u.rows, u.cols);
                transpose_t(t_u, u);
                if (verbose >= 3){
                    mxLog("u inside nc <=0 is: \n");
                    for (int ilog = 0; ilog < t_u.cols * t_u.rows; ilog++) mxLog("%f",t_u.t[ilog]);
                }
            }
            else{
                //y = qr.solve(t(cz) %*% t(a), yg)
                if (aTranspose.t == NULL) aTranspose = new_matrix(a.rows, a.cols);
                transpose_t(aTranspose, a);
                if (verbose >= 3){
                    mxLog("aTranspose is: \n");
                    for (int ilog = 0; ilog < aTranspose.cols; ilog++) mxLog("%f",aTranspose.t[ilog]);
                }
                if (t_cz.t == NULL) t_cz = new_matrix(cz.rows, cz.cols);
                transpose_t(t_cz, cz);
                if (firstMatrix.t == NULL) firstMatrix = new_matrix(aTranspose.cols, t_cz.rows);
                timessEigen(firstMatrix, t_cz, aTranspose);
                if (verbose >= 3){
                    mxLog("firstMatrix is: \n");
                    for (int ilog = 0; ilog < firstMatrix.cols; ilog++) mxLog("%f",firstMatrix.t[ilog]);
                }
                if (verbose >= 3){
                    mxLog("yg2 is: \n");
                    for (int ilog = 0; ilog < yg2.cols; ilog++) mxLog("%f",yg2.t[ilog]);
                }
                if (solution.t == NULL) solution = new_matrix(yg2.cols, firstMatrix.cols);
                if (input.t == NULL) input = new_matrix(firstMatrix.cols, firstMatrix.rows);
                duplicateIt_t(input, firstMatrix);
                if (rhs.t == NULL) rhs = fill(yg2.cols, max(firstMatrix.cols, firstMatrix.rows), (double)0.0);
                for (int i = 0; i < yg2.rows; i++)
                    for (int j = 0; j < yg2.cols; j++)
                        M(rhs, j, i) = M(yg2, j, i);
                QRdsolve_t(solution, input, rhs);
                if (verbose >= 3){
                    mxLog("solution is: \n");
                    for (int ilog = 0; ilog < solution.cols; ilog++) mxLog("%f",solution.t[ilog]);
                }
                //Matrix solution = QRd(firstMatrix, secondMatrix);
                if (t_sol.t == NULL) t_sol = new_matrix(solution.rows, solution.cols);
                transpose_t(t_sol, solution);
                if (verbose >= 3){
                    mxLog("y is: \n");
                    for (int ilog = 0; ilog < t_sol.cols * t_sol.rows; ilog++) mxLog("%f",t_sol.t[ilog]);
                }
                //yMatrix = subset(solution, 0, 0, nc-1);
                
                //u = -cz %*% (yg - ( t(cz) %*% t(a) ) %*% y)
                multiplyByScalar2D(cz, -1.0);
                if (toSubtract.t == NULL) toSubtract = new_matrix(solution.cols, firstMatrix.rows);
                if (verbose >= 3){
                    mxLog("solution is: \n");
                    for (int ilog = 0; ilog < solution.cols; ilog++) mxLog("%f",solution.t[ilog]);
                }
                if (verbose >= 3){
                    mxLog("firstMatrix is: \n");
                    for (int ilog = 0; ilog < firstMatrix.cols; ilog++) mxLog("%f",firstMatrix.t[ilog]);
                }
                
                timessEigen(toSubtract, firstMatrix, solution);
                if (verbose >= 3){
                    mxLog("toSubtract.rows: %d", toSubtract.rows);
                    mxLog("toSubtract.cols: %d", toSubtract.cols);
                    
                    mxLog("toSubtract is: \n");
                    for (int ilog = 0; ilog < toSubtract.cols * toSubtract.rows; ilog++) mxLog("%f",toSubtract.t[ilog]);
                }
                subtractEigen(yg2, toSubtract);
                if (verbose >= 3){
                    mxLog("yg2.rows: %d", yg2.rows);
                    mxLog("yg2.cols: %d", yg2.cols);
                    mxLog("yg2 is: \n");
                    for (int ilog = 0; ilog < yg2.cols * yg2.rows; ilog++) mxLog("%f",yg2.t[ilog]);
                }
                if (uu.t == NULL) uu = new_matrix(yg2.cols, cz.rows);
                timessEigen(uu, cz, yg2);
                
                if (t_u.t == NULL) t_u = new_matrix(uu.rows, uu.cols);
                transpose_t(t_u, uu);
            }
            
            if (p0_1.t == NULL) p0_1 = new_matrix(npic, 1);
            subsetEigen(p0_1, t_u, 0, 0, npic-1);
            
            addEigen(p0_1, p);
            if (verbose >= 3){
                mxLog("p0_1 is: \n");
                for (int ilog = 0; ilog < p0_1.cols; ilog++) mxLog("%f",p0_1.t[ilog]);
            }
            
            if (M(ind, 10, 0) <= 0.5){
                go = 1;
            } else{
                if (listPartOne.t == NULL) listPartOne = new_matrix(mm, 1);
                subsetEigen(listPartOne, p0_1, 0, 0, mm-1);
                if (t16.t == NULL) t16 = new_matrix(pb.rows, 1);
                getColumn_t(t16, pb, 0);
                subtractEigen(listPartOne, t16);
                if (listPartTwo.t == NULL) listPartTwo = new_matrix(pb.rows, 1);
                getColumn_t(listPartTwo, pb, 1);
                if (t17.t == NULL) t17 = new_matrix(mm, 1);
                subsetEigen(t17, p0_1, 0, 0, mm-1);
                subtractEigen(listPartTwo, t17);
                if (llist.t == NULL) llist = new_matrix(listPartOne.cols + listPartTwo.cols, listPartOne.rows);
                copyEigen(llist, listPartOne, listPartTwo);
                
                go = findMin(llist);
                lambdaValue = 3 * lambdaValue;
                if (verbose >= 3){
                    mxLog("go is: \n");
                    mxLog("%.20f", go);
                    mxLog("lambdaValue is: \n");
                    mxLog("%.20f", lambdaValue);
                    
                }
            }
        } // end while(go <= 0){
        
        
        alp[0] = 0;
        if (ob1.t == NULL) ob1 = new_matrix(ob.cols, ob.rows);
        duplicateIt_t(ob1, ob);
        if (ob2.t == NULL) ob2 = new_matrix(ob1.cols, ob1.rows);
        duplicateIt_t(ob2, ob1);
        
        sob[0] = j;
        sob[1] = j;
        
        
        if (verbose >= 3){
            mxLog("sob is is: \n");
            for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f", sob[ilog]);
        }
        
        if (t18.t == NULL) t18 = new_matrix(p.rows, p.cols);
        transpose_t(t18, p);
        if (ptt.t == NULL) ptt = new_matrix(t18.cols + t18.cols, t18.rows);
        copyEigen(ptt, t18, t18);
        
        alp[2] = 1.0;
        if (verbose >= 3){
            mxLog("ptt is: \n");
            for (int ilog = 0; ilog < ptt.cols*ptt.rows; ilog++) mxLog("%f", ptt.t[ilog]);
        }
        
        if (p0_1_t.t == NULL) p0_1_t = new_matrix(p0_1.rows, p0_1.cols);
        transpose_t(p0_1_t, p0_1);
        if (ptt2.t == NULL) ptt2 = new_matrix(ptt.cols + p0_1_t.cols, ptt.rows);
        copyEigen(ptt2, ptt, p0_1_t);
        
        if (verbose >= 3){
            mxLog("ptt2 is: \n");
            for (int ilog = 0; ilog < ptt2.cols*ptt2.rows; ilog++) mxLog("%f",ptt2.t[ilog]);
        }
        
        if (pttCol.t == NULL) pttCol = new_matrix(ptt2.rows, 1);
        getColumn_t(pttCol, ptt2, 2);
        if (verbose >= 3){
            mxLog("pttCol is: \n");
            for (int ilog = 0; ilog < pttCol.cols; ilog++) mxLog("%f",pttCol.t[ilog]);
        }
        if (tmpv.t == NULL) tmpv = new_matrix(npic - nineq, 1);
        subsetEigen(tmpv, pttCol, 0, nineq, (npic-1));
        if (verbose >= 3){
            mxLog("tmpv is: \n");
            for (int ilog = 0; ilog < tmpv.cols; ilog++) mxLog("%f",tmpv.t[ilog]);
        }
        
        if (t19.t == NULL) t19 = new_matrix(nc+np-nc, 1);
        subsetEigen(t19, vscale, 0, (nc+1), (nc+np));
        multiplyEigen(tmpv, t19);
        
        if (verbose >= 3){
            mxLog("tmpv is: \n");
            for (int ilog = 0; ilog < tmpv.cols; ilog++) mxLog("%f",tmpv.t[ilog]);
        }
        if (verbose >= 2){
            //printf("10th call is \n");
        }
        funv = solFun(tmpv, mode, verbose);
        if (verbose >= 3){
            mxLog("hessv is: \n");
            for (int ilog = 0; ilog < hessv.cols; ilog++) mxLog("%f",hessv.t[ilog]);
            
            mxLog("g is: \n");
            for (int ilog = 0; ilog < g.cols; ilog++) mxLog("%f",g.t[ilog]);
            
            mxLog("funv is: \n");
            mxLog("%.20f", funv);
        }
        
        if (*mode == -1)
        {
            funv = 1e24;
            *mode = 0;
        }
        
        eqv = solEqBFun(verbose);
        
        ineqv = myineqFun(verbose);
        
        solnp_nfn = solnp_nfn + 1;
        
        if (M(ineqv, 0, 0) != EMPTY){
            if(M(eqv,0,0) != EMPTY)
            {
                if (t20.t == NULL) t20 = new_matrix(1, 1);
                fill_t(t20, 1, 1, funv);
                if (firstPartt.t == NULL) firstPartt = new_matrix(t20.cols + eqv.cols, t20.rows);
                copyEigen(firstPartt, t20, eqv);
                
                if (firstPart.t == NULL) firstPart = new_matrix(firstPartt.cols + ineqv.cols, firstPartt.rows);
                copyEigen(firstPart, firstPartt, ineqv);
            }
            else{
                if (t20.t == NULL) t20 = new_matrix(1, 1);
                fill_t(t20, 1, 1, funv);
                if (firstPart.t == NULL) firstPart = new_matrix(t20.cols + ineqv.cols, t20.rows);
                copyEigen(firstPart, t20, ineqv);
            }
        }
        else if (M(eqv, 0, 0) != EMPTY){
            if (t20.t == NULL) t20 = new_matrix(1, 1);
            fill_t(t20, 1, 1, funv);
            if (firstPart.t == NULL) firstPart = new_matrix(t20.cols + eqv.cols, t20.rows);
            copyEigen(firstPart, t20, eqv);
        }
        else {
            if (firstPart.t == NULL) firstPart = new_matrix(1, 1);
            fill_t(firstPart, 1, 1, funv);
            
        }
        
        
        if (secondPart.t == NULL) secondPart = new_matrix(nc+1, 1);
        subsetEigen(secondPart, vscale, 0, 0, nc);
        if (verbose >= 3){
            mxLog("firstPart is: \n");
            for (int ilog = 0; ilog < firstPart.cols * firstPart.rows; ilog++) mxLog("%f",firstPart.t[ilog]);
            mxLog("secondPart is: \n");
            for (int ilog = 0; ilog < secondPart.cols * secondPart.rows; ilog++) mxLog("%f",secondPart.t[ilog]);
        }
        divideEigen(firstPart, secondPart);
        if (ob3.t == NULL) ob3 = new_matrix(firstPart.cols, firstPart.rows);
        duplicateIt_t(ob3, firstPart);
        
        if (verbose >= 3){
            mxLog("ob3 is: \n");
            for (int ilog = 0; ilog < ob3.cols; ilog++) mxLog("%f",ob3.t[ilog]);
        }
        
        sob[2] = M(ob3, 0, 0);
        
        if (M(ind, 3, 0) > 0.5){
            // ob3[ (neq + 2):(nc + 1) ] = ob3[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq, 3 ]
            if (partOne.t == NULL) partOne = new_matrix(nc-neq, 1);
            subsetEigen(partOne, ob3, 0, neq+1, nc);
            if (tempPttCol.t == NULL) tempPttCol = new_matrix(ptt2.rows, 1);
            getColumn_t(tempPttCol, ptt2, 2);
            if (partTwo.t == NULL) partTwo = new_matrix(nineq, 1);
            subsetEigen(partTwo, tempPttCol, 0, 0, nineq-1);
            subtractEigen(partOne, partTwo);
            copyIntoInplace(ob3, partOne, 0, neq+1, nc);
        }
        
        if (nc > 0){
            //sob[ 3 ] = ob3[ 1 ] - t(yy) %*% ob3[ 2:(nc + 1) ] + rho * .vnorm(ob3[ 2:(nc + 1) ]) ^ 2
            if (firstp.t == NULL) firstp = new_matrix(nc, 1);
            subsetEigen(firstp, ob3, 0, 1, nc);
            if (t21.t == NULL) t21 = new_matrix(ptt2.rows, 1);
            getColumn_t(t21, ptt2, 2);
            if (t22.t == NULL) t22 = new_matrix(t21.rows, t21.cols);
            transpose_t(t22, t21);
            if (t23.t == NULL) t23 = new_matrix(t22.cols, a.rows);
            timessEigen(t23, a, t22);
            if (t24.t == NULL) t24 = new_matrix(t23.rows, t23.cols);
            transpose_t(t24, t23);
            subtractEigen(firstp, t24);
            addEigen(firstp, b);
            copyIntoInplace(ob3, firstp, 0, 1, nc);
            if (t25.t == NULL) t25 = new_matrix(nc, 1);
            subsetEigen(t25, ob3, 0, 1, nc);
            
            double vnormTerm = vnorm(t25) * vnorm(t25);
            if (yyTerm.t == NULL) yyTerm = new_matrix(yy.rows, yy.cols);
            transpose_t(yyTerm, yy);
            if (t26.t == NULL) t26 = new_matrix(yyTerm.cols, 1);
            getRow_t(t26, yyTerm, 0);
            if (t27.t == NULL) t27 = new_matrix(t25.cols, 1);
            getRow_t(t27, t25, 0);
            double dotProductTerm = dotProduct(t26, t27);
            
            sob[2] = M(ob3, 0, 0) - dotProductTerm + (rho * vnormTerm);
        }
        
        go = 1;
        
        while(go > tol){
            alp[1] = (alp[0] + alp[2]) / 2.0;
            if (p_copy.t == NULL) p_copy = new_matrix(p.cols, p.rows);
            duplicateIt_t(p_copy, p);
            multiplyByScalar2D(p_copy, (1.0 - alp[1]));
            if (p0_copy.t == NULL) p0_copy = new_matrix(p0_1.cols, p0_1.rows);
            duplicateIt_t(p0_copy, p0_1);
            multiplyByScalar2D(p0_copy, alp[1]);
            addEigen(p_copy, p0_copy);
            setColumnInplace(ptt2, p_copy, 1);
            if (pttColOne.t == NULL) pttColOne = new_matrix(ptt2.rows, 1);
            getColumn_t(pttColOne, ptt2, 1);
            if (tmpv.t == NULL) tmpv = new_matrix(npic-nineq, 1);
            subsetEigen(tmpv, pttColOne, 0, nineq, (npic-1));
            if (t29.t == NULL) t29 = new_matrix(nc+np-nc, 1);
            subsetEigen(t29, vscale, 0, (nc+1), (nc+np));
            multiplyEigen(tmpv, t29);
            if (verbose >= 3){
                mxLog("tmpv is: \n");
                for (int ilog = 0; ilog < tmpv.cols; ilog++) mxLog("%f",tmpv.t[ilog]);
            }
            
            if (verbose >= 3){
                mxLog("11th call is \n");
            }
            
            funv = solFun(tmpv, mode, verbose);
            if (verbose >= 3){
                mxLog("funv is: \n");
                mxLog("%2f", funv);
            }
            
            if (*mode == -1)
            {
                funv = 1e24;
                *mode = 0;
            }
            
            eqv = solEqBFun(verbose);
            
            ineqv = myineqFun(verbose);
            solnp_nfn = solnp_nfn + 1;
            Matrix firstPart, secondPart, firstPartt;
            
            if (M(ineqv, 0, 0) != EMPTY){
                if(M(eqv,0,0) != EMPTY)
                {
                    if (t20.t == NULL) t20 = new_matrix(1, 1);
                    fill_t(t20, 1, 1, funv);
                    if (firstPartt.t == NULL) firstPartt = new_matrix(t20.cols + eqv.cols, t20.rows);
                    copyEigen(firstPartt, t20, eqv);
                    
                    if (firstPart.t == NULL) firstPart = new_matrix(firstPartt.cols + ineqv.cols, firstPartt.rows);
                    copyEigen(firstPart, firstPartt, ineqv);
                }
                else{
                    if (t20.t == NULL) t20 = new_matrix(1, 1);
                    fill_t(t20, 1, 1, funv);
                    if (firstPart.t == NULL) firstPart = new_matrix(t20.cols + ineqv.cols, t20.rows);
                    copyEigen(firstPart, t20, ineqv);
                }
            }
            else if (M(eqv, 0, 0) != EMPTY){
                if (t20.t == NULL) t20 = new_matrix(1, 1);
                fill_t(t20, 1, 1, funv);
                if (firstPart.t == NULL) firstPart = new_matrix(t20.cols + eqv.cols, t20.rows);
                copyEigen(firstPart, t20, eqv);
            }
            else {
                if (firstPart.t == NULL) firstPart = new_matrix(1, 1);
                fill_t(firstPart, 1, 1, funv);
                
            }
            if (secondPart.t == NULL) secondPart = new_matrix(nc+1, 1);
            subsetEigen(secondPart, vscale, 0, 0, nc);
            divideEigen(firstPart, secondPart);
            duplicateIt_t(ob2, firstPart);
            
            if (verbose >= 3){
                mxLog("ob2 is: \n");
                for (int ilog = 0; ilog < ob2.cols; ilog++) mxLog("%f",ob2.t[ilog]);
            }
            
            sob[1] = M(ob2, 0, 0);
            if (verbose >= 3){
                mxLog("sob is: \n");
                for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
            }
            if (M(ind, 3, 0) > 0.5){
                if (partOne.t == NULL) partOne = new_matrix(nc-neq, 1);
                subsetEigen(partOne, ob2, 0, neq+1, nc);
                if (tempPttCol.t == NULL) tempPttCol = new_matrix(ptt2.rows, 1);
                getColumn_t(tempPttCol, ptt2, 1);
                if (partTwo.t == NULL) partTwo = new_matrix(nineq, 1);
                subsetEigen(partTwo, tempPttCol, 0, 0, nineq-1);
                subtractEigen(partOne, partTwo);
                copyIntoInplace(ob2, partOne, 0, neq+1, nc);
            }
            
            if (nc > 0){
                if (t30.t == NULL) t30 = new_matrix(nc, 1);
                subsetEigen(t30, ob2, 0, 1, nc);
                if (t31.t == NULL) t31 = new_matrix(ptt2.rows, 1);
                getColumn_t(t31, ptt2, 1);
                if (t32.t == NULL) t32 = new_matrix(t31.rows, t31.cols);
                transpose_t(t32, t31);
                if (t33.t == NULL) t33 = new_matrix(t32.cols, a.rows);
                timessEigen(t33, a, t32);
                if (t34.t == NULL) t34 = new_matrix(t33.rows, t33.cols);
                transpose_t(t34, t33);
                subtractEigen(t30, t34);
                addEigen(t30, b);
                copyIntoInplace(ob2, t30, 0, 1, nc);
                if (temp.t == NULL) temp = new_matrix(nc, 1);
                subsetEigen(temp, ob2, 0, 1, nc);
                double vnormTerm = vnorm(temp) * vnorm(temp);
                if (yyTerm.t == NULL) yyTerm = new_matrix(yy.rows, yy.cols);
                transpose_t(yyTerm, yy);
                if (t26.t == NULL) t26 = new_matrix(yyTerm.cols, 1);
                getRow_t(t26, yyTerm, 0);
                if (t27.t == NULL) t27 = new_matrix(t25.cols, 1);
                getRow_t(t27, temp, 0);
                double dotProductTerm = dotProduct(t26, t27);
                sob[1] = M(ob2, 0, 0) - dotProductTerm + rho * vnormTerm;
            }
            if (verbose >= 3){
                mxLog("sob is: \n");
                for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
            }
            if (verbose >= 3){
                mxLog("obm is: \n");
                for (int ilog = 0; ilog < obm.cols; ilog++) mxLog("%f",obm.t[ilog]);
            }
            if (obm.t == NULL) obm = new_matrix(1, 1);
            M(obm, 0, 0) = sob.maxCoeff();
            if (verbose >= 3){
                mxLog("obm is: \n");
                for (int ilog = 0; ilog < obm.cols; ilog++) mxLog("%f",obm.t[ilog]);
            }
            if (M(obm, 0, 0) < j){
                double obn = sob.minCoeff();
                go = tol * (M(obm, 0, 0) - obn) / (j - M(obm, 0, 0));
            }
            
            const bool condif1 = (sob[1] >= sob[0]);
            const bool condif2 = (sob[0] <= sob[2]) && (sob[1] < sob[0]);
            const bool condif3 = (sob[1] <  sob[0]) && (sob[0] > sob[2]);
            
            if (condif1){
                sob[2] = sob[1];
                if (ob3.t == NULL) ob3 = new_matrix(ob2.cols, ob2.rows);
                duplicateIt_t(ob3, ob2);
                alp[2] = alp[1];
                if (tempCol.t == NULL) tempCol = new_matrix(ptt2.rows, 1);
                getColumn_t(tempCol, ptt2, 1);
                setColumnInplace(ptt2, tempCol, 2);
                
                if (verbose >= 3){
                    mxLog("sob is: \n");
                    for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
                    mxLog("ob3 is: \n");
                    for (int ilog = 0; ilog < ob3.cols; ilog++) mxLog("%f",ob3.t[ilog]);
                    mxLog("alp is: \n");
                    for (int ilog = 0; ilog < alp.size(); ilog++) mxLog("%f",alp[ilog]);
                    mxLog("ptt2 is: \n");
                    for (int ilog = 0; ilog < ptt2.cols; ilog++) mxLog("%f",ptt2.t[ilog]);
                }
            }
            
            if (condif2){
                sob[2] = sob[1];
                if (ob3.t == NULL) ob3 = new_matrix(ob2.cols, ob2.rows);
                duplicateIt_t(ob3, ob2);
                alp[2] = alp[1];
                if (tempCol.t == NULL) tempCol = new_matrix(ptt2.rows, 1);
                getColumn_t(tempCol, ptt2, 1);
                setColumnInplace(ptt2, tempCol, 2);
                
                if (verbose >= 3){
                    mxLog("sob is: \n");
                    for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
                    mxLog("ob3 is: \n");
                    for (int ilog = 0; ilog < ob3.cols; ilog++) mxLog("%f",ob3.t[ilog]);
                    mxLog("alp is: \n");
                    for (int ilog = 0; ilog < alp.size(); ilog++) mxLog("%f",alp[ilog]);
                    mxLog("ptt2 is: \n");
                    for (int ilog = 0; ilog < ptt2.cols; ilog++) mxLog("%f",ptt2.t[ilog]);				}
            }
            
            if (condif3){
                sob[0] = sob[1];
                if (ob1.t == NULL) ob1 = new_matrix(ob2.cols, ob2.rows);
                duplicateIt_t(ob1, ob2);
                alp[0] = alp[1];
                if (tempCol.t == NULL) tempCol = new_matrix(ptt2.rows, 1);
                getColumn_t(tempCol, ptt2, 1);
                setColumnInplace(ptt2, tempCol, 0);
                if (verbose >= 3){
                    mxLog("sob is: \n");
                    for (int ilog = 0; ilog < sob.size(); ilog++) mxLog("%f",sob[ilog]);
                    mxLog("ob3 is: \n");
                    for (int ilog = 0; ilog < ob3.cols; ilog++) mxLog("%f",ob3.t[ilog]);
                    mxLog("alp is: \n");
                    for (int ilog = 0; ilog < alp.size(); ilog++) mxLog("%f",alp[ilog]);
                    mxLog("ptt2 is: \n");
                    for (int ilog = 0; ilog < ptt2.cols; ilog++) mxLog("%f",ptt2.t[ilog]);				}
            }
            
            if (go >= tol){
                go = alp[2] - alp[0];
                if (verbose >= 3){
                    mxLog("go is: \n");
                    mxLog("%.20f", go);
                }
            }
            
        } // 	while(go > tol){
        
        if (verbose >= 3){
            mxLog("go is: \n");
            mxLog("%.16f", go);
        }
        //mxLog("sx_Matrix is: \n");
        //for (i = 0; i < sx_Matrix.cols * sx_Matrix.rows; i++) mxLog("%f",sx_Matrix.t[i]);
        //sx_Matrix = duplicateIt(sx);
        // mxLog("sx is:");
        //for (int ilog = 0; ilog < sx.cols * sx.rows; ilog++) mxLog("%f",sx.t[ilog]);
        //mxLog("sx.rows: %d", sx.rows);
        //mxLog("sx.cols: %d", sx.cols);
        if (sx.t == NULL) sx = fill(p.cols, p.rows, (double)0.0);
        duplicateIt_t(sx_Matrix, sx);
        duplicateIt_t(sx, p);
        if (yg.t == NULL) yg = new_matrix(g.cols, g.rows);
        duplicateIt_t(yg, g);
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
            mxLog("%.20f", obn);
        }
        if (j <= obn){
            maxit = minit;
        }
        if (verbose >= 3){
            mxLog("j is: \n");
            mxLog("%.20f", j);
        }
        double reduce = (j - obn) / ((double)1.0 + (double)fabs(j));
        if (verbose >= 3){
            mxLog("reduce is: \n");
            mxLog("%.20f", reduce);
        }
        if (reduce < tol){
            maxit = minit;
        }
        
        const bool condif1 = (sob[0] <  sob[1]);
        const bool condif2 = (sob[2] <  sob[1]) && (sob[0] >=  sob[1]);
        const bool condif3 = (sob[0] >= sob[1]) && (sob[2] >= sob[1]);
        
        if (condif1){
            j = sob[0];
            if (p.t == NULL) p = new_matrix(ptt2.rows, 1);
            getColumn_t(p, ptt2, 0);
            if (ob.t == NULL) ob = new_matrix(ob1.cols, ob1.rows);
            duplicateIt_t(ob, ob1);
            if (verbose >= 3){
                mxLog("condif1\n");
                mxLog("j is: \n");
                mxLog("%2f", j);
                mxLog("p is: \n");
                for (int ilog = 0; ilog < p.cols; ilog++) mxLog("%f",p.t[ilog]);
                mxLog("ob is: \n");
                for (int ilog = 0; ilog < ob.cols; ilog++) mxLog("%f",ob.t[ilog]);
            }
        }
        
        if (condif2){
            
            j = sob[2];
            if (p.t == NULL) p = new_matrix(ptt2.rows, 1);
            getColumn_t(p, ptt2, 2);
            if (ob.t == NULL) ob = new_matrix(ob3.cols, ob3.rows);
            duplicateIt_t(ob, ob3);
            if (verbose >= 3){
                mxLog("condif2\n");
                mxLog("j is: \n");
                mxLog("%2f", j);
                mxLog("p is: \n");
                for (int ilog = 0; ilog < p.cols; ilog++) mxLog("%f",p.t[ilog]);
                mxLog("ob is: \n");
                for (int ilog = 0; ilog < ob.cols; ilog++) mxLog("%f",ob.t[ilog]);
            }
            
        }
        
        if (condif3){
            j = sob[1];
            if (p.t == NULL) p = new_matrix(ptt2.rows, 1);
            getColumn_t(p, ptt2, 1);
            if (ob.t == NULL) ob = new_matrix(ob2.cols, ob2.rows);
            duplicateIt_t(ob, ob2);
            if (verbose >= 3){
                mxLog("condif3\n");
                mxLog("j is: \n");
                mxLog("%2f", j);
                mxLog("p is: \n");
                for (int ilog = 0; ilog < p.cols; ilog++) mxLog("%f",p.t[ilog]);
                mxLog("ob is: \n");
                for (int ilog = 0; ilog < ob.cols; ilog++) mxLog("%f",ob.t[ilog]);
            }
        }
        if (verbose >= 3){
            mxLog("yg\n");
            for (int ilog = 0; ilog < yg.cols; ilog++) mxLog("%f",yg.t[ilog]);
        }
        //printSize();
    } // end while (minit < maxit){
    
    M(yg_rec, 1, 0) = vnorm(yg);
    if(M(yg_rec, 0, 0) / M(yg_rec, 1, 0) > 1000)  flag_NormgZ = 1;
    
    minr_rec = minit;
    Matrix result2 = getColumn(ptt2, 1);
    subtractEigen(result2, getColumn(ptt2, 0));
    Matrix result3 = getColumn(ptt2, 1);
    subtractEigen(result3, getColumn(ptt2, 2));
    if (all(result2) || all(result3)) flag_step = 1;
    //p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
    Matrix vscalePart = subset(vscale, 0, (neq+1), (nc+np));
    // I need vscale, p, y, hessv
    multiplyEigen(p, vscalePart);
    
    if (nc > 0){
        //y = vscale[ 0 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
        multiplyByScalar2D(t_sol, M(vscale,0,0));
        divideEigen(t_sol, subset(vscale, 0, 1, nc));
    }
    
    // hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*%
    //                                t(vscale[ (neq + 2):(nc + np + 1) ]) )
    
    Matrix transposePart = transpose2D(subset(vscale, 0, (neq+1), (nc+np)));
    divideEigen(hessv, transposePart);
    multiplyByScalar2D(hessv, M(vscale,0,0));
    
    if (verbose >= 1 && reduce > tol) {
        mxLog("m3 solnp Rf_error message being reported.");
    }
    
    resP = duplicateIt(p);
    resY = transpose(subset(t_sol, 0, 0, (yyRows-1)));
    resHessv = duplicateIt(hessv);
    resLambda = lambdaValue;
    
    if (verbose >= 3){
        mxLog("------------------------RETURNING FROM SUBNP------------------------");
        mxLog("p information: ");
        for (int ilog = 0; ilog < resP.cols; ilog++) mxLog("%f",resP.t[ilog]);
        mxLog("y information: ");
        for (int ilog = 0; ilog < resY.cols; ilog++) mxLog("%f",resY.t[ilog]);
        mxLog("hessv information: ");
        for (int ilog = 0; ilog < resHessv.cols; ilog++) mxLog("%f",resHessv.t[ilog]);
        mxLog("lambda information: ");
        mxLog("%f", resLambda);
        mxLog("minit information: ");
        mxLog("%d", minit);
        mxLog("------------------------END RETURN FROM SUBNP------------------------");
    }
    
    return g;
    
} // end subnp
