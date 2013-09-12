#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>

#include "matrix.h"
#include "omxCsolnp.h"


double EMPTY;

bool DEBUG;

Matrix ineqLB;
Matrix ineqUB;

Matrix LB;

Matrix UB;

Matrix control;

Matrix ind;

Matrix pars;

Matrix eqB;

Matrix resP;
double resLambda;
Matrix resMu;
Matrix resHessv;
Matrix resY;

int ineqLBLength;
int ineqUBLength;
int LBLength;
int UBLength;
int parsLength;
int eqBLength;
double eps;
int outerIter;

Matrix subnp(Matrix pars, solFun_t solFun, solEqBFun_t solEqBFun, solIneqFun_t myineqFun,
	     Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,  Matrix ctrl, int verbose);

Param_Obj solnp(Matrix solPars, solFun_t solFun,
		Matrix solEqB, solEqBFun_t solEqBFun, solIneqFun_t myineqFun,
		Matrix solLB,  Matrix solUB,  Matrix solIneqUB,  Matrix solIneqLB,
		Matrix solctrl, bool debugToggle, int verbose)
{
	if(verbose >= 3){
		printf("solPars is: \n");
		print(solPars); putchar('\n');
		printf("4th call is: \n");
		printf("%2f", solFun(solPars, 0)); putchar('\n');
		printf("solEqB is: \n");
		print(solEqB); putchar('\n');
		printf("solEqBFun is: \n");
		print(solEqBFun(solPars, 0)); putchar('\n');
		printf("myineqFun is: \n");
		print(myineqFun(solPars, 0)); putchar('\n');
		printf("solLB is: \n");
		print(solLB); putchar('\n');
		printf("solUB is: \n");
		print(solUB); putchar('\n');
		printf("solIneqUB is: \n");
		print(solIneqUB); putchar('\n');
		printf("solIneqLB is: \n");
		print(solIneqLB); putchar('\n');
	}

	double solnp_nfn = 0;
	eps = 2.220446e-16;
	time_t sec;
	sec = time (NULL);
	ind = fill(11, 1, (double) 0.0);
	DEBUG = debugToggle;
	EMPTY = -999999.0;
	
	ineqLBLength = solIneqLB.cols;
	ineqUBLength = solIneqUB.cols;
	LBLength = solLB.cols;
	UBLength = solUB.cols;
	parsLength = solPars.cols;
	eqBLength = solEqB.cols;
	ineqLB.cols = solIneqLB.cols;
	ineqUB.cols = solIneqUB.cols;
    
	Matrix grad = fill(solPars.cols, 1, (double)0.0);
    
    
	pars = duplicateIt(solPars);
    
	eqB = duplicateIt(solEqB);
	
	control = duplicateIt(solctrl);

	if(verbose >= 2){
		printf("control is: \n");
		print(control); putchar('\n');
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
    
	if (LBLength > 1){
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
	
	if (UBLength > 1){
		UB = duplicateIt(solUB);
		if (LB.cols < 1){
			LB = fill(UB.cols, 1, (double) -DBL_MAX/2);
		}
        
	}
	else{
		UB = fill(1, 1, (double) 0.0);
		M(UB, 0, 0) = EMPTY;
	}
    
	if (verbose >= 2){
		printf("LB is: \n");
		print(LB); putchar('\n');

		printf("UB is: \n");
		print(UB); putchar('\n');
    	}
    
   
	int np = pars.cols;
    
	// [0] length of pars
	// [1] has function gradient?
	// [2] has hessian?
	// [3] has ineq?
	// [4] ineq length
	// [5] has jacobian (inequality)
	// [6] has eq?
	// [7] eq length
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
	double funv = solFun(pars, verbose);
    
    
	// does not have a hessian (currently not supported in Rsolnp)
	M(ind, 2, 0) = 0;
	
	// do inequality checks and return starting values
	int nineq;
	Matrix ineqx0 = fill(ineqLB.cols, 1, (double)0.0);
    
	Matrix ineqv = myineqFun(pars, verbose);
    
	if ( M(ineqv, 0, 0) != EMPTY){
		
		M(ind, 3, 0) = 1;
		nineq = ineqLB.cols;
		
		M(ind, 4, 0) = nineq;
		
		// check for infitnites/nans
        
		Matrix ineqLBx = ineqLB;
		Matrix ineqUBx = ineqUB;
        
		int i;
		for (i = 0; i<ineqLBx.cols; i++)
			{
				if (M(ineqLBx,i,0) <= -99999999.0){ //-99999999.0
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
	Matrix eqv = solEqBFun(pars, verbose);
    
	if( M(eqv, 0, 0) != EMPTY){
		M(ind, 6, 0) = 1;
		neq = eqB.cols;
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
		printf("ind is: \n");
		print(ind);
	}

	Matrix pb;
    
    
	if(M(ind, 10, 0))
		{   if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) != EMPTY)) 
				{   pb = fill(2, nineq, (double)0.0);
					pb = setColumn(pb, ineqLB, 0);
					pb = setColumn(pb, ineqUB, 1);
					Matrix pb_cont = fill(2, np, (double)0.0);
					pb_cont = setColumn(pb_cont, LB, 0);
					pb_cont = setColumn(pb_cont, UB, 1);
					pb = transpose(copy(transpose(pb), transpose(pb_cont)));
				}
			else if((M(LB, 0, 0) == EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
				{
					pb = fill(2, nineq, (double)0.0);
					pb = setColumn(pb, ineqLB, 0);
					pb = setColumn(pb, ineqUB, 1);
				}
			else if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) == EMPTY))
				{
					pb = fill(2, np, (double)0.0);
					pb = setColumn(pb, LB, 0);
					pb = setColumn(pb, UB, 1);
				}
		}
    
    
	//  if (M(LB, 0, 0) != EMPTY)
	else{
		pb = fill(1, 1, EMPTY);
	}
    
        int pbRows = pb.rows;
	int pbCols = pb.cols;
	double rho   = M(control, 0, 0);
	int maxit = M(control, 1, 0);
	int minit = M(control, 2, 0);
	double delta = M(control, 3, 0);
	double tol   = M(control, 4, 0);
	double trace = M(control, 5, 0);
	
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
		else{
			constraint = duplicateIt(eqv);
		}
        
		if( M(ind, 3, 0) > 0 ) {
			
			// 	tmpv = cbind(constraint[ (neq[0]):(tc[0]-1) ] - .ineqLB, .ineqUB - constraint[ (neq + 1):tc ] )
			Matrix difference1 = subtract(subset(constraint, 0, neq, tc-1), ineqLB);
			Matrix difference2 = subtract(ineqUB, subset(constraint, 0, neq, tc-1));
			Matrix tmpv = fill(2, nineq, (double)0.0);
			tmpv = setColumn(tmpv, difference1, 0);
			tmpv = setColumn(tmpv, difference2, 1);
			Matrix testMin = rowWiseMin(tmpv);
            
			if( allGreaterThan(testMin, 0) ) {
				ineqx0 = subset(constraint, 0, neq, tc-1);
			}
            
			constraint = copyInto(constraint, subtract(subset(constraint, 0, neq, tc-1), ineqx0), 0, neq, tc-1);
            
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
	else ob = funvMatrix;
    
	if(verbose >= 3){
		printf("ob is: \n");
		print(ob); putchar('\n');
	}

	Matrix vscale;
    
	while(solnp_iter < 2){
		solnp_iter = solnp_iter + 1;
		outerIter = solnp_iter;
		Matrix subnp_ctrl = fill(5, 1, (double)0.0);
		M(subnp_ctrl, 0, 0) = rho;
		M(subnp_ctrl, 1, 0) = minit;
		M(subnp_ctrl, 2, 0) = delta;
		M(subnp_ctrl, 3, 0) = tol;
		M(subnp_ctrl, 4, 0) = trace;
        
		if ( M(ind, 6, 0) > 0){
			Matrix subsetMat = subset(ob, 0, 1, neq);
			double max = findMax(matrixAbs(subsetMat));
            
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
		vscale = minMaxAbs(vscale, tol);
        
		if (verbose >= 1){
			printf("------------------------CALLING SUBNP------------------------"); putchar('\n');
			printf("p information: "); putchar('\n');
			print(p); putchar('\n');
			printf("lambda information: "); putchar('\n');
			print(lambda); putchar('\n');
			printf("ob information: "); putchar('\n');
			print(ob); putchar('\n');
			printf("hessv information: "); putchar('\n');
			print(hessv); putchar('\n');
			printf("mu information: "); putchar('\n');
			printf("%2d", mu); putchar('\n');
			printf("vscale information: "); putchar('\n');
			print(vscale); putchar('\n');
			printf("subnp_ctrl information: "); putchar('\n');
			print(subnp_ctrl); putchar('\n');
			printf("------------------------END CALLING SUBNP------------------------"); putchar('\n');
		}
        
		grad = subnp(p,solFun, solEqBFun, myineqFun, lambda, ob, hessv, mu, vscale, subnp_ctrl, verbose);
        
		p = duplicateIt(resP);
        
		lambda = duplicateIt(resY);
        
		hessv = duplicateIt(resHessv);
        
		mu = resLambda;
        
        
		Matrix temp = subset(p, 0, nineq, (nineq+np-1));
        
        if (verbose >= 2){

        printf("6th call is \n");
        }
        funv = solFun(temp, verbose);
        
		solnp_nfn = solnp_nfn + 1;
        
		Matrix funv_mat = fill(1, 1, funv);
		Matrix tempdf = copy(temp, funv_mat);
        
		eqv = solEqBFun(temp, verbose);
        
		ineqv = myineqFun(temp, verbose);
        
        
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
        
		double resultForTT = (j - M(ob, 0, 0)) / max(fabs(M(ob, 0, 0)), 1.0);
		M(tt, 0, 0) = resultForTT;
		j = M(ob, 0, 0);
        
		if (tc > 0){
			// constraint = ob[ 2:(tc + 1) ]
			constraint = subset(ob, 0, 1, tc);
            
			if ( M(ind, 3, 0) > 0.5){
				//tempv = rbind( constraint[ (neq + 1):tc ] - pb[ 1:nineq, 1 ], pb[ 1:nineq, 2 ] - constraint[ (neq + 1):tc ] )
				Matrix subsetOne = subset(constraint, 0, neq, tc-1);
				Matrix subsetTwo = subset(getColumn(pb, 0), 0, 0, nineq-1);
				Matrix subsetThree = subset(getColumn(pb, 1), 0, 0, nineq-1);
				Matrix diff1 = subtract(subsetOne, subsetTwo);
				Matrix diff2 = subtract(subsetThree, subsetOne);
				Matrix tempv = fill(nineq, 2, (double)0.0);
				tempv = setRow(tempv, 0, diff1);
				tempv = setRow(tempv, 1, diff2);
                
				if (findMin(tempv) > 0){
					Matrix copyValues = subset(constraint, 0, neq, tc-1);
					p = copyInto(p, copyValues, 0, 0, nineq-1);
				}
			} // end if (ind[0][3] > 0.5){
            
			Matrix diff = subtract(subset(constraint, 0, neq, tc-1),
					       subset(p, 0, 0, nineq-1));
            
			constraint = copyInto(constraint, diff, 0, neq, tc-1);
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
            
			Matrix list = fill(2, 1, (double)0.0);
            
			M(list, 0, 0) = tol + M(tt, 0, 0);
			M(list, 1, 0) = M(tt, 0, 1) - M(tt, 0, 2);
            
			if (findMax(list) <= 0){
				//hessv = diag( diag ( hessv ) )
				/** DOESN'T AFFECT US NOW EVENTUALLY IT WILL **/
				lambda = fill(1, 1, (double)0.0);
				hessv = diag(diag(hessv));
			}
            
			M(tt, 0, 1) = M(tt, 0, 2);
            
		} // end if (tc > 0){
        
		double vnormValue;
		Matrix tempTTVals = fill(2, 1, (double)0.0);
		M(tempTTVals, 0, 0) = M(tt, 0, 0);
		M(tempTTVals, 1, 0) = M(tt, 0, 1);
        
		vnormValue = vnorm(tempTTVals);
        
		if (vnormValue <= tol){
			maxit = solnp_iter;
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
		Matrix tempTTVals = fill(2, 1, (double) 0.0);
		M(tempTTVals, 0, 0) = M(tt, 0, 0);
		M(tempTTVals, 1, 0) = M(tt, 0, 1);
		vnormValue = vnorm(tempTTVals);
        
		if (verbose >= 1) {
			if (vnormValue <= tol){
				mxLog("The solution converged in %d iterations. It is:", solnp_iter);
				print(p);
			}
			else{
				mxLog("Solution failed to converge. Final parameters are:");
				print(p);
			}
		}
	}
    
	struct Param_Obj pfunv;
	Matrix hessi = fill(hessv.cols*hessv.cols,1 , (double)0.0);
    
	int i, ii;
	int ind_hess = 0;
	for (i = 0; i < hessv.cols; i++)
		{
			for (ii = 0; ii < hessv.rows; ii++)
				{
					M(hessi, ind_hess, 0) = M(hessv, i, ii);
					ind_hess = ind_hess + 1;
				}
		}
    

	Matrix p_hess = copy(p, hessi);
	pfunv.parameter = copy(p_hess, grad);
	pfunv.objValue = funv;
	return pfunv;
    
}

Matrix subnp(Matrix pars, solFun_t solFun, solEqBFun_t solEqBFun, solIneqFun_t myineqFun,
	     Matrix yy,  Matrix ob,  Matrix hessv, double lambda,  Matrix vscale,  Matrix ctrl, int verbose)
{
	
	int yyRows = yy.rows;
	int yyCols = yy.cols;
    
	double rho   = M(ctrl, 0, 0);
	int maxit = M(ctrl, 1, 0);
	double delta = M(ctrl, 2, 0);
	double tol =   M(ctrl, 3, 0);
	double trace =  M(ctrl, 4, 0);
    
	int neq =  (int)M(ind, 7, 0);
	int nineq = (int)M(ind, 4, 0);
	int np = (int)M(ind, 0, 0);

	double ch = 1;
	Matrix argum;
	Matrix y;

	if (verbose >= 2){
		printf("ind inside subnp is: \n");
		print(ind); putchar('\n');
	}

	Matrix alp = fill(3, 1, (double)0.0);
    
	int nc = neq + nineq;
	int npic = np + nineq;
	
	Matrix p0 = duplicateIt(pars);
	Matrix pb;
	Matrix col1_pb;
    
	if(M(ind, 10, 0))
		{
			if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
				{   pb = fill(2, nineq, (double)0.0);
					pb = setColumn(pb, ineqLB, 0);
					pb = setColumn(pb, ineqUB, 1);
            
					Matrix pb_cont = fill(2, np, (double)0.0);
					pb_cont = setColumn(pb_cont, LB, 0);
					pb_cont = setColumn(pb_cont, UB, 1);
            
					pb = transpose(copy(transpose(pb), transpose(pb_cont)));
            
				}
			else if((M(LB, 0, 0) == EMPTY) && (M(ineqLB, 0, 0) != EMPTY))
				{
					pb = fill(2, nineq, (double)0.0);
					pb = setColumn(pb, ineqLB, 0);
					pb = setColumn(pb, ineqUB, 1);
            
				}
			else if((M(LB, 0, 0) != EMPTY) && (M(ineqLB, 0, 0) == EMPTY))
				{
					pb = fill(2, np, (double)0.0);
					pb = setColumn(pb, LB, 0);
					pb = setColumn(pb, UB, 1);
            
				}
		}
	else{
		pb = fill(1,1,EMPTY);
	}
    
	if (verbose >= 3){
		printf("pb is: \n");
		print(pb); putchar('\n');
	}

	Matrix sob = fill(3, 1, (double)0.0);
	Matrix ptt;
	
	Matrix yyMatrix = duplicateIt(yy);
    
	ob = divide(ob, subset(vscale, 0, 0, nc));
	p0 = divide(p0, subset(vscale, 0, (neq+1), (nc + np)));
    
	int mm;
	if (M(ind, 10, 0) > 0){
		if (M(ind, 9, 0) <= 0){
			mm = nineq;
		}
		else{
			mm=npic;
		}
		Matrix vscaleSubset = subset(vscale, 0, neq+1, neq+mm);
		double vscaleSubsetLength = (neq+mm) - (neq+1) + 1;
		Matrix vscaleTwice = fill(pb.cols, pb.rows, (double)0.0);
		vscaleTwice = setColumn(vscaleTwice, vscaleSubset, 0);
		vscaleTwice = setColumn(vscaleTwice, vscaleSubset, 1);
        
		if (M(pb, 0, 0) != EMPTY){
			pb = divide(pb, vscaleTwice);
		}
	} // end if (ind [0][10] > 0)
    
	if (verbose >= 3){
		printf("pb is: \n");
		print(pb); putchar('\n');
	}

	// scale the lagrange multipliers and the Hessian
	if( nc > 0) {
		// yy [total constraints = nineq + neq]
		// scale here is [tc] and dot multiplied by yy
		//yy = vscale[ 2:(nc + 1) ] * yy / vscale[ 1 ]
		yy = divideByScalar2D(yy, M(vscale,0,0));
		yy = multiply(transpose(subset(vscale, 0, 1, nc)), yy);
	}
    
	// hessv [ (np+nineq) x (np+nineq) ]
	// hessv = hessv * (vscale[ (neq + 2):(nc + np + 1) ] %*% t(vscale[ (neq + 2):(nc + np + 1)]) ) / vscale[ 1 ]
    
	Matrix vscaleSubset = subset(vscale, 0, (neq+1), (nc + np));
	Matrix transDotProduct = transposeDP(vscaleSubset);//Mahsa: THIS IS WRONG. WRONG. TOTALLY WRONG. why transposeDP? let's try transposeDotProduct.
	hessv = divideByScalar2D(multiply(hessv, transDotProduct), M(vscale, 0, 0));
    
	double j = M(ob, 0, 0);

	Matrix a;
    
	if( M(ind, 3, 0) > 0){
		if ( M(ind, 6, 0) <= 0){
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
	double go, reduce;
	int minit;
	double lambdaValue = lambda;
    
    
	if (nc > 0) {
		constraint = subset(ob, 0, 1, nc);
        
		int i;
        
		for (i=0; i<np; i++){
			int index = nineq + i;
			M(p0, index, 0) = M(p0, index, 0) + delta;
			tmpv = multiply(subset(p0, 0, nineq, (npic-1)), subset(vscale, 0, (nc+1), (nc+np)));
            
			if (verbose >= 2){
				printf("7th call is \n");
			}
			funv = solFun(tmpv, verbose);
            
			eqv = solEqBFun(tmpv, verbose);
            
			ineqv = myineqFun(tmpv, verbose);
            
			//exit(0);
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
            
			ob = divide(firstPart, secondPart);
            
			M(g, index, 0) = (M(ob, 0, 0)-j) / delta;

			if (verbose >= 3){
				printf("g is: \n");
				print(g); putchar('\n');
			}
			Matrix colValues = subtract(subset(ob, 0, 1, nc), constraint);
			colValues = divideByScalar2D(colValues, delta);
			a = setColumn(a, colValues, index);
            
			M(p0, index, 0) = M(p0, index, 0) - delta;
		} // end for (int i=0; i<np, i++){
        
		if(M(ind, 3, 0) > 0){
			//constraint[ (neq + 1):(neq + nineq) ] = constraint[ (neq + 1):(neq + nineq) ] - p0[ 1:nineq ]
			Matrix firstPart, secondPart;
			firstPart  = subset(constraint, 0, neq, (neq+nineq-1));
			secondPart = subset(p0, 0, 0, (nineq-1));
			Matrix values = subtract(firstPart, secondPart);
            
			constraint = copyInto(constraint, values, 0, neq, (neq+nineq-1));
            
		}
        
		b = fill(nc, 1, (double)0.0);
        
        Matrix timess_a_p0 = timess(a, transpose(p0));
        //  b [nc,1]
		b = subtract(transpose(timess_a_p0), constraint);

        
		ch = -1;
		M(alp, 0, 0) = tol - findMax(matrixAbs(constraint));
		if ( M(alp, 0, 0) <= 0){
            
			ch = 1;
            
			if ( M(ind, 10, 0) < 0.5){
				Matrix dotProd = transposeDotProduct(a); //Mahsa: this is equal to "a %*% t(a)"
				Matrix solution = solve(dotProd, constraint);
				p0 = subtract(p0, matrixDotProduct(transpose(a), solution));
				M(alp, 0, 0) = 1;
			}
            
		} // end if (alp[0][0] <= 0){
        
        
		if (M(alp, 0, 0) <= 0){
            
			int npic_int = npic;
			p0 = copy(p0, fill(1, 1, (double)1.0));
			a = copy(a, transpose(multiplyByScalar2D(constraint, -1.0)));
            
			Matrix cx = copy(fill(npic, 1, (double)0.0), fill(1, 1, (double)1.0));
            
			dx = fill(1, npic+1, (double)1.0);
            
			go = 1;
			minit = 0;
            
			while(go >= tol)
				{
					minit = minit + 1;
					gap = fill(2, mm, (double)0.0);
					gap = setColumn(gap, subtract(subset(p0, 0, 0, mm-1), getColumn(pb, 0)), 0);
					gap = setColumn(gap, subtract(getColumn(pb, 1), subset(p0, 0, 0, mm-1)), 1);
					gap = rowSort(gap);
					dx = copyInto(transpose(dx), getColumn(gap,0), 0, 0, mm-1);
                
					M(dx, npic_int, 0) = M(p0, npic_int, 0);
                
					if(M(ind, 9, 0) <= 0)
						{
							argum = multiplyByScalar2D(fill(1, npic-mm, (double)1.0) , max(findMax(subset(dx, 0, 0, mm-1)), 100));
                    
							dx = copyInto(dx, argum, 0, mm, npic-1);
                    
						}
                
					 y = qrSolve(transpose(timess(a, transpose(diag(dx)))) , transpose(multiply(dx, transpose(cx))));
                
					//Matrix y = QRd(transpose(timess(a, transpose(diag(dx)))) , transpose(multiply(dx, transpose(cx))));
					y = subset(y, 0, 0, nc - 1);
                
					Matrix v = multiply(dx, multiply(dx, subtract(transpose(cx),timess(transpose(a),y))));
                
					int indexx = npic;
					int i;
                
                
					if (M(v, indexx, 0) > 0)
						{
							double z = M(p0, indexx, 0)/M(v, indexx, 0);
                    
							for (i=0; i<mm; i++)
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
                    
							if(z >= (M(p0, indexx, 0)/M(v, indexx, 0)))
								{
									p0 = subtract(p0, multiplyByScalar2D(v, z));
                        
								}
							else{
								p0 = subtract(p0, multiplyByScalar2D(v, 0.9 * z));
                        
							}
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
            
			//exit(0);
			if (minit >= 10){
				printf("m2 solnp error message being reported.");
				putchar('\n');
			}
			int h;
			Matrix aMatrix = fill(npic, nc, (double)0.0);
			for (h = 0; h<a.rows; h++)
				{
					aMatrix = setRow(aMatrix, h, subset(getRow(a, h), 0, 0, npic-1));
				}
			a = aMatrix;
            
			b = timess(a, transpose(subset(p0, 0, 0, npic-1)));

            
		}// end if(M(alp, 0, 0) <= 0)
	} // end if (nc > 0){
    
	p = subset(p0, 0, 0, npic-1);
    y = fill(1, 1, (double)0.0);
    
	if (verbose >= 3){
		printf("p is: \n");
		print(p); putchar('\n');
	}
    
	if (ch > 0){
		tmpv = multiply(subset(p, 0, nineq, (npic-1)), subset(vscale, 0, (nc+1), (nc+np)));
		if (verbose >= 2){
			printf("tmpv is: \n");
			print(tmpv); putchar('\n');
			printf("8th call is \n");
		}
		funv = solFun(tmpv, verbose);
		if (verbose >= 3){
			printf("funv is: \n");
			printf("%2f", funv); putchar('\n');
		}
		eqv = solEqBFun(tmpv, verbose);
		if (verbose >= 3){
			printf("eqv is: \n");
			print(eqv); putchar('\n');
		}
		ineqv = myineqFun(tmpv, verbose);
		if (verbose >= 3){
			printf("ineqv is: \n");
			print(ineqv); putchar('\n');
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
		ob = divide(firstPart, secondPart);
        
	} // end of if (ch>0)
    
	if (verbose >= 3){
		printf("ob is: \n");
		print(ob); putchar('\n');
	}
   
	j = M(ob, 0, 0);

	if (M(ind, 3, 0) > 0){
		ob = copyInto(ob, subtract(subset(ob, 0, neq+1, nc), subset(p, 0, 0, nineq-1)), 0, neq+1, nc);

	}

	if (nc > 0){

		ob = copyInto(ob, add(subtract(subset(ob, 0, 1, nc), matrixDotProduct(a, p)), b), 0, 1, nc);

        //exit(0);
		Matrix temp = subset(ob, 0, 1, nc);

		double vnormTerm = vnorm(temp) * vnorm(temp);
		Matrix yyTerm = transpose(yy);

		double dotProductTerm = dotProduct(getRow(yyTerm, 0), getRow(temp, 0));

		j = M(ob, 0, 0) - dotProductTerm + rho * vnormTerm;

        } // end if (nc > 0)
    
	minit = 0;
	Matrix obm = fill(1, 1, (double)0.0);
	Matrix yg = fill(npic, 1, (double)0.0);
	Matrix sx = fill(p.cols, 1, (double)0.0);
	Matrix sc = fill(2, 1, (double)0.0);
	Matrix cz;
	Matrix czSolution;
	Matrix u;
    
	int i;
	while (minit < maxit){
		minit = minit + 1;
        
		if (ch > 0){
            
			for (i=0; i<np; i++){
				int index = nineq+i;
				M(p, index, 0) = M(p, index, 0) + delta;
				tmpv = multiply(subset(p, 0, nineq, (npic-1)), subset(vscale, 0, (nc+1), (nc+np)));
				if (verbose >= 2){
					printf("9th call is \n");

				}
				funv = solFun(tmpv, verbose);
				if (verbose >= 3){
					printf("funv is: \n");
					printf("%2f", funv); putchar('\n');
				}
				eqv = solEqBFun(tmpv, verbose);
				ineqv = myineqFun(tmpv, verbose);
				solnp_nfn = solnp_nfn + 1;
				Matrix firstPart, secondPart, firstPartt;
                
				if (M(ineqv, 0, 0) != EMPTY){
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
				else {firstPart = fill(1, 1, funv);}
                
				secondPart = subset(vscale, 0, 0, nc);
				obm = divide(firstPart, secondPart);
				if (verbose >= 3){
					printf("obm is: \n");
					print(obm); putchar('\n');
					printf("j is: \n");
					printf("%.20f", j); putchar('\n');
				}
                
				if (M(ind, 3, 0) > 0.5){
					obm = copyInto(obm, subtract(subset(obm, 0, neq+1, nc), subset(p, 0, 0, nineq-1)), 0, neq+1, nc);
				}
                
				if (nc > 0){
                    
					Matrix first_part = subtract(subset(obm, 0, 1, nc),matrixDotProduct(a, p));
					obm = copyInto(obm, add(first_part, b), 0, 1, nc);
					Matrix temp = subset(obm, 0, 1, nc);
					double vnormTerm = vnorm(temp) * vnorm(temp);
					Matrix yyTerm = transpose(yy);
					double dotProductTerm = dotProduct(getRow(yyTerm, 0), getRow(temp, 0));
					double newOBMValue = M(obm, 0, 0) - dotProductTerm + rho * vnormTerm;
					obm = fill(1, 1, newOBMValue);    
				}
                
				M(g, index, 0) = (M(obm, 0, 0) - j)/delta;
				M(p, index, 0) = M(p, index, 0) - delta;
                
				if (verbose >= 3){
					printf("g is: \n");
					print(g); putchar('\n');
					printf("p is: \n");
					print(p); putchar('\n');
				}
			} // end for (i=0; i<np; i++){
            
			if (M(ind, 3, 0) > 0.5){
				g = copyInto(g, fill(nineq, 1, (double)0.0), 0, 0, (nineq-1));
			}
            
		} // end if (ch > 0){
        
		if (verbose >= 3){
			printf("yg is: \n");
			print(yg); putchar('\n');
		}

		if (minit > 1){
			yg = subtract(g, yg);
			sx = subtract(p, sx);
			M(sc, 0, 0) = dotProduct(getRow(matrixDotProduct(hessv, sx), 0), getRow(sx, 0));
			M(sc, 1, 0) = dotProduct(getRow(sx, 0), getRow(yg, 0));
			if ((M(sc, 0, 0) * M(sc, 1, 0)) > 0){
				//hessv  = hessv - ( sx %*% t(sx) ) / sc[ 1 ] + ( yg %*% t(yg) ) / sc[ 2 ]
				sx = matrixDotProduct(hessv, sx);
                
				Matrix sxMatrix = divideByScalar2D(transpose2D(sx), M(sc, 0, 0));
				Matrix ygMatrix = divideByScalar2D(transpose2D(yg), M(sc, 1, 0));
				
				hessv = subtract(hessv, sxMatrix);
				hessv = add(hessv, ygMatrix);
			}

		}
		if (verbose >= 3){ 
			printf("yg is: \n");
			print(yg); putchar('\n');
			printf("sx is: \n");
			print(sx); putchar('\n');
			printf("sc is: \n");
			print(sc); putchar('\n');
			printf("hessv is: \n");
			print(hessv); putchar('\n');
		}
        
		dx = fill(npic, 1, 0.01);
        
        
		if (M(ind, 10, 0) > 0.5){
			/** LOTS HERE BUT IT DOESN'T AFFECT US **/
            
			gap = fill(pb.cols, pb.rows, (double)0.0);
            
			gap = setColumn(gap, subtract(subset(p, 0, 0, mm-1), getColumn(pb, 0)), 0);
            
			gap = setColumn(gap, subtract(getColumn(pb, 1), subset(p, 0, 0, mm-1)), 1);
            
			gap = rowSort(gap);
            
			gap = add(getColumn(gap, 0), multiplyByScalar2D(fill(1, mm,(double)1.0),sqrt(eps)));
                        
			dx = copyInto(dx, divide(fill(mm, 1, (double)1.0), gap), 0, 0, mm-1);
            
			if (verbose >= 3){
				printf("dx is: \n");
				print(dx); putchar('\n');
			}

			if(M(ind, 9, 0) <= 0)
				{
					argum = multiplyByScalar2D(fill(1, npic-mm, (double)1.0) , min(findMin(subset(dx, 0, 0, mm-1)), 0.01));
					dx = copyInto(dx, argum, 0, mm, npic-1);
                
				}
            
		}
        
		if (verbose >= 3){
			printf("dx is: \n");
			print(dx); putchar('\n');
		}
		go = -1;
		lambdaValue = lambdaValue/10.0;
		int mahsa_count = 0;
		Matrix yMatrix;
        
		if (verbose >= 3){
			printf("lambdaValue is: \n");
			printf("%.20f", lambdaValue); putchar('\n');
			printf("hessv is: \n");
			print(hessv); putchar('\n');
		}

		while(go <= 0){
			Matrix dxDiagValues = multiply(dx, dx);
			if (verbose >= 3){
				printf("dxDiagValues is: \n");
				print(dxDiagValues); putchar('\n');
			}
			Matrix dxDiag = diag(dxDiagValues);
			if (verbose >= 3){
				printf("dxDiag is: \n");
				print(dxDiag); putchar('\n');
			}
			Matrix dxMultbyLambda = multiplyByScalar2D(dxDiag, lambdaValue);
			if (verbose >= 3){
				printf("dxMultbyLambda is: \n");
				print(dxMultbyLambda); putchar('\n');
			}
			Matrix addMatrices = add(hessv, dxMultbyLambda);
			if (verbose >= 3){
				printf("addMatrices is: \n");
				print(addMatrices); putchar('\n');
			}
			cz = cholesky(addMatrices);
			if (verbose >= 3){
				printf("cz is: \n");
				print(cz); putchar('\n');
			}
			Matrix identityMatrix = diag(fill(hessv.cols, 1, (double)1.0));
			cz = MatrixInvert(cz);
			//cz = luSolve(cz, identityMatrix);
			if (verbose >= 3){
				printf("cz is: \n");
				print(cz); putchar('\n');
			}
            
			Matrix getRowed = getRow(cz, 0);
			if (verbose >= 3){
				printf("g is: \n");
				print(g); putchar('\n');
			}
			Matrix getRowedtwo = getRow(g, 0);
			//double rr = dotProduct(getRowed, getRowedtwo);
			yg = matrixDotProduct(cz, g);
			if (verbose >= 3){
				printf("yg is: \n");
				print(yg); putchar('\n');
			}
			if (nc <= 0){
				u = matrixDotProduct(divideByScalar2D(transpose(cz), -1.0), yg);
				if (verbose >= 3){
					printf("u inside nc <=0 is: \n");
					print(u); putchar('\n');
				}
			}
			else{
				//y = qr.solve(t(cz) %*% t(a), yg)
                
				Matrix aTranspose = transpose(a);
				if (verbose >= 3){
					printf("aTranspose is: \n");
					print(aTranspose); putchar('\n');
				}
				Matrix firstMatrix = timess(cz, aTranspose);
				if (verbose >= 3){
					printf("firstMatrix is: \n");
					print(firstMatrix); putchar('\n');
				}
				Matrix secondMatrix = transpose(yg);
				if (verbose >= 3){
					printf("secondMatrix is: \n");
					print(secondMatrix); putchar('\n');
				}
				Matrix solution = qrSolve(firstMatrix, secondMatrix);
				if (verbose >= 3){
					printf("solution is: \n");
					print(solution); putchar('\n');
				}
				//Matrix solution = QRd(firstMatrix, secondMatrix);
				y = transpose(solution);
				if (verbose >= 3){
					printf("y is: \n");
					print(y); putchar('\n');
				}
				//y = subset(solution, 0, 0, nc-1);
				yMatrix = duplicateIt(solution);
				if (verbose >= 3){
					printf("yMatrix is: \n");
					print(yMatrix); putchar('\n');
				}
				//u = -cz %*% (yg - ( t(cz) %*% t(a) ) %*% y)
				Matrix minuscz = multiplyByScalar2D(transpose(cz), -1.0);
				Matrix toSubtract = timess(firstMatrix, yMatrix);
				Matrix partU = subtract(secondMatrix, toSubtract);
                
				u = timess(minuscz, partU);
				if (verbose >= 3){
					printf("u is: \n");
					print(u); putchar('\n');
				}
				u = transpose(u);
			}
            
			p0 = add(subset(u, 0, 0, npic-1), p);
			if (verbose >= 3){
				printf("p0 is: \n");
				print(p0); putchar('\n');
			}
            
			if (M(ind, 10, 0) <= 0.5){
				go = 1;
			} else{
				Matrix listPartOne = subtract(subset(p0, 0, 0, mm-1), getColumn(pb, 0));
                
				Matrix listPartTwo = subtract(getColumn(pb, 1), subset(p0, 0, 0, mm-1));
				Matrix list = copy(listPartOne, listPartTwo);
				go = findMin(list);
				lambdaValue = 3 * lambdaValue;
				if (verbose >= 3){
					printf("go is: \n");
					printf("%2f", go); putchar('\n');
				}
			}
		} // end while(go <= 0){
        
        
		M(alp, 0, 0) = 0;
		Matrix ob1 = duplicateIt(ob);
		Matrix ob2 = duplicateIt(ob1);
        
		M(sob, 0, 0) = j;
		M(sob, 1, 0) = j;
        
        
		if (verbose >= 3){
			printf("sob is is: \n");
			print(sob); putchar('\n');
		}
		ptt = copy(transpose(p), transpose(p));
        
		M(alp, 2, 0) = 1.0;
		if (verbose >= 3){
			printf("ptt is: \n");
			print(ptt); putchar('\n');
			printf("transpose(p0) is: \n");
			print(transpose(p0)); putchar('\n');
		}
        
		ptt = copy(ptt, transpose(p0));
       
		if (verbose >= 3){
			printf("ptt2 is: \n");
			print(ptt); putchar('\n');
		}
        
		bool condif1, condif2, condif3;
        
		Matrix pttCol = getColumn(ptt, 2);
		if (verbose >= 3){
			printf("pttCol is: \n");
			print(subset(pttCol, 0, nineq, (npic-1))); putchar('\n');
        
			printf("vscale is: \n");
			print(subset(vscale, 0, (nc+1), (nc+np))); putchar('\n');
		}
        
		tmpv = multiply(subset(pttCol, 0, nineq, (npic-1)), subset(vscale, 0, (nc+1), (nc+np)));
       
		if (verbose >= 3){
			printf("tmpv is: \n");
			print(tmpv); putchar('\n');
		}
		if (verbose >= 2){
			printf("10th call is \n");
		}
		funv = solFun(tmpv, verbose);
		if (verbose >= 3){
			printf("funv is: \n");
			printf("%.20f", funv); putchar('\n');
		}
		eqv = solEqBFun(tmpv, verbose);
        
		ineqv = myineqFun(tmpv, verbose);
        
		solnp_nfn = solnp_nfn + 1;
		Matrix firstPart, secondPart, firstPartt;
        
		if (M(ineqv, 0, 0) != EMPTY){
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
		else {firstPart = fill(1, 1, funv);}
        
		secondPart = subset(vscale, 0, 0, nc);
        
		Matrix ob3 = divide(firstPart, secondPart);
		if (verbose >= 3){
			printf("ob3 is: \n");
			print(ob3); putchar('\n');
		}
        
		M(sob, 2, 0) = M(ob3, 0, 0);
        
		if (M(ind, 3, 0) > 0.5){
			// ob3[ (neq + 2):(nc + 1) ] = ob3[ (neq + 2):(nc + 1) ] - ptt[ 1:nineq, 3 ]
			Matrix diff = fill((nineq+1), 1, (double)0.0);
			Matrix partOne = subset(ob3, 0, neq+1, nc);
			Matrix tempPttCol = getColumn(ptt, 2);
			Matrix partTwo = subset(tempPttCol, 0, 0, nineq);
			diff = subtract(partOne, partTwo);
			ob3 = copyInto(ob3, diff, 0, neq+1, nc);
		}
        
		if (nc > 0){
			//sob[ 3 ] = ob3[ 1 ] - t(yy) %*% ob3[ 2:(nc + 1) ] + rho * .vnorm(ob3[ 2:(nc + 1) ]) ^ 2
			Matrix firstp = subtract(subset(ob3, 0, 1, nc), matrixDotProduct(a, getColumn(ptt, 2)));
            
            ob3 = copyInto(ob3, add(firstp, b), 0, 1, nc);            
            Matrix temp = subset(ob3, 0, 1, nc);
            
            double vnormTerm = vnorm(temp) * vnorm(temp);
            
            Matrix yyTerm = transpose(yy);
            
            double dotProductTerm = dotProduct(getRow(yyTerm, 0), getRow(temp, 0));
            
            M(sob, 2, 0) = M(ob3, 0, 0) - dotProductTerm + (rho * vnormTerm);            
		}
                
		go = 1;
        
		while(go > tol){
            
			M(alp, 1, 0) = (M(alp, 0, 0) + M(alp, 2, 0)) / 2.0;
            
			Matrix colValues = add(multiplyByScalar2D(p, (1.0 - M(alp, 1, 0))), multiplyByScalar2D(p0, (M(alp, 1, 0))));
			ptt = setColumn(ptt, colValues, 1);
            
			Matrix pttColOne = getColumn(ptt, 1);
            
			tmpv = multiply(subset(pttColOne, 0, nineq, (npic-1)),
					subset(vscale, 0, (nc+1), (nc+np)));
			if (verbose >= 3){
				printf("tmpv is: \n");
				print(tmpv); putchar('\n');
			}
            
			if (verbose >= 2){ 
				printf("11th call is \n");
			}
         
			funv = solFun(tmpv, verbose);
			if (verbose >= 3){
				printf("funv is: \n");
				printf("%2f", funv); putchar('\n');
			}
            
			eqv = solEqBFun(tmpv, verbose);
            
			ineqv = myineqFun(tmpv, verbose);
			solnp_nfn = solnp_nfn + 1;
			Matrix firstPart, secondPart, firstPartt;
			if (M(ineqv, 0, 0) != EMPTY){
				if (M(eqv,0,0) != EMPTY)
					{
						firstPartt = copy(fill(1, 1, funv), eqv);
						firstPart = copy( firstPartt, ineqv);
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
            
			Matrix ob2 = divide(firstPart, secondPart);
			if (verbose >= 3){
				printf("ob2 is: \n");
				print(ob2); putchar('\n');
			}
            
			M(sob, 1, 0) = M(ob2, 0, 0);
			if (verbose >= 3){
				printf("sob is: \n");
				print(sob); putchar('\n');
			}
            //exit(0);
			if (M(ind, 3, 0) > 0.5){
				Matrix diff = fill(nineq+1, 1, (double)0.0);
				Matrix partOne = subset(ob2, 0, neq+1, nc);
				Matrix tempPttCol = getColumn(ptt, 1);
				Matrix partTwo = subset(tempPttCol, 0, 0, nineq-1);
				diff = subtract(partOne, partTwo);
				ob2 = copyInto(ob2, diff, 0, neq+1, nc);
                
			}
            
			if (nc > 0){
				ob2 = copyInto(ob2, add(subtract(subset(ob2, 0, 1, nc), matrixDotProduct(a, getColumn(ptt, 1))), b), 0, 1, nc);
				Matrix temp = subset(ob2, 0, 1, nc);
				double vnormTerm = vnorm(temp) * vnorm(temp);
				Matrix yyTerm = transpose(yy);
				double dotProductTerm = dotProduct(getRow(yyTerm, 0), getRow(temp, 0));
				M(sob, 1, 0) = M(ob2, 0, 0) - dotProductTerm + rho * vnormTerm;
			}
			if (verbose >= 3){
				printf("sob is: \n");
				print(sob); putchar('\n');
			}
			M(obm, 0, 0) = findMax(sob);
			if (verbose >= 3){
				printf("obm is: \n");
				print(obm); putchar('\n');
			}
			if (M(obm, 0, 0) < j){
				double obn = findMin(sob);
				go = tol * (M(obm, 0, 0) - obn) / (j - M(obm, 0, 0));
			}
            
            
			condif1 = (M(sob, 1, 0) >= M(sob, 0, 0));
			condif2 = (M(sob, 0, 0) <= M(sob, 2, 0)) && (M(sob, 1, 0) < M(sob, 0, 0));
			condif3 = (M(sob, 1, 0) <  M(sob, 0, 0)) && (M(sob, 0, 0) > M(sob, 2, 0));
            
			if (condif1){
				M(sob, 2, 0) = M(sob, 1, 0);
                
				ob3 = duplicateIt(ob2);
                
				M(alp, 2, 0) = M(alp, 1, 0);
                
				Matrix tempCol = getColumn(ptt, 1);
                
				ptt = setColumn(ptt, tempCol, 2);
                
				if (verbose >= 3){
					printf("sob is: \n");
					print(sob); putchar('\n');
					printf("ob3 is: \n");
					print(ob3); putchar('\n');
					printf("alp is: \n");
					print(alp); putchar('\n');
					printf("ptt is: \n");
					print(ptt); putchar('\n');
				}
                
			}
            
			if (condif2){
				M(sob, 2, 0) = M(sob, 1, 0);
				ob3 = duplicateIt(ob2);
				M(alp, 2, 0) = M(alp, 1, 0);
				Matrix tempCol = getColumn(ptt, 1);
                ptt = setColumn(ptt, tempCol, 2);

				if (verbose >= 3){
                    printf("sob is: \n");
					print(sob); putchar('\n');
					printf("ob3 is: \n");
					print(ob3); putchar('\n');
					printf("alp is: \n");
					print(alp); putchar('\n');
					printf("ptt is: \n");
					print(ptt); putchar('\n');
				}
			}
            
			if (condif3){
				M(sob, 0, 0) = M(sob, 1, 0);
				ob1 = duplicateIt(ob2);
				M(alp, 0, 0) = M(alp, 1, 0);
				Matrix tempCol = getColumn(ptt, 1);
				ptt = setColumn(ptt, tempCol, 0);
				if (verbose >= 3){
					printf("sob is: \n");
					print(sob); putchar('\n');
					printf("ob3 is: \n");
					print(ob3); putchar('\n');
					printf("alp is: \n");
					print(alp); putchar('\n');
					printf("ptt is: \n");
					print(ptt); putchar('\n');
				}
			}
			if (go >= tol){
				go = M(alp, 2, 0) - M(alp, 0, 0);
				if (verbose >= 3){
					printf("go is: \n");
					printf("%2f", go); putchar('\n');
				}
			}
		} // 	while(go > tol){
        //exit(0);
		if (verbose >= 3){
			printf("go is: \n");
			printf("%.16f", go); putchar('\n');
		}
		sx = duplicateIt(p);
		yg = duplicateIt(g);
		if (verbose >= 3){
			printf("sx is: \n");
			print(sx); putchar('\n');
			printf("yg is: \n");
			print(yg); putchar('\n');
		}
		ch = 1;
        
		double obn = findMin(sob);
		if (verbose >= 3){
			printf("obn is: \n");
			printf("%2f", obn); putchar('\n');
		}
		if (j <= obn){
			maxit = minit;
		}
		if (verbose >= 3){
			printf("j is: \n");
			printf("%2f", j); putchar('\n');
		}
		double reduce = (j - obn) / ((double)1.0 + (double)fabs(j));
		if (verbose >= 3){
			printf("reduce is: \n");
			printf("%.20f", reduce); putchar('\n');
		}
		if (reduce < tol){
			maxit = minit;
		}
        
		condif1 = (M(sob, 0, 0) <  M(sob, 1, 0));
		condif2 = (M(sob, 2, 0) <  M(sob, 1, 0)) && (M(sob, 0, 0) >=  M(sob, 1, 0));
		condif3 = (M(sob, 0, 0) >= M(sob, 1, 0)) && (M(sob, 2, 0) >= M(sob, 1, 0));
        
		if (condif1){
			j = M(sob, 0, 0);
			p = getColumn(ptt, 0);
			ob = duplicateIt(ob1);
			if (verbose >= 3){
				printf("j is: \n");
				printf("%2f", j); putchar('\n');
				printf("p is: \n");
				print(p); putchar('\n');
				printf("ob is: \n");
				print(ob); putchar('\n');
			}
		}
        
		if (condif2){
            
			j = M(sob, 2, 0);
			p = getColumn(ptt, 2);
			ob = duplicateIt(ob3);
			if (verbose >= 3){
				printf("j is: \n");
				printf("%2f", j); putchar('\n');
				printf("p is: \n");
				print(p); putchar('\n');
				printf("ob is: \n");
				print(ob); putchar('\n');
			}
            
		}
        
		if (condif3){
			j = M(sob, 1, 0);
			p = getColumn(ptt, 1);
			ob = duplicateIt(ob2);
			if (verbose >= 3){
				printf("j is: \n");
				printf("%2f", j); putchar('\n');
				printf("p is: \n");
				print(p); putchar('\n');
				printf("ob is: \n");
				print(ob); putchar('\n');
			}
		}
		if (verbose >= 3){
			printf("yg\n");
			print(yg); putchar('\n');
		}
	} // end while (minit < maxit){
    
	
	//p = p * vscale[ (neq + 2):(nc + np + 1) ]  # unscale the parameter vector
	Matrix vscalePart = subset(vscale, 0, (neq+1), (nc+np));
    
	p = multiply(p, vscalePart);
    
	if (nc > 0){
		//y = vscale[ 0 ] * y / vscale[ 2:(nc + 1) ] # unscale the lagrange multipliers
		y = multiplyByScalar2D(y, M(vscale,0,0));
		y = divide(y, subset(vscale, 0, 1, nc));
	}
    
	// hessv = vscale[ 1 ] * hessv / (vscale[ (neq + 2):(nc + np + 1) ] %*%
	//                                t(vscale[ (neq + 2):(nc + np + 1) ]) )
    
	Matrix transposePart = transpose2D(subset(vscale, 0, (neq+1), (nc+np)));
	hessv = divide(hessv, transposePart);
	hessv = multiplyByScalar2D(hessv, M(vscale,0,0));
    
	if (verbose >= 1 && reduce > tol) {
		printf("m3 solnp error message being reported.");
		putchar('\n');
	}
	
	resP = duplicateIt(p);
	resY = transpose(subset(y, 0, 0, (yyRows-1)));

	resHessv = duplicateIt(hessv);

	resLambda = lambdaValue;
    
	if (DEBUG && outerIter==4){
		printf("------------------------RETURNING FROM SUBNP------------------------"); putchar('\n');
		printf("p information: "); putchar('\n');
		print(resP); putchar('\n');
		printf("y information: "); putchar('\n');
		print(resY); putchar('\n');
		printf("hessv information: "); putchar('\n');
		print(resHessv); putchar('\n');
		printf("lambda information: "); putchar('\n');
		print(fill(1, 1, resLambda)); putchar('\n');
		printf("minit information: "); putchar('\n');
		print(fill(1, 1, minit)); putchar('\n');
		printf("------------------------END RETURN FROM SUBNP------------------------"); putchar('\n');
	}
    

	return g;
	
} // end subnp
















