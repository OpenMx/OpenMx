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

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include <sys/stat.h>
#include <errno.h>

#include "omxDefines.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxImportFrontendState.h"
#include "matrix.h"

Matrix result_csolnpEqBStartFun;
Matrix result_csolnpEqB;
Matrix result_ineqLB;
Matrix result_ineqUB;
Matrix result_ineqVal;
//double EMPTY = -999999.0;

void setupIneqLess(struct Matrix *bl_ineqless, struct Matrix *bu_ineqless, int size)
    {
        int index = 0;
        for(int offset = 0; offset < size; offset++) {
            bl_ineqless->t[index] = NEG_INF;
            bu_ineqless->t[index] = -0.0;
            index++;
        }
    }

void setupIneqGreater(struct Matrix *lb, struct Matrix *ub, int size)
{
        int index = 0;
        for(int offset = 0; offset < size; offset++) {
            lb->t[index] = 0.0;
            ub->t[index] = INF;
            printf("%2f", lb->t[index]); putchar('\n');
            printf("%2f", ub->t[index]); putchar('\n');
            index++;
        }
}
    
void setupEqB(struct Matrix *eqPointer, int size)
    {
        int index = 0;
        for(int offset = 0; offset < size; offset++) {
            eqPointer->t[index] = 0.0;
            index++;
        }
    }
    
    
void omxProcessConstraintsCsolnp(struct Matrix *lb_ineq, struct Matrix *ub_ineq, struct Matrix *eqb)  {
    if (globalState->numConstraints == 0){
        return;
    }
    
    int i;
    int size;
	double EMPTY = -999999.0;
    Matrix lb_ineqless = fill(1, 1, EMPTY);
    Matrix ub_ineqless = fill(1, 1, EMPTY);
    Matrix lb_ineqmore = fill(1, 1, EMPTY);
    Matrix ub_ineqmore = fill(1, 1, EMPTY);
    Matrix eqbound = fill(1, 1, EMPTY);
    Matrix ineq_lb = fill(1, 1, EMPTY);
    Matrix ineq_ub = fill(1, 1, EMPTY);
    Matrix myeqb = fill(1, 1, EMPTY);
    
	for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {

		size = globalState->conList[constraintIndex].size;
        
            if (globalState->conList[constraintIndex].opCode == 0)
            {
                Matrix lb_ineqless = fill(size, 1, EMPTY);
                Matrix ub_ineqless = fill(size, 1, EMPTY);
                setupIneqLess(&lb_ineqless, &ub_ineqless, size);
            }
            
            if (globalState->conList[constraintIndex].opCode == 1)
            {
                eqbound = fill(size, 1, EMPTY);
                setupEqB(&eqbound, size);
            }
            
            if (globalState->conList[constraintIndex].opCode == 2)
            {
                lb_ineqmore = fill(size, 1, EMPTY);
                ub_ineqmore = fill(size, 1, EMPTY);
                setupIneqGreater(&lb_ineqmore, &ub_ineqmore, size);
            }
            
        printf("lb_ineqmore is: \n");
        print(lb_ineqmore); putchar('\n');
        printf("ub_ineqmore is: \n");
        print(ub_ineqmore); putchar('\n');
        
        if (M(lb_ineqless, 0, 0) != EMPTY)
        {
            printf("inside lb_ineqless \n");
            ineq_lb = copy(ineq_lb, lb_ineqless);
            ineq_ub = copy(ineq_ub, ub_ineqless);
        }
        else if (M(lb_ineqmore, 0, 0) != EMPTY)
        {
            printf("inside lb_ineqmore \n");
            ineq_lb = copy(ineq_lb, lb_ineqmore);
            ineq_ub = copy(ineq_ub, ub_ineqmore);
        }
        else if (M(eqbound, 0, 0) != EMPTY)
        {
            printf("inside lb_eqb \n");
            myeqb = copy(myeqb, eqbound);
        }
    }

    if (ineq_lb.cols > 1)
    {
        ineq_lb = subset(ineq_lb, 0, 1, ineq_lb.cols - 1);
        ineq_ub = subset(ineq_ub, 0, 1, ineq_ub.cols - 1);
    }
    
    if (myeqb.cols > 1)
    {
        myeqb = subset(myeqb, 0, 1, myeqb.cols - 1);
    }
    
    for (i = 0; i < ineq_lb.cols; i++)
    {
        lb_ineq->t[i] = M(ineq_lb, i, 0);
        ub_ineq->t[i] = M(ineq_ub, i, 0);
        printf("pointer check \n");
        printf("%2f", lb_ineq->t[i]); putchar('\n');
        printf("%2f", ub_ineq->t[i]); putchar('\n');
    }
    printf("%d", myeqb.cols); putchar('\n');
    
    for (i = 0; i < myeqb.cols; i++)
    {
        eqb->t[i] = M(myeqb, i, 0);
        printf("%2f", eqb->t[i]); putchar('\n');
    }
    
    lb_ineq->cols = ineq_lb.cols;
    ub_ineq->cols = ineq_ub.cols;
    eqb->cols = myeqb.cols;
    //printf("eqb is: \n");
    //print(eqb); putchar('\n');

}
