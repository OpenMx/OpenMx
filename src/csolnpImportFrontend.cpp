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

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include <sys/stat.h>
#include <errno.h>

#include "omxDefines.h"
#include "omxState.h"
#include "omxNPSOLSpecific.h"
#include "omxImportFrontendState.h"
#include "matrix.h"

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


void omxProcessConstraintsCsolnp(FitContext *fc, struct Matrix lb_ineq, struct Matrix ub_ineq, struct Matrix eqb)  {
    omxState *globalState = fc->state;
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
    Matrix ineq_lb = fill(0, 1, EMPTY);
    Matrix ineq_ub = fill(0, 1, EMPTY);
    Matrix myeqb = fill(0, 1, EMPTY);
    
    for(int constraintIndex = 0; constraintIndex < globalState->numConstraints; constraintIndex++) {
        
        size = globalState->conList[constraintIndex].size;
	if (size == 0) continue;
        
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
        
        if (M(lb_ineqless, 0, 0) != EMPTY)
        {
            ineq_lb = copy(ineq_lb, lb_ineqless);
            ineq_ub = copy(ineq_ub, ub_ineqless);
        }
        else if (M(lb_ineqmore, 0, 0) != EMPTY)
        {
            ineq_lb = copy(ineq_lb, lb_ineqmore);
            ineq_ub = copy(ineq_ub, ub_ineqmore);
        }
        else if (M(eqbound, 0, 0) != EMPTY)
        {
            myeqb = copy(myeqb, eqbound);
        }
    }
    
    if (lb_ineq.cols < ineq_lb.cols ||
        ub_ineq.cols < ineq_ub.cols ||
        eqb.cols < myeqb.cols) {
        Rf_error("Failed to unpack constraint bounds");
    }
    
    for (i = 0; i < ineq_lb.cols; i++)
    {
        lb_ineq.t[i] = M(ineq_lb, i, 0);
        ub_ineq.t[i] = M(ineq_ub, i, 0);
    }
    
    for (i = 0; i < myeqb.cols; i++)
    {
        eqb.t[i] = M(myeqb, i, 0);
    }
    
}
