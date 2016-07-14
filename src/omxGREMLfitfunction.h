/*
 *  Copyright 2007-2016 The OpenMx Project
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

struct omxGREMLFitState { 
	//TODO(?): Some of these members might be redundant with what's stored in the FitContext, 
	//and could therefore be cut
	omxMatrix *y, *X, *cov, *invcov, *means, *origVdim_om;
	std::vector< omxMatrix* > dV;
	std::vector< const char* > dVnames;
	std::vector<int> indyAlg; //will keep track of which algebras don't get marked dirty after dropping cases
	std::vector<int> origdVdim;
	void dVupdate(FitContext *fc);
	void dVupdate_final();
	int dVlength, usingGREMLExpectation, parallelDerivScheme;
	double nll, REMLcorrection;
	Eigen::VectorXd gradient;
	Eigen::MatrixXd avgInfo; //the Average Information matrix
	FreeVarGroup *varGroup;
	std::vector<int> gradMap;
	void buildParamMap(FreeVarGroup *newVarGroup);
	std::vector< Eigen::VectorXi > rowbins, AIMelembins;
	void planParallelDerivs(int nThreadz, int wantHess, int Vrows);
	omxMatrix *aug, *augGrad, *augHess;
	std::vector<int> dAugMap;
	double pullAugVal(int thing, int row, int col);
	void recomputeAug(int thing, FitContext *fc);
}; 

void omxDestroyGREMLFitFunction(omxFitFunction *oo);

void omxInitGREMLFitFunction(omxFitFunction *oo);

void omxCallGREMLFitFunction(omxFitFunction *oo, int want, FitContext *fc);

static void omxPopulateGREMLAttributes(omxFitFunction *oo, SEXP algebra);

