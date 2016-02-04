/*
  Copyright 2015 Joshua Nathaniel Pritikin and contributors

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

// Named in honor of Fellner (1987) "Sparse matrices, and the
// estimation of variance components by likelihood methods"
// Fellner was probably the first to apply sparse matrix algorithms
// to this kind of problem.

#include "glue.h"
#include <iterator>
#include <exception>
#include <stdexcept>
#include <Rconfig.h>
#include <Rmath.h>
#include "omxFitFunction.h"
#include "RAMInternal.h"
#include <Eigen/Cholesky>

namespace FellnerFitFunction {
	struct state {
		int verbose;

		int numProfiledOut;
		std::vector<int> olsVarNum;     // index into fc->est
		Eigen::MatrixXd olsDesign;      // a.k.a "X"

		void compute(omxFitFunction *oo, int want, FitContext *fc);
		void setupProfiledParam(omxFitFunction *oo, FitContext *fc);
	};

	static void compute(omxFitFunction *oo, int want, FitContext *fc)
	{
		state *st = (state *) oo->argStruct;
		st->compute(oo, want, fc);
	}

	void state::setupProfiledParam(omxFitFunction *oo, FitContext *fc)
	{
		omxExpectation *expectation             = oo->expectation;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
		omxExpectationCompute(fc, expectation, "nothing", "flat");

		ProtectedSEXP Rprofile(R_do_slot(oo->rObj, Rf_install("profileOut")));
		numProfiledOut = Rf_length(Rprofile);
		if (numProfiledOut == 0) return;
		
		RelationalRAMExpectation::state *rram = ram->rram;
		if (rram->group.size() > 1) {
			Rf_error("Cannot profile out parameters when problem is split into independent groups");
		}

		RelationalRAMExpectation::independentGroup &ig = *rram->group[0];
		omxData *data               = expectation->data;
		fc->profiledOut.assign(fc->numParam, false);

		olsVarNum.reserve(numProfiledOut);
		olsDesign.resize(ig.dataVec.size(), numProfiledOut);

		for (int px=0; px < numProfiledOut; ++px) {
			const char *pname = CHAR(STRING_ELT(Rprofile, px));
			int vx = fc->varGroup->lookupVar(pname);
			if (vx < 0) {
				mxLog("Parameter [%s] not found", pname);
				continue;
			}

			omxFreeVar &fv = *fc->varGroup->vars[vx];
			olsVarNum.push_back(vx);
			bool found = false;
			bool moreThanOne;
			const omxFreeVarLocation *loc =
				fv.getOnlyOneLocation(ram->M, moreThanOne);
			if (loc) {
				if (moreThanOne) {
					mxLog("Parameter [%s] appears in more than one spot in %s",
					      pname, ram->M->name());
					continue;
				}
				found = true;
				int vnum = loc->row + loc->col;
				// Should ensure the loading is fixed and not a defvar TODO
				// Should ensure zero variance & no cross-level links TODO
				olsDesign.col(px) = (ig.dataColumn.array() == vnum).cast<double>();
			}
			loc = fv.getOnlyOneLocation(ram->A, moreThanOne);
			if (loc) {
				if (moreThanOne) {
					mxLog("Parameter [%s] appears in more than one spot in %s",
					      pname, ram->A->name());
					continue;
				}
				found = true;
				int vnum = loc->col;
				EigenMatrixAdaptor eA(ram->A);
				int rnum;
				eA.col(vnum).array().abs().maxCoeff(&rnum);
				// ensure only 1 nonzero in column TODO
				for (size_t ax=0; ax < ig.placements.size(); ++ax) {
					RelationalRAMExpectation::placement &pl = ig.placements[ax];
					RelationalRAMExpectation::addr &a1 = rram->layout[ pl.aIndex ];
					if (a1.model != expectation) continue;
					data->handleDefinitionVarList(ram->M->currentState, a1.row);
					double weight = omxVectorElement(ram->M, vnum);
					olsDesign.col(px).segment(pl.obsStart, a1.numObs()) =
						weight * (ig.dataColumn.segment(pl.obsStart, a1.numObs()) == rnum).cast<double>();
				}
			}
			if (!found) Rf_error("oops");

			fc->profiledOut[vx] = true;
		}
	}

	void state::compute(omxFitFunction *oo, int want, FitContext *fc)
	{
		omxExpectation *expectation             = oo->expectation;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		if (want & (FF_COMPUTE_PREOPTIMIZE)) {
			setupProfiledParam(oo, fc);
			return;
		}

		if (!(want & (FF_COMPUTE_FIT | FF_COMPUTE_INITIAL_FIT))) Rf_error("Not implemented");

		double lpOut = NA_REAL;
		try {
			double lp = 0.0;
			omxExpectationCompute(fc, expectation, "covariance", "flat");

			double remlAdj = 0.0;
			RelationalRAMExpectation::state *rram   = ram->rram;
			for (size_t gx=0; gx < rram->group.size(); ++gx) {
				RelationalRAMExpectation::independentGroup &ig = *rram->group[gx];
				if (0 == ig.dataVec.size()) continue;

				if (!ig.covDecomp.analyzedPattern()) {
					ig.fullCov.makeCompressed();
					ig.covDecomp.analyzePattern(ig.fullCov);
				}

				ig.covDecomp.factorize(ig.fullCov);
				ig.covDecomp.refreshInverse();
				Eigen::Map< Eigen::MatrixXd > iV(ig.covDecomp.getInverseData(),
								 ig.fullCov.rows(), ig.fullCov.rows());

				if (numProfiledOut) {
					Eigen::MatrixXd constCov =
						olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() * olsDesign;
					Eigen::LLT< Eigen::MatrixXd > cholConstCov;
					cholConstCov.compute(constCov);
					if(cholConstCov.info() != Eigen::Success){
						// ought to report error detail TODO
						throw std::exception();
					}
					remlAdj = 2*Eigen::MatrixXd(cholConstCov.matrixL()).diagonal().array().log().sum();

					Eigen::MatrixXd ident = Eigen::MatrixXd::Identity(numProfiledOut, numProfiledOut);
					Eigen::MatrixXd cholConstPrec = cholConstCov.solve(ident).triangularView<Eigen::Lower>();
					Eigen::VectorXd param =
						(cholConstPrec.selfadjointView<Eigen::Lower>() *
						 olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() * ig.dataVec);

					for (int px=0; px < numProfiledOut; ++px) {
						fc->est[ olsVarNum[px] ] = param[px];
						fc->varGroup->vars[ olsVarNum[px] ]->copyToState(ram->M->currentState, param[px]);
					}
				}
			}

			omxExpectationCompute(fc, expectation, "mean", "flat");

			for (size_t gx=0; gx < rram->group.size(); ++gx) {
				RelationalRAMExpectation::independentGroup &ig = *rram->group[gx];
				if (0 == ig.dataVec.size()) continue;

				//mxPrintMat("dataVec", ig.dataVec);
				//mxPrintMat("fullMeans", ig.fullMeans);
				//ig.applyRotationPlan(ig.expectedMean);
				//mxPrintMat("expectedMean", ig.expectedMean);

				Eigen::VectorXd resid = ig.dataVec - ig.expectedMean;
				//mxPrintMat("resid", resid);

				double logDet = ig.covDecomp.log_determinant();
				Eigen::Map< Eigen::MatrixXd > iV(ig.covDecomp.getInverseData(),
								 ig.fullCov.rows(), ig.fullCov.rows());
				double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
				double cterm = M_LN_2PI * (ig.dataVec.size() - numProfiledOut);
				if (verbose >= 2) mxLog("log det %f iqf %f cterm %f remlAdj %f", logDet, iqf, cterm, remlAdj);
				lp += logDet + iqf + cterm + remlAdj;
			}
			lpOut = lp;
		} catch (const std::exception& e) {
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		}
		oo->matrix->data[0] = lpOut;
	}

	static void popAttr(omxFitFunction *oo, SEXP algebra)
	{
		// use Eigen_cholmod_wrap to return a sparse matrix? TODO
		// always return it?

		/*
		state *st                               = (state *) oo->argStruct;
		SEXP expCovExt, expMeanExt;
		if (st->fullCov.rows() > 0) {
			Rf_protect(expCovExt = Rf_allocMatrix(REALSXP, expCovInt->rows, expCovInt->cols));
			memcpy(REAL(expCovExt), expCovInt->data, sizeof(double) * expCovInt->rows * expCovInt->cols);
			Rf_setAttrib(algebra, Rf_install("expCov"), expCovExt);
		}

		if (expMeanInt && expMeanInt->rows > 0) {
			Rf_protect(expMeanExt = Rf_allocMatrix(REALSXP, expMeanInt->rows, expMeanInt->cols));
			memcpy(REAL(expMeanExt), expMeanInt->data, sizeof(double) * expMeanInt->rows * expMeanInt->cols);
			Rf_setAttrib(algebra, Rf_install("expMean"), expMeanExt);
			}   */
	}

	static void destroy(omxFitFunction *oo)
	{
		state *st = (state*) oo->argStruct;
		delete st;
	}

	static void init(omxFitFunction *oo)
	{
		omxExpectation* expectation = oo->expectation;
		if(expectation == NULL) {
			omxRaiseErrorf("%s cannot fit without a model expectation", oo->fitType);
			return;
		}
		if (!strEQ(expectation->expType, "MxExpectationRAM")) {
			Rf_error("%s: only MxExpectationRAM is implemented", oo->matrix->name());
		}

		// prohibit ordinal for now TODO
		if (expectation->numOrdinal != 0) {
			Rf_error("%s cannot handle ordinal data yet", oo->fitType);
		}

		oo->computeFun = FellnerFitFunction::compute;
		oo->destructFun = FellnerFitFunction::destroy;
		oo->populateAttrFun = FellnerFitFunction::popAttr;
		FellnerFitFunction::state *st = new FellnerFitFunction::state;
		oo->argStruct = st;

		st->numProfiledOut = 0;
		{
			SEXP tmp;
			ScopedProtect p1(tmp, R_do_slot(oo->rObj, Rf_install("verbose")));
			st->verbose = Rf_asInteger(tmp) + OMX_DEBUG;
		}
	}
};

void InitFellnerFitFunction(omxFitFunction *oo)
{
	FellnerFitFunction::init(oo);
}
