/*
 * Copyright 2015-2016 Joshua Nathaniel Pritikin and contributors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
#include "Compute.h"

namespace FellnerFitFunction {
	struct state : omxFitFunction {
		int verbose;
		int numProfiledOut;
		std::vector<int> olsVarNum;     // index into fc->est
		Eigen::MatrixXd olsDesign;      // a.k.a "X"

		virtual void init();
		template <typename T1>
		int computeCov(RelationalRAMExpectation::independentGroup &ig, FitContext *fc, T1 &covDecomp);
		virtual void compute(int want, FitContext *fc);
		void setupProfiledParam(FitContext *fc);
	};

	void state::setupProfiledParam(FitContext *fc)
	{
		auto *oo = this;
		omxExpectation *expectation             = oo->expectation;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

		if (numProfiledOut) ram->forceSingleGroup = true;
		omxExpectationCompute(fc, expectation, "nothing", "flat");
		
		if (numProfiledOut == 0) return;

		RelationalRAMExpectation::state &rram = ram->rram->getParent();
		if (rram.group.size() > 1) {
			mxThrow("Cannot profile out parameters when problem is split into independent groups");
		}

		RelationalRAMExpectation::independentGroup &ig = *rram.group[0];
		fc->profiledOut.assign(fc->numParam, false);

		olsVarNum.reserve(numProfiledOut);
		olsDesign.resize(ig.getParent().dataVec.size(), numProfiledOut);

		ProtectedSEXP Rprofile(R_do_slot(oo->rObj, Rf_install("profileOut")));
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
					RelationalRAMExpectation::addr &a1 = rram.layout[ ig.gMap[ax] ];
					if (a1.getExpNum() != expectation->expNum) continue;
					expectation->loadDefVars(a1.row);
					double weight = omxVectorElement(ram->M, vnum);
					olsDesign.col(px).segment(pl.obsStart, a1.numObs()) =
						weight * (ig.dataColumn.segment(pl.obsStart, a1.numObs()) == rnum).cast<double>();
				}
			}
			if (!found) mxThrow("oops");

			fc->profiledOut[vx] = true;
		}
	}

	template <typename T1>
	int state::computeCov(RelationalRAMExpectation::independentGroup &ig, FitContext *fc, T1 &covDecomp)
	{
		if (0 == ig.getParent().dataVec.size()) return 0;

		ig.computeCov(fc);

		/*
		if (!ig.analyzedCov) {
			ig.fullCov.makeCompressed();
			ig.covDecomp.analyzePattern(ig.fullCov);
			ig.analyzedCov = true;
		}
		ig.covDecomp.factorize(ig.fullCov);
		*/

		Eigen::MatrixXd &denseCov = ig.fullCov;
		covDecomp.compute(denseCov);

		if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) return 1;

		covDecomp.refreshInverse();
		return 0;
	}

	void state::compute(int want, FitContext *fc)
	{
		auto *oo = this;
		omxExpectation *expectation             = oo->expectation;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

		if (want & (FF_COMPUTE_INITIAL_FIT | FF_COMPUTE_PREOPTIMIZE)) {
			if (fc->isClone()) return;
			
			setupProfiledParam(fc);

			RelationalRAMExpectation::state *rram   = ram->rram;
			if (verbose >= 1) {
				mxLog("%s: %d groups", oo->name(),
				      int(rram->group.size()));
			}
			return;
		}

		if (!(want & (FF_COMPUTE_FIT))) mxThrow("Not implemented");

		double lpOut = NA_REAL;
		try {
			if (!ram->rram) {
				// it is possible to skip FF_COMPUTE_PREOPTIMIZE (e.g. omxRecompute)
				omxExpectationCompute(fc, expectation, "nothing", "flat");
			}

			RelationalRAMExpectation::state &rram   = *ram->rram;
			state *parent = this; // better to cache it TODO
			if (fc->isClone()) {
				omxMatrix *pfitMat = fc->getParentState()->getMatrixFromIndex(oo->matrix);
				parent = (state*) pfitMat->fitFunction;
			}

			//mxLog("%s: compute fit", oo->name());

			SimpCholesky< Eigen::MatrixXd > covDecomp;
			bool haveMean = false;
			double lp = 0.0;
			for (size_t gx=0; gx < rram.group.size(); ++gx) {
				RelationalRAMExpectation::independentGroup &ig = *rram.group[gx];
				double lp1 = 0.0;

				if (computeCov(*rram.group[gx], fc, covDecomp)) {
					throw std::runtime_error("Cholesky decomposition failed");
				}

				if (rram.group.size() == 1 && parent->numProfiledOut) {
					double remlAdj = 0.0;
					const Eigen::MatrixXd &iV = covDecomp.getInverse();
					Eigen::MatrixXd constCov =
						parent->olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() * parent->olsDesign;
					Eigen::LLT< Eigen::MatrixXd > cholConstCov;
					cholConstCov.compute(constCov);
					if(cholConstCov.info() != Eigen::Success){
						// ought to report error detail TODO
						throw std::exception();
					}
					remlAdj = 2*Eigen::MatrixXd(cholConstCov.matrixL()).diagonal().array().log().sum();

					Eigen::MatrixXd ident =
						Eigen::MatrixXd::Identity(parent->numProfiledOut, parent->numProfiledOut);
					Eigen::MatrixXd cholConstPrec = cholConstCov.solve(ident).triangularView<Eigen::Lower>();
					Eigen::VectorXd param =
						(cholConstPrec.selfadjointView<Eigen::Lower>() *
						 parent->olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() *
						 ig.getParent().dataVec);

					for (int px=0; px < parent->numProfiledOut; ++px) {
						fc->est[ parent->olsVarNum[px] ] = param[px];
						fc->varGroup->vars[ parent->olsVarNum[px] ]->copyToState(ram->M->currentState, param[px]);
					}
					lp1 += remlAdj - M_LN_2PI * parent->numProfiledOut;
				}

				if (!haveMean) {
					// delay until after remlAdj done
					omxExpectationCompute(fc, expectation, "mean", "flat");
					haveMean = true;
				}

				if (0 == ig.getParent().dataVec.size()) continue;

				//mxPrintMat("dataVec", ig.dataVec);
				//mxPrintMat("fullMeans", ig.fullMeans);
				//ig.applyRotationPlan(ig.expectedVec);
				//mxPrintMat("expectedVec", ig.expectedVec);

				const Eigen::MatrixXd &iV = covDecomp.getInverse();
				double logDet = 2.0 * covDecomp.log_determinant();
				//mxPrintMat("iV", iV);
				int clumps = ig.numLooseClumps();
				if (clumps) {
					int residLen = clumps * ig.clumpObs;
					Eigen::VectorXd resid =
						ig.getParent().dataVec.segment(0,residLen) - ig.expectedVec.segment(0,residLen);
					//mxPrintMat("resid", resid);

				// Eigen::Map< Eigen::MatrixXd > iV(covDecomp.getInverseData(),
				// 				 ig.fullCov.rows(), ig.fullCov.rows());
					double iqf = 0.0;
					for (int cx=0; cx < clumps; ++cx) {
						iqf += (resid.segment(cx*ig.clumpObs, ig.clumpObs).transpose() *
							iV.selfadjointView<Eigen::Lower>() *
							resid.segment(cx*ig.clumpObs, ig.clumpObs));
					}
					double cterm = M_LN_2PI * residLen;
					if (verbose + !std::isfinite(iqf) >= 2) {
						mxLog("group[%d] log det %f iqf %f cterm %f",
						      int(1+gx), clumps * logDet, iqf, cterm);
					}
					lp1 += clumps * logDet + iqf + cterm;
				}
				if (ig.getParent().sufficientSets.size()) {
					double cterm = M_LN_2PI * ig.clumpObs;
					for (int sx=0; sx < (int)ig.getParent().sufficientSets.size(); ++sx) {
						RelationalRAMExpectation::sufficientSet &ss = ig.getParent().sufficientSets[sx];
						Eigen::VectorXd resid =
							ss.dataMean - ig.expectedVec.segment(ss.start * ig.clumpObs, ig.clumpObs);
						//mxPrintMat("resid", resid);
						double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
						double tr1 = trace_prod(iV, ss.dataCov);
						if (verbose + !std::isfinite(iqf) >= 2) {
							mxLog("group[%d] ss[%d] iqf %f tr1 %f logDet %f cterm %f",
							      int(1+gx), (1+sx), iqf, tr1, logDet, cterm);
						}
						lp1 += ss.length * (iqf + logDet + cterm) + (ss.length-1) * tr1;
					}
				}
				ig.fit = lp1;
				lp += lp1;
				//mxLog("%s: group[%d] lp1 %.6g lp %.6g", oo->name(), int(1+gx), lp1, lp);
			}
			lpOut = lp;
			//mxLog("%s: total lp %.7g", oo->name(), lpOut);
		} catch (const std::exception& e) {
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		} catch (...) {
			if (fc) fc->recordIterationError("%s: unknown error", oo->name());
		}
		oo->matrix->data[0] = lpOut;
	}

	void state::init()
	{
		auto *oo = this;
		auto *st = this;

		omxExpectation* expectation = oo->expectation;
		if(expectation == NULL) {
			omxRaiseErrorf("%s cannot fit without a model expectation", oo->fitType);
			return;
		}
		if (!strEQ(expectation->name, "MxExpectationRAM")) {
			mxThrow("%s: only MxExpectationRAM is implemented", oo->matrix->name());
		}

		oo->canDuplicate = true;

		ProtectedSEXP Rprofile(R_do_slot(oo->rObj, Rf_install("profileOut")));
		st->numProfiledOut = Rf_length(Rprofile);

		{
			SEXP tmp;
			ScopedProtect p1(tmp, R_do_slot(oo->rObj, Rf_install("verbose")));
			st->verbose = Rf_asInteger(tmp) + OMX_DEBUG;
		}
	}
};

omxFitFunction *InitFellnerFitFunction()
{
	return new FellnerFitFunction::state;
}
