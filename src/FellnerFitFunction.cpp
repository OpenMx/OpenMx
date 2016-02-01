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
	// Based on lme4CholmodDecomposition.h from lme4
	template<typename _MatrixType, int _UpLo = Eigen::Lower>
	class Cholmod : public Eigen::CholmodDecomposition<_MatrixType, _UpLo> {
	private:
		Eigen::MatrixXd ident;
		cholmod_dense* x_cd;

	protected:
		typedef Eigen::CholmodDecomposition<_MatrixType, _UpLo> Base;
		using Base::m_factorizationIsOk;
		typedef void (*cholmod_error_type)(int status, const char *file, int line, const char *message);

	        cholmod_factor* factor() const { return Base::m_cholmodFactor; }
		cholmod_common& cholmod() const {
			return const_cast<Cholmod<_MatrixType, _UpLo>*>(this)->Base::cholmod();
		}

		cholmod_error_type oldHandler;
		static void cholmod_error(int status, const char *file, int line, const char *message) {
			throw std::runtime_error(message);
		}

     // * If you are going to factorize hundreds or more matrices with the same
     // * nonzero pattern, you may wish to spend a great deal of time finding a
     // * good permutation.  In this case, try setting Common->nmethods to CHOLMOD_MAXMETHODS
     // * The time spent in cholmod_analysis will be very high, but you need to
     // * call it only once. TODO

	public:
		Cholmod() : x_cd(NULL) {
			oldHandler = cholmod().error_handler;
			cholmod().error_handler = cholmod_error;
			cholmod().supernodal = CHOLMOD_AUTO;
			cholmod().nmethods = CHOLMOD_MAXMETHODS;
			if (0) {
				cholmod().nmethods = 2;
				cholmod().method[0].ordering = CHOLMOD_NESDIS;
				cholmod().method[1].ordering = CHOLMOD_AMD;
			}
			m_factorizationIsOk = false; // base class should take care of it TODO
		};
		~Cholmod() {
			if (x_cd) cholmod_free_dense(&x_cd, &cholmod());
			cholmod().error_handler = oldHandler;
		};

		bool analyzedPattern() const { return m_factorizationIsOk; };

		void analyzePattern(const typename Base::MatrixType& matrix)
		{
			if (ident.rows() != matrix.rows()) {
				ident.setIdentity(matrix.rows(), matrix.rows());
			}
			Base::analyzePattern(matrix);
			cholmod_common &cm = cholmod();
			if (OMX_DEBUG) {
				mxLog("Cholmod: selected ordering %d lnz=%f fl=%f super=%d",
				      cm.method[cm.selected].ordering,
				      cm.method[cm.selected].lnz, cm.method[cm.selected].fl, factor()->is_super);
			}
		}

		double log_determinant() const {
			// Based on https://github.com/njsmith/scikits-sparse/blob/master/scikits/sparse/cholmod.pyx
			cholmod_factor *cf = factor();
			if (cf->xtype == CHOLMOD_PATTERN) Rf_error("Cannot extract diagonal from symbolic factor");
			double logDet = 0;
			double *x = (double*) cf->x;
			if (cf->is_super) {
				// This is a supernodal factorization, which is stored as a bunch
				// of dense, lower-triangular, column-major arrays packed into the
				// x vector. This is not documented in the CHOLMOD user-guide, or
				// anywhere else as far as I can tell; I got the details from
				// CVXOPT's C/cholmod.c.

				int *super = (int*) cf->super;
				int *pi = (int*) cf->pi;
				int *px = (int*) cf->px;
				for (size_t sx=0; sx < cf->nsuper; ++sx) {
					int ncols = super[sx + 1] - super[sx];
					int nrows = pi[sx + 1] - pi[sx];

					Eigen::Map<const Eigen::Array<double,1,Eigen::Dynamic>, 0, Eigen::InnerStride<> >
						s1(x + px[sx], ncols, Eigen::InnerStride<>(nrows+1));
					logDet += s1.real().log().sum();
				}
			} else {
				// This is a simplicial factorization, which is simply stored as a
				// sparse CSC matrix in x, p, i. We want the diagonal, which is
				// just the first entry in each column; p gives the offsets in x to
				// the beginning of each column.
				//
				// The ->p array actually has n+1 entries, but only the first n
				// entries actually point to real columns (the last entry is a
				// sentinel)
				int *p = (int*) cf->p;
				for (size_t ex=0; ex < cf->n; ++ex) {
					logDet += log( x[p[ex]] );
				}
			}
			if (cf->is_ll) {
				logDet *= 2.0;
			}
			return logDet;
		};

		double *getInverseData() const
		{
			return (double*) x_cd->x;
		}

		void refreshInverse()
		{
			eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
			cholmod_dense b_cd(viewAsCholmod(ident));
			if (x_cd) cholmod_free_dense(&x_cd, &cholmod());
			x_cd = cholmod_solve(CHOLMOD_A, factor(), &b_cd, &cholmod());
			if(!x_cd) throw std::runtime_error("cholmod_solve failed");
		};
	};

	struct state {
		int verbose;
		Cholmod< Eigen::SparseMatrix<double> > covDecomp;

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
		RelationalRAMExpectation::state *rram = ram->rram;
		omxData *data               = expectation->data;
		fc->profiledOut.assign(fc->numParam, false);

		ProtectedSEXP Rprofile(R_do_slot(oo->rObj, Rf_install("profileOut")));
		numProfiledOut = Rf_length(Rprofile);

		olsVarNum.reserve(numProfiledOut);
		olsDesign.resize(rram->dataVec.size(), numProfiledOut);

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
				olsDesign.col(px) = (rram->dataColumn.array() == vnum).cast<double>();
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
				for (size_t ax=0; ax < rram->layout.size(); ++ax) {
					RelationalRAMExpectation::addr &a1 = rram->layout[ax];
					if (a1.model != expectation) continue;
					data->handleDefinitionVarList(ram->M->currentState, a1.row);
					double weight = omxVectorElement(ram->M, vnum);
					olsDesign.col(px).segment(a1.obsStart, a1.numObs()) =
						weight * (rram->dataColumn.segment(a1.obsStart, a1.numObs()) == rnum).cast<double>();
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

		double lp = NA_REAL;
		try {
			omxExpectationCompute(fc, expectation, "covariance", "flat");
			RelationalRAMExpectation::state *rram   = ram->rram;

			if (!covDecomp.analyzedPattern()) {
				rram->fullCov.makeCompressed();
				covDecomp.analyzePattern(rram->fullCov);
			}

			covDecomp.factorize(rram->fullCov);
			covDecomp.refreshInverse();
			Eigen::Map< Eigen::MatrixXd > iV(covDecomp.getInverseData(),
							 rram->fullCov.rows(), rram->fullCov.rows());

			double remlAdj = 0.0;
			if (numProfiledOut) {
				Eigen::MatrixXd constCov = olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() * olsDesign;
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
					 olsDesign.transpose() * iV.selfadjointView<Eigen::Lower>() * rram->dataVec);

				for (int px=0; px < numProfiledOut; ++px) {
					fc->est[ olsVarNum[px] ] = param[px];
					fc->varGroup->vars[ olsVarNum[px] ]->copyToState(ram->M->currentState, param[px]);
				}
			}

			omxExpectationCompute(fc, expectation, "mean", "flat");
			//mxPrintMat("dataVec", rram->dataVec);
			//mxPrintMat("fullMeans", rram->fullMeans);
			Eigen::VectorXd resid = rram->dataVec - rram->expectedMean;
			rram->applyRotationPlan(resid);
			//mxPrintMat("resid", resid);

			lp = covDecomp.log_determinant();
			double iqf = resid.transpose() * iV.selfadjointView<Eigen::Lower>() * resid;
			double cterm = M_LN_2PI * (rram->dataVec.size() - numProfiledOut);
			if (verbose >= 2) mxLog("log det %f iqf %f cterm %f remlAdj %f", lp, iqf, cterm, remlAdj);
			lp += iqf + cterm + remlAdj;
		} catch (const std::exception& e) {
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		}
		oo->matrix->data[0] = lp;
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
