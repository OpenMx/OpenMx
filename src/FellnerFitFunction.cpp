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
#include <Rconfig.h>
#include <Rmath.h>
#include <RcppEigenCholmod.h>
#include <RcppEigenStubs.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>
#include "omxFitFunction.h"

// New Expectation API : per case to full distribution adapter TODO

namespace FellnerFitFunction {
	// Based on lme4CholmodDecomposition.h from lme4
	template<typename _MatrixType, int _UpLo = Eigen::Lower>
	class Cholmod : public Eigen::CholmodDecomposition<_MatrixType, _UpLo> {
	private:
		Eigen::MatrixXd ident;

	protected:
		typedef Eigen::CholmodDecomposition<_MatrixType, _UpLo> Base;
		using Base::m_factorizationIsOk;

	        cholmod_factor* factor() const { return Base::m_cholmodFactor; }
		cholmod_common& cholmod() const {
			return const_cast<Cholmod<_MatrixType, _UpLo>*>(this)->Base::cholmod();
		}

     // * If you are going to factorize hundreds or more matrices with the same
     // * nonzero pattern, you may wish to spend a great deal of time finding a
     // * good permutation.  In this case, try setting Common->nmethods to 9.
     // * The time spent in cholmod_analysis will be very high, but you need to
     // * call it only once. TODO

	public:
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
					for (int cx=px[sx]; cx < px[sx] + nrows * ncols; cx += nrows+1) {
						logDet += log(x[cx]);
					}
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

		template<typename MB>
		double inv_quad_form(const Eigen::MatrixBase<MB> &vec) {
			eigen_assert(m_factorizationIsOk && "The decomposition is not in a valid state for solving, you must first call either compute() or symbolic()/numeric()");
			if (ident.rows() != vec.rows()) {
				ident.setIdentity(vec.rows(), vec.rows());
			}
			cholmod_dense b_cd(viewAsCholmod(ident));
			cholmod_dense* x_cd = cholmod_solve(CHOLMOD_A, factor(), &b_cd, &cholmod());
			if(!x_cd) throw std::runtime_error("cholmod_solve failed");
			Eigen::Map< Eigen::MatrixXd > iA((double*) x_cd->x, vec.rows(), vec.rows());
			double ans = vec.transpose() * iA.selfadjointView<Eigen::Lower>() * vec;
			cholmod_free_dense(&x_cd, &cholmod());
			return ans;
		};
	};

	struct state {
		omxMatrix *smallRow;
		int totalNotMissing;
		std::vector<bool> notMissing;
		Eigen::VectorXd data;
		omxMatrix *cov;
		omxMatrix *means;
		omxMatrix *smallCov;
		omxMatrix *smallMeans;
		Eigen::SparseMatrix<double> fullCov;
		Cholmod< Eigen::SparseMatrix<double> > covDecomp;
		Eigen::VectorXd fullMeans;
	};
	
	static void compute(omxFitFunction *oo, int want, FitContext *fc)
	{
		if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
		if (!(want & (FF_COMPUTE_FIT | FF_COMPUTE_INITIAL_FIT))) Rf_error("Not implemented");

		state *st                               = (state *) oo->argStruct;
		omxExpectation *expectation             = oo->expectation;
		omxData *data                           = expectation->data;
		omxMatrix *cov                          = st->cov;
		omxMatrix *means                        = st->means;
		Eigen::SparseMatrix<double> &fullCov    = st->fullCov;
		Eigen::VectorXd &fullMeans              = st->fullMeans;

		Eigen::VectorXi contRemove(cov->cols);
		Eigen::VectorXd oldDefs;
		oldDefs.resize(data->defVars.size());
		oldDefs.setConstant(NA_REAL);

		fullMeans.resize(st->totalNotMissing);
		fullCov.resize(st->totalNotMissing, st->totalNotMissing);
		fullCov.setZero();

		int filteredPos = 0;
		for (int row=0; row < data->rows; ++row) {
			int numVarsFilled = data->handleDefinitionVarList(oo->matrix->currentState, row, oldDefs.data());
			if (row == 0 || numVarsFilled) {
				omxExpectationCompute(expectation, NULL);
			}

			int fullPos = row * cov->rows;
			for (int dx=0; dx < cov->rows; ++dx) contRemove[dx] = !st->notMissing[fullPos + dx];
			omxCopyMatrix(st->smallMeans, means);
			omxRemoveElements(st->smallMeans, contRemove.data());
			omxCopyMatrix(st->smallCov, cov);
			omxRemoveRowsAndColumns(st->smallCov, contRemove.data(), contRemove.data());

			EigenVectorAdaptor smallMeans(st->smallMeans);
			EigenMatrixAdaptor smallCov(st->smallCov);

			fullMeans.segment(filteredPos, st->smallCov->rows) = smallMeans;
			for (int sr=0; sr < smallCov.rows(); ++sr) {
				for (int sc=0; sc <= sr; ++sc) {
					fullCov.coeffRef(filteredPos + sr, filteredPos + sc) = smallCov(sr,sc);
				}
			}
			filteredPos += st->smallCov->rows;
		}

		if(0){Eigen::MatrixXd tmp = fullCov.block(0,0,8,8);
			mxPrintMat("R", tmp);}

		for (size_t vx=0; vx < expectation->varying.size(); ++vx) {
			varyBy &vb = expectation->varying[vx];
			if (!omxDataColumnIsFactor(data, vb.factorCol)) {
				Rf_error("Column %s (%d) is not a factor",
					 omxDataColumnName(data, vb.factorCol), vb.factorCol);
			}
			omxMatrix *cov = omxGetExpectationComponent(vb.model, "unfilteredCov");
			omxMatrix *Zspec = vb.model->Zmatrix;

			int levels = omxDataGetNumFactorLevels(data, vb.factorCol);
			std::vector<bool> curLevelMask;
			curLevelMask.resize(data->rows);
			for (int lx=1; lx <= levels; ++lx) {
				int numAtLevel = 0;
				for (int rx=0; rx < data->rows; ++rx) {
					bool yes = omxIntDataElement(data, rx, vb.factorCol) == lx;
					curLevelMask[rx] = yes;
					numAtLevel += yes;
				}
				Eigen::MatrixXd Zmat(Zspec->rows * numAtLevel, Zspec->cols);
				Eigen::VectorXd voldDefs;
				voldDefs.resize(vb.model->data->defVars.size());
				voldDefs.setConstant(NA_REAL);
				int zr=0;
				for (int rx=0; rx < data->rows; ++rx) {
					if (!curLevelMask[rx]) continue;
					vb.model->data->handleDefinitionVarList(oo->matrix->currentState, rx, voldDefs.data());
					omxRecompute(Zspec, fc);
					EigenMatrixAdaptor eZspec(Zspec);
					//mxPrintMat("Zspec", eZspec);
					Zmat.block(zr, 0, Zspec->rows, Zspec->cols).array() = eZspec.array();
					zr += Zspec->rows;
				}
				EigenMatrixAdaptor ecov(cov);
				//mxPrintMat("Z", Zmat);
				//mxPrintMat("G", ecov);
				Eigen::MatrixXd ZGZ = Zmat * ecov.selfadjointView<Eigen::Lower>() * Zmat.transpose();
				//mxPrintMat("ZGZ", ZGZ.block(0,0,8,8));
				for (int r1=0,v1=0; r1 < data->rows; ++r1) {
					if (!curLevelMask[r1]) continue;
					for (int r2=0,v2=0; r2 <= r1; ++r2) {
						if (!curLevelMask[r2]) continue;
						for (int b1=0; b1 < Zspec->rows; ++b1) {
							for (int b2=0; b2 < Zspec->rows; ++b2) {
								double val = ZGZ(Zspec->rows * v1 + b1, Zspec->rows * v2 + b2);
								if (val == 0) continue;
								fullCov.coeffRef(Zspec->rows * r1 + b1,
										 Zspec->rows * r2 + b2) += val;
							}
						}
						++v2;
					}
					++v1;
				}
			}
		}

		if(0){Eigen::MatrixXd tmp = fullCov.block(0,0,8,8);
			mxPrintMat("V", tmp);}

		double lp = NA_REAL;
		try {
			st->covDecomp.analyzePattern(fullCov);
			st->covDecomp.factorize(fullCov);
			lp = st->covDecomp.log_determinant();
			Eigen::VectorXd resid = st->data - fullMeans;
			double iqf = st->covDecomp.inv_quad_form(resid);
			lp += iqf;
			lp += M_LN_2PI * st->totalNotMissing;
		} catch (const std::exception& e) {
			if (fc) fc->recordIterationError("%s: %s", oo->name(), e.what());
		}
		//mxLog("%f", lp);
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
		omxFreeMatrix(st->smallMeans);
		omxFreeMatrix(st->smallCov);
		omxFreeMatrix(st->smallRow);
		delete st;
	}
};

void InitFellnerFitFunction(omxFitFunction *oo)
{
	omxExpectation* expectation = oo->expectation;
	if(expectation == NULL) {
		omxRaiseErrorf("%s cannot fit without a model expectation", oo->fitType);
		return;
	}
	omxMatrix *cov = omxGetExpectationComponent(expectation, "cov");
	if(cov == NULL) { 
		omxRaiseError("No covariance expectation in FIML evaluation.");
		return;
	}

	omxMatrix *means = omxGetExpectationComponent(expectation, "means");
	if(means == NULL) { 
		omxRaiseError("No means model in FIML evaluation.");
		return;
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

	st->cov = cov;
	st->means = means;
	st->smallCov   = omxInitMatrix(1, 1, TRUE, oo->matrix->currentState);
	st->smallMeans = omxInitMatrix(1, 1, TRUE, oo->matrix->currentState);
	st->smallRow = omxInitMatrix(1, cov->cols, TRUE, oo->matrix->currentState);
	omxData *data               = expectation->data;
	omxMatrix *dataColumns	    = expectation->dataColumns;

	st->totalNotMissing = 0;
	st->notMissing.reserve(data->rows * cov->cols);
	for (int row=0; row < data->rows; ++row) {
		omxDataRow(data, row, dataColumns, st->smallRow);
		for (int col=0; col < cov->cols; ++col) {
			double val = omxMatrixElement(st->smallRow, 0, col);
			bool yes = std::isfinite(val);
			st->notMissing.push_back(yes);
			if (yes) ++st->totalNotMissing;
		}
	}

	st->data.resize(st->totalNotMissing);
	for (int row=0, dx=0; row < data->rows; ++row) {
		omxDataRow(data, row, dataColumns, st->smallRow);
		for (int col=0; col < cov->cols; ++col) {
			double val = omxMatrixElement(st->smallRow, 0, col);
			if (!std::isfinite(val)) continue;
			st->data[ dx++ ] = val;
		}
	}
}
