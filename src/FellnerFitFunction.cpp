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
#include <RcppEigenCholmod.h>
#include <RcppEigenStubs.h>
#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/CholmodSupport>
#include <Eigen/SparseLU>
//#include <Eigen/UmfPackSupport>
#include "omxFitFunction.h"
#include "RAMInternal.h"

namespace FellnerFitFunction {
	// Based on lme4CholmodDecomposition.h from lme4
	template<typename _MatrixType, int _UpLo = Eigen::Lower>
	class Cholmod : public Eigen::CholmodDecomposition<_MatrixType, _UpLo> {
	private:
		Eigen::MatrixXd ident;

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
     // * good permutation.  In this case, try setting Common->nmethods to 9.
     // * The time spent in cholmod_analysis will be very high, but you need to
     // * call it only once. TODO

	public:
		Cholmod() {
			oldHandler = cholmod().error_handler;
			cholmod().error_handler = cholmod_error;
		};
		~Cholmod() {
			cholmod().error_handler = oldHandler;
		};
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
		omxMatrix *smallCol;
		std::vector<bool> latentFilter; // use to reduce the A matrix
		bool AmatDependsOnParameters;
		bool haveFilteredAmat;
		Eigen::VectorXd dataVec;
		Eigen::SparseMatrix<double>      fullA;
		Eigen::SparseLU< Eigen::SparseMatrix<double>,
				 Eigen::COLAMDOrdering<Eigen::SparseMatrix<double>::Index> > Asolver;
		//Eigen::UmfPackLU< Eigen::SparseMatrix<double> > Asolver;
		Eigen::SparseMatrix<double>      filteredA;
		Eigen::SparseMatrix<double>      fullS;
		Eigen::SparseMatrix<double>      fullCov;
		Cholmod< Eigen::SparseMatrix<double> > covDecomp;
		Eigen::VectorXd fullMeans;

		void loadOneRow(omxExpectation *expectation, FitContext *fc, int row, int &lx);
		void placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize);
		void prepOneRow(omxExpectation *expectation, int row_or_key, int &lx, int &dx);
	};
	
	// verify whether sparse can deal with parameters set to exactly zero TODO

	void state::loadOneRow(omxExpectation *expectation, FitContext *fc, int key_or_row, int &lx)
	{
		omxData *data = expectation->data;

		int row;
		if (!data->hasPrimaryKey()) {
			row = key_or_row;
		} else {
			row = data->lookupRowOfKey(key_or_row);
			if (data->rowToOffsetMap[row] != lx) return;
		}

		data->handleDefinitionVarList(expectation->currentState, row);
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
		omxRecompute(ram->A, fc);
		omxRecompute(ram->S, fc);
		if (ram->M) omxRecompute(ram->M, fc);

		for (size_t jx=0; jx < ram->joins.size(); ++jx) {
			join &j1 = ram->joins[jx];
			int key = omxKeyDataElement(data, row, j1.foreignKey);
			if (key == NA_INTEGER) continue;
			loadOneRow(j1.ex, fc, key, lx);
		}
		for (size_t jx=0; jx < ram->joins.size(); ++jx) {
			join &j1 = ram->joins[jx];
			int key = omxKeyDataElement(data, row, j1.foreignKey);
			if (key == NA_INTEGER) continue;
			int frow = j1.data->lookupRowOfKey(key);
			int jOffset = j1.data->rowToOffsetMap[frow];
			omxMatrix *betA = j1.regression;
			omxRecompute(betA, fc);
			omxRAMExpectation *ram2 = (omxRAMExpectation*) j1.ex->argStruct;
			for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
				for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
					fullA.coeffRef(lx + rx, jOffset + cx) -=
						omxMatrixElement(betA, rx, cx);
				}
			}
		}

		if (!haveFilteredAmat) {
			EigenMatrixAdaptor eA(ram->A);
			for (int cx=0; cx < eA.cols(); ++cx) {
				for (int rx=0; rx < eA.rows(); ++rx) {
					if (rx != cx && eA(rx,cx) != 0) {
						// can't use eA.block(..) -= because fullA must remain sparse
						fullA.coeffRef(lx + rx, lx + cx) -= eA(rx, cx);
					}
				}
			}
		}

		EigenMatrixAdaptor eS(ram->S);
		for (int cx=0; cx < eS.cols(); ++cx) {
			for (int rx=0; rx < eS.rows(); ++rx) {
				if (rx >= cx && eS(rx,cx) != 0) {
					fullS.coeffRef(lx + rx, lx + cx) += eS(rx, cx);
				}
			}
		}

		if (ram->M) {
			EigenVectorAdaptor eM(ram->M);
			for (int mx=0; mx < eM.size(); ++mx) {
				fullMeans[lx + mx] = eM[mx];
			}
		} else {
			fullMeans.segment(lx, ram->A->cols).setZero();
		}

		lx += ram->A->cols;
	}

	static void compute(omxFitFunction *oo, int want, FitContext *fc)
	{
		if (want & (FF_COMPUTE_PREOPTIMIZE)) return;
		if (!(want & (FF_COMPUTE_FIT | FF_COMPUTE_INITIAL_FIT))) Rf_error("Not implemented");

		state *st                               = (state *) oo->argStruct;
		omxExpectation *expectation             = oo->expectation;
		omxData *data                           = expectation->data;
		Eigen::SparseMatrix<double> &fullA      = st->fullA;
		Eigen::SparseMatrix<double> &filteredA  = st->filteredA;
		Eigen::SparseMatrix<double> &fullS      = st->fullS;
		Eigen::SparseMatrix<double> &fullCov    = st->fullCov;
		Eigen::VectorXd &fullMeans              = st->fullMeans;

		fullMeans.resize(st->latentFilter.size());

		if (fullA.nonZeros() == 0) {
			fullA.resize(st->latentFilter.size(), st->latentFilter.size());
			fullA.setIdentity();
		} else {
			for (int k=0; k<fullA.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(fullA, k); it; ++it) {
					it.valueRef() = it.row() == it.col()? 1 : 0;
				}
			}
		}
		if (fullS.nonZeros() == 0) {
			fullS.resize(st->latentFilter.size(), st->latentFilter.size());
		} else {
			for (int k=0; k<fullS.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(fullS, k); it; ++it) {
					it.valueRef() = 0;
				}
			}
		}

		for (int lx=0, row=0; row < data->rows; ++row) {
			int key_or_row = data->hasPrimaryKey()? data->primaryKeyOfRow(row) : row;
			st->loadOneRow(expectation, fc, key_or_row, lx);
		}

		//{ Eigen::MatrixXd tmp = fullA; mxPrintMat("fullA", tmp); }
		//{ Eigen::MatrixXd tmp = fullS; mxPrintMat("fullS", tmp); }

		double lp = NA_REAL;
		try {
			if (!fullA.isCompressed()) {
				fullA.makeCompressed();
				st->Asolver.analyzePattern(fullA);
			}
			st->Asolver.factorize(fullA);
			if (!st->haveFilteredAmat) {
				// consider http://users.clas.ufl.edu/hager/papers/Lightning/update.pdf ?

				if (1) {
					Eigen::SparseMatrix<double> ident;
					ident.resize(fullA.rows(), fullA.rows());
					ident.setIdentity();
					Eigen::SparseMatrix<double> result = st->Asolver.solve(ident);
					//{ Eigen::MatrixXd tmp = result; mxPrintMat("result", tmp); }
					filteredA.resize(st->dataVec.size(), fullA.rows());
					filteredA.setZero();
					for (int rx=0, dx=-1; rx < fullA.rows(); ++rx) {
						if (!st->latentFilter[rx]) continue;
						++dx;
						for (int cx=0; cx < fullA.cols(); ++cx) {
							if (result.coeff(rx, cx) == 0) continue;
							filteredA.coeffRef(dx, cx) = result.coeff(rx, cx);
						}
					}
				} else {
					// at least 3x slower than inverting the whole matrix in one go
					filteredA.resize(st->dataVec.size(), fullA.rows());
					filteredA.reserve(st->dataVec.size());
					Eigen::VectorXd a1(fullA.rows());
					a1.setZero();
					Eigen::VectorXd result(fullA.rows());
					for (int ax=0; ax < fullA.rows(); ++ax) {
						// is there a faster way to do this? TODO
						a1[ax] = 1.0;
						result = st->Asolver.solve(a1);
						for (int lx=0, ox=-1; lx < fullA.rows(); ++lx) {
							if (!st->latentFilter[lx]) continue;
							++ox;
							if (result[lx] == 0) continue;
							filteredA.coeffRef(ox, ax) = result[lx];
						}
						a1[ax] = 0.0;
					}
				}
				filteredA.makeCompressed();
				st->haveFilteredAmat = !st->AmatDependsOnParameters;
				//{ Eigen::MatrixXd tmp = filteredA; mxPrintMat("filteredA", tmp); }
			}
			//mxPrintMat("S", fullS);
			//Eigen::MatrixXd fullCovDense =
			bool firstTime = fullCov.nonZeros() == 0;
			fullCov = (filteredA * fullS.selfadjointView<Eigen::Lower>() * filteredA.transpose());
			//{ Eigen::MatrixXd tmp = fullCov; mxPrintMat("fullcov", tmp); }

			if (firstTime) {
				fullCov.makeCompressed();
				st->covDecomp.analyzePattern(fullCov);
			}

			st->covDecomp.factorize(fullCov);
			lp = st->covDecomp.log_determinant();
			//mxPrintMat("dataVec", st->dataVec);
			//mxPrintMat("fullMeans", fullMeans);
			Eigen::VectorXd resid = st->dataVec - filteredA * fullMeans;
			double iqf = st->covDecomp.inv_quad_form(resid);
			lp += iqf;
			lp += M_LN_2PI * st->dataVec.size();
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
		omxFreeMatrix(st->smallCol);
		delete st;
	}

	void state::placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize)
	{
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		for (size_t jx=0; jx < ram->joins.size(); ++jx) {
			join &j1 = ram->joins[jx];
			int key = omxKeyDataElement(data, frow, j1.foreignKey);
			if (key == NA_INTEGER) continue;
			placeOneRow(j1.ex, j1.data->lookupRowOfKey(key), totalObserved, maxSize);
		}
		if (data->hasPrimaryKey()) {
			// insert_or_assign would be nice here
			std::map<int,int>::const_iterator it = data->rowToOffsetMap.find(frow);
			if (it != data->rowToOffsetMap.end()) return;

			if (OMX_DEBUG) {
				mxLog("%s: place row %d at %d", expectation->name, frow, maxSize);
			}
			data->rowToOffsetMap[frow] = maxSize;
		}
		int jCols = expectation->dataColumns->cols;
		if (jCols) {
			if (!ram->M) {
				Rf_error("'%s' has manifest observations but '%s' has no mean model",
					 data->name, expectation->name);
			}
			if (smallCol->cols < jCols) {
				omxResizeMatrix(smallCol, 1, jCols);
			}
			omxDataRow(expectation, frow, smallCol);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(smallCol, 0, col);
				bool yes = std::isfinite(val);
				if (yes) ++totalObserved;
			}
		}
		maxSize += ram->F->cols;
		AmatDependsOnParameters |= ram->A->dependsOnParameters();
	}

	void state::prepOneRow(omxExpectation *expectation, int row_or_key, int &lx, int &dx)
	{
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		int frow;
		if (!data->hasPrimaryKey()) {
			frow = row_or_key;
		} else {
			frow = data->lookupRowOfKey(row_or_key);
			if (data->rowToOffsetMap[frow] != lx) return;
		}

		for (size_t jx=0; jx < ram->joins.size(); ++jx) {
			join &j1 = ram->joins[jx];
			int key = omxKeyDataElement(data, frow, j1.foreignKey);
			if (key == NA_INTEGER) continue;
			prepOneRow(j1.ex, key, lx, dx);
		}

		int jCols = expectation->dataColumns->cols;
		if (jCols) {
			omxDataRow(expectation, frow, smallCol);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(smallCol, 0, col);
				bool yes = std::isfinite(val);
				if (!yes) continue;
				latentFilter[ lx + col ] = true;
				dataVec[ dx++ ] = val;
			}
		}
		lx += ram->F->cols;
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

		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
		ram->ensureTrivialF();
		int numManifest = ram->F->rows;

		st->AmatDependsOnParameters = ram->A->dependsOnParameters();
		st->haveFilteredAmat = false;
		st->smallCol = omxInitMatrix(1, numManifest, TRUE, oo->matrix->currentState);
		omxData *data               = expectation->data;

		// don't permit reuse of our expectations by some other fit function TODO

		int totalObserved = 0;
		int maxSize = 0;
		for (int row=0; row < data->rows; ++row) {
			st->placeOneRow(expectation, row, totalObserved, maxSize);
		}
		//mxLog("AmatDependsOnParameters=%d", st->AmatDependsOnParameters);

		//mxLog("total observations %d", totalObserved);
		st->latentFilter.assign(maxSize, false); // will have totalObserved true entries
		st->dataVec.resize(totalObserved);
	
		for (int row=0, dx=0, lx=0; row < data->rows; ++row) {
			int key_or_row = data->hasPrimaryKey()? data->primaryKeyOfRow(row) : row;
			st->prepOneRow(expectation, key_or_row, lx, dx);
		}
	}
};

void InitFellnerFitFunction(omxFitFunction *oo)
{
	FellnerFitFunction::init(oo);
}
