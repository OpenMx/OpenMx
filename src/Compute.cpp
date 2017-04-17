/*
 *  Copyright 2007-2017 The OpenMx Project
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

/*File created 2013*/

#include <algorithm>
#include <stdarg.h>
#include <limits>
#include <random>
//#include <iostream>

#include "glue.h"
#include "Compute.h"
#include "omxState.h"
#include "omxExportBackendState.h"
#include "omxRFitFunction.h"
#include "matrix.h"
#include "omxBuffer.h"
#include "omxState.h"
#include <Eigen/Cholesky>
#include "finiteDifferences.h"
#include "EnableWarnings.h"

void pda(const double *ar, int rows, int cols);

SEXP allocInformVector(int size)
{
	const char *statusCodeLabels[] = { // see ComputeInform type
		"OK", "OK/green",
		"infeasible linear constraint",
		"infeasible non-linear constraint",
		"iteration limit",
		"not convex",
		"nonzero gradient",
		"bad deriv",
		"?",
		"internal error",
		"infeasible start"
	};
	return makeFactor(Rf_allocVector(INTSXP, size),
			  OMX_STATIC_ARRAY_SIZE(statusCodeLabels), statusCodeLabels);
}

void FitContext::queue(HessianBlock *hb)
{
	if (hb->vars.size() == 0) {
		delete hb;
		return;
	}

	// std::cout << "queue " << allBlocks.size() << " var map: ";
	// for (int vm=0; vm < (int) hb->vars.size(); ++vm) {
	// 	std::cout << vm << ":" << hb->vars[vm] << ",";
	// }
	// std::cout << "\n" << hb->mat << "\n";

	if (OMX_DEBUG) {
		std::vector<int>::iterator it = std::unique(hb->vars.begin(), hb->vars.end());
		if (std::distance(hb->vars.end(),it) != 0) {
			Rf_error("HessianBlock var mapping is not 1-to-1");
		}
		int prev = hb->vars[0];
		for (int vx=1; vx < int(hb->vars.size()); vx++) {
			if (prev > hb->vars[vx]) {
				Rf_error("hb->vars must be sorted");
			}
			prev = hb->vars[vx];
		}
	}

	minBlockSize = std::max(int(hb->vars.size()), minBlockSize);
	allBlocks.push_back(hb);
}

void FitContext::analyzeHessianBlock(HessianBlock *hb)
{
	for (size_t vx=0; vx < hb->vars.size(); ++vx) {
		HessianBlock *hb2 = blockByVar[ hb->vars[vx] ];
		if (!hb2) continue;

		for (size_t vx2=0; vx2 < hb2->vars.size(); ++vx2) {
			blockByVar[ hb2->vars[vx2] ] = NULL;
		}
		estNonZero -= hb2->estNonZero();

		if (hb->vars.size() < hb2->vars.size()) std::swap(hb, hb2);

		// std::cout << "hb: ";
		// for (int i : hb->vars) std::cout << i << ' ';
		// std::cout << "\nhb2: ";
		// for (int i : hb2->vars) std::cout << i << ' ';
		// std::cout << "\n";
		bool inc = std::includes(hb->vars.begin(), hb->vars.end(), hb2->vars.begin(), hb2->vars.end());
		//std::cout << "includes " << inc << "\n";
		if (inc) {
			hb->subBlocks.push_back(hb2);
		} else {
			HessianBlock *par2 = new HessianBlock;
			mergeBlocks.push_back(par2);
			par2->merge = true;
			par2->vars = hb->vars;
			par2->vars.insert(par2->vars.end(), hb2->vars.begin(), hb2->vars.end());
			std::inplace_merge(par2->vars.begin(), par2->vars.begin() + hb->vars.size(), par2->vars.end());
			int nsize = std::unique(par2->vars.begin(), par2->vars.end()) - par2->vars.begin();
			par2->vars.resize(nsize);
			par2->mat.resize(nsize, nsize);
			par2->mat.triangularView<Eigen::Upper>().setZero();
			par2->subBlocks.push_back(hb);
			par2->subBlocks.push_back(hb2);
			analyzeHessianBlock(par2);
			return;
		}
	}

	for (size_t vx=0; vx < hb->vars.size(); ++vx) {
		blockByVar[ hb->vars[vx] ] = hb;
	}
	estNonZero += hb->estNonZero();
	maxBlockSize = std::max(int(hb->vars.size()), maxBlockSize);
	//mxLog("queue maxBlocksize %d sparse %f", maxBlockSize, estNonZero / double(numParam * numParam));
}

void FitContext::analyzeHessian()
{
	// If we knew the minBlockSize was large then we wouldn't even
	// try to build merge blocks. 
	// If maxBlockSize is greater than some threshold then we should
	// also give up.

	if (blockByVar.size()) return;

	blockByVar.assign(numParam, NULL);

	for (size_t bx=0; bx < allBlocks.size(); ++bx) {
		analyzeHessianBlock(allBlocks[bx]);
	}
}

void FitContext::negateHessian()
{
	// this assumes that we haven't processed the blocks yet, need to check TODO
	for (size_t bx=0; bx < allBlocks.size(); ++bx) {
		allBlocks[bx]->mat *= -1.0;
	}
}

void FitContext::refreshDenseHess()
{
	if (haveDenseHess) return;

	hess.triangularView<Eigen::Upper>().setZero();

	for (size_t bx=0; bx < allBlocks.size(); ++bx) {
		HessianBlock *hb = allBlocks[bx];

		std::vector<int> &map = hb->vars;
		size_t bsize = map.size();

		for (size_t v1=0; v1 < bsize; ++v1) {
			for (size_t v2=0; v2 <= v1; ++v2) {
				hess(map[v2], map[v1]) += hb->mat(v2, v1);
			}
		}
	}

	haveDenseHess = true;
}

void FitContext::copyDenseHess(double *dest)
{
	refreshDenseHess();
	for (size_t v1=0; v1 < numParam; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			double coef = hess.selfadjointView<Eigen::Upper>()(v2,v1);
			if (v1==v2) {
				dest[v1 * numParam + v2] = coef;
			} else {
				dest[v1 * numParam + v2] = coef;
				dest[v2 * numParam + v1] = coef;
			}
		}
	}
}

double *FitContext::getDenseHessUninitialized()
{
	// Assume the caller is going to fill it out
	haveDenseHess = true;
	haveDenseIHess = false;
	return hess.data();
}

struct allFiniteHelper {
	typedef double Scalar;
	typedef Eigen::DenseIndex  Index;
	bool finite;
	void init(const Scalar& value, Index i, Index j) {
		finite = true;
		if (i <= j) finite = std::isfinite(value);
	};
	void operator() (const Scalar& value, Index i, Index j) {
		if (i <= j) finite &= std::isfinite(value);
	};
};

static void InvertSymmetricNR(Eigen::MatrixXd &hess, Eigen::MatrixXd &ihess)
{
	ihess = hess;
	Matrix ihessMat(ihess.data(), ihess.rows(), ihess.cols());

	if (0) {
		// for comparison
		InvertSymmetricIndef(ihessMat, 'U');
		return;
	}

	if (!InvertSymmetricPosDef(ihessMat, 'U')) return;

	int numParams = hess.rows();
	allFiniteHelper afh;
	hess.visit(afh);
	if (!afh.finite) {
		ihess = Eigen::MatrixXd::Zero(numParams, numParams);
		return;
	}

	omxBuffer<double> hessWork(numParams * numParams);
	memcpy(hessWork.data(), hess.data(), sizeof(double) * numParams * numParams);
	//pda(hess.data(), numParams, numParams);

	char jobz = 'V';
	char range = 'A';
	char uplo = 'U';
	double abstol = 0;
	int m;
	omxBuffer<double> w(numParams);
	omxBuffer<double> z(numParams * numParams);
	double optWork;
	int optIwork;
	int lwork = -1;
	int liwork = -1;
	int info;
	double realIgn = 0;
	int intIgn = 0;
	omxBuffer<int> isuppz(numParams * 2);

	F77_CALL(dsyevr)(&jobz, &range, &uplo, &numParams, hessWork.data(),
			 &numParams, &realIgn, &realIgn, &intIgn, &intIgn, &abstol, &m, w.data(),
			 z.data(), &numParams, isuppz.data(), &optWork, &lwork, &optIwork, &liwork, &info);

	lwork = optWork;
	omxBuffer<double> work(lwork);
	liwork = optIwork;
	omxBuffer<int> iwork(liwork);
	F77_CALL(dsyevr)(&jobz, &range, &uplo, &numParams, hessWork.data(),
			 &numParams, &realIgn, &realIgn, &intIgn, &intIgn, &abstol, &m, w.data(),
			 z.data(), &numParams, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info);
	if (info < 0) {
		Rf_error("dsyevr %d", info);
	} else if (info) {
		mxLog("Eigen decomposition failed %d", info);
		ihess = Eigen::MatrixXd::Zero(numParams, numParams);
		return;
	}

	std::vector<double> evalDiag(numParams * numParams);
	for (int px=0; px < numParams; ++px) {
		evalDiag[px * numParams + px] = 1/fabs(w[px]);
	}

	Matrix evMat(z.data(), numParams, numParams);
	Matrix edMat(evalDiag.data(), numParams, numParams);
	omxBuffer<double> prod1(numParams * numParams);
	Matrix p1Mat(prod1.data(), numParams, numParams);
	SymMatrixMultiply('R', 'U', 1.0, 0, edMat, evMat, p1Mat);
	char transa = 'N';
	char transb = 'T';
	double alpha = 1.0;
	double beta = 0;
	F77_CALL(dgemm)(&transa, &transb, &numParams, &numParams, &numParams, &alpha,
			prod1.data(), &numParams, z.data(), &numParams, &beta, ihess.data(), &numParams);
	//pda(ihess.data(), numParams, numParams);
}

void FitContext::refreshDenseIHess()
{
	if (haveDenseIHess) return;

	refreshDenseHess();

	ihess = hess;
	Matrix ihessMat(ihess.data(), ihess.rows(), ihess.cols());
	InvertSymmetricIndef(ihessMat, 'U');

	haveDenseIHess = true;
}

void FitContext::refreshSparseHess()
{
	if (haveSparseHess) return;

	sparseHess.resize(numParam, numParam);
	sparseHess.setZero();

	// need to sort allBlocks for performance TODO
	for (size_t bx=0; bx < allBlocks.size(); ++bx) {
		HessianBlock *hb = allBlocks[bx];

		std::vector<int> &map = hb->vars;
		size_t bsize = map.size();

		for (size_t v1=0; v1 < bsize; ++v1) {
			for (size_t v2=0; v2 <= v1; ++v2) {
				sparseHess.coeffRef(map[v2], map[v1]) += hb->mat(v2, v1);
			}
		}
	}

	haveSparseHess = true;
}

static double norm1(const Eigen::SparseMatrix<double> &mat)
{
	// norm_1, maximum column sum
	double maxCol = 0;
	for (int k=0; k < mat.outerSize(); ++k) {
		double col = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
			col += fabs(it.value());
		}
		maxCol = std::max(maxCol, col);
	}
	return maxCol;
}

static double frobeniusNorm(const Eigen::SparseMatrix<double> &mat)
{
	double total=0;
	for (int k=0; k < mat.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it) {
			total += it.value() * it.value();
		}
	}
	return total / mat.nonZeros();
}

static bool soleymani2013(const Eigen::SparseMatrix<double> &mat, Eigen::SparseMatrix<double> &imat)
{
       imat.setIdentity();
       imat /= frobeniusNorm(mat);

       Eigen::SparseMatrix<double> I1(mat.rows(), mat.cols());
       I1.setIdentity();
       Eigen::SparseMatrix<double> I7(mat.rows(), mat.cols());
       I7 = I1 * 7;
       Eigen::SparseMatrix<double> I21(mat.rows(), mat.cols());
       I21 = I1 * 21;
       Eigen::SparseMatrix<double> I35(mat.rows(), mat.cols());
       I35 = I1 * 35;

       double maxCol;
       int iter = 0;
       Eigen::SparseMatrix<double> AV = mat.selfadjointView<Eigen::Upper>() * imat;
       while (++iter < 100) {
               // try prune() on the multiplications TODO
               // try selfadjointView TODO
               Eigen::SparseMatrix<double> rpart =
                       (AV * (AV * (AV * (AV * (AV * (AV - I7) + I21) -I35) + I35) -I21) + I7);

               Eigen::SparseMatrix<double> prev = imat;
               imat = prev * rpart;
               AV = mat.selfadjointView<Eigen::Upper>() * imat;

               Eigen::SparseMatrix<double> diff = I1 - AV;
               maxCol = norm1(diff);
               if (maxCol < 1e-6) {
                       break;
               }
       }
       //mxLog("soleymani2013 invert in %d iter, maxCol %.2g", iter, maxCol);
       return false;
}

SEXP sparseInvert_wrapper(SEXP Rmat)
{
	omxManageProtectInsanity mpi;

	SEXP matrixDims;
	Rf_protect(matrixDims = Rf_getAttrib(Rmat, R_DimSymbol));
	int *dimList = INTEGER(matrixDims);
	int rows = dimList[0];
	int cols = dimList[1];
	if (rows != cols) Rf_error("Must be square");

	double *matData = REAL(Rmat);

	Eigen::SparseMatrix<double> mat(rows,cols);
	for (int cx=0; cx < cols; ++cx) {
		for (int rx=0; rx < rows; ++rx) {
			double val = matData[cx * rows + rx];
			if (val == 0) continue;
			mat.coeffRef(rx, cx) = val;
		}
	}

	Eigen::SparseMatrix<double> imat(rows,cols);
	if (soleymani2013(mat, imat)) Rf_error("Invert failed");
	
	SEXP ret;
	Rf_protect(ret = Rf_allocMatrix(REALSXP, rows, cols));
	double *retData = REAL(ret);
	for (int cx=0; cx < cols; ++cx) {
		for (int rx=0; rx < rows; ++rx) {
			retData[cx * rows + rx] = imat.coeff(rx, cx);
		}
	}

	return ret;
}

void FitContext::testMerge()
{
	const int UseId = 2;
	
	analyzeHessian();

	//std::cout << "block count " << allBlocks.size() << std::endl;

	sparseHess.resize(numParam, numParam);
	sparseHess.setZero();

	for (size_t vx=0; vx < numParam; ++vx) {
		HessianBlock *hb = blockByVar[vx];
		hb->addSubBlocks();
	}

	for (size_t vx=0; vx < numParam; ++vx) {
		HessianBlock *hb = blockByVar[vx];
		if (hb->useId == UseId) continue;
		hb->useId = UseId;

		// std::cout << "add block " << vx << "\nvar map: ";
		// for (int vm=0; vm < (int) hb->vars.size(); ++vm) {
		// 	std::cout << vm << ":" << hb->vars[vm] << ",";
		// }
		// std::cout << "\n" << hb->mmat << "\n";

		size_t size = hb->mmat.rows();
		for (size_t col=0; col < size; ++col) {
			for (size_t row=0; row <= col; ++row) {
				int vr = hb->vars[row];
				int vc = hb->vars[col];
				sparseHess.coeffRef(vr,vc) = hb->mmat(row,col);
			}
		}
	}

	refreshDenseHess();
	Eigen::MatrixXd dense = sparseHess;
	Eigen::MatrixXd diff = (dense - hess).selfadjointView<Eigen::Upper>();
	// std::cout << "difference\n" << diff << std::endl;
	// std::cout << "sparse\n" << dense << std::endl;
	// std::cout << "dense\n" << hess << std::endl;
	double bad = diff.cwiseAbs().maxCoeff();
	if (bad > .0001) Rf_error("Hess: dense sparse mismatch %f", bad);
}

bool FitContext::refreshSparseIHess()
{
	if (haveSparseIHess) return true;

	const int AcceptableDenseInvertSize = 100;
	const bool checkResult = OMX_DEBUG;

	//testMerge();

	sparseIHess.resize(numParam, numParam);
	sparseIHess.setZero();

	// sparseness by size simulation
	//mxLog("minBlockSize %d maxBlocksize %d sparse %f", maxBlockSize, estNonZero / double(numParam * numParam));

	if (minBlockSize < AcceptableDenseInvertSize) {
		analyzeHessian();
	}
	if (maxBlockSize < std::min(int(numParam), AcceptableDenseInvertSize)) {
		const int UseId = 1;
		for (size_t vx=0; vx < numParam; ++vx) {
			HessianBlock *hb = blockByVar[vx];
			if (!hb) Rf_error("Attempting to invert Hessian, "
					  "but no Hessian information for '%s'", varGroup->vars[vx]->name);
			if (hb->useId == UseId) continue;
			hb->useId = UseId;

			hb->addSubBlocks();
			size_t size = hb->mmat.rows();

			InvertSymmetricNR(hb->mmat, hb->imat);
		
			for (size_t col=0; col < size; ++col) {
				for (size_t row=0; row <= col; ++row) {
					int vr = hb->vars[row];
					int vc = hb->vars[col];
					sparseIHess.coeffRef(vr,vc) = hb->imat(row,col);
				}
			}
		}
	} else {
		return false;

		// Needs more work TODO
		// if (estNonZero / double(numParam * numParam) < .2)
		refreshSparseHess();
		if (soleymani2013(sparseHess, sparseIHess)) {
			sparseIHess.setZero();  // NR will try steepest descent
		}
	}

	if (checkResult) {
		refreshDenseHess();
		InvertSymmetricNR(hess, ihess);
		ihess.triangularView<Eigen::Lower>() = ihess.transpose().triangularView<Eigen::Lower>();
		Eigen::MatrixXd denseI = sparseIHess;
		denseI.triangularView<Eigen::Lower>() = denseI.transpose().triangularView<Eigen::Lower>();
		Eigen::MatrixXd resid = ihess - denseI;
		double bad = resid.cwiseAbs().maxCoeff();
		if (bad > .01) {
			// std::cout << "dense\n" << ihess << std::endl;
			// std::cout << "sparse\n" << denseI << std::endl;
			Rf_error("IHess: dense sparse mismatch %f", bad);
		}

		if (0) {
			int nonZero = 0;
			for (size_t cx=0; cx < numParam; ++cx) {
				for (size_t rx=0; rx <= cx; ++rx) {
					if (hess(rx,cx) != 0) {
						if (rx == cx) nonZero += 1;
						else nonZero += 2;
					}
				}
			}
			mxLog("done %f maxBlocksize %d est sparse %f actual sparse %f", bad, maxBlockSize,
			      estNonZero / double(numParam * numParam),
			      nonZero / double(numParam * numParam));
		}
	}

        haveSparseIHess = true;
	return true;
}

Eigen::VectorXd FitContext::ihessGradProd()
{
	if (refreshSparseIHess()) {
		return sparseIHess.selfadjointView<Eigen::Upper>() * grad;
	} else {
		refreshDenseHess();
		InvertSymmetricNR(hess, ihess);
		return ihess.selfadjointView<Eigen::Upper>() * grad;
	}
}

Eigen::VectorXd FitContext::ihessDiag()
{
	refreshDenseIHess();
	return ihess.diagonal();
}

double *FitContext::getDenseIHessUninitialized()
{
	// Assume the caller is going to fill it out
	haveDenseIHess = true;
	haveDenseHess = false;
	return ihess.data();
}

void FitContext::copyDenseIHess(double *dest)
{
	refreshDenseIHess();

	for (size_t v1=0; v1 < numParam; ++v1) {
		for (size_t v2=0; v2 <= v1; ++v2) {
			double coef = ihess.selfadjointView<Eigen::Upper>()(v2,v1);
			if (v1==v2) {
				dest[v1 * numParam + v2] = coef;
			} else {
				dest[v1 * numParam + v2] = coef;
				dest[v2 * numParam + v1] = coef;
			}
		}
	}
}

double *FitContext::getDenseHessianish()
{
	if (haveDenseHess) return hess.data();
	if (haveDenseIHess) return ihess.data();
	return NULL;
}

HessianBlock *HessianBlock::clone()
{
	HessianBlock *hb = new HessianBlock;
	hb->vars = vars;
	hb->mat.resize(vars.size(), vars.size());
	return hb;
}

void HessianBlock::addSubBlocks()
{
	if (mmat.rows()) return;

	mmat = mat;

	std::vector<int> vmap; // we could localize the vars map instead of computing the mapping each time TODO

	for (size_t bx=0; bx < subBlocks.size(); ++bx) {
		HessianBlock *sb = subBlocks[bx];
		sb->addSubBlocks();
	}

	//std::cout << "initial " << id << "\n" << mmat << std::endl;

	for (size_t bx=0; bx < subBlocks.size(); ++bx) {
		HessianBlock *sb = subBlocks[bx];
		
		//std::cout << "subblock " << sb->id << "\n" << sb->mmat << std::endl;

		size_t numVars = sb->vars.size();
		vmap.resize(numVars);
		for (size_t vx=0; vx < numVars; ++vx) {
			// std::set has special method for this TODO
			int px = std::lower_bound(vars.begin(), vars.end(), sb->vars[vx]) - vars.begin();
			vmap[vx] = px;
		}

		for (size_t v1=0; v1 < numVars; ++v1) {
			for (size_t v2=0; v2 <= v1; ++v2) {
				mmat(vmap[v2], vmap[v1]) += sb->mmat(v2, v1);
			}
		}
	}
	
	//std::cout << "result " << id << "\n" << mmat << std::endl;
}

int HessianBlock::estNonZero() const
{
	if (!merge) {
		return int(vars.size() * vars.size());
	} else {
		int total = 0;
		for (size_t bx=0; bx < subBlocks.size(); ++bx) {
			HessianBlock *sb = subBlocks[bx];
			total += sb->estNonZero();
		}
		return std::min(total, int(vars.size() * vars.size()));
	}
}

void FitContext::init()
{
	numParam = varGroup->vars.size();
	wanted = 0;
	mac = parent? parent->mac : 0;
	fit = parent? parent->fit : NA_REAL;
	fitUnits = parent? parent->fitUnits : FIT_UNITS_UNINITIALIZED;
	skippedRows = 0;
	est = new double[numParam];
	infoDefinite = NA_LOGICAL;
	infoCondNum = NA_REAL;
	infoA = NULL;
	infoB = NULL;
	stderrs = NULL;
	inform = INFORM_UNINITIALIZED;
	iterations = 0;
	ciobj = 0;
	openmpUser = false;
	computeCount = 0;

	hess.resize(numParam, numParam);
	ihess.resize(numParam, numParam);
	clearHessian();
}

void FitContext::clearHessian()
{
	for (size_t bx=0; bx < mergeBlocks.size(); ++bx) {
		delete mergeBlocks[bx];
	}
	for (size_t bx=0; bx < allBlocks.size(); ++bx) {
		delete allBlocks[bx];
	}

	allBlocks.clear();
	mergeBlocks.clear();
	blockByVar.clear();
	haveSparseHess = false;
	haveSparseIHess = false;
	haveDenseHess = false;
	haveDenseIHess = false;
	estNonZero = 0;
	minBlockSize = 0;
	maxBlockSize = 0;
}

void FitContext::allocStderrs()
{
	if (stderrs) return;

	stderrs = new double[numParam];

	for (size_t px=0; px < numParam; ++px) {
		stderrs[px] = NA_REAL;
	}
}

FitContext::FitContext(omxState *_state, std::vector<double> &startingValues)
{
	parent = NULL;
	varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
	init();
	profiledOut.assign(numParam, false);

	state = _state;
	if (numParam) {
		if (startingValues.size() != numParam) {
			Rf_error("Got %d starting values for %d parameters",
				 startingValues.size(), numParam);
		}
		memcpy(est, startingValues.data(), sizeof(double) * numParam);
	}
}

FitContext::FitContext(FitContext *_parent, FreeVarGroup *_varGroup)
{
	this->parent = _parent;
	this->varGroup = _varGroup;
	init();

	state = parent->state;
	FreeVarGroup *src = parent->varGroup;
	FreeVarGroup *dest = varGroup;
	size_t dvars = varGroup->vars.size();
	if (dvars == 0) return;
	mapToParent.resize(dvars);
	profiledOut.resize(numParam);

	size_t d1 = 0;
	for (size_t s1=0; s1 < src->vars.size(); ++s1) {
		if (src->vars[s1] != dest->vars[d1]) continue;
		mapToParent[d1] = s1;
		est[d1] = parent->est[s1];
		profiledOut[d1] = parent->profiledOut[s1];
		if (++d1 == dvars) break;
	}
	if (d1 != dvars) Rf_error("Parent free parameter group (id=%d) is not a superset of %d",
			       src->id[0], dest->id[0]);


	wanted = parent->wanted;
	infoDefinite = parent->infoDefinite;
	infoCondNum = parent->infoCondNum;
	iterations = parent->iterations;
	ciobj = parent->ciobj;
}

void FitContext::updateParent()
{
	FreeVarGroup *src = varGroup;
	FreeVarGroup *dest = parent->varGroup;
	size_t svars = varGroup->vars.size();

	parent->wanted |= wanted;
	parent->fit = fit;
	parent->fitUnits = fitUnits;
	parent->skippedRows = skippedRows;
	parent->mac = mac;
	parent->infoDefinite = infoDefinite;
	parent->infoCondNum = infoCondNum;
	parent->iterations = iterations;

	// rewrite using mapToParent TODO

	if (svars > 0) {
		for (size_t d1=0, s1 = 0; d1 < dest->vars.size(); ++d1) {
			if (dest->vars[d1] != src->vars[s1]) continue;
			parent->est[d1] = est[s1];
			if (++s1 == svars) break;
		}
		if (stderrs) {
			parent->allocStderrs();
			for (size_t s1=0; s1 < src->vars.size(); ++s1) {
				parent->stderrs[mapToParent[s1]] = stderrs[s1];
			}
		}
	}
	
	// pda(est, 1, svars);
	// pda(parent->est, 1, dvars);
}

int FitContext::getGlobalComputeCount()
{
	FitContext *fc = this;
	if (fc->parent && fc->parent->childList.size()) {
		// Kids for multithreading only appear as leaves of tree
		fc = fc->parent;
	}
	int cc = fc->getLocalComputeCount();
	while (fc->parent) {
		fc = fc->parent;
		cc += fc->getLocalComputeCount();
	}
	return cc;
}

int FitContext::getLocalComputeCount()
{
	int cc = computeCount;
	for (size_t cx=0; cx < childList.size(); ++cx) {
		cc += childList[cx]->getLocalComputeCount();
	}
	return cc;
}

void FitContext::updateParentAndFree()
{
	updateParent();
	delete this;
}

void FitContext::log(int what)
{
	size_t count = varGroup->vars.size();
	std::string buf;
	if (what & FF_COMPUTE_MAXABSCHANGE) buf += string_snprintf("MAC: %.5f\n", mac);
	if (what & FF_COMPUTE_FIT) buf += string_snprintf("fit: %.5f (scale %f)\n", fit, Global->llScale);
	if (what & FF_COMPUTE_ESTIMATE) {
		buf += string_snprintf("est %d: c(", (int) count);
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.5f", est[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_GRADIENT) {
		buf += string_snprintf("gradient %d: c(", (int) count);
		for (size_t vx=0; vx < count; ++vx) {
			buf += string_snprintf("%.5f", grad[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_HESSIAN) {
		refreshDenseHess();
		buf += string_snprintf("hessian %d x %d: c(\n", (int) count, (int) count);
		for (size_t v1=0; v1 < count; ++v1) {
			for (size_t v2=0; v2 < count; ++v2) {
				double coef;
				if (v1 > v2) {
					coef = hess.selfadjointView<Eigen::Upper>()(v2,v1);
				} else {
					coef = hess.selfadjointView<Eigen::Upper>()(v1,v2);
				}
				buf += string_snprintf("%.5f", coef);
				if (v1 < count - 1 || v2 < count - 1) buf += ", ";
			}
			buf += "\n";
		}
		buf += ")\n";
	}
	mxLogBig(buf);
}

void FitContext::resetIterationError()
{
	IterationError.clear();
}

void FitContext::recordIterationError(const char* msg, ...)
{
	// Can avoid overhead of setting error if one is already set TODO
	va_list ap;
	va_start(ap, msg);
	IterationError = string_vsnprintf(msg, ap);
	va_end(ap);
}

std::string FitContext::getIterationError()
{
	if (childList.size()) {
		size_t tlen = 0;
		for (size_t cx=0; cx < childList.size(); ++cx) {
			tlen += childList[cx]->IterationError.size();
		}
		// Fit function may not have used parallel processing
		if (tlen == 0) return IterationError;

		std::string result;
		for (size_t cx=0; cx < childList.size(); ++cx) {
			auto &str = childList[cx]->IterationError;
			if (!str.size()) continue;
			result += string_snprintf("%d: %s\n", int(cx), str.c_str());
		}
		return result;
	} else {
		return IterationError;
	}
}

static void _fixSymmetry(const char *name, double *mat, size_t numParam, bool force)
{
	for (size_t h1=1; h1 < numParam; h1++) {
		for (size_t h2=0; h2 < h1; h2++) {
			if (!force && mat[h2 * numParam + h1] != 0) {
				omxRaiseErrorf("%s is not upper triangular", name);
				break;
			}
			mat[h2 * numParam + h1] = mat[h1 * numParam + h2];
		}
	}
}

static void omxRepopulateRFitFunction(omxFitFunction* oo, double* x, int n)
{
	omxRFitFunction* rFitFunction = (omxRFitFunction*)oo;

	SEXP theCall, estimate;

	ScopedProtect(estimate, Rf_allocVector(REALSXP, n));
	double *est = REAL(estimate);
	for(int i = 0; i < n ; i++) {
		est[i] = x[i];
	}

	{ScopedProtect p1(theCall, Rf_allocVector(LANGSXP, 4));

		// imxUpdateModelValues does not handle parameters with equality
		// constraints. This is a bug.
	SETCAR(theCall, Rf_install("imxUpdateModelValues"));
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->flatModel);
	SETCADDDR(theCall, estimate);

	R_Reprotect(rFitFunction->model = Rf_eval(theCall, R_GlobalEnv), rFitFunction->modelIndex);
	}

	omxMarkDirty(oo->matrix);
}

void FitContext::ensureParamWithinBox(bool nudge)
{
	if (OMX_DEBUG) mxLog("FitContext::ensureParamWithinBox(nudge=%d)", nudge);
	for (size_t px = 0; px < varGroup->vars.size(); ++px) {
		omxFreeVar *fv = varGroup->vars[px];
		if (nudge && est[px] == 0.0) {
			est[px] += 0.1;
		}
		if (fv->lbound > est[px]) {
			est[px] = fv->lbound + 1.0e-6;
		}
		if (fv->ubound < est[px]) {
			est[px] = fv->ubound - 1.0e-6;
		}
        }
}

void FitContext::copyParamToModel()
{
	copyParamToModelClean();
	varGroup->markDirty(state);
}

void copyParamToModelInternal(FreeVarGroup *varGroup, omxState *os, double *at)
{
	size_t numParam = varGroup->vars.size();

	for(size_t k = 0; k < numParam; k++) {
		omxFreeVar* freeVar = varGroup->vars[k];
		freeVar->copyToState(os, at[k]);
	}
}

void FitContext::copyParamToModelClean()
{
	if(numParam == 0) return;

	copyParamToModelInternal(varGroup, state, est);

	if (RFitFunction) omxRepopulateRFitFunction(RFitFunction, est, numParam);

	if (childList.size() == 0 || !openmpUser) return;

	for(size_t i = 0; i < childList.size(); i++) {
		memcpy(childList[i]->est, est, sizeof(double) * numParam);
		childList[i]->copyParamToModel();
	}
}

omxMatrix *FitContext::lookupDuplicate(omxMatrix* element)
{
	if (element == NULL) return NULL;
	return state->lookupDuplicate(element);
}
	
double *FitContext::take(int want)
{
	if (!(want & (wanted | FF_COMPUTE_ESTIMATE))) {
		Rf_error("Attempt to take %d but not available", want);
	}

	double *ret = NULL;
	switch(want) {
	case FF_COMPUTE_ESTIMATE:
		ret = est;
		est = NULL;
		break;
	default:
		Rf_error("Taking of %d is not implemented", want);
	}
	if (!ret) Rf_error("Attempt to take %d, already taken", want);
	return ret;
}

// Rethink this whole design TODO
// If we compute things blockwise then most of this can go away?
void FitContext::preInfo()
{
	size_t npsq = numParam * numParam;

	if (!infoA) infoA = new double[npsq];
	if (!infoB) infoB = new double[npsq];

	switch (infoMethod) {
	case INFO_METHOD_SANDWICH:
	case INFO_METHOD_MEAT:
		OMXZERO(infoB, npsq);
	case INFO_METHOD_BREAD:
		OMXZERO(infoA, npsq);
		break;
	case INFO_METHOD_HESSIAN:
		clearHessian();
		break;
	default:
		Rf_error("Unknown information matrix estimation method %d", infoMethod);
	}
}

void FitContext::postInfo()
{
	switch (infoMethod) {
	case INFO_METHOD_SANDWICH:{
		// move into FCDeriv TODO
		omxBuffer<double> work(numParam * numParam);
		Matrix amat(infoA, numParam, numParam);
		InvertSymmetricIndef(amat, 'U');
		_fixSymmetry("InfoB", infoB, numParam, false);
		Matrix bmat(infoB, numParam, numParam);
		Matrix wmat(work.data(), numParam, numParam);
		Matrix hmat(getDenseIHessUninitialized(), numParam, numParam);
		SymMatrixMultiply('L', 'U', 1, 0, amat, bmat, wmat);
		SymMatrixMultiply('R', 'U', 1, 0, amat, wmat, hmat);
		wanted |= FF_COMPUTE_IHESSIAN;
		break;}
	case INFO_METHOD_MEAT:{
		memcpy(getDenseHessUninitialized(), infoB, sizeof(double) * numParam * numParam); // avoid copy TODO
		wanted |= FF_COMPUTE_HESSIAN;
		break;}
	case INFO_METHOD_BREAD:{
		memcpy(getDenseHessUninitialized(), infoA, sizeof(double) * numParam * numParam); // avoid copy TODO
		wanted |= FF_COMPUTE_HESSIAN;
		break;}
	case INFO_METHOD_HESSIAN:
		if (Global->llScale > 0) negateHessian();
		wanted |= FF_COMPUTE_HESSIAN;
		break;
	default:
		Rf_error("Unknown information matrix estimation method %d", infoMethod);
	}
}

bool FitContext::isClone() const
{
	return state->isClone();
}

void FitContext::createChildren(omxMatrix *alg)
{
	if (Global->numThreads <= 1) {
		diagParallel(OMX_DEBUG, "FitContext::createChildren: max threads set to 1");
		return;
	}
	if (childList.size()) return;

	for(size_t j = 0; j < state->expectationList.size(); j++) {
		if (!state->expectationList[j]->canDuplicate) {
			diagParallel(OMX_DEBUG, "FitContext::createChildren: %s cannot be duplicated",
				     state->expectationList[j]->name);
			return;
		}
	}
	for(size_t j = 0; j < state->algebraList.size(); j++) {
		omxFitFunction *ff = state->algebraList[j]->fitFunction;
		if (!ff) continue;
		if (!ff->canDuplicate) {
			diagParallel(OMX_DEBUG, "FitContext::createChildren: %s cannot be duplicated",
				     state->algebraList[j]->name());
			return;
		}
		if (ff->openmpUser) {
			diagParallel(OMX_DEBUG, "FitContext::createChildren: %s is an OpenMP user",
				     state->algebraList[j]->name());
		}
		openmpUser |= ff->openmpUser;
	}

	diagParallel(OMX_DEBUG, "FitContext::createChildren: create %d FitContext for parallel processing; OpenMP user=%d",
		     Global->numThreads, openmpUser);

	int numThreads = Global->numThreads;

	childList.reserve(numThreads);

	for(int ii = 0; ii < numThreads; ii++) {
		//omxManageProtectInsanity mpi;
		FitContext *kid = new FitContext(this, varGroup);
		kid->state = new omxState(state);
		kid->state->initialRecalc(kid);
		omxAlgebraPreeval(alg, kid);
		childList.push_back(kid);
		//if (OMX_DEBUG) mxLog("Protect depth at line %d: %d", __LINE__, mpi.getDepth());
	}

	if (OMX_DEBUG) mxLog("FitContext::createChildren: done creating %d omxState", Global->numThreads);
}

void FitContext::destroyChildren()
{
	if (0 == childList.size()) return;
	IterationError = getIterationError();
	for (int cx=0; cx < int(childList.size()); ++cx) {
		delete childList[cx];
	}
	childList.clear();
}

FitContext::~FitContext()
{
	destroyChildren();

	if (parent) {
		parent->computeCount += computeCount;
		computeCount = 0;

		if (parent->state != state) delete state;
	}

	clearHessian();
	if (est) delete [] est;
	if (stderrs) delete [] stderrs;
	if (infoA) delete [] infoA;
	if (infoB) delete [] infoB;
}

omxFitFunction *FitContext::RFitFunction = NULL;

void FitContext::setRFitFunction(omxFitFunction *rff)
{
	if (rff) {
		Global->numThreads = 1;
		if (RFitFunction) {
			Rf_error("You can only create 1 MxRFitFunction per independent model");
		}
	}
	RFitFunction = rff;
}

void CIobjective::evalFit(omxFitFunction *ff, int want, FitContext *fc)
{
	omxFitFunctionCompute(ff, want, fc);
}

void CIobjective::checkSolution(FitContext *fc)
{
	if (fc->getInform() > INFORM_UNCONVERGED_OPTIMUM) return;

	if (getDiag() != DIAG_SUCCESS) {
		fc->setInform(INFORM_NONLINEAR_CONSTRAINTS_INFEASIBLE);
	}
}

class EMAccel {
protected:
	FitContext *fc;
	int numParam;
	std::vector<double> prevAdj1;
	std::vector<double> prevAdj2;
	int verbose;

public:
	EMAccel(FitContext *_fc, int _verbose) : fc(_fc), verbose(_verbose) {
		numParam = fc->varGroup->vars.size();
		prevAdj1.assign(numParam, 0);
		prevAdj2.resize(numParam);
		dir.resize(numParam);
	};
	virtual ~EMAccel() {};
	template <typename T> void recordTrajectory(std::vector<T> &prevEst) {
		prevAdj2 = prevAdj1;
		for (int px=0; px < numParam; ++px) {
			prevAdj1[px] = fc->est[px] - prevEst[px];
		}
	};

	Eigen::VectorXd dir;
	virtual bool calcDirection(bool major) = 0;
	virtual void recalibrate() = 0;
	virtual bool retry() { return false; };
};

class Ramsay1975 : public EMAccel {
	// Ramsay, J. O. (1975). Solving Implicit Equations in
	// Psychometric Data Analysis.  Psychometrika, 40(3), 337-360.

	double minCaution;
	double highWatermark;
	bool goingWild;
	void restart(bool myFault);  // no longer in use

public:
	double maxCaution;
	double caution;

	Ramsay1975(FitContext *fc, int verbose, double minCaution);
	virtual void recalibrate();
	virtual bool calcDirection(bool major);
};

Ramsay1975::Ramsay1975(FitContext *_fc, int _verbose, double _minCaution) :
	EMAccel(_fc, _verbose), minCaution(_minCaution)
{
	maxCaution = 0.0;
	caution = 0;
	highWatermark = std::max(0.5, caution);  // arbitrary guess

	if (verbose >= 2) {
		mxLog("Ramsay: %d parameters, caution %f, min caution %f",
		      numParam, caution, minCaution);
	}
}

bool Ramsay1975::calcDirection(bool major)
{
	for (int vx=0; vx < numParam; ++vx) {
		const double prevEst = fc->est[vx] - prevAdj1[vx];
		dir[vx] = ((1 - caution) * fc->est[vx] + caution * prevEst) - fc->est[vx];
	}
	return true;
}

void Ramsay1975::recalibrate()
{
	if (numParam == 0) return;

	double normPrevAdj2 = 0;
	double normAdjDiff = 0;
	std::vector<double> adjDiff(numParam);

	// The choice of norm is also arbitrary. Other norms might work better.
	for (int px=0; px < numParam; ++px) {
		adjDiff[px] = prevAdj1[px] - prevAdj2[px];
		normPrevAdj2 += prevAdj2[px] * prevAdj2[px];
	}

	for (int px=0; px < numParam; ++px) {
		normAdjDiff += adjDiff[px] * adjDiff[px];
	}
	if (normAdjDiff == 0) {
		return;
	}

	double ratio = sqrt(normPrevAdj2 / normAdjDiff);
	//if (verbose >= 3) mxLog("Ramsay[%d]: sqrt(%.5f/%.5f) = %.5f",
	// flavor, normPrevAdj2, normAdjDiff, ratio);

	double newCaution = 1 - (1-caution) * ratio;
	if (newCaution > .95) newCaution = .95;  // arbitrary guess
	if (newCaution < 0) newCaution /= 2;     // don't get overconfident
	if (newCaution < minCaution) newCaution = minCaution;
	if (newCaution < caution) {
		caution = newCaution/3 + 2*caution/3;  // don't speed up too fast, arbitrary ratio
	} else {
		caution = newCaution;
	}
	maxCaution = std::max(maxCaution, caution);
	goingWild = false;
	if (caution < highWatermark || (normPrevAdj2 < 1e-3 && normAdjDiff < 1e-3)) {
		if (verbose >= 3) mxLog("Ramsay: %.2f caution", caution);
	} else {
		if (verbose >= 3) {
			mxLog("Ramsay: caution %.2f > %.2f, extreme oscillation, restart recommended",
			      caution, highWatermark);
		}
		goingWild = true;
	}
	highWatermark += .02; // arbitrary guess
}

void Ramsay1975::restart(bool myFault)
{
	prevAdj1.assign(numParam, 0);
	prevAdj2.assign(numParam, 0);
	myFault |= goingWild;
	if (myFault) {
		highWatermark = 1 - (1 - highWatermark) * .5; // arbitrary guess
		caution = std::max(caution, highWatermark);   // arbitrary guess
		maxCaution = std::max(maxCaution, caution);
		highWatermark = caution;
	}
	if (numParam && verbose >= 3) {
		mxLog("Ramsay: restart%s with %.2f caution %.2f highWatermark",
		      myFault? " (my fault)":"", caution, highWatermark);
	}
}

class Varadhan2008 : public EMAccel {
	bool retried;
	double maxAlpha;
	double alpha;
	Eigen::Map< Eigen::VectorXd > rr;
	Eigen::VectorXd vv;

public:
	Varadhan2008(FitContext *_fc, int _verbose) :
		EMAccel(_fc, _verbose), rr(&prevAdj2[0], numParam), vv(numParam)
	{
		alpha = 0;
		maxAlpha = 0;
		retried = false;
	};
	virtual void recalibrate();
	virtual bool retry();
	virtual bool calcDirection(bool major);
};

bool Varadhan2008::calcDirection(bool major)
{
	if (!major) return false;
	if (verbose >= 3) mxLog("Varadhan: alpha = %.2f", alpha);

	for (int vx=0; vx < numParam; ++vx) {
		double adj2 = prevAdj1[vx] + prevAdj2[vx];
		double t0 = fc->est[vx] - adj2;
		dir[vx] = (t0 + 2 * alpha * rr[vx] + alpha * alpha * vv[vx]) - fc->est[vx];
	}
	return true;
}

void Varadhan2008::recalibrate()
{
	if (numParam == 0) return;

	memcpy(vv.data(), &prevAdj1[0], sizeof(double) * numParam);
	vv -= rr;

	if (maxAlpha && !retried && alpha > 0) maxAlpha = alpha*2;
	double newAlpha = rr.norm() / vv.norm();
	alpha = newAlpha - 0.5;     // slightly more conservative seems to help
	if (!std::isfinite(alpha) || alpha < 1) alpha = 1;
	if (maxAlpha && alpha > maxAlpha) alpha = maxAlpha;

	retried = false;
}

bool Varadhan2008::retry()
{
	retried = true;
	if (alpha == 1) return false;

	alpha = alpha / 4;
	if (alpha < 1.5) alpha = 1;
	maxAlpha = alpha;
	return true;
}

omxCompute::omxCompute()
{
	varGroup = NULL;
}

void omxCompute::collectResultsHelper(FitContext *fc, std::vector< omxCompute* > &clist,
				      LocalComputeResult *lcr, MxRList *out)
{
	for (std::vector< omxCompute* >::iterator it = clist.begin(); it != clist.end(); ++it) {
		omxCompute *c1 = *it;
		c1->collectResults(fc, lcr, out);
	}
}

void omxCompute::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	MxRList *slots = new MxRList();
        reportResults(fc, slots, out);
	if (slots->size()) {
		lcr->push_back(std::make_pair(computeId, slots));
	} else {
		delete slots;
	}
}

omxCompute::~omxCompute()
{}

void omxCompute::initFromFrontend(omxState *globalState, SEXP rObj)
{
	SEXP slotValue;
	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("id")));
	if (Rf_length(slotValue) != 1) Rf_error("MxCompute has no ID");
	computeId = INTEGER(slotValue)[0];
	}

	varGroup = Global->findVarGroup(computeId);

	if (!varGroup) {
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("freeSet")));
		if (Rf_length(slotValue) == 0) {
			varGroup = Global->findVarGroup(FREEVARGROUP_NONE);
		} else if (strcmp(CHAR(STRING_ELT(slotValue, 0)), ".")==0) {
			varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
		} else {
			Rf_warning("MxCompute ID %d references matrix '%s' in its freeSet "
				"but this matrix contains no free parameters",
				computeId, CHAR(STRING_ELT(slotValue, 0)));
			varGroup = Global->findVarGroup(FREEVARGROUP_NONE);
		}
	}
	if (OMX_DEBUG) {
		mxLog("MxCompute id %d assigned to var group %d", computeId, varGroup->id[0]);
	}
}

void omxCompute::compute(FitContext *fc)
{
	FitContext *narrow = fc;
	if (fc->varGroup != varGroup) narrow = new FitContext(fc, varGroup);
	computeWithVarGroup(narrow);
	if (fc->varGroup != varGroup) narrow->updateParentAndFree();
}

void omxCompute::computeWithVarGroup(FitContext *fc)
{
	ComputeInform origInform = fc->getInform();
	if (OMX_DEBUG) { mxLog("enter %s varGroup %d", name, varGroup->id[0]); }
	fc->setInform(INFORM_UNINITIALIZED);
	computeImpl(fc);
	fc->setInform(std::max(origInform, fc->getInform()));
	if (OMX_DEBUG) { mxLog("exit %s varGroup %d inform %d", name, varGroup->id[0], fc->getInform()); }
	fc->destroyChildren();
	Global->checkpointMessage(fc, fc->est, "%s", name);
}

class ComputeContainer : public omxCompute {
	typedef omxCompute super;
protected:
	std::vector< omxCompute* > clist;
public:
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
};

void ComputeContainer::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);
	collectResultsHelper(fc, clist, lcr, out);
}

class omxComputeSequence : public ComputeContainer {
	typedef ComputeContainer super;
	bool independent;

 public:
	virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
	virtual ~omxComputeSequence();
};

class omxComputeIterate : public ComputeContainer {
	typedef ComputeContainer super;
	int maxIter;
	double maxDuration;
	double tolerance;
	int iterations;
	int verbose;

 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
	virtual ~omxComputeIterate();
};

class omxComputeOnce : public omxCompute {
	typedef omxCompute super;
	std::vector< omxMatrix* > algebras;
	std::vector< omxExpectation* > expectations;
	std::vector< const char* > predict;
	const char *how;
	int verbose;
	bool mac;
	bool starting;
	bool fit;
	bool gradient;
	bool hessian;
	bool ihessian;
	bool infoMat;
	enum ComputeInfoMethod infoMethod;
	bool hgprod;
	bool isBestFit; // for backward compatibility

 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeEM : public omxCompute {
	typedef omxCompute super;
	omxCompute *estep;
	omxCompute *mstep;
	omxMatrix *fit3;   // rename to observedFit
	int EMcycles;
	int maxIter;
	int mstepIter;
	int totalMstepIter;
	double tolerance;
	double semTolerance;
	int verbose;
	Eigen::VectorXd lbound;
	Eigen::VectorXd ubound;
	const char *accelName;
	bool useRamsay;
	bool useVaradhan;
	EMAccel *accel;
	enum EMInfoMethod {
		EMInfoNone,
		EMInfoMengRubinFamily,
		EMInfoOakes
	} information;
	std::vector< omxMatrix * > infoFitFunction;
	enum ComputeInfoMethod infoMethod;
	enum SEMMethod { ClassicSEM, TianSEM, GridSEM, AgileSEM } semMethod;
	double *semMethodData;
	int semMethodLen;
	bool semDebug;
	bool semFixSymmetry;
	bool semForcePD;
	int agileMaxIter;
	SEXP rateMatrix; //debug
	SEXP inputInfoMatrix; //debug
	SEXP outputInfoMatrix; //debug
	SEXP origEigenvalues; //debug
	double noiseTarget;
	double noiseTolerance;
	std::vector<double*> estHistory;
	Eigen::MatrixXd probeOffset;
	std::vector<double> diffWork;
	std::vector<int> paramProbeCount;
	Eigen::VectorXd optimum;
	double bestFit;
 	static const double MIDDLE_START;
	static const double MIDDLE_END;
	size_t maxHistLen;
	int semProbeCount;

	void observedFit(FitContext *fc);
	template <typename T1>
	void accelLineSearch(bool major, FitContext *fc, Eigen::MatrixBase<T1> &preAccel);
	template <typename T>
	bool probeEM(FitContext *fc, int vx, double offset, Eigen::MatrixBase<T> &rijWork);

	template <typename T>
	void recordDiff(FitContext *fc, int v1, Eigen::MatrixBase<T> &rijWork,
			double *stdDiff, bool *mengOK);
	void MengRubinFamily(FitContext *fc);

	void Oakes(FitContext *fc);

 public:
	template <typename T1, typename T2>
	void dEstep(FitContext *fc, Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result);

        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
	virtual ~ComputeEM();
};

const double ComputeEM::MIDDLE_START = 0.105360515657826281366; // -log(.9) constexpr
const double ComputeEM::MIDDLE_END = 0.001000500333583534363566; // -log(.999) constexpr

class ComputeStandardError : public omxCompute {
	typedef omxCompute super;
 public:
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeHessianQuality : public omxCompute {
	typedef omxCompute super;
	int verbose;
 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeReportDeriv : public omxCompute {
	typedef omxCompute super;
 public:
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeReportExpectation : public omxCompute {
	typedef omxCompute super;
 public:
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeBootstrap : public omxCompute {
	typedef omxCompute super;
	
	struct context {
		omxData *data;
		double *origRowWeights;
		std::vector<double> origCumSum;
		std::vector<double> resample;
	};
	std::vector< context > contexts;
	omxCompute *plan;
	int verbose;
	int numReplications;
	int seed;
	//std::vector<double> quantile;
	bool parallel;
	int only;

	int previousNumParam;
	SEXP previousData;
	SEXP rawOutput;
	MxRList onlyWeight;

 public:
	ComputeBootstrap() : plan(0) {};
	virtual ~ComputeBootstrap();
	virtual void initFromFrontend(omxState *, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
	virtual void reportResults(FitContext *fc, MxRList *, MxRList *result);
};

static class omxCompute *newComputeSequence()
{ return new omxComputeSequence(); }

static class omxCompute *newComputeIterate()
{ return new omxComputeIterate(); }

static class omxCompute *newComputeOnce()
{ return new omxComputeOnce(); }

static class omxCompute *newComputeEM()
{ return new ComputeEM(); }

static class omxCompute *newComputeStandardError()
{ return new ComputeStandardError(); }

static class omxCompute *newComputeHessianQuality()
{ return new ComputeHessianQuality(); }

static class omxCompute *newComputeReportDeriv()
{ return new ComputeReportDeriv(); }

static class omxCompute *newComputeReportExpectation()
{ return new ComputeReportExpectation(); }

static class omxCompute *newComputeBootstrap()
{ return new ComputeBootstrap(); }

struct omxComputeTableEntry {
        char name[32];
        omxCompute *(*ctor)();
};

static const struct omxComputeTableEntry omxComputeTable[] = {
        {"MxComputeNumericDeriv", &newComputeNumericDeriv},
        {"MxComputeGradientDescent", &newComputeGradientDescent},
	{"MxComputeSequence", &newComputeSequence },
	{"MxComputeIterate", &newComputeIterate },
	{"MxComputeOnce", &newComputeOnce },
        {"MxComputeNewtonRaphson", &newComputeNewtonRaphson},
        {"MxComputeEM", &newComputeEM },
	{"MxComputeStandardError", &newComputeStandardError},
	{"MxComputeHessianQuality", &newComputeHessianQuality},
	{"MxComputeReportDeriv", &newComputeReportDeriv},
	{"MxComputeConfidenceInterval", &newComputeConfidenceInterval},
	{"MxComputeReportExpectation", &newComputeReportExpectation},
	{"MxComputeTryHard", &newComputeTryHard},
	{"MxComputeNelderMead", &newComputeNelderMead},
	{"MxComputeBootstrap", &newComputeBootstrap},
};

omxCompute *omxNewCompute(omxState* os, const char *type)
{
        omxCompute *got = NULL;

        for (size_t fx=0; fx < OMX_STATIC_ARRAY_SIZE(omxComputeTable); fx++) {
                const struct omxComputeTableEntry *entry = omxComputeTable + fx;
                if(strcmp(type, entry->name) == 0) {
                        got = entry->ctor();
			got->name = entry->name;
                        break;
                }
        }

        if (!got) Rf_error("Compute plan step '%s' is not implemented", type);

        return got;
}

void omxComputeSequence::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("independent")));
	independent = Rf_asLogical(slotValue);
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("steps")));

	for (int cx = 0; cx < Rf_length(slotValue); cx++) {
		SEXP step = VECTOR_ELT(slotValue, cx);
		SEXP s4class;
		const char *s4name;
		{
			ScopedProtect p1(s4class, STRING_ELT(Rf_getAttrib(step, R_ClassSymbol), 0));
			s4name = CHAR(s4class);
		}
		omxCompute *compute = omxNewCompute(globalState, s4name);
		compute->initFromFrontend(globalState, step);
		if (isErrorRaised()) break;
		clist.push_back(compute);
	}

	if (independent) {
		bool fail = false;
		for (int c1 = 1; c1 < (int) clist.size(); ++c1) {
			for (int c2 = 0; c2 < c1; ++c2) {
				if (clist[c1]->varGroup->isDisjoint(clist[c2]->varGroup)) continue;
				omxRaiseErrorf("mxComputeSequence(independent=TRUE) but steps "
					       "%d and %d contain some of the same "
					       "free parameters", 1+c1, 1+c2);
				fail = true;
				break;
			}
			if (fail) break;
		}
	}
}

void omxComputeSequence::computeImpl(FitContext *fc)
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		clist[cx]->compute(fc);
		if (isErrorRaised()) break;
	}
}

omxComputeSequence::~omxComputeSequence()
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		delete clist[cx];
	}
}

void omxComputeIterate::initFromFrontend(omxState *globalState, SEXP rObj)
{
	SEXP slotValue;

	super::initFromFrontend(globalState, rObj);

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];
	}

	{
		ProtectedSEXP RmaxDur(R_do_slot(rObj, Rf_install("maxDuration")));
		maxDuration = Rf_asReal(RmaxDur);
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (std::isfinite(tolerance) && tolerance <= 0) Rf_error("tolerance must be positive");
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("steps")));

	for (int cx = 0; cx < Rf_length(slotValue); cx++) {
		SEXP step = VECTOR_ELT(slotValue, cx);
		SEXP s4class;
		const char *s4name;
		{
			ScopedProtect p1(s4class, STRING_ELT(Rf_getAttrib(step, R_ClassSymbol), 0));
			s4name = CHAR(s4class);
		}
		omxCompute *compute = omxNewCompute(globalState, s4name);
		compute->initFromFrontend(globalState, step);
		if (isErrorRaised()) break;
		clist.push_back(compute);
	}

	{
		ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(slotValue);
	}
	iterations = 0;
}

void omxComputeIterate::computeImpl(FitContext *fc)
{
	double prevFit = 0;
	double mac = tolerance * 10;
	time_t startTime = time(0);
	while (1) {
		++iterations;
		++fc->iterations;
		for (size_t cx=0; cx < clist.size(); ++cx) {
			clist[cx]->compute(fc);
			if (isErrorRaised()) break;
		}
		if (fc->wanted & FF_COMPUTE_MAXABSCHANGE) {
			if (fc->mac < 0) {
				Rf_warning("MAC estimated at %.4f; something is wrong", fc->mac);
				break;
			} else {
				mac = fc->mac;
				if (verbose) mxLog("ComputeIterate: mac %.9g", mac);
			}
		}
		if (fc->wanted & FF_COMPUTE_FIT) {
			if (fc->fit == 0) {
				Rf_warning("Fit estimated at 0; something is wrong");
				break;
			}
			if (prevFit != 0) {
				double change = (prevFit - fc->fit) / fc->fit;
				if (verbose) mxLog("ComputeIterate: fit %.9g rel change %.9g", fc->fit, change);
				mac = fabs(change);
			} else {
				if (verbose) mxLog("ComputeIterate: initial fit %.9g", fc->fit);
			}
			prevFit = fc->fit;
		}
		if (std::isfinite(tolerance)) {
			if (!(fc->wanted & (FF_COMPUTE_MAXABSCHANGE | FF_COMPUTE_FIT))) {
				omxRaiseErrorf("ComputeIterate: neither MAC nor fit available");
			}
			if (mac < tolerance) break;
		}
		if (std::isfinite(maxDuration) && time(0) - startTime > maxDuration) break;
		if (isErrorRaised() || iterations >= maxIter) break;
	}
}

void omxComputeIterate::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	MxRList output;
	output.add("iterations", Rf_ScalarInteger(iterations));
	slots->add("output", output.asR());
}

omxComputeIterate::~omxComputeIterate()
{
	for (size_t cx=0; cx < clist.size(); ++cx) {
		delete clist[cx];
	}
}

void ComputeEM::initFromFrontend(omxState *globalState, SEXP rObj)
{
	SEXP slotValue;
	SEXP s4class;

	super::initFromFrontend(globalState, rObj);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("estep")));
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	estep = omxNewCompute(globalState, CHAR(s4class));
	estep->initFromFrontend(globalState, slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("mstep")));
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	mstep = omxNewCompute(globalState, CHAR(s4class));
	mstep->initFromFrontend(globalState, slotValue);

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("observedFit")));
	fit3 = globalState->algebraList[ INTEGER(slotValue)[0] ];
	if (fit3->fitFunction) {
		omxCompleteFitFunction(fit3);
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];
	if (maxIter < 0) Rf_error("maxIter must be non-negative");
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) Rf_error("tolerance must be positive");
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
	}

	useRamsay = false;
	useVaradhan = false;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("accel")));
	accelName = CHAR(STRING_ELT(slotValue, 0));
	if (strEQ(accelName, "ramsay1975")) {
		useRamsay = true;
	} else if (strEQ(accelName, "varadhan2008")) {
		useVaradhan = true;
	} else if (STRING_ELT(slotValue, 0) == NA_STRING || strEQ(accelName, "")) {
		accelName = "none";
	} else {
		Rf_warning("%s: unknown acceleration method %s ignored", name, accelName);
		accelName = "none";
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("information")));
	const char *infoName = CHAR(STRING_ELT(slotValue, 0));
	information = EMInfoNone;
	if (STRING_ELT(slotValue, 0) == NA_STRING) {
		// ok
	} else if (strEQ(infoName, "mr1991")) {
		information = EMInfoMengRubinFamily;
	} else if (strEQ(infoName, "oakes1999")) {
		information = EMInfoOakes;
	} else {
		Rf_warning("%s: unknown information method %s ignored", name, infoName);
	}

	if (information == EMInfoOakes) {
		infoMethod = INFO_METHOD_HESSIAN;

		SEXP infoArgs, argNames;
		Rf_protect(infoArgs = R_do_slot(rObj, Rf_install("infoArgs")));
		Rf_protect(argNames = Rf_getAttrib(infoArgs, R_NamesSymbol));

		for (int ax=0; ax < Rf_length(infoArgs); ++ax) {
			const char *key = R_CHAR(STRING_ELT(argNames, ax));
			slotValue = VECTOR_ELT(infoArgs, ax);
			if (strEQ(key, "fitfunction")) {
				for (int fx=0; fx < Rf_length(slotValue); ++fx) {
					omxMatrix *ff = globalState->algebraList[INTEGER(slotValue)[fx]];
					if (!ff->fitFunction) Rf_error("infoArgs$fitfunction is %s, not a fitfunction", ff->name());
					infoFitFunction.push_back(ff);
				}
			} else if (strEQ(key, "inputInfo")) {
				infoMethod = stringToInfoMethod(CHAR(slotValue));
			} else {
				mxLog("%s: unknown key %s", name, key);
			}
		}
	}

	if (information == EMInfoMengRubinFamily) {
		// TODO remove SEM prefix everywhere
		infoMethod = INFO_METHOD_HESSIAN;
		semMethod = AgileSEM;
		agileMaxIter = 1;
		semDebug = false;
		semFixSymmetry = true;
		semForcePD = false;
		noiseTarget = exp(-5.2); //constexpr
		noiseTolerance = exp(2.0); //constexpr

		// Meng & Rubin set this parameter in terms of the absolute tolerance
		// instead of the relative tolerance. semTolerance is used to compare
		// parameter estimates, not log-likelihood, so we cannot provide a
		// sane default based on the tolerance.
		semTolerance = 0.01; // sqrt(1e-4)

		SEXP infoArgs, argNames;
		Rf_protect(infoArgs = R_do_slot(rObj, Rf_install("infoArgs")));
		Rf_protect(argNames = Rf_getAttrib(infoArgs, R_NamesSymbol));

		for (int ax=0; ax < Rf_length(infoArgs); ++ax) {
			const char *key = R_CHAR(STRING_ELT(argNames, ax));
			slotValue = VECTOR_ELT(infoArgs, ax);
			if (strEQ(key, "fitfunction")) {
				for (int fx=0; fx < Rf_length(slotValue); ++fx) {
					omxMatrix *ff = globalState->algebraList[INTEGER(slotValue)[fx]];
					if (!ff->fitFunction) Rf_error("infoArgs$fitfunction is %s, not a fitfunction", ff->name());
					infoFitFunction.push_back(ff);
				}
			} else if (strEQ(key, "inputInfo")) {
				infoMethod = stringToInfoMethod(CHAR(slotValue));
			} else if (strEQ(key, "semMethod")) {
				semMethodLen = Rf_length(slotValue);
				if (semMethodLen == 0) {
					semMethod = AgileSEM;
					semMethodData = NULL;
				} else if (Rf_isReal(slotValue)) {
					semMethodData = REAL(slotValue);
					semMethod = GridSEM;
				} else if (Rf_isString(slotValue)) {
					const char *methodName = CHAR(Rf_asChar(slotValue));
					if (strEQ(methodName, "mr")) {
						semMethod = ClassicSEM;
					} else if (strEQ(methodName, "tian")) {
						semMethod = TianSEM;
					} else if (strEQ(methodName, "agile")) {
						semMethod = AgileSEM;
					} else {
						Rf_error("Unknown SEM method '%s'", methodName);
					}
				}
			} else if (strEQ(key, "agileMaxIter")) {
				agileMaxIter = INTEGER(slotValue)[0];
			} else if (strEQ(key, "semDebug")) {
				semDebug = Rf_asLogical(slotValue);
			} else if (strEQ(key, "semFixSymmetry")) {
				semFixSymmetry = Rf_asLogical(slotValue);
			} else if (strEQ(key, "semForcePD")) {
				semForcePD = Rf_asLogical(slotValue);
			} else if (strEQ(key, "semTolerance")) {
				semTolerance = Rf_asReal(slotValue);
				if (semTolerance <= 0) Rf_error("semTolerance must be positive");
			} else if (strEQ(key, "noiseTarget")) {
				noiseTarget = REAL(slotValue)[0];
				if (noiseTarget <= 0) Rf_error("noiseTarget must be positive");
			} else if (strEQ(key, "noiseTolerance")) {
				noiseTolerance = REAL(slotValue)[0];
				if (noiseTolerance < 1) Rf_error("noiseTolerance must be >=1");
			} else {
				mxLog("%s: unknown key %s", name, key);
			}
		}
		if (!semFixSymmetry && semForcePD) {
			Rf_warning("%s: semFixSymmetry must be enabled for semForcePD", name);
			semForcePD = false;
		}
	}

	if (information != EMInfoNone && infoFitFunction.size() == 0) {
		omxRaiseErrorf("%s: at least one fitfunction is required to estimate the information matrix. "
			       "Add something like infoArgs=list(fitfunction='fitfunction')", name);
	}

	inputInfoMatrix = NULL;
	outputInfoMatrix = NULL;
	rateMatrix = NULL;
	origEigenvalues = NULL;
}

template <typename T>
bool ComputeEM::probeEM(FitContext *fc, int vx, double offset, Eigen::MatrixBase<T> &rijWork)
{
	bool failed = false;
	probeOffset(paramProbeCount[vx], vx) = offset;

	Eigen::Map< Eigen::VectorXd > Est(fc->est, (int) fc->varGroup->vars.size());
	Est = optimum;
	fc->est[vx] += offset;
	fc->copyParamToModel();

	if (verbose >= 3) mxLog("ComputeEM: probe %d of %s offset %.6f",
				1+paramProbeCount[vx], fc->varGroup->vars[vx]->name, offset);

	estep->compute(fc);
	fc->wanted &= ~FF_COMPUTE_HESSIAN;  // discard garbage
	int informSave = fc->getInform();  // not sure if we want to hide inform here TODO
	mstep->compute(fc);
	if (fc->getInform() > INFORM_UNCONVERGED_OPTIMUM) {
		if (verbose >= 3) mxLog("ComputeEM: probe failed with code %d", fc->getInform());
		failed = true;
	}
	fc->setInform(informSave);

	rijWork.col(paramProbeCount[vx]) = (Est - optimum) / offset;

	paramProbeCount[vx] += 1;
	++semProbeCount;
	return failed;
}

template <typename T>
void ComputeEM::recordDiff(FitContext *fc, int v1, Eigen::MatrixBase<T> &rijWork,
			   double *stdDiff, bool *mengOK)
{
	const int h1 = paramProbeCount[v1]-2;
	const int h2 = h1+1;

	Eigen::ArrayXd diff = (rijWork.col(h1) - rijWork.col(h2)).array().abs();
	*mengOK = (diff < semTolerance).all();

	double dist = fabs(probeOffset(h1, v1) - probeOffset(h2, v1));
	if (dist < tolerance/4) Rf_error("SEM: invalid probe offset distance %.9f", dist);
	*stdDiff = diff.sum() / (diff.size() * dist);
	diffWork[v1 * maxHistLen + h1] = *stdDiff;
	if (verbose >= 2) mxLog("ComputeEM: (%f,%f) mengOK %d diff %f stdDiff %f",
				probeOffset(h1, v1), probeOffset(h2, v1),
				*mengOK, diff.sum() / diff.size(), *stdDiff);
}

void ComputeEM::observedFit(FitContext *fc)
{
	fc->copyParamToModel();
	ComputeFit("EM", fit3, FF_COMPUTE_FIT, fc);
	if (verbose >= 4) mxLog("ComputeEM[%d]: observed fit = %f", EMcycles, fc->fit);

	if (!(fc->wanted & FF_COMPUTE_FIT)) {
		omxRaiseErrorf("ComputeEM: fit not available");
	}
	if (fc->fit == 0) {
		omxRaiseErrorf("Fit estimated at 0; something is wrong");
	}
}

template <typename T1>
void ComputeEM::accelLineSearch(bool major, FitContext *fc, Eigen::MatrixBase<T1> &preAccel)
{
	Eigen::Map< Eigen::VectorXd > pVec(fc->est, fc->numParam);
	if (!accel->calcDirection(major)) {
		observedFit(fc);
		return;
	}
	if (verbose >= 4) mxPrintMat("accelDir", accel->dir);
	double speed = 1.0;
	int retry = 3;
	while (--retry) {
		pVec = (accel->dir * speed + preAccel).cwiseMax(lbound).cwiseMin(ubound);
		observedFit(fc);
		if (std::isfinite(fc->fit)) return;
		speed *= .3;
		if (verbose >= 3) mxLog("%s: fit NaN; reduce accel speed to %f", name, speed);
	}
	// give up
	pVec = preAccel;
	observedFit(fc);
}

void ComputeEM::computeImpl(FitContext *fc)
{
	const double Scale = fabs(Global->llScale);
	double prevFit = 0;
	double mac = tolerance * 10;
	bool converged = false;
	const int freeVars = int(fc->varGroup->vars.size());
	bool in_middle = false;
	maxHistLen = 0;
	EMcycles = 0;
	semProbeCount = 0;
	lbound.resize(0);

	if (verbose >= 1) mxLog("ComputeEM: Welcome, tolerance=%g accel=%s info=%d",
				tolerance, accelName, information);

	if (useRamsay) accel = new Ramsay1975(fc, verbose, -1.25);
	if (useVaradhan) accel = new Varadhan2008(fc, verbose);

	int mstepInform = INFORM_UNINITIALIZED;
	std::vector<double> prevEst(fc->numParam);
	while (EMcycles < maxIter) {
		++ EMcycles;
		memcpy(&prevEst[0], fc->est, sizeof(double) * fc->numParam);
		if (verbose >= 4) mxLog("ComputeEM[%d]: E-step", EMcycles);
		estep->compute(fc);
		fc->wanted &= ~FF_COMPUTE_DERIV;

		{
			if (verbose >= 4) mxLog("ComputeEM[%d]: M-step", EMcycles);
			FitContext *fc1 = new FitContext(fc, mstep->varGroup);
			int startIter = fc1->iterations;
			mstep->compute(fc1);
			fc1->wanted &= ~FF_COMPUTE_HESSIAN;  // discard garbage
			mstepIter = fc1->iterations - startIter;
			totalMstepIter += mstepIter;
			mstepInform = fc1->getInform();
			fc1->updateParentAndFree();
		}

		if (accel) {
			if (!lbound.size()) {
				// bounds might have changed
				lbound.resize(freeVars);
				ubound.resize(freeVars);
				for(int px = 0; px < int(freeVars); px++) {
					lbound[px] = varGroup->vars[px]->lbound;
					ubound[px] = varGroup->vars[px]->ubound;
				}
				if (verbose >= 3) {
					mxPrintMat("lbound", lbound);
					mxPrintMat("ubound", ubound);
				}
			}

			Eigen::Map< Eigen::VectorXd > pVec(fc->est, fc->numParam);
			Eigen::VectorXd preAccel = pVec;
			accel->recordTrajectory(prevEst);

			if (EMcycles > 3 && (EMcycles + 1) % 3 == 0) {
				accel->recalibrate();
				accelLineSearch(true, fc, preAccel);
				while (prevFit < fc->fit) {
					if (!accel->retry()) break;
					accelLineSearch(true, fc, preAccel);
				}
			} else {
				accelLineSearch(false, fc, preAccel);
			}
		} else {
			observedFit(fc);
		}

		if (!std::isfinite(fc->fit)) {
			omxRaiseErrorf("%s: fit not finite in iteration %d", name, EMcycles);
		}

		double change = 0;
		if (prevFit != 0) {
			if (verbose >= 5) {
				for (int px=0; px < freeVars; ++px) {
					omxFreeVar *fv = fc->varGroup->vars[px];
					mxLog("%d~%s %.4f -> %.4f",
					      px, fv->name, prevEst[px], fc->est[px]);
				}
			}

			change = (prevFit - fc->fit) / fc->fit;
			if (verbose >= 2) mxLog("ComputeEM[%d]: msteps %d fit %.9g rel change %.9g",
						EMcycles, mstepIter, fc->fit, change);
			mac = fabs(change);

			// For Tian, in_middle depends on the absolute (not relative) change in LL!
			const double absMac = fabs(prevFit - fc->fit);
			if (absMac < MIDDLE_START * Scale) in_middle = true;
			if (absMac < MIDDLE_END * Scale) in_middle = false;
		} else {
			if (verbose >= 2) mxLog("ComputeEM: msteps %d initial fit %.9g",
						mstepIter, fc->fit);
		}

		prevFit = fc->fit;
		converged = mac < tolerance;
		++fc->iterations;
		if (isErrorRaised() || converged) break;

		if (semMethod == ClassicSEM || ((semMethod == TianSEM || semMethod == AgileSEM) && in_middle)) {
			double *estCopy = new double[freeVars];
			memcpy(estCopy, fc->est, sizeof(double) * freeVars);
			estHistory.push_back(estCopy);
		}
	}

	int wanted = FF_COMPUTE_FIT | FF_COMPUTE_BESTFIT | FF_COMPUTE_ESTIMATE;
	fc->wanted = wanted;
	fc->setInform(converged? mstepInform : INFORM_ITERATION_LIMIT);
	bestFit = fc->fit;
	if (verbose >= 1) mxLog("ComputeEM: cycles %d/%d total mstep %d fit %f inform %d",
				EMcycles, maxIter, totalMstepIter, bestFit, fc->getInform());

	if (!converged || fc->skippedRows || information == EMInfoNone) return;

	optimum.resize(freeVars);
	memcpy(optimum.data(), fc->est, sizeof(double) * freeVars);

	if (information == EMInfoMengRubinFamily) {
		MengRubinFamily(fc);
	} else if (information == EMInfoOakes) {
		Oakes(fc);
	} else {
		Rf_error("Unknown information method %d", information);
	}

	fc->fit = bestFit;
	memcpy(fc->est, optimum.data(), sizeof(double) * freeVars);
	fc->copyParamToModel();
}

struct estep_jacobian_functional {
	ComputeEM *em;
	FitContext *fc;

	estep_jacobian_functional(ComputeEM *_em, FitContext *_fc) : em(_em), fc(_fc) {};

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result) const {
		em->dEstep(fc, x, result);
	};
};

template <typename T1, typename T2>
void ComputeEM::dEstep(FitContext *fc, Eigen::MatrixBase<T1> &x, Eigen::MatrixBase<T2> &result)
{
	Eigen::Map< Eigen::ArrayXd > Est(fc->est, fc->numParam);
	Est = x;
	fc->copyParamToModel();

	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
	}

	Est = optimum;
	fc->copyParamToModelClean();

	fc->grad = Eigen::VectorXd::Zero(fc->numParam);
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_GRADIENT, fc);
	}
	result = fc->grad;
	reportProgress(fc);
}

void ComputeEM::Oakes(FitContext *fc)
{
	if (verbose >= 1) mxLog("ComputeEM: Oakes1999 method=simple");

	int wanted = fc->wanted;
	const int freeVars = (int) fc->varGroup->vars.size();

	estep->compute(fc);
	fc->wanted &= ~FF_COMPUTE_HESSIAN;  // discard garbage

	fc->grad = Eigen::VectorXd::Zero(fc->numParam);
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_GRADIENT, fc);
	}

	Eigen::VectorXd optimumCopy = optimum;  // will be modified
	Eigen::VectorXd refGrad(freeVars);
	refGrad = fc->grad;
	//mxPrintMat("refGrad", refGrad);

	Eigen::MatrixXd jacobian(freeVars, freeVars);
	estep_jacobian_functional ejf(this, fc);
	fd_jacobian<false>(GradientAlgorithm_Forward, 1, 1e-5, ejf, refGrad, optimumCopy, jacobian);

	fc->infoMethod = infoMethod;
	fc->preInfo();
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_INFO, fc);
	}
	fc->postInfo();

	fc->refreshDenseHess();
	double *hess = fc->getDenseHessUninitialized();
	Eigen::Map< Eigen::MatrixXd > hessMat(hess, freeVars, freeVars);
	hessMat += (jacobian + jacobian.transpose()) * 0.5;  //only need upper triangle TODO

	fc->wanted = wanted | FF_COMPUTE_HESSIAN;
}

void ComputeEM::MengRubinFamily(FitContext *fc)
{
	const int freeVars = (int) fc->varGroup->vars.size();

	if (verbose >= 1) mxLog("ComputeEM: MengRubinFamily tolerance=%f semMethod=%d, semTolerance=%f ideal noise=[%f,%f]",
				tolerance, semMethod, semTolerance,
				noiseTarget/noiseTolerance, noiseTarget*noiseTolerance);

	if (semMethod == AgileSEM) {
		maxHistLen = 2 + agileMaxIter * 2;
	} else if (semMethod == ClassicSEM || semMethod == TianSEM) {
		maxHistLen = estHistory.size();
	} else {
		maxHistLen = semMethodLen;
	}

	probeOffset.resize(maxHistLen, freeVars);
	diffWork.resize(maxHistLen * freeVars);
	paramProbeCount.assign(freeVars, 0);
	Eigen::MatrixXd rij(freeVars, freeVars);

	int semConverged=0;
	for (int v1=0; v1 < freeVars; ++v1) {
		Eigen::MatrixXd rijWork(freeVars, maxHistLen);
		int pick = 0;
		bool paramConverged = false;
		if (semMethod == AgileSEM) {
			double offset1 = .001;
			const double stepSize = offset1 * .01;

			if (probeEM(fc, v1, offset1, rijWork)) break;
			double offset2 = offset1 + stepSize;
			if (probeEM(fc, v1, offset2, rijWork)) break;
			double diff;
			bool mengOK;
			recordDiff(fc, v1, rijWork, &diff, &mengOK);
			double midOffset = (offset1 + offset2) / 2;

			paramConverged = true;
			int iter = 0;
			omxBuffer<double> coefHist(agileMaxIter);
			while (++iter <= agileMaxIter &&
			       !(noiseTarget/noiseTolerance < diff && diff < noiseTarget*noiseTolerance)) {
				coefHist[iter-1] = diff * midOffset * midOffset;
				double coef = 0;
				for (int cx=0; cx < iter; ++cx) coef += coefHist[cx];
				coef /= iter;
				if (verbose >= 4) mxLog("ComputeEM: agile iter[%d] coef=%.6g", iter, coef);
				offset1 = sqrt(coef/noiseTarget);
				if (probeEM(fc, v1, offset1, rijWork)) {
					paramConverged = false;
					break;
				}
				if (iter < agileMaxIter || semDebug) {
					offset2 = offset1 + stepSize;
					if (probeEM(fc, v1, offset2, rijWork)) {
						paramConverged = false;
						break;
					}
					midOffset = (offset1 + offset2) / 2;
					recordDiff(fc, v1, rijWork, &diff, &mengOK);
				}
				pick += 2;
			}
		} else if (semMethod == ClassicSEM || semMethod == TianSEM) {
			if (!estHistory.size()) {
				if (verbose >= 1) mxLog("ComputeEM: no history available;"
							" Classic or Tian SEM require convergence history");
				return;
			}
			for (size_t hx=0; hx < estHistory.size(); ++hx) {
				double popt = optimum[v1];
				double offset1 = estHistory[hx][v1] - popt;
				if (paramProbeCount[v1] && fabs(probeOffset(paramProbeCount[v1]-1, v1) -
							     offset1) < tolerance) continue;
				if (fabs(offset1) < tolerance) continue;
				if (probeEM(fc, v1, offset1, rijWork)) break;
				if (hx == 0) continue;
				pick = hx;
				double diff;
				bool mengOK;
				recordDiff(fc, v1, rijWork, &diff, &mengOK);
				if (mengOK) {
					paramConverged = true;
					break;
				}
			}
		} else {
			for (int hx=0; hx < semMethodLen; ++hx) {
				probeEM(fc, v1, semMethodData[hx], rijWork); // ignore errors
				if (hx == 0) continue;
				double diff;
				bool mengOK;
				recordDiff(fc, v1, rijWork, &diff, &mengOK);
			}
			paramConverged = true;
		}

		const char *pname = fc->varGroup->vars[v1]->name;
		if (paramConverged) {
			++semConverged;
			rij.col(v1) = rijWork.col(pick);
			if (verbose >= 2) mxLog("ComputeEM: %s converged in %d probes",
						pname, paramProbeCount[v1]);
		} else {
			if (verbose >= 2) mxLog("ComputeEM: %s failed to converge after %d probes",
						pname, paramProbeCount[v1]);
			break;
		}
	}

	if (verbose >= 1) {
		if (semConverged == freeVars) {
			mxLog("ComputeEM: %d probes used to estimate Hessian", semProbeCount);
		} else {
			mxLog("ComputeEM: %d probes used for SEM but failed to converge", semProbeCount);
		}
	}
	if (semConverged < freeVars) return;

	if (semDebug) {
		Rf_protect(rateMatrix = Rf_allocMatrix(REALSXP, freeVars, freeVars));
		memcpy(REAL(rateMatrix), rij.data(), sizeof(double) * freeVars * freeVars);
	}

	rij = Eigen::MatrixXd::Identity(freeVars, freeVars) - rij;

	//mxLog("rij symm");
	//pda(rij.data(), freeVars, freeVars);

	estep->compute(fc);
	fc->wanted &= ~FF_COMPUTE_HESSIAN;  // discard garbage

	int wanted = fc->wanted;
	fc->wanted = 0;
	fc->infoMethod = infoMethod;
	fc->preInfo();
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_INFO, fc);
	}
	fc->postInfo();

	Rf_protect(inputInfoMatrix = Rf_allocMatrix(REALSXP, freeVars, freeVars));
	double *hess = REAL(inputInfoMatrix);
	fc->copyDenseHess(hess);

	Matrix rijMat(rij.data(), freeVars, freeVars);
	Matrix hessMat(hess, freeVars, freeVars);
	omxBuffer<double> infoBuf(freeVars * freeVars);
	Matrix infoMat(infoBuf.data(), freeVars, freeVars);

	SymMatrixMultiply('L', 'U', 1, 0, hessMat, rijMat, infoMat);  // result not symmetric!

	if (semDebug) {
		// ihess is always symmetric, this could be asymmetric
		Rf_protect(outputInfoMatrix = Rf_allocMatrix(REALSXP, freeVars, freeVars));
		memcpy(REAL(outputInfoMatrix), infoBuf.data(), sizeof(double) * freeVars * freeVars);
	}

	double *ihess = fc->getDenseIHessUninitialized();
	int singular;
	if (semFixSymmetry) {
		MeanSymmetric(infoMat);
		singular = InvertSymmetricIndef(infoMat, 'U');
		memcpy(ihess, infoBuf.data(), sizeof(double) * freeVars * freeVars);
	} else {
		Matrix ihessMat(ihess, freeVars, freeVars);
		singular = MatrixSolve(infoMat, ihessMat, true);
	}
	if (singular) {
		if (verbose >= 1) mxLog("ComputeEM: SEM Hessian is singular %d", singular);
		return;
	}

	if (semForcePD) {
		double *oev = NULL;
		if (semDebug) {
			Rf_protect(origEigenvalues = Rf_allocVector(REALSXP, freeVars));
			oev = REAL(origEigenvalues);
		}
		Matrix mat(ihess, freeVars, freeVars);
		InplaceForcePosSemiDef(mat, oev, &fc->infoCondNum);
	}

	fc->wanted = wanted | FF_COMPUTE_IHESSIAN;
	//pda(fc->ihess, freeVars, freeVars);
}

void ComputeEM::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);

	std::vector< omxCompute* > clist(1);
	clist[0] = mstep;

	collectResultsHelper(fc, clist, lcr, out);
}

void ComputeEM::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	size_t numFree = fc->varGroup->vars.size();
	if (!numFree) return;

	MxRList out;
	out.add("EMcycles", Rf_ScalarInteger(EMcycles));
	out.add("totalMstep", Rf_ScalarInteger(totalMstepIter));
	out.add("semProbeCount", Rf_ScalarInteger(semProbeCount));
	slots->add("output", out.asR());

	if (semDebug) {
		const int freeVars = (int) fc->varGroup->vars.size();
		MxRList dbg;

		if (probeOffset.size()) {
			SEXP Rpo;
			Rf_protect(Rpo = Rf_allocMatrix(REALSXP, maxHistLen, freeVars));
			memcpy(REAL(Rpo), probeOffset.data(), sizeof(double) * maxHistLen * freeVars);
			dbg.add("probeOffset", Rpo);
		}

		if (diffWork.size()) {
			SEXP Rdiff;
			Rf_protect(Rdiff = Rf_allocMatrix(REALSXP, maxHistLen, freeVars));
			memcpy(REAL(Rdiff), diffWork.data(), sizeof(double) * maxHistLen * freeVars);
			dbg.add("semDiff", Rdiff);
		}

		if (paramProbeCount.size()) {
			SEXP Rphl;
			Rf_protect(Rphl = Rf_allocVector(INTSXP, freeVars));
			memcpy(INTEGER(Rphl), paramProbeCount.data(), sizeof(int) * freeVars);
			dbg.add("paramHistLen", Rphl);
		}

		if (inputInfoMatrix) dbg.add("inputInfo", inputInfoMatrix);
		if (outputInfoMatrix) dbg.add("outputInfo", outputInfoMatrix);
		if (rateMatrix) dbg.add("rateMatrix", rateMatrix);
		if (origEigenvalues) dbg.add("origEigenvalues", origEigenvalues);

		slots->add("debug", dbg.asR());
	}
}

ComputeEM::~ComputeEM()
{
	if (accel) delete accel;

	delete estep;
	delete mstep;

	for (size_t hx=0; hx < estHistory.size(); ++hx) {
		delete [] estHistory[hx];
	}
	estHistory.clear();
}

enum ComputeInfoMethod omxCompute::stringToInfoMethod(const char *iMethod)
{
	enum ComputeInfoMethod infoMethod = INFO_METHOD_DEFAULT; // to avoid gcc warning
	if (strcmp(iMethod, "sandwich")==0) {
		infoMethod = INFO_METHOD_SANDWICH;
	} else if (strcmp(iMethod, "meat")==0) {
		infoMethod = INFO_METHOD_MEAT;
	} else if (strcmp(iMethod, "bread")==0) {
		infoMethod = INFO_METHOD_BREAD;
	} else if (strcmp(iMethod, "hessian")==0) {
		infoMethod = INFO_METHOD_HESSIAN;
	} else {
		Rf_error("Unknown information matrix estimation method '%s'", iMethod);
	}
	return infoMethod;
}

void omxComputeOnce::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("from")));
	for (int wx=0; wx < Rf_length(slotValue); ++wx) {
		if (isErrorRaised()) return;
		int objNum = INTEGER(slotValue)[wx];
		if (objNum >= 0) {
			omxMatrix *algebra = globalState->algebraList[objNum];
			if (algebra->fitFunction) {
				omxCompleteFitFunction(algebra);
			}
			algebras.push_back(algebra);
		} else {
			omxExpectation *expectation = globalState->expectationList[~objNum];
			omxCompleteExpectation(expectation);
			expectations.push_back(expectation);
		}
	}
	if (algebras.size() && expectations.size()) {
		Rf_error("MxComputeOnce cannot evaluate expectations and fitfunctions at the same time");
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("what")));
	int whatLen = Rf_length(slotValue);
	if (algebras.size()) {
		if (whatLen == 0) {
			fit = true;
		}
		for (int wx=0; wx < whatLen; ++wx) {
			SEXP elem;
			Rf_protect(elem = STRING_ELT(slotValue, wx));
			const char *what = CHAR(elem);
			if      (strcmp(what, "maxAbsChange")==0) mac = true;
			else if (strcmp(what, "set-starting-values")==0) starting = true;
			else if (strcmp(what, "fit")         ==0) fit = true;
			else if (strcmp(what, "gradient")    ==0) gradient = true;
			else if (strcmp(what, "hessian")     ==0) hessian = true;
			else if (strcmp(what, "information") ==0) infoMat = true;
			else if (strcmp(what, "ihessian")    ==0) ihessian = true;
			else omxRaiseErrorf("mxComputeOnce: don't know how to compute %s", what);
		}

		if (hessian && infoMat) Rf_error("Cannot compute the Hessian and Fisher Information matrix simultaneously");
	} else {
		for (int wx=0; wx < whatLen; ++wx) {
			SEXP elem;
			ScopedProtect p1(elem, STRING_ELT(slotValue, wx));
			predict.push_back(CHAR(elem));
		}
	}

	{ ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install(".is.bestfit")));
	isBestFit = Rf_asLogical(slotValue);
	}

	bool howConflict = false;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("how")));
	if (Rf_length(slotValue) > 1) {
		omxRaiseErrorf("mxComputeOnce: more than one method specified");
	} else if (Rf_length(slotValue) == 1) {
		SEXP elem;
		Rf_protect(elem = STRING_ELT(slotValue, 0));
		if (algebras.size()) {
			const char *iMethod = CHAR(elem);
			if (infoMat) {
				infoMethod = stringToInfoMethod(iMethod);
				if (infoMethod == INFO_METHOD_MEAT && gradient && whatLen == 2) {
					//OK
				} else if (whatLen > 1) {
					howConflict = true;
				}
			} else {
				omxRaiseErrorf("mxComputeOnce: unknown method %s requested", iMethod);
			}
		} else {
			how = CHAR(elem);
			if (whatLen > 1) howConflict = true;
		}
	}
	if (howConflict) {
		omxRaiseErrorf("mxComputeOnce: when how is specified, you can only compute one thing at a time");
	}

	for (int ax=0; ax < (int) algebras.size(); ++ax) {
		omxFitFunction *ff = algebras[ax]->fitFunction;
		if (gradient && (!ff || !ff->gradientAvailable)) {
			Rf_error("Gradient requested but not available");
		}
		if ((hessian || ihessian || hgprod) && (!ff || !ff->hessianAvailable)) {
			// add a separate flag for hgprod TODO
			Rf_error("Hessian requested but not available");
		}
		// add check for information TODO
	}
}

void omxComputeOnce::computeImpl(FitContext *fc)
{
	if (algebras.size()) {
		int want = 0;
		if (starting) {
			want |= FF_COMPUTE_STARTING;
		}
		if (mac) {
			want |= FF_COMPUTE_MAXABSCHANGE;
			fc->mac = 0;
		}
		if (fit) {
			want |= FF_COMPUTE_FIT;
			if (isBestFit) want |= FF_COMPUTE_BESTFIT;
			fc->fit = 0;
		}
		if (gradient) {
			want |= FF_COMPUTE_GRADIENT;
			fc->grad = Eigen::VectorXd::Zero(fc->numParam);
		}
		if (hessian) {
			want |= FF_COMPUTE_HESSIAN;
			fc->clearHessian();
		}
		if (infoMat) {
			want |= FF_COMPUTE_INFO;
			fc->infoMethod = infoMethod;
			fc->grad = Eigen::VectorXd::Zero(fc->numParam);
			fc->clearHessian();
			fc->preInfo();
		}
		if (ihessian) {
			want |= FF_COMPUTE_IHESSIAN;
			fc->clearHessian();
		}
		if (!want) return;

		for (size_t wx=0; wx < algebras.size(); ++wx) {
			omxMatrix *algebra = algebras[wx];
			if (algebra->fitFunction) {
				omxAlgebraPreeval(algebra, fc);
				ComputeFit("Once", algebra, want, fc);
				if (infoMat) {
					fc->postInfo();
				}
			} else {
				omxMarkDirty(algebra);
				omxRecompute(algebra, fc);
			}
		}
	} else if (expectations.size()) {
		if (predict.size() > 1) Rf_error("Not implemented");
		const char *pr1 = "nothing"; // better to default to 0 ?
		if (predict.size()) pr1 = predict[0];
		for (size_t wx=0; wx < expectations.size(); ++wx) {
			omxExpectation *expectation = expectations[wx];
			omxExpectationCompute(fc, expectation, pr1, how);
		}
	}
}

void omxComputeOnce::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	for (size_t ax=0; ax < algebras.size(); ++ax) {
		omxFitFunction *ff = algebras[ax]->fitFunction;
		if (!ff) continue;
		omxPopulateFitFunction(algebras[ax], out);
	}
}

void ComputeStandardError::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	if (!(fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN))) return;

	fc->allocStderrs();  // at least report NAs

	const size_t numParams = fc->numParam;

	Eigen::VectorXd ihessDiag(fc->ihessDiag());

	const double scale = fabs(Global->llScale);

	// This function calculates the standard errors from the Hessian matrix
	// sqrt(scale * diag(solve(hessian)))

	for(size_t i = 0; i < numParams; i++) {
		double got = ihessDiag[i];
		if (got <= 0) continue;
		fc->stderrs[i] = sqrt(scale * got);
	}
}

void ComputeHessianQuality::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;

	ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(slotValue);
}

void ComputeHessianQuality::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	if (!(fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN))) return;

	// See Luenberger & Ye (2008) Second Order Test (p. 190) and Condition Number (p. 239)

	if (fc->infoDefinite != NA_LOGICAL || !doubleEQ(fc->infoCondNum, NA_REAL)) {
		if (verbose >= 1) {
			mxLog("%s: information matrix already determined to be non positive definite; skipping", name);
		}
		return; // already set elsewhere
	}

	int numParams = int(fc->varGroup->vars.size());

	fc->infoDefinite = false;
	double *hessMem = fc->getDenseHessianish();
	if (!hessMem) {
		if (verbose >= 1) {
			mxLog("%s: information matrix not available; skipping", name);
		}
		return;
	}
	Eigen::Map< Eigen::MatrixXd > hess(hessMem, numParams, numParams);
	hess.triangularView<Eigen::Lower>() = hess.transpose().triangularView<Eigen::Lower>();
	Eigen::LDLT< Eigen::MatrixXd > cholH(hess);
	if (cholH.info() != Eigen::Success || !(cholH.vectorD().array() > 0.0).all()) {
		if (verbose >= 1) {
			mxLog("%s: Cholesky decomposition failed", name);
		}
		return;
	}

	Eigen::MatrixXd ihess(numParams, numParams);
	ihess.setIdentity();
	ihess = cholH.solve(ihess);

	// see LAPACK's dtrcon
	double cn = ihess.colwise().sum().maxCoeff() * hess.colwise().sum().maxCoeff();
	if (!std::isfinite(cn)) {
		if (verbose >= 1) {
			mxLog("%s: result is not finite", name);
		}
		return;
	}

	if (cn < 1.0) cn = 1/cn;

	fc->infoDefinite = true;
	fc->infoCondNum = cn;
}

void ComputeReportDeriv::reportResults(FitContext *fc, MxRList *, MxRList *result)
{
	size_t numFree = fc->numParam;

	if (fc->wanted & FF_COMPUTE_GRADIENT) {
		if (!fc->grad.data()) {
			// oh well
		} else {
			SEXP Rgradient;
			Rgradient = Rf_allocVector(REALSXP, numFree);
			result->add("gradient", Rgradient);
			memcpy(REAL(Rgradient), fc->grad.data(), sizeof(double) * numFree);
		}
	}
	if (fc->wanted & FF_COMPUTE_HESSIAN) {
		SEXP Rhessian;
		Rhessian = Rf_allocMatrix(REALSXP, numFree, numFree);
		result->add("hessian", Rhessian);
		fc->copyDenseHess(REAL(Rhessian));
	}
	if (fc->wanted & FF_COMPUTE_IHESSIAN) {
		SEXP Rihessian;
		Rihessian = Rf_allocMatrix(REALSXP, numFree, numFree);
		result->add("ihessian", Rihessian);
		fc->copyDenseIHess(REAL(Rihessian));
	}
}

void ComputeReportExpectation::reportResults(FitContext *fc, MxRList *, MxRList *result)
{
	std::vector< omxExpectation* > &expectationList = fc->state->expectationList;

	SEXP expectations;
	Rf_protect(expectations = Rf_allocVector(VECSXP, expectationList.size()));

	for(size_t index = 0; index < expectationList.size(); index++) {
		if(OMX_DEBUG) { mxLog("Final Calculation of Expectation %d.", (int) index); }
		omxExpectation *curExpectation = expectationList[index];
		omxExpectationRecompute(fc, curExpectation);
		SEXP rExpect;
		Rf_protect(rExpect = Rf_allocVector(LGLSXP, 1)); // placeholder to attach attributes
		if(OMX_DEBUG) { mxLog("Expectation %d has attribute population.", (int) index); }
		curExpectation->populateAttr(rExpect);
		SET_VECTOR_ELT(expectations, index, rExpect);
	}

	result->add("expectations", expectations);
}

ComputeBootstrap::~ComputeBootstrap()
{
	if (plan) delete plan;
}

void ComputeBootstrap::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	SEXP slotValue;
	SEXP s4class;

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = omxNewCompute(globalState, CHAR(s4class));
	plan->initFromFrontend(globalState, slotValue);

	ProtectedSEXP Rdata(R_do_slot(rObj, Rf_install("data")));
	for (int wx=0; wx < Rf_length(Rdata); ++wx) {
		if (isErrorRaised()) return;
		int objNum = INTEGER(Rdata)[wx];
		context ctx;
		ctx.data = globalState->dataList[objNum];
		int numRows = ctx.data->numRawRows();
		if (!numRows) {
			Rf_error("%s: '%s' cannot have row weights",
				 name, ctx.data->name);
		}
		ctx.origRowWeights = ctx.data->getWeightColumn();
		ctx.origCumSum.resize(numRows);
		ctx.resample.resize(ctx.origCumSum.size());
		if (ctx.origRowWeights) {
			std::partial_sum(ctx.origRowWeights, ctx.origRowWeights + ctx.origCumSum.size(),
					 ctx.origCumSum.begin());
		} else {
			for (int rx=0; rx < numRows; ++rx) ctx.origCumSum[rx] = 1+rx;
		}
		contexts.push_back(ctx);
	}

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose);

	ProtectedSEXP Rrepl(R_do_slot(rObj, Rf_install("replications")));
	numReplications = Rf_asInteger(Rrepl);

	ProtectedSEXP Rparallel(R_do_slot(rObj, Rf_install("parallel")));
	parallel = Rf_asLogical(Rparallel);

	ProtectedSEXP Ronly(R_do_slot(rObj, Rf_install("only")));
	only = Rf_asInteger(Ronly);
	if (only != NA_INTEGER) {
		numReplications = 1;
	}

	previousNumParam = -1;
	previousData = 0;

	ProtectedSEXP Routput(R_do_slot(rObj, Rf_install("output")));
	ProtectedSEXP RoutputNames(Rf_getAttrib(Routput, R_NamesSymbol));
	for (int ax=0; ax < Rf_length(Routput); ++ax) {
		const char *key = R_CHAR(STRING_ELT(RoutputNames, ax));
		SEXP val = VECTOR_ELT(Routput, ax);
		if (strEQ(key, "raw")) {
			previousData = val;
		} else if (strEQ(key, "numParam")) {
			previousNumParam = Rf_asInteger(val);
		}
	}
}

void ComputeBootstrap::computeImpl(FitContext *fc)
{
	if (verbose >= 1) mxLog("%s: %d replications seed=%d parallel=%d",
				name, numReplications, seed, int(parallel));

	int numCols = fc->numParam + 3;
	Rf_protect(rawOutput = Rf_allocVector(VECSXP, numCols));
	SEXP colNames = Rf_allocVector(STRSXP, numCols);
	Rf_setAttrib(rawOutput, R_NamesSymbol, colNames);

	SET_STRING_ELT(colNames, 0, Rf_mkChar("seed"));
	SET_VECTOR_ELT(rawOutput, 0, Rf_allocVector(INTSXP, numReplications));
	SET_STRING_ELT(colNames, 1, Rf_mkChar("fit"));
	SET_VECTOR_ELT(rawOutput, 1, Rf_allocVector(REALSXP, numReplications));
	for (int px=0; px < int(fc->numParam); ++px) {
		SET_STRING_ELT(colNames, 2+px, Rf_mkChar(varGroup->vars[px]->name));
		SET_VECTOR_ELT(rawOutput, 2+px, Rf_allocVector(REALSXP, numReplications));
	}
	SET_STRING_ELT(colNames, 2+fc->numParam, Rf_mkChar("statusCode"));
	SET_VECTOR_ELT(rawOutput, 2+fc->numParam, allocInformVector(numReplications));
	markAsDataFrame(rawOutput, numReplications);

	if (previousData && (previousNumParam != int(fc->numParam) ||
			     Rf_length(previousData) != Rf_length(rawOutput))) {
		if (verbose >= 1) mxLog("%s: discarded mismatching previous data (%d/%d %d/%d)",
					name, previousNumParam, int(fc->numParam),
					Rf_length(previousData), Rf_length(rawOutput));
		previousData = 0;
	}

	for (int repl=0; repl < numReplications; ++repl) {
		INTEGER(VECTOR_ELT(rawOutput, 0))[repl] = NA_INTEGER;
		for (int cx=0; cx <= int(fc->numParam); ++cx) {
			REAL(VECTOR_ELT(rawOutput, 1 + cx))[repl] = NA_REAL;
		}
		INTEGER(VECTOR_ELT(rawOutput, 2 + fc->numParam))[repl] = NA_INTEGER;
	}
	if (only == NA_INTEGER && previousData) {
		int toCopy = std::min(Rf_length(VECTOR_ELT(previousData, 0)),
				      numReplications);
		if (verbose >= 1) mxLog("%s: copying %d rows from previous run", name, toCopy);
		memcpy(INTEGER(VECTOR_ELT(rawOutput, 0)),
		       INTEGER(VECTOR_ELT(previousData, 0)),
		       toCopy * sizeof(int));
		for (int cx=0; cx <= int(fc->numParam); ++cx) {
			memcpy(REAL(VECTOR_ELT(rawOutput, 1+cx)),
			       REAL(VECTOR_ELT(previousData, 1+cx)),
			       toCopy * sizeof(double));
		}
		memcpy(INTEGER(VECTOR_ELT(rawOutput, 2 + fc->numParam)),
		       INTEGER(VECTOR_ELT(previousData, 2 + fc->numParam)),
		       toCopy * sizeof(int));
	}

	// implement parallel TODO

	auto *seedVec = INTEGER(VECTOR_ELT(rawOutput, 0));
	if (only == NA_INTEGER || !previousData) {
		GetRNGstate();
		for (int repl=0; repl < numReplications; ++repl) {
			if (seedVec[repl] != NA_INTEGER) continue;
			int seed1 = unif_rand() * std::numeric_limits<int>::max();
			if (seed1 == NA_INTEGER) seed1 = 0; // maybe impossible
			seedVec[repl] = seed1;
		}
		PutRNGstate();
	} else {
		if (only <= Rf_length(VECTOR_ELT(previousData, 0))) {
			if (verbose >= 1) mxLog("%s: using only=%d", name, only);
			seedVec[0] = INTEGER(VECTOR_ELT(previousData, 0))[only - 1];
		} else {
			Rf_error("%s: only=%d but previous data has just %d replications",
				 name, only, Rf_length(VECTOR_ELT(previousData, 0)));
		}
	}

	Eigen::VectorXd origEst = fc->getEst();

	for (int repl=0; repl < numReplications && !isErrorRaised(); ++repl) {
		std::mt19937 generator(seedVec[repl]);
		if (INTEGER(VECTOR_ELT(rawOutput, 2 + fc->numParam))[repl] != NA_INTEGER) continue;
		if (verbose >= 2) mxLog("%s: replication %d", name, repl);
		for (auto &ctx : contexts) {
			ctx.resample.assign(ctx.origCumSum.size(), 0.0);
			int last = ctx.origCumSum.size() - 1;
			int total = ctx.origCumSum[last];
			std::uniform_int_distribution<int> dist(1, total);
			for (int sx=0; sx < total; ++sx) {
				int pick = dist(generator);
				auto rowPick = std::lower_bound(ctx.origCumSum.begin(), ctx.origCumSum.end(), pick);
				int row = rowPick - ctx.origCumSum.begin();
				ctx.resample[row] += 1.0;
			}
			if (verbose >= 4) {
				EigenStdVectorAdaptor<double> rs(ctx.resample);
				mxPrintMat(ctx.data->name, rs);
			}
			ctx.data->setWeightColumn(ctx.resample.data());
			if (only != NA_INTEGER) {
				onlyWeight.add(ctx.data->name, Rcpp::wrap(ctx.resample));
			}
		}
		fc->state->invalidateCache();
		fc->getEst() = origEst;
		plan->compute(fc);
		if (only == NA_INTEGER) {
			fc->wanted &= ~FF_COMPUTE_DERIV;  // discard garbage
		}
		if (verbose >= 3) {
			auto est = fc->getEst();
			mxPrintMat("est", est);
		}
		REAL(VECTOR_ELT(rawOutput, 1))[repl] = fc->fit;
		for (int px=0; px < int(fc->numParam); ++px) {
			REAL(VECTOR_ELT(rawOutput, 2 + px))[repl] = fc->est[px];
		}
		INTEGER(VECTOR_ELT(rawOutput, 2 + fc->numParam))[repl] = fc->wrapInform();
		reportProgress(fc);
	}

	for (auto &ctx : contexts) {
		ctx.data->setWeightColumn(ctx.origRowWeights);
	}

	if (only == NA_INTEGER) {
		fc->setInform(INFORM_UNINITIALIZED);
		fc->getEst() = origEst;
		fc->copyParamToModel();
	}
}

void ComputeBootstrap::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);
	std::vector< omxCompute* > clist(1);
	clist[0] = plan;
	collectResultsHelper(fc, clist, lcr, out);
}

void ComputeBootstrap::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	MxRList output;
	output.add("numParam", Rcpp::wrap(int(fc->numParam)));
	output.add("raw", rawOutput);
	if (only != NA_INTEGER) {
		output.add("weight", onlyWeight.asR());
	}
	slots->add("output", output.asR());
}
