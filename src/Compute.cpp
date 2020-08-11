/*
 *  Copyright 2013-2019 by the individuals mentioned in the source code history
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

#include <algorithm>
#include <stdarg.h>
#include <limits>
#include <random>
#include <fstream>
#include <forward_list>
#include <map>

#include "glue.h"
#include "Compute.h"
#include "omxState.h"
#include "omxExportBackendState.h"
#include "omxRFitFunction.h"
#include "matrix.h"
#include "omxState.h"
#include "omxData.h"
#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <RcppEigenWrap.h>
#include "finiteDifferences.h"
#include "minicsv.h"
#include "LoadDataAPI.h"
#include "EnableWarnings.h"

void pda(const double *ar, int rows, int cols);

static const char *statusCodeLabels[] = { // see ComputeInform type
		"OK", "OK/green",
		"infeasible linear constraint",
		"infeasible non-linear constraint",
		"iteration limit/blue",
		"not convex/red",
		"nonzero gradient/red",
		"bad deriv",
		"?",
		"internal error",
		"infeasible start"
};

SEXP allocInformVector(int size)
{
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
			mxThrow("HessianBlock var mapping is not 1-to-1");
		}
		int prev = hb->vars[0];
		for (int vx=1; vx < int(hb->vars.size()); vx++) {
			if (prev > hb->vars[vx]) {
				mxThrow("hb->vars must be sorted");
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

	int numFree = calcNumFree();
	hess.resize(numFree, numFree);

	hess.triangularView<Eigen::Upper>().setZero();

	// if allBlocks.size() == 0 then what? TODO
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
	int np = hess.rows();
	for (int v1=0; v1 < np; ++v1) {
		for (int v2=0; v2 <= v1; ++v2) {
			double coef = hess.selfadjointView<Eigen::Upper>()(v2,v1);
			if (v1==v2) {
				dest[v1 * np + v2] = coef;
			} else {
				dest[v1 * np + v2] = coef;
				dest[v2 * np + v1] = coef;
			}
		}
	}
}

double *FitContext::getDenseHessUninitialized()
{
	int numFree = calcNumFree();
	hess.resize(numFree, numFree);

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
	ThinMatrix ihessMat(ihess.data(), ihess.rows(), ihess.cols());

	if (!InvertSymmetricPosDef(ihessMat, 'U')) return;

	int numParams = hess.rows();
	allFiniteHelper afh;
	hess.visit(afh);
	if (!afh.finite) {
		ihess = Eigen::MatrixXd::Zero(numParams, numParams);
		return;
	}

	Eigen::Map< Eigen::MatrixXd > EihessMat(ihess.data(), ihess.rows(), ihess.cols());
	ForceInvertSymmetricPosDef(EihessMat);
}

void FitContext::refreshDenseIHess()
{
	if (haveDenseIHess) return;

	refreshDenseHess();

	ihess = hess;
	ThinMatrix ihessMat(ihess.data(), ihess.rows(), ihess.cols());
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
	ProtectAutoBalanceDoodad mpi;

	SEXP matrixDims;
	Rf_protect(matrixDims = Rf_getAttrib(Rmat, R_DimSymbol));
	int *dimList = INTEGER(matrixDims);
	int rows = dimList[0];
	int cols = dimList[1];
	if (rows != cols) mxThrow("Must be square");

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
	if (soleymani2013(mat, imat)) mxThrow("Invert failed");

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
	if (bad > .0001) mxThrow("Hess: dense sparse mismatch %f", bad);
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
			if (!hb) mxThrow("Attempting to invert Hessian, "
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
			mxThrow("IHess: dense sparse mismatch %f", bad);
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
	for (int px=0; px < int(haveGrad.size()); ++px) {
		if (haveGrad[px]) continue;
		mxLog("FitContext::ihessGradProd grad[%d/%s] missing",
		      px, varGroup->vars[px]->name);
	}
	if (refreshSparseIHess()) {
		return sparseIHess.selfadjointView<Eigen::Upper>() * gradZ;
	} else {
		refreshDenseHess();
		InvertSymmetricNR(hess, ihess);
		return ihess.selfadjointView<Eigen::Upper>() * gradZ;
	}
}

Eigen::VectorXd FitContext::ihessDiag()
{
	refreshDenseIHess();
	return ihess.diagonal();
}

double *FitContext::getDenseIHessUninitialized()
{
	int numFree = calcNumFree();
	ihess.resize(numFree, numFree);

	// Assume the caller is going to fill it out
	haveDenseIHess = true;
	haveDenseHess = false;
	return ihess.data();
}

void FitContext::copyDenseIHess(double *dest)
{
	refreshDenseIHess();

	int np = ihess.rows();
	for (int v1=0; v1 < np; ++v1) {
		for (int v2=0; v2 <= v1; ++v2) {
			double coef = ihess.selfadjointView<Eigen::Upper>()(v2,v1);
			if (v1==v2) {
				dest[v1 * np + v2] = coef;
			} else {
				dest[v1 * np + v2] = coef;
				dest[v2 * np + v1] = coef;
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

int FitContext::getDenseHessianishSize()
{
	if (haveDenseHess) return hess.rows();
	if (haveDenseIHess) return ihess.rows();
	return 0;
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
	previousReportFit = nan("uninit");
	fitUnits = parent? parent->fitUnits : FIT_UNITS_UNINITIALIZED;
	skippedRows = 0;
	est = new double[numParam];
	infoDefinite = NA_LOGICAL;
	infoCondNum = NA_REAL;
	infoA = NULL;
	infoB = NULL;
	inform = INFORM_UNINITIALIZED;
	iterations = 0;
	ciobj = 0;
	openmpUser = false;
	ordinalRelativeError = 0;
	computeCount = 0;

	hess.resize(numParam, numParam);  // TODO why needed?
	ihess.resize(numParam, numParam);  // TODO why needed?

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

void FitContext::calcStderrs()
{
	int numFree = calcNumFree();
	stderrs.resize(numFree);
	if (vcov.rows() != numFree || vcov.cols() != numFree) {
		mxThrow("FitContext::calcStderrs vcov size wrong %d vs %d",
          vcov.rows(), numFree);
	}
	const double scale = fabs(Global->llScale);

	if(constraintJacobian.rows()){
		Eigen::MatrixXd hesstmp(numFree, numFree);
		if(fitUnits == FIT_UNITS_MINUS2LL){
			copyDenseHess(hesstmp.data());
			hesstmp = hesstmp/scale;
		}
		else{ //WLS case--'hesstmp' is actually not a Hessian matrix, but the inverse of the WLS sampling-covariance matrix.
			Eigen::LLT< Eigen::MatrixXd > cholVcov;
			cholVcov.compute(vcov.selfadjointView<Eigen::Lower>());
			if(cholVcov.info() != Eigen::Success){
				Eigen::FullPivLU< Eigen::MatrixXd > luVcov(vcov.selfadjointView<Eigen::Lower>());
				if(luVcov.isInvertible()){
					hesstmp = luVcov.inverse();
				}
				else{ //This is supposed to never happen with WLS...
					Rf_warning(
						"constraint-adjusted standard errors could not be calculated because the sampling covariance matrix was uninvertible");
					return;
				}
			}
			else{
				hesstmp = cholVcov.solve(Eigen::MatrixXd::Identity( vcov.rows(), vcov.cols() ));
			}
		}
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qrj(constraintJacobian.transpose());
		Eigen::MatrixXd Q = qrj.householderQ();
		Eigen::MatrixXd U = Q.block(0, qrj.rank(), Q.rows(), Q.cols()-qrj.rank());
		if(U.rows()==0 || U.cols()==0){
			Rf_warning(
				"standard errors could not be calculated because no basis could be found for the nullspace of the constraint Jacobian");
			return;
		}
		if(OMX_DEBUG){mxPrintMat("basis",U);}
		Eigen::MatrixXd centr = U.transpose() * hesstmp * U;
		//centr should be symmetric, will almost always be invertible, may sometimes not be PD:
		Eigen::LLT< Eigen::MatrixXd > cholCentr;
		//Ff center is PD, we'd rather calculate its inverse via its Cholesky factorization:
		cholCentr.compute(centr.selfadjointView<Eigen::Lower>());
		if(cholCentr.info() != Eigen::Success){ //<--Will be true if centr is not PD.
			Eigen::FullPivLU< Eigen::MatrixXd > luCentr(centr.selfadjointView<Eigen::Lower>());
			if(luCentr.isInvertible()){
				centr = luCentr.inverse();
			}
			else{
				Rf_warning(
					"constraint-adjusted standard errors could not be calculated because the coefficient matrix of the quadratic form was uninvertible");
				return;
			}
		}
		else{
			centr = cholCentr.solve(Eigen::MatrixXd::Identity( centr.rows(), centr.cols() ));
		}
		vcov = U * centr * U.transpose();
	}

	for(int i = 0; i < numFree; i++) {
		double got = vcov(i,i);
		if (got <= 0) {
			stderrs[i] = NA_REAL;
			continue;
		}
		stderrs[i] = sqrt(got);
	}
}

FitContext::FitContext(omxState *_state)
{
	parent = NULL;
	varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
	init();
	profiledOut.assign(numParam, false);

	auto &startingValues = Global->startingValues;
	state = _state;
	if (numParam) {
		if (startingValues.size() != numParam) {
			mxThrow("Got %d starting values for %d parameters",
				int(startingValues.size()), int(numParam));
		}
		memcpy(est, startingValues.data(), sizeof(double) * numParam);
	}
	//equality.resize(state->numEqC);
	//inequality.resize(state->numIneqC);
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
	if (d1 != dvars) mxThrow("Parent free parameter group (id=%d) is not a superset of %d",
			       src->id[0], dest->id[0]);


	wanted = parent->wanted;
	infoDefinite = parent->infoDefinite;
	infoCondNum = parent->infoCondNum;
	iterations = parent->iterations;
	ciobj = parent->ciobj;
	//equality.resize(state->numEqC);
	//inequality.resize(state->numIneqC);
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
	parent->recordOrdinalRelativeError(getOrdinalRelativeError());

	// rewrite using mapToParent TODO

	if (svars > 0) {
		for (size_t d1=0, s1 = 0; d1 < dest->vars.size(); ++d1) {
			if (dest->vars[d1] != src->vars[s1]) continue;
			parent->est[d1] = est[s1];
			if (++s1 == svars) break;
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
			buf += string_snprintf("%.16f", est[vx]);
			if (vx < count - 1) buf += ", ";
		}
		buf += ")\n";
	}
	if (what & FF_COMPUTE_GRADIENT) {
		buf += string_snprintf("gradient %d: c(", (int) count);
		for (size_t vx=0; vx < count; ++vx) {
			if (haveGrad[vx]) {
				buf += string_snprintf("%.5f", gradZ[vx]);
			} else {
				buf += '-';
			}
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

std::string FitContext::asProgressReport()
{
	std::string str;
	if (!std::isfinite(previousReportFit) || !std::isfinite(fit)) {
		str = string_snprintf("evaluations %d fit %.6g", getGlobalComputeCount(), fit);
	} else {
		str = string_snprintf("evaluations %d fit %.6g change %.4g",
													getGlobalComputeCount(), fit, fit - previousReportFit);
	}
	previousReportFit = fit;
	return str;
}

void FitContext::resetIterationError()
{
	IterationError.clear();
}

void FitContext::resetOrdinalRelativeError()
{
	if (childList.size()) {
		for (size_t cx=0; cx < childList.size(); ++cx) {
			childList[cx]->resetOrdinalRelativeError();
		}
	}
	if (OMX_DEBUG && ordinalRelativeError > .01) {
		mxLog("reset ordinalRelativeError %g back to zero", ordinalRelativeError);
	}
	ordinalRelativeError = 0;
}

void FitContext::recordIterationError(const char* msg, ...)
{
	// Can avoid overhead of setting error if one is already set TODO
	va_list ap;
	va_start(ap, msg);
	string_vsnprintf(msg, ap, IterationError);
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

	ProtectedSEXP estimate(Rf_allocVector(REALSXP, n));
	double *est = REAL(estimate);
	for(int i = 0; i < n ; i++) {
		est[i] = x[i];
	}

	ProtectedSEXP theCall(Rf_allocVector(LANGSXP, 4));

		// imxUpdateModelValues does not handle parameters with equality
		// constraints. This is a bug.
	SETCAR(theCall, Rf_install("imxUpdateModelValues"));
	SETCADR(theCall, rFitFunction->model);
	SETCADDR(theCall, rFitFunction->flatModel);
	SETCADDDR(theCall, estimate);

	rFitFunction->model = Rf_eval(theCall, R_GlobalEnv);

	Rf_setAttrib(rFitFunction->rObj, Rf_install("model"), rFitFunction->model);

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

// NOTE: All non-linear constraints are applied regardless of free
// variable group.
void FitContext::solEqBFun(bool wantAJ, int verbose) //<--"want analytic Jacobian"
{
	const int eq_n = (int) equality.size();

	if (!eq_n) return;

	/*Note that this needs to happen even if no equality constraints have analytic Jacobians, because
	 analyticEqJacTmp is copied to the Jacobian matrix the elements of which are populated by code in
	 finiteDifferences.h, which knows to numerically populate an element if it's NA:*/
	analyticEqJacTmp.setConstant(NA_REAL);

	int cur=0, j=0, c=0, roffset=0;
	for(j = 0; j < int(state->conListX.size()); j++) {
		omxConstraint &con = *state->conListX[j];
		if (con.opCode != omxConstraint::EQUALITY) continue;

		con.refreshAndGrab(this, &equality(cur));
		if(wantAJ && isUsingAnalyticJacobian() && con.jacobian != NULL){
			omxRecompute(con.jacobian, this);
			for(c=0; c<con.jacobian->cols; c++){
				if(con.jacMap[c]<0){continue;}
				for(roffset=0; roffset<con.size; roffset++){
					analyticEqJacTmp(cur+roffset,con.jacMap[c]) = con.jacobian->data[c * con.size + roffset];
				}
			}
		}
		cur += con.size;
	}

	if (verbose >= 3) {
		mxPrintMat("equality", equality);
	}
}

void FitContext::myineqFun(bool wantAJ, int verbose, int ineqType, bool CSOLNP_HACK)
{
	const int ineq_n = (int) inequality.size();

	if (!ineq_n) return;

	analyticIneqJacTmp.setConstant(NA_REAL);

	int cur=0, j=0, c=0, roffset=0;
	for (j=0; j < int(state->conListX.size()); j++) {
		omxConstraint &con = *state->conListX[j];
		if (con.opCode == omxConstraint::EQUALITY) continue;

		con.refreshAndGrab(this, (omxConstraint::Type) ineqType, &inequality(cur));
		if(wantAJ && isUsingAnalyticJacobian() && con.jacobian != NULL){
			omxRecompute(con.jacobian, this);
			for(c=0; c<con.jacobian->cols; c++){
				if(con.jacMap[c]<0){continue;}
				for(roffset=0; roffset<con.size; roffset++){
					analyticIneqJacTmp(cur+roffset,con.jacMap[c]) = con.jacobian->data[c * con.size + roffset];
				}
			}
		}
		cur += con.size;
	}

	if (CSOLNP_HACK) {
		// CSOLNP doesn't know that inequality constraints can be inactive (by design, since it's an interior-point algorithm)
	} else {
		//SLSQP seems to require inactive inequality constraint functions to be held constant at zero:
		inequality = inequality.array().max(0.0);
		if(wantAJ && isUsingAnalyticJacobian()){
			for(int i=0; i<analyticIneqJacTmp.rows(); i++){
				/*The Jacobians of each inactive constraint are set to zero here;
				 as their elements will be zero rather than NaN, the code in finiteDifferences.h will leave them alone:*/
				if(!inequality[i]){analyticIneqJacTmp.row(i).setZero();}
			}
		}
	}

	if (verbose >= 3) {
		mxPrintMat("inequality", inequality);
	}
}

//Optimizers care about separating equality and inequality constraints, but the ComputeNumericDeriv step doesn't:
void FitContext::allConstraintsF(bool wantAJ, int verbose, int ineqType, bool CSOLNP_HACK, bool maskInactive){
	int c_n = state->numEqC + state->numIneqC;
	if(!c_n){return;}
	std::vector<bool> is_inactive_ineq(c_n);

	constraintJacobian.setConstant(NA_REAL);

	int cur=0;
	for (int j=0; j < int(state->conListX.size()); j++) {
		omxConstraint &con = *state->conListX[j];
		if (con.opCode == omxConstraint::EQUALITY) {
			con.refreshAndGrab(this, &constraintFunVals(cur));
			for(int i=0; i < con.size; i++){
				is_inactive_ineq[cur+i] = false;
			}
		} else{
			con.refreshAndGrab(this, (omxConstraint::Type) ineqType, &constraintFunVals(cur));
			for(int i=0; i < con.size; i++){
				if(constraintFunVals(cur+i) < 0 && maskInactive){
					constraintFunVals(cur+i) = 0;
					is_inactive_ineq[cur+i] = true;
				} else{
					is_inactive_ineq[cur+i] = false;
				}
			}
		}
		if(wantAJ && isUsingAnalyticJacobian() && con.jacobian != NULL){
			omxRecompute(con.jacobian, this);
			for(int c=0; c<con.jacobian->cols; c++){
				if(con.jacMap[c]<0){continue;}
				for(int roffset=0; roffset<con.size; roffset++){
					constraintJacobian(cur+roffset,con.jacMap[c]) = con.jacobian->data[c * con.size + roffset];
				}
			}
		}
		cur += con.size;
	}

	if (CSOLNP_HACK) {

	} else {
		if(wantAJ && isUsingAnalyticJacobian() && maskInactive){
			for(int i=0; i<constraintJacobian.rows(); i++){
				/*The Jacobians of each inactive constraint are set to zero here;
				 as their elements will be zero rather than NaN, the code in finiteDifferences.h will leave them alone:*/
				if(is_inactive_ineq[i]){constraintJacobian.row(i).setZero();}
			}
		}
	}

	if (verbose >= 3) {
		mxPrintMat("constraint Jacobian", constraintJacobian);
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
		mxThrow("Attempt to take %d but not available", want);
	}

	double *ret = NULL;
	switch(want) {
	case FF_COMPUTE_ESTIMATE:
		ret = est;
		est = NULL;
		break;
	default:
		mxThrow("Taking of %d is not implemented", want);
	}
	if (!ret) mxThrow("Attempt to take %d, already taken", want);
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
		OMXZERO(infoA, npsq);
		OMXZERO(infoB, npsq);
		break;
	case INFO_METHOD_BREAD:
		OMXZERO(infoA, npsq);
		break;
	case INFO_METHOD_HESSIAN:
		clearHessian();
		break;
	default:
		mxThrow("Unknown information matrix estimation method %d", infoMethod);
	}
}

void FitContext::postInfo()
{
	switch (infoMethod) {
	case INFO_METHOD_SANDWICH:{
		// move into FCDeriv TODO
		std::vector<double> work(numParam * numParam);
		ThinMatrix amat(infoA, numParam, numParam);
		InvertSymmetricIndef(amat, 'U');
		_fixSymmetry("InfoB", infoB, numParam, false);
		ThinMatrix bmat(infoB, numParam, numParam);
		ThinMatrix wmat(work.data(), numParam, numParam);
		ThinMatrix hmat(getDenseIHessUninitialized(), numParam, numParam);
		SymMatrixMultiply('L', amat, bmat, wmat);
		SymMatrixMultiply('R', amat, wmat, hmat);
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
		mxThrow("Unknown information matrix estimation method %d", infoMethod);
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
		FitContext *kid = new FitContext(this, varGroup);
		kid->state = new omxState(state);
		kid->state->initialRecalc(kid);
		omxAlgebraPreeval(alg, kid);
		childList.push_back(kid);
	}

	if (OMX_DEBUG) mxLog("FitContext::createChildren: done creating %d omxState", Global->numThreads);
}

void FitContext::destroyChildren()
{
	if (0 == childList.size()) return;
	IterationError = getIterationError();
	for (int cx=0; cx < int(childList.size()); ++cx) {
		recordOrdinalRelativeError(childList[cx]->getOrdinalRelativeError());
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
	if (infoA) delete [] infoA;
	if (infoB) delete [] infoB;
}

omxFitFunction *FitContext::RFitFunction = NULL;

void FitContext::setRFitFunction(omxFitFunction *rff)
{
	if (rff) {
		Global->numThreads = 1;
		if (RFitFunction) {
			mxThrow("You can only create 1 MxRFitFunction per independent model");
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
	}

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
	ProtectedSEXP Rid(R_do_slot(rObj, Rf_install("id")));
	if (Rf_length(Rid) != 1) mxThrow("MxCompute has no ID");
	computeId = INTEGER(Rid)[0];

	ProtectedSEXP Rpersist(R_do_slot(rObj, Rf_install(".persist")));
	dotPersist = Rf_asLogical(Rpersist);

	varGroup = Global->findVarGroup(computeId);

	if (!varGroup) {
		ProtectedSEXP Rfreeset(R_do_slot(rObj, Rf_install("freeSet")));
		if (Rf_length(Rfreeset) == 0) {
			varGroup = Global->findVarGroup(FREEVARGROUP_NONE);
		} else if (strcmp(CHAR(STRING_ELT(Rfreeset, 0)), ".")==0) {
			varGroup = Global->findVarGroup(FREEVARGROUP_ALL);
		} else {
			Rf_warning("MxCompute ID %d references matrix '%s' in its freeSet "
				"but this matrix contains no free parameters",
				computeId, CHAR(STRING_ELT(Rfreeset, 0)));
			varGroup = Global->findVarGroup(FREEVARGROUP_NONE);
		}
	}
	if (OMX_DEBUG) {
		mxLog("MxCompute id %d assigned to var group %d", computeId, varGroup->id[0]);
	}
}

void omxCompute::complainNoFreeParam()
{
	if (Global->ComputePersist) {
		omxRaiseErrorf("%s: model has no free parameters; "
			       "You may want to reset your model's "
			       "compute plan with model$compute <- mxComputeDefault() "
			       "and try again", name);
	} else {
		// should never happen, see MxRun.R
		omxRaiseErrorf("%s: model has no free parameters", name);
	}
}

class EnterVarGroup {
	bool narrowed;
 public:
	FitContext *narrow;
	EnterVarGroup(FitContext *fc, FreeVarGroup *varGroup) {
		narrow = fc;
		narrowed = fc->varGroup != varGroup;
		if (narrowed) narrow = new FitContext(fc, varGroup);
	}
	~EnterVarGroup() {
		if (narrowed) narrow->updateParentAndFree();
	}
};

void omxCompute::compute(FitContext *fc)
{
	EnterVarGroup evg(fc, varGroup);
	computeWithVarGroup(evg.narrow);
}

struct LeaveComputeWithVarGroup {
	FitContext *fc;
	bool toResetInform;
	ComputeInform origInform;
	const char *name;

	LeaveComputeWithVarGroup(FitContext *_fc, struct omxCompute *compute) : fc(_fc), name(compute->name) {
		origInform = fc->getInform();
		toResetInform = compute->accumulateInform();
		if (toResetInform) fc->setInform(INFORM_UNINITIALIZED);
		if (Global->debugProtectStack) {
			mxLog("enter %s: protect depth %d", name, Global->mpi->getDepth());
		}
	};
	~LeaveComputeWithVarGroup() {
		fc->destroyChildren();
		if (toResetInform) fc->setInform(std::max(origInform, fc->getInform()));
		Global->checkpointMessage(fc, fc->est, "%s", name);
		if (Global->debugProtectStack) {
			mxLog("exit %s: protect depth %d", name, Global->mpi->getDepth());
		}

	}
};

void omxCompute::computeWithVarGroup(FitContext *fc)
{
	LeaveComputeWithVarGroup lcwvg(fc, this);
	computeImpl(fc);
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

class ComputeLoop : public ComputeContainer {
	typedef ComputeContainer super;
	int verbose;
	int indicesLength;
	int *indices;
	int maxIter;
	double maxDuration;
	int iterations;
	int startFrom;

 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
	virtual ~ComputeLoop();
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

struct ComputeSetOriginalStarts : public omxCompute {
	virtual bool accumulateInform() { return false; };
        virtual void computeImpl(FitContext *fc);
};

class ComputeStandardError : public omxCompute {
	typedef omxCompute super;
	omxMatrix *fitMat;
	std::vector<omxExpectation *> exList;
	std::vector<int> numStats;
	bool wlsStats;
	double x2;
	int df;
	double x2m, x2mv;
	double madj, mvadj, dstar;

	struct visitEx {
		ComputeStandardError &top;
		visitEx(ComputeStandardError *cse) : top(*cse) {};
		void operator()(omxMatrix *mat) const
		{
			if (!mat->fitFunction) {
				omxRaiseErrorf("%s: Cannot compute SEs when '%s' included in fit",
					 top.name, mat->name());
				return;
			}
			omxExpectation *e1 = mat->fitFunction->expectation;
			if (!e1) return;
			if (!e1->data) {
				omxRaiseErrorf("%s: expectation '%s' does not have data",
					 top.name, e1->name);
				return;
			}
			omxData *d1 = e1->data;
			d1->visitObsStats([this, d1](obsSummaryStats &o1) {
					if (o1.fullWeight) return;
					mxThrow("%s: terribly sorry, master, but '%s' does not "
						"include the full weight matrix hence "
						"standard errors cannot be computed",
						top.name, d1->name);
				});
			top.exList.push_back(e1);
		}
	};

 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

struct ParJacobianSense {
	FitContext *fc;
	std::vector<omxExpectation *> *exList;
	std::vector<omxMatrix *> *alList;
	int numOf;
	std::vector<int> numStats;
	int maxNumStats;
	int totalNumStats;
	int defvar_row;
	Eigen::VectorXd ref;
	Eigen::MatrixXd result;

	// default defvar_row to 1 or 0? TODO
	ParJacobianSense() : exList(0), alList(0), defvar_row(1) {};

	void attach(std::vector<omxExpectation *> *_exList, std::vector<omxMatrix *> *_alList) {
		if (_exList && _alList) mxThrow("_exList && _alList");
		exList = _exList;
		alList = _alList;
		numOf = exList? exList->size() : alList->size();
		numStats.reserve(numOf);
		maxNumStats = 0;
		totalNumStats = 0;
		for (int ex=0; ex < numOf; ++ex) {
			int nss = exList? (*exList)[ex]->numSummaryStats() : (*alList)[ex]->size();
			numStats.push_back(nss);
			totalNumStats += nss;
			maxNumStats = std::max(nss, maxNumStats);
		}
	};

	void measureRef(FitContext *_fc) {
		using Eigen::Map;
		using Eigen::VectorXd;
		fc = _fc;
		int numFree = fc->calcNumFree();
		Map< VectorXd > curEst(fc->est, numFree);
		result.resize(totalNumStats, numFree);
		ref.resize(totalNumStats);
		(*this)(curEst, ref);
	}

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &, Eigen::MatrixBase<T2> &result1) const {
		fc->copyParamToModel();
		Eigen::VectorXd tmp(maxNumStats);
		for (int ex=0, offset=0; ex < numOf; offset += numStats[ex++]) {
			if (exList) {
				(*exList)[ex]->asVector(fc, defvar_row, tmp);
				result1.block(offset, 0, numStats[ex], 1) =
					tmp.segment(0, numStats[ex]);
			} else {
				omxMatrix *mat = (*alList)[ex];
				omxRecompute(mat, fc);
				EigenVectorAdaptor vec(mat);
				if (numStats[ex] != vec.size()) {
					mxThrow("Algebra '%s' changed size during Jacobian", mat->name());
				}
				result1.block(offset, 0, numStats[ex], 1) = vec;
			}
		}
	}
};

struct ParJacobianSense1 {
	FitContext *fc;
	omxMatrix *alg;
	Eigen::MatrixXd ref;
	Eigen::MatrixXd result;

	ParJacobianSense1() : alg(0) {};

	void attach(omxMatrix * _alg) {
		if (_alg) mxThrow("_alg");
		alg = _alg;
	};

	void measureRef(FitContext *_fc) {
		using Eigen::Map;
		using Eigen::VectorXd;
		fc = _fc;
		int numFree = fc->calcNumFree();
		Map< VectorXd > curEst(fc->est, numFree);
		result.resize(alg->rows, alg->cols);
		ref.resize(alg->rows, alg->cols);
		(*this)(curEst, ref);
	}

	template <typename T1, typename T2>
	void operator()(Eigen::MatrixBase<T1> &, Eigen::MatrixBase<T2> &result1) const {
		fc->copyParamToModel();
		omxRecompute(alg, fc);
		EigenMatrixAdaptor Ealg(alg);
		result1 = Ealg;
	}
};

// usage:
//
// ParJacobianSense1 sense;
// sense.attach(myAlgebra);
// sense.measureRef(fc);
// fd_jacobian1<false>(GradientAlgorithm_Forward, 2, 1e-4, sense, sense.ref, curEst, px, sense.result);

class ComputeJacobian : public omxCompute {
	typedef omxCompute super;
	std::vector<omxExpectation *> exList;
	std::vector<omxMatrix *> alList;
	omxData *data;
	ParJacobianSense sense;

 public:
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
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

class ComputeCheckpoint : public omxCompute {
	typedef omxCompute super;

	struct snap {
		int evaluations;
		int iterations;
		time_t timestamp;
		std::vector<int> computeLoopIndex;
		Eigen::VectorXd est;
		double fit;
		FitStatisticUnits fitUnits;
		int inform;
		Eigen::VectorXd stderrs;
		Eigen::VectorXd gradient;
		Eigen::VectorXd vcov;
		Eigen::VectorXd algebraEnt;
		std::vector< std::string > extra;
	};

	const char *path;
	std::ofstream ofs;
	bool toReturn;
	std::vector<omxMatrix*> algebras;
	int numAlgebraEnt;
	bool wroteHeader;
	std::vector<std::string> colnames;
	std::forward_list<snap> snaps;
	int numSnaps;
	int numParam;
	bool inclPar, inclLoop, inclFit, inclCounters, inclStatus;
	bool inclSEs;
	bool inclGradient;
	bool inclVcov;
	bool badSEWarning;
	bool firstTime;
	size_t numExtraCols;

 public:
	virtual bool accumulateInform() { return false; };
        virtual void initFromFrontend(omxState *, SEXP rObj);
        virtual void computeImpl(FitContext *fc);
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
		int *origRowFreq;
		std::vector<int> origCumSum;
		std::vector<int> resample;
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
	virtual bool accumulateInform() { return false; };
	virtual void initFromFrontend(omxState *, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
	virtual void reportResults(FitContext *fc, MxRList *, MxRList *result);
};

class ComputeGenerateData : public omxCompute {
	typedef omxCompute super;
	std::vector< omxExpectation* > expectations;
	MxRList simData;

 public:
	virtual void initFromFrontend(omxState *globalState, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
        virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *out);
};

class ComputeLoadData : public omxCompute {
	typedef omxCompute super;

	static std::vector<LoadDataProviderBase2*> Providers;
	std::unique_ptr<LoadDataProviderBase2> provider;

	omxData *data;
	bool useOriginalData;

	struct ColumnInvalidator : StateInvalidator {
		typedef StateInvalidator super;
		omxData *data;
		const std::vector< int > &columns;
		ColumnInvalidator(omxState &_st, omxData *_data,
				  const std::vector< int > &_columns) :
			super(_st), data(_data), columns(_columns) {};
		virtual void doData() { data->invalidateColumnsCache(columns); };
	};

 public:
	static void loadedHook();
	static void addProvider(LoadDataProviderBase2 *ldp) { Providers.push_back(ldp); }
	virtual void initFromFrontend(omxState *globalState, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *);
};

std::vector<LoadDataProviderBase2*> ComputeLoadData::Providers;

void ComputeLoadDataLoadedHook()
{ ComputeLoadData::loadedHook(); }

class ComputeLoadMatrix : public omxCompute {
	typedef omxCompute super;

	enum LoadMethod {
		LoadCSV,
		LoadDataFrame
	} loadMethod;

	std::vector< omxMatrix* > mat;
	std::vector< mini::csv::ifstream* > streams;
	std::vector<bool> hasRowNames;
	bool useOriginalData;
	// origData is not really needed until it is possible to seek backwards
	std::vector<Eigen::MatrixXd> origData;
	int line;
	Rcpp::DataFrame observed;

	void loadFromCSV(FitContext *fc, int index);
	void loadDataFrame(FitContext *fc, int index);

 public:
	virtual ~ComputeLoadMatrix();
	virtual void initFromFrontend(omxState *globalState, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
};

class ComputeLoadContext : public omxCompute {
	typedef omxCompute super;

	int verbose;
	int loadCounter;
	char sep;
	bool header;
	std::vector<const char *> colnames;
	std::string path;
	std::unique_ptr< mini::csv::ifstream > st;
	int cpIndex;
	int numColumns;
	int *columnPtr;
	int maxColumn;
	int curLine;

	void reopen();

 public:
	virtual void initFromFrontend(omxState *globalState, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void reportResults(FitContext *fc, MxRList *slots, MxRList *);
};

class ComputeTryCatch : public omxCompute {
	typedef omxCompute super;
	std::unique_ptr< omxCompute > plan;
	int cpIndex;

public:
	virtual void initFromFrontend(omxState *, SEXP rObj);
	virtual void computeImpl(FitContext *fc);
	virtual void collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out);
};

static class omxCompute *newComputeSequence()
{ return new omxComputeSequence(); }

static class omxCompute *newComputeIterate()
{ return new omxComputeIterate(); }

static class omxCompute *newComputeLoop()
{ return new ComputeLoop(); }

static class omxCompute *newComputeOnce()
{ return new omxComputeOnce(); }

static class omxCompute *newComputeEM()
{ return new ComputeEM(); }

static class omxCompute *newComputeHessianQuality()
{ return new ComputeHessianQuality(); }

static class omxCompute *newComputeReportDeriv()
{ return new ComputeReportDeriv(); }

static class omxCompute *newComputeReportExpectation()
{ return new ComputeReportExpectation(); }

static class omxCompute *newComputeBootstrap()
{ return new ComputeBootstrap(); }

static class omxCompute *newComputeGenerateData()
{ return new ComputeGenerateData(); }

static class omxCompute *newComputeLoadMatrix()
{ return new ComputeLoadMatrix(); }

static class omxCompute *newComputeCheckpoint()
{ return new ComputeCheckpoint(); }

struct omxComputeTableEntry {
        char name[32];
        omxCompute *(*ctor)();
};

static const struct omxComputeTableEntry omxComputeTable[] = {
        {"MxComputeNumericDeriv", &newComputeNumericDeriv},
        {"MxComputeGradientDescent", &newComputeGradientDescent},
	{"MxComputeSequence", &newComputeSequence },
	{"MxComputeIterate", &newComputeIterate },
	{"MxComputeLoop", &newComputeLoop },
	{"MxComputeOnce", &newComputeOnce },
        {"MxComputeNewtonRaphson", &newComputeNewtonRaphson},
        {"MxComputeEM", &newComputeEM },
	{"MxComputeStandardError",
	 []()->omxCompute* { return new ComputeStandardError; }},
	{"MxComputeHessianQuality", &newComputeHessianQuality},
	{"MxComputeReportDeriv", &newComputeReportDeriv},
	{"MxComputeConfidenceInterval", &newComputeConfidenceInterval},
	{"MxComputeReportExpectation", &newComputeReportExpectation},
	{"MxComputeTryHard", &newComputeTryHard},
	{"MxComputeNelderMead", &newComputeNelderMead},
	{"MxComputeBootstrap", &newComputeBootstrap},
	{"MxComputeGenerateData", &newComputeGenerateData},
	{"MxComputeLoadData",
	 []()->omxCompute* { return new ComputeLoadData; }},
	{"MxComputeLoadMatrix", &newComputeLoadMatrix},
	{"MxComputeCheckpoint", newComputeCheckpoint},
	{"MxComputeSimAnnealing", &newComputeGenSA},
	{"MxComputeJacobian",
	 []()->omxCompute* { return new ComputeJacobian; }},
	{"MxComputeSetOriginalStarts",
	 []()->omxCompute* { return new ComputeSetOriginalStarts; }},
	{"MxComputeTryCatch",
	 []()->omxCompute* { return new ComputeTryCatch; }},
	{"MxComputeLoadContext",
	 []()->omxCompute* { return new ComputeLoadContext; }}
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

        if (!got) mxThrow("Compute plan step '%s' is not implemented", type);

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

void ComputeTryCatch::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	auto &cp = Global->checkpointColnames;
	cpIndex = cp.size();
	std::string errCol(string_snprintf("catch%d", int(Global->computeLoopIndex.size())));
	cp.push_back(errCol);

	SEXP slotValue;
	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("plan")));
	SEXP s4class;
	Rf_protect(s4class = STRING_ELT(Rf_getAttrib(slotValue, R_ClassSymbol), 0));
	plan = std::unique_ptr< omxCompute >(omxNewCompute(globalState, CHAR(s4class)));
	plan->initFromFrontend(globalState, slotValue);
}

void ComputeTryCatch::computeImpl(FitContext *fc)
{
	auto &cv = Global->checkpointValues;
	cv[cpIndex] = "";
	try {
		plan->compute(fc);
	} catch( std::exception &ex ) {
		cv[cpIndex] = ex.what();
	} catch(...) {
		cv[cpIndex] = "c++ exception (unknown reason)";
	}
	if (isErrorRaisedIgnTime()) {
		cv[cpIndex] = Global->getBads();
		Global->bads.clear();
	}
	Global->throwOnUserInterrupted();
}

void ComputeTryCatch::collectResults(FitContext *fc, LocalComputeResult *lcr, MxRList *out)
{
	super::collectResults(fc, lcr, out);
	std::vector< omxCompute* > clist(1);
	clist[0] = plan.get();
	collectResultsHelper(fc, clist, lcr, out);
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
	if (std::isfinite(tolerance) && tolerance <= 0) mxThrow("tolerance must be positive");
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

void ComputeLoop::initFromFrontend(omxState *globalState, SEXP rObj)
{
	SEXP slotValue;

	super::initFromFrontend(globalState, rObj);

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("maxIter")));
	maxIter = INTEGER(slotValue)[0];
	}

	{
		ProtectedSEXP RstartFrom(R_do_slot(rObj, Rf_install("startFrom")));
		startFrom = Rf_asInteger(RstartFrom);
		ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
		verbose = Rf_asInteger(Rverbose);
		ProtectedSEXP Rindices(R_do_slot(rObj, Rf_install("indices")));
		indicesLength = Rf_length(Rindices);
		indices = INTEGER(Rindices);
		ProtectedSEXP RmaxDur(R_do_slot(rObj, Rf_install("maxDuration")));
		maxDuration = Rf_asReal(RmaxDur);
	}

	Rf_protect(slotValue = R_do_slot(rObj, Rf_install("steps")));

	PushLoopIndex pli(name);

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

	iterations = 0;
}

void ComputeLoop::computeImpl(FitContext *fc)
{
	bool hasIndices = indicesLength != 0;
	bool hasMaxIter = maxIter != NA_INTEGER;
	time_t startTime = time(0);
	int lastIndex = indicesLength;
	if (hasMaxIter) lastIndex = std::min(lastIndex, maxIter);
	while (1) {
		PushLoopIndex pli(name, hasIndices? indices[iterations] : startFrom+iterations,
											iterations, lastIndex);
		++iterations;
		++fc->iterations;
		for (size_t cx=0; cx < clist.size(); ++cx) {
			clist[cx]->compute(fc);
			if (isErrorRaised()) {
				if (verbose) mxLog("%s: error raised at step %d", name, int(cx));
				break;
			}
		}
		if (std::isfinite(maxDuration) && time(0) - startTime > maxDuration) {
			if (verbose) mxLog("%s: maximum duration", name);
			break;
		}
		if (hasMaxIter && iterations >= maxIter) {
			if (verbose) mxLog("%s: maximum iterations", name);
			break;
		}
		if (hasIndices && iterations >= indicesLength) {
			if (verbose) mxLog("%s: completed todo list", name);
			break;
		}
		if (isErrorRaised()) {
			if (verbose) mxLog("%s: error raised", name);
			break;
		}
		if (!hasMaxIter) {
			auto &clm = Global->computeLoopMax;
			int m1 = clm[clm.size() - 1];
			if (m1) {
				// filled in by sub compute plan
				maxIter = m1;
				hasMaxIter = true;
			}
		}
	}
}

void ComputeLoop::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	MxRList output;
	output.add("iterations", Rf_ScalarInteger(iterations));
	slots->add("output", output.asR());
}

ComputeLoop::~ComputeLoop()
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
	if (maxIter < 0) mxThrow("maxIter must be non-negative");
	}

	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("tolerance")));
	tolerance = REAL(slotValue)[0];
	if (tolerance <= 0) mxThrow("tolerance must be positive");
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
					if (!ff->fitFunction) mxThrow("infoArgs$fitfunction is %s, not a fitfunction", ff->name());
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
					if (!ff->fitFunction) mxThrow("infoArgs$fitfunction is %s, not a fitfunction", ff->name());
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
						mxThrow("Unknown SEM method '%s'", methodName);
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
				if (semTolerance <= 0) mxThrow("semTolerance must be positive");
			} else if (strEQ(key, "noiseTarget")) {
				noiseTarget = REAL(slotValue)[0];
				if (noiseTarget <= 0) mxThrow("noiseTarget must be positive");
			} else if (strEQ(key, "noiseTolerance")) {
				noiseTolerance = REAL(slotValue)[0];
				if (noiseTolerance < 1) mxThrow("noiseTolerance must be >=1");
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
	if (dist < tolerance/4) mxThrow("SEM: invalid probe offset distance %.9f", dist);
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
		mxThrow("Unknown information method %d", information);
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
	}
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

	fc->gradZ = Eigen::VectorXd::Zero(fc->numParam);
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_GRADIENT, fc);
	}
	result = fc->gradZ;
	reportProgress(fc);
}

void ComputeEM::Oakes(FitContext *fc)
{
	if (verbose >= 1) mxLog("ComputeEM: Oakes1999 method=simple");

	int wanted = fc->wanted;
	const int freeVars = (int) fc->varGroup->vars.size();

	estep->compute(fc);
	fc->wanted &= ~FF_COMPUTE_HESSIAN;  // discard garbage

	fc->initGrad();
	for (size_t fx=0; fx < infoFitFunction.size(); ++fx) {
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_PREOPTIMIZE, fc);
		omxFitFunctionCompute(infoFitFunction[fx]->fitFunction, FF_COMPUTE_GRADIENT, fc);
	}

	Eigen::VectorXd optimumCopy = optimum;  // will be modified
	Eigen::VectorXd refGrad(freeVars);
	refGrad = fc->gradZ;
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
			std::vector<double> coefHist(agileMaxIter);
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

	ThinMatrix rijMat(rij.data(), freeVars, freeVars);
	ThinMatrix hessMat(hess, freeVars, freeVars);
	std::vector<double> infoBuf(freeVars * freeVars);
	ThinMatrix infoMat(infoBuf.data(), freeVars, freeVars);

	SymMatrixMultiply('L', hessMat, rijMat, infoMat);  // result not symmetric!

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
		ThinMatrix ihessMat(ihess, freeVars, freeVars);
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
		Eigen::Map< Eigen::MatrixXd > mat(ihess, freeVars, freeVars);
		ForceInvertSymmetricPosDef(mat, oev, &fc->infoCondNum);
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
		mxThrow("Unknown information matrix estimation method '%s'", iMethod);
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
		mxThrow("MxComputeOnce cannot evaluate expectations and fitfunctions at the same time");
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

		if (hessian && infoMat) mxThrow("Cannot compute the Hessian and Fisher Information matrix simultaneously");
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
			mxThrow("Gradient requested but not available");
		}
		if ((hessian || ihessian || hgprod) && (!ff || !ff->hessianAvailable)) {
			// add a separate flag for hgprod TODO
			mxThrow("Hessian requested but not available");
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
			fc->initGrad();
		}
		if (hessian) {
			want |= FF_COMPUTE_HESSIAN;
			fc->clearHessian();
		}
		if (infoMat) {
			want |= FF_COMPUTE_INFO;
			fc->infoMethod = infoMethod;
			fc->initGrad();
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
		if (predict.size() > 1) mxThrow("Not implemented");
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

void ComputeJacobian::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	ProtectedSEXP Rof(R_do_slot(rObj, Rf_install("of")));
	int numOf = Rf_length(Rof);
	if (!numOf) mxThrow("%s: must provide at least one expectation", name);
	exList.reserve(numOf);
	for (int ex=0; ex < numOf; ++ex) {
		int objNum = INTEGER(Rof)[ex];
		if (objNum < 0) {
			omxExpectation *e1 = state->expectationList[~objNum];
			omxCompleteExpectation(e1);
			exList.push_back(e1);
		} else {
			omxMatrix *algebra = state->algebraList[objNum];
			if (algebra->fitFunction) {
				omxCompleteFitFunction(algebra);
			}
			alList.push_back(algebra);
		}
	}

	if (exList.size()) {
		sense.attach(&exList, 0);
	} else {
		sense.attach(0, &alList);
	}

	data = 0;
	ProtectedSEXP Rdata(R_do_slot(rObj, Rf_install("data")));
	int objNum = Rf_asInteger(Rdata);
	if (objNum != NA_INTEGER) data = state->dataList[objNum];

	ProtectedSEXP Rdefvar_row(R_do_slot(rObj, Rf_install("defvar.row")));
	sense.defvar_row = Rf_asInteger(Rdefvar_row);
}

void ComputeJacobian::computeImpl(FitContext *fc)
{
	int numFree = fc->calcNumFree();
	Eigen::Map< Eigen::VectorXd > curEst(fc->est, numFree);
	if (sense.defvar_row != NA_INTEGER) {
		data->loadDefVars(fc->state, sense.defvar_row - 1);
	}
	sense.measureRef(fc);
	fd_jacobian<false>(GradientAlgorithm_Forward, 2, 1e-4, sense, sense.ref, curEst, sense.result);
}

void ComputeJacobian::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	MxRList output;
	output.add("jacobian", Rcpp::wrap(sense.result));
	slots->add("output", output.asR());
}

void FitContext::resetToOriginalStarts()
{
	setInform(INFORM_UNINITIALIZED);
	auto &startingValues = Global->startingValues;
	auto &vars = varGroup->vars;
	for (int vx=0; vx < int(vars.size()); ++vx) {
		auto *fv = vars[vx];
		est[vx] = startingValues[fv->id];
	}
	fit = NA_REAL;
	mac = NA_REAL;
	fitUnits = FIT_UNITS_UNINITIALIZED;
	skippedRows = 0;
	vcov.resize(0,0);
	stderrs.resize(0);
	clearHessian();
	resetIterationError();
}

void ComputeSetOriginalStarts::computeImpl(FitContext *fc)
{
	fc->resetToOriginalStarts();
}

void ComputeStandardError::initFromFrontend(omxState *state, SEXP rObj)
{
	super::initFromFrontend(state, rObj);

	wlsStats = false;
	fitMat = omxNewMatrixFromSlot(rObj, state, "fitfunction");
}

template <typename T1> bool isULS(const Eigen::MatrixBase<T1> &acov)
{
	for (int cx=0; cx < acov.cols(); ++cx) {
		for (int rx=cx; rx < acov.rows(); ++rx) {
			if (acov(rx,cx) == 0.0) continue;
			if (cx == rx && acov(rx,cx) == 1.0) continue;
			return false;
		}
	}
	return true;
}

void ComputeStandardError::computeImpl(FitContext *fc)
{
	if (fc->fitUnits == FIT_UNITS_UNINITIALIZED) return;
	int numFree = fc->calcNumFree();
	if (fc->fitUnits == FIT_UNITS_MINUS2LL) {
		if (!fc->vcov.size()) {
			fc->vcov.resize(numFree, numFree);
			const double Scale = fabs(Global->llScale);
			fc->refreshDenseIHess();
			//fc->ihess is not actually the inverted Hessian,
			//but there's a method for constructing the inverted Hessian from it:
			fc->copyDenseIHess(fc->vcov.data());
			fc->vcov = Scale * fc->vcov;
		}
		fc->calcStderrs();
		return;
	}

	if (!(fc->fitUnits == FIT_UNITS_SQUARED_RESIDUAL ||
	      fc->fitUnits == FIT_UNITS_SQUARED_RESIDUAL_CHISQ)) return;
	if (!fitMat) return;

	exList.clear();
	std::function<void(omxMatrix*)> ve = visitEx(this);
	fitMat->fitFunction->traverse(ve);

	double totalWeight = 0;
	int numOrdinal = 0;
	int totalStats = 0;
	numStats.reserve(exList.size());
	for (auto &e1 : exList) {
		e1->data->visitObsStats([&](obsSummaryStats &o1){
				numOrdinal += o1.numOrdinal;
				totalWeight += o1.totalWeight;
				int sz = o1.fullWeight->rows;
				numStats.push_back(sz);
				totalStats += sz;
			});
	}

	Eigen::VectorXd obStats(totalStats);
	Eigen::VectorXd exStats(totalStats);
	Eigen::MatrixXd Vmat(totalStats,totalStats);
	Eigen::MatrixXd Wmat(totalStats,totalStats);
	Vmat.setZero();
	Wmat.setZero();
	for (int ex=0, sx=0, offset=0; ex < int(exList.size()); ++ex) {
		omxData *d1 = exList[ex]->data;

		d1->visitObsStats([&](obsSummaryStats &o1){
				int sz = numStats[sx];

				// filter out missing variables TODO
				Eigen::VectorXd vec1(sz);
				exList[ex]->asVector(fc, 0, vec1);
				exStats.segment(offset, sz) = vec1;
				if (o1.thresholdMat) {
					EigenMatrixAdaptor oTh(o1.thresholdMat);
          auto &oThresh = o1.thresholdCols;
					normalToStdVector(o1.covMat, o1.meansMat, o1.slopeMat,
                            [&oThresh, &oTh](int r,int c)->double{ return oTh(r, oThresh[c].column); },
														oThresh, vec1);
				} else {
					normalToStdVector(o1.covMat, o1.meansMat, o1.slopeMat,
														[](int r,int c)->double{ return 0; },
														o1.thresholdCols, vec1);
				}
				obStats.segment(offset, sz) = vec1;

				if (o1.acovMat) {
					EigenMatrixAdaptor acov(o1.acovMat);
					Vmat.block(offset,offset,sz,sz) = acov;
				} else {
					Vmat.block(offset,offset,sz,sz).setIdentity();
				}

				EigenMatrixAdaptor fw(o1.fullWeight);
				int nonZeroDims = (fw.diagonal().array() != 0.0).count();
				if (nonZeroDims == 0) {
					omxRaiseErrorf("%s: fullWeight matrix is missing", exList[ex]->name);
					return;
				}
				Eigen::MatrixXd ifw(fw.rows(), fw.cols());
				ifw.setZero();
				Eigen::MatrixXd dense;
				subsetCovariance(fw, [&](int xx){ return fw.diagonal().coeff(xx,xx) != 0.0; },
						 nonZeroDims, dense);
				if (InvertSymmetricIndef(dense, 'L') > 0) {
					MoorePenroseInverse(dense);
				}
				Eigen::MatrixXd idense = dense.selfadjointView<Eigen::Lower>();
				subsetCovarianceStore(ifw,
						      [&](int xx){ return fw.diagonal().coeff(xx,xx) != 0.0; },
						      idense);
				Wmat.block(offset,offset,sz,sz) = ifw;
				offset += sz;
				sx += 1;
			});
	}

	Eigen::Map< Eigen::VectorXd > curEst(fc->est, numFree);
	ParJacobianSense sense;
	sense.attach(&exList, 0);
	sense.measureRef(fc);
	fd_jacobian<false>(GradientAlgorithm_Forward, 2, 1e-4, sense, sense.ref, curEst, sense.result);

	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr1(sense.result);
	Eigen::MatrixXd q1 = qr1.householderQ();
	Eigen::MatrixXd jacOC = q1.block(0, qr1.rank(), q1.rows(), q1.cols() - qr1.rank());
	if (jacOC.cols()) {
		Eigen::MatrixXd zqb = jacOC.transpose() * Wmat * jacOC;
		MoorePenroseInverse(zqb);

		Eigen::VectorXd diff = obStats - exStats;
		x2 = diff.transpose() * jacOC * zqb * jacOC.transpose() * diff;
		Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr2(jacOC);
		df = qr2.rank();
	} else {
		x2 = 0;
		df = 0;
	}

	Eigen::MatrixXd dvd = sense.result.transpose() * Vmat * sense.result;
	if (InvertSymmetricIndef(dvd, 'L') > 0) return;

	fc->vcov = dvd.selfadjointView<Eigen::Lower>() * sense.result.transpose() *
		Vmat * Wmat * Vmat * sense.result * dvd.selfadjointView<Eigen::Lower>();
	fc->calcStderrs();

	Eigen::MatrixXd Umat = Vmat - Vmat * sense.result * dvd.selfadjointView<Eigen::Lower>() *
		sense.result.transpose() * Vmat;
	Eigen::MatrixXd UW = Umat * Wmat;
	Eigen::MatrixXd UW2 = UW * UW; // unclear if this should be UW^2 i.e. elementwise power
	double trUW = UW.diagonal().array().sum();
	madj = trUW / df;
	x2m = fc->fit / madj;
	dstar = round((trUW * trUW) / UW2.diagonal().array().sum());
	mvadj = (trUW*trUW) / dstar;
	x2mv = fc->fit / mvadj;
	// N.B. x2mv is off by a factor of N where N is the total number of rows in all data sets for the ULS case.
	if (isULS(Vmat)) x2mv /= totalWeight;
	wlsStats = true;
}

void ComputeStandardError::reportResults(FitContext *fc, MxRList *slots, MxRList *out)
{
	int numFree=0;
	SEXP parNames=0, dimnames=0;

	if (fc->vcov.size() || fc->stderrs.size()) {
		numFree = fc->calcNumFree();
		if (numFree != fc->stderrs.size()) {
			mxThrow("%s: numFree != fc->stderrs.size() %d != %d",
				name, numFree, fc->stderrs.size());
		}

		parNames = Rf_allocVector(STRSXP, numFree);
		Rf_protect(parNames);
		for (int vx=0, px=0; vx < int(fc->numParam) && px < numFree; ++vx) {
			if (fc->profiledOut[vx]) continue;
			SET_STRING_ELT(parNames, px++, Rf_mkChar(varGroup->vars[vx]->name));
		}

		dimnames = Rf_allocVector(VECSXP, 2);
		Rf_protect(dimnames);
		SET_VECTOR_ELT(dimnames, 0, parNames);
	}

	if (fc->vcov.size()) {
		SEXP Vcov;
		Rf_protect(Vcov = Rf_allocMatrix(REALSXP, fc->vcov.rows(), fc->vcov.cols()));
		memcpy(REAL(Vcov), fc->vcov.data(), sizeof(double) * fc->vcov.rows() * fc->vcov.cols());
		Rf_setAttrib(Vcov, R_DimNamesSymbol, dimnames);
		out->add("vcov", Vcov);
	}

	if (fc->stderrs.size()) {
		SEXP stdErrors;
		Rf_protect(stdErrors = Rf_allocMatrix(REALSXP, numFree, 1));
		memcpy(REAL(stdErrors), fc->stderrs.data(), sizeof(double) * numFree);
		Rf_setAttrib(stdErrors, R_DimNamesSymbol, dimnames);
		out->add("standardErrors", stdErrors);
	}

	if (wlsStats) {
		out->add("chi", Rf_ScalarReal(x2));
		out->add("chiDoF", Rf_ScalarInteger(df));
		out->add("chiM", Rf_ScalarReal(x2m));
		out->add("chiMV", Rf_ScalarReal(x2mv));
		out->add("chiMadjust", Rf_ScalarReal(madj));
		out->add("chiMVadjust", Rf_ScalarReal(mvadj));
		out->add("chiDoFstar", Rf_ScalarReal(dstar));
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

	/*
	 * If there are equality MxConstraints, then the quality of the generically calculated Hessian will not be informative;
	 * we must rely upon the optimizer's status code.
	 * (Completely boneheaded equality MxConstraints, like constraining a free parameter equal to itself, are assumed to be the user's problem.)
	 * The quality of the generically calculated Hessian is likewise uninformative if there are active inequality MxConstraints.
	 * But, if there are no equalities, and all of the inequalities are inactive, then we might as well proceed (corner case as that may be):
	*/
	if(fc->state->conListX.size()){
		//Any compute step that cares about how many constraints there are needs to ask the omxState to recount:
		fc->state->countNonlinearConstraints(fc->state->numEqC, fc->state->numIneqC, false);
		int nf = fc->calcNumFree();
		fc->inequality.resize(fc->state->numIneqC);
		fc->analyticIneqJacTmp.resize(fc->state->numIneqC, nf);
		fc->myineqFun(true, verbose, omxConstraint::LESS_THAN, false);
		if(fc->state->numEqC || fc->inequality.array().sum()){return;}
	}

	// See Luenberger & Ye (2008) Second Order Test (p. 190) and Condition Number (p. 239)

	if (fc->infoDefinite != NA_LOGICAL || !doubleEQ(fc->infoCondNum, NA_REAL)) {
		if (verbose >= 1) {
			mxLog("%s: information matrix already determined to be non positive definite; skipping", name);
		}
		return; // already set elsewhere
	}

	int numParams = fc->getDenseHessianishSize();
	if (numParams == 0) return;

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
	if( fc->state->conListX.size() ){
		/* After the call to the backend,
		 * frontend function nameGenericConstraintOutput(), in R/MxRunHelperFunctions.R, uses 'constraintNames',
		 * 'constraintRows', and 'constraintCols' to populate the dimnames of 'constraintFunctionValues' and
		 * 'constraintJacobian'.
		 */
		SEXP cn, cr, cc, cv, cjac;
		size_t i=0;
		{
			Rf_protect(cn = Rf_allocVector( STRSXP, fc->state->conListX.size() ));
			Rf_protect(cr = Rf_allocVector( INTSXP, fc->state->conListX.size() ));
			Rf_protect(cc = Rf_allocVector( INTSXP, fc->state->conListX.size() ));
			for(i=0; i < fc->state->conListX.size(); i++){
				SET_STRING_ELT( cn, i, Rf_mkChar(fc->state->conListX[i]->name) );
				INTEGER(cr)[i] = fc->state->conListX[i]->nrows;
				INTEGER(cc)[i] = fc->state->conListX[i]->ncols;
			}
			result->add("constraintNames", cn);
			result->add("constraintRows", cr);
			result->add("constraintCols", cc);
		}
		if( fc->constraintFunVals.size() ){
			Rf_protect(cv = Rf_allocVector( REALSXP, fc->constraintFunVals.size() ));
			memcpy( REAL(cv), fc->constraintFunVals.data(), sizeof(double) * fc->constraintFunVals.size() );
			result->add("constraintFunctionValues", cv);
		}
		if( fc->constraintJacobian.rows() ){
			Rf_protect(cjac = Rf_allocMatrix( REALSXP, fc->constraintJacobian.rows(), fc->constraintJacobian.cols() ));
			memcpy( REAL(cjac), fc->constraintJacobian.data(), sizeof(double) * fc->constraintJacobian.rows() * fc->constraintJacobian.cols() );
			result->add("constraintJacobian", cjac);
		}
	}

	if (!(fc->wanted & (FF_COMPUTE_GRADIENT|FF_COMPUTE_HESSIAN|FF_COMPUTE_IHESSIAN))) return;

	int numFree = fc->calcNumFree();

	SEXP parNames = Rf_allocVector(STRSXP, numFree);
	Rf_protect(parNames);
	for (int vx=0, px=0; vx < int(fc->numParam) && px < numFree; ++vx) {
		if (fc->profiledOut[vx]) continue;
		SET_STRING_ELT(parNames, px++, Rf_mkChar(varGroup->vars[vx]->name));
	}

	if (fc->wanted & FF_COMPUTE_GRADIENT) {
		SEXP Rgradient = Rf_allocVector(REALSXP, numFree);
		result->add("gradient", Rgradient);
		double *gp = REAL(Rgradient);
		fc->copyGradToOptimizer(gp);
		Rf_setAttrib(Rgradient, R_NamesSymbol, parNames);
	}

	if (fc->wanted & (FF_COMPUTE_HESSIAN | FF_COMPUTE_IHESSIAN)) {
		SEXP dimnames = Rf_allocVector(VECSXP, 2);
		Rf_protect(dimnames);
		for (int dx=0; dx < 2; ++dx) SET_VECTOR_ELT(dimnames, dx, parNames);

		if (numFree != fc->hess.rows()) {
			if (OMX_DEBUG) mxLog("free %d hessian %d", numFree, fc->hess.rows());
			return;
		}
		if (fc->wanted & FF_COMPUTE_HESSIAN) {
			SEXP Rhessian;
			Rhessian = Rf_allocMatrix(REALSXP, numFree, numFree);
			result->add("hessian", Rhessian);
			fc->copyDenseHess(REAL(Rhessian));
			Rf_setAttrib(Rhessian, R_DimNamesSymbol, dimnames);
		}
		if (numFree != fc->ihess.rows()) {
			if (OMX_DEBUG) mxLog("free %d ihessian %d", numFree, fc->ihess.rows());
			return;
		}
		if (fc->wanted & FF_COMPUTE_IHESSIAN) {
			SEXP Rihessian;
			Rihessian = Rf_allocMatrix(REALSXP, numFree, numFree);
			result->add("ihessian", Rihessian);
			fc->copyDenseIHess(REAL(Rihessian));
			Rf_setAttrib(Rihessian, R_DimNamesSymbol, dimnames);
		}
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
		omxExpectationCompute(fc, curExpectation);
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
			mxThrow("%s: data '%s' of type '%s' cannot have row weights",
				 name, ctx.data->name, ctx.data->getType());
		}
		ctx.origRowFreq = ctx.data->getFreqColumn();
		ctx.origCumSum.resize(numRows);
		ctx.resample.resize(ctx.origCumSum.size());
		if (ctx.origRowFreq) {
			std::partial_sum(ctx.origRowFreq, ctx.origRowFreq + ctx.origCumSum.size(),
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
	SEXP colNames;
	{
		ProtectedSEXP colNamesP(Rf_allocVector(STRSXP, numCols));
		Rf_setAttrib(rawOutput, R_NamesSymbol, colNamesP);
		colNames = colNamesP;
	}

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
		BorrowRNGState grs;
		for (int repl=0; repl < numReplications; ++repl) {
			if (seedVec[repl] != NA_INTEGER) continue;
			int seed1 = unif_rand() * std::numeric_limits<int>::max();
			if (seed1 == NA_INTEGER) seed1 = 0; // maybe impossible
			seedVec[repl] = seed1;
		}
	} else {
		if (only <= Rf_length(VECTOR_ELT(previousData, 0))) {
			if (verbose >= 1) mxLog("%s: using only=%d", name, only);
			seedVec[0] = INTEGER(VECTOR_ELT(previousData, 0))[only - 1];
		} else {
			mxThrow("%s: only=%d but previous data has just %d replications",
				 name, only, Rf_length(VECTOR_ELT(previousData, 0)));
		}
	}

	for (int repl=0; repl < numReplications && !isErrorRaised(); ++repl) {
		std::mt19937 generator(seedVec[repl]);
		if (INTEGER(VECTOR_ELT(rawOutput, 2 + fc->numParam))[repl] != NA_INTEGER) continue;
		if (verbose >= 2) mxLog("%s: replication %d", name, repl);
		fc->state->invalidateCache();
		for (auto &ctx : contexts) {
			ctx.resample.assign(ctx.origCumSum.size(), 0);
			int last = ctx.origCumSum.size() - 1;
			int total = ctx.origCumSum[last];
			std::uniform_int_distribution<int> dist(1, total);
			for (int sx=0; sx < total; ++sx) {
				int pick = dist(generator);
				auto rowPick = std::lower_bound(ctx.origCumSum.begin(), ctx.origCumSum.end(), pick);
				int row = rowPick - ctx.origCumSum.begin();
				ctx.resample[row] += 1;
			}
			ctx.data->setFreqColumn(ctx.resample.data());
			if (only != NA_INTEGER) {
				onlyWeight.add(ctx.data->name, Rcpp::wrap(ctx.resample));
			}
		}
		fc->state->connectToData();
		fc->resetToOriginalStarts();
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
		ctx.data->setFreqColumn(ctx.origRowFreq);
	}

	if (only == NA_INTEGER) {
		fc->resetToOriginalStarts();
		fc->copyParamToModel();
		fc->wanted &= ~FF_COMPUTE_FIT;  // discard garbage
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
		output.add("frequency", onlyWeight.asR());
	}
	slots->add("output", output.asR());
}

void ComputeGenerateData::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	ProtectedSEXP Rexp(R_do_slot(rObj, Rf_install("expectation")));
	for (int wx=0; wx < Rf_length(Rexp); ++wx) {
		if (isErrorRaised()) return;
		int objNum = INTEGER(Rexp)[wx];
		omxExpectation *expectation = globalState->expectationList[objNum];
		expectations.push_back(expectation);
	}
}

void ComputeGenerateData::computeImpl(FitContext *fc)
{
	if (simData.size()) mxThrow("Cannot generate data more than once");

	BorrowRNGState grs;

	for (auto ex : expectations) {
		ex->generateData(fc, simData);
	}
}

void ComputeGenerateData::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	slots->add("output", simData.asR());
}

class LoadDataCSVProvider : public LoadDataProvider<LoadDataCSVProvider> {
	std::unique_ptr< mini::csv::ifstream > icsv;
	int cpIndex;
	bool byrow;

	virtual const char *getName() { return "csv"; };
	virtual void init(SEXP rObj) {
		ProtectedSEXP Rbyrow(R_do_slot(rObj, Rf_install("byrow")));
		byrow = Rf_asLogical(Rbyrow);
		ProtectedSEXP Rcs(R_do_slot(rObj, Rf_install("cacheSize")));

		if (!byrow) stripeSize = std::max(Rf_asInteger(Rcs), 1);
		requireFile(rObj);
	}
	virtual void addCheckpointColumns(std::vector< std::string > &cp)
	{
		if (rowNames == 0 || !byrow) return;
		cpIndex = cp.size();
		auto rc = *rawCols;
		for (int cx=0; cx < int(columns.size()); ++cx) {
			std::string c1 = fileName + ":" + rc[ columns[cx] ].name;
			cp.push_back(c1);
		}
	}
	void loadByCol(int index);
	void loadByRow(int index);
	virtual void loadRowImpl(int index)
	{
		if (!byrow) {
			loadByCol(index);
		} else {
			loadByRow(index);
		}
	}
	void mxScanInt(mini::csv::ifstream &st, ColumnData &rc, int *out);
};

void LoadDataCSVProvider::mxScanInt(mini::csv::ifstream &st, ColumnData &rc, int *out)
{
	const std::string &rn = st.get_delimited_str();
	if (isNA(rn)) {
		*out = NA_INTEGER;
		return;
	}
	if (rc.levels.size()) {
		bool found = false; // maybe use a map for better performance?
		for (int lx=0; lx < int(rc.levels.size()); ++lx) {
			if (rn == rc.levels[lx]) {
				found = true;
				*out = 1+lx;
				break;
			}
		}
		if (!found) mxThrow("%s: factor level '%s' unrecognized in column '%s'",
				    name, rn.c_str(), rc.name);
	} else {
		std::istringstream is(rn);
		is >> *out;
	}
}

void LoadDataCSVProvider::loadByCol(int index)
{
	if (stripeStart == -1 ||
	    index < stripeStart || index >= stripeEnd) {
		bool backward = index < stripeStart;
		stripeStart = std::max(0, index - backward * (stripeSize - 1));

		loadCounter += 1;
		mini::csv::ifstream st(filePath);
		st.set_delimiter(' ', "##");
		for (int rx=0; rx < skipRows; ++rx) st.skip_line();
		int stripeAvail = stripeSize;
		for (int rx=0,dr=0; rx < srcRows; ++rx) {
			bool gotLine = false;
			try {
				gotLine = st.read_line();
			} catch (...) {
				// ignore
			}
			if (!gotLine) {
				mxThrow("%s: ran out of data for '%s' (need %d rows but only found %d)",
					 name, dataName, srcRows, 1+rx);
			}
			if (skipRow(rx)) continue;
			int toSkip = stripeStart * columns.size() + skipCols;
			for (int jx=0; jx < toSkip; ++jx) {
				std::string rn;
				st >> rn;
			}
			for (int sx=0,dx=0; sx < stripeAvail; ++sx) {
				auto rc = *rawCols;
				try {
					for (int cx=0; cx < int(columns.size()); ++cx) {
						if (colTypes[cx] == COLUMNDATA_NUMERIC) {
							st >> stripeData[dx].realData[dr];
						} else {
							mxScanInt(st, rc[ columns[cx] ],
								  &stripeData[dx].intData[dr]);
						}
						dx += 1;
					}
				} catch (...) {
					// assume we tried to read off the end of the line
					stripeAvail = sx;
				}
			}
			dr += 1;
		}
		stripeEnd = stripeStart + stripeAvail;
		if (verbose >= 2) {
			mxLog("%s: loaded stripes [%d,%d) of %d columns each",
			      name, stripeStart, stripeEnd, int(columns.size()));
		}
	}

	if (index < stripeStart || index >= stripeEnd) {
		mxThrow("%s: no data available for %d", name, index);
	}

	int offset = (index - stripeStart) * columns.size();
	auto &rc = *rawCols;
	for (int cx=0; cx < int(columns.size()); ++cx) {
		rc[ columns[cx] ].setBorrow(stripeData[offset + cx]);
	}
}

void LoadDataCSVProvider::loadByRow(int index)
{
	auto &rc = *rawCols;
	if (!icsv || index < curRecord) {
		icsv = std::unique_ptr< mini::csv::ifstream >(new mini::csv::ifstream(filePath));
		icsv->set_delimiter(' ', "##");
		for (int rx=0; rx < skipRows; ++rx) {
			if (1+rx == colNames) {
				// TODO
			}
			icsv->skip_line();
		}
		curRecord = 0;
		loadCounter += 1;
	}

	while (index > curRecord) {
		icsv->skip_line();
		curRecord += 1;
	}
	for (int cx=0; cx < int(columns.size()); ++cx) {
		if (!icsv->read_line()) {
			mxThrow("%s: ran out of data for '%s' at record %d",
				name, dataName, 1+index);
		}
		for (int sx=0; sx < skipCols; ++sx) {
			std::string rn;
			*icsv >> rn;
			if (checkpoint && 1+sx == rowNames) {
				auto &cv = *checkpointValues;
				cv[cpIndex + cx] = rn;
			}
		}
		if (colTypes[cx] == COLUMNDATA_NUMERIC) {
			for (int rx=0, dr=0; rx < srcRows; ++rx) {
				const std::string& str = icsv->get_delimited_str();
				if (skipRow(rx)) continue;
				if (isNA(str)) {
					stripeData[cx].realData[dr] = NA_REAL;
				} else {
					std::istringstream is(str);
					is >> stripeData[cx].realData[dr];
				}
				dr += 1;
			}
		} else {
			for (int rx=0,dr=0; rx < srcRows; ++rx) {
				if (skipRow(rx)) {
					icsv->get_delimited_str();
					continue;
				}
				mxScanInt(*icsv, rc[ columns[cx] ],
					  &stripeData[cx].intData[dr]);
				dr += 1;
			}
		}
	}
	curRecord += 1;
	for (int cx=0; cx < int(columns.size()); ++cx) {
		rc[ columns[cx] ].setBorrow(stripeData[cx]);
	}
}

class LoadDataDFProvider : public LoadDataProvider<LoadDataDFProvider> {
	bool byrow;
	Rcpp::DataFrame observed;

	virtual const char *getName() { return "data.frame"; };
	virtual int getNumVariants() {
		return (observed.nrows() / srcRows) * (observed.ncol() / columns.size());
	}
	virtual void init(SEXP rObj) {
		ProtectedSEXP Rbyrow(R_do_slot(rObj, Rf_install("byrow")));
		byrow = Rf_asLogical(Rbyrow);
		if (byrow) mxThrow("byrow=TRUE not implemented for data.frame data");

		ProtectedSEXP Robs(R_do_slot(rObj, Rf_install("observed")));
		observed = Robs;
		if (int(observed.size()) < int(colTypes.size())) {
			mxThrow("%s: provided observed data only has %d columns but %d requested",
				name, int(observed.size()), int(colTypes.size()));
		}
		if (observed.nrows() % srcRows != 0) {
			mxThrow("%s: original data has %d rows, "
					 "does not divide the number of observed rows %d evenly (remainder %d)",
					 name, srcRows, observed.nrows(), observed.nrows() % srcRows);
		}
		if (observed.ncol() % columns.size() != 0) {
			mxThrow("%s: original data has %d columns, "
					 "does not divide the number of observed columns %d evenly (remainder %d)",
					 name, int(columns.size()), observed.ncol(), observed.ncol() % columns.size());
		}
		if (observed.nrows() != srcRows && observed.ncol() != int(columns.size())) {
			mxThrow("%s: additional data must be in rows or columns, but not both");
		}
		CharacterVector obNames = observed.attr("names");
		for (int cx=0; cx < int(colTypes.size()); ++cx) {
			if (colTypes[cx] == COLUMNDATA_NUMERIC) {
				// OK
			} else {
				IntegerVector vec = observed[cx];
				if (vec.hasAttribute("levels")) {
					CharacterVector lev = vec.attr("levels");
					auto &rc = (*rawCols)[ columns[cx] ];
					if (int(rc.levels.size()) != int(lev.size())) {
						mxThrow("%s: observed column %d (%s) has a different number"
								 "of factor levels %d compare to the original data %d",
								 name, 1+cx, as<const char *>(obNames[cx]),
								 int(lev.size()), int(rc.levels.size()));
					}
				}
			}
		}
	}
	virtual void loadRowImpl(int index)
	{
		auto &rc = *rawCols;
		if (observed.nrows() != srcRows) {
			int rowBase = index * srcRows;
			if (observed.nrows() < rowBase + srcRows) {
				mxThrow("%s: index %d requested but observed data only has %d sets of rows",
						 name, index, observed.nrows() / srcRows);
			}
			for (int cx=0; cx < int(columns.size()); ++cx) {
				RObject vec = observed[cx];
				if (colTypes[cx] == COLUMNDATA_NUMERIC) {
					NumericVector avec(vec);
					for (int rx=0,dr=0; rx < srcRows; ++rx) {
						if (skipRow(rx)) continue;
						stripeData[cx].realData[dr] = avec[rowBase + rx];
						dr += 1;
					}
				} else {
					IntegerVector ivec(vec);
					for (int rx=0,dr=0; rx < srcRows; ++rx) {
						if (skipRow(rx)) continue;
						stripeData[cx].intData[dr] = ivec[rowBase + rx];
						dr += 1;
					}
				}
				rc[ columns[cx] ].setBorrow(stripeData[cx]);
			}
		} else {
			int colBase = index * columns.size();
			if (observed.ncol() < int(colBase + columns.size())) {
				mxThrow("%s: index %d requested but observed data only has %d sets of columns",
						 name, index, observed.ncol() / columns.size());
			}
			for (int cx=0; cx < int(columns.size()); ++cx) {
				RObject vec = observed[colBase + cx];
				if (colTypes[cx] == COLUMNDATA_NUMERIC) {
					NumericVector avec(vec);
					rc[ columns[cx] ].setBorrow(avec.begin());
				} else {
					IntegerVector ivec(vec);
					rc[ columns[cx] ].setBorrow(ivec.begin());
				}
			}
		}
	}
};

void LoadDataProviderBase2::commonInit(SEXP rObj, const char *_name,
                                       const char *_dataName, int _rows,
                                       std::vector<ColumnData> &_rawCols,
                                       ColMapType &_rawColMap,
                                       std::vector< std::string > &_checkpointValues,
                                       bool useOriginalData)
{
	name = _name;
	dataName = _dataName;
	destRows = _rows;
	srcRows = _rows;
	rawCols = &_rawCols;
	rawColMap = &_rawColMap;
	checkpointValues = &_checkpointValues;

	curRecord = -1;
	loadCounter = 0;
	stripeSize = 1;
	stripeStart = -1;
	stripeEnd = -1;

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose);

	rowNames = NA_INTEGER;
	colNames = NA_INTEGER;
	ProtectedSEXP Rrownames(R_do_slot(rObj, Rf_install("row.names")));
	if (Rf_length(Rrownames)) rowNames = Rf_asInteger(Rrownames);
	ProtectedSEXP Rcolnames(R_do_slot(rObj, Rf_install("col.names")));
	if (Rf_length(Rcolnames)) colNames = Rf_asInteger(Rcolnames);

	ProtectedSEXP Rskiprows(R_do_slot(rObj, Rf_install("skip.rows")));
	skipRows = Rf_asInteger(Rskiprows);
	ProtectedSEXP Rskipcols(R_do_slot(rObj, Rf_install("skip.cols")));
	skipCols = Rf_asInteger(Rskipcols);

	ProtectedSEXP RnaStr(R_do_slot(rObj, Rf_install("na.strings")));
	for (int x1=0; x1 < Rf_length(RnaStr); ++x1) {
		naStrings.push_back(R_CHAR(STRING_ELT(RnaStr, x1)));
	}

	ProtectedSEXP Rcol(R_do_slot(rObj, Rf_install("column")));
	for (int cx=0; cx < Rf_length(Rcol); ++cx) {
		auto cn = R_CHAR(STRING_ELT(Rcol, cx));
		auto &rcm = *rawColMap;
		auto rci = rcm.find(cn);
		if (rci == rcm.end()) {
			omxRaiseErrorf("%s: column '%s' not found in '%s'",
				       name, cn, dataName);
			continue;
		}
		columns.push_back(rci->second);
		auto &rc = _rawCols[rci->second];
		colTypes.push_back(rc.type);
		if (useOriginalData) origData.emplace_back(rc.steal());
	}

	ProtectedSEXP Rcheckpoint(R_do_slot(rObj, Rf_install("checkpointMetadata")));
	checkpoint = Rf_asLogical(Rcheckpoint);

	ProtectedSEXP RrowFilter(R_do_slot(rObj, Rf_install("rowFilter")));
	rowFilter = 0;
	if (Rf_length(RrowFilter)) {
		rowFilter = INTEGER(RrowFilter);
		srcRows = Rf_length(RrowFilter);
		int numSkip = std::accumulate(rowFilter, rowFilter + srcRows, 0);
		if (destRows != srcRows - numSkip) {
			mxThrow("rowFilter skips %d rows but %d-%d does not match the number of "
							"rows of observed data %d", numSkip, srcRows, numSkip, destRows);
		}
	}
}

void ComputeLoadData::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	ProtectedSEXP RoriginalData(R_do_slot(rObj, Rf_install("originalDataIsIndexOne")));
	useOriginalData = Rf_asLogical(RoriginalData);

	ProtectedSEXP Rmethod(R_do_slot(rObj, Rf_install("method")));
	const char *methodName = R_CHAR(STRING_ELT(Rmethod, 0));

	ProtectedSEXP Rdata(R_do_slot(rObj, Rf_install("dest")));
	if (Rf_length(Rdata) != 1)
		mxThrow("%s: can only handle 1 destination MxData", name);
	int objNum = Rf_asInteger(Rdata);
	data = globalState->dataList[objNum];
	auto &rd = data->getUnfilteredRawData();

	for (auto pr : Providers) {
		if (strEQ(methodName, pr->getName())) {
			provider = pr->clone();
			provider->commonInit(rObj, name, data->name, rd.rows, rd.rawCols,
                           data->rawColMap, Global->checkpointValues, useOriginalData);
			provider->init(rObj);
			break;
		}
	}
	if (!provider) {
		std::string avail;
		for (auto pr : Providers) {
			avail += " ";
			avail += pr->getName();
		}
		mxThrow("%s: unknown provider '%s'; available providers are:%s",
						name, methodName, avail.c_str());
	}

	if (provider->wantCheckpoint()) {
		auto &cp = Global->checkpointColnames;
		provider->addCheckpointColumns(cp);
	}
}

void ComputeLoadData::computeImpl(FitContext *fc)
{
	std::vector<int> &clc = Global->computeLoopIndex;
	if (clc.size() == 0) mxThrow("%s: must be used within a loop", name);
	int index = clc[clc.size()-1] - 1;  // innermost loop index

	data->setModified();
	if (useOriginalData && index == 0) {
		provider->loadOrigRow();
	} else {
		index -= useOriginalData; // 0 == the first record
		provider->loadRow(index);
		auto &clm = Global->computeLoopMax;
		auto &max = clm.at(clm.size()-1);
		if (max == 0) max = provider->getNumVariants();
	}

	auto &columns = provider->getColumns();
	ColumnInvalidator ci(*fc->state, data, columns);
	ci();
	data->evalAlgebras(fc);
  fc->state->connectToData();
}

void ComputeLoadData::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	MxRList dbg;
	dbg.add("loadCounter", Rf_ScalarInteger(provider->getLoadCounter()));
	slots->add("debug", dbg.asR());
}

void ComputeLoadData::loadedHook()
{
	Providers.clear();
	Providers.push_back(new LoadDataCSVProvider());
	Providers.push_back(new LoadDataDFProvider());
}

void AddLoadDataProvider(double version, int ldpbSz, LoadDataProviderBase2 *ldp)
{
	if (version == OPENMX_LOAD_DATA_API_VERSION) {
		if (ldpbSz != sizeof(LoadDataProviderBase2)) {
			mxThrow("Cannot add mxComputeLoadData provider, version matches "
							"but OpenMx is compiled with different compiler options (%d != %d)",
							ldpbSz, int(sizeof(LoadDataProviderBase2)));
		}
	} else {
		mxThrow("Cannot add mxComputeLoadData provider, version mismatch");
	}
	ComputeLoadData::addProvider(ldp);
}

void ComputeLoadContext::reopen()
{
	loadCounter += 1;
	try {
		st = std::unique_ptr<mini::csv::ifstream>(new mini::csv::ifstream(path));
	} catch (...) {
		mxThrow("%s: failed to open '%s'", name, path.c_str());
	}
	st->set_delimiter(sep, "##");
}

void ComputeLoadContext::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	loadCounter = 0;

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	verbose = Rf_asInteger(Rverbose);
	ProtectedSEXP Rheader(R_do_slot(rObj, Rf_install("header")));
	header = Rf_asInteger(Rheader);
	ProtectedSEXP Rcolnames(R_do_slot(rObj, Rf_install("col.names")));
	for (int cx=0; cx < Rf_length(Rcolnames); ++cx) {
		colnames.push_back(R_CHAR(STRING_ELT(Rcolnames, cx)));
	}

	ProtectedSEXP Rcol(R_do_slot(rObj, Rf_install("column")));
	numColumns = Rf_length(Rcol);
	columnPtr = INTEGER(Rcol);
	if (numColumns == 0) return;

	ProtectedSEXP Rsep(R_do_slot(rObj, Rf_install("sep")));
	const char *sepStr = R_CHAR(STRING_ELT(Rsep, 0));
	if (strlen(sepStr) != 1) mxThrow("%s: sep must be a single character, not '%s'", name, sepStr);
	sep = sepStr[0];

	auto &cp = Global->checkpointColnames;
	cpIndex = cp.size();
	Eigen::Map< Eigen::ArrayXi > col(columnPtr, numColumns);
	if (col.minCoeff() < 1) mxThrow("%s: the first column is 1, not %d",
				   name, col.minCoeff());
	maxColumn = col.maxCoeff();
	//mxLog("%s: len %d max %d", name, numColumns, maxColumn);

	ProtectedSEXP Rpath(R_do_slot(rObj, Rf_install("path")));
	path = R_CHAR(STRING_ELT(Rpath, 0));
	reopen();
	if (header) {
		if (!st->read_line()) mxThrow("%s: cannot read header of '%s'", name, path.c_str());
	}
	if (colnames.size()) {
		for (int cx=0; cx < numColumns; ++cx) cp.push_back(colnames[cx]);
	} else if (header) {
		int xx=0;
		for (int cx=0; cx < maxColumn; ++cx) {
			std::string c1;
			*st >> c1;
			if (cx == col[xx]-1) {
				if (verbose) mxLog("cx %d xx %d %s", cx, xx, c1.c_str());
				cp.push_back(c1);
				if (++xx == numColumns) break;
			}
		}
		if (xx != numColumns) mxThrow("%s: columns must be ordered from first to last", name);
	} else {
		for (int cx=0; cx < numColumns; ++cx) {
			std::string c1 = string_snprintf(":%d", col[cx]);
			cp.push_back(path + c1);
		}
	}
	curLine = 0;
}

void ComputeLoadContext::computeImpl(FitContext *fc)
{
	if (numColumns == 0) return;

	std::vector<int> &clc = Global->computeLoopIndex;
	if (clc.size() == 0) mxThrow("%s: must be used within a loop", name);
	int index = clc[clc.size()-1] - 1;  // innermost loop index

	if (index < curLine) {
		reopen();
		st->skip_line(); // header
		curLine = 0;
	}
	while (index > curLine) {
		st->skip_line();
		curLine += 1;
	}

	bool fail = false; // not sure how failure is reported
	try {
		fail = !st->read_line();
	} catch (...) {
		fail = true;
	}
	if (fail) mxThrow("%s: '%s' ran out of data at record %d",
			  name, path.c_str(), 1+index);

	Eigen::Map< Eigen::ArrayXi > col(columnPtr, numColumns);
	auto &cv = Global->checkpointValues;

	for (int cx=0, xx=0; cx < maxColumn; ++cx) {
		std::string c1;
		*st >> c1;
		if (cx == col[xx]-1) {
			cv[cpIndex + xx] = c1;
			if (++xx == numColumns) break;
		}
	}
	curLine += 1;
}

void ComputeLoadContext::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	MxRList dbg;
	dbg.add("loadCounter", Rf_ScalarInteger(loadCounter));
	slots->add("debug", dbg.asR());
}

void ComputeLoadMatrix::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);

	ProtectedSEXP RoriginalData(R_do_slot(rObj, Rf_install("originalDataIsIndexOne")));
	useOriginalData = Rf_asLogical(RoriginalData);

	ProtectedSEXP Rdata(R_do_slot(rObj, Rf_install("dest")));
	ProtectedSEXP Rpath(R_do_slot(rObj, Rf_install("path")));
	hasRowNames.resize(Rf_length(Rdata));
	if (useOriginalData) origData.resize(Rf_length(Rdata));
	ProtectedSEXP Rrownames(R_do_slot(rObj, Rf_install("row.names")));
	ProtectedSEXP Rcolnames(R_do_slot(rObj, Rf_install("col.names")));

	ProtectedSEXP Rmethod(R_do_slot(rObj, Rf_install("method")));
	const char *methodName = R_CHAR(STRING_ELT(Rmethod, 0));
	if (strEQ(methodName, "csv")) {
		loadMethod = LoadCSV;
		if (Rf_length(Rpath) != Rf_length(Rdata)) {
			mxThrow("%s: %d matrices to load but only %d paths",
				name, Rf_length(Rdata), Rf_length(Rpath));
		}
	} else if (strEQ(methodName, "data.frame")) {
		loadMethod = LoadDataFrame;
	} else {
		mxThrow("%s: unknown method '%s'", name, methodName);
	}

	int numElem = 0;
	for (int wx=0; wx < Rf_length(Rdata); ++wx) {
		if (isErrorRaised()) return;
		int objNum = ~INTEGER(Rdata)[wx];
		omxMatrix *m1 = globalState->matrixList[objNum];
		if (m1->hasPopulateSubstitutions()) {
			omxRaiseErrorf("%s: matrix '%s' has populate substitutions",
				       name, m1->name());
		}
		numElem += m1->numNonConstElements();
		EigenMatrixAdaptor Em1(m1);
		if (useOriginalData) origData[wx] = Em1;
		mat.push_back(m1);

		if (loadMethod != LoadDataFrame) {
			const char *p1 = R_CHAR(STRING_ELT(Rpath, wx));
			// convert to std::string to work around mini::csv constructor bug
			// https://github.com/shaovoon/minicsv/issues/8
			streams.push_back(new mini::csv::ifstream(std::string(p1)));
			mini::csv::ifstream &st = *streams[wx];
			st.set_delimiter(' ', "##");
			if (INTEGER(Rcolnames)[wx % Rf_length(Rcolnames)]) {
				st.skip_line();
			}
		}

		hasRowNames[wx] = INTEGER(Rrownames)[wx % Rf_length(Rrownames)];
		//mxLog("ld %s %s", d1->name, p1);
	}
	line = 1;

	if (loadMethod == LoadDataFrame) {
		ProtectedSEXP Robs(R_do_slot(rObj, Rf_install("observed")));
		observed = Robs;
		if (int(observed.size()) < numElem) {
			mxThrow("%s: provided observed data only has %d columns but %d requested",
				name, int(observed.size()), numElem);
		}
		Rcpp::CharacterVector obNames = observed.attr("names");
		for (int cx=0; cx < numElem; ++cx) {
			if (Rcpp::is<Rcpp::NumericVector>(observed[cx])) continue;

			mxThrow("%s: observed column %d (%s) is not type 'numeric'",
				name, 1+cx, Rcpp::as<const char *>(obNames[cx]));
		}
	}

}

ComputeLoadMatrix::~ComputeLoadMatrix()
{
	for (auto st : streams) delete st;
	streams.clear();
}

void ComputeLoadMatrix::computeImpl(FitContext *fc)
{
	std::vector<int> &clc = Global->computeLoopIndex;
	if (clc.size() == 0) mxThrow("%s: must be used within a loop", name);
	int index = clc[clc.size()-1];  // innermost loop index
	if (useOriginalData && index == 1) {
		for (int dx=0; dx < int(mat.size()); ++dx) {
			EigenMatrixAdaptor Em(mat[dx]);
			Em.derived() = origData[dx];
		}
		return;
	} else {
		index -= useOriginalData; // 1 == the first record
		switch (loadMethod) {
		case LoadCSV:
			loadFromCSV(fc, index);
			break;
		case LoadDataFrame:
			loadDataFrame(fc, index);
			break;
		default:
			mxThrow("%s: unknown load method %d", name, loadMethod);
		}
	}

	fc->state->invalidateCache();
	fc->state->connectToData();
	fc->state->omxInitialMatrixAlgebraCompute(fc);
	if (isErrorRaised()) mxThrow("%s", Global->getBads()); // ?still necessary?
}

struct clmStream {
	Rcpp::DataFrame &observed;
	const int row;
	int curCol;
	clmStream(Rcpp::DataFrame &_ob, int _row) : observed(_ob), row(_row) { curCol = 0; };

	void operator >> (double& val)
	{
		auto vec = observed[curCol];
		val = REAL(vec)[row];
		curCol += 1;
	}
};

void ComputeLoadMatrix::loadDataFrame(FitContext *fc, int index)
{
	if (observed.nrows() < index) {
		mxThrow("%s: index %d requested but observed data only has %d rows",
			name, index, observed.nrows());
	}

	clmStream st(observed, index - 1);
	for (int dx=0; dx < int(mat.size()); ++dx) {
		mat[dx]->loadFromStream(st);
	}
}

void ComputeLoadMatrix::loadFromCSV(FitContext *fc, int index)
{
	if (line > index) {
		mxThrow("%s: at line %d, cannot seek backwards to line %d",
			 name, line, index);
	}
	while (line < index) {
		for (int dx=0; dx < int(mat.size()); ++dx) {
			mini::csv::ifstream &st = *streams[dx];
			st.skip_line();
		}
		line += 1;
	}
	for (int dx=0; dx < int(mat.size()); ++dx) {
		mini::csv::ifstream &st = *streams[dx];
		if (!st.read_line()) {
			mxThrow("%s: ran out of data for matrix '%s'",
				 name, mat[dx]->name());
		}
		if (hasRowNames[dx]) {
			std::string rn;
			st >> rn;  // discard it
		}
		mat[dx]->loadFromStream(st);
	}
	line += 1;
}

void ComputeCheckpoint::initFromFrontend(omxState *globalState, SEXP rObj)
{
	super::initFromFrontend(globalState, rObj);
	badSEWarning = false;

	ProtectedSEXP Rappend(R_do_slot(rObj, Rf_install("append")));
	bool append = Rf_asLogical(Rappend);

	ProtectedSEXP Rheader(R_do_slot(rObj, Rf_install("header")));
	wroteHeader = !Rf_asLogical(Rheader);

	firstTime = true;
	numExtraCols = 0;

	path = 0;
	ProtectedSEXP Rpath(R_do_slot(rObj, Rf_install("path")));
	if (Rf_length(Rpath) == 1) {
		path = R_CHAR(STRING_ELT(Rpath, 0));
		ofs.open(path, append? std::ofstream::app : std::ofstream::trunc);
		if (!ofs.is_open()) {
			mxThrow("Failed to open '%s' for writing", path);
		}
	}

	ProtectedSEXP RtoReturn(R_do_slot(rObj, Rf_install("toReturn")));
	toReturn = Rf_asLogical(RtoReturn);

	ProtectedSEXP Rpar(R_do_slot(rObj, Rf_install("parameters")));
	inclPar = Rf_asLogical(Rpar);

	ProtectedSEXP RloopInd(R_do_slot(rObj, Rf_install("loopIndices")));
	inclLoop = Rf_asLogical(RloopInd);

	ProtectedSEXP Rfit(R_do_slot(rObj, Rf_install("fit")));
	inclFit = Rf_asLogical(Rfit);

	ProtectedSEXP Rcounters(R_do_slot(rObj, Rf_install("counters")));
	inclCounters = Rf_asLogical(Rcounters);

	ProtectedSEXP Rstatus(R_do_slot(rObj, Rf_install("status")));
	inclStatus = Rf_asLogical(Rstatus);

	ProtectedSEXP Rse(R_do_slot(rObj, Rf_install("standardErrors")));
	inclSEs = Rf_asLogical(Rse);

	ProtectedSEXP Rgradient(R_do_slot(rObj, Rf_install("gradient")));
	inclGradient = Rf_asLogical(Rgradient);

	ProtectedSEXP Rvcov(R_do_slot(rObj, Rf_install("vcov")));
	inclVcov = Rf_asLogical(Rvcov);

	ProtectedSEXP Rwhat(R_do_slot(rObj, Rf_install("what")));
	for (int wx=0; wx < Rf_length(Rwhat); ++wx) {
		if (isErrorRaised()) return;
		int objNum = INTEGER(Rwhat)[wx];
		omxMatrix *algebra = globalState->algebraList[objNum];
		if (algebra->fitFunction) {
			omxCompleteFitFunction(algebra);
		}
		algebras.push_back(algebra);
	}

	if (inclCounters) {
		colnames.push_back("OpenMxEvals");
		colnames.push_back("iterations");
	}
	colnames.push_back("timestamp");
	if (inclLoop) {
		auto &clc = Global->computeLoopContext;
		for (int lx=0; lx < int(clc.size()); ++lx) {
			colnames.push_back(string_snprintf("%s%d", clc[lx], 1+lx));
		}
	}

	std::vector< omxFreeVar* > &vars = Global->findVarGroup(FREEVARGROUP_ALL)->vars;
	numParam = vars.size();

	if (inclPar) {
		for(int j = 0; j < numParam; j++) {
			colnames.push_back(vars[j]->name);
		}
	}

	if (inclFit) {
		colnames.push_back("objective");
		colnames.push_back("fitUnits");
	}
	if (inclStatus) colnames.push_back("statusCode");
	if (inclSEs) {
		for(int j = 0; j < numParam; j++) {
			std::string c1 = vars[j]->name;
			c1 += "SE";
			colnames.push_back(c1);
		}
	}
	if (inclGradient) {
		for (auto v1 : vars) {
			std::string c1 = v1->name;
			c1 += "Grad";
			colnames.push_back(c1);
		}
	}
	if (inclVcov) {
		for (int cx=0; cx < numParam; ++cx) {
			for (int rx=cx; rx < numParam; ++rx) {
				std::string c1 = "V";
				c1 += vars[rx]->name;
				c1 += ":";
				c1 += vars[cx]->name;
				colnames.push_back(c1);
			}
		}
	}

	numAlgebraEnt = 0;
	for (auto &mat : algebras) {
		for (int cx=0; cx < mat->cols; ++cx) {
			for (int rx=0; rx < mat->rows; ++rx) {
				colnames.push_back(string_snprintf("%s[%d,%d]", mat->name(), 1+rx, 1+cx));
			}
		}
		numAlgebraEnt += mat->cols * mat->rows;
	}
	// TODO: confidence intervals, hessian, constraint algebras
	// what does Eric Schmidt include in per-voxel output?
	// remove old checkpoint code?
	numSnaps = 0;
}

void ComputeCheckpoint::computeImpl(FitContext *fc)
{
	// if ((timePerCheckpoint && timePerCheckpoint <= now - lastCheckpoint) ||
	//     (iterPerCheckpoint && iterPerCheckpoint <= fc->iterations - lastIterations) ||
	//     (evalsPerCheckpoint && evalsPerCheckpoint <= curEval - lastEvaluation)) {

	snap s1;
	s1.evaluations = fc->getGlobalComputeCount();
	s1.iterations = fc->iterations;
	s1.timestamp = time(0);
	if (inclLoop) s1.computeLoopIndex = Global->computeLoopIndex;
	if (inclPar) {
		Eigen::Map< Eigen::VectorXd > Eest(fc->est, fc->numParam);
		s1.est = Eest;
	}
	s1.fit = fc->fit;
	s1.fitUnits = fc->fitUnits;
	s1.inform = fc->wrapInform();
	if (inclSEs) {
		if (fc->stderrs.size() && fc->stderrs.size() != int(fc->numParam)) {
			if (!badSEWarning) {
				Rf_warning("%s: there are %d standard errors but %d parameters",
					   name, fc->stderrs.size(), int(fc->numParam));
				badSEWarning = true;
			}
		}
		s1.stderrs = fc->stderrs;
	}
	if (inclGradient) s1.gradient = fc->gradZ;
	if (inclVcov && fc->vcov.rows() == numParam) {
		s1.vcov.resize(triangleLoc1(numParam));
		int lx=0;
		for (int cx=0; cx < numParam; ++cx) {
			for (int rx=cx; rx < numParam; ++rx) {
				s1.vcov[lx++] = fc->vcov(rx, cx);
			}
		}
	}
	s1.algebraEnt.resize(numAlgebraEnt);

	{int xx=0;
	for (auto &mat : algebras) {
		for (int cx=0; cx < mat->cols; ++cx) {
			for (int rx=0; rx < mat->rows; ++rx) {
				EigenMatrixAdaptor eM(mat);
				s1.algebraEnt[xx] = eM(rx, cx);
				xx += 1;
			}
		}
	}}

	if (firstTime) {
		auto &xcn = Global->checkpointColnames;
		numExtraCols = xcn.size();
		colnames.insert(colnames.end(), xcn.begin(), xcn.end());
		firstTime = false;
	}
	s1.extra = Global->checkpointValues;

	if (ofs.is_open()) {
		const int digits = std::numeric_limits<double>::digits10 + 1;
		if (!wroteHeader) {
			bool first = true;
			for (auto &cn : colnames) {
				if (first) { first=false; }
				else       { ofs << '\t'; }
				ofs << cn;
			}
			ofs << std::endl;
			wroteHeader = true;
		}

		bool first = true;
		if (inclCounters) {
			if (first) { first=false; }
			else       { ofs << '\t'; }
			ofs << s1.evaluations;
			ofs << '\t' << s1.iterations;
		}
		if (1) {
			if (first) { first=false; }
			else       { ofs << '\t'; }
			const int timeBufSize = 32;
			char timeBuf[timeBufSize];
			struct tm *nowTime = localtime(&s1.timestamp);
			strftime(timeBuf, timeBufSize, "%b %d %Y %I:%M:%S %p", nowTime);
			ofs << timeBuf;
		}
		if (inclLoop) {
			auto &clc = s1.computeLoopIndex;
			for (int lx=0; lx < int(clc.size()); ++lx) {
				ofs << '\t' << clc[lx];
			}
		}
		if (inclPar) {
			for (int x1=0; x1 < int(s1.est.size()); ++x1) {
				ofs << '\t' << std::setprecision(digits) << s1.est[x1];
			}
		}
		if (inclFit) {
			ofs << '\t' << std::setprecision(digits) << s1.fit;
			ofs << '\t' << fitUnitsToName(s1.fitUnits);
		}
		if (inclStatus) {
			if (s1.inform == NA_INTEGER) {
				ofs << '\t' << "NA";
			} else {
				ofs << '\t' << statusCodeLabels[s1.inform - 1];
			}
		}
		if (inclSEs) {
			if (s1.stderrs.size()) {
				for (int x1=0; x1 < int(s1.est.size()); ++x1) {
					ofs << '\t' << std::setprecision(digits) << s1.stderrs[x1];
				}
			} else {
				for (int x1=0; x1 < int(s1.est.size()); ++x1) {
					ofs << '\t' << NA_REAL;
				}
			}
		}
		if (inclGradient) {
			if (s1.gradient.size()) {
				for (int x1=0; x1 < int(s1.est.size()); ++x1) {
					ofs << '\t' << std::setprecision(digits) << s1.gradient[x1];
				}
			} else {
				for (int x1=0; x1 < int(s1.est.size()); ++x1) {
					ofs << '\t' << NA_REAL;
				}
			}
		}
		if (inclVcov) {
			int numVcov = triangleLoc1(numParam);
			if (s1.vcov.size()) {
				for (int x1=0; x1 < numVcov; ++x1) {
					ofs << '\t' << std::setprecision(digits) << s1.vcov[x1];
				}
			} else {
				for (int x1=0; x1 < numVcov; ++x1) {
					ofs << '\t' << NA_REAL;
				}
			}
		}
		for (int x1=0; x1 < int(s1.algebraEnt.size()); ++x1) {
			ofs << '\t' << std::setprecision(digits) << s1.algebraEnt[x1];
		}
		for (auto &x1 : s1.extra) ofs << '\t' << x1;
		ofs << std::endl;
		ofs.flush();
	}

	if (toReturn) {
		snaps.push_front(s1);
		numSnaps += 1;
	}
}

void ComputeCheckpoint::reportResults(FitContext *fc, MxRList *slots, MxRList *)
{
	const bool debug = false;
	if (ofs.is_open()) {
		ofs.close();
	}
	if (!toReturn) return;

	snaps.reverse();

	SEXP log;
	Rf_protect(log = Rf_allocVector(VECSXP, colnames.size()));
	int curCol=0;
	if (inclCounters) {
		{
			SEXP col = Rf_allocVector(INTSXP, numSnaps);
			if (debug) mxLog("log[%d] = %s (evaluations)", curCol, colnames[curCol].c_str());
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = INTEGER(col);
			int sx=0;
			for (auto &s1 : snaps) v[sx++] = s1.evaluations;
		}
		{
			SEXP col = Rf_allocVector(INTSXP, numSnaps);
			if (debug) mxLog("log[%d] = %s (iterations)", curCol, colnames[curCol].c_str());
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = INTEGER(col);
			int sx=0;
			for (auto &s1 : snaps) v[sx++] = s1.iterations;
		}
	}
	{
		SEXP POSIXct;
		Rf_protect(POSIXct = Rf_allocVector(STRSXP, 2));
		SET_STRING_ELT(POSIXct, 0, Rf_mkChar("POSIXct"));
		SET_STRING_ELT(POSIXct, 1, Rf_mkChar("POSIXt"));
		SEXP col = Rf_allocVector(REALSXP, numSnaps);
		Rf_setAttrib(col, R_ClassSymbol, POSIXct);
		if (debug) mxLog("log[%d] = %s (timestamp)", curCol, colnames[curCol].c_str());
		SET_VECTOR_ELT(log, curCol++, col);
		auto *v = REAL(col);
		int sx=0;
		for (auto &s1 : snaps) v[sx++] = s1.timestamp;
	}
	if (inclLoop) {
		auto &clc = snaps.front().computeLoopIndex;
		for (int lx=0; lx < int(clc.size()); ++lx) {
			SEXP col = Rf_allocVector(INTSXP, numSnaps);
			if (debug) mxLog("log[%d] = %s (loop index %d)", curCol, colnames[curCol].c_str(), lx);
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = INTEGER(col);
			int sx=0;
			for (auto &s1 : snaps) v[sx++] = s1.computeLoopIndex[lx];
		}
	}
	if (inclPar) {
		auto numEst = int(snaps.front().est.size());
		for (int x1=0; x1 < numEst; ++x1) {
			SEXP col = Rf_allocVector(REALSXP, numSnaps);
			if (debug) mxLog("log[%d] = %s (parameter %d/%d)", curCol, colnames[curCol].c_str(),
			      x1, numEst);
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = REAL(col);
			int sx=0;
			for (auto &s1 : snaps) v[sx++] = s1.est[x1];
		}
	}
	if (inclFit) {
		SEXP col1 = Rf_allocVector(REALSXP, numSnaps);
		if (debug) mxLog("log[%d] = %s (fit)", curCol, colnames[curCol].c_str());
		SET_VECTOR_ELT(log, curCol++, col1);
		SEXP col2 = makeFitUnitsFactor(Rf_allocVector(INTSXP, numSnaps));
		SET_VECTOR_ELT(log, curCol++, col2);
		auto *v = REAL(col1);
		auto *vu = INTEGER(col2);
		int sx=0;
		for (auto &s1 : snaps) {
			v[sx] = s1.fit;
			if (s1.fitUnits == FIT_UNITS_UNINITIALIZED) {
				vu[sx] = NA_INTEGER;
			} else {
				vu[sx] = s1.fitUnits;
			}
			sx += 1;
		}
	}
	if (inclStatus) {
		SEXP col = allocInformVector(numSnaps);
		if (debug) mxLog("log[%d] = %s (inform)", curCol, colnames[curCol].c_str());
		SET_VECTOR_ELT(log, curCol++, col);
		auto *v = INTEGER(col);
		int sx=0;
		for (auto &s1 : snaps) v[sx++] = s1.inform;
	}
	if (inclSEs) {
		auto numEst = int(snaps.front().est.size());
		for (int x1=0; x1 < numEst; ++x1) {
			SEXP col = Rf_allocVector(REALSXP, numSnaps);
			if (debug) mxLog("log[%d] = %s (SE %d/%d)", curCol, colnames[curCol].c_str(), x1, numEst);
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = REAL(col);
			int sx=0;
			for (auto &s1 : snaps) {
				if (s1.stderrs.size()) {
					v[sx++] = s1.stderrs[x1];
				} else {
					v[sx++] = NA_REAL;
				}
			}
		}
	}
	if (inclGradient) {
		auto numEst = int(snaps.front().est.size());
		for (int x1=0; x1 < numEst; ++x1) {
			SEXP col = Rf_allocVector(REALSXP, numSnaps);
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = REAL(col);
			int sx=0;
			for (auto &s1 : snaps) {
				v[sx++] = s1.gradient.size()? s1.gradient[x1] : NA_REAL;
			}
		}
	}
	if (inclVcov) {
		int numVcov = triangleLoc1(numParam);
		for (int x1=0; x1 < numVcov; ++x1) {
			SEXP col = Rf_allocVector(REALSXP, numSnaps);
			SET_VECTOR_ELT(log, curCol++, col);
			auto *v = REAL(col);
			int sx=0;
			for (auto &s1 : snaps) {
				v[sx++] = s1.vcov.size()? s1.vcov[x1] : NA_REAL;
			}
		}
	}
	for (int x1=0; x1 < numAlgebraEnt; ++x1) {
		SEXP col = Rf_allocVector(REALSXP, numSnaps);
		if (debug) mxLog("log[%d] = %s (alge %d/%d)", curCol, colnames[curCol].c_str(), x1, numAlgebraEnt);
		SET_VECTOR_ELT(log, curCol++, col);
		auto *v = REAL(col);
		int sx=0;
		for (auto &s1 : snaps) v[sx++] = s1.algebraEnt[x1];
	}
	auto &xcn = Global->checkpointColnames;
	if (xcn.size() != numExtraCols) {
		mxThrow("%s: xcn.size() != numExtraCols; %d != %d",
			name, int(xcn.size()), int(numExtraCols));
	}
	for (int x1=0; x1 < int(xcn.size()); ++x1) {
		SEXP col = Rf_allocVector(STRSXP, numSnaps);
		if (debug) mxLog("log[%d] = %s (extra %d/%d)", curCol, colnames[curCol].c_str(), x1, int(xcn.size()));
		SET_VECTOR_ELT(log, curCol++, col);
		int sx=0;
		for (auto &s1 : snaps) SET_STRING_ELT(col, sx++, Rf_mkChar(s1.extra[x1].c_str()));
	}

	markAsDataFrame(log, numSnaps);
	SEXP logCols;
	Rf_protect(logCols = Rf_allocVector(STRSXP, colnames.size()));
	Rf_setAttrib(log, R_NamesSymbol, logCols);
	for (int cx=0; cx < int(colnames.size()); ++cx) {
		SET_STRING_ELT(logCols, cx, Rf_mkChar(colnames[cx].c_str()));
	}

	slots->add("log", log);
}
