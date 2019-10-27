/*
 *  Copyright 2007-2019 by the individuals mentioned in the source code history
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
 *
 */

/***********************************************************
 *
 *  omxMatrix.h
 *
 *  Created: Timothy R. Brick 	Date: 2008-11-13 12:33:06
 *
 *	Contains header information for the omxMatrix class
 *   omxDataMatrices hold necessary information to simplify
 * 	dealings between the OpenMX back end and BLAS.
 *
 **********************************************************/

#ifndef _OMXMATRIX_H_
#define _OMXMATRIX_H_

#include <R_ext/Arith.h>

#include "omxDefines.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include "minicsv.h"

struct populateLocation {
	int from;
	int srcRow, srcCol;
	int destRow, destCol;

	populateLocation() {};
	populateLocation(int _from, int _srcRow, int _srcCol, int _destRow, int _destCol)
	: from(_from), srcRow(_srcRow), srcCol(_srcCol), destRow(_destRow), destCol(_destCol) {};
	void transpose() { std::swap(destRow, destCol); }
};

//
// Ways that values can be loaded into plain matrices (not algebras)
// and how to ensure that non-zero values are loaded into these cells,
//
// + Free parameters, omxState->setFakeParam(estSave);
// + Populated labels, omxMatrix->markPopulatedEntries();
// + Definition variables, omxData->loadFakeData(currentState, 1.0);
//

class omxMatrix {
	std::vector< populateLocation > populate;  // For inclusion of values from other matrices
	// Note: some overlap with FreeVarGroup::cacheDependencies
	bool dependsOnParametersCache;    // Ignores free variable groups
	bool dependsOnDefVarCache;
	int joinKey;
	class omxExpectation *joinModel;
	int shape;
	bool allocationLock;   // whether the data can be moved
 public:
 	omxMatrix() : dependsOnParametersCache(false), dependsOnDefVarCache(false), joinKey(-1),
								joinModel(0), shape(0), allocationLock(false),
								freeRownames(false), freeColnames(false)
		{};
	void setDependsOnParameters() { dependsOnParametersCache = true; };
	void setDependsOnDefinitionVariables() { dependsOnDefVarCache = true; };
	bool dependsOnParameters() const { return dependsOnParametersCache; };
	bool dependsOnDefinitionVariables() const { return dependsOnDefVarCache; };
	bool hasPopulateSubstitutions() const { return populate.size(); };
	void addPopulate(omxMatrix *from, int srcRow, int srcCol, int destRow, int destCol);
	void transposePopulate();
	void lockAllocation() { allocationLock = true; }
	void setData(double *ptr);
	void setJoinInfo(SEXP Rmodel, SEXP Rkey);
	void omxProcessMatrixPopulationList(SEXP matStruct);
	void omxPopulateSubstitutions(int want, FitContext *fc);
	void markPopulatedEntries();
	int getJoinKey() const { return joinKey; }
	class omxExpectation *getJoinModel() const { return joinModel; }
										//TODO: Improve encapsulation
/* Actually Useful Members */
	int rows, cols;						// Matrix size  (specifically, its leading edge)
	double* data;						// Actual Data Pointer
	unsigned short colMajor;			// used for quick transpose
	unsigned short hasMatrixNumber;		// is this object in the matrix or algebra arrays?
	int matrixNumber;					// the offset into the matrices or algebras arrays

	SEXP owner;	// The R object owning data or NULL if we own it.

/* For BLAS Multiplication Speedup */ 	// TODO: Replace some of these with inlines or macros.
	const char* majority;				// Filled by compute(), included for speed
	const char* minority;				// Filled by compute(), included for speed
	int leading;						// Leading edge; depends on original majority
	int lagging;						// Non-leading edge.

/* Curent State */
	omxState* currentState;				// Optimizer State
	unsigned cleanVersion;
	unsigned version;

/* For Algebra Functions */				// At most, one of these may be non-NULL.
	bool canDiscard();
	omxAlgebra* algebra;				// If it's not an algebra, this is NULL.
	omxFitFunction* fitFunction;		// If it's not a fit function, this is NULL.

	std::string nameStr;
	const char *name() const { return nameStr.c_str(); }

	bool freeRownames, freeColnames;
	std::vector<const char *> rownames;
	std::vector<const char *> colnames;
	int lookupColumnByName(const char *target);

	friend void omxCopyMatrix(omxMatrix *dest, omxMatrix *src);  // turn into method later TODO
	void take(omxMatrix *orig);

	void unshareMemoryWithR();
	void loadDimnames(SEXP dimnames);
	const char *getType() const {
		const char *what = "matrix";
		if (algebra) what = "algebra";
		if (fitFunction) what = "fitfunction";
		return what;
	}
	void copyAttr(omxMatrix *src);
	bool isSimple() const { return !algebra && !fitFunction && populate.size()==0; };
	bool isAlgebra() const { return algebra != 0; }
	int numNonConstElements() const;
	template <typename T> void loadFromStream(T &st);
	int size() const { return rows * cols; }
	SEXP asR();
	bool isValidElem(int row, int col)
	{ return row >= 0 && col >= 0 && row < rows && col < cols; };
};

void omxEnsureColumnMajor(omxMatrix *mat);

inline double *omxMatrixDataColumnMajor(omxMatrix *mat)
{
	omxEnsureColumnMajor(mat);
	return mat->data;
}

// NOTE: These Eigen wrapper are not thread safe when mixed with omxCopyMatrix

struct EigenMatrixAdaptor : Eigen::Map< Eigen::MatrixXd > {
	EigenMatrixAdaptor(omxMatrix *mat) :
	  Eigen::Map< Eigen::MatrixXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols) {}
};

struct EigenMatrixAdaptor0 : Eigen::Map< Eigen::MatrixXd > {
	EigenMatrixAdaptor0(omxMatrix *mat) :
	Eigen::Map< Eigen::MatrixXd >(mat? omxMatrixDataColumnMajor(mat) : 0,
				      mat? mat->rows : 0,
				      mat? mat->cols : 0) {}
};

struct EigenArrayAdaptor : Eigen::Map< Eigen::ArrayXXd > {
	EigenArrayAdaptor(omxMatrix *mat) :
	  Eigen::Map< Eigen::ArrayXXd >(omxMatrixDataColumnMajor(mat), mat->rows, mat->cols) {}
};

struct EigenVectorAdaptor : Eigen::Map< Eigen::VectorXd > {
	EigenVectorAdaptor(omxMatrix *mat) :
	  Eigen::Map< Eigen::VectorXd >(mat->data, mat->rows * mat->cols) {}
};

template <typename Scalar>
struct EigenStdVectorAdaptor : Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> > {
	EigenStdVectorAdaptor(std::vector<Scalar> &vec) :
	Eigen::Map< Eigen::Matrix<Scalar, Eigen::Dynamic, 1> >(vec.data(), vec.size()) {};
};

// If you call these functions directly then you need to free the memory with omxFreeMatrix.
// If you obtain a matrix from omxNewMatrixFromSlot then you must NOT free it.
omxMatrix* omxInitMatrix(int nrows, int ncols, unsigned short colMajor, omxState* os);
inline omxMatrix* omxInitMatrix(int nrows, int ncols, omxState* os)
{ return omxInitMatrix(nrows, ncols, 1, os); }
omxMatrix *omxCreateCopyOfMatrix(omxMatrix *orig, omxState *os);

	void omxFreeMatrix(omxMatrix* om);						// Ditto, traversing argument trees

/* Matrix Creation Functions */
omxMatrix* omxNewMatrixFromRPrimitive0(SEXP rObject, omxState* state, 
       unsigned short hasMatrixNumber, int matrixNumber);
	omxMatrix* omxNewMatrixFromRPrimitive(SEXP rObject, omxState *state,
	unsigned short hasMatrixNumber, int matrixNumber); 							// Create an omxMatrix from an R object
	omxMatrix* omxNewIdentityMatrix(int nrows, omxState* state);				// Creates an Identity Matrix of a given size
	extern omxMatrix* omxMatrixLookupFromState1(SEXP matrix, omxState* os);	// Create a matrix/algebra from a matrix pointer

	omxMatrix* omxDuplicateMatrix(omxMatrix* src, omxState* newState);
	SEXP omxExportMatrix(omxMatrix *om);

/* Getters 'n Setters (static functions declared below) */
	// static OMXINLINE double omxMatrixElement(omxMatrix *om, int row, int col);
	// static OMXINLINE double omxVectorElement(omxMatrix *om, int index);
	// static OMXINLINE void omxSetMatrixElement(omxMatrix *om, int row, int col, double value);
	// static OMXINLINE void omxSetVectorElement(omxMatrix *om, int index, double value);

	double* omxLocationOfMatrixElement(omxMatrix *om, int row, int col);
	void omxMarkDirty(omxMatrix *om);
void omxMarkClean(omxMatrix *om);

/* Matrix Modification Functions */
	void omxZeroByZeroMatrix(omxMatrix *source);
void omxResizeMatrix(omxMatrix *source, int nrows, int ncols);
	omxMatrix* omxFillMatrixFromRPrimitive(omxMatrix* om, SEXP rObject, omxState *state,
		unsigned short hasMatrixNumber, int matrixNumber); 								// Populate an omxMatrix from an R object
	void omxTransposeMatrix(omxMatrix *mat);												// Transpose a matrix in place.
	void omxToggleRowColumnMajor(omxMatrix *mat);										// Transform row-major into col-major and vice versa 

/* Function wrappers that switch based on inclusion of algebras */
	void omxPrint(omxMatrix *source, const char* d);

void omxRecompute(omxMatrix *matrix, FitContext *fc);
void CheckAST(omxMatrix *matrix, FitContext *fc);

void omxRemoveElements(omxMatrix *om, int removed[]);
void omxRemoveRowsAndColumns(omxMatrix *om, int rowsRemoved[], int colsRemoved[]);

/* Matrix-Internal Helper functions */
	void omxMatrixLeadingLagging(omxMatrix *matrix);
void omxPrintMatrix(omxMatrix *source, const char* header);  // deprecated, use omxPrint

/* OMXINLINE functions and helper functions */

void setMatrixError(omxMatrix *om, int row, int col, int numrow, int numcol);
void setVectorError(int index, int numrow, int numcol);
void matrixElementError(int row, int col, omxMatrix *om);
void vectorElementError(int index, int numrow, int numcol);

bool omxNeedsUpdate(omxMatrix *matrix);

OMXINLINE static bool omxMatrixIsDirty(omxMatrix *om) { return om->cleanVersion != om->version; }
OMXINLINE static bool omxMatrixIsClean(omxMatrix *om) { return om->cleanVersion == om->version; }
OMXINLINE static unsigned omxGetMatrixVersion(omxMatrix *om) { return om->version; }

static OMXINLINE int omxIsMatrix(omxMatrix *mat) {
    return (mat->algebra == NULL && mat->fitFunction == NULL);
}

/* BLAS Wrappers */

static OMXINLINE void omxSetMatrixElement(omxMatrix *om, int row, int col, double value) {
	if (!om->isValidElem(row, col)) {
		setMatrixError(om, row + 1, col + 1, om->rows, om->cols);
		return;
	}
	int index = 0;
	if(om->colMajor) {
		index = col * om->rows + row;
	} else {
		index = row * om->cols + col;
	}
	om->data[index] = value;
}

static OMXINLINE void omxAccumulateMatrixElement(omxMatrix *om, int row, int col, double value) {
	if (!om->isValidElem(row, col)) {
                setMatrixError(om, row + 1, col + 1, om->rows, om->cols);
                return;
        }
        int index = 0;
        if(om->colMajor) {
                index = col * om->rows + row;
        } else {
                index = row * om->cols + col;
        }
        om->data[index] += value;
}

static OMXINLINE double omxMatrixElement(omxMatrix *om, int row, int col) {
	int index = 0;
	if (!om->isValidElem(row, col)) {
		matrixElementError(row + 1, col + 1, om);
        return (NA_REAL);
	}
	if(om->colMajor) {
		index = col * om->rows + row;
	} else {
		index = row * om->cols + col;
	}
	return om->data[index];
}

static OMXINLINE double *omxMatrixColumn(omxMatrix *om, int col) {
  if (!om->colMajor) mxThrow("omxMatrixColumn requires colMajor order");
  if (col < 0 || col >= om->cols) mxThrow("omxMatrixColumn(%d) but only %d columns", col, om->cols);
  return om->data + col * om->rows;
}

static OMXINLINE void omxAccumulateVectorElement(omxMatrix *om, int index, double value) {
	if (index < 0 || index >= (om->rows * om->cols)) {
		setVectorError(index + 1, om->rows, om->cols);
		return;
	} else {
		om->data[index] += value;
    }
}

static OMXINLINE void omxSetVectorElement(omxMatrix *om, int index, double value) {
	if (index < 0 || index >= (om->rows * om->cols)) {
		setVectorError(index + 1, om->rows, om->cols);
		return;
	} else {
		om->data[index] = value;
    }
}

static OMXINLINE double omxVectorElement(omxMatrix *om, int index) {
	if (index < 0 || index >= (om->rows * om->cols)) {
		vectorElementError(index + 1, om->rows, om->cols);
        return (NA_REAL);
	} else {
		return om->data[index];
    }
}

static OMXINLINE void omxUnsafeSetVectorElement(omxMatrix *om, int index, double value) {
	om->data[index] = value;
}

static OMXINLINE double omxUnsafeVectorElement(omxMatrix *om, int index) {
	return om->data[index];
}

// In the process of debugging some other problem, I noticed that
// DDEBUGMX_ALGEBRA causes models/passing/SimpleAlgebraCIs.R to fail.
//
// I used GIT to bisect the problem. It seems b1e8136b758c838 broke
// DDEBUGMX_ALGEBRA. In this commit, I changed the Fortran argument
// passing from character to integer at request of some CRAN maintainers.
// However, the change also appears to have broken something subtle in
// Fortran/C communication.
//
// I am really not sure how to identify the underlying issue.
//
// Acting on the belief that it's better to have correct code than fast
// code, I decided to replace omxBLAS.f with equivalent Eigen code.
// omxBLAS was (supposedly) fast because it omitted the code to
// properly handle NaN (in some cases).

static OMXINLINE void omxDGEMM(unsigned short int transposeA, unsigned short int transposeB,		// result <- alpha * A %*% B + beta * C
				double alpha, omxMatrix* a, omxMatrix *b, double beta, omxMatrix* result) {

	EigenMatrixAdaptor eA(a);
	EigenMatrixAdaptor eB(b);
	EigenMatrixAdaptor eResult(result);
	Eigen::MatrixXd ccopy;

	if (beta) {
		ccopy = eResult * beta;
	}
	if (!transposeA && !transposeB) {
		eResult.derived() = alpha * eA * eB;
	} else if (transposeA && !transposeB) {
		eResult.derived() = alpha * eA.transpose() * eB;
	} else if (!transposeA && transposeB) {
		eResult.derived() = alpha * eA * eB.transpose();
	} else {
		eResult.derived() = alpha * eA.transpose() * eB.transpose();
	}
	if (beta) {
		eResult.derived() += ccopy;
	}
	result->colMajor = true;
	omxMatrixLeadingLagging(result);
}

static OMXINLINE void omxDGEMV(bool transpose, double alpha, omxMatrix* mat,	// result <- alpha * A %*% B + beta * C
				omxMatrix* vec, double beta, omxMatrix*result) {							// where B is treated as a vector
	EigenMatrixAdaptor eMat(mat);
	EigenVectorAdaptor eVec(vec);
	EigenMatrixAdaptor eResult(result);
	Eigen::VectorXd vcopy;

	//mxLog("%dx%d = %dx%d * %d (t=%d)", eResult.rows(), eResult.cols(),
	//eMat.rows(), eMat.cols(), eVec.rows(), transpose);

	if (beta) {
		vcopy = eResult * beta;
	}
	int rows = transpose? eMat.cols() : eMat.rows();
	if (eResult.cols() == rows) {
		if (transpose) {
			eResult.derived() = (alpha * eMat.transpose() * eVec).transpose();
		} else {
			eResult.derived() = (alpha * eMat * eVec).transpose();
		}
		if (beta) {
			eResult.derived() += vcopy.transpose();
		}
	} else {
		if (transpose) {
			eResult.derived() = alpha * eMat.transpose() * eVec;
		} else {
			eResult.derived() = alpha * eMat * eVec;
		}
		if (beta) {
			eResult.derived() += vcopy;
		}
	}
	result->colMajor = true;
	omxMatrixLeadingLagging(result);
}

static OMXINLINE void omxDSYMV(double alpha, omxMatrix* mat,            // result <- alpha * A %*% B + beta * C
				omxMatrix* vec, double beta, omxMatrix* result) {       // only A is symmetric, and B is a vector
	if(OMX_DEBUG) {
		int nVecEl = vec->rows * vec->cols;
		// mxLog("DSYMV: %c, %d, %f, 0x%x, %d, 0x%x, %d, %f, 0x%x, %d\n", u, (mat->cols),alpha, mat->data, (mat->leading), 
	                    // vec->data, onei, beta, result->data, onei); //:::DEBUG:::
		if(mat->cols != nVecEl) {
			mxThrow("Mismatch in symmetric vector/matrix multiply: %s (%d x %d) * (%d x 1).\n", "symmetric", mat->rows, mat->cols, nVecEl); // :::DEBUG:::
		}
	}

	EigenMatrixAdaptor Emat(mat);
	EigenVectorAdaptor Evec(vec);
	EigenVectorAdaptor Eresult(result);
	
	Eresult.derived() = alpha * (Emat.selfadjointView<Eigen::Upper>() * Evec).array() + beta;
}

static OMXINLINE void omxDSYMM(unsigned short int symmOnLeft, omxMatrix* symmetric,
			       omxMatrix *other, omxMatrix* result) {
	EigenMatrixAdaptor Es(symmetric);
	EigenMatrixAdaptor Eo(other);
	EigenMatrixAdaptor Eresult(result);

	if (symmOnLeft) {
		Eresult.derived() = Es.template selfadjointView<Eigen::Upper>() * Eo;
	} else {
		Eresult.derived() = Eo * Es.template selfadjointView<Eigen::Upper>();
	}
}

static OMXINLINE void omxDAXPY(double alpha, omxMatrix* lhs, omxMatrix* rhs)
{
	EigenVectorAdaptor El(lhs);
	EigenVectorAdaptor Er(rhs);
	Er += alpha * El;
}

static OMXINLINE double omxDDOT(omxMatrix* lhs, omxMatrix* rhs)
{
	EigenVectorAdaptor El(lhs);
	EigenVectorAdaptor Er(rhs);
	return El.transpose() * Er;
}

void omxShallowInverse(int numIters, omxMatrix* A, omxMatrix* Z, omxMatrix* Ax, omxMatrix* I );

double omxMaxAbsDiff(omxMatrix *m1, omxMatrix *m2);

bool thresholdsIncreasing(omxMatrix* om, int column, int count, FitContext *fc);

void omxMatrixHorizCat(omxMatrix** matList, int numArgs, omxMatrix* result);

void omxMatrixVertCat(omxMatrix** matList, int numArgs, omxMatrix* result);

void omxMatrixTrace(omxMatrix** matList, int numArgs, omxMatrix* result);

void expm_eigen(int n, double *rz, double *out);
void logm_eigen(int n, double *rz, double *out);

template <typename T>
std::string mxStringifyMatrix(const char *name, const Eigen::DenseBase<T> &mat, std::string &xtra,
			      bool debug=false)
{
	std::string buf;

	if (!debug && mat.rows() * mat.cols() > 1000) {
		buf = string_snprintf("%s is too large to print # %dx%d\n",
				name, mat.rows(), mat.cols());
		return buf;
	}

	bool transpose = false; // maybe easier for debugging use: mat.rows() > mat.cols();
	buf += string_snprintf("%s = %s matrix(c(    # %dx%d",
			       name, transpose? "t(" : "", mat.rows(), mat.cols());

	bool first=true;
	int rr = mat.rows();
	int cc = mat.cols();
	if (transpose) std::swap(rr,cc);
	if (!mat.derived().data()) {
		buf += "\nNULL";
	} else {
		for(int j = 0; j < rr; j++) {
			buf += "\n";
			for(int k = 0; k < cc; k++) {
				if (first) first=false;
				else buf += ",";
				double val;
				if (transpose) {
					val = mat(k,j);
				} else {
					val = mat(j,k);
				}
				buf += string_snprintf(" %3.6g", val);
			}
		}
	}

	int rows = mat.rows();
	int cols = mat.cols();
	if (transpose) std::swap(rows, cols);
	buf += string_snprintf("), byrow=TRUE, nrow=%d, ncol=%d",
			       rows, cols);
	buf += xtra;
	buf += ")";
	if (transpose) buf += ")";
	buf += "\n";
	return buf;
}

template <typename T>
void mxPrintMatX(const char *name, const Eigen::DenseBase<T> &mat, std::string &xtra)
{
	std::string buf = mxStringifyMatrix(name, mat, xtra);
	mxLogBig(buf);
}

template <typename T>
void mxPrintMat(const char *name, const Eigen::DenseBase<T> &mat)
{
	std::string xtra;
	mxPrintMatX(name, mat, xtra);
}

template <typename T1, typename T2, typename T3>
void computeMeanCov(const Eigen::MatrixBase<T1> &dataVec, int stride,
		    Eigen::MatrixBase<T2> &meanOut, Eigen::MatrixBase<T3> &covOut)
{
	if (stride == 0) return;
	int units = dataVec.size() / stride;
	meanOut.derived().resize(stride);
	meanOut.setZero();
	covOut.derived().resize(stride, stride);
	covOut.setZero();
	// read the data only once to minimize memory thrashing
	for (int ux=0; ux < units; ++ux) {
		meanOut += dataVec.segment(ux * stride, stride);
		covOut += (dataVec.segment(ux * stride, stride) *
			   dataVec.segment(ux * stride, stride).transpose());
	}
	meanOut /= units;
	covOut -= units * meanOut * meanOut.transpose();
	covOut /= units-1;
}

template <typename T1, typename T2>
double trace_prod(const Eigen::MatrixBase<T1> &t1, const Eigen::MatrixBase<T2> &t2)
{
	double sum = 0.0;
	for (int rx=0; rx < t1.rows(); ++rx) {
		sum += t1.row(rx) * t2.col(rx);
	}
	return sum;
}

void MoorePenroseInverse(Eigen::Ref<Eigen::MatrixXd> mat);

// https://forum.kde.org/viewtopic.php?f=74&t=96706
// https://forum.kde.org/viewtopic.php?f=74&t=124421
// https://forum.kde.org/viewtopic.php?f=74&t=91271
template <typename T1>
void filterJacobianRows(Eigen::MatrixBase<T1>& A, int& rankA){
	//TODO: check for conformability
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qra(A.transpose());
	rankA = qra.rank();
	Eigen::MatrixXd Q(A.cols(), A.rows());
	Q.setIdentity(A.cols(), A.rows());
	qra.householderQ().applyThisOnTheLeft(Q);
	Eigen::MatrixXd R = qra.matrixR().triangularView<Eigen::Upper>();
	R.conservativeResize(A.rows(), rankA);
	//mxLog("rank: %d",rankA[0]);
	//mxPrintMat("Q ",Q);
	//mxPrintMat("R ",R);
	A = (Q * R).transpose();
}

template <typename T> void omxMatrix::loadFromStream(T &st)
{
	using namespace mini::csv;
	EigenMatrixAdaptor v(this);

	switch(shape) {
	case 1: //Diag
		for (int rx=0; rx < rows; ++rx) {
			st >> v(rx, rx);
		}
		break;

	case 2: //Full
		for (int cx=0; cx < cols; ++cx) {
			for (int rx=0; rx < rows; ++rx) {
				st >> v(rx,cx);
			}
		}
		break;
		
	case 4: //Lower
		for (int cx=0; cx < cols; ++cx) {
			for (int rx=cx; rx < rows; ++rx) {
				st >> v(rx,cx);
			}
		}
		break;

	case 5: //Sdiag
		for (int cx=0; cx < cols-1; ++cx) {
			for (int rx=cx+1; rx < rows; ++rx) {
				st >> v(rx,cx);
			}
		}
		break;

	case 6: //Stand
		for (int cx=0; cx < cols-1; ++cx) {
			for (int rx=cx+1; rx < rows; ++rx) {
				double tmp;
				st >> tmp;
				v(rx,cx) = tmp;
				v(cx,rx) = tmp;
			}
		}
		break;

	case 7: //Symm
		for (int cx=0; cx < cols; ++cx) {
			for (int rx=cx; rx < rows; ++rx) {
				double tmp;
				st >> tmp;
				v(rx,cx) = tmp;
				v(cx,rx) = tmp;
			}
		}
		break;

	case 8: //Unit
	case 9: //Zero
	case 3: //Iden
		mxThrow("loadFromStream: matrix '%s' is constant (type %d);"
			 " use a Full matrix if you wish to update it", name(), shape);
		break;

	default:
		mxThrow("loadFromStream: matrix '%s' with shape %d is unimplemented",
			 name(), shape);
		break;
	}
}

void MatrixInvert1(omxMatrix *target);

#endif /* _OMXMATRIX_H_ */
