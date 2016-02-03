/*
  Copyright 2007-2016 The OpenMx Project

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

#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "omxBLAS.h"
#include "omxRAMExpectation.h"
#include "RAMInternal.h"
#include <Eigen/LU>

static omxMatrix* omxGetRAMExpectationComponent(omxExpectation* ox, const char* component);

void omxRAMExpectation::ensureTrivialF() // move to R side? TODO
{
	if (trivialF) return;

	omxRecompute(F, NULL);  // should not do anything
	EigenMatrixAdaptor eF(F);
	Eigen::MatrixXd ident(F->rows, F->rows);
	ident.setIdentity();
	if (ident != eF.block(0, 0, F->rows, F->rows)) {
		Rf_error("Square part of F matrix is not trivial");
	}
	trivialF = true;
}

static void omxCallRAMExpectation(omxExpectation* oo, FitContext *fc, const char *what, const char *how)
{
	omxRAMExpectation* oro = (omxRAMExpectation*)(oo->argStruct);

	if (what && how && strEQ(how, "flat")) {
		bool wantCov = false;
		bool wantMean = false;
		if (strEQ(what, "distribution")) { wantCov = true; wantMean = true; }
		if (strEQ(what, "covariance")) wantCov = true;
		if (strEQ(what, "mean")) wantMean = true;
		if (!oro->rram) {
			oro->rram = new RelationalRAMExpectation::state;
			oro->rram->init(oo, fc);
		}
		if (wantCov)  oro->rram->computeCov(fc);
		if (wantMean) oro->rram->computeMean(fc);
		return;
	}

	omxRecompute(oro->A, fc);
	omxRecompute(oro->S, fc);
	omxRecompute(oro->F, fc);
	if(oro->M != NULL)
	    omxRecompute(oro->M, fc);
	    
	oro->CalculateRAMCovarianceAndMeans();
}

static void omxDestroyRAMExpectation(omxExpectation* oo) {

	if(OMX_DEBUG) { mxLog("Destroying RAM Expectation."); }
	
	omxRAMExpectation* argStruct = (omxRAMExpectation*)(oo->argStruct);

	if (argStruct->rram) delete argStruct->rram;

	omxFreeMatrix(argStruct->cov);

	if(argStruct->means != NULL) {
		omxFreeMatrix(argStruct->means);
	}

	omxFreeMatrix(argStruct->I);
	omxFreeMatrix(argStruct->X);
	omxFreeMatrix(argStruct->Y);
	omxFreeMatrix(argStruct->Ax);
	delete argStruct;
}

static void refreshUnfilteredCov(omxExpectation *oo)
{
	// Ax = ZSZ' = Covariance matrix including latent variables
	omxRAMExpectation* oro = (omxRAMExpectation*) (oo->argStruct);
	omxMatrix* A = oro->A;
	omxMatrix* S = oro->S;
	omxMatrix* Ax= oro->Ax;
    
    omxRecompute(A, NULL);
    omxRecompute(S, NULL);
	
    omxMatrix* Z = oro->getZ(NULL);

    EigenMatrixAdaptor eZ(Z);
    EigenMatrixAdaptor eS(S);
    EigenMatrixAdaptor eAx(Ax);

    eAx.block(0, 0, eAx.rows(), eAx.cols()) = eZ * eS * eZ.transpose();
}

static void omxPopulateRAMAttributes(omxExpectation *oo, SEXP robj) {
    if(OMX_DEBUG) { mxLog("Populating RAM Attributes."); }

    refreshUnfilteredCov(oo);
	omxRAMExpectation* oro = (omxRAMExpectation*) (oo->argStruct);
	omxMatrix* Ax= oro->Ax;
	
	{
		ProtectedSEXP expCovExt(Rf_allocMatrix(REALSXP, Ax->rows, Ax->cols));
		memcpy(REAL(expCovExt), Ax->data, sizeof(double) * Ax->rows * Ax->cols);
		Rf_setAttrib(robj, Rf_install("UnfilteredExpCov"), expCovExt);
	}
	Rf_setAttrib(robj, Rf_install("numStats"), Rf_ScalarReal(omxDataDF(oo->data)));

	MxRList out;
	MxRList dbg;

	if (oro->rram) {
		RelationalRAMExpectation::state *rram = oro->rram;
		RelationalRAMExpectation::independentGroup &ig = rram->group;
		rram->exportInternalState(dbg);
		ig.exportInternalState(out, dbg);
	}

	Rf_setAttrib(robj, Rf_install("output"), out.asR());
	Rf_setAttrib(robj, Rf_install("debug"), dbg.asR());
}

// reimplement inverse using eigen::sparsematrix TODO
omxMatrix *omxRAMExpectation::getZ(FitContext *fc)
{
	if (Zversion != omxGetMatrixVersion(A)) {
		omxShallowInverse(fc, numIters, A, _Z, Ax, I);
		Zversion = omxGetMatrixVersion(A);
	}
	return _Z;
}

/*
 * omxCalculateRAMCovarianceAndMeans
 * 			Just like it says on the tin.  Calculates the mean and covariance matrices
 * for a RAM model.  M is the number of total variables, latent and manifest. N is
 * the number of manifest variables.
 *
 * params:
 * omxMatrix *A, *S, *F 	: matrices as specified in the RAM model.  MxM, MxM, and NxM
 * omxMatrix *M				: vector containing model implied means. 1xM
 * omxMatrix *Cov			: On output: model-implied manifest covariance.  NxN.
 * omxMatrix *Means			: On output: model-implied manifest means.  1xN.
 * int numIterations		: Precomputed number of iterations of taylor series expansion.
 * omxMatrix *I				: Identity matrix.  If left NULL, will be populated.  MxM.
 * omxMatrix *Z				: On output: Computed (I-A)^-1. MxM.
 * omxMatrix *Y, *X, *Ax	: Space for computation. NxM, NxM, MxM.  On exit, populated.
 */

void omxRAMExpectation::CalculateRAMCovarianceAndMeans()
{
	if (F->rows == 0) return;

	if(OMX_DEBUG) { mxLog("Running RAM computation with numIters is %d\n.", numIters); }
		
	if(Ax == NULL || I == NULL || Y == NULL || X == NULL) {
		Rf_error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
	}
		
	if(cov == NULL && means == NULL) {
		return; // We're not populating anything, so why bother running the calculation?
	}
	
	omxMatrix *Z = getZ(NULL);
	
	/* Cov = FZSZ'F' */
	omxDGEMM(FALSE, FALSE, 1.0, F, Z, 0.0, Y);

	omxDGEMM(FALSE, FALSE, 1.0, Y, S, 0.0, X);

	omxDGEMM(FALSE, TRUE, 1.0, X, Y, 0.0, cov);
	 // Cov = FZSZ'F' (Because (FZ)' = Z'F')
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(cov, "....RAM: Model-implied Covariance Matrix:");}
	
	if(M != NULL && means != NULL) {
		// F77_CALL(omxunsafedgemv)(Y->majority, &(Y->rows), &(Y->cols), &oned, Y->data, &(Y->leading), M->data, &onei, &zerod, means->data, &onei);
		omxDGEMV(FALSE, 1.0, Y, M, 0.0, means);
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(means, "....RAM: Model-implied Means Vector:");}
	}
}

struct IsZeroManifestCovariance {
	typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> ManifestMaskType;
	typedef Eigen::MatrixXd::Scalar Scalar;
	typedef Eigen::MatrixXd::Index Index;
	bool ok;
	ManifestMaskType manifestMask;

	IsZeroManifestCovariance(ManifestMaskType manifestMask) : ok(true), manifestMask(manifestMask) {};
	void init(const Scalar& value, Index i, Index j) { check(value,i,j); };
	void operator() (const Scalar& value, Index i, Index j) { check(value,i,j); };
	void check(const Scalar& value, Index i, Index j) {
		if (i != j && value != 0.0 && (manifestMask[i] || manifestMask[j])) {
			ok = false;
		}
	};
};
		
void omxInitRAMExpectation(omxExpectation* oo) {
	
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;

	if(OMX_DEBUG) { mxLog("Initializing RAM expectation."); }
	
	int l, k;

	SEXP slotValue;
	
	omxMatrix *Zmat = omxInitMatrix(0, 0, TRUE, currentState);
	omxRAMExpectation *RAMexp = new omxRAMExpectation(Zmat);
	RAMexp->rram = 0;
	
	/* Set Expectation Calls and Structures */
	oo->computeFun = omxCallRAMExpectation;
	oo->destructFun = omxDestroyRAMExpectation;
	oo->componentFun = omxGetRAMExpectationComponent;
	oo->populateAttrFun = omxPopulateRAMAttributes;
	oo->argStruct = (void*) RAMexp;
	
	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	RAMexp->verbose = Rf_asInteger(Rverbose) + OMX_DEBUG;

	/* Set up expectation structures */
	if(OMX_DEBUG) { mxLog("Initializing RAM expectation."); }

	if(OMX_DEBUG) { mxLog("Processing M."); }
	RAMexp->M = omxNewMatrixFromSlot(rObj, currentState, "M");

	if(OMX_DEBUG) { mxLog("Processing A."); }
	RAMexp->A = omxNewMatrixFromSlot(rObj, currentState, "A");

	if(OMX_DEBUG) { mxLog("Processing S."); }
	RAMexp->S = omxNewMatrixFromSlot(rObj, currentState, "S");

	if(OMX_DEBUG) { mxLog("Processing F."); }
	RAMexp->F = omxNewMatrixFromSlot(rObj, currentState, "F");

	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { mxLog("Generating I."); }
	RAMexp->I = omxNewIdentityMatrix(RAMexp->A->rows, currentState);

	if(OMX_DEBUG) { mxLog("Processing expansion iteration depth."); }
	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("depth")));
	RAMexp->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { mxLog("Using %d iterations.", RAMexp->numIters); }
	}

	ProtectedSEXP Rrampart(R_do_slot(rObj, Rf_install("rampart")));
	RAMexp->rampart = Rf_asInteger(Rrampart);

	ProtectedSEXP Rbetween(R_do_slot(rObj, Rf_install("between")));
	if (Rf_length(Rbetween)) {
		if (!oo->data) Rf_error("%s: data is required for joins", oo->name);
		if (!Rf_isInteger(Rbetween)) Rf_error("%s: between must be an integer vector", oo->name);
		RAMexp->between.reserve(Rf_length(Rbetween));
		int *bnumber = INTEGER(Rbetween);
		for (int jx=0; jx < Rf_length(Rbetween); ++jx) {
			omxMatrix *bmat = currentState->getMatrixFromIndex(bnumber[jx]);
			int foreignKey = bmat->getJoinKey();
			omxExpectation *fex = bmat->getJoinModel();
			omxCompleteExpectation(fex);
			if (!strEQ(fex->expType, "MxExpectationRAM")) {
				Rf_error("%s: only MxExpectationRAM can be joined with MxExpectationRAM", oo->name);
			}
			if (omxDataIsSorted(fex->data)) {
				Rf_error("%s join with %s but observed data is sorted",
					 oo->name, fex->name);
			}
			if (!omxDataColumnIsKey(oo->data, foreignKey)) {
				Rf_error("Cannot join using non-integer type column '%s' in '%s'. "
					 "Did you forget to use mxData(..., sort=FALSE)?",
					 omxDataColumnName(oo->data, foreignKey),
					 oo->data->name);
			}

			if (OMX_DEBUG) {
				mxLog("%s: join col %d against %s using between matrix %s",
				      oo->name, foreignKey, fex->name, bmat->name());
			}
				
			RAMexp->between.push_back(bmat);
		}
	}

	l = RAMexp->F->rows;
	k = RAMexp->A->cols;

	if (k != RAMexp->S->cols || k != RAMexp->S->rows || k != RAMexp->A->rows) {
		Rf_error("RAM matrices '%s' and '%s' must have the same dimensions",
			 RAMexp->S->name(), RAMexp->A->name());
	}

	if(OMX_DEBUG) { mxLog("Generating internals for computation."); }

	omxResizeMatrix(Zmat, k, k);

	RAMexp->Ax = 	omxInitMatrix(k, k, TRUE, currentState);
	RAMexp->Ax->rownames = RAMexp->S->rownames;
	RAMexp->Ax->colnames = RAMexp->S->colnames;
	RAMexp->Y = 	omxInitMatrix(l, k, TRUE, currentState);
	RAMexp->X = 	omxInitMatrix(l, k, TRUE, currentState);
	
	RAMexp->cov = 		omxInitMatrix(l, l, TRUE, currentState);

	if(RAMexp->M != NULL) {
		RAMexp->means = 	omxInitMatrix(1, l, TRUE, currentState);
	} else {
	    RAMexp->means  = 	NULL;
    }
}

static omxMatrix* omxGetRAMExpectationComponent(omxExpectation* ox, const char* component) {
	
	if(OMX_DEBUG) { mxLog("RAM expectation: %s requested--", component); }

	omxRAMExpectation* ore = (omxRAMExpectation*)(ox->argStruct);
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = ore->cov;
	} else if(strEQ("means", component)) {
		retval = ore->means;
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}
	
	return retval;
}

namespace RelationalRAMExpectation {

	int addr::numVars() const
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) model->argStruct;
		return ram->F->cols;
	}

	// verify whether sparse can deal with parameters set to exactly zero TODO

	void independentGroup::refreshUnitA(FitContext *fc, int px)
	{
		struct placement &pl = placements[px];
		addr &a1 = st.layout[pl.aIndex];
		omxExpectation *expectation = a1.model;
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		EigenMatrixAdaptor eA(ram->A);
		for (int cx=0; cx < eA.cols(); ++cx) {
			for (int rx=0; rx < eA.rows(); ++rx) {
				double val = eA(rx, cx);
				if (val != 0) {
					if (rx == cx) {
						Rf_error("%s: nonzero diagonal entry in A matrix at %d",
							 st.homeEx->name, 1+pl.modelStart+cx);
					}
					// can't use eA.block(..) -= because fullA must remain sparse
					fullA.coeffRef(pl.modelStart + cx, pl.modelStart + rx) = signA * val;
				}
			}
		}

		const double scale = a1.rampartScale;
		if (scale == 0.0) return;

		for (size_t jx=0; jx < ram->between.size(); ++jx) {
			omxMatrix *betA = ram->between[jx];
			int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
			if (key == NA_INTEGER) continue;
			omxData *data1 = betA->getJoinModel()->data;
			int frow = data1->lookupRowOfKey(key);
			placement &p2 = placements[ rowToPlacementMap[std::make_pair(data1, frow)] ];
			omxRecompute(betA, fc);
			omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel()->argStruct;
			for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
				for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
					double val = omxMatrixElement(betA, rx, cx);
					if (val == 0.0) continue;
					fullA.coeffRef(p2.modelStart + cx, pl.modelStart + rx) =
						signA * val * scale;
				}
			}
		}
	}

	void independentGroup::refreshModel(FitContext *fc)
	{
		for (size_t ax=0; ax < placements.size(); ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[pl.aIndex];
			omxExpectation *expectation = a1.model;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
			omxData *data = expectation->data;
			data->handleDefinitionVarList(expectation->currentState, a1.row);
			omxRecompute(ram->A, fc);
			omxRecompute(ram->S, fc);

			refreshUnitA(fc, ax);

			EigenMatrixAdaptor eS(ram->S);
			for (int cx=0; cx < eS.cols(); ++cx) {
				for (int rx=cx; rx < eS.rows(); ++rx) {
					if (eS(rx,cx) != 0) {
						fullS.coeffRef(pl.modelStart + rx, pl.modelStart + cx) = eS(rx, cx);
					}
				}
			}
		}
	}

	// 3rd visitor
	void state::refreshModel(FitContext *fc)
	{
		group.refreshModel(fc);
	}

	// 1st visitor
	int state::flattenOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize, int level)
	{
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		struct addr a1;
		a1.group = 0;
		a1.copy = 0;
		a1.level = level;
		a1.clumped = false;
		a1.rampartScale = 1.0;
		a1.parent1 = NA_INTEGER;
		a1.fk1 = NA_INTEGER;
		a1.numJoins = 0;
		int parent1Pos = -1;
		for (size_t jx=0; jx < ram->between.size(); ++jx) {
			omxMatrix *b1 = ram->between[jx];
			int key = omxKeyDataElement(data, frow, b1->getJoinKey());
			if (key == NA_INTEGER) continue;
			a1.numJoins += 1;
			if (a1.numJoins == 1) a1.fk1 = key;
			omxExpectation *e1 = b1->getJoinModel();
			int parentPos = flattenOneRow(e1, e1->data->lookupRowOfKey(key), totalObserved, maxSize, level+1);
			if (a1.numJoins == 1) parent1Pos = parentPos;
		}

		a1.row = frow;
		a1.key = frow;
		if (data->hasPrimaryKey()) {
			ram->ensureTrivialF();

			// insert_or_assign would be nice here
			RowToLayoutMapType::const_iterator it = rowToLayoutMap.find(std::make_pair(data, frow));
			if (it != rowToLayoutMap.end()) return it->second;

			rowToLayoutMap[ std::make_pair(data, frow) ] = layout.size();
			a1.key = data->primaryKeyOfRow(frow);
		}

		if (parent1Pos >= 0) {
			addr &pop = layout[parent1Pos];
			pop.numKids += 1;
			a1.parent1 = parent1Pos;
		}
		a1.model = expectation;
		a1.numKids = 0;
		int obsStart = totalObserved;

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

		a1.numObsCache = totalObserved - obsStart;
		layout.push_back(a1);

		maxSize += ram->F->cols;
		return layout.size()-1;
	}

	void independentGroup::determineShallowDepth(FitContext *fc)
	{
		if (!Global->RAMInverseOpt) return;

		Eigen::VectorXd vec(fc->varGroup->vars.size());
		vec.setConstant(1);
		copyParamToModelInternal(fc->varGroup, st.homeEx->currentState, vec.data());

		for (size_t ax=0; ax < placements.size(); ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[ax];
			omxExpectation *expectation = a1.model;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
			omxData *data = expectation->data;

			data->handleDefinitionVarList(expectation->currentState, a1.row);
			omxRecompute(ram->A, NULL);

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxData *data1 = betA->getJoinModel()->data;
				int frow = data1->lookupRowOfKey(key);
				placement &p2 = placements[ rowToPlacementMap[std::make_pair(data1, frow)] ];
				omxRecompute(betA, NULL);
				betA->markPopulatedEntries();
				omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel()->argStruct;
				for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
					for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
						double val = omxMatrixElement(betA, rx, cx);
						if (val == 0.0) continue;
						fullA.coeffRef(p2.modelStart + cx, pl.modelStart + rx) = 1;
					}
				}
			}
			ram->A->markPopulatedEntries();
			EigenMatrixAdaptor eA(ram->A);
			for (int cx=0; cx < eA.cols(); ++cx) {
				for (int rx=0; rx < eA.rows(); ++rx) {
					if (rx != cx && eA(rx,cx) != 0) {
						fullA.coeffRef(pl.modelStart + cx, pl.modelStart + rx) = 1;
					}
				}
			}
		}

		int maxDepth = std::min(fullA.cols(), 30);
		if (Global->RAMMaxDepth != NA_INTEGER) maxDepth = Global->RAMMaxDepth;
		Eigen::SparseMatrix<double> curProd = fullA;
		for (int tx=1; tx < maxDepth; ++tx) {
			if (verbose() >= 3) { Eigen::MatrixXd tmp = curProd; mxPrintMat("curProd", tmp); }
			curProd = (curProd * fullA.transpose()).eval();
			bool allZero = true;
			for (int k=0; k < curProd.outerSize(); ++k) {
				for (Eigen::SparseMatrix<double>::InnerIterator it(curProd, k); it; ++it) {
					if (it.value() != 0.0) {
						allZero = false;
						break;
					}
				}
			}
			if (allZero) {
				AshallowDepth = tx;
				break;
			}
		}
		fullA.setZero();
		fc->copyParamToModelClean();

		if (AshallowDepth >= 0) signA = 1.0;
		if (verbose() >= 1) {
			mxLog("%s: RAM shallow inverse depth = %d", st.homeEx->name, AshallowDepth);
		}
	}

	struct CompareLib {
		state &st;
		CompareLib(state *st) : st(*st) {};

		// actually stores !missingness
		template <typename T>
		void getMissingnessPattern(const addr &a1, std::vector<T> &out) const
		{
			omxDataRow(a1.model, a1.row, st.smallCol);
			int jCols = a1.model->dataColumns->cols;
			out.reserve(jCols);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(st.smallCol, 0, col);
				out.push_back(std::isfinite(val));
			}
		}

		bool compareModelAndMissingness(addr &la, addr &ra, bool &mismatch) const
		{
			mismatch = true;
			if (la.model != ra.model)
				return strcmp(la.model->name, ra.model->name) < 0;

			if (la.numVars() != ra.numVars())
				return la.numVars() < ra.numVars();

			std::vector<bool> lmp;
			getMissingnessPattern(la, lmp);
			std::vector<bool> rmp;
			getMissingnessPattern(ra, rmp);

			if (lmp.size() != rmp.size())
				return lmp.size() < rmp.size();

			for (size_t lx=0; lx < lmp.size(); ++lx) {
				if (lmp[lx] == rmp[lx]) continue;
				return int(lmp[lx]) < int(rmp[lx]);
			}
			mismatch = false;
			return false;
		}
	};

	struct WithinClumpCompare : CompareLib {
		WithinClumpCompare(state *st) : CompareLib(st) {};

		bool operator() (const int lhs, const int rhs) const
		{
			addr &la = st.layout[lhs];
			addr &ra = st.layout[rhs];

			if (la.level != ra.level)
				return la.level < ra.level;

			bool mismatch;
			return compareModelAndMissingness(la, ra, mismatch);
		}
	};

	struct CompatibleGroupCompare : CompareLib {
		CompatibleGroupCompare(state *st) : CompareLib(st) {};

		bool operator() (const std::vector<int> &lhs, const std::vector<int> &rhs) const
		{
			if (lhs.size() != rhs.size()) return lhs.size() < rhs.size();
			for (size_t ux=0; ux < lhs.size(); ++ux) {
				addr &la = st.layout[lhs[ux]];
				addr &ra = st.layout[rhs[ux]];
				bool mismatch;
				bool got = compareModelAndMissingness(la, ra, mismatch);
				if (mismatch) return got;
				// defvars TODO
			}
			return false;
		}
	};

	// 2nd visitor
	void state::planModelEval(int maxSize, int totalObserved, FitContext *fc)
	{
		typedef std::map< std::vector<int>,
				  std::set<std::vector<int> >,
				  CompatibleGroupCompare> CompatibleGroupMapType;

		macroA.resize(layout.size(), layout.size());
		macroA.setIdentity();
		group.placements.reserve(layout.size());

		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			if (a1.rampartScale == 0.0) continue;
			omxRAMExpectation *ram = (omxRAMExpectation*) a1.model->argStruct;
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *b1 = ram->between[jx];
				int key = omxKeyDataElement(a1.model->data, a1.row, b1->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxExpectation *e1 = b1->getJoinModel();
				int row = e1->data->lookupRowOfKey(key);
				RowToLayoutMapType::const_iterator it =
					rowToLayoutMap.find(std::make_pair(e1->data, row));
				if (it == rowToLayoutMap.end())
					Rf_error("Cannot find row %d in %s", row, e1->data->name);
				macroA(it->second, ax) = -1;
			}
		}

		Eigen::FullPivLU<Eigen::MatrixXf> lu(macroA);
		Eigen::MatrixXf macroAi = lu.inverse();

		CompatibleGroupMapType cgm(this);
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			if (a1.numJoins) continue;
			int clumpSize = (macroAi.row(ax).array() != 0).count();
			std::vector<int> clump;
			clump.reserve(clumpSize);
			for (size_t mx=0; mx < layout.size(); ++mx) {
				if (macroAi(ax, mx)) clump.push_back(mx);
			}
			std::sort(clump.begin(), clump.end(), WithinClumpCompare(this));
			cgm[ clump ].insert(clump);
		}

		int groupNum = 1;
		for (CompatibleGroupMapType::iterator it = cgm.begin();
		     it != cgm.end(); ++it, ++groupNum) {
			int copyNum = 1;
			for (std::set<std::vector<int> >::iterator px = it->second.begin();
			     px != it->second.end(); ++px, ++copyNum) {
				const std::vector<int> &clump = *px;
				for (size_t cx=0; cx < clump.size(); ++cx) {
					addr &a1 = layout[ clump[cx] ];
					if (a1.group) Rf_error("Unit[%d] already assigned a group; this is a bug", clump[cx]);
					a1.group = groupNum;
					a1.copy = copyNum;
				}
			}
		}
 
		int dx = 0, mx = 0;
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];

			placement pl;
			pl.aIndex = ax;
			pl.modelStart = mx;
			pl.obsStart = dx;
			group.placements.push_back(pl);
			mx += a1.numVars();
			dx += a1.numObs();
		}

		group.prep(maxSize, totalObserved, fc);

		// after rampart then models are either clumped or independent?

		// numJoins == 0
		// whole clump:
		//   no defvars or defvars match
		//   missingness matches

		rotationPlan.clear();
	}

	void independentGroup::prep(int maxSize, int totalObserved, FitContext *fc)
	{
		ident.resize(maxSize, maxSize);
		ident.setIdentity();
		fullA.resize(maxSize, maxSize);
		latentFilter.assign(maxSize, false); // will have totalObserved true entries
		obsNameVec = Rf_protect(Rf_allocVector(STRSXP, totalObserved));
		varNameVec = Rf_protect(Rf_allocVector(STRSXP, maxSize));
		dataVec.resize(totalObserved);
		dataColumn.resize(totalObserved);
		dataColumn.setConstant(-1);
		fullMean.resize(maxSize);
		fullMean.setZero();
		fullS.conservativeResize(maxSize, maxSize);

		int dx=0;
		for (size_t ax=0; ax < placements.size(); ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[ pl.aIndex ];
			// also build index TODO
			if (verbose() >= 2) {
				int modelEnd = pl.modelStart + a1.numVars() - 1;
				if (a1.numObs()) {
					mxLog("place %s[%d] at %d %d obs %d %d", a1.modelName().c_str(),
					      a1.row, pl.modelStart, modelEnd, pl.obsStart, pl.obsStart + a1.numObs() - 1);
				} else {
					mxLog("place latent %s[%d] at %d %d", a1.modelName().c_str(),
					      a1.row, pl.modelStart, modelEnd);
				}
			}

			omxData *data = a1.model->data;
			rowToPlacementMap[ std::make_pair(data, a1.row) ] = ax;

			omxRAMExpectation *ram = (omxRAMExpectation*) a1.model->argStruct;

			std::string modelName(data->name);
			modelName = modelName.substr(0, modelName.size() - 4); // remove "data" suffix

			int jCols = a1.model->dataColumns->cols;
			if (jCols) {
				omxDataRow(a1.model, a1.row, st.smallCol);
				omxMatrix *colList = a1.model->dataColumns;
				for (int col=0; col < jCols; ++col) {
					double val = omxMatrixElement(st.smallCol, 0, col);
					bool yes = std::isfinite(val);
					if (!yes) continue;
					latentFilter[ pl.modelStart + col ] = true;
					std::string dname =
						modelName + omxDataColumnName(data, omxVectorElement(colList, col));
					SET_STRING_ELT(obsNameVec, dx, Rf_mkChar(dname.c_str()));
					dataVec[ dx ] = val;
					if (a1.model == st.homeEx) dataColumn[ dx ] = col;
					dx += 1;
				}
			}
			for (int vx=0; vx < ram->F->cols; ++vx) {
				std::string dname = modelName + ram->F->colnames[vx];
				SET_STRING_ELT(varNameVec, pl.modelStart + vx, Rf_mkChar(dname.c_str()));
			}
		}

		determineShallowDepth(fc);

		// need to remap in opposite direction of aIndex TODO
		rotationPlan = st.rotationPlan;
	}

	struct RampartCompareLib {
		state *st;
		RampartCompareLib(state *st) : st(st) {};

		omxExpectation *getJoinModel(const addr *a1) const {
			omxRAMExpectation *ram = (omxRAMExpectation*) a1->model->argStruct;
			omxMatrix *b1 = ram->between[0];
			return b1->getJoinModel();
		};

		// actually stores !missingness
		template <typename T>
		void getMissingnessPattern(const addr *a1, std::vector<T> &out) const
		{
			omxDataRow(a1->model, a1->row, st->smallCol);
			int jCols = a1->model->dataColumns->cols;
			out.reserve(jCols);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(st->smallCol, 0, col);
				out.push_back(std::isfinite(val));
			}
		}

		bool cmpUpper(const addr *lhs, const addr *rhs, bool &result) const
		{
			if (getJoinModel(lhs) != getJoinModel(rhs)) {
				result = strcmp(getJoinModel(lhs)->name, getJoinModel(rhs)->name) < 0;
				return true;
			}
			return false;
		}

		bool cmp1(const addr *lhs, const addr *rhs, bool &result) const
		{
			if (lhs->model != rhs->model) {
				result = strcmp(lhs->model->name, rhs->model->name) < 0;
				return true;
			}
			if (lhs->numVars() != rhs->numVars()) {
				result = lhs->numVars() < rhs->numVars();
				return true;
			}

			std::vector<bool> lmp;
			getMissingnessPattern(lhs, lmp);
			std::vector<bool> rmp;
			getMissingnessPattern(rhs, rmp);

			if (lmp.size() != rmp.size()) {
				result = lmp.size() < rmp.size();
				return true;
			}
			for (size_t lx=0; lx < lmp.size(); ++lx) {
				if (lmp[lx] == rmp[lx]) continue;
				result = int(lmp[lx]) < int(rmp[lx]);
				return true;
			}
			if (lhs->clump.size() != rhs->clump.size()) {
				result = lhs->clump.size() < rhs->clump.size();
				return true;
			}
			for (size_t cx=0; cx < lhs->clump.size(); ++cx) {
				if (cmp1(&st->layout[lhs->clump[cx]], &st->layout[rhs->clump[cx]], result))
					return true;
			}
			return false;
		}
	};

	struct RampartTodoCompare : RampartCompareLib {
		RampartTodoCompare(state *st) : RampartCompareLib(st) {};

		bool operator() (const addr *lhs, const addr *rhs) const
		{
			bool result = false;
			if (cmpUpper(lhs, rhs, result)) return result;
			if (lhs->fk1 != rhs->fk1)
				return lhs->fk1 < rhs->fk1;
			cmp1(lhs, rhs, result);
			return result;
		}
	};

	struct RampartClumpCompare : RampartCompareLib {
		RampartClumpCompare(state *st) : RampartCompareLib(st) {};

		bool clumpCmp(const int lhs, const int rhs) const {
			const addr *lhsObj = &st->layout[lhs];
			const addr *rhsObj = &st->layout[rhs];
			bool result = false;
			if (cmp1(lhsObj, rhsObj, result)) return result;
			return lhs < rhs;
		}

		bool operator() (const int lhs, const int rhs) const
		{
			bool got = clumpCmp(lhs, rhs);
			//mxLog("cmp %d %d -> %d", lhs, rhs, got);
			return got;
		}
	};

	template <typename T>
	void state::oertzenRotate(std::vector<T> &t1)
	{
		// get covariance and sort by mahalanobis distance TODO
		rotationPlan.push_back(t1);
		addr &specimen = layout[ t1[0] ];
		for (size_t cx=0; cx < specimen.clump.size(); ++cx) {
			std::vector<int> t2;
			t2.reserve(t1.size());
			for (size_t tx=0; tx < t1.size(); ++tx) {
				addr &a1 = layout[ t1[tx] ];
				t2.push_back(a1.clump[cx]);
			}
			oertzenRotate(t2);
		}
	}

	int state::rampartRotate(int level)
	{
		typedef std::map< addr*, std::vector<int>, RampartTodoCompare > RampartMap;
		RampartMap todo(RampartTodoCompare(this));
		int unlinked = 0;

		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr& a1 = layout[ax];
			if (a1.numKids != 0 || a1.numJoins != 1 || a1.clumped || a1.rampartScale != 1.0) continue;

			omxRAMExpectation *ram = (omxRAMExpectation*) a1.model->argStruct;
			omxMatrix *b1 = ram->between[0];
			if (b1->dependsOnDefinitionVariables()) continue;
			std::vector<int> &t1 = todo[&a1];
			t1.push_back(int(ax));
		}
		for (RampartMap::iterator it = todo.begin(); it != todo.end(); ++it) {
			std::vector<int> &t1 = it->second;

			if (t1.size() >= 2) {
				//std::string buf = "rotate units";
				oertzenRotate(t1);
				addr &specimen = layout[ t1[0] ];
				specimen.rampartScale = sqrt(double(t1.size()));
				if (false) {
				for (size_t cx=0; cx < specimen.clump.size(); ++cx) {
					addr &a1 = layout[ specimen.clump[cx] ];
					double orig = a1.rampartScale*a1.rampartScale;
					a1.rampartScale = sqrt(orig / double(t1.size())) * sqrt(orig);
					a1.clumped = true;
				}
				}
				layout[specimen.parent1].numKids -= t1.size();
				layout[specimen.parent1].clump.push_back(t1[0]);
				for (size_t ux=1; ux < t1.size(); ++ux) {
					addr &a1 = layout[ t1[ux] ];
					a1.rampartScale = 0;
					a1.numJoins = 0;
					//buf += string_snprintf(" %d", 1+t1[ux]);
				}
				//buf += string_snprintf(" -> %d\n", 1+t1[0]);
				//mxLogBig(buf);
			} else {
				// Don't rotate, just clump units together with parent.
				addr &specimen = layout[ t1[0] ];
				layout[specimen.parent1].numKids -= t1.size();
				layout[specimen.parent1].clump.insert(layout[specimen.parent1].clump.end(),
								      t1.begin(), t1.end());
				for (size_t ux=0; ux < t1.size(); ++ux) {
					addr &a1 = layout[ t1[ux] ];
					a1.clumped = true;
				}
			}
			// not really unlinked in clumped case, but layout is changed; fix reporting TODO
			unlinked += t1.size() - 1;
		}
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr& a1 = layout[ax];
			if (a1.clump.size() <= 1) continue;
			std::sort(a1.clump.begin(), a1.clump.end(), RampartClumpCompare(this));
			if (false) {
				// if Compare is screwed up then it should show up here
				std::vector<int> dc(a1.clump);
				std::sort(dc.begin(), dc.end(), RampartClumpCompare(this));
				for (size_t cx=0; cx < dc.size(); ++cx) {
					if (dc[cx] != a1.clump[cx]) Rf_error("oops");
				}
			}
		}
		return unlinked;
	}

	void state::init(omxExpectation *expectation, FitContext *fc)
	{
		homeEx = expectation;

		omxRAMExpectation *ram = (omxRAMExpectation*) homeEx->argStruct;
		ram->ensureTrivialF();

		int numManifest = ram->F->rows;

		smallCol = omxInitMatrix(1, numManifest, TRUE, homeEx->currentState);
		omxData *data = homeEx->data;

		int totalObserved = 0;
		int maxSize = 0;
		for (int row=0; row < data->rows; ++row) {
			flattenOneRow(homeEx, row, totalObserved, maxSize, 1);
		}

		if (verbose() >= 1) {
			mxLog("%s: total observations %d", homeEx->name, totalObserved);
		}

		if (ram->rampartEnabled()) {
			int maxIter = ram->rampart;
			int unlinked = 0;
			int level = -1; // mainly for debugging
			while (int more = rampartRotate(++level)) {
				rampartUsage.push_back(more);
				unlinked += more;
				if (--maxIter == 0) break;
			}
			if (verbose() >= 1) {
				mxLog("%s: rampart unlinked %d units", homeEx->name, unlinked);
			}
		}

		planModelEval(maxSize, totalObserved, fc);
		rowToLayoutMap.clear();
	}

	state::~state()
	{
		omxFreeMatrix(smallCol);
	}

	void independentGroup::invertAndFilterA()
	{
		// consider http://users.clas.ufl.edu/hager/papers/Lightning/update.pdf ?
		if (AshallowDepth >= 0) {
			fullA.makeCompressed();
			IAF = fullA + ident;
			for (int iter=1; iter <= AshallowDepth; ++iter) {
				IAF = (IAF * fullA + ident).eval();
				//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
			}
		} else {
			fullA += ident;
			if (!analyzed) {
				analyzed = true;
				fullA.makeCompressed();
				Asolver.analyzePattern(fullA);
			}
			Asolver.factorize(fullA);
			if (Asolver.info() != Eigen::Success) {
				Rf_error("%s: failed to invert flattened A matrix; %s",
					 st.homeEx->name, Asolver.lastErrorMessage().c_str());
			}

			IAF = Asolver.solve(ident);
			fullA -= ident;  // leave unchanged
			//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
		}

		const bool doubleCheck = true;
		Eigen::MatrixXd denseA;
		if (doubleCheck) {
			denseA = IAF;
		}

		// We built A transposed so we can quickly filter columns
		// Switch to filterOuter http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1130 TODO
		IAF.uncompress();
		Eigen::SparseMatrix<double>::Index *op = IAF.outerIndexPtr();
		Eigen::SparseMatrix<double>::Index *nzp = IAF.innerNonZeroPtr();
		int dx = 0;
		for (int cx=0; cx < fullA.cols(); ++cx) {
			if (!latentFilter[cx]) continue;
			op[dx] = op[cx];
			nzp[dx] = nzp[cx];
			++dx;
		}
		op[dx] = op[fullA.cols()];
		IAF.conservativeResize(fullA.rows(), dataVec.size());

		if (doubleCheck) {
			Eigen::MatrixXd denseAF;
			denseAF.resize(fullA.rows(), dataVec.size());
			int dx=0;
			for (int cx=0; cx < fullA.cols(); ++cx) {
				if (!latentFilter[cx]) continue;
				denseAF.col(dx) = denseA.col(cx);
				++dx;
			}
			if (dx != dataVec.size()) Rf_error("latentFilter has wrong count %d != %d",
							   dx, dataVec.size());
			Eigen::MatrixXd denseFilteredA = IAF;
			if ((denseAF.array() != denseFilteredA.array()).any()) {
				for (int rx=0; rx<denseAF.rows(); ++rx) {
					for (int cx=0; cx<denseAF.cols(); ++cx) {
						if (denseAF.coeff(rx,cx) != denseFilteredA.coeff(rx,cx)) {
							mxLog("[%d,%d] %f != %f",
							      rx, cx, denseAF.coeff(rx,cx), denseFilteredA.coeff(rx,cx));
						}
					}
				}
				Rf_error("stop");
			}
		}
		//{ Eigen::MatrixXd tmp = out; mxPrintMat("out", tmp); }
	}

	void independentGroup::computeCov(FitContext *fc)
	{
		if (fullA.nonZeros() == 0) {
			fullA.resize(latentFilter.size(), latentFilter.size());
		}
		fullS.conservativeResize(latentFilter.size(), latentFilter.size());

		refreshModel(fc);

		//{ Eigen::MatrixXd tmp = fullA; mxPrintMat("fullA", tmp); }
		//{ Eigen::MatrixXd tmp = fullS; mxPrintMat("fullS", tmp); }

		invertAndFilterA();

		//mxPrintMat("S", fullS);
		//Eigen::MatrixXd fullCovDense =

		fullCov = (IAF.transpose() * fullS.selfadjointView<Eigen::Lower>() * IAF);
		//mxLog("fullCov %d%% nonzero", int(fullCov.nonZeros() * 100.0 / (fullCov.rows() * fullCov.cols())));
		//{ Eigen::MatrixXd tmp = fullCov; mxPrintMat("fullcov", tmp); }
	}

	void state::computeCov(FitContext *fc)
	{
		group.computeCov(fc);
	}

	void independentGroup::computeMean(FitContext *fc)
	{
		expectedMean.conservativeResize(dataVec.size());

		for (size_t ax=0; ax < placements.size(); ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[pl.aIndex];
			omxExpectation *expectation = a1.model;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
			if (!ram->M) continue;

			omxData *data = expectation->data;
			data->handleDefinitionVarList(expectation->currentState, a1.row);
			omxRecompute(ram->A, fc);
			omxRecompute(ram->M, fc);
			EigenMatrixAdaptor eZ(ram->getZ(fc));
			EigenVectorAdaptor eM(ram->M);
			fullMean.segment(pl.modelStart, a1.numVars()) = eZ * eM;

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxData *data1 = betA->getJoinModel()->data;
				int frow = data1->lookupRowOfKey(key);
				placement &p2 = placements[ rowToPlacementMap[std::make_pair(data1, frow)] ];
				omxRecompute(betA, fc);
				EigenMatrixAdaptor eBA(betA);
				fullMean.segment(pl.modelStart, a1.numVars()) +=
					eBA * fullMean.segment(p2.modelStart, eBA.cols());
			}
		}
		int ox = 0;
		for (size_t lx=0; lx < latentFilter.size(); ++lx) {
			if (!latentFilter[lx]) continue;
			expectedMean[ox++] = fullMean[lx];
		}
	}

	void state::computeMean(FitContext *fc)
	{
		group.computeMean(fc);
	}

	static int plusOne(int val) {
		if (val == NA_INTEGER) return val;
		else return val + 1;
	}

	Eigen::SparseMatrix<double> independentGroup::getInputMatrix() const
	{
		return signA * fullA.transpose();
	}

	void independentGroup::exportInternalState(MxRList &out, MxRList &dbg)
	{
		if (expectedMean.size()) {
			SEXP m1 = Rcpp::wrap(expectedMean);
			Rf_protect(m1);
			Rf_setAttrib(m1, R_NamesSymbol, obsNameVec);
			out.add("mean", m1);
		}
		if (fullCov.nonZeros()) {
			out.add("covariance", Rcpp::wrap(fullCov));
		}

		SEXP fmean = Rcpp::wrap(fullMean);
		dbg.add("mean", fmean);
		Rf_setAttrib(fmean, R_NamesSymbol, varNameVec);
		Eigen::SparseMatrix<double> A = getInputMatrix();
		dbg.add("A", Rcpp::wrap(A));
		if (0) {
			// regularize internal representation
			Eigen::SparseMatrix<double> fAcopy = IAF.transpose();
			dbg.add("filteredA", Rcpp::wrap(fAcopy));
		}
		Eigen::SparseMatrix<double> fullSymS = fullS.selfadjointView<Eigen::Lower>();
		dbg.add("S", Rcpp::wrap(fullSymS));
		dbg.add("latentFilter", Rcpp::wrap(latentFilter));
		SEXP dv = Rcpp::wrap(dataVec);
		Rf_protect(dv);
		Rf_setAttrib(dv, R_NamesSymbol, obsNameVec);
		dbg.add("dataVec", dv);
	}

	void state::exportInternalState(MxRList &dbg)
	{
		dbg.add("rampartUsage", Rcpp::wrap(rampartUsage));
		dbg.add("macroA", Rcpp::wrap(macroA));

		SEXP modelName, key, numJoins, numKids, parent1, fk1, rscale, group, copy;
		//SEXP startLoc, endLoc, obsStart, obsEnd;
		Rf_protect(modelName = Rf_allocVector(STRSXP, layout.size()));
		Rf_protect(key = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numKids = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numJoins = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(parent1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(fk1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(rscale = Rf_allocVector(REALSXP, layout.size()));
		Rf_protect(group = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(copy = Rf_allocVector(INTSXP, layout.size()));
		//Rf_protect(startLoc = Rf_allocVector(INTSXP, layout.size()));
		//Rf_protect(endLoc = Rf_allocVector(INTSXP, layout.size()));
		//Rf_protect(obsStart = Rf_allocVector(INTSXP, layout.size()));
		//Rf_protect(obsEnd = Rf_allocVector(INTSXP, layout.size()));
		for (size_t mx=0; mx < layout.size(); ++mx) {
			SET_STRING_ELT(modelName, mx, Rf_mkChar(layout[mx].modelName().c_str()));
			INTEGER(key)[mx] = layout[mx].key;
			INTEGER(numKids)[mx] = layout[mx].numKids;
			INTEGER(numJoins)[mx] = layout[mx].numJoins;
			INTEGER(parent1)[mx] = plusOne(layout[mx].parent1);
			INTEGER(fk1)[mx] = layout[mx].fk1;
			REAL(rscale)[mx] = layout[mx].rampartScale;
			INTEGER(group)[mx] = layout[mx].group? layout[mx].group : NA_INTEGER;
			INTEGER(copy)[mx] = layout[mx].copy? layout[mx].copy : NA_INTEGER;
			/*
			INTEGER(startLoc)[mx] = 1 + layout[mx].modelStart;
			INTEGER(endLoc)[mx] = layout[mx].modelStart + layout[mx].numVars();
			if (layout[mx].numObs()) {
				INTEGER(obsStart)[mx] = 1+layout[mx].obsStart;
				INTEGER(obsEnd)[mx] = layout[mx].obsStart+layout[mx].numObs();
			} else {
				INTEGER(obsStart)[mx] = NA_INTEGER;
				INTEGER(obsEnd)[mx] = NA_INTEGER;
			}
			*/
		}
		dbg.add("layout", Rcpp::DataFrame::create(Rcpp::Named("model")=modelName,
							  Rcpp::Named("key")=key,
							  Rcpp::Named("numKids")=numKids,
							  Rcpp::Named("numJoins")=numJoins,
							  Rcpp::Named("parent1")=parent1,
							  Rcpp::Named("fk1")=fk1,
							  Rcpp::Named("rampartScale")=rscale,
							  Rcpp::Named("group")=group,
							  Rcpp::Named("copy")=copy));
							  //Rcpp::Named("modelStart")=startLoc,
							  //Rcpp::Named("modelEnd")=endLoc,
							  //Rcpp::Named("obsStart")=obsStart,
							  //Rcpp::Named("obsEnd")=obsEnd));
	}

};
