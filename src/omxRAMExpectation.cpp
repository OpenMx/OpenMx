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

static void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, 
    omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, 
    omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax);
static omxMatrix* omxGetRAMExpectationComponent(omxExpectation* ox, const char* component);

void omxRAMExpectation::ensureTrivialF() // move to R side? TODO
{
	omxRecompute(F, NULL);  // should not do anything
	EigenMatrixAdaptor eF(F);
	Eigen::MatrixXd ident(F->rows, F->rows);
	ident.setIdentity();
	if (ident != eF.block(0, 0, F->rows, F->rows)) {
		Rf_error("Square part of F matrix is not trivial");
	}
}

static void omxCallRAMExpectation(omxExpectation* oo, FitContext *fc, const char *what, const char *how)
{
	omxRAMExpectation* oro = (omxRAMExpectation*)(oo->argStruct);

	if (what && strEQ(what, "distribution") && how && strEQ(how, "flat")) {
		if (!oro->rram) {
			oro->rram = new RelationalRAMExpectation::state;
			oro->rram->init(oo, fc);
		}
		oro->rram->compute(fc);
		return;
	}

	omxRecompute(oro->A, fc);
	omxRecompute(oro->S, fc);
	omxRecompute(oro->F, fc);
	if(oro->M != NULL)
	    omxRecompute(oro->M, fc);
	    
	omxCalculateRAMCovarianceAndMeans(oro->A, oro->S, oro->F, oro->M, oro->cov, 
		oro->means, oro->numIters, oro->I, oro->Z, oro->Y, oro->X, oro->Ax);
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
	omxFreeMatrix(argStruct->Z);
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
	omxMatrix* Z = oro->Z;
	omxMatrix* I = oro->I;
    int numIters = oro->numIters;
    
    omxRecompute(A, NULL);
    omxRecompute(S, NULL);
	
    omxShallowInverse(NULL, numIters, A, Z, Ax, I ); // Z = (I-A)^-1

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

	if (oro->rram && oro->rram->fullCov.nonZeros()) {
		RelationalRAMExpectation::state *rram = oro->rram;
		SEXP m1 = Rcpp::wrap(rram->regularA.out.transpose() * rram->fullMeans);
		Rf_protect(m1);
		Rf_setAttrib(m1, R_NamesSymbol, rram->obsNameVec);
		out.add("mean", m1);
		out.add("covariance", Rcpp::wrap(rram->fullCov));
		rram->exportInternalState(dbg);
	}

	Rf_setAttrib(robj, Rf_install("output"), out.asR());
	Rf_setAttrib(robj, Rf_install("debug"), dbg.asR());
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

static void omxCalculateRAMCovarianceAndMeans(omxMatrix* A, omxMatrix* S, omxMatrix* F, 
	omxMatrix* M, omxMatrix* Cov, omxMatrix* Means, int numIters, omxMatrix* I, 
	omxMatrix* Z, omxMatrix* Y, omxMatrix* X, omxMatrix* Ax) {
	
	if (F->rows == 0) return;

	if(OMX_DEBUG) { mxLog("Running RAM computation with numIters is %d\n.", numIters); }
		
	if(Ax == NULL || I == NULL || Z == NULL || Y == NULL || X == NULL) {
		Rf_error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
	}
		
	if(Cov == NULL && Means == NULL) {
		return; // We're not populating anything, so why bother running the calculation?
	}
	
	omxShallowInverse(NULL, numIters, A, Z, Ax, I );
	
	/* Cov = FZSZ'F' */
	omxDGEMM(FALSE, FALSE, 1.0, F, Z, 0.0, Y);

	omxDGEMM(FALSE, FALSE, 1.0, Y, S, 0.0, X);

	omxDGEMM(FALSE, TRUE, 1.0, X, Y, 0.0, Cov);
	 // Cov = FZSZ'F' (Because (FZ)' = Z'F')
	
	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Cov, "....RAM: Model-implied Covariance Matrix:");}
	
	if(M != NULL && Means != NULL) {
		// F77_CALL(omxunsafedgemv)(Y->majority, &(Y->rows), &(Y->cols), &oned, Y->data, &(Y->leading), M->data, &onei, &zerod, Means->data, &onei);
		omxDGEMV(FALSE, 1.0, Y, M, 0.0, Means);
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(Means, "....RAM: Model-implied Means Vector:");}
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
	
	omxRAMExpectation *RAMexp = new omxRAMExpectation;
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

	RAMexp->Z = 	omxInitMatrix(k, k, TRUE, currentState);
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

	// verify whether sparse can deal with parameters set to exactly zero TODO

	void state::refreshLevelTransitions(FitContext *fc, addr &a1, Amatrix &dest, double scale)
	{
		if (scale == 0.0) return;
		omxExpectation *expectation = a1.model;
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		for (size_t jx=0; jx < ram->between.size(); ++jx) {
			omxMatrix *betA = ram->between[jx];
			int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
			if (key == NA_INTEGER) continue;
			omxData *data1 = betA->getJoinModel()->data;
			int frow = data1->lookupRowOfKey(key);
			int jOffset = data1->rowToOffsetMap[frow];
			omxRecompute(betA, fc);
			omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel()->argStruct;
			for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
				for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
					double val = omxMatrixElement(betA, rx, cx);
					if (val == 0.0) continue;
					dest.in.coeffRef(jOffset + cx, a1.modelStart + rx) =
						signA * val * scale;
					//mxLog("A(%d,%d) = %f", jOffset + cx, lx + rx, val);
				}
			}
		}
	}

	void state::refreshUnitA(FitContext *fc, addr &a1, Amatrix &dest)
	{
		omxExpectation *expectation = a1.model;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		EigenMatrixAdaptor eA(ram->A);
		for (int cx=0; cx < eA.cols(); ++cx) {
			for (int rx=0; rx < eA.rows(); ++rx) {
				double val = eA(rx, cx);
				if (val != 0) {
					if (rx == cx) {
						Rf_error("%s: nonzero diagonal entry in A matrix at %d",
							 homeEx->name, 1+a1.modelStart+cx);
					}
					// can't use eA.block(..) -= because fullA must remain sparse
					dest.in.coeffRef(a1.modelStart + cx, a1.modelStart + rx) = signA * val;
				}
			}
		}
	}

	// 3rd visitor
	void state::refreshModel(FitContext *fc)
	{
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			omxExpectation *expectation = a1.model;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
			omxData *data = expectation->data;
			data->handleDefinitionVarList(expectation->currentState, a1.row);
			omxRecompute(ram->A, fc);
			omxRecompute(ram->S, fc);

			refreshLevelTransitions(fc, a1, regularA, 1.0);
			if (rampartUsage.size()) {
				refreshLevelTransitions(fc, a1, rampartA, a1.rampartScale);
			}

			if (!haveFilteredAmat) {
				refreshUnitA(fc, a1, regularA);
				if (rampartUsage.size()) {
					refreshUnitA(fc, a1, rampartA);
				}
			}

			EigenMatrixAdaptor eS(ram->S);
			for (int cx=0; cx < eS.cols(); ++cx) {
				for (int rx=cx; rx < eS.rows(); ++rx) {
					if (eS(rx,cx) != 0) {
						fullS.coeffRef(a1.modelStart + rx, a1.modelStart + cx) = eS(rx, cx);
					}
				}
			}

			if (ram->M) {
				omxRecompute(ram->M, fc);
				EigenVectorAdaptor eM(ram->M);
				for (int mx=0; mx < eM.size(); ++mx) {
					fullMeans[a1.modelStart + mx] = eM[mx];
				}
			} else {
				fullMeans.segment(a1.modelStart, ram->A->cols).setZero();
			}
		}
	}

	// 1st visitor
	int state::placeOneRow(omxExpectation *expectation, int frow, int &totalObserved, int &maxSize)
	{
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;

		struct addr a1;
		a1.rampartScale = 1.0;
		a1.parent1 = NA_INTEGER;
		a1.fk1 = NA_INTEGER;
		int parent1Pos = -1;
		for (size_t jx=0; jx < ram->between.size(); ++jx) {
			omxMatrix *b1 = ram->between[jx];
			int key = omxKeyDataElement(data, frow, b1->getJoinKey());
			if (key == NA_INTEGER) continue;
			if (jx==0) a1.fk1 = key;
			AmatDependsOnParameters |= b1->dependsOnParameters();
			omxExpectation *e1 = b1->getJoinModel();
			int parentPos = placeOneRow(e1, e1->data->lookupRowOfKey(key), totalObserved, maxSize);
			if (jx==0) parent1Pos = parentPos;
		}

		a1.row = frow;
		a1.key = frow;
		if (data->hasPrimaryKey()) {
			if (data->rowToOffsetMap.size() == 0) {
				ram->ensureTrivialF();
			}

			// insert_or_assign would be nice here
			std::map<int,int>::const_iterator it = data->rowToOffsetMap.find(frow);
			if (it != data->rowToOffsetMap.end()) return it->second;

			data->rowToOffsetMap[frow] = maxSize;
			a1.key = data->primaryKeyOfRow(frow);
		}

		if (parent1Pos >= 0) {
			std::vector<addr>::iterator low =
				std::lower_bound(layout.begin(), layout.end(), parent1Pos, addr::CompareWithModelStart);
			if (parent1Pos != low->modelStart) Rf_error("Parent search failed"); //impossible
			low->numKids += 1;
			a1.parent1 = low - layout.begin();
		}
		a1.model = expectation;
		a1.numKids = 0;
		a1.numJoins = ram->between.size();
		a1.modelStart = maxSize;
		a1.modelEnd = maxSize + ram->F->cols - 1;
		a1.obsStart = totalObserved;

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

		a1.obsEnd = totalObserved - 1;
		layout.push_back(a1);
		if (verbose() >= 2) {
			if (a1.obsStart <= a1.obsEnd) {
				mxLog("place %s[%d] at %d %d obs %d %d", a1.modelName().c_str(),
				      frow, a1.modelStart, a1.modelEnd, a1.obsStart, a1.obsEnd);
			} else {
				mxLog("place latent %s[%d] at %d %d", a1.modelName().c_str(),
				      frow, a1.modelStart, a1.modelEnd);
			}
		}

		maxSize += ram->F->cols;
		AmatDependsOnParameters |= ram->A->dependsOnParameters();
		return a1.modelStart;
	}

	// 2nd visitor
	void state::examineModel()
	{
		int dx = 0;
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			omxExpectation *expectation = a1.model;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
			omxData *data = expectation->data;

			if (Global->RAMInverseOpt || ram->rampartEnabled()) {
				data->handleDefinitionVarList(expectation->currentState, a1.row);
				omxRAMExpectation *ram = (omxRAMExpectation*) expectation->argStruct;
				omxRecompute(ram->A, NULL);

				for (size_t jx=0; jx < ram->between.size(); ++jx) {
					omxMatrix *betA = ram->between[jx];
					int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
					if (key == NA_INTEGER) continue;
					omxData *data1 = betA->getJoinModel()->data;
					int frow = data1->lookupRowOfKey(key);
					int jOffset = data1->rowToOffsetMap[frow];
					omxRecompute(betA, NULL);
					betA->markPopulatedEntries();
					omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel()->argStruct;
					for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
						for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
							double val = omxMatrixElement(betA, rx, cx);
							if (val == 0.0) continue;
							testA.in.coeffRef(jOffset + cx, a1.modelStart + rx) = 1;
						}
					}
				}
				ram->A->markPopulatedEntries();
				EigenMatrixAdaptor eA(ram->A);
				for (int cx=0; cx < eA.cols(); ++cx) {
					for (int rx=0; rx < eA.rows(); ++rx) {
						if (rx != cx && eA(rx,cx) != 0) {
							testA.in.coeffRef(a1.modelStart + cx, a1.modelStart + rx) = 1;
						}
					}
				}
				// only care about diagonal for rampart homogeneous error variance test TODO
				EigenMatrixAdaptor eS(ram->S);
				for (int cx=0; cx < eS.cols(); ++cx) {
					for (int rx=cx; rx < eS.rows(); ++rx) {
						if (eS(rx,cx) != 0) {
							fullS.coeffRef(a1.modelStart + rx, a1.modelStart + cx) = 1;
						}
					}
				}
			}

			std::string modelName(expectation->data->name);
			modelName = modelName.substr(0, modelName.size() - 4); // remove "data" suffix

			int jCols = expectation->dataColumns->cols;
			if (jCols) {
				omxDataRow(expectation, a1.row, smallCol);
				omxMatrix *colList = expectation->dataColumns;
				for (int col=0; col < jCols; ++col) {
					double val = omxMatrixElement(smallCol, 0, col);
					bool yes = std::isfinite(val);
					if (!yes) continue;
					latentFilter[ a1.modelStart + col ] = true;
					std::string dname = modelName + omxDataColumnName(expectation->data,
											  omxVectorElement(colList, col));
					SET_STRING_ELT(obsNameVec, dx, Rf_mkChar(dname.c_str()));
					dataVec[ dx++ ] = val;
				}
			}
			for (int vx=0; vx < ram->F->cols; ++vx) {
				std::string dname = modelName + ram->F->colnames[vx];
				SET_STRING_ELT(varNameVec, a1.modelStart + vx, Rf_mkChar(dname.c_str()));
			}
		}
	}

	struct RampartCompareLib {
		state *st;
		RampartCompareLib(state *st) : st(st) {};

		omxExpectation *getJoinModel(const addr *a1) const {
			omxRAMExpectation *ram = (omxRAMExpectation*) a1->model->argStruct;
			omxMatrix *b1 = ram->between[0];
			return b1->getJoinModel();
		};

		bool cmpUpper(const addr *lhs, const addr *rhs) const
		{
			if (getJoinModel(lhs) != getJoinModel(rhs))
				return strcmp(getJoinModel(lhs)->name, getJoinModel(rhs)->name) < 0;
			return false;
		}

		bool cmp1(const addr *lhs, const addr *rhs) const
		{
			if (lhs->model != rhs->model)
				return strcmp(lhs->model->name, rhs->model->name) < 0;
			if (lhs->numObs() != rhs->numObs())
				return lhs->numObs() < rhs->numObs();
			if (lhs->numVars() != rhs->numVars())
				return lhs->numVars() < rhs->numVars();

			for (int lx=0; lx < lhs->numVars(); ++lx) {
				bool p1 = st->latentFilter[ lhs->modelStart + lx ];
				bool p2 = st->latentFilter[ rhs->modelStart + lx ];
				if (p1 == p2) continue;
				return int(p1) < int(p2);
			}
			if (lhs->clump.size() != rhs->clump.size())
				return lhs->clump.size() < rhs->clump.size();
			for (size_t cx=0; cx < lhs->clump.size(); ++cx) {
				return cmp1(&st->layout[lhs->clump[cx]], &st->layout[rhs->clump[cx]]);
			}
			return false;
		}
	};

	struct RampartTodoCompare : RampartCompareLib {
		RampartTodoCompare(state *st) : RampartCompareLib(st) {};

		bool operator() (const addr *lhs, const addr *rhs) const
		{
			if (cmpUpper(lhs, rhs)) return true;
			if (lhs->fk1 != rhs->fk1)
				return lhs->fk1 < rhs->fk1;
			return cmp1(lhs, rhs);
		}
	};

	struct RampartClumpCompare : RampartCompareLib {
		RampartClumpCompare(state *st) : RampartCompareLib(st) {};

		bool operator() (const int lhs, const int rhs) const
		{
			const addr *lhsObj = &st->layout[lhs];
			const addr *rhsObj = &st->layout[rhs];
			if (cmp1(lhsObj, rhsObj)) return true;
			return lhs < rhs;
		}
	};

	template <typename T>
	void state::oertzenRotate(std::vector<T> &t1)
	{
		// get covariance and sort by mahalanobis distance TODO
		addr &specimen = layout[ t1[0] ];
		for (int ox=0; ox < specimen.numObs(); ++ox) {
			std::vector<int> tmp;
			for (size_t ux=0; ux < t1.size(); ++ux) {
				addr &a1 = layout[ t1[ux] ];
				tmp.push_back(a1.obsStart + ox);
			}
			rotationPlan.push_back(tmp);
		}
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
			if (a1.numKids != 0 || a1.numJoins != 1 || a1.rampartScale != 1.0) continue;

			omxRAMExpectation *ram = (omxRAMExpectation*) a1.model->argStruct;
			omxMatrix *b1 = ram->between[0];
			if (b1->dependsOnDefinitionVariables()) continue;
			std::vector<int> &t1 = todo[&a1];
			t1.push_back(int(ax));
		}
		for (RampartMap::iterator it = todo.begin(); it != todo.end(); ++it) {
			std::vector<int> &t1 = it->second;
			if (t1.size() <= 1) continue;

			if (true) {
				if (maxRotationUnits < t1.size())
					maxRotationUnits = t1.size();
				std::string buf = "rotate units";
				oertzenRotate(t1);
				addr &specimen = layout[ t1[0] ];
				specimen.rampartScale = sqrt(double(t1.size()));
				if (false) {
				for (size_t cx=0; cx < specimen.clump.size(); ++cx) {
					addr &a1 = layout[ specimen.clump[cx] ];
					double orig = a1.rampartScale*a1.rampartScale;
					a1.rampartScale = sqrt(orig / double(t1.size())) * sqrt(orig);
				}
				}
				layout[specimen.parent1].numKids -= t1.size();
				layout[specimen.parent1].clump.push_back(t1[0]);
				for (size_t ux=1; ux < t1.size(); ++ux) {
					addr &a1 = layout[ t1[ux] ];
					a1.rampartScale = 0;
					a1.numJoins = 0;
					buf += string_snprintf(" %d", 1+t1[ux]);
				}
				buf += string_snprintf(" -> %d\n", 1+t1[0]);
				//mxLogBig(buf);
			} else {
				// Don't rotate, just clump units together with parent.
				addr &specimen = layout[ t1[0] ];
				layout[specimen.parent1].numKids -= t1.size();
				layout[specimen.parent1].clump.insert(layout[specimen.parent1].clump.end(),
								      t1.begin(), t1.end());
				for (size_t ux=0; ux < t1.size(); ++ux) {
					addr &a1 = layout[ t1[ux] ];
					a1.numJoins = 0;
				}
			}
			unlinked += t1.size() - 1;
		}
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr& a1 = layout[ax];
			if (a1.clump.size() <= 1) continue;
			std::sort(a1.clump.begin(), a1.clump.end(), RampartClumpCompare(this));
		}
		return unlinked;
	}

	void state::init(omxExpectation *expectation, FitContext *fc)
	{
		homeEx = expectation;

		omxRAMExpectation *ram = (omxRAMExpectation*) homeEx->argStruct;
		ram->ensureTrivialF();

		int numManifest = ram->F->rows;

		AmatDependsOnParameters = ram->A->dependsOnParameters();
		haveFilteredAmat = false;
		smallCol = omxInitMatrix(1, numManifest, TRUE, homeEx->currentState);
		omxData *data               = homeEx->data;

		// don't permit reuse of our expectations by some other fit function TODO

		int totalObserved = 0;
		int maxSize = 0;
		for (int row=0; row < data->rows; ++row) {
			placeOneRow(homeEx, row, totalObserved, maxSize);
		}

		if (verbose() >= 1) {
			mxLog("%s: total observations %d AmatDependsOnParameters=%d",
			      homeEx->name, totalObserved, AmatDependsOnParameters);
		}
		latentFilter.assign(maxSize, false); // will have totalObserved true entries
		obsNameVec = Rf_protect(Rf_allocVector(STRSXP, totalObserved));
		varNameVec = Rf_protect(Rf_allocVector(STRSXP, maxSize));
		dataVec.resize(totalObserved);
		testA.in.resize(maxSize, maxSize);
		ident.resize(maxSize, maxSize);
		ident.setIdentity();

		FreeVarGroup *varGroup = Global->findVarGroup(FREEVARGROUP_ALL); // ignore freeSet

		if (Global->RAMInverseOpt) {
			Eigen::VectorXd vec(varGroup->vars.size());
			vec.setConstant(1);
			copyParamToModelInternal(varGroup, homeEx->currentState, vec.data());
		}

		fullS.conservativeResize(maxSize, maxSize);
		examineModel();

		if (OMX_DEBUG) {
			int lfCount = std::count(latentFilter.begin(), latentFilter.end(), true);
			if (lfCount != totalObserved) {
				Rf_error("lfCount %d != totalObserved %d", lfCount, totalObserved);
			}
		}

		AshallowDepth = -1;

		if (Global->RAMInverseOpt) {
			int maxDepth = std::min(maxSize, 30);
			if (Global->RAMMaxDepth != NA_INTEGER) maxDepth = Global->RAMMaxDepth;
			Eigen::SparseMatrix<double> curProd = testA.in;
			for (int tx=1; tx < maxDepth; ++tx) {
				if (verbose() >= 3) { Eigen::MatrixXd tmp = curProd; mxPrintMat("curProd", tmp); }
				curProd = (curProd * testA.in.transpose()).eval();
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
			fc->copyParamToModelClean();
		}
		signA = AshallowDepth >= 0 ? 1.0 : -1.0;
		if (verbose() >= 1) {
			mxLog("%s: RAM shallow inverse depth = %d", homeEx->name, AshallowDepth);
		}

		maxRotationUnits = 0;
		if (ram->rampartEnabled()) {
			invertAndFilterA(testA);
			for (size_t ax=0; ax < layout.size(); ++ax) {
				addr &a1 = layout[ax];
				omxRAMExpectation *ram = (omxRAMExpectation*) a1.model->argStruct;

				if (0) {
					Eigen::ArrayXd sdiag = fullS.diagonal();
					for (int cx=a1.obsStart; cx <= a1.obsEnd; ++cx) {
						Eigen::VectorXd loadings = testA.out.col(cx);
						Eigen::ArrayXd contrib = loadings.array() * sdiag;
						for (int vx=0; vx < maxSize; ++vx) {
							if (vx == cx || contrib[vx] == 0.0) continue;
							mxLog("%s contributes variance to %s",
							      CHAR(STRING_ELT(varNameVec, vx)),
							      CHAR(STRING_ELT(obsNameVec, cx)));
						}
					}
				}
			}

			int maxIter = ram->rampart;
			int unlinked = 0;
			int level = -1;
			while (int more = rampartRotate(++level)) {
				rampartUsage.push_back(more);
				unlinked += more;
				if (--maxIter == 0) break;
			}
			if (verbose() >= 1) {
				mxLog("%s: rampart unlinked %d units", homeEx->name, unlinked);
			}
		}
		fullS.setZero();

		ProtectedSEXP RscaleOverride(R_do_slot(expectation->rObj, Rf_install("scaleOverride")));
		if (Rf_length(RscaleOverride)) {
			double *override = REAL(RscaleOverride);
			for (int ox=0; ox < Rf_length(RscaleOverride); ox += 2) {
				layout[ override[ox] - 1 ].rampartScale = override[ox+1];
			}
		}
	}

	state::~state()
	{
		omxFreeMatrix(smallCol);
	}

	void state::invertAndFilterA(Amatrix &Amat)
	{
		// consider http://users.clas.ufl.edu/hager/papers/Lightning/update.pdf ?
		if (AshallowDepth >= 0) {
			Amat.in.makeCompressed();
			Amat.out = Amat.in + ident;
			for (int iter=1; iter <= AshallowDepth; ++iter) {
				Amat.out = (Amat.out * Amat.in + ident).eval();
				//{ Eigen::MatrixXd tmp = Amat.out; mxPrintMat("Amat.out", tmp); }
			}
		} else {
			Amat.in += ident;
			if (!Amat.analyzed) {
				Amat.analyzed = true;
				Amat.in.makeCompressed();
				Amat.solver.analyzePattern(Amat.in);
			}
			Amat.solver.factorize(Amat.in);
			if (Amat.solver.info() != Eigen::Success) {
				Rf_error("%s: failed to invert flattened A matrix; %s",
					 homeEx->name, Amat.solver.lastErrorMessage().c_str());
			}

			Amat.out = Amat.solver.solve(ident);
			Amat.in -= ident;  // leave unchanged
			//{ Eigen::MatrixXd tmp = Amat.out; mxPrintMat("Amat.out", tmp); }
		}

		const bool doubleCheck = true;
		Eigen::MatrixXd denseA;
		if (doubleCheck) {
			denseA = Amat.out;
		}

		// We built A transposed so we can quickly filter columns
		// Switch to filterOuter http://eigen.tuxfamily.org/bz/show_bug.cgi?id=1130 TODO
		Amat.out.uncompress();
		Eigen::SparseMatrix<double>::Index *op = Amat.out.outerIndexPtr();
		Eigen::SparseMatrix<double>::Index *nzp = Amat.out.innerNonZeroPtr();
		int dx = 0;
		for (int cx=0; cx < Amat.in.cols(); ++cx) {
			if (!latentFilter[cx]) continue;
			op[dx] = op[cx];
			nzp[dx] = nzp[cx];
			++dx;
		}
		op[dx] = op[Amat.in.cols()];
		Amat.out.conservativeResize(Amat.in.rows(), dataVec.size());

		if (doubleCheck) {
			Eigen::MatrixXd denseAF;
			denseAF.resize(Amat.in.rows(), dataVec.size());
			int dx=0;
			for (int cx=0; cx < Amat.in.cols(); ++cx) {
				if (!latentFilter[cx]) continue;
				denseAF.col(dx) = denseA.col(cx);
				++dx;
			}
			if (dx != dataVec.size()) Rf_error("latentFilter has wrong count %d != %d",
							   dx, dataVec.size());
			Eigen::MatrixXd denseFilteredA = Amat.out;
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
		//{ Eigen::MatrixXd tmp = Amat.out; mxPrintMat("Amat.out", tmp); }
	}

	void state::compute(FitContext *fc)
	{
		fullMeans.conservativeResize(latentFilter.size());

		if (regularA.in.nonZeros() == 0) {
			regularA.in.resize(latentFilter.size(), latentFilter.size());
			if (rampartUsage.size()) {
				rampartA.in.resize(latentFilter.size(), latentFilter.size());
			}
		}
		fullS.conservativeResize(latentFilter.size(), latentFilter.size());

		refreshModel(fc);

		//{ Eigen::MatrixXd tmp = fullA; mxPrintMat("fullA", tmp); }
		//{ Eigen::MatrixXd tmp = fullS; mxPrintMat("fullS", tmp); }

		if (!haveFilteredAmat) {
			invertAndFilterA(regularA);
			if (rampartUsage.size()) {
				invertAndFilterA(rampartA);
			}
			haveFilteredAmat = !AmatDependsOnParameters;
		}
		//mxPrintMat("S", fullS);
		//Eigen::MatrixXd fullCovDense =

		if (!rampartUsage.size()) {
			fullCov = (regularA.out.transpose() * fullS.selfadjointView<Eigen::Lower>() * regularA.out);
		} else {
			fullCov = (rampartA.out.transpose() * fullS.selfadjointView<Eigen::Lower>() * rampartA.out);
		}
		//mxLog("fullCov %d%% nonzero", int(fullCov.nonZeros() * 100.0 / (fullCov.rows() * fullCov.cols())));
		//{ Eigen::MatrixXd tmp = fullCov; mxPrintMat("fullcov", tmp); }
	}

	static int plusOne(int val) {
		if (val == NA_INTEGER) return val;
		else return val + 1;
	}

	void state::exportInternalState(MxRList &dbg)
	{
		SEXP fmean = Rcpp::wrap(fullMeans);
		dbg.add("mean", fmean);
		Rf_setAttrib(fmean, R_NamesSymbol, varNameVec);
		// output rampartA TODO
		Eigen::SparseMatrix<double> A = signA * regularA.in.transpose();
		dbg.add("A", Rcpp::wrap(A));
		if (rampartUsage.size()) {
			Eigen::SparseMatrix<double> rA = signA * rampartA.in.transpose();
			dbg.add("rA", Rcpp::wrap(rA));
		}
		if (0) {
			// regularize internal representation
			Eigen::SparseMatrix<double> fAcopy = regularA.out.transpose();
			dbg.add("filteredA", Rcpp::wrap(fAcopy));
		}
		Eigen::SparseMatrix<double> fullSymS = fullS.selfadjointView<Eigen::Lower>();
		dbg.add("rampartUsage", Rcpp::wrap(rampartUsage));
		dbg.add("S", Rcpp::wrap(fullSymS));
		dbg.add("latentFilter", Rcpp::wrap(latentFilter));
		SEXP dv = Rcpp::wrap(dataVec);
		Rf_protect(dv);
		Rf_setAttrib(dv, R_NamesSymbol, obsNameVec);
		dbg.add("dataVec", dv);

		SEXP modelName, key, numJoins, numKids, parent1, fk1,
			startLoc, endLoc, obsStart, obsEnd, rscale;
		Rf_protect(modelName = Rf_allocVector(STRSXP, layout.size()));
		Rf_protect(key = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numKids = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numJoins = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(parent1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(fk1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(rscale = Rf_allocVector(REALSXP, layout.size()));
		Rf_protect(startLoc = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(endLoc = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(obsStart = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(obsEnd = Rf_allocVector(INTSXP, layout.size()));
		for (size_t mx=0; mx < layout.size(); ++mx) {
			SET_STRING_ELT(modelName, mx, Rf_mkChar(layout[mx].modelName().c_str()));
			INTEGER(key)[mx] = layout[mx].key;
			INTEGER(numKids)[mx] = layout[mx].numKids;
			INTEGER(numJoins)[mx] = layout[mx].numJoins;
			INTEGER(parent1)[mx] = plusOne(layout[mx].parent1);
			INTEGER(fk1)[mx] = layout[mx].fk1;
			REAL(rscale)[mx] = layout[mx].rampartScale;
			INTEGER(startLoc)[mx] = 1 + layout[mx].modelStart;
			INTEGER(endLoc)[mx] = 1 + layout[mx].modelEnd;
			if (layout[mx].obsStart <= layout[mx].obsEnd) {
				INTEGER(obsStart)[mx] = 1+layout[mx].obsStart;
				INTEGER(obsEnd)[mx] = 1+layout[mx].obsEnd;
			} else {
				INTEGER(obsStart)[mx] = NA_INTEGER;
				INTEGER(obsEnd)[mx] = NA_INTEGER;
			}
		}
		dbg.add("layout", Rcpp::DataFrame::create(Rcpp::Named("model")=modelName,
							  Rcpp::Named("key")=key,
							  Rcpp::Named("numKids")=numKids,
							  Rcpp::Named("numJoins")=numJoins,
							  Rcpp::Named("parent1")=parent1,
							  Rcpp::Named("fk1")=fk1,
							  Rcpp::Named("rampartScale")=rscale,
							  Rcpp::Named("modelStart")=startLoc,
							  Rcpp::Named("modelEnd")=endLoc,
							  Rcpp::Named("obsStart")=obsStart,
							  Rcpp::Named("obsEnd")=obsEnd));
	}

};
