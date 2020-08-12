/*
 * Copyright 2007-2019 by the individuals mentioned in the source code history
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

#include "omxExpectation.h"
#include "omxFitFunction.h"
#include "RAMInternal.h"
#include "Compute.h"
//#include <Eigen/LU>
#include "EnableWarnings.h"

void omxRAMExpectation::MpcIO::recompute(FitContext *fc)
{
	omxMatrix *M = !fc? M0 : fc->state->lookupDuplicate(M0);
	omxRecompute(M, fc);
}

unsigned omxRAMExpectation::MpcIO::getVersion(FitContext *fc)
{
	omxMatrix *M = !fc? M0 : fc->state->lookupDuplicate(M0);
	return omxGetMatrixVersion(M);
}

void omxRAMExpectation::MpcIO::refresh(FitContext *fc)
{
	omxMatrix *M = !fc? M0 : fc->state->lookupDuplicate(M0);
	EigenVectorAdaptor Emean(M);
	full = Emean;
}

void omxRAMExpectation::ApcIO::recompute(FitContext *fc)
{
	omxMatrix *A = !fc? A0 : fc->state->lookupDuplicate(A0);
	omxRecompute(A, fc);
}

unsigned omxRAMExpectation::ApcIO::getVersion(FitContext *fc)
{
	omxMatrix *A = !fc? A0 : fc->state->lookupDuplicate(A0);
	return omxGetMatrixVersion(A);
}

template <typename T>
void omxRAMExpectation::ApcIO::
_refresh(FitContext *fc, T &mat, double sign)
{
	omxMatrix *A = !fc? A0 : fc->state->lookupDuplicate(A0);
	if (sign == 1) {
		for (auto &v : vec) mat.coeffRef(v.c, v.r) = A->data[v.off];
	} else {
		for (auto &v : vec) mat.coeffRef(v.c, v.r) = -A->data[v.off];
	}
}

void omxRAMExpectation::SpcIO::recompute(FitContext *fc)
{
	omxMatrix *S = !fc? S0 : fc->state->lookupDuplicate(S0);
	omxRecompute(S, fc);
}

unsigned omxRAMExpectation::SpcIO::getVersion(FitContext *fc)
{
	omxMatrix *S = !fc? S0 : fc->state->lookupDuplicate(S0);
	return omxGetMatrixVersion(S);
}

template <typename T>
void omxRAMExpectation::SpcIO::_refresh(FitContext *fc, T &mat)
{
	omxMatrix *S = !fc? S0 : fc->state->lookupDuplicate(S0);
	omxRecompute(S, fc);
	for (auto &v : vec) mat.coeffRef(v.r, v.c) = S->data[v.off];
}

void omxRAMExpectation::flatten(FitContext *fc)
{
	if (rram) return;
	rram = new RelationalRAMExpectation::state;
	rram->init(this, fc);
}

void omxRAMExpectation::getExogenousPredictors(std::vector<int> &out)
{
	out = exoDataColumns;
}

void omxRAMExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	omxRAMExpectation* oro = this;

	if (what && how && strEQ(how, "flat")) {
		bool wantCov = false;
		bool wantMean = false;
		if (strEQ(what, "distribution")) { wantCov = true; wantMean = true; }
		if (strEQ(what, "covariance")) wantCov = true;
		if (strEQ(what, "mean")) wantMean = true;
		flatten(fc);
		if (wantCov)  oro->rram->computeCov(fc);
		if (wantMean) oro->rram->computeMean(fc);
		return;
	}

	oro->CalculateRAMCovarianceAndMeans(fc);

	super::compute(fc, what, how);
}

void omxRAMExpectation::invalidateCache()
{
  super::invalidateCache();
	if (rram) {
		delete rram;
		rram = 0;
	}
}

omxRAMExpectation::~omxRAMExpectation()
{
	if (rram) delete rram;
}

void omxRAMExpectation::populateAttr(SEXP robj)
{
	{
		ProtectedSEXP expCovExt(Rf_allocMatrix(REALSXP, A->rows, A->cols));
		Eigen::Map<Eigen::MatrixXd> covMem(REAL(expCovExt), A->rows, A->cols);
		pcalc.fullCov(0, covMem);

		Rf_setAttrib(robj, Rf_install("UnfilteredExpCov"), expCovExt);
		if (F->colnames.size()) {
			ProtectedSEXP names(Rf_allocVector(STRSXP, A->rows));
			for (int px=0; px < A->rows; ++px) {
				SET_STRING_ELT(names, px, Rf_mkChar(F->colnames[px]));
			}
			ProtectedSEXP dimnames(Rf_allocVector(VECSXP, 2));
			for (int dx=0; dx < 2; ++dx) SET_VECTOR_ELT(dimnames, dx, names);
			Rf_setAttrib(expCovExt, R_DimNamesSymbol, dimnames);
		}

		ProtectedSEXP RnumStats(Rf_ScalarReal(omxDataDF(data)));
		Rf_setAttrib(robj, Rf_install("numStats"), RnumStats);
	}

	MxRList out;
	MxRList dbg;

	omxRAMExpectation* oro = this;
	if (oro->rram) {
		rram->exportInternalState(dbg);
	} else {
    // This code should be shared with all normal expectations TODO
    compute(0, 0, 0);
		EigenMatrixAdaptor Ecov(oro->cov);
		out.add("covariance", Rcpp::wrap(Ecov));
		if (oro->means) {
			EigenVectorAdaptor Emean(oro->means);
			out.add("mean", Rcpp::wrap(Emean));
		}
		if (hasProductNodes) {
			dbg.add("polyRep", Rcpp::wrap(pcalc.getPolyRep()));
		}
    populateNormalAttr(robj, out);
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

void omxRAMExpectation::CalculateRAMCovarianceAndMeans(FitContext *fc)
{
	if (F->rows == 0) return;

	if (slope) omxRecompute(slope, fc);

	if(cov == NULL && means == NULL) {
		return; // We're not populating anything, so why bother running the calculation?
	}

	EigenMatrixAdaptor Ecov(cov);
	pcalc.cov(fc, Ecov);

	if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(cov, "....RAM: Model-implied Covariance Matrix:");}

	if (M) {
		EigenVectorAdaptor Emean(means);
		pcalc.mean(fc, Emean);
		if (slope) {
			EigenMatrixAdaptor Eslope(slope);
			Emean += Eslope * exoPredMean;
		}
		if(OMX_DEBUG_ALGEBRA) {omxPrintMatrix(means, "....RAM: Model-implied Means Vector:");}
	}
}

omxExpectation *omxInitRAMExpectation(omxState *st, int num)
{ return new omxRAMExpectation(st, num); }

static void recordNonzeroCoeff(omxMatrix *m, std::vector<coeffLoc> &vec, bool lowerTri)
{
	omxRecompute(m, 0);
	m->markPopulatedEntries();
	EigenMatrixAdaptor em(m);

	for (int cx=0; cx < m->cols; ++cx) {
		int rStart = lowerTri? cx : 0;
		for (int rx=rStart; rx < m->rows; ++rx) {
			if (em(rx,cx) == 0.0) continue;
			vec.emplace_back(&em.coeffRef(rx,cx) - &em.coeffRef(0,0), rx, cx);
		}
	}
}

omxRAMExpectation::omxRAMExpectation(omxState *st, int num)
	: super(st, num), studiedF(false), numExoPred(0), slope(0), rram(0)
{
	if (st->isClone()) {
		omxRAMExpectation *ram = (omxRAMExpectation *)st->getParent(this);
		Scoeff = &ram->ScoeffStorage;
		Acoeff = &ram->AcoeffStorage;
		pcalc.clone(ram->pcalc);
	} else {
		Scoeff = &ScoeffStorage;
		Acoeff = &AcoeffStorage;
	}
}

void omxRAMExpectation::init()
{
	loadDataColFromR();
	loadThresholdFromR();

	int l, k;

	omxRAMExpectation *RAMexp = this;
	RAMexp->rram = 0;

	auto oo=this;

	oo->canDuplicate = true;

	ProtectedSEXP Rverbose(R_do_slot(rObj, Rf_install("verbose")));
	RAMexp->verbose = Rf_asInteger(Rverbose) + OMX_DEBUG;

	ProtectedSEXP RsingleGroup(R_do_slot(rObj, Rf_install(".forceSingleGroup")));
	RAMexp->forceSingleGroup = Rf_asLogical(RsingleGroup);

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

	ProtectedSEXP Rrampart(R_do_slot(rObj, Rf_install(".rampartCycleLimit")));
	RAMexp->rampartCycleLimit = Rf_asInteger(Rrampart);

	ProtectedSEXP RrampartLimit(R_do_slot(rObj, Rf_install(".rampartUnitLimit")));
	RAMexp->rampartUnitLimit = Rf_asInteger(RrampartLimit);

	RAMexp->maxDebugGroups = 0;
	if (R_has_slot(rObj, Rf_install(".maxDebugGroups"))) {
		ProtectedSEXP Rmdg(R_do_slot(rObj, Rf_install(".maxDebugGroups")));
		RAMexp->maxDebugGroups = Rf_asInteger(Rmdg);
	}

	RAMexp->useSufficientSets = true;
	if (R_has_slot(rObj, Rf_install(".useSufficientSets"))) {
		ProtectedSEXP Rss(R_do_slot(rObj, Rf_install(".useSufficientSets")));
		RAMexp->useSufficientSets = Rf_asLogical(Rss);
	}

	RAMexp->optimizeMean = 0;
	if (R_has_slot(rObj, Rf_install(".optimizeMean"))) {
		ProtectedSEXP Rom(R_do_slot(rObj, Rf_install(".optimizeMean")));
		RAMexp->optimizeMean = Rf_asInteger(Rom);
	}

	useSparse = NA_INTEGER;
	if (R_has_slot(rObj, Rf_install(".useSparse"))) {
		ProtectedSEXP Rus(R_do_slot(rObj, Rf_install(".useSparse")));
		useSparse = Rf_asInteger(Rus);
	}

	ProtectedSEXP Rbetween(R_do_slot(rObj, Rf_install("between")));
	if (Rf_length(Rbetween)) {
		if (!oo->data) mxThrow("%s: data is required for joins", oo->name);
		if (!Rf_isInteger(Rbetween)) mxThrow("%s: between must be an integer vector", oo->name);
		RAMexp->between.reserve(Rf_length(Rbetween));
		int *bnumber = INTEGER(Rbetween);
		for (int jx=0; jx < Rf_length(Rbetween); ++jx) {
			omxMatrix *bmat = currentState->getMatrixFromIndex(bnumber[jx]);
			int foreignKey = bmat->getJoinKey();
			omxExpectation *fex = bmat->getJoinModel();
			if (!fex) mxThrow("%s: level transition matrix '%s' does not reference the upper level model",
					   oo->name, bmat->name());
			omxCompleteExpectation(fex);
			if (!strEQ(fex->name, "MxExpectationRAM")) {
				mxThrow("%s: only MxExpectationRAM can be joined with MxExpectationRAM (not %s)",
								oo->name, fex->name);
			}
			omxDataKeysCompatible(fex->data, oo->data, foreignKey);
			if (!omxDataColumnIsKey(oo->data, foreignKey)) {
				mxThrow("Cannot join using non-integer type column '%s' in '%s'. "
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
		mxThrow("RAM matrices '%s' and '%s' must have the same dimensions",
			 RAMexp->S->name(), RAMexp->A->name());
	}

	hasProductNodes = false;
	isProductNode.assign(A->cols, false);
	if (R_has_slot(rObj, Rf_install("isProductNode"))) {
		ProtectedSEXP RprodNode(R_do_slot(rObj, Rf_install("isProductNode")));
    // PPML can transform a model such that RprodNode length does
    // not match A->cols
    if (Rf_length(RprodNode) == A->cols) {
      for (int px = 0; px < A->cols; ++px) {
        if (INTEGER(RprodNode)[px]) {
          isProductNode[px] = true;
          hasProductNodes = true;
        }
      }
    }
	}

	if(OMX_DEBUG) { mxLog("Generating internals for computation."); }

	cov = omxNewMatrixFromSlotOrAnon(rObj, currentState, "expectedCovariance", l, l);
	if (!cov->hasMatrixNumber) covOwner = omxMatrixPtr(cov);
	else connectMatrixToExpectation(cov, this, "covariance");

	if (M) {
		means =	omxNewMatrixFromSlotOrAnon(rObj, currentState, "expectedMean", 1, l);
		if (!means->hasMatrixNumber) meanOwner = omxMatrixPtr(means);
		else connectMatrixToExpectation(means, this, "mean");
	} else {
		RAMexp->means  = 	NULL;
	}

	RAMexp->studyF();
	//mxPrintMat("RAM corrected dc", oo->getDataColumns());

	if (!currentState->isClone()) {
		if (hasProductNodes) rampartCycleLimit = 0;

		Eigen::VectorXd estSave;
		currentState->setFakeParam(estSave);
		loadFakeDefVars();

		MpcIO *mio = 0;
		if (M) {
			mio = new MpcIO;
			mio->M0 = M;
		}
		recordNonzeroCoeff(A, *Acoeff, false);
		ApcIO *aio = new ApcIO(*Acoeff);
		aio->A0 = A;
		recordNonzeroCoeff(S, *Scoeff, true);
		SpcIO *sio = new SpcIO(*Scoeff);
		sio->S0 = S;

		pcalc.attach(k, l, latentFilter, isProductNode, mio, aio, sio);
		pcalc.setAlgo(0, hasProductNodes, useSparse);

		currentState->restoreParam(estSave);
	} else {
    auto pex = (omxRAMExpectation*) currentState->getParent(this);
    if (pex->slope) {
      numExoPred = pex->numExoPred;
      exoDataColIndex = pex->exoDataColIndex;
      addSlopeMatrix();
    }
  }
}

void omxRAMExpectation::studyExoPred()
{
	if (data->defVars.size() == 0 || !M || !M->isSimple() || !S->isSimple()) return;

	Eigen::VectorXd estSave;
	currentState->setFakeParam(estSave);
	omxRecompute(A, 0);

	EigenMatrixAdaptor eA(A);
	EigenMatrixAdaptor eS(S);
	hasVariance = eS.diagonal().array().abs().matrix();

  exoDataColIndex.resize(S->rows, -1);
	int mNum = ~M->matrixNumber;
	for (int k=0; k < int(data->defVars.size()); ++k) {
		omxDefinitionVar &dv = data->defVars[k];
		if (dv.matrix == mNum && hasVariance[ dv.col ] == 0.0) {
			bool toManifest = false;
			const char *latentName=0;
			for (int cx=0; cx < eA.cols(); ++cx) {
				if (eA(cx, dv.col) == 0.0) continue;
				if (latentFilter[cx]) toManifest = true;
				else latentName = S->colnames[cx];
			}
			if (!toManifest && !latentName) continue;
			if (latentName) mxThrow("%s: latent exogenous variables are not supported (%s -> %s)", name,
						 S->colnames[dv.col], latentName);
			exoDataColIndex[dv.col] = dv.column;
			numExoPred += 1;
			dv.loadData(currentState, 0.);
			if (OMX_DEBUG + verbose >= 1) {
				mxLog("%s: set defvar '%s' for latent '%s' to exogenous mode",
				      name, data->columnName(dv.column), S->colnames[dv.col]);
			}
			data->defVars.erase(data->defVars.begin() + k--);
		}
	}

	currentState->restoreParam(estSave);

  addSlopeMatrix();
}

void omxRAMExpectation::addSlopeMatrix()
{
	if (!numExoPred) return;

	slope = omxInitMatrix(F->rows, numExoPred, currentState);
	EigenMatrixAdaptor eSl(slope);
	eSl.setZero();

	for (int cx=0, ex=0; cx < S->rows; ++cx) {
		if (exoDataColIndex[cx] == -1) continue;
		auto &rc = data->rawCol( exoDataColIndex[cx] );
		if (rc.type != COLUMNDATA_NUMERIC) {
			omxRaiseErrorf("%s: exogenous predictor '%s' must be type numeric (not '%s')",
				       name, rc.name, rc.typeName());
			continue;
		}
		exoDataColumns.push_back(exoDataColIndex[cx]);
		for (int rx=0, dx=0; rx < S->rows; ++rx) {
			if (!latentFilter[rx]) continue;
			slope->addPopulate(A, rx, cx, dx, ex);
			dx += 1;
		}
		ex += 1;
	}
}

void omxRAMExpectation::connectToData()
{
  super::connectToData();

	exoPredMean.resize(exoDataColumns.size());
	for (int cx=0; cx < int(exoDataColumns.size()); ++cx) {
		auto &e1 = data->rawCol( exoDataColumns[cx] );
		Eigen::Map< Eigen::VectorXd > vec(e1.d(), data->numRawRows());
		exoPredMean[cx] = vec.mean();
	}
}

void omxRAMExpectation::studyF()
{
	// Permute the data columns such that the manifest
	// part of F is a diagonal matrix. This permits
	// trivial filtering of the latent variables.
	auto dataColumns = super::getDataColumns();
	auto origDataColumnNames = super::getDataColumnNames();
	auto origThresholdInfo = super::getThresholdInfo();
	EigenMatrixAdaptor eF(F);
	latentFilter.assign(eF.cols(), false);
	dataCols.resize(eF.rows());
	dataColNames.resize(eF.rows(), 0);
	if (!eF.rows()) return;  // no manifests
	for (int cx =0, dx=0; cx < eF.cols(); ++cx) {
		int dest;
		double isManifest = eF.col(cx).maxCoeff(&dest);
		latentFilter[cx] = isManifest;
		if (isManifest) {
			dataColNames[dx] = origDataColumnNames[dest];
			int newDest = dataColumns.size()? dataColumns[dest] : dest;
			dataCols[dx] = newDest;
			if (origThresholdInfo.size()) {
				omxThresholdColumn adj = origThresholdInfo[dest];
				adj.dataColumn = newDest;
				thresholds.push_back(adj);
			}
			dx += 1;
		}
	}
  // ensure we don't try to use the wrong stuff
  origDataColumnNames.clear();
  origThresholdInfo.clear();
  studiedF = true;
  //mxLog("studyF: dataColumns permuted to (%d):", int(dataColNames.size()));
  //for (auto &dc : dataColNames) mxLog("%s", dc);
}

omxMatrix* omxRAMExpectation::getComponent(const char* component)
{
	if(OMX_DEBUG) { mxLog("RAM expectation: %s requested--", component); }

	omxRAMExpectation* ore = this;
	omxMatrix* retval = NULL;

	if(strEQ("cov", component)) {
		retval = ore->cov;
	} else if(strEQ("means", component)) {
		retval = ore->means;
	} else if(strEQ("slope", component)) {
		if (!slope) studyExoPred();
		retval = slope;
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}

	return retval;
}


void omxRAMExpectation::generateData(FitContext *fc, MxRList &out)
{
	if (!between.size()) {
		super::generateData(fc, out);
	}

	flatten(fc);

	rram->simulate(fc, out);
}

void omxRAMExpectation::analyzeDefVars(FitContext *fc)
{
	int numDefVars = data->defVars.size();

	data->loadFakeData(currentState, 1.0);

	hasMean.resize(S->rows);
	if (M && M->isSimple()) {
		omxRecompute(M, fc);
		EigenVectorAdaptor eM(M);
		hasMean = eM.array().abs();
		dvInfluenceMean.assign(numDefVars, false);
	} else {
		hasMean.setConstant(M? 1. : 0.);
		dvInfluenceMean.assign(numDefVars, M? true : false);
		if (verbose >= 1) mxLog("%s: defvar effect on mean unknown", name);
	}

	hasVariance.resize(S->rows);
	if (S->isSimple()) {
		omxRecompute(S, fc);
		EigenMatrixAdaptor eS(S);
		hasVariance = eS.diagonal().array().abs().matrix();
		dvInfluenceVar.assign(numDefVars, false);
	} else {
		hasVariance.setConstant(1.0);
		dvInfluenceVar.assign(numDefVars, true);
		if (verbose >= 1) mxLog("%s: defvar effect on variance unknown", name);
	}

	dvContribution.resize(S->rows);

	int sNum = ~S->matrixNumber;
	int mNum = M? ~M->matrixNumber : 0;

	std::set<int> tracked;
	if (!A->algebra) tracked.insert(~A->matrixNumber);
	for (auto mat : between) {
		if (!mat->algebra) tracked.insert(~mat->matrixNumber);
	}

	for (int k=0; k < numDefVars; ++k) {
		omxDefinitionVar &dv = data->defVars[k];

		if (M && dv.matrix == mNum) {
			dvInfluenceMean[k] = true;
			dvInfluenceVar[k] = dvInfluenceVar[k] | (hasVariance[ dv.col ] != 0.0);
			dvContribution[dv.col].insert(std::make_pair(this, k));
			continue;
		}
		if (dv.matrix == sNum) {
			dvInfluenceMean[k] = dvInfluenceMean[k] | (hasMean[ dv.col ] != 0.0);
			dvInfluenceVar[k] = true;
			dvContribution[dv.col].insert(std::make_pair(this, k));
			continue;
		}
		if (tracked.find(dv.matrix) == tracked.end()) {
			auto mat = currentState->matrixList[dv.matrix];
			if (verbose >= 1) mxLog("%s: %s at %s[%d,%d] tracking not implemented",
						name, omxDataColumnName(data, dv.column),
						mat->name(), 1+dv.row, 1+dv.col);
			dvInfluenceMean[k] = true;
			dvInfluenceVar[k] = true;
		}
	}
}

void omxRAMExpectation::logDefVarsInfluence()
{
	int numDefVars = data->defVars.size();
	for (int k=0; k < numDefVars; ++k) {
		omxDefinitionVar &dv = data->defVars[k];
		auto mat = currentState->matrixList[dv.matrix];
		mxLog("%s: %s->%s[%d,%d] affects mean=%d var=%d", name,
		      omxDataColumnName(data, dv.column),
		      mat->name(), 1+dv.row, 1+dv.col,
		      int(dvInfluenceMean[k]), int(dvInfluenceVar[k]));
	}
}

namespace RelationalRAMExpectation {

	omxExpectation *addr::getModel(FitContext *fc)
	{
		return omxExpectationFromIndex(model->expNum, fc->state);
	}

	omxRAMExpectation *addr::getRAMExpectation(FitContext *fc)
	{
		return (omxRAMExpectation*) getModel(fc);
	}

	std::vector< omxMatrix* > &addr::getBetween() const
	{
		return getRAMExpectationReadOnly()->between;
	}

	void omxDataRow(omxExpectation *model, int frow, omxMatrix *smallCol)
	{
		omxDataRow(model->data, frow, model->getDataColumns(), smallCol);
	}

	int addr::numVars() const
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) model;
		return ram->F->cols;
	}

	// verify whether sparse can deal with parameters set to exactly zero TODO

	independentGroup &independentGroup::getParent()
	{
		return *st.getParent().group[arrayIndex];
	}

	int state::getOptimizeMean()
	{
		return ((omxRAMExpectation*)homeEx)->optimizeMean;
	}

	// Similar to connectedness of an undirected graph
	void state::computeConnected(std::vector<int> &region, SubgraphType &connected)
	{
		Connectedness cc(region, connected, layout.size(), verbose() >= 3);

		for (int ax=int(layout.size())-1; ax >= 0; --ax) {
			cc.log();
			addr &a1 = layout[ax];
			std::vector< omxMatrix* > &between = a1.getBetween();
			if (a1.rampartScale == 0.0 || !between.size()) continue;
			for (size_t jx=0; jx < between.size(); ++jx) {
				omxMatrix *b1 = between[jx];
				int key = omxKeyDataElement(a1.getData(), a1.row, b1->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxExpectation *e1 = b1->getJoinModel();
				int row = e1->data->lookupRowOfKey(key);
				RowToLayoutMapType::const_iterator it =
					rowToLayoutMap.find(std::make_pair(e1->data, row));
				if (it == rowToLayoutMap.end())
					mxThrow("Cannot find row %d in %s", row, e1->data->name);
				int bx = it->second;
				cc.connect(ax, bx);
			}
		}
	}

	void independentGroup::MpcIO::recompute(FitContext *fc)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			if (ram->M) omxRecompute(ram->M, fc);
		}
	}

	unsigned independentGroup::MpcIO::getVersion(FitContext *fc)
	{
		unsigned v = 0;
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			if (ram->M) v += omxGetMatrixVersion(ram->M);
		}
		return v;
	}

	void independentGroup::MpcIO::refresh(FitContext *fc)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			placement &pl = par.placements[ax];
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			expectation->loadDefVars(a1.row);
			if (ram->M) {
				omxRecompute(ram->M, fc);
				EigenVectorAdaptor eM(ram->M);
				full.col(0).segment(pl.modelStart, eM.size()) = eM;
			}
		}
	}

	void independentGroup::ApcIO::recompute(FitContext *fc)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			omxRecompute(ram->A, fc);

			const double scale = a1.rampartScale;
			if (scale == 0.0) continue;

			omxData *data = expectation->data;
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxRecompute(betA, fc);
			}
		}
	}

	unsigned independentGroup::ApcIO::getVersion(FitContext *fc)
	{
		unsigned v = useRampart * 100000;
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			v += omxGetMatrixVersion(ram->A);

			const double scale = a1.rampartScale;
			if (scale == 0.0) continue;

			omxData *data = expectation->data;
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				v += omxGetMatrixVersion(betA);
			}
		}
		return v;
	}

	template <typename T>
	void independentGroup::ApcIO::_refresh(FitContext *fc, T &mat, double sign)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			placement &pl = par.placements[ax];
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			expectation->loadDefVars(a1.row);
			omxRecompute(ram->A, fc);

			if (sign == 1) {
				double *dat = ram->A->data;
				for (auto &v : *ram->Acoeff) mat.coeffRef(pl.modelStart + v.c, pl.modelStart + v.r) = dat[v.off];
			} else {
				double *dat = ram->A->data;
				for (auto &v : *ram->Acoeff) mat.coeffRef(pl.modelStart + v.c, pl.modelStart + v.r) = -dat[v.off];
			}

			double scale = 1.0;
			if (useRampart) {
				scale = a1.rampartScale;
				if (scale == 0.0) continue;
			}

			omxData *data = expectation->data;
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxData *data1 = betA->getJoinModel()->data;
				int frow = data1->lookupRowOfKey(key);
				RowToPlacementMapType::iterator plIndex =
					par.rowToPlacementMap.find(std::make_pair(data1, frow));
				placement &p2 = par.placements[ plIndex->second ];
				omxRecompute(betA, fc);
				EigenMatrixAdaptor eBetA(betA);
				omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel();
				for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
					for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
						double val = eBetA(rx,cx);
						if (val == 0.0) continue;
						mat.coeffRef(p2.modelStart + cx, pl.modelStart + rx) = sign * val * scale;
					}
				}
			}
		}
	}

	void independentGroup::SpcIO::recompute(FitContext *fc)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			omxRecompute(ram->S, fc);
		}
	}

	unsigned independentGroup::SpcIO::getVersion(FitContext *fc)
	{
		unsigned v = 0;
		for (int ax=0; ax < clumpSize; ++ax) {
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			v += omxGetMatrixVersion(ram->S);
		}
		return v;
	}

	template <typename T>
	void independentGroup::SpcIO::_refresh(FitContext *fc, T &mat)
	{
		for (int ax=0; ax < clumpSize; ++ax) {
			placement &pl = par.placements[ax];
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			expectation->loadDefVars(a1.row);
			omxRecompute(ram->S, fc);
			double *dat = ram->S->data;
			for (auto &v : *ram->Scoeff) {
				mat.coeffRef(pl.modelStart + v.r, pl.modelStart + v.c) = dat[v.off];
			}
		}
	}

	// 1st visitor
	int state::flattenOneRow(omxExpectation *expectation, int frow, int &maxSize)
	{
		allEx.insert(expectation);
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

		if (data->hasPrimaryKey()) {
			// insert_or_assign would be nice here
			RowToLayoutMapType::const_iterator it = rowToLayoutMap.find(std::make_pair(data, frow));
			if (it != rowToLayoutMap.end()) {
				if (it->second < 0) mxThrow("%s cycle detected: '%s' row %d joins against itself",
							     homeEx->name, data->name, 1+frow);
				return it->second;
			}

			rowToLayoutMap[ std::make_pair(data, frow) ] = -1;
		}

		struct addr a1;
		struct addrSetup as1;
		as1.clumped = false;
		a1.ig = 0;
		as1.parent1 = NA_INTEGER;
		as1.fk1 = NA_INTEGER;
		as1.numJoins = 0;
		as1.numKids = 0;
		as1.heterogenousMean = false;
		as1.rset = NA_INTEGER;
		a1.rampartScale = 1.0;
		a1.quickRotationFactor = 1.0;
		as1.skipMean = NA_INTEGER;
		a1.row = frow;
		a1.nextMean = 1;
		a1.setModel(expectation);
		//as1.key = frow;

		std::vector<int> parents;
		parents.reserve(ram->between.size());

		for (size_t jx=0; jx < ram->between.size(); ++jx) {
			omxMatrix *b1 = ram->between[jx];
			int key = omxKeyDataElement(data, frow, b1->getJoinKey());
			if (key == NA_INTEGER) continue;
			omxExpectation *e1 = b1->getJoinModel();
			int parentPos = flattenOneRow(e1, e1->data->lookupRowOfKey(key), maxSize);
			if (jx == 0) {
				as1.fk1 = key;
				as1.parent1 = parentPos;
			}
			parents.push_back(parentPos);
		}

		for (size_t jx=0; jx < parents.size(); ++jx) {
			addrSetup &pop = layoutSetup[ parents[jx] ];
			pop.numKids += 1;
			as1.numJoins += 1;
		}

		a1.numObsCache = 0;
		int jCols = expectation->getDataColumns().size();
		if (jCols) {
			if (!ram->M) {
				complainAboutMissingMeans(expectation);
				return 0;
			}
			if (smallCol->cols < jCols) {
				omxResizeMatrix(smallCol, 1, jCols);
			}
			omxDataRow(expectation, frow, smallCol);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(smallCol, 0, col);
				bool yes = std::isfinite(val);
				if (yes) ++a1.numObsCache;
			}
		}

		layout.push_back(a1);
		layoutSetup.push_back(as1);

		if (data->hasPrimaryKey()) {
			rowToLayoutMap[ std::make_pair(data, frow) ] = layout.size() - 1;
			//a1.key = data->primaryKeyOfRow(frow);
		}

		maxSize += ram->F->cols;
		return layout.size()-1;
	}

	const std::vector<bool> &addr::getDefVarInfluenceMean() const
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) model;
		return ram->dvInfluenceMean;
	}

	const std::vector<bool> &addr::getDefVarInfluenceVar() const
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) model;
		return ram->dvInfluenceVar;
	}

	void addr::dataRow(omxMatrix *out) const
	{
		omxDataRow(model, row, out);
	}

	struct CompareLib {
		state &st;
		CompareLib(state *_st) : st(_st->getParent()) {};

		// actually stores !missingness
		template <typename T>
		void getMissingnessPattern(const addr &a1, std::vector<T> &out) const
		{
			a1.dataRow(st.smallCol);
			int jCols = a1.getDataColumns().size();
			out.reserve(jCols);
			for (int col=0; col < jCols; ++col) {
				double val = omxMatrixElement(st.smallCol, 0, col);
				out.push_back(std::isfinite(val));
			}
		}

		bool compareMissingnessAndCov(const addr &la, const addr &ra, bool &mismatch) const
		{
			mismatch = la.getExpNum() != ra.getExpNum();
			if (mismatch) return la.getExpNum() < ra.getExpNum();

			mismatch = la.numVars() != ra.numVars();
			if (mismatch) return la.numVars() < ra.numVars();

			std::vector<bool> lmp;
			getMissingnessPattern(la, lmp);
			std::vector<bool> rmp;
			getMissingnessPattern(ra, rmp);

			mismatch = lmp.size() != rmp.size();
			if (mismatch) return lmp.size() < rmp.size();

			for (size_t lx=0; lx < lmp.size(); ++lx) {
				mismatch = lmp[lx] != rmp[lx];
				if (!mismatch) continue;
				return int(lmp[lx]) < int(rmp[lx]);
			}

			mismatch = la.rampartScale != ra.rampartScale;
			if (mismatch) return la.rampartScale < ra.rampartScale;

			bool got = compareCovDefVars(la, ra, mismatch);
			if (mismatch) return got;

			return false;
		}

		bool cmpCovClump(const addr &la, const addr &ra, bool &mismatch) const
		{
			bool got;

			got = compareMissingnessAndCov(la, ra, mismatch);
			if (mismatch) return got;

			const addrSetup &lhss = st.layoutSetup[&la - &st.layout[0]];
			const addrSetup &rhss = st.layoutSetup[&ra - &st.layout[0]];
			mismatch = lhss.clump.size() != rhss.clump.size();
			if (mismatch) return lhss.clump.size() < rhss.clump.size();

			for (size_t cx=0; cx < lhss.clump.size(); ++cx) {
				got = cmpCovClump(st.layout[lhss.clump[cx]],
						   st.layout[rhss.clump[cx]], mismatch);
				if (mismatch) return got;
			}

			return false;
		}

		bool compareCovDefVars(const addr &la, const addr &ra, bool &mismatch) const
		{
			mismatch = true;

			auto &dvInfluenceVar = la.getDefVarInfluenceVar();
			omxData *data = la.getData();  // both la & ra have same data
			for (size_t k=0; k < data->defVars.size(); ++k) {
				if (!dvInfluenceVar[k]) continue;
				int col = data->defVars[k].column;
				double lv = omxDoubleDataElement(data, la.row, col);
				double rv = omxDoubleDataElement(data, ra.row, col);
				if (lv == rv) continue;
				return lv < rv;
			}

			mismatch = false;
			return false;
		}

		bool compareMeanDefVars(const addr &la, const addr &ra, bool &mismatch) const
		{
			mismatch = true;

			auto &dvInfluenceMean = la.getDefVarInfluenceMean();
			omxData *data = la.getData();  // both la & ra have same data
			for (size_t k=0; k < data->defVars.size(); ++k) {
				if (!dvInfluenceMean[k]) continue;
				int col = data->defVars[k].column;
				double lv = omxDoubleDataElement(data, la.row, col);
				double rv = omxDoubleDataElement(data, ra.row, col);
				if (lv == rv) continue;
				return lv < rv;
			}

			mismatch = false;
			return false;
		}

	};

	struct CompatibleCovCompare : CompareLib {
		CompatibleCovCompare(state *_st) : CompareLib(_st) {};

		bool operator() (const std::vector<int> &lhs, const std::vector<int> &rhs) const
		{
			if (lhs.size() != rhs.size()) return lhs.size() < rhs.size();
			bool mismatch;
			for (size_t ux=0; ux < lhs.size(); ++ux) {
				addr &la = st.layout[lhs[ux]];
				addr &ra = st.layout[rhs[ux]];
				bool got = compareMissingnessAndCov(la, ra, mismatch);
				if (mismatch) return got;
			}
			return false;
		}
	};

	struct CompatibleMeanCompare : CompareLib {
		CompatibleMeanCompare(state *_st) : CompareLib(_st) {};

		addr *joinedWith(const addr &la, int jx) const
		{
			omxRAMExpectation *ram = la.getRAMExpectationReadOnly();
			omxData *data = la.getData();
			omxMatrix *betA = ram->between[jx];
			int key = omxKeyDataElement(data, la.row, betA->getJoinKey());
			if (key == NA_INTEGER) return 0;
			omxData *data1 = betA->getJoinModel()->data;
			omxExpectation *e1 = betA->getJoinModel();
			int row = data1->lookupRowOfKey(key);
			state::RowToLayoutMapType::const_iterator it =
				st.rowToLayoutMap.find(std::make_pair(e1->data, row));
			if (it == st.rowToLayoutMap.end())
				mxThrow("Cannot find row %d in %s", row, e1->data->name);
			return &st.layout[it->second];
		}

		bool compareMeanDeep(const addr &la, const addr &ra, bool &mismatch) const
		{
			omxRAMExpectation *ram = la.getRAMExpectationReadOnly();
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				addr *lp = joinedWith(la, jx);
				addr *rp = joinedWith(ra, jx);
				if (!lp && !rp) continue;
				if (!lp || !rp) {
					mismatch = true;
					return !lp;
				}
				bool got = compareMeanDeep(*lp, *rp, mismatch);
				if (mismatch) return got;
			}

			bool got = compareMeanDefVars(la, ra, mismatch);
			if (mismatch) return got;

			return false;
		}

		bool operator() (const std::vector<int> &lhs, const std::vector<int> &rhs) const
		{
			for (size_t ux=0; ux < lhs.size(); ++ux) {
				addr &la = st.layout[lhs[ux]];
				addr &ra = st.layout[rhs[ux]];
				bool mismatch;
				bool got = compareMeanDeep(la, ra, mismatch);
				if (mismatch) return got;
			}
			return false;
		}
	};

	struct RampartTodoCompare : CompareLib {
		RampartTodoCompare(state *_st) : CompareLib(_st) {};

		bool operator() (const addr *lhs, const addr *rhs) const
		{
			const addrSetup *lhss = &st.layoutSetup[lhs - &st.layout[0]];
			const addrSetup *rhss = &st.layoutSetup[rhs - &st.layout[0]];

			bool mismatch = lhss->fk1 != rhss->fk1;
			if (mismatch) return lhss->fk1 < rhss->fk1;

			return cmpCovClump(*lhs, *rhs, mismatch);
		}
	};

	struct RampartClumpCompare : CompareLib {
		RampartClumpCompare(state *_st) : CompareLib(_st) {};

		bool clumpCmp(const int lhs, const int rhs) const {
			bool mismatch;
			bool got = cmpCovClump(st.layout[lhs], st.layout[rhs], mismatch);
			if (mismatch) return got;
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
	void state::appendClump(int ax, std::vector<T> &clump)
	{
		clump.push_back(ax);
		addrSetup &a1 = layoutSetup[ax];
		for (size_t cx = 0; cx < a1.clump.size(); ++cx) {
			appendClump(a1.clump[cx], clump);
		}
	}

	template <typename T>
	void state::propagateDefVar(omxRAMExpectation *to, Eigen::MatrixBase<T> &transition,
				    omxRAMExpectation *from)
	{
		bool within = to == from;
		to->hasMean     += (transition * from->hasMean    ).array().abs().matrix();
		to->hasVariance += (transition * from->hasVariance).array().abs().matrix();

		for (int rx=0; rx < transition.rows(); ++rx) {
			for (int cx=0; cx < transition.cols(); ++cx) {
				if (within && rx == cx) continue;
				if (transition(rx,cx) == 0) continue;
				auto &fromDv = from->dvContribution[cx];

				bool hasMean = false;
				bool hasVar = false;
				if (!within) {
					hasMean = from->hasMean[cx] != 0.;
					hasVar  = from->hasVariance[cx] != 0.;
				}
				if (verbose() >= 1) {
					for (auto it : fromDv) {
						omxRAMExpectation *dvHome =
							(omxRAMExpectation*) it.first;
						auto &dv = dvHome->data->defVars[ it.second ];
						auto mat = dvHome->currentState->matrixList[dv.matrix];
						mxLog("%s at %s[%d,%d] goes from %s to %s w/ mean %d var %d (0=no effect)",
						      omxDataColumnName(dvHome->data, dv.column),
						      mat->name(), 1+dv.row, 1+dv.col,
						      from->S->rownames[cx], to->S->rownames[rx],
						      hasMean, hasVar);
					}
				}

				if (hasMean || hasVar) {
					for (auto it : fromDv) {
						omxRAMExpectation *dvHome =
							(omxRAMExpectation*) it.first;
						if (hasMean)
							dvHome->dvInfluenceMean[ it.second ] = true;
						if (hasVar)
							dvHome->dvInfluenceVar[ it.second ] = true;
					}
				}
				auto &dv1 = to->dvContribution[rx];
				dv1.insert(fromDv.begin(), fromDv.end());
			}
		}
	}

	void state::analyzeModel1(FitContext *fc)
	{
		std::set<omxExpectation*> checkedEx;

		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			omxExpectation *expectation = a1.getModel(fc);

			if (checkedEx.find(expectation) != checkedEx.end()) continue;

			omxData *data = expectation->data;
			// def vars are addressed in analyzeModel2.
			// Not sure if we need to consider them here, but it is
			// the conservative approach.
			data->loadFakeData(expectation->currentState, 1.0);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxExpectation *ex2 = betA->getJoinModel();
				omxRAMExpectation *ram2 = (omxRAMExpectation*) ex2;
				omxRecompute(betA, fc);
				betA->markPopulatedEntries();
				EigenMatrixAdaptor eBA(betA);
				propagateDefVar(ram, eBA, ram2);
			}

			omxRecompute(ram->A, fc);
			ram->A->markPopulatedEntries();

			// Can ignore latent node type. Compute the inverse since it's a one-time cost.
			EigenMatrixAdaptor eA(ram->A);
			Eigen::MatrixXd fullA(eA);
			fullA.diagonal().array() = 1;
			Eigen::FullPivLU< Eigen::MatrixXd > lu(fullA);
			Eigen::MatrixXd IA(fullA.rows(), fullA.cols());
			IA.setIdentity();
			IA = lu.solve(IA);

			propagateDefVar(ram, IA, ram);

			checkedEx.insert(expectation);
			if (checkedEx.size() == allEx.size()) break;
		}
	}

	void state::analyzeModel2(FitContext *fc)
	{
		for (auto it : allEx) {
			omxRAMExpectation *ram = (omxRAMExpectation*) it;
			if (ram->getThresholdInfo().size()) {
				mxThrow("%s: Ordinal indicators are not supported in multilevel models", ram->name);
			}
			omxData *data = ram->data;
			int numDefVars = data->defVars.size();
			if (numDefVars == 0) continue;

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				auto *from = (omxRAMExpectation*) betA->getJoinModel();
				int bNum = ~betA->matrixNumber;

				for (int k=0; k < numDefVars; ++k) {
					omxDefinitionVar &dv = data->defVars[k];
					if (dv.matrix == bNum) {
						ram->dvInfluenceMean[k] = from->hasMean[dv.col] != 0.;
						ram->dvInfluenceVar[k] = from->hasVariance[dv.col] != 0.;
					}
				}
			}

			int aNum = ~ram->A->matrixNumber;
			for (int k=0; k < numDefVars; ++k) {
				omxDefinitionVar &dv = data->defVars[k];
				if (dv.matrix == aNum) {
					ram->dvInfluenceMean[k] = ram->hasMean[dv.col] != 0.;
					ram->dvInfluenceVar[k] = ram->hasVariance[dv.col] != 0.;
				}
			}
		}
	}

	template <typename T>
	bool state::placeSet(std::set<std::vector<T> > &toPlace, independentGroup *ig)
	{
		bool heterogenousMean = false;
		for (std::set<std::vector<int> >::iterator px = toPlace.begin();
		     px != toPlace.end(); ++px) {
			const std::vector<int> &clump = *px;
			for (size_t cx=0; cx < clump.size(); ++cx) {
				heterogenousMean |= layoutSetup[ clump[cx] ].heterogenousMean;
				ig->place(clump[cx]);
			}
		}
		return heterogenousMean;
	}

	// 2nd visitor
	void state::planModelEval(int maxSize, FitContext *fc)
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) homeEx;
		if (ram->forceSingleGroup) {
			independentGroup *ig = new independentGroup(this, layout.size(), layout.size());
			for (size_t ax=0; ax < layout.size(); ++ax) ig->place(ax);
			ig->prep(fc);
			group.push_back(ig);
			return;
		}

		if (verbose() >= 2) {
			mxLog("%s: analyzing unit dependencies", homeEx->name);
		}

		// Here we sort clumps into independent groups of identical covariance

		std::vector<int> region;
		SubgraphType connected;
		computeConnected(region, connected);

		// connected gives the complete dependency information,
		// but we already have partial dependency information
		// from Rampart clumping. We need to preserve the
		// Rampart clumping order when we determine the
		// grouping. Otherwise we can get more groups (and
		// fewer copies) than ideal.

		typedef std::map< std::vector<int>,
				  std::set<std::vector<int> >,
				  CompatibleCovCompare> CompatibleCovMapType;
		CompatibleCovMapType cgm(this);
		for (size_t ax=0; ax < layout.size(); ++ax) {
			if (region[ax] == -1) {
				std::vector<int> clump;
				clump.push_back(ax);
				cgm[ clump ].insert(clump);
				continue;
			}
			std::set<int> &unsortedClump = connected[ region[ax] ];
			if (!unsortedClump.size()) continue;  //already done
			std::vector<int> clump;
			clump.reserve(unsortedClump.size());
			while (unsortedClump.size()) {
				bool movedSome = false;
				for (std::set<int>::iterator it=unsortedClump.begin(); it!=unsortedClump.end(); ++it) {
					addrSetup &as2 = layoutSetup[*it];
					if (as2.clumped) continue;
					int beforeSize = clump.size();
					appendClump(*it, clump);
					for (size_t cx=beforeSize; cx < clump.size(); ++cx) {
						unsortedClump.erase(clump[cx]);
					}
					movedSome = true;
					break;
				}
				if (!movedSome) break;
			}
			// Not sure if order matters here TODO
			clump.insert(clump.end(), unsortedClump.begin(), unsortedClump.end());
			cgm[ clump ].insert(clump);
		}

		// if lots of copies==1 then we need a different strategy TODO

		if (verbose() >= 1) {
			mxLog("%s: will create %d independent groups", homeEx->name, int(cgm.size()));
		}
		group.reserve(cgm.size());

		for (CompatibleCovMapType::iterator it = cgm.begin();
		     it != cgm.end(); ++it) {
			independentGroup *ig = new independentGroup(this, it->second.size(),
								    it->second.begin()->size());
			typedef std::map< std::vector<int>,
					  std::set<std::vector<int> >,
					  CompatibleMeanCompare> CompatibleMeanMapType;
			CompatibleMeanMapType cmm(this);
			for (std::set<std::vector<int> >::iterator px = it->second.begin();
			     px != it->second.end(); ++px) {
				const std::vector<int> &clump = *px;
				cmm[ clump ].insert(clump);
			}

			int ssCount = 0;
			for (CompatibleMeanMapType::iterator mit = cmm.begin();
			     mit != cmm.end(); ++mit) {
				if (mit->second.size() > 1) {
					++ssCount;
					continue;
				}
				placeSet(mit->second, ig);
			}
			ig->sufficientSets.reserve(ssCount);
			for (CompatibleMeanMapType::iterator mit = cmm.begin();
			     mit != cmm.end(); ++mit) {
				if (mit->second.size() == 1) continue;
				int from = ig->placements.size();
				if (placeSet(mit->second, ig)) continue;
				if (verbose() >= 3) {
					mxLog("group %d same mean %d -> %d clumpsize %d",
					      int(group.size()), from, int(ig->placements.size() - 1),
					      int(it->second.begin()->size()));
				}
				sufficientSet ss;
				ss.start = from / ig->clumpSize;
				ss.length = (ig->placements.size() - from) / ig->clumpSize;
				ig->sufficientSets.push_back(ss);
			}
			if (!ram->useSufficientSets) ig->sufficientSets.clear();
			ig->prep(fc);
			group.push_back(ig);
		}
	}

	void independentGroup::place(int ax)
	{
		if (st.layout[ax].ig) {
			mxThrow("Unit[%d] already assigned; this is a bug", ax);
		}
		st.layout[ax].ig = this;
		int mx = 0;
		int dx = 0;
		if (placements.size()) {
			int last = placements.size()-1;
			placement &prev = placements[last];
			addr &a1 = st.layout[ gMap[last] ];
			mx = prev.modelStart + a1.numVars();
			dx = prev.obsStart + a1.numObs();
		}
		placement pl;
		pl.modelStart = mx;
		pl.obsStart = dx;
		placements.push_back(pl);
		gMap.push_back(ax);
	}

	independentGroup::independentGroup(class state *_st, int size, int _clumpSize)
		: st(*_st), clumpSize(_clumpSize)
	{
		placements.reserve(size);
	}

	independentGroup::independentGroup(independentGroup *ig)
		: st(ig->st), clumpSize(ig->clumpSize)
	{
		arrayIndex = ig->arrayIndex;
		obsNameVec = 0;
		varNameVec = 0;
		skipMean = ig->skipMean;
		expectedVec.resize(ig->expectedVec.size());
		if (skipMean) expectedVec.setZero();
		fullMean.resize(ig->fullMean.size());
		fullMean.setZero();
		clumpVars = ig->clumpVars;
		clumpObs = ig->clumpObs;
		pcalc.clone(ig->pcalc);
	}

	void independentGroup::prep(FitContext *fc)
	{
		skipMean = NA_INTEGER;
		int totalObserved = 0;
		int maxSize = 0;
		if (placements.size()) {
			int last = placements.size()-1;
			placement &prev = placements[last];
			addr &a1 = st.layout[ gMap[last] ];
			totalObserved = prev.obsStart + a1.numObs();
			maxSize = prev.modelStart + a1.numVars();
		}
		if (verbose() >= 2) {
			mxLog("%s: create independentGroup[%d] maxSize=%d totalObserved=%d",
			      st.homeEx->name, (int)st.group.size(), maxSize, totalObserved);
		}
		latentFilter.assign(maxSize, false); // will have totalObserved true entries
		isProductNode.assign(maxSize, false);
		varNameVec = Rf_protect(Rf_allocVector(STRSXP, maxSize));
		dataVec.resize(totalObserved);
		dataColumn.resize(totalObserved);
		dataColumn.setConstant(-1);
		fullMean.resize(maxSize);
		fullMean.setZero();
		if (0) {
			rawFullMean.resize(maxSize);
			rawFullMean.setZero();
		}

		{
			int last = clumpSize-1;
			placement &end = placements[last];
			addr &a1 = st.layout[ gMap[last] ];
			clumpVars = end.modelStart + a1.numVars();
			clumpObs = end.obsStart + a1.numObs();
		}
		int tmpObs = st.hasProductNodes? clumpObs : totalObserved;
		expectedVec.resize(tmpObs);
		obsNameVec = Rf_protect(Rf_allocVector(STRSXP, tmpObs));

		int dx=0;
		for (size_t ax=0; ax < placements.size(); ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[ gMap[ax] ];
			a1.igIndex = ax;

			if (verbose() >= 3) {
				// useless diagnostic?
				int modelEnd = pl.modelStart + a1.numVars() - 1;
				if (a1.numObs()) {
					mxLog("place %s[%d] at %d %d obs %d %d", a1.modelName().c_str(),
					      a1.row, pl.modelStart, modelEnd, pl.obsStart, pl.obsStart + a1.numObs() - 1);
				} else {
					mxLog("place latent %s[%d] at %d %d", a1.modelName().c_str(),
					      a1.row, pl.modelStart, modelEnd);
				}
			}

			omxData *data = a1.getData();
			rowToPlacementMap[ std::make_pair(data, a1.row) ] = ax;

			omxRAMExpectation *ram = a1.getRAMExpectation(fc);

			std::string modelName(data->name);
			modelName = modelName.substr(0, modelName.size() - 4); // remove "data" suffix

			auto dc = a1.getDataColumns();
			if (dc.size()) {
				int prevDx = dx;
				a1.dataRow(st.smallCol);
				for (int vx=0, ncol=0; vx < ram->F->cols; ++vx) {
					isProductNode[ pl.modelStart + vx ] = ram->isProductNode[vx];
					if (!ram->latentFilter[vx]) continue;
					int col = ncol++;
					double val = omxMatrixElement(st.smallCol, 0, col);
					bool yes = std::isfinite(val);
					if (!yes) continue;
					latentFilter[ pl.modelStart + vx ] = true;
					std::string dname =
						modelName + omxDataColumnName(data, dc[col]);
					if (dx < tmpObs) SET_STRING_ELT(obsNameVec, dx, Rf_mkChar(dname.c_str()));
					dataVec[ dx ] = val;
					if (a1.getExpNum() == st.homeEx->expNum) dataColumn[ dx ] = col;
					dx += 1;
				}
				if (a1.numObs() != dx - prevDx) {
					mxThrow("numObs() %d != %d", a1.numObs(), dx - prevDx);
				}
			}
			for (int vx=0; vx < ram->F->cols; ++vx) {
				std::string dname = modelName + ram->F->colnames[vx];
				SET_STRING_ELT(varNameVec, pl.modelStart + vx, Rf_mkChar(dname.c_str()));
			}
		}

		pcalc.attach(clumpVars, clumpObs, latentFilter, isProductNode, new MpcIO(*this),
								 new ApcIO(*this), new SpcIO(*this));
		pcalc.setAlgo(fc, st.hasProductNodes, NA_INTEGER);
	}

	template <typename T>
	void state::oertzenRotate(std::vector<T> &t1, bool canOptimize)
	{
		addrSetup &specimen = layoutSetup[ t1[0] ];
		CompatibleMeanCompare cmp(this);

		bool mismatch = false;
		for (int cx=1; cx < int(t1.size()); ++cx) {
			cmp.compareMeanDeep(layout[ t1[0] ], layout[ t1[cx] ], mismatch);
			if (mismatch) break;
		}
		bool keep = true;
		if (mismatch) {
			for (int cx=0; cx < int(t1.size()); ++cx) {
				layoutSetup[ t1[cx] ].heterogenousMean = true;
			}
		} else if (canOptimize) {
			layout[ t1[0] ].quickRotationFactor *= sqrt(double(t1.size()));
			for (int vx=1; vx < int(t1.size()); ++vx) {
				layoutSetup[ t1[vx] ].skipMean = 1;
				layout[ t1[vx] ].quickRotationFactor = 0.;
			}
			keep = false;
		}

		modelRotationPlanFilter.push_back(keep);
		rotationPlan.push_back(t1);

		for (size_t cx=0; cx < specimen.clump.size(); ++cx) {
			std::vector<int> t2;
			t2.reserve(t1.size());
			for (size_t tx=0; tx < t1.size(); ++tx) {
				addrSetup &a1 = layoutSetup[ t1[tx] ];
				if (specimen.clump.size() != a1.clump.size()) {
					mxThrow("BUG: clump size mismatch; alert developers");
				}
				t2.push_back(a1.clump[cx]);
			}
			oertzenRotate(t2, canOptimize);
		}
	}

	int state::rampartRotate(int level)
	{
		typedef std::map< addr*, std::vector<int>, RampartTodoCompare > RampartMap;
		RampartMap todo(RampartTodoCompare(this));
		int unlinked = 0;

		int loopTo = layout.size();
		int rampartUnitLimit = ((omxRAMExpectation*) homeEx)->rampartUnitLimit;
		if (rampartUnitLimit != NA_INTEGER) loopTo = std::min(rampartUnitLimit, loopTo);
		for (int ax=0; ax < loopTo; ++ax) {
			addr &a1 = layout[ax];
			addrSetup &as1 = layoutSetup[ax];
			if (as1.numKids != 0 || as1.numJoins != 1 || as1.clumped) continue;
			// numKids == 0 && numJoints == 1 && !as1.clumped
			std::vector<int> &t1 = todo[&a1];
			t1.push_back(int(ax));
		}
		for (RampartMap::iterator it = todo.begin(); it != todo.end(); ++it) {
			std::vector<int> &t1 = it->second;

			if (t1.size() >= 2) {
				//std::string buf = "rotate units";

				oertzenRotate(t1, getOptimizeMean() >= 1);
				addr &specimen = layout[ t1[0] ];
				specimen.rampartScale = sqrt(double(t1.size()));
				int parent1 = layoutSetup[ t1[0] ].parent1;
				layoutSetup[parent1].numKids -= t1.size();
				clumpWith(parent1, t1[0]);
				for (size_t ux=1; ux < t1.size(); ++ux) {
					layout[ t1[ux] ].rampartScale = 0;
					layoutSetup[ t1[ux] ].numJoins = 0;
					//buf += string_snprintf(" %d", 1+t1[ux]);
				}
				//buf += string_snprintf(" -> %d\n", 1+t1[0]);
				//mxLogBig(buf);
			} else {
				// Don't rotate, just clump units together with parent.
				int parent1 = layoutSetup[ t1[0] ].parent1;
				layoutSetup[parent1].numKids -= t1.size();
				for (size_t ux=0; ux < t1.size(); ++ux) {
					clumpWith(parent1, t1[ux]);
				}
			}
			// not really unlinked in clumped case, but layout is changed; fix reporting TODO
			unlinked += t1.size() - 1;
		}
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addrSetup &a1 = layoutSetup[ax];
			if (a1.clump.size() <= 1) continue;
			std::sort(a1.clump.begin(), a1.clump.end(), RampartClumpCompare(this));
			if (false) {
				// if Compare is screwed up then it should show up here
				std::vector<int> dc(a1.clump);
				std::sort(dc.begin(), dc.end(), RampartClumpCompare(this));
				for (size_t cx=0; cx < dc.size(); ++cx) {
					if (dc[cx] != a1.clump[cx]) mxThrow("oops");
				}
			}
		}
		return unlinked;
	}

	template <bool model>
	struct UnitAccessor {
		state &st;
		UnitAccessor(state *_st) : st(*_st) {};
		bool isModel() const { return model; };

		// split into coeff & coeffRef versions TODO
		double &operator() (const int unit, const int obs)
		{
			addr &ad = st.getParent().layout[unit];
			independentGroup &ig = *ad.ig;
			independentGroup &tig = *st.group[ad.ig->arrayIndex];
			int obsStart = ig.placements[ad.igIndex].obsStart;
			return (model? tig.expectedVec : ig.dataVec).coeffRef(obsStart + obs);
		};
	};

	template <typename T>
	void state::unapplyRotationPlan(T accessor)
	{
		for (size_t rx=0; rx < rotationPlan.size(); ++rx) {
			const std::vector<int> &units = rotationPlan[rx];
			int numUnits = units.size();
			const addr &specimen = layout[units[0]];
			for (int ox=0; ox < specimen.numObs(); ++ox) {
				double p1 = sqrt(1.0/numUnits) * accessor(units[0], ox);
				for (int ii = 0; ii < numUnits; ii++) {
					double k = numUnits-ii;
					if (ii >= 1 && ii < numUnits-1) {
						p1 += sqrt(1.0/(k*(k+1))) * accessor(units[ii], ox);
					}
					double p2;
					if (ii >= numUnits-2) {
						p2 = M_SQRT1_2;
						if (ii == numUnits - 1) p2 = -p2;
					} else {
						p2 = -sqrt((k-1.0)/k);
					}
					accessor(units[ii], ox) =
						p1 + p2 * accessor(units[std::min(ii+1, numUnits-1)], ox);
				}
			}
		}
	}

	template <typename T>
	void state::applyRotationPlan(T accessor)
	{
		// maybe faster to do all observations in parallel
		// to allow more possibility of instruction reordering TODO
		const bool debug = false;
		std::string buf;
		for (size_t rx=0; rx < rotationPlan.size(); ++rx) {
			if (debug) {
				buf += string_snprintf("rotate<model=%d> step[%d]",
						       accessor.isModel(), int(rx));
			}
			const std::vector<int> &units = rotationPlan[rx];

			const addr &specimen = layout[units[0]];
			for (int ox=0; ox < specimen.numObs(); ++ox) {
				if (debug) buf += string_snprintf(" obs[%d]", ox);
				double partialSum = 0.0;
				for (size_t ux=0; ux < units.size(); ++ux) {
					partialSum += accessor(units[ux], ox);
					if (debug) buf += string_snprintf(" %d", int(1+ units[ux]));
				}

				double prev = accessor(units[0], ox);
				accessor(units[0], ox) = partialSum / sqrt(double(units.size()));
				if (debug) buf += string_snprintf(": %f", accessor(units[0], ox));

				for (size_t i=1; i < units.size(); i++) {
					double k=units.size()-i;
					partialSum -= prev;
					double prevContrib = sqrt(k / (k+1)) * prev;
					prev = accessor(units[i], ox);
					accessor(units[i], ox) =
						partialSum * sqrt(1.0 / (k*(k+1))) - prevContrib;
					if (debug) buf += string_snprintf(" %f", accessor(units[i], ox));
				}
			}
			if (debug) buf += "\n";
		}
		if (debug && buf.size()) mxLogBig(buf);
	}

	void state::optimizeModelRotation()
	{
		std::vector< std::vector<int> > origRotationPlan = rotationPlan;
		rotationPlan.clear();

		for (int px=0; px < int(origRotationPlan.size()); ++px) {
			if (modelRotationPlanFilter[px])
				rotationPlan.push_back(origRotationPlan[px]);
		}

		if (getOptimizeMean() <= 1) return;

		origRotationPlan = rotationPlan;
		rotationPlan.clear();

		// Transitive closure on reverse dependencies
		for (int px=origRotationPlan.size()-1; px >= 0; --px) {
			auto &vec = origRotationPlan[px];
			bool skip = true;
			for (int vx=0; vx < int(vec.size()); ++vx) {
				if (layoutSetup[vec[vx]].skipMean != 1) {
					skip = false;
					break;
				}
			}
			if (skip == false) {
				for (int vx=0; vx < int(vec.size()); ++vx) {
					layoutSetup[vec[vx]].skipMean = 0;
				}
			}
		}
		for (auto &vec : origRotationPlan) {
			if (layoutSetup[vec[0]].skipMean == 0) {
				rotationPlan.push_back(vec);
			}
		}
	}

	void state::init(omxExpectation *expectation, FitContext *fc)
	{
		parent = this;
		homeEx = expectation;

		omxRAMExpectation *ram = (omxRAMExpectation*) homeEx;
		int numManifest = ram->F->rows;
		smallCol = omxInitMatrix(1, numManifest, TRUE, homeEx->currentState);

		if (fc->isClone()) {
			omxExpectation *phomeEx = omxExpectationFromIndex(homeEx->expNum, fc->getParentState());
			omxRAMExpectation *pram = (omxRAMExpectation*) phomeEx;
			parent = pram->rram;
			hasProductNodes = parent->hasProductNodes;
			group.reserve(parent->group.size());
			for (size_t gx=0; gx < parent->group.size(); ++gx) {
				independentGroup *ig = parent->group[gx];
				group.push_back(new independentGroup(ig));
			}
			return;
		}

		{
			ProtectedSEXP Rdvhack(R_do_slot(expectation->rObj, Rf_install(".analyzeDefVars")));
			doAnalyzeDefVars = Rf_asLogical(Rdvhack);
		}

		int maxSize = 0;
		int homeDataRows = homeEx->data->nrows();
		for (int row=0; row < homeDataRows; ++row) {
			flattenOneRow(homeEx, row, maxSize);
			if (isErrorRaised()) return;
		}

		hasProductNodes = false;
		for (auto *ex : allEx) {
			omxRAMExpectation *ram1 = (omxRAMExpectation*) ex;
			hasProductNodes |= ram1->getHasProductNodes();
			if (!ex->data->hasWeight() && !ex->data->hasFreq()) continue;
			mxThrow("%s: row frequencies or weights provided in '%s' are not compatible with joins",
				 expectation->name, ex->data->name);
		}
		if (hasProductNodes && getOptimizeMean() < 2) {
			mxThrow("%s: .optimizeMean=2L is required for product nodes", homeEx->name);
		}

		Eigen::VectorXd paramSave;
		homeEx->currentState->setFakeParam(paramSave);

		for (auto it : allEx) {
			omxRAMExpectation *ram2 = (omxRAMExpectation*) it;
			ram2->analyzeDefVars(fc);
			if (verbose() >= 1) ram2->logDefVarsInfluence();
		}

		if (doAnalyzeDefVars) {
			analyzeModel1(fc);
			analyzeModel2(fc);

			for (auto it : allEx) {
				omxRAMExpectation *ram2 = (omxRAMExpectation*) it;
				if (verbose() >= 1) ram2->logDefVarsInfluence();
			}
		}

		for (auto it : allEx) {
			omxRAMExpectation *ram2 = (omxRAMExpectation*) it;
			ram2->dvContribution.clear();
		}

		if (ram->rampartEnabled()) {
			int maxIter = ram->rampartCycleLimit;
			int unlinked = 0;
			int level = -1; // mainly for debugging
			while (int more = rampartRotate(++level)) {
				rampartUsage.push_back(more);
				unlinked += more;
				if (maxIter != NA_INTEGER && --maxIter == 0) break;
			}
			if (verbose() >= 1) {
				mxLog("%s: rampart unlinked %d units", homeEx->name, unlinked);
			}
		}

		planModelEval(maxSize, fc);

		homeEx->currentState->restoreParam(paramSave);

		for (std::vector<independentGroup*>::iterator it = group.begin() ; it != group.end(); ++it) {
			(*it)->arrayIndex = it - group.begin();
		}

		applyRotationPlan(UnitAccessor<false>(this));
		//unapplyRotationPlan(UnitAccessor<false>(this));
		//applyRotationPlan(UnitAccessor<false>(this));

		for (std::vector<independentGroup*>::iterator it = group.begin() ; it != group.end(); ++it) {
			(*it)->finalizeData();
		}

		if (getOptimizeMean() >= 1) optimizeModelRotation();

		for (int r1=0; r1 < int(rotationPlan.size()); ++r1) {
			auto &vec = rotationPlan[r1];
			for (int r2=0; r2 < int(vec.size()); ++r2) {
				int rset = layoutSetup[vec[r2]].rset;
				if (rset == NA_INTEGER) {
					layoutSetup[vec[r2]].rset = r1;
				} else {
					// Can be subject to multiple rotations
					layoutSetup[vec[r2]].rset += 1000 + r1;
				}
			}
		}

		rotationCount = 0;
		for (auto &vec : rotationPlan) {
			rotationCount += vec.size();
		}

		// skipMean for layout[0] is always false
		for (int ax=0; ax < int(layout.size()); ax += layout[ax].nextMean) {
			int incr = 1;
			while (ax+incr < int(layout.size()) &&
			       layoutSetup[ax+incr].skipMean == 1) incr += 1;
			layout[ax].nextMean = incr;
		}
	}

	state::~state()
	{
		for (std::vector<independentGroup*>::iterator it = group.begin() ; it != group.end(); ++it) {
			delete *it;
		}

		omxFreeMatrix(smallCol);
	}

	void independentGroup::finalizeData()
	{
		if (clumpObs == 0) return;
		for (int sx=0; sx < int(sufficientSets.size()); ++sx) {
			RelationalRAMExpectation::sufficientSet &ss = sufficientSets[sx];
			placement &first = placements[ss.start * clumpSize];
			computeMeanCov(dataVec.segment(first.obsStart, ss.length * clumpObs),
				       clumpObs, ss.dataMean, ss.dataCov);
			if (st.getOptimizeMean() < 2) continue;
			for (int cx=0; cx < clumpSize; ++cx) {
				int gx = ss.start * clumpSize + cx;
				auto &sm = st.layoutSetup[gMap[gx]].skipMean;
				if (sm == NA_INTEGER) sm = 0;
			}
			for (int px=1; px < ss.length; ++px) {
				for (int cx=0; cx < clumpSize; ++cx) {
					int gx = (ss.start + px) * clumpSize + cx;
					auto &sm = st.layoutSetup[gMap[gx]].skipMean;
					if (sm == NA_INTEGER) sm = 1;
				}
			}
		}
		if (st.hasProductNodes) {
			skipMean = st.layoutSetup[gMap[0]].skipMean;
			if (skipMean) expectedVec.setZero();
		}
	}

	void independentGroup::computeCov(FitContext *fc)
	{
		if (0 == getParent().dataVec.size()) return;

		((ApcIO&)*pcalc.aio).useRampart = true;
		pcalc.cov(fc, fullCov);
	}

	void independentGroup::computeMean(FitContext *fc)
	{
		if (0 == getParent().dataVec.size()) return;

		if (skipMean) return;

		((ApcIO&)*pcalc.aio).useRampart = false;
		pcalc.mean(fc, expectedVec);
	}

	void independentGroup::simulate()
	{
		if (!dataVec.size()) return;

		simDataVec = expectedVec;

		SimpCholesky< Eigen::MatrixXd > covDecomp;
		Eigen::MatrixXd &denseCov = fullCov;
		covDecomp.compute(denseCov);
		if (covDecomp.info() != Eigen::Success || !(covDecomp.vectorD().array() > 0.0).all()) {
			omxRaiseErrorf("%s: covariance is non-positive definite", st.homeEx->name);
			return;
		}

		Eigen::MatrixXd res(fullCov.rows(), fullCov.cols());
		res.setIdentity();
		res = covDecomp.transpositionsP() * res;
		// L^* P
		res = covDecomp.matrixU() * res;
		// D(L^*P)
		res = covDecomp.vectorD().array().sqrt().matrix().asDiagonal() * res;

		Eigen::VectorXd sim1(clumpObs);
		int clumps = placements.size() / clumpSize;
		for (int cx=0; cx < clumps; ++cx) {
			for (int ob=0; ob < clumpObs; ++ob) {
				sim1[ob] = Rf_rnorm(0, 1.0);
			}
			simDataVec.segment(cx*clumpObs, clumpObs) += sim1.transpose() * res;
		}
	}

	struct SimUnitAccessor {
		state &st;
		SimUnitAccessor(state *_st) : st(*_st) {};

		double &operator() (const int unit, const int obs)
		{
			addr &ad = st.getParent().layout[unit];
			independentGroup &ig = *ad.ig;
			int obsStart = ig.placements[ad.igIndex].obsStart;
			return ig.simDataVec.coeffRef(obsStart + obs);
		};
	};

	void state::simulate(FitContext *fc, MxRList &out)
	{
		computeMean(fc);

		for (auto &ig : group) {
			ig->computeCov(fc);
			ig->simulate();
		}

		unapplyRotationPlan(SimUnitAccessor(this));

		std::map<omxExpectation*, SEXP> DataMap;
		for (auto &ex1 : allEx) {
			auto &dc = ex1->getDataColumns();
			if (dc.size() == 0) continue;
			omxData *data = ex1->data;
			int rows = data->nrows();

			SEXP df;
			Rf_protect(df = Rf_allocVector(VECSXP, dc.size()));
			SEXP colnames = Rf_allocVector(STRSXP, dc.size());
			Rf_setAttrib(df, R_NamesSymbol, colnames);
			for (int col=0; col < int(dc.size()); ++col) {
				SEXP colData = Rf_allocVector(REALSXP, rows);
				SET_VECTOR_ELT(df, col, colData);
				double *colPtr = REAL(colData);
				for (int rx=0; rx < rows; ++rx) {
					colPtr[rx] = NA_REAL;
				}
				SET_STRING_ELT(colnames, col, Rf_mkChar(omxDataColumnName(data, dc[col])));
			}
			markAsDataFrame(df, rows);

			DataMap[ex1] = df;
			out.add(data->name, df);
		}

		// NOTE: Does not copy foreign and primaryKeys

		for (auto &ig : group) {
			if (0 == ig->dataVec.size()) continue;
			for (int px=0, dx=0; px < int(ig->gMap.size()); ++px) {
				addr &a1 = layout[ ig->gMap[px] ];
				omxRAMExpectation *ram = a1.getRAMExpectationReadOnly();
				SEXP df = DataMap[ram];
				auto &pl = ig->placements[px];
				for (int vx=0, ncol=0; vx < ram->F->cols; ++vx) {
					if (!ram->latentFilter[vx]) continue;
					int col = ncol++;
					if (!ig->latentFilter[ pl.modelStart + vx ]) continue;

					REAL(VECTOR_ELT(df, col))[a1.row] = ig->simDataVec[dx];
					dx += 1;
				}
			}
		}
	}

	void state::computeCov(FitContext *fc)
	{
		for (size_t gx=0; gx < group.size(); ++gx) {
			independentGroup *ig = group[gx];
			ig->computeCov(fc);
		}
	}

	void independentGroup::filterFullMean()
	{
		// With optimizeMean, copies lots of extra data TODO
		independentGroup &pig = getParent();
		if (0 == pig.dataVec.size()) return;
		int ox = 0;
		for (size_t lx=0; lx < pig.latentFilter.size(); ++lx) {
			if (!pig.latentFilter[lx]) continue;
			expectedVec[ox++] = fullMean[lx];
		}
	}

	void state::computeMean(FitContext *fc)
	{
		if (hasProductNodes) {
			computeMeanByGroup(fc);
		} else {
			computeMeanByModel(fc);
		}

		state &pst = getParent();
		const int layoutSize = pst.layout.size();
		if (pst.getOptimizeMean() >= 1) {
			for (int ax=0; ax < layoutSize; ax += pst.layout[ax].nextMean) {
				addr &a1 = pst.layout[ax];
				int a1Start = a1.ig->placements[a1.igIndex].obsStart;
				independentGroup &tig1 = *group[a1.ig->arrayIndex];
				tig1.expectedVec.segment(a1Start, a1.numObs())
					*= a1.quickRotationFactor;
			}
		}
	}

	void state::computeMeanByGroup(FitContext *fc)
	{
		for (size_t gx=0; gx < group.size(); ++gx) {
			group[gx]->computeMean(fc);
		}
	}

	void state::computeMeanByModel(FitContext *fc)
	{
		// maybe there is a way to sort by dependency
		// so this loop can be parallelized

		state &pst = getParent();
		const int layoutSize = pst.layout.size();

		// can detect whether all units within an independent group are self contained
		for (int ax=0; ax < layoutSize; ax += pst.layout[ax].nextMean) {
			addr &a1 = pst.layout[ax];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

			omxData *data = expectation->data;
			expectation->loadDefVars(a1.row);
			int a1Start = a1.ig->placements[a1.igIndex].modelStart;
			independentGroup &tig1 = *group[a1.ig->arrayIndex];
			if (ram->M) {
				omxRecompute(ram->M, fc);
				EigenVectorAdaptor eM(ram->M);
				tig1.fullMean.segment(a1Start, a1.numVars()) = eM;
				if (0) {
					tig1.rawFullMean.segment(a1Start, a1.numVars()) = eM;
				}
			} else {
				tig1.fullMean.segment(a1Start, a1.numVars()).setZero();
			}

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxData *data1 = betA->getJoinModel()->data;
				int frow = data1->lookupRowOfKey(key);
				int a2Offset = pst.rowToLayoutMap[std::make_pair(data1, frow)];
				if (ax < a2Offset) mxThrow("Not in topological order");
				addr &a2 = pst.layout[a2Offset];
				independentGroup &tig2 = *group[a2.ig->arrayIndex];
				omxRecompute(betA, fc);
				EigenMatrixAdaptor eBA(betA);
				tig1.fullMean.segment(a1Start, a1.numVars()) +=
					eBA * tig2.fullMean.segment(a2.ig->placements[a2.igIndex].modelStart, eBA.cols());
			}

			// Cannot filter latents here because this could be an upper level
			// for subsequent units.

			expectation->loadDefVars(a1.row);
			omxRecompute(ram->A, fc);
			Eigen::VectorXd meanIn(tig1.fullMean.segment(a1Start, a1.numVars()));
			tig1.fullMean.segment(a1Start, a1.numVars()) = ram->pcalc.fullMean(fc, meanIn);
		}

		for (size_t gx=0; gx < group.size(); ++gx) {
			group[gx]->filterFullMean();
		}

		if (false) {
			size_t totalObserved = 0;
			for (size_t gx=0; gx < group.size(); ++gx) {
				totalObserved += group[gx]->dataVec.size();
			}
			Eigen::VectorXd expectedVec(totalObserved);
			int ox=0;
			for (size_t ax=0; ax < pst.layout.size(); ++ax) {
				addr &a1 = pst.layout[ax];
				int a1Start = a1.ig->placements[a1.igIndex].obsStart;
				independentGroup &tig1 = *group[a1.ig->arrayIndex];
				expectedVec.segment(ox, a1.numObs()) =
					tig1.expectedVec.segment(a1Start, a1.numObs());
			}
		}

		pst.applyRotationPlan(UnitAccessor<true>(this));
	}

	static int plusOne(int val) {
		if (val == NA_INTEGER) return val;
		else return val + 1;
	}

	void independentGroup::exportInternalState(MxRList &out, MxRList &dbg)
	{
		dbg.add("clumpSize", Rf_ScalarInteger(clumpSize));
		dbg.add("clumpObs", Rf_ScalarInteger(clumpObs));
		dbg.add("numLooseClumps", Rf_ScalarInteger(numLooseClumps()));
		dbg.add("skipMean", Rf_ScalarInteger(skipMean));
		if (st.hasProductNodes) {
			dbg.add("polyRep", Rcpp::wrap(pcalc.getPolyRep()));
		}

		if (clumpObs < 500) {
			// Can crash R because vectors are too long.
			// Maybe could allow more, but clumpObs==4600 is too much.
			if (expectedVec.size()) {
				SEXP m1 = Rcpp::wrap(expectedVec);
				Rf_protect(m1);
				Rf_setAttrib(m1, R_NamesSymbol, obsNameVec);
				out.add("mean", m1);
			}
			if (fullCov.nonZeros()) {
				out.add("covariance", Rcpp::wrap(fullCov));
			}
			SEXP fmean = Rcpp::wrap(fullMean);
			dbg.add("fullMean", fmean);
			Rf_setAttrib(fmean, R_NamesSymbol, varNameVec);
			if (0) {
				fmean = Rcpp::wrap(rawFullMean);
				dbg.add("rawFullMean", fmean);
				Rf_setAttrib(fmean, R_NamesSymbol, varNameVec);
			}
			//Eigen::SparseMatrix<double> A = getInputMatrix();
			//dbg.add("A", Rcpp::wrap(A));
			if (0) {
				// regularize internal representation
				//Eigen::SparseMatrix<double> fAcopy = asymT.IAF.transpose();
				//dbg.add("filteredA", Rcpp::wrap(fAcopy));
			}
			//Eigen::SparseMatrix<double> fullSymS = fullS.selfadjointView<Eigen::Lower>();
			//dbg.add("S", Rcpp::wrap(fullSymS));
			dbg.add("latentFilter", Rcpp::wrap(latentFilter));
			SEXP dv = Rcpp::wrap(dataVec);
			Rf_protect(dv);
			//Rf_setAttrib(dv, R_NamesSymbol, obsNameVec);
			dbg.add("dataVec", dv);
		} else {
			Rf_warning("%s: group %d too large to transfer to back to R",
				   st.homeEx->name, arrayIndex+1);
		}

		SEXP aIndex, modelStart, obsStart;
		Rf_protect(aIndex = Rf_allocVector(INTSXP, placements.size()));
		Rf_protect(modelStart = Rf_allocVector(INTSXP, placements.size()));
		Rf_protect(obsStart = Rf_allocVector(INTSXP, placements.size()));
		for (size_t mx=0; mx < placements.size(); ++mx) {
			INTEGER(aIndex)[mx] = 1 + gMap[mx];
			INTEGER(modelStart)[mx] = 1 + placements[mx].modelStart;
			INTEGER(obsStart)[mx] = 1 + placements[mx].obsStart;
		}
		SEXP layoutColNames, layoutDF;
		int numLayoutCols = 3;
		Rf_protect(layoutColNames = Rf_allocVector(STRSXP, numLayoutCols));
		SET_STRING_ELT(layoutColNames, 0, Rf_mkChar("aIndex"));
		SET_STRING_ELT(layoutColNames, 1, Rf_mkChar("modelStart"));
		SET_STRING_ELT(layoutColNames, 2, Rf_mkChar("obsStart"));
		Rf_protect(layoutDF = Rf_allocVector(VECSXP, numLayoutCols));
		Rf_setAttrib(layoutDF, R_NamesSymbol, layoutColNames);
		SET_VECTOR_ELT(layoutDF, 0, aIndex);
		SET_VECTOR_ELT(layoutDF, 1, modelStart);
		SET_VECTOR_ELT(layoutDF, 2, obsStart);
		markAsDataFrame(layoutDF, placements.size());
		dbg.add("layout", layoutDF);

		dbg.add("numSufficientSets", Rcpp::wrap(int(sufficientSets.size())));
		dbg.add("fit", Rcpp::wrap(fit));

		int digits = ceilf(log10f(sufficientSets.size()));
		std::string fmt = string_snprintf("ss%%0%dd", digits);
		for (size_t gx=0; gx < sufficientSets.size(); ++gx) {
			sufficientSet &ss = sufficientSets[gx];
			MxRList info;
			info.add("start", Rcpp::wrap(1 + ss.start));
			info.add("length", Rcpp::wrap(ss.length));
			info.add("mean", Rcpp::wrap(ss.dataMean));
			info.add("covariance", Rcpp::wrap(ss.dataCov));
			std::string name = string_snprintf(fmt.c_str(), int(1+gx));
			dbg.add(name.c_str(), info.asR());
		}
	}

	void state::exportInternalState(MxRList &dbg)
	{
		dbg.add("rampartUsage", Rcpp::wrap(rampartUsage));
		dbg.add("rotationCount", Rcpp::wrap(rotationCount));
		dbg.add("numGroups", Rcpp::wrap(int(group.size())));
		dbg.add("hasProductNodes", Rcpp::wrap(hasProductNodes));

		SEXP modelName, row, numJoins, numKids, parent1, fk1, rscale, qrf,
			hmean, skipMean, rset, ugroup;
		Rf_protect(modelName = Rf_allocVector(STRSXP, layout.size()));
		Rf_protect(row = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numKids = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numJoins = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(parent1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(fk1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(rscale = Rf_allocVector(REALSXP, layout.size()));
		Rf_protect(qrf = Rf_allocVector(REALSXP, layout.size()));
		Rf_protect(hmean = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(skipMean = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(rset = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(ugroup = Rf_allocVector(INTSXP, layout.size()));
		for (size_t mx=0; mx < layout.size(); ++mx) {
			SET_STRING_ELT(modelName, mx, Rf_mkChar(layout[mx].modelName().c_str()));
			INTEGER(row)[mx] = 1+layout[mx].row;
			INTEGER(numKids)[mx] = layoutSetup[mx].numKids;
			INTEGER(numJoins)[mx] = layoutSetup[mx].numJoins;
			INTEGER(parent1)[mx] = plusOne(layoutSetup[mx].parent1);
			INTEGER(fk1)[mx] = layoutSetup[mx].fk1;
			REAL(rscale)[mx] = layout[mx].rampartScale;
			REAL(qrf)[mx] = layout[mx].quickRotationFactor;
			INTEGER(hmean)[mx] = layoutSetup[mx].heterogenousMean;
			INTEGER(skipMean)[mx] = layoutSetup[mx].skipMean;
			INTEGER(rset)[mx] = layoutSetup[mx].rset;
			INTEGER(ugroup)[mx] = layout[mx].ig? 1+layout[mx].ig->arrayIndex : NA_INTEGER;
		}
		dbg.add("layout", Rcpp::DataFrame::create(Rcpp::Named("model")=modelName,
							  Rcpp::Named("row")=row,
							  Rcpp::Named("numKids")=numKids,
							  Rcpp::Named("numJoins")=numJoins,
							  Rcpp::Named("parent1")=parent1,
							  Rcpp::Named("fk1")=fk1,
							  Rcpp::Named("rampartScale")=rscale,
							  Rcpp::Named("qrf")=qrf,
							  Rcpp::Named("hmean")=hmean,
							  Rcpp::Named("skip")=skipMean,
							  Rcpp::Named("rset")=rset,
							  Rcpp::Named("group")=ugroup));

		int digits = ceilf(log10f(1.0+group.size()));
		std::string fmt = string_snprintf("g%%0%dd", digits);
		size_t maxIndex = std::min(group.size(),
					size_t(((omxRAMExpectation*)homeEx)->maxDebugGroups));
		for (size_t gx=0; gx < maxIndex; ++gx) {
			independentGroup &ig = *group[gx];
			MxRList info;
			ig.exportInternalState(info, info);
			std::string name = string_snprintf(fmt.c_str(), int(1+gx));
			dbg.add(name.c_str(), info.asR());
		}
	}
}
