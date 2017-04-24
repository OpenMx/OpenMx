/*
 * Copyright 2007-2017 The OpenMx Project
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

void omxRAMExpectation::compute(FitContext *fc, const char *what, const char *how)
{
	omxRAMExpectation* oro = this;

	if (what && how && strEQ(how, "flat")) {
		bool wantCov = false;
		bool wantMean = false;
		if (strEQ(what, "distribution")) { wantCov = true; wantMean = true; }
		if (strEQ(what, "covariance")) wantCov = true;
		if (strEQ(what, "mean")) wantMean = true;
		if (!oro->rram) {
			oro->rram = new RelationalRAMExpectation::state;
			oro->rram->init(this, fc);
		}
		if (wantCov)  oro->rram->computeCov(fc);
		if (wantMean) oro->rram->computeMean(fc);
		return;
	}

	oro->CalculateRAMCovarianceAndMeans(fc);
}

omxRAMExpectation::~omxRAMExpectation()
{
	if(OMX_DEBUG) { mxLog("Destroying RAM Expectation."); }
	
	omxRAMExpectation* argStruct = this;

	if (argStruct->rram) delete argStruct->rram;

	omxFreeMatrix(argStruct->cov);

	if(argStruct->means != NULL) {
		omxFreeMatrix(argStruct->means);
	}

	omxFreeMatrix(argStruct->I);
	omxFreeMatrix(argStruct->X);
	omxFreeMatrix(argStruct->Y);
	omxFreeMatrix(argStruct->Ax);
	omxFreeMatrix(_Z);
}

static void refreshUnfilteredCov(omxExpectation *oo)
{
	// Ax = ZSZ' = Covariance matrix including latent variables
	omxRAMExpectation* oro = (omxRAMExpectation*) oo;
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

void omxRAMExpectation::populateAttr(SEXP robj)
{
    refreshUnfilteredCov(this);
    omxRAMExpectation* oro = this;
	
	{
		ProtectedSEXP expCovExt(Rf_allocMatrix(REALSXP, Ax->rows, Ax->cols));
		memcpy(REAL(expCovExt), Ax->data, sizeof(double) * Ax->rows * Ax->cols);
		Rf_setAttrib(robj, Rf_install("UnfilteredExpCov"), expCovExt);
	}
	Rf_setAttrib(robj, Rf_install("numStats"), Rf_ScalarReal(omxDataDF(data)));

	MxRList out;
	MxRList dbg;

	if (oro->rram) {
		rram->exportInternalState(dbg);
	} else {
		oro->CalculateRAMCovarianceAndMeans(0);
		EigenMatrixAdaptor Ecov(oro->cov);
		out.add("covariance", Rcpp::wrap(Ecov));
		if (oro->means) {
			EigenVectorAdaptor Emean(oro->means);
			out.add("mean", Rcpp::wrap(Emean));
		}
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

void omxRAMExpectation::CalculateRAMCovarianceAndMeans(FitContext *fc)
{
	if (F->rows == 0) return;

	omxRecompute(A, fc);
	omxRecompute(S, fc);
	omxRecompute(F, fc);
	if (M) omxRecompute(M, fc);
	    
	if(OMX_DEBUG) { mxLog("Running RAM computation with numIters is %d\n.", numIters); }
		
	if(Ax == NULL || I == NULL || Y == NULL || X == NULL) {
		Rf_error("Internal Error: RAM Metadata improperly populated.  Please report this to the OpenMx development team.");
	}
		
	if(cov == NULL && means == NULL) {
		return; // We're not populating anything, so why bother running the calculation?
	}
	
	omxMatrix *Z = getZ(NULL);
	EigenMatrixAdaptor eZ(Z);
	EigenMatrixAdaptor eY(Y);
	for (int rx=0, dx=0; rx < eZ.rows(); ++rx) {
		if (!latentFilter[rx]) continue;
		eY.row(dx) = eZ.row(rx);
		dx += 1;
	}

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

omxExpectation *omxInitRAMExpectation() { return new omxRAMExpectation; }

void omxRAMExpectation::init() {
	if(OMX_DEBUG) { mxLog("Initializing RAM expectation."); }
	
	int l, k;

	SEXP slotValue;
	
	omxMatrix *Zmat = omxInitMatrix(0, 0, TRUE, currentState);
	_Z = Zmat;
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

	/* Identity Matrix, Size Of A */
	if(OMX_DEBUG) { mxLog("Generating I."); }
	RAMexp->I = omxNewIdentityMatrix(RAMexp->A->rows, currentState);

	if(OMX_DEBUG) { mxLog("Processing expansion iteration depth."); }
	{ScopedProtect p1(slotValue, R_do_slot(rObj, Rf_install("depth")));
	RAMexp->numIters = INTEGER(slotValue)[0];
	if(OMX_DEBUG) { mxLog("Using %d iterations.", RAMexp->numIters); }
	}

	ProtectedSEXP Rrampart(R_do_slot(rObj, Rf_install(".rampartCycleLimit")));
	RAMexp->rampartCycleLimit = Rf_asInteger(Rrampart);

	ProtectedSEXP RrampartLimit(R_do_slot(rObj, Rf_install(".rampartUnitLimit")));
	RAMexp->rampartUnitLimit = Rf_asInteger(RrampartLimit);

	RAMexp->useSufficientSets = true;
	if (R_has_slot(rObj, Rf_install(".useSufficientSets"))) {
		ProtectedSEXP Rss(R_do_slot(rObj, Rf_install(".useSufficientSets")));
		RAMexp->useSufficientSets = Rf_asLogical(Rss);
	}

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
			if (!fex) Rf_error("%s: level transition matrix '%s' does not reference the upper level model",
					   oo->name, bmat->name());
			omxCompleteExpectation(fex);
			if (!strEQ(fex->expType, "MxExpectationRAM")) {
				Rf_error("%s: only MxExpectationRAM can be joined with MxExpectationRAM", oo->name);
			}
			omxDataKeysCompatible(fex->data, oo->data, foreignKey);
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

	RAMexp->studyF();
	//mxPrintMat("RAM corrected dc", oo->getDataColumns());
}

void omxRAMExpectation::studyF()
{
	auto dataColumns = super::getDataColumns();
	auto origThresholdInfo = super::getThresholdInfo();
	EigenMatrixAdaptor eF(F);
	latentFilter.assign(eF.cols(), false);
	dataCols.resize(eF.rows());
	if (!eF.rows()) return;  // no manifests
	for (int cx =0, dx=0; cx < eF.cols(); ++cx) {
		int dest;
		double isManifest = eF.col(cx).maxCoeff(&dest);
		latentFilter[cx] = isManifest;
		if (isManifest) {
			int newDest = dataColumns.size()? dataColumns[dest] : dest;
			dataCols[dx] = newDest;
			if (origThresholdInfo.size()) {
				omxThresholdColumn adj = origThresholdInfo[dest];
				adj.dColumn = dx;
				thresholds.push_back(adj);
			}
			dx += 1;
		}
	}
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
	} else if(strEQ("pvec", component)) {
		// Once implemented, change compute function and return pvec
	}
	
	return retval;
}

namespace RelationalRAMExpectation {

	omxExpectation *addr::getModel(FitContext *fc)
	{
		return omxExpectationFromIndex(model->expNum, fc->state);
	}

	omxRAMExpectation *addr::getRAMExpectation(FitContext *fc)
	{
		return (omxRAMExpectation*) getModel(fc);
	};

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
					Rf_error("Cannot find row %d in %s", row, e1->data->name);
				int bx = it->second;
				cc.connect(ax, bx);
			}
		}
	}

	void independentGroup::refreshUnitA(FitContext *fc, int px)
	{
		independentGroup &par = getParent();
		struct placement &pl = par.placements[px];
		addr &a1 = par.st.layout[ par.gMap[px] ];
		omxExpectation *expectation = a1.getModel(fc);
		omxData *data = expectation->data;
		omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

		EigenMatrixAdaptor eA(ram->A);
		for (int cx=0; cx < eA.cols(); ++cx) {
			for (int rx=0; rx < eA.rows(); ++rx) {
				double val = eA(rx, cx);
				if (val != 0) {
					if (rx == cx) {
						Rf_error("%s: nonzero diagonal entry in A matrix at %d",
							 st.homeEx->name, 1+pl.modelStart+cx);
					}
					asymT.fullA.coeffRef(pl.modelStart + cx, pl.modelStart + rx) =
						asymT.getSign() * val;
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
			RowToPlacementMapType::iterator plIndex =
				par.rowToPlacementMap.find(std::make_pair(data1, frow));
			placement &p2 = par.placements[ plIndex->second ];
			omxRecompute(betA, fc);
			omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel();
			for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
				for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
					double val = omxMatrixElement(betA, rx, cx);
					if (val == 0.0) continue;
					asymT.fullA.coeffRef(p2.modelStart + cx, pl.modelStart + rx) =
						asymT.getSign() * val * scale;
				}
			}
		}
	}

	void independentGroup::refreshModel(FitContext *fc)
	{
		independentGroup &par = getParent();
		for (int ax=0; ax < clumpSize; ++ax) {
			placement &pl = par.placements[ax];
			addr &a1 = par.st.layout[ par.gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			expectation->loadDefVars(a1.row);
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
				if (it->second < 0) Rf_error("%s cycle detected: '%s' row %d joins against itself",
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
		a1.rampartScale = 1.0;
		a1.row = frow;
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

	void independentGroup::determineShallowDepth(FitContext *fc)
	{
		if (!Global->RAMInverseOpt) return;

		for (int ax=0; ax < clumpSize; ++ax) {
			placement &pl = placements[ax];
			addr &a1 = st.layout[ gMap[ax] ];
			omxExpectation *expectation = a1.getModel(fc);
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;
			omxData *data = expectation->data;

			expectation->loadDefVars(a1.row);
			omxRecompute(ram->A, fc);

			if (a1.rampartScale != 0.0) {
				for (size_t jx=0; jx < ram->between.size(); ++jx) {
					omxMatrix *betA = ram->between[jx];
					int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
					if (key == NA_INTEGER) continue;
					omxData *data1 = betA->getJoinModel()->data;
					int frow = data1->lookupRowOfKey(key);
					RowToPlacementMapType::iterator plIndex =
						rowToPlacementMap.find(std::make_pair(data1, frow));
					if (plIndex == rowToPlacementMap.end()) Rf_error("Cannot find row %d in %s",
											 frow, data1->name);
					placement &p2 = placements[ plIndex->second ];
					omxRecompute(betA, fc);
					betA->markPopulatedEntries();
					omxRAMExpectation *ram2 = (omxRAMExpectation*) betA->getJoinModel();
					for (int rx=0; rx < ram->A->rows; ++rx) {  //lower
						for (int cx=0; cx < ram2->A->rows; ++cx) {  //upper
							double val = omxMatrixElement(betA, rx, cx);
							if (val == 0.0) continue;
							asymT.fullA.coeffRef(p2.modelStart + cx, pl.modelStart + rx) = 1;
						}
					}
				}
			}

			ram->A->markPopulatedEntries();
			EigenMatrixAdaptor eA(ram->A);
			for (int cx=0; cx < eA.cols(); ++cx) {
				for (int rx=0; rx < eA.rows(); ++rx) {
					if (rx != cx && eA(rx,cx) != 0) {
						asymT.fullA.coeffRef(pl.modelStart + cx, pl.modelStart + rx) = 1;
					}
				}
			}
		}

		asymT.determineShallowDepth(fc);

		if (st.verbose() >= 1) {
			mxLog("%s: RAM shallow inverse depth = %d", st.homeEx->name, asymT.getDepth());
		}
	}

	std::vector<bool> &addr::getIgnoreDefVar()
	{
		omxRAMExpectation *ram = (omxRAMExpectation*) model;
		return ram->ignoreDefVar;
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

		bool compareModelAndMissingness(addr &la, addr &ra, bool &mismatch) const
		{
			mismatch = true;
			if (la.getExpNum() != ra.getExpNum())
				return la.getExpNum() < ra.getExpNum();

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

			mismatch = la.rampartScale != ra.rampartScale;
			if (mismatch) return la.rampartScale < ra.rampartScale;

			mismatch = false;
			return false;
		}

		bool compareAllDefVars(const addr &la, const addr &ra, bool &mismatch) const
		{
			mismatch = true;

			omxData *data = la.getData();  // both la & ra have same data
			for (size_t k=0; k < data->defVars.size(); ++k) {
				int col = data->defVars[k].column;
				double lv = omxDoubleDataElement(data, la.row, col);
				double rv = omxDoubleDataElement(data, ra.row, col);
				if (lv == rv) continue;
				return lv < rv;
			}

			mismatch = false;
			return false;
		}

		bool compareDefVars(addr &la, addr &ra, bool &mismatch) const
		{
			mismatch = true;

			std::vector<bool> &ignoreDefVar = la.getIgnoreDefVar();
			omxData *data = la.getData();  // both la & ra have same data
			for (size_t k=0; k < data->defVars.size(); ++k) {
				if (ignoreDefVar[k]) continue;
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
			for (size_t ux=0; ux < lhs.size(); ++ux) {
				addr &la = st.layout[lhs[ux]];
				addr &ra = st.layout[rhs[ux]];
				bool mismatch;
				bool got = compareModelAndMissingness(la, ra, mismatch);
				if (mismatch) return got;

				got = compareDefVars(la, ra, mismatch);
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
				Rf_error("Cannot find row %d in %s", row, e1->data->name);
			return &st.layout[it->second];
		}

		bool compareDeep(const addr &la, const addr &ra, bool &mismatch) const
		{
			omxRAMExpectation *ram = la.getRAMExpectationReadOnly();
			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				addr *lp = joinedWith(la, jx);
				addr *rp = joinedWith(ra, jx);
				if (!lp && !rp) continue;
				bool got = compareDeep(*lp, *rp, mismatch);
				if (mismatch) return got;
			}

			// rampartScale also affects covariance so no need to check here
			if (la.rampartScale != ra.rampartScale) {
				mismatch = true;
				return la.rampartScale < ra.rampartScale;
			}

			bool got = compareAllDefVars(la, ra, mismatch);
			if (mismatch) return got;

			return false;
		}

		bool operator() (const std::vector<int> &lhs, const std::vector<int> &rhs) const
		{
			bool mismatch = false;
			for (size_t ux=0; ux < lhs.size(); ++ux) {
				addr &la = st.layout[lhs[ux]];
				addr &ra = st.layout[rhs[ux]];
				bool got = compareDeep(la, ra, mismatch);
				if (mismatch) return got;
			}
			return false;
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
				    omxRAMExpectation *from, bool within)
	{
		to->hasVariance += (transition * from->hasVariance).array().abs().matrix();
		for (int rx=0; rx < transition.rows(); ++rx) {
			for (int cx=0; cx < transition.cols(); ++cx) {
				if (within && rx == cx) continue;
				if (transition(rx,cx) == 0) continue;
				dvScoreboardSetType &dv2 = from->dvScoreboard[cx];

				bool notIgnoreable = from->hasVariance[cx] && !within;
				if (verbose() >= 1) {
					for (dvScoreboardSetType::iterator it = dv2.begin();
					     it != dv2.end(); ++it) {
						omxExpectation *ex3 = it->first;
						omxRAMExpectation *ram3 = (omxRAMExpectation*) ex3;
						omxDefinitionVar &dv = ex3->data->defVars[ it->second ];
						mxLog("%s at %s[%d,%d] goes from %s to %s => %d (0=no cov effect)",
						      omxDataColumnName(ex3->data, dv.column),
						      ram3->M->name(), 1+dv.row, 1+dv.col,
						      from->S->rownames[cx], to->S->rownames[rx],
						      notIgnoreable);
					}
				}
				
				if (notIgnoreable) {
					for (dvScoreboardSetType::iterator it = dv2.begin();
					     it != dv2.end(); ++it) {
						omxExpectation *ex3 = it->first;
						omxRAMExpectation *ram3 = (omxRAMExpectation*) ex3;
						ram3->ignoreDefVar[ it->second ] = false;
					}
				} else {
					dvScoreboardSetType &dv1 = to->dvScoreboard[rx];
					dv1.insert(dv2.begin(), dv2.end());
				}
			}
		}
	}

	void state::identifyZeroVarPred(FitContext *fc)
	{
		for (std::set<omxExpectation*>::iterator it = allEx.begin() ; it != allEx.end(); ++it) {
			omxRAMExpectation *ram = (omxRAMExpectation*) (*it);
			omxData *data = (*it)->data;
			data->loadFakeData((*it)->currentState, 1.0);
			ram->S->markPopulatedEntries();
			omxRecompute(ram->S, fc);
			EigenMatrixAdaptor eS(ram->S);
			ram->hasVariance = eS.diagonal().array().abs().matrix();
			ram->dvScoreboard.resize(eS.rows());
		}

		std::set<omxExpectation*> checkedEx;
		for (size_t ax=0; ax < layout.size(); ++ax) {
			addr &a1 = layout[ax];
			omxExpectation *expectation = a1.getModel(fc);
			omxData *data = expectation->data;
			omxRAMExpectation *ram = (omxRAMExpectation*) expectation;

			if (ram->M) {
				int mNum = ~ram->M->matrixNumber;
				for (size_t k=0; k < data->defVars.size(); ++k) {
					omxDefinitionVar &dv = data->defVars[k];
					if (dv.matrix != mNum || ram->hasVariance[ dv.col ]) {
						continue;
					}
					ram->dvScoreboard[dv.col].insert(std::make_pair(expectation, k));
					ram->ignoreDefVar[k] = true; // tentative
					if (verbose() >= 1) {
						mxLog("consider def var %s at %s[%d,%d]",
						      omxDataColumnName(data, dv.column),
						      ram->M->name(), 1+dv.row, 1+dv.col);
					}
				}
			}

			for (size_t jx=0; jx < ram->between.size(); ++jx) {
				omxMatrix *betA = ram->between[jx];
				int key = omxKeyDataElement(data, a1.row, betA->getJoinKey());
				if (key == NA_INTEGER) continue;
				omxExpectation *ex2 = betA->getJoinModel();
				omxRAMExpectation *ram2 = (omxRAMExpectation*) ex2;
				omxRecompute(betA, fc);
				EigenMatrixAdaptor eBA(betA);
				propagateDefVar(ram, eBA, ram2, false);
			}

			ram->A->markPopulatedEntries();
			omxRecompute(ram->A, fc);
			EigenMatrixAdaptor Zmat(ram->getZ(fc));
			propagateDefVar(ram, Zmat, ram, true);

			checkedEx.insert(expectation);
			if (checkedEx.size() == allEx.size()) break;
		}

		for (std::set<omxExpectation*>::iterator it = allEx.begin() ; it != allEx.end(); ++it) {
			omxRAMExpectation *ram = (omxRAMExpectation*) (*it);
			ram->dvScoreboard.clear();
		}
	}

	template <typename T>
	void state::placeSet(std::set<std::vector<T> > &toPlace, independentGroup *ig)
	{
		for (std::set<std::vector<int> >::iterator px = toPlace.begin();
		     px != toPlace.end(); ++px) {
			const std::vector<int> &clump = *px;
			for (size_t cx=0; cx < clump.size(); ++cx) {
				ig->place(clump[cx]);
			}
		}
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

		if (verbose() >= 2) {
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
			ig->sufficientSets.resize(ssCount);
			int ssIndex = 0;
			for (CompatibleMeanMapType::iterator mit = cmm.begin();
			     mit != cmm.end(); ++mit) {
				if (mit->second.size() == 1) continue;
				int from = ig->placements.size();
				placeSet(mit->second, ig);
				if (verbose() >= 3) {
					mxLog("group %d same mean %d -> %d clumpsize %d",
					      int(group.size()), from, int(ig->placements.size() - 1),
					      int(it->second.begin()->size()));
				}
				ig->sufficientSets[ssIndex].start = from / ig->clumpSize;
				ig->sufficientSets[ssIndex].length = (ig->placements.size() - from) / ig->clumpSize;
				++ssIndex;
			}
			if (!ram->useSufficientSets) ig->sufficientSets.clear();
			ig->prep(fc);
			group.push_back(ig);
		}
	}

	void independentGroup::place(int ax)
	{
		if (st.layout[ax].ig) {
			Rf_error("Unit[%d] already assigned; this is a bug", ax);
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

	independentGroup::independentGroup(independentGroup *ig)
		: st(ig->st), clumpSize(ig->clumpSize),
		  analyzedCov(false), asymT(ig->latentFilter)
	{
		arrayIndex = ig->arrayIndex;
		obsNameVec = 0;
		varNameVec = 0;
		expectedVec.resize(ig->expectedVec.size());
		fullMean.resize(ig->fullMean.size());
		clumpVars = ig->clumpVars;
		clumpObs = ig->clumpObs;
		asymT.resize(clumpVars, clumpObs);
		asymT.setDepth(ig->asymT.getDepth());
	}

	void independentGroup::prep(FitContext *fc)
	{
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
		obsNameVec = Rf_protect(Rf_allocVector(STRSXP, totalObserved));
		varNameVec = Rf_protect(Rf_allocVector(STRSXP, maxSize));
		expectedVec.resize(totalObserved);
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
					if (!ram->latentFilter[vx]) continue;
					int col = ncol++;
					double val = omxMatrixElement(st.smallCol, 0, col);
					bool yes = std::isfinite(val);
					if (!yes) continue;
					latentFilter[ pl.modelStart + vx ] = true;
					std::string dname =
						modelName + omxDataColumnName(data, dc[col]);
					SET_STRING_ELT(obsNameVec, dx, Rf_mkChar(dname.c_str()));
					dataVec[ dx ] = val;
					if (a1.getExpNum() == st.homeEx->expNum) dataColumn[ dx ] = col;
					dx += 1;
				}
				if (a1.numObs() != dx - prevDx) {
					Rf_error("numObs() %d != %d", a1.numObs(), dx - prevDx);
				}
			}
			for (int vx=0; vx < ram->F->cols; ++vx) {
				std::string dname = modelName + ram->F->colnames[vx];
				SET_STRING_ELT(varNameVec, pl.modelStart + vx, Rf_mkChar(dname.c_str()));
			}
		}

		asymT.resize(clumpVars, clumpObs);
		determineShallowDepth(fc);
	}

	struct RampartCompareLib {
		state *st;
		RampartCompareLib(state *_st) : st(_st) {};

		omxExpectation *getJoinModel(const addr *a1) const {
			omxRAMExpectation *ram = a1->getRAMExpectationReadOnly();
			omxMatrix *b1 = ram->between[0];
			return b1->getJoinModel();
		};

		// actually stores !missingness
		template <typename T>
		void getMissingnessPattern(const addr *a1, std::vector<T> &out) const
		{
			a1->dataRow(st->smallCol);
			int jCols = a1->getDataColumns().size();
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

			omxRAMExpectation *ram = lhs->getRAMExpectationReadOnly();
			omxMatrix *b1 = ram->between[0];

			bool mismatch;
			result = lhs->getData()->CompareDefVarInMatrix(lhs->row, rhs->row, b1, mismatch);
			if (mismatch) return true;

			return false;
		}

		bool cmpRecursive(const addr *lhs, const addr *rhs, bool &result) const
		{
			if (lhs->getExpNum() != rhs->getExpNum()) {
				result = lhs->getExpNum() < rhs->getExpNum();
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

			bool mismatch = lhs->rampartScale != rhs->rampartScale;
			if (mismatch) {
				result = lhs->rampartScale < rhs->rampartScale;
				return true;
			}

			const addrSetup *lhss = &st->layoutSetup[lhs - &st->layout[0]];
			const addrSetup *rhss = &st->layoutSetup[rhs - &st->layout[0]];
			if (lhss->clump.size() != rhss->clump.size()) {
				result = lhss->clump.size() < rhss->clump.size();
				return true;
			}
			for (size_t cx=0; cx < lhss->clump.size(); ++cx) {
				if (cmpRecursive(&st->layout[lhss->clump[cx]], &st->layout[rhss->clump[cx]], result))
					return true;
			}
			return false;
		}
	};

	struct RampartTodoCompare : RampartCompareLib {
		RampartTodoCompare(state *_st) : RampartCompareLib(_st) {};

		bool operator() (const addr *lhs, const addr *rhs) const
		{
			bool result = false;
			if (cmpUpper(lhs, rhs, result)) return result;

			const addrSetup *lhss = &st->layoutSetup[lhs - &st->layout[0]];
			const addrSetup *rhss = &st->layoutSetup[rhs - &st->layout[0]];

			if (lhss->fk1 != rhss->fk1)
				return lhss->fk1 < rhss->fk1;
			cmpRecursive(lhs, rhs, result);
			return result;
		}
	};

	struct RampartClumpCompare : RampartCompareLib {
		RampartClumpCompare(state *_st) : RampartCompareLib(_st) {};

		bool clumpCmp(const int lhs, const int rhs) const {
			const addr *lhsObj = &st->layout[lhs];
			const addr *rhsObj = &st->layout[rhs];
			bool result = false;
			if (cmpRecursive(lhsObj, rhsObj, result)) return result;
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
		addrSetup &specimen = layoutSetup[ t1[0] ];
		for (size_t cx=0; cx < specimen.clump.size(); ++cx) {
			std::vector<int> t2;
			t2.reserve(t1.size());
			for (size_t tx=0; tx < t1.size(); ++tx) {
				addrSetup &a1 = layoutSetup[ t1[tx] ];
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

		int loopTo = layout.size();
		int rampartUnitLimit = ((omxRAMExpectation*) homeEx)->rampartUnitLimit;
		if (rampartUnitLimit != NA_INTEGER) loopTo = std::min(rampartUnitLimit, loopTo);
		for (int ax=0; ax < loopTo; ++ax) {
			addr &a1 = layout[ax];
			addrSetup &as1 = layoutSetup[ax];
			if (as1.numKids != 0 || as1.numJoins != 1 || as1.clumped) continue;
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
					if (dc[cx] != a1.clump[cx]) Rf_error("oops");
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
			group.reserve(parent->group.size());
			for (size_t gx=0; gx < parent->group.size(); ++gx) {
				independentGroup *ig = parent->group[gx];
				group.push_back(new independentGroup(ig));
			}
			return;
		}

		{
			ProtectedSEXP Rdvhack(R_do_slot(expectation->rObj, Rf_install(".identifyZeroVarPred")));
			doIdentifyZeroVarPred = Rf_asLogical(Rdvhack);
			if (verbose()) {
				mxLog("%s: identifyZeroVarPred=%d", homeEx->name, doIdentifyZeroVarPred);
			}
		}

		int maxSize = 0;
		int homeDataRows = homeEx->data->rows;
		for (int row=0; row < homeDataRows; ++row) {
			flattenOneRow(homeEx, row, maxSize);
			if (isErrorRaised()) return;
		}
		for (auto *ex : allEx) {
			if (!ex->data->hasWeight()) continue;
			Rf_error("%s: row weights provided in '%s' are not compatible with joins",
				 expectation->name, ex->data->name);
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

		Eigen::VectorXd vec(fc->varGroup->vars.size());
		vec.setConstant(1);
		copyParamToModelInternal(fc->varGroup, homeEx->currentState, vec.data());

		for (std::set<omxExpectation*>::iterator it = allEx.begin() ; it != allEx.end(); ++it) {
			omxRAMExpectation *ram2 = (omxRAMExpectation*) (*it);
			ram2->ignoreDefVar.assign((*it)->data->defVars.size(), false);
		}

		if (doIdentifyZeroVarPred) identifyZeroVarPred(fc);

		planModelEval(maxSize, fc);

		fc->copyParamToModelClean();

		for (std::vector<independentGroup*>::iterator it = group.begin() ; it != group.end(); ++it) {
			(*it)->arrayIndex = it - group.begin();
		}

		applyRotationPlan(UnitAccessor<false>(this));

		for (std::vector<independentGroup*>::iterator it = group.begin() ; it != group.end(); ++it) {
			(*it)->finalizeData();
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
		}
	}

	void independentGroup::computeCov1(FitContext *fc)
	{
		if (0 == getParent().dataVec.size()) return;

		fullS.conservativeResize(clumpVars, clumpVars);

		refreshModel(fc);
	}

	void independentGroup::computeCov2()
	{
		//{ Eigen::MatrixXd tmp = fullA; mxPrintMat("fullA", tmp); }
		//{ Eigen::MatrixXd tmp = fullS; mxPrintMat("fullS", tmp); }

		asymT.invert();
		asymT.filter();

		//mxPrintMat("S", fullS);
		//Eigen::MatrixXd fullCovDense =

		// IAF tends to be very sparse so we want to do this quadratic
		// product using sparse matrices. However, the result is typically
		// fairly dense.

		fullCov = (asymT.IAF.transpose() * fullS.selfadjointView<Eigen::Lower>() * asymT.IAF);
		//mxLog("fullCov %d%% nonzero", int(fullCov.nonZeros() * 100.0 / (fullCov.rows() * fullCov.cols())));
		//{ Eigen::MatrixXd tmp = fullCov; mxPrintMat("fullcov", tmp); }
	}

	void state::computeCov(FitContext *fc)
	{
		for (size_t gx=0; gx < group.size(); ++gx) {
			independentGroup *ig = group[gx];
			ig->computeCov1(fc);
			ig->computeCov2();
		}
	}

	void independentGroup::filterFullMean()
	{
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
		// maybe there is a way to sort by dependency
		// so this loop can be parallelized

		state &pst = getParent();

		// can detect whether all units within an independent group are self contained
		for (size_t ax=0; ax < pst.layout.size(); ++ax) {
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
				size_t a2Offset = pst.rowToLayoutMap[std::make_pair(data1, frow)];
				if (ax < a2Offset) Rf_error("Not in topological order");
				addr &a2 = pst.layout[a2Offset];
				independentGroup &tig2 = *group[a2.ig->arrayIndex];
				omxRecompute(betA, fc);
				EigenMatrixAdaptor eBA(betA);
				tig1.fullMean.segment(a1Start, a1.numVars()) +=
					eBA * tig2.fullMean.segment(a2.ig->placements[a2.igIndex].modelStart, eBA.cols());
			}

			expectation->loadDefVars(a1.row);
			omxRecompute(ram->A, fc);
			EigenMatrixAdaptor eZ(ram->getZ(fc));
			tig1.fullMean.segment(a1Start, a1.numVars()) =
				eZ * tig1.fullMean.segment(a1Start, a1.numVars());
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
				expectedVec.segment(ox, a1.numObs()) =
					a1.ig->expectedVec.segment(a1Start, a1.numObs());
			}
		}

		pst.applyRotationPlan(UnitAccessor<true>(this));
	}

	static int plusOne(int val) {
		if (val == NA_INTEGER) return val;
		else return val + 1;
	}

	Eigen::SparseMatrix<double> independentGroup::getInputMatrix() const
	{
		return asymT.getSign() * asymT.fullA.transpose();
	}

	void independentGroup::exportInternalState(MxRList &out, MxRList &dbg)
	{
		dbg.add("clumpSize", Rf_ScalarInteger(clumpSize));
		dbg.add("clumpObs", Rf_ScalarInteger(clumpObs));
		dbg.add("numLooseClumps", Rf_ScalarInteger(numLooseClumps()));

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
			Eigen::SparseMatrix<double> A = getInputMatrix();
			dbg.add("A", Rcpp::wrap(A));
			if (0) {
				// regularize internal representation
				Eigen::SparseMatrix<double> fAcopy = asymT.IAF.transpose();
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
		dbg.add("numGroups", Rcpp::wrap(int(group.size())));

		SEXP modelName, row, numJoins, numKids, parent1, fk1, rscale, ugroup;
		Rf_protect(modelName = Rf_allocVector(STRSXP, layout.size()));
		Rf_protect(row = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numKids = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(numJoins = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(parent1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(fk1 = Rf_allocVector(INTSXP, layout.size()));
		Rf_protect(rscale = Rf_allocVector(REALSXP, layout.size()));
		Rf_protect(ugroup = Rf_allocVector(INTSXP, layout.size()));
		for (size_t mx=0; mx < layout.size(); ++mx) {
			SET_STRING_ELT(modelName, mx, Rf_mkChar(layout[mx].modelName().c_str()));
			INTEGER(row)[mx] = 1+layout[mx].row;
			INTEGER(numKids)[mx] = layoutSetup[mx].numKids;
			INTEGER(numJoins)[mx] = layoutSetup[mx].numJoins;
			INTEGER(parent1)[mx] = plusOne(layoutSetup[mx].parent1);
			INTEGER(fk1)[mx] = layoutSetup[mx].fk1;
			REAL(rscale)[mx] = layout[mx].rampartScale;
			INTEGER(ugroup)[mx] = layout[mx].ig? 1+layout[mx].ig->arrayIndex : NA_INTEGER;
		}
		dbg.add("layout", Rcpp::DataFrame::create(Rcpp::Named("model")=modelName,
							  Rcpp::Named("row")=row,
							  Rcpp::Named("numKids")=numKids,
							  Rcpp::Named("numJoins")=numJoins,
							  Rcpp::Named("parent1")=parent1,
							  Rcpp::Named("fk1")=fk1,
							  Rcpp::Named("rampartScale")=rscale,
							  Rcpp::Named("group")=ugroup));

		int digits = ceilf(log10f(group.size()));
		std::string fmt = string_snprintf("g%%0%dd", digits);
		for (size_t gx=0; gx < std::min(group.size(),size_t(64)); ++gx) {
			independentGroup &ig = *group[gx];
			MxRList info;
			ig.exportInternalState(info, info);
			std::string name = string_snprintf(fmt.c_str(), int(1+gx));
			dbg.add(name.c_str(), info.asR());
		}
	}

};
