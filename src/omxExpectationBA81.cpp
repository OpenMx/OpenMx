/*
  Copyright 2012-2013 Joshua Nathaniel Pritikin and contributors

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

#include <limits>
#include <Rmath.h>

#include "omxExpectationBA81.h"
#include "glue.h"
#include "libifa-rpf.h"
#include "dmvnorm.h"
#include "omxBuffer.h"

const struct rpf *rpf_model = NULL;
int rpf_numModels;

void pda(const double *ar, int rows, int cols)
{
	if (rows == 0 || cols == 0) return;
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%.6g, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

void pia(const int *ar, int rows, int cols)
{
	if (rows == 0 || cols == 0) return;
	std::string buf;
	for (int rx=0; rx < rows; rx++) {   // column major order
		for (int cx=0; cx < cols; cx++) {
			buf += string_snprintf("%d, ", ar[cx * rows + rx]);
		}
		buf += "\n";
	}
	mxLogBig(buf);
}

void ba81LikelihoodSlow2(BA81Expect *state, int px, double *out)
{
	const long totalQuadPoints = state->totalQuadPoints;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	const size_t numItems = state->itemSpec.size();
	omxData *data = state->data;
	const int *rowMap = state->rowMap;
	double *oProb = state->outcomeProb;
	std::vector<double> &priQarea = state->priQarea;

	for (long qx=0; qx < totalQuadPoints; ++qx) {
		out[qx] = priQarea[qx];
	}

	for (size_t ix=0; ix < numItems; ix++) {
		int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;

		for (long qx=0; qx < totalQuadPoints; ++qx) {
			out[qx] *= oProb[pick];
			oProb += itemOutcomes[ix];
		}
	}
}

void cai2010EiEis(BA81Expect *state, int px, double *lxk, double *Eis, double *Ei)
{
	const int numSpecific = state->numSpecific;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	double *oProb = state->outcomeProb;
	const long totalQuadPoints = state->totalQuadPoints;
	const long totalPrimaryPoints = state->totalPrimaryPoints;
	const long specificPoints = state->quadGridSize;
	const size_t numItems = state->itemSpec.size();
	const double OneOverLargest = state->OneOverLargestDouble;
	omxData *data = state->data;
	const int *rowMap = state->rowMap;
	std::vector<double> &speQarea = state->speQarea;
	std::vector<double> &priQarea = state->priQarea;

	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (long sx=0; sx < specificPoints * numSpecific; sx++) {
			lxk[qloc] = speQarea[sx];
			++qloc;
		}
	}

	for (size_t ix=0; ix < numItems; ix++) {
		int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
		if (pick == NA_INTEGER) {
			oProb += itemOutcomes[ix] * totalQuadPoints;
			continue;
		}
		pick -= 1;
		int Sgroup = state->Sgroup[ix];
		double *out1 = lxk;
		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			out1[Sgroup] *= oProb[pick];
			oProb += itemOutcomes[ix];
			out1 += numSpecific;
		}
	}

	for (long qx=0; qx < totalPrimaryPoints * numSpecific; ++qx) Eis[qx] = 0;
	for (long qx=0; qx < totalPrimaryPoints; ++qx) Ei[qx] = priQarea[qx];

	long eisloc = 0;
	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (long sx=0; sx < specificPoints; sx++) {
			for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
				double piece = lxk[qloc];
				Eis[eisloc + sgroup] += piece;
				++qloc;
			}
		}
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Ei[qx] *= Eis[eisloc + sgroup] * OneOverLargest;
		}
		eisloc += numSpecific;
	}

	for (long qx=0, qloc = 0; qx < totalPrimaryPoints; qx++) {
		for (int sgroup=0; sgroup < numSpecific; ++sgroup) {
			Eis[qloc] = Ei[qx] / Eis[qloc];
			++qloc;
		}
	}
}

OMXINLINE static void
mapLatentSpace(BA81Expect *state, int sgroup, double piece, const double *where,
	       const double *whereGram, double *latentDist)
{
	// could specialize this for regular and cai2010 for a small gain TODO
	// also, don't need to estimate latent distribution if no free parameters TODO
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;
	int pmax = maxDims;
	if (state->numSpecific) pmax -= 1;

	if (sgroup == 0) {
		int gx = 0;
		int cx = maxAbilities;
		for (int d1=0; d1 < pmax; d1++) {
			double piece_w1 = piece * where[d1];
			latentDist[d1] += piece_w1;
			for (int d2=0; d2 <= d1; d2++) {
				double piece_cov = piece * whereGram[gx];
				latentDist[cx] += piece_cov;
				++cx; ++gx;
			}
		}
	}

	if (state->numSpecific) {
		int sdim = pmax + sgroup;
		double piece_w1 = piece * where[pmax];
		latentDist[sdim] += piece_w1;

		double piece_var = piece * whereGram[triangleLoc0(pmax)];
		int to = maxAbilities + triangleLoc0(sdim);
		latentDist[to] += piece_var;
	}
}

// Depends on item parameters, but not latent distribution
void ba81OutcomeProb(BA81Expect *state, bool estep, bool wantLog)
{
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	std::vector<int> &cumItemOutcomes = state->cumItemOutcomes;
	omxMatrix *itemParam = state->itemParam;
	omxMatrix *design = state->design;
	const int maxDims = state->maxDims;
	const size_t numItems = state->itemSpec.size();
	state->outcomeProb = Realloc(state->outcomeProb, state->totalOutcomes * state->totalQuadPoints, double);
	double *param = (estep && state->EitemParam)? state->EitemParam : itemParam->data;

#pragma omp parallel for num_threads(Global->numThreads)
	for (size_t ix=0; ix < numItems; ix++) {
		double *qProb = state->outcomeProb + cumItemOutcomes[ix] * state->totalQuadPoints;
		const double *spec = state->itemSpec[ix];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		double *iparam = param + ix * itemParam->rows;
		rpf_prob_t prob_fn = wantLog? rpf_model[id].logprob : rpf_model[id].prob;

		for (long qx=0; qx < state->totalQuadPoints; qx++) {
			double *where = state->wherePrep.data() + qx * maxDims;

			double ptheta[dims];
			for (int dx=0; dx < dims; dx++) {
				int ability = (int)omxMatrixElement(design, dx, ix) - 1;
				if (ability >= maxDims) ability = maxDims-1;
				ptheta[dx] = where[ability];
			}

			(*prob_fn)(spec, iparam, ptheta, qProb);

			qProb += itemOutcomes[ix];
		}
	}
}

static void ba81Estep1(omxExpectation *oo)
{
	if(OMX_DEBUG) {mxLog("Beginning %s Computation.", oo->name);}

	BA81Expect *state = (BA81Expect*) oo->argStruct;
	const int numThreads = Global->numThreads;
	const int numUnique = state->numUnique;
	const int numSpecific = state->numSpecific;
	const int maxDims = state->maxDims;
	const int whereGramSize = triangleLoc1(maxDims);
	const int maxAbilities = state->maxAbilities;
	const int primaryDims = numSpecific? maxDims-1 : maxDims;
	omxData *data = state->data;
	int *numIdentical = state->numIdentical;
	const long totalQuadPoints = state->totalQuadPoints;

	state->excludedPatterns = 0;
	state->patternLik = Realloc(state->patternLik, numUnique, double);
	double *patternLik = state->patternLik;

	const int numLatents = maxAbilities + triangleLoc1(maxAbilities);
	std::vector<double> latentDist(numLatents, 0.0);

	const size_t numItems = state->itemSpec.size();
	const int totalOutcomes = state->totalOutcomes;
	std::vector<int> &itemOutcomes = state->itemOutcomes;
	const int *rowMap = state->rowMap;
	std::vector<double> thrExpected(totalOutcomes * totalQuadPoints * numThreads, 0.0);
	double *wherePrep = state->wherePrep.data();
	double *whereGram = state->whereGram.data();

	if (numSpecific == 0) {
		omxBuffer<double> thrLxk(totalQuadPoints * numThreads);
		omxBuffer<double> thrQweight(totalQuadPoints * numThreads);
		std::vector<double> thrDweight(totalQuadPoints * numThreads, 0.0);

#pragma omp parallel for num_threads(numThreads)
		for (int px=0; px < numUnique; px++) {
			int thrId = omx_absolute_thread_num();
			double *Qweight = thrQweight.data() + totalQuadPoints * thrId;
			double *Dweight = thrDweight.data() + totalQuadPoints * thrId;
			double *lxk = thrLxk.data() + thrId * totalQuadPoints;
			ba81LikelihoodSlow2(state, px, lxk);

			double patternLik1 = 0;
			for (long qx=0; qx < totalQuadPoints; qx++) {
				double tmp = lxk[qx];
				Qweight[qx] = tmp;
				patternLik1 += tmp;
			}

			patternLik[px] = patternLik1;

			// This uses the previous iteration's latent distribution.
			// If we recompute patternLikelihood to get the current
			// iteration's expected scores then it speeds up convergence.
			// However, recomputing patternLikelihood and dependent
			// math takes much longer than simply using the data
			// we have available here. This is even more true for the
			// two-tier model.
			if (!validPatternLik(state, patternLik1)) {
#pragma omp atomic
				state->excludedPatterns += 1;
				continue;
			}

			double weight = numIdentical[px] / patternLik1;
			for (long qx=0; qx < totalQuadPoints; ++qx) {
				double tmp = Qweight[qx] * weight;
				Dweight[qx] += tmp;
				Qweight[qx] = tmp;
			}

			double *myExpected = thrExpected.data() + thrId * totalOutcomes * totalQuadPoints;
			double *out = myExpected;
			for (size_t ix=0; ix < numItems; ++ix) {
				int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
				if (pick == NA_INTEGER) {
					out += itemOutcomes[ix] * totalQuadPoints;
					continue;
				}
				pick -= 1;

				for (long qx=0; qx < totalQuadPoints; ++qx) {
					out[pick] += Qweight[qx];
					out += itemOutcomes[ix];
				}
			}
		}

		for (int tx=1; tx < numThreads; ++tx) {
			double *Dweight = thrDweight.data() + totalQuadPoints * tx;
			double *dest = thrDweight.data();
			for (long qx=0; qx < totalQuadPoints; ++qx) {
				dest[qx] += Dweight[qx];
			}
		}

		for (long qx=0; qx < totalQuadPoints; ++qx) {
			mapLatentSpace(state, 0, thrDweight[qx], wherePrep + qx * maxDims,
				       whereGram + qx * whereGramSize, latentDist.data());
		}
	} else {
		omxBuffer<double> thrLxk(totalQuadPoints * numSpecific * numThreads);
		const long totalPrimaryPoints = state->totalPrimaryPoints;
		const long specificPoints = state->quadGridSize;
		omxBuffer<double> thrEi(totalPrimaryPoints * numThreads);
		omxBuffer<double> thrEis(totalPrimaryPoints * numSpecific * numThreads);
		omxBuffer<double> thrQweight(totalQuadPoints * numSpecific * numThreads);
		std::vector<double> thrDweight(totalQuadPoints * numSpecific * numThreads);

#pragma omp parallel for num_threads(numThreads)
		for (int px=0; px < numUnique; px++) {
			int thrId = omx_absolute_thread_num();
			double *Qweight = thrQweight.data() + totalQuadPoints * numSpecific * thrId;
			double *Dweight = thrDweight.data() + totalQuadPoints * numSpecific * thrId;

			double *lxk = thrLxk.data() + totalQuadPoints * numSpecific * thrId;
			double *Ei = thrEi.data() + totalPrimaryPoints * thrId;
			double *Eis = thrEis.data() + totalPrimaryPoints * numSpecific * thrId;
			cai2010EiEis(state, px, lxk, Eis, Ei);

			for (long qloc=0, eisloc=0; eisloc < totalPrimaryPoints * numSpecific; eisloc += numSpecific) {
				for (long sx=0; sx < specificPoints; sx++) {
					for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
						double lxk1 = lxk[qloc];
						double Eis1 = Eis[eisloc + Sgroup];
						double tmp = Eis1 * lxk1;
						Qweight[qloc] = tmp;
						++qloc;
					}
				}
			}

			double patternLik1 = 0;
			for (long qx=0; qx < totalPrimaryPoints; ++qx) {
				patternLik1 += Ei[qx];
			}
			patternLik[px] = patternLik1;

			if (!validPatternLik(state, patternLik1)) {
#pragma omp atomic
				state->excludedPatterns += 1;
				continue;
			}

			double *myExpected = thrExpected.data() + thrId * totalOutcomes * totalQuadPoints;
			double weight = numIdentical[px] / patternLik1;
			for (long qx=0; qx < totalQuadPoints * numSpecific; qx++) {
				double tmp = Qweight[qx] * weight;
				Qweight[qx] = tmp;
				Dweight[qx] += tmp;
			}

			double *out = myExpected;
			for (size_t ix=0; ix < numItems; ++ix) {
				int pick = omxIntDataElementUnsafe(data, rowMap[px], ix);
				if (pick == NA_INTEGER) {
					out += itemOutcomes[ix] * totalQuadPoints;
					continue;
				}
				pick -= 1;

				int Sgroup = state->Sgroup[ix];
				double *Qw = Qweight;
				for (long qx=0; qx < totalQuadPoints; ++qx) {
					out[pick] += Qw[Sgroup];
					out += itemOutcomes[ix];
					Qw += numSpecific;
				}
			}
		}

		for (int tx=1; tx < numThreads; ++tx) {
			double *Dweight = thrDweight.data() + totalQuadPoints * numSpecific * tx;
			double *dest = thrDweight.data();
			for (long qx=0; qx < totalQuadPoints * numSpecific; ++qx) {
				dest[qx] += Dweight[qx];
			}
		}

		long qloc=0;
		for (long qx=0; qx < totalQuadPoints; qx++) {
			double *whPrep = wherePrep + qx * maxDims;
			double *whGram = whereGram + qx * whereGramSize;
			for (int Sgroup=0; Sgroup < numSpecific; Sgroup++) {
				mapLatentSpace(state, Sgroup, thrDweight[qloc], whPrep, whGram, latentDist.data());
				++qloc;
			}
		}
	}

	const long expectedSize = totalQuadPoints * totalOutcomes;
	OMXZERO(state->expected, expectedSize);

	double *e1 = thrExpected.data();
	for (int tx=0; tx < numThreads; ++tx) {
		for (long ex=0; ex < expectedSize; ++ex) {
			state->expected[ex] += *e1;
			++e1;
		}
	}

	//mxLog("raw latent after weighting");
	//pda(latentDist, numLatents, numUnique);

	std::vector<double> &ElatentMean = state->ElatentMean;
	std::vector<double> &ElatentCov = state->ElatentCov;
	
	{
		double *latentDist1 = latentDist.data();
		for (int d1=0; d1 < maxAbilities; d1++) {
			int cx = maxAbilities + triangleLoc1(d1);
			if (d1 < primaryDims) {
				ElatentMean[d1] = latentDist1[d1];
				for (int d2=0; d2 <= d1; d2++) {
					int cell = d2 * maxAbilities + d1;
					ElatentCov[cell] = latentDist1[cx];
					++cx;
				}
			} else {
				ElatentMean[d1] = latentDist1[d1];
				int cell = d1 * maxAbilities + d1;
				int loc = maxAbilities + triangleLoc0(d1);
				ElatentCov[cell] = latentDist1[loc];
			}
		}
	}

	//pda(ElatentMean.data(), 1, state->maxAbilities);
	//pda(ElatentCov.data(), state->maxAbilities, state->maxAbilities);

	for (int d1=0; d1 < maxAbilities; d1++) {
		ElatentMean[d1] /= data->rows;
	}

	for (int d1=0; d1 < primaryDims; d1++) {
		for (int d2=0; d2 <= d1; d2++) {
			int cell = d2 * maxAbilities + d1;
			int tcell = d1 * maxAbilities + d2;
			ElatentCov[tcell] = ElatentCov[cell] =
				ElatentCov[cell] / data->rows - ElatentMean[d1] * ElatentMean[d2];
		}
	}
	for (int sdim=primaryDims; sdim < maxAbilities; sdim++) {
		int cell = sdim * maxAbilities + sdim;
		ElatentCov[cell] = ElatentCov[cell] / data->rows - ElatentMean[sdim] * ElatentMean[sdim];
	}

	// (mirt) need to check for negative eigenvalues in the cov matrix when more than 1 factor? TODO

	if (state->verbose) {
		mxLog("%s: lxk(%d) patternLik (%d/%d excluded) ElatentMean ElatentCov",
		      oo->name, omxGetMatrixVersion(state->itemParam),
		      state->excludedPatterns, numUnique);
		//pda(ElatentMean.data(), 1, state->maxAbilities);
		//pda(ElatentCov.data(), state->maxAbilities, state->maxAbilities);
	}

	++state->ElatentVersion;
}

static int getLatentVersion(BA81Expect *state)
{
	return omxGetMatrixVersion(state->latentMeanOut) + omxGetMatrixVersion(state->latentCovOut);
}

// Attempt G-H grid? http://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
void ba81SetupQuadrature(omxExpectation* oo)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	if (state->Qpoint.size() == 0 && latentClean) return;

	if (state->verbose) {
		mxLog("%s: quadrature(%d)", oo->name, getLatentVersion(state));
	}

	int gridsize = state->targetQpoints;
	const int maxDims = state->maxDims;
	double Qwidth = state->Qwidth;
	int numSpecific = state->numSpecific;
	int priDims = maxDims - (numSpecific? 1 : 0);

	state->totalQuadPoints = 1;
	for (int dx=0; dx < maxDims; dx++) {
		state->totalQuadPoints *= gridsize;
	}

	state->Qpoint.resize(gridsize);
	double qgs = gridsize-1;
	for (int px=0; px < gridsize; ++px) {
		state->Qpoint[px] = Qwidth - px * 2 * Qwidth / qgs;
	}

	if (state->quadGridSize != gridsize) {
		const long totalQuadPoints = state->totalQuadPoints;
		std::vector<double> &wherePrep = state->wherePrep;
		wherePrep.resize(totalQuadPoints * maxDims);
		std::vector<double> &whereGram = state->whereGram;
		whereGram.resize(totalQuadPoints * triangleLoc1(maxDims));
		
		for (long qx=0; qx < totalQuadPoints; qx++) {
			double *wh = wherePrep.data() + qx * maxDims;
			int quad[maxDims];
			decodeLocation(qx, maxDims, gridsize, quad);
			pointToWhere(state, quad, wh, maxDims);
			gramProduct(wh, maxDims, whereGram.data() + qx * triangleLoc1(maxDims));
		}

		// try starting small and increasing to the cap? TODO
		state->quadGridSize = gridsize;
	}

	state->totalPrimaryPoints = state->totalQuadPoints;

	if (numSpecific) {
		state->totalPrimaryPoints /= state->quadGridSize;
		state->speQarea.resize(gridsize * numSpecific);
	}

	state->priQarea.resize(state->totalPrimaryPoints);

	//pda(state->latentMeanOut->data, 1, state->maxAbilities);
	//pda(state->latentCovOut->data, state->maxAbilities, state->maxAbilities);

	omxBuffer<double> priCovData(priDims * priDims);
	for (int d1=0; d1 < priDims; ++d1) {
		for (int d2=0; d2 < priDims; ++d2) {
			priCovData[d1 * priDims + d2] = omxMatrixElement(state->latentCovOut, d1, d2);
		}
	}

	const double Largest = state->LargestDouble;
	double totalArea = 0;
	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		int quad[priDims];
		decodeLocation(qx, priDims, state->quadGridSize, quad);
		double where[priDims];
		pointToWhere(state, quad, where, priDims);
		state->priQarea[qx] = exp(dmvnorm(priDims, where,
						  state->latentMeanOut->data,
						  priCovData.data()));
		totalArea += state->priQarea[qx];
	}
	for (int qx=0; qx < state->totalPrimaryPoints; qx++) {
		state->priQarea[qx] *= Largest;
		state->priQarea[qx] /= totalArea;
		//mxLog("%.5g,", state->priQarea[qx]);
	}

	for (int sgroup=0; sgroup < numSpecific; sgroup++) {
		totalArea = 0;
		int covCell = (priDims + sgroup) * state->maxAbilities + priDims + sgroup;
		double mean = state->latentMeanOut->data[priDims + sgroup];
		double var = state->latentCovOut->data[covCell];
		//mxLog("setup[%d] %.2f %.2f", sx, mean, var);
		for (int qx=0; qx < state->quadGridSize; qx++) {
			double den = dnorm(state->Qpoint[qx], mean, sqrt(var), FALSE);
			state->speQarea[sIndex(state, sgroup, qx)] = den;
			totalArea += den;
		}
		for (int qx=0; qx < state->quadGridSize; qx++) {
			state->speQarea[sIndex(state, sgroup, qx)] *= Largest;
			state->speQarea[sIndex(state, sgroup, qx)] /= totalArea;
		}
		//pda(state->speQarea.data() + sIndex(state, sgroup, 0), 1, state->quadGridSize);
	}

	// The idea here is to avoid denormalized values if they are
	// enabled (5e-324 vs 2e-308).  It would be bad if results
	// changed depending on the denormalization setting.
	// Moreover, we don't lose too much even if denormalized
	// values are disabled.

	state->SmallestPatternLik = 1e16 * std::numeric_limits<double>::min();

	state->expected = Realloc(state->expected, state->totalOutcomes * state->totalQuadPoints, double);
	state->latentParamVersion = getLatentVersion(state);
}

OMXINLINE static void
accumulateScores(BA81Expect *state, int px, int sgroup, double piece, const double *where,
		 int primaryDims, int covEntries, std::vector<double> *mean, std::vector<double> *cov)
{
	int maxDims = state->maxDims;
	int maxAbilities = state->maxAbilities;

	if (sgroup == 0) {
		int cx=0;
		for (int d1=0; d1 < primaryDims; d1++) {
			double piece_w1 = piece * where[d1];
			double &dest1 = (*mean)[px * maxAbilities + d1];
			dest1 += piece_w1;
			for (int d2=0; d2 <= d1; d2++) {
				double &dest2 = (*cov)[px * covEntries + cx];
				dest2 += where[d2] * piece_w1;
				++cx;
			}
		}
	}

	if (state->numSpecific) {
		int sdim = maxDims + sgroup - 1;
		double piece_w1 = piece * where[primaryDims];
		double &dest3 = (*mean)[px * maxAbilities + sdim];
		dest3 += piece_w1;

		double &dest4 = (*cov)[px * covEntries + triangleLoc0(sdim)];
		dest4 += piece_w1 * where[primaryDims];
	}
}

// re-factor to share code with E-step TODO
static void
EAPinternalFast(omxExpectation *oo, std::vector<double> *mean, std::vector<double> *cov)
{
	BA81Expect *state = (BA81Expect*) oo->argStruct;
	if (state->verbose) mxLog("%s: EAP", oo->name);

	const int numUnique = state->numUnique;
	const int numSpecific = state->numSpecific;
	const int maxDims = state->maxDims;
	const int maxAbilities = state->maxAbilities;
	const int primaryDims = numSpecific? maxDims-1 : maxDims;
	const int covEntries = triangleLoc1(maxAbilities);
	double *patternLik = state->patternLik;
	const long totalQuadPoints = state->totalQuadPoints;
	const long totalPrimaryPoints = state->totalPrimaryPoints;

	mean->assign(numUnique * maxAbilities, 0);
	cov->assign(numUnique * covEntries, 0);

	if (numSpecific == 0) {
		omxBuffer<double> thrLxk(totalQuadPoints * Global->numThreads);

#pragma omp parallel for num_threads(Global->numThreads)
		for (int px=0; px < numUnique; px++) {
			if (!validPatternLik(state, patternLik[px])) {
				continue;
			}

			int thrId = omx_absolute_thread_num();
			double *lxk = thrLxk.data() + thrId * totalQuadPoints;
			ba81LikelihoodSlow2(state, px, lxk);

			for (long qx=0; qx < state->totalQuadPoints; qx++) {
				int quad[maxDims];
				decodeLocation(qx, maxDims, state->quadGridSize, quad);
				double where[maxDims];
				pointToWhere(state, quad, where, maxDims);

				double tmp = lxk[qx];
				accumulateScores(state, px, 0, tmp, where, primaryDims, covEntries, mean, cov);
			}
		}
	} else {
		int sDim = primaryDims;
		const long specificPoints = state->quadGridSize;
		omxBuffer<double> thrLxk(totalQuadPoints * numSpecific * Global->numThreads);
		omxBuffer<double> thrEi(totalPrimaryPoints * Global->numThreads);
		omxBuffer<double> thrEis(totalPrimaryPoints * numSpecific * Global->numThreads);

#pragma omp parallel for num_threads(Global->numThreads)
		for (int px=0; px < numUnique; px++) {
			if (!validPatternLik(state, patternLik[px])) {
				continue;
			}

			int thrId = omx_absolute_thread_num();
			double *lxk = thrLxk.data() + totalQuadPoints * numSpecific * thrId;
			double *Eis = thrEis.data() + totalPrimaryPoints * numSpecific * thrId;
			double *Ei = thrEi.data() + totalPrimaryPoints * thrId;
			cai2010EiEis(state, px, lxk, Eis, Ei);

			long qloc = 0;
			long eisloc = 0;
			for (long qx=0; qx < totalPrimaryPoints; qx++) {
				int quad[maxDims];
				decodeLocation(qx, primaryDims, state->quadGridSize, quad);
				for (long sx=0; sx < specificPoints; sx++) {
					for (int Sgroup=0; Sgroup < numSpecific; ++Sgroup) {
						quad[sDim] = sx;
						double where[maxDims];
						pointToWhere(state, quad, where, maxDims);
						double lxk1 = lxk[qloc];
						double Eis1 = Eis[eisloc + Sgroup];
						double tmp = Eis1 * lxk1;
						accumulateScores(state, px, Sgroup, tmp, where, primaryDims,
								 covEntries, mean, cov);
						++qloc;
					}
				}
				eisloc += numSpecific;
			}
		}
	}

	for (int px=0; px < numUnique; px++) {
		double denom = patternLik[px];
		if (!validPatternLik(state, denom)) {
			for (int ax=0; ax < maxAbilities; ++ax) {
				(*mean)[px * maxAbilities + ax] = NA_REAL;
			}
			for (int cx=0; cx < covEntries; ++cx) {
				(*cov)[px * covEntries + cx] = NA_REAL;
			}
			continue;
		}
		for (int ax=0; ax < maxAbilities; ax++) {
			(*mean)[px * maxAbilities + ax] /= denom;
		}
		for (int cx=0; cx < triangleLoc1(primaryDims); ++cx) {
			(*cov)[px * covEntries + cx] /= denom;
		}
		for (int sx=0; sx < numSpecific; sx++) {
			(*cov)[px * covEntries + triangleLoc0(primaryDims + sx)] /= denom;
		}
		int cx=0;
		for (int a1=0; a1 < primaryDims; ++a1) {
			for (int a2=0; a2 <= a1; ++a2) {
				double ma1 = (*mean)[px * maxAbilities + a1];
				double ma2 = (*mean)[px * maxAbilities + a2];
				(*cov)[px * covEntries + cx] -= ma1 * ma2;
				++cx;
			}
		}
		for (int sx=0; sx < numSpecific; sx++) {
			int sdim = primaryDims + sx;
			double ma1 = (*mean)[px * maxAbilities + sdim];
			(*cov)[px * covEntries + triangleLoc0(sdim)] -= ma1 * ma1;
		}
        }
}

static void
ba81compute(omxExpectation *oo, const char *context)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;

	if (context) {
		if (strcmp(context, "scores")==0) {
			state->type = EXPECTATION_AUGMENTED;
		} else if (strcmp(context, "nothing")==0) {
			state->type = EXPECTATION_OBSERVED;
		} else {
			omxRaiseErrorf(globalState, "%s: don't know how to predict '%s'",
				       oo->name, context);
		}
		return;
	}

	bool latentClean = state->latentParamVersion == getLatentVersion(state);
	bool itemClean = state->itemParamVersion == omxGetMatrixVersion(state->itemParam) && latentClean;

	if (state->verbose) {
		mxLog("%s: Qinit %d itemClean %d latentClean %d (1=clean)",
		      oo->name, state->Qpoint.size() != 0, itemClean, latentClean);
	}

	ba81SetupQuadrature(oo);

	if (!itemClean) {
		ba81OutcomeProb(state, TRUE, FALSE);
		ba81Estep1(oo);
	}

	state->itemParamVersion = omxGetMatrixVersion(state->itemParam);
}

static void
copyScore(int rows, int maxAbilities, std::vector<double> &mean,
	  std::vector<double> &cov, const int rx, double *scores, const int dest)
{
	for (int ax=0; ax < maxAbilities; ++ax) {
		scores[rows * ax + dest] = mean[maxAbilities * rx + ax];
	}
	for (int ax=0; ax < maxAbilities; ++ax) {
		scores[rows * (maxAbilities + ax) + dest] =
			sqrt(cov[triangleLoc1(maxAbilities) * rx + triangleLoc0(ax)]);
	}
	for (int ax=0; ax < triangleLoc1(maxAbilities); ++ax) {
		scores[rows * (2*maxAbilities + ax) + dest] =
			cov[triangleLoc1(maxAbilities) * rx + ax];
	}
}

/**
 * MAP is not affected by the number of items. EAP is. Likelihood can
 * get concentrated in a single quadrature ordinate. For 3PL, response
 * patterns can have a bimodal likelihood. This will confuse MAP and
 * is a key advantage of EAP (Thissen & Orlando, 2001, p. 136).
 *
 * Thissen, D. & Orlando, M. (2001). IRT for items scored in two
 * categories. In D. Thissen & H. Wainer (Eds.), \emph{Test scoring}
 * (pp 73-140). Lawrence Erlbaum Associates, Inc.
 */
static void
ba81PopulateAttributes(omxExpectation *oo, SEXP robj)
{
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	int maxAbilities = state->maxAbilities;

	setAttrib(robj, install("numStats"), ScalarReal(state->numUnique - 1)); // missingness? latent params? TODO

	if (state->debugInternal) {
		const double LogLargest = state->LogLargestDouble;
		int numUnique = state->numUnique;
		int totalOutcomes = state->totalOutcomes;
		SEXP Rlik;
		SEXP Rexpected;

		PROTECT(Rlik = allocVector(REALSXP, numUnique));
		memcpy(REAL(Rlik), state->patternLik, sizeof(double) * numUnique);
		double *lik_out = REAL(Rlik);
		for (int px=0; px < numUnique; ++px) {
			// Must return value in log units because it may not be representable otherwise
			lik_out[px] = log(lik_out[px]) - LogLargest;
		}

		PROTECT(Rexpected = allocVector(REALSXP, state->totalQuadPoints * totalOutcomes));
		memcpy(REAL(Rexpected), state->expected, sizeof(double) * totalOutcomes * state->totalQuadPoints);

		MxRList dbg;
		dbg.push_back(std::make_pair(mkChar("patternLikelihood"), Rlik));
		dbg.push_back(std::make_pair(mkChar("em.expected"), Rexpected));

		SEXP Rmean, Rcov;
		PROTECT(Rmean = allocVector(REALSXP, maxAbilities));
		memcpy(REAL(Rmean), state->ElatentMean.data(), maxAbilities * sizeof(double));

		PROTECT(Rcov = allocMatrix(REALSXP, maxAbilities, maxAbilities));
		memcpy(REAL(Rcov), state->ElatentCov.data(), maxAbilities * maxAbilities * sizeof(double));

		dbg.push_back(std::make_pair(mkChar("mean"), Rmean));
		dbg.push_back(std::make_pair(mkChar("cov"), Rcov));

		setAttrib(robj, install("debug"), dbg.asR());
	}

	if (state->scores == SCORES_OMIT || state->type == EXPECTATION_UNINITIALIZED) return;

	// TODO Wainer & Thissen. (1987). Estimating ability with the wrong
	// model. Journal of Educational Statistics, 12, 339-368.

	/*
	int numQpoints = state->targetQpoints * 2;  // make configurable TODO

	if (numQpoints < 1 + 2.0 * sqrt(state->itemSpec->cols)) {
		// Thissen & Orlando (2001, p. 136)
		warning("EAP requires at least 2*sqrt(items) quadrature points");
	}

	ba81SetupQuadrature(oo, numQpoints, 0);
	ba81Estep1(oo);
	*/

	std::vector<double> mean;
	std::vector<double> cov;
	EAPinternalFast(oo, &mean, &cov);

	int numUnique = state->numUnique;
	omxData *data = state->data;
	int rows = state->scores == SCORES_FULL? data->rows : numUnique;
	int cols = 2 * maxAbilities + triangleLoc1(maxAbilities);
	SEXP Rscores;
	PROTECT(Rscores = allocMatrix(REALSXP, rows, cols));
	double *scores = REAL(Rscores);

	const int SMALLBUF = 10;
	char buf[SMALLBUF];
	SEXP names;
	PROTECT(names = allocVector(STRSXP, cols));
	for (int nx=0; nx < maxAbilities; ++nx) {
		snprintf(buf, SMALLBUF, "s%d", nx+1);
		SET_STRING_ELT(names, nx, mkChar(buf));
		snprintf(buf, SMALLBUF, "se%d", nx+1);
		SET_STRING_ELT(names, maxAbilities + nx, mkChar(buf));
	}
	for (int nx=0; nx < triangleLoc1(maxAbilities); ++nx) {
		snprintf(buf, SMALLBUF, "cov%d", nx+1);
		SET_STRING_ELT(names, maxAbilities*2 + nx, mkChar(buf));
	}
	SEXP dimnames;
	PROTECT(dimnames = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dimnames, 1, names);
	setAttrib(Rscores, R_DimNamesSymbol, dimnames);

	if (state->scores == SCORES_FULL) {
#pragma omp parallel for num_threads(Global->numThreads)
		for (int rx=0; rx < numUnique; rx++) {
			int dups = omxDataNumIdenticalRows(state->data, state->rowMap[rx]);
			for (int dup=0; dup < dups; dup++) {
				int dest = omxDataIndex(data, state->rowMap[rx]+dup);
				copyScore(rows, maxAbilities, mean, cov, rx, scores, dest);
			}
		}
	} else {
#pragma omp parallel for num_threads(Global->numThreads)
		for (int rx=0; rx < numUnique; rx++) {
			copyScore(rows, maxAbilities, mean, cov, rx, scores, rx);
		}
	}

	MxRList out;
	out.push_back(std::make_pair(mkChar("scores"), Rscores));
	setAttrib(robj, install("output"), out.asR());
}

static void ba81Destroy(omxExpectation *oo) {
	if(OMX_DEBUG) {
		mxLog("Freeing %s function.", oo->name);
	}
	BA81Expect *state = (BA81Expect *) oo->argStruct;
	omxFreeAllMatrixData(state->design);
	omxFreeAllMatrixData(state->latentMeanOut);
	omxFreeAllMatrixData(state->latentCovOut);
	omxFreeAllMatrixData(state->itemParam);
	Free(state->numIdentical);
	Free(state->rowMap);
	Free(state->patternLik);
	Free(state->Sgroup);
	Free(state->expected);
	Free(state->outcomeProb);
	delete state;
}

void getMatrixDims(SEXP r_theta, int *rows, int *cols)
{
    SEXP matrixDims;
    PROTECT(matrixDims = getAttrib(r_theta, R_DimSymbol));
    int *dimList = INTEGER(matrixDims);
    *rows = dimList[0];
    *cols = dimList[1];
    UNPROTECT(1);
}

static void ignoreSetVarGroup(omxExpectation*, FreeVarGroup *)
{}

void omxInitExpectationBA81(omxExpectation* oo) {
	omxState* currentState = oo->currentState;	
	SEXP rObj = oo->rObj;
	SEXP tmp;
	
	if(OMX_DEBUG) {
		mxLog("Initializing %s.", oo->name);
	}
	if (!rpf_model) {
		if (0) {
			const int wantVersion = 3;
			int version;
			get_librpf_t get_librpf = (get_librpf_t) R_GetCCallable("rpf", "get_librpf_model_GPL");
			(*get_librpf)(&version, &rpf_numModels, &rpf_model);
			if (version < wantVersion) error("librpf binary API %d installed, at least %d is required",
							 version, wantVersion);
		} else {
			rpf_numModels = librpf_numModels;
			rpf_model = librpf_model;
		}
	}
	
	BA81Expect *state = new BA81Expect;

	// These two constants should be as identical as possible
	state->LogLargestDouble = log(std::numeric_limits<double>::max()) - 1;
	state->LargestDouble = exp(state->LogLargestDouble);
	state->OneOverLargestDouble = 1/state->LargestDouble;

	state->numSpecific = 0;
	state->excludedPatterns = 0;
	state->numIdentical = NULL;
	state->rowMap = NULL;
	state->design = NULL;
	state->patternLik = NULL;
	state->outcomeProb = NULL;
	state->expected = NULL;
	state->type = EXPECTATION_UNINITIALIZED;
	state->scores = SCORES_OMIT;
	state->itemParam = NULL;
	state->EitemParam = NULL;
	state->itemParamVersion = 0;
	state->latentParamVersion = 0;
	state->quadGridSize = 0;
	oo->argStruct = (void*) state;

	PROTECT(tmp = GET_SLOT(rObj, install("data")));
	state->data = omxDataLookupFromState(tmp, currentState);

	if (strcmp(omxDataType(state->data), "raw") != 0) {
		omxRaiseErrorf(currentState, "%s unable to handle data type %s", oo->name, omxDataType(state->data));
		return;
	}

	PROTECT(tmp = GET_SLOT(rObj, install("ItemSpec")));
	for (int sx=0; sx < length(tmp); ++sx) {
		SEXP model = VECTOR_ELT(tmp, sx);
		if (!OBJECT(model)) {
			error("Item models must inherit rpf.base");
		}
		SEXP spec;
		PROTECT(spec = GET_SLOT(model, install("spec")));
		state->itemSpec.push_back(REAL(spec));
	}

	PROTECT(tmp = GET_SLOT(rObj, install("design")));
	if (!isNull(tmp)) {
		// better to demand integers and not coerce to real TODO
		state->design = omxNewMatrixFromRPrimitive(tmp, globalState, FALSE, 0);
	}

	state->latentMeanOut = omxNewMatrixFromSlot(rObj, currentState, "mean");
	if (!state->latentMeanOut) error("Failed to retrieve mean matrix");
	state->latentCovOut  = omxNewMatrixFromSlot(rObj, currentState, "cov");
	if (!state->latentCovOut) error("Failed to retrieve cov matrix");

	state->itemParam =
		omxNewMatrixFromSlot(rObj, globalState, "ItemParam");

	PROTECT(tmp = GET_SLOT(rObj, install("EItemParam")));
	if (!isNull(tmp)) {
		int rows, cols;
		getMatrixDims(tmp, &rows, &cols);
		if (rows != state->itemParam->rows || cols != state->itemParam->cols) {
			error("EItemParam must have same dimensions as ItemParam");
		}
		state->EitemParam = REAL(tmp);
	}

	oo->computeFun = ba81compute;
	oo->setVarGroup = ignoreSetVarGroup;
	oo->destructFun = ba81Destroy;
	oo->populateAttrFun = ba81PopulateAttributes;
	
	// TODO: Exactly identical rows do not contribute any information.
	// The sorting algorithm ought to remove them so we don't waste RAM.
	// The following summary stats would be cheaper to calculate too.

	int numUnique = 0;
	omxData *data = state->data;
	if (omxDataNumFactor(data) != data->cols) {
		// verify they are ordered factors TODO
		omxRaiseErrorf(currentState, "%s: all columns must be factors", oo->name);
		omxPrintData(data, "data", 5);
		return;
	}

	for (int rx=0; rx < data->rows;) {
		rx += omxDataNumIdenticalRows(state->data, rx);
		++numUnique;
	}
	state->numUnique = numUnique;

	state->rowMap = Realloc(NULL, numUnique, int);
	state->numIdentical = Realloc(NULL, numUnique, int);

	const int numItems = state->itemParam->cols;
	if (data->cols != numItems) {
		error("Data has %d columns for %d items", data->cols, numItems);
	}

	int maxSpec = 0;
	int maxParam = 0;
	int maxItemDims = 0;

	std::vector<int> &itemOutcomes = state->itemOutcomes;
	std::vector<int> &cumItemOutcomes = state->cumItemOutcomes;
	itemOutcomes.resize(numItems);
	cumItemOutcomes.resize(numItems);
	int totalOutcomes = 0;
	for (int cx = 0; cx < data->cols; cx++) {
		const double *spec = state->itemSpec[cx];
		int id = spec[RPF_ISpecID];
		int dims = spec[RPF_ISpecDims];
		if (maxItemDims < dims)
			maxItemDims = dims;

		int no = spec[RPF_ISpecOutcomes];
		itemOutcomes[cx] = no;
		cumItemOutcomes[cx] = totalOutcomes;
		totalOutcomes += no;

		// TODO this summary stat should be available from omxData
		int dataMax=0;
		for (int rx=0; rx < data->rows; rx++) {
			int pick = omxIntDataElementUnsafe(data, rx, cx);
			if (dataMax < pick)
				dataMax = pick;
		}
		if (dataMax > no) {
			error("Data for item %d has %d outcomes, not %d", cx+1, dataMax, no);
		}

		int numSpec = (*rpf_model[id].numSpec)(spec);
		if (maxSpec < numSpec)
			maxSpec = numSpec;

		int numParam = (*rpf_model[id].numParam)(spec);
		if (maxParam < numParam)
			maxParam = numParam;
	}

	state->totalOutcomes = totalOutcomes;

	if (int(state->itemSpec.size()) != data->cols) {
		omxRaiseErrorf(currentState, "ItemSpec must contain %d item model specifications",
			       data->cols);
		return;
	}

	if (state->design == NULL) {
		state->maxDims = maxItemDims;
		state->maxAbilities = maxItemDims;
		state->design = omxInitTemporaryMatrix(NULL, state->maxDims, numItems,
				       TRUE, currentState);
		for (int ix=0; ix < numItems; ix++) {
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			for (int dx=0; dx < state->maxDims; dx++) {
				omxSetMatrixElement(state->design, dx, ix, dx < dims? (double)dx+1 : nan(""));
			}
		}
	} else {
		omxMatrix *design = state->design;
		if (design->cols != numItems) {
			omxRaiseErrorf(currentState, "Design matrix should have %d columns", numItems);
			return;
		}

		state->maxAbilities = 0;
		for (int ix=0; ix < design->rows * design->cols; ix++) {
			double got = design->data[ix];
			if (!R_FINITE(got)) continue;
			if (round(got) != (int)got) error("Design matrix can only contain integers"); // TODO better way?
			if (state->maxAbilities < got)
				state->maxAbilities = got;
		}
		maxItemDims = 0;
		for (int ix=0; ix < design->cols; ix++) {
			const double *idesign = omxMatrixColumn(design, ix);
			int ddim = 0;
			for (int rx=0; rx < design->rows; rx++) {
				if (isfinite(idesign[rx])) ddim += 1;
			}
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			if (ddim > dims) error("Item %d has %d dims but design assigns %d dims", ix, dims, ddim);
			if (maxItemDims < ddim) {
				maxItemDims = ddim;
			}
		}
		state->maxDims = maxItemDims;
	}
	if (state->maxAbilities <= state->maxDims) {
		state->Sgroup = Calloc(numItems, int);
	} else {
		// Not sure if this is correct, revisit TODO
		int Sgroup0 = -1;
		state->Sgroup = Realloc(NULL, numItems, int);
		for (int dx=0; dx < state->maxDims; dx++) {
			for (int ix=0; ix < numItems; ix++) {
				int ability = omxMatrixElement(state->design, dx, ix);
				if (dx < state->maxDims - 1) {
					if (Sgroup0 <= ability)
						Sgroup0 = ability+1;
					continue;
				}
				int ss=-1;
				if (ability >= Sgroup0) {
					if (ss == -1) {
						ss = ability;
					} else {
						omxRaiseErrorf(currentState, "Item %d cannot belong to more than "
							       "1 specific dimension (both %d and %d)",
							       ix, ss, ability);
						return;
					}
				}
				if (ss == -1) ss = Sgroup0;
				state->Sgroup[ix] = ss - Sgroup0;
			}
		}
		state->numSpecific = state->maxAbilities - state->maxDims + 1;
	}

	// Rows with no information about an ability will obtain the
	// prior distribution as an ability estimate. This will
	// throw off multigroup latent distribution estimates.
	for (int rx=0, ux=0; rx < data->rows; ux++) {
		int dups = omxDataNumIdenticalRows(state->data, rx);
		state->numIdentical[ux] = dups;
		state->rowMap[ux] = rx;

		std::vector<bool> hasScore(state->maxAbilities);
		for (int ix=0; ix < numItems; ix++) {
			int pick = omxIntDataElementUnsafe(data, rx, ix);
			if (pick == NA_INTEGER) continue;
			const double *spec = state->itemSpec[ix];
			int dims = spec[RPF_ISpecDims];
			int dr = 0;
			for (int dx=0; dx < dims; dx++) {
				int ability = (int)omxMatrixElement(state->design, dr + dx, ix);
				while (ability == NA_INTEGER) {
					++dr;
					ability = (int)omxMatrixElement(state->design, dr + dx, ix);
				}
				// assume factor loadings are the first item parameters
				if (omxMatrixElement(state->itemParam, dx, ix) == 0) continue;
				hasScore[ability - 1] = true;
			}
		}
		for (int ax=0; ax < state->maxAbilities; ++ax) {
			if (!hasScore[ax]) {
				int dest = omxDataIndex(data, ux);
				omxRaiseErrorf(currentState,
					       "Data row %d has no information about ability %d", 1+dest, 1+ax);
				return;
			}
		}
		rx += dups;
	}

	if (state->latentMeanOut->rows * state->latentMeanOut->cols != state->maxAbilities) {
		error("The mean matrix '%s' must be 1x%d or %dx1", state->latentMeanOut->name,
		      state->maxAbilities, state->maxAbilities);
	}
	if (state->latentCovOut->rows != state->maxAbilities ||
	    state->latentCovOut->cols != state->maxAbilities) {
		error("The cov matrix '%s' must be %dx%d",
		      state->latentCovOut->name, state->maxAbilities, state->maxAbilities);
	}

	PROTECT(tmp = GET_SLOT(rObj, install("verbose")));
	state->verbose = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("debugInternal")));
	state->debugInternal = asLogical(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qpoints")));
	state->targetQpoints = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("qwidth")));
	state->Qwidth = asReal(tmp);

	PROTECT(tmp = GET_SLOT(rObj, install("scores")));
	const char *score_option = CHAR(asChar(tmp));
	if (strcmp(score_option, "omit")==0) state->scores = SCORES_OMIT;
	if (strcmp(score_option, "unique")==0) state->scores = SCORES_UNIQUE;
	if (strcmp(score_option, "full")==0) state->scores = SCORES_FULL;

	state->ElatentVersion = 0;
	state->ElatentMean.resize(state->maxAbilities);
	state->ElatentCov.resize(state->maxAbilities * state->maxAbilities);

	// verify data bounded between 1 and numOutcomes TODO
	// hm, looks like something could be added to omxData for column summary stats?
}
