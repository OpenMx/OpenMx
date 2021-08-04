/*
 *  Copyright 2021 by the individuals mentioned in the source code history
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

#include "omxDefines.h"
#include "Compute.h"
#include "penalty.h"
#include "EnableWarnings.h"

Penalty::Penalty(S4 _obj, omxMatrix *mat) : matrix(mat)
{
  robj = _obj;
  params = robj.slot("params");
  epsilon = robj.slot("epsilon");
  scale = robj.slot("scale");
  smoothProportion = as<double>(robj.slot("smoothProportion"));
}

Penalty::~Penalty() {}

const char *Penalty::name() const
{ return matrix->name(); }

int Penalty::countNumZero(FitContext *fc) const
{
  int count = 0;
  for (int px = 0; px < params.size(); ++px) {
    if (fabs(fc->est[params[px]] / scale[px % scale.size()]) <=
        epsilon[px % epsilon.size()]) ++count;
  }
  return count;
}

double Penalty::penaltyStrength(double absPar, int px) const
{
  double active = epsilon[px % epsilon.size()];
  double width = active * smoothProportion;
  double inactive = active - width;
  if (absPar > active) return 1;
  if (absPar < inactive) return 0;
  return (absPar - inactive) / width;
}

void Penalty::copyFrom(const Penalty *pen)
{
  params = pen->params;
  epsilon = pen->epsilon;
  scale = pen->scale;
  smoothProportion = pen->smoothProportion;
}

double Penalty::getValue() const { return matrix->data[0]; }

LassoPenalty::LassoPenalty(S4 obj, omxMatrix *mat) : Penalty(obj, mat)
{
  IntegerVector pv = obj.slot("hyperparameters");
  if (pv.size() != 1) mxThrow("LassoPenalty should have 1 hyperparameter (not %d)", pv.size());
  lambda = pv[0];
}

std::unique_ptr<Penalty> LassoPenalty::clone(omxMatrix *mat) const
{
  auto pen = std::make_unique<LassoPenalty>(robj, mat);
  pen->copyFrom(this);
  pen->lambda = lambda;
  return pen;
}

void LassoPenalty::compute(int want, FitContext *fc)
{
	if (want & FF_COMPUTE_FIT) {
    double tmp = 0;
    for (int px = 0; px < params.size(); ++px) {
      double par = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      tmp += penaltyStrength(par, px) * par;
    }
    matrix->data[0] = tmp * fc->est[lambda];
  }
  if (want & FF_COMPUTE_GRADIENT) {
    double la = fc->est[lambda];
    for (int px = 0; px < params.size(); ++px) {
      double par = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      fc->gradZ[ params[px] ] += penaltyStrength(par, px) * copysign(la, fc->est[ params[px] ]);
    }
  }
}

RidgePenalty::RidgePenalty(S4 obj, omxMatrix *mat) : Penalty(obj, mat)
{
  IntegerVector pv = obj.slot("hyperparameters");
  if (pv.size() != 1) mxThrow("RidgePenalty should have 1 hyperparameter (not %d)", pv.size());
  lambda = pv[0];
}

std::unique_ptr<Penalty> RidgePenalty::clone(omxMatrix *mat) const
{
  auto pen = std::make_unique<RidgePenalty>(robj, mat);
  pen->copyFrom(this);
  pen->lambda = lambda;
  return pen;
}

void RidgePenalty::compute(int want, FitContext *fc)
{
	if (want & FF_COMPUTE_FIT) {
    double tmp = 0;
    for (int px = 0; px < params.size(); ++px) {
      double p1 = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      tmp += penaltyStrength(p1, px) * p1 * p1;
    }
    matrix->data[0] = tmp * fc->est[lambda];
  }
  if (want & FF_COMPUTE_GRADIENT) {
    double la2 = 2 * fc->est[lambda];
    for (int px = 0; px < params.size(); ++px) {
      double p1 = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      fc->gradZ[ params[px] ] += la2 * penaltyStrength(p1, px) * p1;
    }
  }
}

ElasticNetPenalty::ElasticNetPenalty(S4 obj, omxMatrix *mat) : Penalty(obj, mat)
{
  IntegerVector pv = obj.slot("hyperparameters");
  if (pv.size() != 2) mxThrow("ElasticNetPenalty should have 2 hyperparameters (not %d)", pv.size());
  alpha = pv[0];
  lambda = pv[1];
}

std::unique_ptr<Penalty> ElasticNetPenalty::clone(omxMatrix *mat) const
{
  auto pen = std::make_unique<ElasticNetPenalty>(robj, mat);
  pen->copyFrom(this);
  pen->alpha = alpha;
  pen->lambda = lambda;
  return pen;
}

void ElasticNetPenalty::compute(int want, FitContext *fc)
{
	if (want & FF_COMPUTE_FIT) {
    double lasso = 0;
    double ridge = 0;
    for (int px = 0; px < params.size(); ++px) {
      double p1 = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      double str = penaltyStrength(p1, px);
      lasso += str * p1;
      ridge += str * p1 * p1;
    }
    double prop = fc->est[ alpha ];
    matrix->data[0] = fc->est[lambda] * ((1-prop) * ridge + prop * lasso);
  }
  if (want & FF_COMPUTE_GRADIENT) {
    double prop = fc->est[ alpha ];
    double la = fc->est[lambda];
    for (int px = 0; px < params.size(); ++px) {
      double p1 = fabs(fc->est[ params[px] ] / scale[px % scale.size()]);
      double str = penaltyStrength(p1, px);
      fc->gradZ[ params[px] ] +=
        prop * str * copysign(la, fc->est[ params[px] ]) + (1-prop) * 2 * la * str * p1;
    }
  }
}

void omxGlobal::importPenalty(omxMatrix *mat, S4 obj, FitContext *fc)
{
  auto &fvg = *findVarGroup(FREEVARGROUP_ALL);

  const char *type = obj.slot("type");
  std::unique_ptr<Penalty> rp;
  if (strEQ(type, "lasso")) rp = std::make_unique<LassoPenalty>(obj, mat);
  else if (strEQ(type, "ridge")) rp = std::make_unique<RidgePenalty>(obj, mat);
  else if (strEQ(type, "elasticNet")) rp = std::make_unique<ElasticNetPenalty>(obj, mat);
  else mxThrow("Unknown type of mxPenalty '%s'", type);
  mat->penalty = std::move(rp);
  omxResizeMatrix(mat, 1, 1);

  List hpr = obj.slot("hpranges");
  for (int rx=0; rx < hpr.size(); ++rx) {
    CharacterVector hprNames = hpr.names();
    const char *varName = hprNames[rx];
    int vx = fvg.lookupVar(varName);
    if (vx == -1) mxThrow("Cannot find free variable named '%s' for penalty '%s'",
                          varName, type);
    auto got = penaltyGrid.find(vx);
    if (got != penaltyGrid.end()) {
      NumericVector g1(got->second);
      NumericVector g2(hpr[rx]);
      if (g1.size() != g2.size()) mxThrow("Different size grids for '%s'", varName);
      for (int gx=0; gx < g1.size(); ++gx) {
        if (g1[gx] != g2[gx]) mxThrow("Different grids for '%s'", varName);
      }
    } else {
      NumericVector g1(hpr[rx]);
      penaltyGrid.emplace(vx, g1);
      if (fc->profiledOutZ[vx]) {
        mxThrow("processPenalties: parameter '%s' is unexpectedly already profiled out",
                varName);
      }
      fc->profiledOutZ[vx] = true;
    }
  }
  fc->calcNumFree();
}
