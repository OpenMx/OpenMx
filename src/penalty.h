/*
 *  Copyright 2021 by the individuals mentioned in the source code history
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

#ifndef u_PENALTY_H_
#define u_PENALTY_H_

class omxMatrix;

class Penalty {
protected:
  S4 robj;
	omxMatrix *matrix;
  IntegerVector params;
  NumericVector epsilon;
  NumericVector scale;
  double smoothProportion;
  void copyFrom(const Penalty *pen);

 public:
  Penalty(S4 _obj, omxMatrix *_mat);
  virtual ~Penalty();
  double penaltyStrength(double absPar, int px) const;
  int countNumZero(FitContext *fc) const;
	virtual void compute(int ffcompute, FitContext *fc)=0;
	const char *name() const;
  virtual std::unique_ptr<Penalty> clone(omxMatrix *mat) const = 0;
  double getValue() const;
};

class LassoPenalty : public Penalty {
  typedef Penalty super;
  int lambda;
public:
  LassoPenalty(S4 obj, omxMatrix *_mat);
	virtual void compute(int ffcompute, FitContext *fc) override;
  virtual std::unique_ptr<Penalty> clone(omxMatrix *mat) const override;
};

class RidgePenalty : public Penalty {
  typedef Penalty super;
  int lambda;
public:
  RidgePenalty(S4 obj, omxMatrix *_mat);
	virtual void compute(int ffcompute, FitContext *fc) override;
  virtual std::unique_ptr<Penalty> clone(omxMatrix *mat) const override;
};

class ElasticNetPenalty : public Penalty {
  typedef Penalty super;
  int lambda, alpha;
public:
  ElasticNetPenalty(S4 obj, omxMatrix *_mat);
	virtual void compute(int ffcompute, FitContext *fc) override;
  virtual std::unique_ptr<Penalty> clone(omxMatrix *mat) const override;
};

#endif
