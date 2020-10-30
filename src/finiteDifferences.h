#ifndef _finiteDifferences_H_
#define _finiteDifferences_H_

#include <limits>

// See http://en.wikipedia.org/wiki/Finite_difference

struct forward_difference_jacobi {
	template <typename T1, typename T2, typename T3, typename T4>
	void operator()(T1 ff, Eigen::MatrixBase<T4> &refFit, Eigen::MatrixBase<T2> &point,
			double offset, int px, int numIter, Eigen::MatrixBase<T3> &Gaprox)
	{
		double orig = point[px];
		Eigen::MatrixXd result(refFit.rows(), refFit.cols());
		for(int k = 0; k < numIter; k++) {
			point[px] = orig + offset;
			ff(point, result);
			Eigen::Map<Eigen::MatrixXd> Gaprox1(&Gaprox.coeffRef(0,k), refFit.rows(), refFit.cols());
			Gaprox1 = (result - refFit) / offset;
			offset *= .5;
		}
		point[px] = orig;
	}
};

struct central_difference_jacobi {
	template <typename T1, typename T2, typename T3, typename T4>
	void operator()(T1 ff, Eigen::MatrixBase<T4> &refFit, Eigen::MatrixBase<T2> &point,
			double offset, int px, int numIter, Eigen::MatrixBase<T3> &Gaprox)
	{
		double orig = point[px];
		Eigen::MatrixXd result1(refFit.rows(), refFit.cols());
		Eigen::MatrixXd result2(refFit.rows(), refFit.cols());
		for(int k = 0; k < numIter; k++) {
			point[px] = orig + offset;
			ff(point, result1);
			point[px] = orig - offset;
			ff(point, result2);
			Eigen::Map<Eigen::MatrixXd> Gaprox1(&Gaprox.coeffRef(0,k), refFit.rows(), refFit.cols());
			Gaprox1 = (result1 - result2) / (2.0 * offset);
			offset *= .5;
		}
		point[px] = orig;
	}
};

template <bool initialized, typename T1, typename T2, typename T3, typename T4, typename T5>
void jacobianImpl(T1 ff,  Eigen::MatrixBase<T2> &ref, Eigen::MatrixBase<T3> &point, int px,
                  int numIter, const double eps, T4 dfn, Eigen::MatrixBase<T5> &jacobiOut, int dest)
{
	Eigen::Map<Eigen::VectorXd> jacobiOut1(&jacobiOut.coeffRef(0,dest), ref.size());
	if (!initialized || std::isnan(jacobiOut1.sum())) {
		double offset = std::max(fabs(point[px] * eps), eps);
		Eigen::MatrixXd Gaprox(ref.size(), numIter);
		dfn(ff, ref, point, offset, px, numIter, Gaprox);
		for(int m = 1; m < numIter; m++) {						// Richardson Step
			for(int k = 0; k < (numIter - m); k++) {
				// NumDeriv Hard-wires 4s for r here. Why?
				Gaprox.col(k) = (Gaprox.col(k+1) * pow(4.0, m) - Gaprox.col(k))/(pow(4.0, m)-1);
			}
		}
		for(int i=0; i<jacobiOut1.rows(); i++){
			if (!initialized || std::isnan(jacobiOut1(i))) {
				jacobiOut1(i) = Gaprox(i,0);
			}
		}
	}
}

// fd_jacobian is deprecated. Should replace with JacobianGadget TODO
// all parameters
template <bool initialized, typename T1, typename T2, typename T3, typename T4>
void fd_jacobian(GradientAlgorithm algo, int numIter, double eps, T1 ff, Eigen::MatrixBase<T2> &ref,
	      Eigen::MatrixBase<T3> &point, Eigen::MatrixBase<T4> &jacobiOut)
{
  //mxLog("fd_jacobian %dx%d", jacobiOut.rows(), jacobiOut.cols());
	// TODO refactor, evaluate jacobian in parallel
	switch (algo) {
	case GradientAlgorithm_Forward:{
		forward_difference_jacobi dfn;
		for (int px=0; px < int(point.size()); ++px) {
			jacobianImpl<initialized>(ff, ref, point, px, numIter, eps, dfn, jacobiOut, px);
		}
		break;}
	case GradientAlgorithm_Central:{
		central_difference_jacobi dfn;
		for (int px=0; px < int(point.size()); ++px) {
			jacobianImpl<initialized>(ff, ref, point, px, numIter, eps, dfn, jacobiOut, px);
		}
		break;}
	default: mxThrow("Unknown gradient algorithm %d", algo);
	}
}

// single parameter
template <bool initialized, typename T1, typename T2, typename T3, typename T4>
void fd_jacobian1(GradientAlgorithm algo, int numIter, double eps, T1 ff, Eigen::MatrixBase<T2> &ref,
		  Eigen::MatrixBase<T3> &point, int px, Eigen::MatrixBase<T4> &jacobiOut)
{
	switch (algo) {
	case GradientAlgorithm_Forward:{
		forward_difference_jacobi dfn;
		jacobianImpl<initialized>(ff, ref, point, px, numIter, eps, dfn, jacobiOut, 0);
		break;}
	case GradientAlgorithm_Central:{
		central_difference_jacobi dfn;
		jacobianImpl<initialized>(ff, ref, point, px, numIter, eps, dfn, jacobiOut, 0);
		break;}
	default: mxThrow("Unknown gradient algorithm %d", algo);
	}
}

template <class Derived>
struct finite_difference_jacobian {
	double *refData;
  int refRows;

  void setRef(Eigen::Ref<Eigen::ArrayXd> ref)
  {
    refData = ref.data();
    refRows = ref.rows();
  }
  Eigen::Map< Eigen::ArrayXd > getRef() {
    Eigen::Map< Eigen::ArrayXd > ref(refData, refRows);
    return ref;
  }
	int thrId;
	double *point;
	double orig;

	template <typename T1>
	void approx(T1 ff, double offset, int px, double *out)
	{ static_cast<Derived *>(this)->approx(ff, offset, px, out); }

	template <typename T1>
	void operator()(T1 ff, int _thrId, double *_point,
			double offset, int px, int numIter, double *grid, int verbose)
	{
		thrId = _thrId;
		point = _point;
		orig = point[px];
    Eigen::Map< Eigen::ArrayXXd > Egrid(grid, refRows, numIter);

		for(int k = 0; k < numIter;) {
      while (1) {
				approx(ff, offset, px, &Egrid.coeffRef(0, k));
        offset *= .5;
        if (k==0 && !Egrid.col(k).isFinite().all()) {
          if (offset > std::numeric_limits<double>::epsilon()) {
            if (verbose + OMX_DEBUG >= 1) {
              mxLog("finite differences[%d]: retry with offset %.4g",
                    px, offset);
            }
            continue;
          }
        }
        break;
      }
			k += 1;
		}
    for(int m = 1; m < numIter; m++) {	// Richardson Step
      for(int k = 0; k < (numIter - m); k++) {
        // NumDeriv Hard-wires 4s for r here. Why?
        Egrid.col(k) = (Egrid.col(k+1) * pow(4.0, m) - Egrid.col(k))/(pow(4.0, m)-1);
      }
    }
		point[px] = orig;
	}
};

struct forward_difference_jacobian : finite_difference_jacobian<forward_difference_jacobian> {
	template <typename T1>
	void approx(T1 ff, double offset, int px, double *out)
	{
    auto ref = getRef();
    Eigen::Map< Eigen::ArrayXd > Eout(out, ref.rows());
		Eigen::ArrayXd result(ref.rows());
		point[px] = orig + offset;
		ff(point, thrId, result);
		Eout.derived() = (result - ref) / offset;
	}
};

struct central_difference_jacobian  : finite_difference_jacobian<central_difference_jacobian> {
	template <typename T1>
	void approx(T1 ff, double offset, int px, double *out)
	{
    auto ref = getRef();
    Eigen::Map< Eigen::ArrayXd > Eout(out, ref.rows());
		Eigen::ArrayXd result1(ref.rows());
		Eigen::ArrayXd result2(ref.rows());
		point[px] = orig + offset;
		ff(point, thrId, result1);
		point[px] = orig - offset;
		ff(point, thrId, result2);
		Eout.derived() = (result1 - result2) / (2.0 * offset);
	}
};

class JacobianGadget {
  const char *name;
	const int ELAPSED_HISTORY_SIZE;
	const int maxAvailThreads;
	const int numFree;
	GradientAlgorithm algo;
	int numIter;
	double eps;
	int verbose;
  bool used;
	int curElapsed;
	int numThreadsBookmark;
	int curNumThreads;
  nanotime_t startTime;
	std::vector<nanotime_t> elapsed0;
	std::vector<nanotime_t> elapsed1;
	Eigen::ArrayXXd grid;
	Eigen::MatrixXd thrPoint;

	template <typename T1, typename T3, typename T4, typename T5>
	void myJacobianImpl(T1 ff, Eigen::MatrixBase<T3> &point,
                      T4 dfn, bool initialized, T5 &out)
	{
		thrPoint.resize(point.size(), curNumThreads);
		thrPoint.colwise() = point;

#pragma omp parallel for num_threads(curNumThreads)
		for (int px=0; px < int(point.size()); ++px) {
			int thrId = omp_get_thread_num();
			int thrSelect = curNumThreads==1? -1 : thrId;
			double offset = std::max(fabs(point[px] * eps), eps);
      if (!initialized || !out.col(px).array().isFinite().all()) {
        try {
          dfn[thrId](ff, thrSelect, &thrPoint.coeffRef(0, thrId), offset, px,
                     numIter, &grid.coeffRef(0,thrId), verbose);
        } catch (const std::exception& e) {
          omxRaiseErrorf("%s", e.what());
        } catch (...) {
          omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
        }
        if (!initialized) {
          out.col(px) = grid.block(0, thrId, out.rows(), 1);
        } else {
          for (int rx=0; rx < out.rows(); ++rx) {
            if (!std::isfinite(out(rx,px))) out(rx,px) = grid(rx,thrId);
          }
        }
      }
		}
	}

  void start()
  {
    used = true;
		startTime = get_nanotime();
		curNumThreads = std::max(1, numThreadsBookmark - curElapsed % 2);
  }

  void finish()
  {
		double el1 = (get_nanotime() - startTime);
    if (curElapsed < ELAPSED_HISTORY_SIZE * 2) {
      if (verbose >= 2) {
        mxLog("%s: test[%d] curNumThreads=%d %fms",
              name, curElapsed, curNumThreads, el1/1000000.0);
      }
			int repl = curElapsed / 2;
			if (curElapsed % 2) {
				elapsed1[repl] = el1;
			} else {
				elapsed0[repl] = el1;
			}
			curElapsed += 1;
			if (curElapsed == ELAPSED_HISTORY_SIZE * 2) {
				std::sort(elapsed0.begin(), elapsed0.end());
				std::sort(elapsed1.begin(), elapsed1.end());
				double e0 = elapsed0[elapsed0.size()/2];
				double e1 = elapsed1[elapsed1.size()/2];
        if (verbose) {
          mxLog("%s: took %fms with %d threads and %fms with %d threads",
                name, e0/1000000.0, numThreadsBookmark, e1/1000000.0, std::max(1, numThreadsBookmark-1));
        }
				if (e0 > e1 && numThreadsBookmark > 1) {
					numThreadsBookmark = std::max(numThreadsBookmark - 1, 1);
					if (numThreadsBookmark > 1) curElapsed = 0;
				}
        if (verbose && curElapsed > 0) {
          mxLog("%s: looks like %d threads offer the best performance",
                name, numThreadsBookmark);
        }
			}
		}
  }

 public:

	JacobianGadget(int numThreads, int _numFree) :
    name("JacobianGadget"), ELAPSED_HISTORY_SIZE(3), maxAvailThreads(numThreads),
    numFree(_numFree), algo(Global->gradientAlgo), numIter(Global->gradientIter),
    eps(Global->gradientStepSize)
	{
		verbose = numThreads>1 && Global->parallelDiag;
    used = false;
		curElapsed = 0;
		numThreadsBookmark = std::min(numThreads, numFree); // could break work into smaller pieces TODO
    if (numThreadsBookmark < 1) numThreadsBookmark = 1;
		if (numThreadsBookmark == 1) {
			curElapsed = ELAPSED_HISTORY_SIZE * 2;
		} else {
			elapsed0.resize(ELAPSED_HISTORY_SIZE);
			elapsed1.resize(ELAPSED_HISTORY_SIZE);
		}
	}
  ~JacobianGadget()
  {
    if (used) {
      diagParallel(OMX_DEBUG, "%s: used %d/%d threads for %d free parameters",
                   name, numThreadsBookmark, maxAvailThreads, numFree);
    } else {
      diagParallel(OMX_DEBUG, "%s: not used (analytic?)", name);
    }
  }

  void setAlgoOptions(GradientAlgorithm _algo,	int _numIter, double _eps)
  {
    algo = _algo;
    numIter = _numIter;
    eps = _eps;
  }

  // ff -- the function to evaluate
  // ref -- reference function value
  // point -- location where to take the jacobian
  // initialized -- whether to do all entries or only NaN entries
  // jacobiOut -- output
template <typename T1, typename T2, typename T3, typename T4>
	void operator()(T1 ff, T2 &ref,
	      Eigen::MatrixBase<T3> &point, bool initialized, T4 &jacobiOut)
	{
		if (point.size() != numFree) mxThrow("%s line %d: expecting %d parameters, got %d",
                                         __FILE__, __LINE__, numFree, point.size());
    if (jacobiOut.cols() != point.size()) mxThrow("%s line %d: non-conformable %d != %d",
                                                  __FILE__, __LINE__, jacobiOut.cols(), point.size());
    if (ref.size() != jacobiOut.rows()) OOPS;

    start();
		grid.resize(numIter * jacobiOut.rows(), curNumThreads);

    switch (algo) {
    case GradientAlgorithm_Forward:{
      std::vector<forward_difference_jacobian> dfn(curNumThreads);
      for (auto &fd : dfn) fd.setRef(ref);
      myJacobianImpl(ff, point, dfn, initialized, jacobiOut);
      break;}
    case GradientAlgorithm_Central:{
      std::vector<central_difference_jacobian> dfn(curNumThreads);
      for (auto &fd : dfn) fd.setRef(ref);
      myJacobianImpl(ff, point, dfn, initialized, jacobiOut);
      break;}
    default: mxThrow("%s: Unknown algorithm %d", name, algo);
    }

    finish();

    //for (int rx=0; rx < jacobiOut.rows(); ++rx) robustify(jacobiOut.row(rx)); // TODO
	}

  bool needRefFit() const { return algo == GradientAlgorithm_Forward; }
};

template <typename T5>
void robustifyInplace(Eigen::MatrixBase<T5> &out)
{
  Eigen::ArrayXXd absOut = out.array().abs();
  std::nth_element(absOut.data(), absOut.data() + absOut.size()/2,
                   absOut.data()+absOut.size());
  double m1 = std::max(absOut.data()[absOut.size()/2], 1.0);
  double big = 1e4 * m1;
  int adj=0;
  for (int cx=0; cx < out.cols(); ++cx) {
    for (int rx=0; rx < out.rows(); ++rx) {
      if (fabs(out(rx,cx)) < big) continue;
      bool neg = out(rx,cx) < 0;
      double gg = m1;
      if (neg) gg = -gg;
      out(rx,cx) = gg;
      ++adj;
    }
  }
  if (false && adj) {
    mxLog("robustify: %d outlier", adj);
    mxPrintMat("robust", out);
  }
}

#endif
