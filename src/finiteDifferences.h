#ifndef _finiteDifferences_H_
#define _finiteDifferences_H_

#include <limits>

// See http://en.wikipedia.org/wiki/Finite_difference

enum GradientAlgorithm {
	GradientAlgorithm_Forward,
	GradientAlgorithm_Central
};

template <class Derived>
struct finite_difference_grad {
	double refFit;
	int thrId;
	double *point;
	double orig;

	template <typename T1>
	double approx(T1 ff, double offset, int px)
	{
		return static_cast<Derived *>(this)->approx(ff, offset, px);
	}

	template <typename T1>
	void operator()(T1 ff, double _refFit, int _thrId, double *_point,
			double offset, int px, int numIter, double *Gaprox, int verbose)
	{
		refFit = _refFit;
		thrId = _thrId;
		point = _point;
		orig = point[px];

		for(int k = 0; k < numIter;) {
			double grad = 0;
			if (offset > std::numeric_limits<double>::epsilon()) {
				grad = approx(ff, offset, px);
			}
			offset *= .5;
			if (!std::isfinite(grad)) {
				if (verbose + OMX_DEBUG >= 1) {
					mxLog("finite differences[%d]: retry with offset %.4g",
					      px, offset);
				}
				continue;
			}
			Gaprox[k] = grad;
			k += 1;
		}
		point[px] = orig;
	}
};

struct forward_difference_grad : finite_difference_grad<forward_difference_grad> {
	template <typename T1>
	double approx(T1 ff, double offset, int px)
	{
		point[px] = orig + offset;
		double f1 = ff(point, thrId);
		return (f1 - refFit) / offset;
	}
};

struct central_difference_grad  : finite_difference_grad<central_difference_grad> {
	template <typename T1>
	double approx(T1 ff, double offset, int px)
	{
		point[px] = orig + offset;
		double f1 = ff(point, thrId);
		point[px] = orig - offset;
		double f2 = ff(point, thrId);
		return (f1 - f2) / (2.0 * offset);
	}
};

class GradientWithRef {
	const int ELAPSED_HISTORY_SIZE;
	const int maxAvailThreads;
	const int numFree;
	const GradientAlgorithm algo;
	const int numIter;
	const double eps;
	int verbose;
	int curElapsed;
	int numGradientThreads;
	int curNumThreads;
	std::vector<nanotime_t> elapsed0;
	std::vector<nanotime_t> elapsed1;
	Eigen::MatrixXd grid;
	Eigen::MatrixXd thrPoint;

	template <typename T1, typename T2, typename T3, typename T4>
	void gradientImpl(T1 ff, double refFit, Eigen::MatrixBase<T2> &point,
			  T4 dfn, Eigen::MatrixBase<T3> &gradOut)
	{
		thrPoint.resize(point.size(), curNumThreads);
		thrPoint.colwise() = point;

#pragma omp parallel for num_threads(curNumThreads)
		for (int px=0; px < int(point.size()); ++px) {
			int thrId = omp_get_thread_num();
			int thrSelect = curNumThreads==1? -1 : thrId;
			double offset = std::max(fabs(point[px] * eps), eps);
			try {
				dfn[thrId](ff, refFit, thrSelect, &thrPoint.coeffRef(0, thrId), offset, px,
					   numIter, &grid.coeffRef(0,px), verbose);
			} catch (const std::exception& e) {
				omxRaiseErrorf("%s", e.what());
			} catch (...) {
				omxRaiseErrorf("%s line %d: unknown exception", __FILE__, __LINE__);
			}
			// push down into per-thread code TODO
			for(int m = 1; m < numIter; m++) {	// Richardson Step
				for(int k = 0; k < (numIter - m); k++) {
					// NumDeriv Hard-wires 4s for r here. Why?
					grid(k,px) = (grid(k+1,px) * pow(4.0, m) - grid(k,px))/(pow(4.0, m)-1);
				}
			}
			gradOut[px] = grid(0,px);
		}
	}

 public:
	GradientWithRef(int numThreads, int _numFree, GradientAlgorithm _algo, int _numIter, double _eps) :
    ELAPSED_HISTORY_SIZE(3), maxAvailThreads(numThreads), numFree(_numFree), algo(_algo),
    numIter(_numIter), eps(_eps)
	{
		verbose = numThreads>1 && Global->parallelDiag;
		curElapsed = 0;
		numGradientThreads = std::min(numThreads, numFree); // could break work into smaller pieces TODO
		if (numGradientThreads == 1) {
			curElapsed = ELAPSED_HISTORY_SIZE * 2;
		} else {
			elapsed0.resize(ELAPSED_HISTORY_SIZE);
			elapsed1.resize(ELAPSED_HISTORY_SIZE);
		}
		grid.resize(numIter, numFree);
	}
  ~GradientWithRef()
  {
    diagParallel(OMX_DEBUG, "Gradient used %d/%d threads for %d free parameters",
                 numGradientThreads, maxAvailThreads, numFree);
  }

	template <typename T1, typename T2, typename T3>
	void operator()(T1 ff, double refFit,
			Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &gradOut)
	{
		if (point.size() != numFree) mxThrow("%s line %d: expecting %d parameters, got %d",
						      __FILE__, __LINE__, numFree, point.size());

		nanotime_t startTime = get_nanotime();
		curNumThreads = std::max(1, numGradientThreads - curElapsed % 2);

		switch (algo) {
		case GradientAlgorithm_Forward:{
			std::vector<forward_difference_grad> dfn(curNumThreads);
			gradientImpl(ff, refFit, point, dfn, gradOut);
			break;}
		case GradientAlgorithm_Central:{
			std::vector<central_difference_grad> dfn(curNumThreads);
			gradientImpl(ff, refFit, point, dfn, gradOut);
			break;}
		default: mxThrow("Unknown gradient algorithm %d", algo);
		}

		double el1 = (get_nanotime() - startTime);
    if (curElapsed < ELAPSED_HISTORY_SIZE * 2) {
      if (verbose >= 2) {
        mxLog("gradient: test[%d] curNumThreads=%d %fms",
              curElapsed, curNumThreads, el1/1000000.0);
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
          mxLog("gradient took %fms with %d threads and %fms with %d threads",
                e0/1000000.0, numGradientThreads, e1/1000000.0, std::max(1, numGradientThreads-1));
        }
				if (e0 > e1 && numGradientThreads > 1) {
					numGradientThreads = std::max(numGradientThreads - 1, 1);
					if (numGradientThreads > 1) curElapsed = 0;
				}
        if (verbose && curElapsed > 0) {
          mxLog("gradient: looks like %d threads offers the best performance",
                numGradientThreads);
        }
			}
		}
	}
};

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

// all parameters
template <bool initialized, typename T1, typename T2, typename T3, typename T4>
void fd_jacobian(GradientAlgorithm algo, int numIter, double eps, T1 ff, Eigen::MatrixBase<T2> &ref,
	      Eigen::MatrixBase<T3> &point, Eigen::MatrixBase<T4> &jacobiOut)
{
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

#endif
