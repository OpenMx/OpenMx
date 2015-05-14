#ifndef _finiteDifferences_H_
#define _finiteDifferences_H_

// See http://en.wikipedia.org/wiki/Finite_difference

struct forward_difference_grad {
	template <typename T1, typename T2>
	void operator()(T1 ff, double refFit, Eigen::MatrixBase<T2> &point,
			double offset, int px, int numIter, std::vector<double> &Gaprox)
	{
		double orig = point[px];
		for(int k = 0; k < numIter; k++) {
			point[px] = orig + offset;
			Gaprox[k] = (ff(point) - refFit) / offset;
			offset *= .5;
		}
		point[px] = orig;
	};
};

struct central_difference_grad {
	template <typename T1, typename T2>
	void operator()(T1 ff, double refFit, Eigen::MatrixBase<T2> &point,
			double offset, int px, int numIter, std::vector<double> &Gaprox)
	{
		double orig = point[px];
		for(int k = 0; k < numIter; k++) {
			point[px] = orig + offset;
			double f1 = ff(point);
			point[px] = orig - offset;
			double f2 = ff(point);
			Gaprox[k] = (f1 - f2) / (2.0 * offset);
			offset *= .5;
		}
		point[px] = orig;
	};
};

template <typename T1, typename T2, typename T3, typename T4>
void gradientImpl(T1 ff, double refFit, Eigen::MatrixBase<T2> &point, int numIter,
		  const double eps, T4 dfn, Eigen::MatrixBase<T3> &gradOut)
{
	for (int px=0; px < int(point.size()); ++px) {
		double offset = std::max(fabs(point[px] * eps), eps);
		std::vector<double> Gaprox(numIter);
		dfn(ff, refFit, point, offset, px, numIter, Gaprox);
		for(int m = 1; m < numIter; m++) {						// Richardson Step
			for(int k = 0; k < (numIter - m); k++) {
				// NumDeriv Hard-wires 4s for r here. Why?
				Gaprox[k] = (Gaprox[k+1] * pow(4.0, m) - Gaprox[k])/(pow(4.0, m)-1);
			}
		}
		gradOut[px] = Gaprox[0];
	}
}

enum GradientAlgorithm {
	GradientAlgorithm_Forward,
	GradientAlgorithm_Central
};

template <typename T1, typename T2, typename T3>
void gradient_with_ref(GradientAlgorithm algo, int order, double eps, T1 ff, double refFit,
		       Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &gradOut)
{
	switch (algo) {
	case GradientAlgorithm_Forward:{
		forward_difference_grad dfn;
		gradientImpl(ff, refFit, point, order, eps, dfn, gradOut);
		break;}
	case GradientAlgorithm_Central:{
		central_difference_grad dfn;
		gradientImpl(ff, refFit, point, order, eps, dfn, gradOut);
		break;}
	default: Rf_error("Unknown gradient algorithm %d", algo);
	}
}

template <typename T1, typename T2, typename T3, typename T4>
void fd_jacobian(T1 ff, Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &ref,
		 bool refOnly, Eigen::MatrixBase<T4> &jacobiOut)
{
        const double eps = 1e-5;
	ff(point, ref);
	if (refOnly) return;
	
	Eigen::VectorXd p2;
	for (int px=0; px < int(point.size()); ++px) {
		p2 = point;
		double offset = std::max(fabs(p2[px] * eps), eps);
		p2[px] += offset;
		Eigen::VectorXd probe(jacobiOut.cols());
		ff(p2, probe);
		jacobiOut.row(px) = (probe - ref) / offset;
	}
}

#endif
