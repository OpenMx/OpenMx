#ifndef _finiteDifferences_H_
#define _finiteDifferences_H_

template <typename T1, typename T2, typename T3>
void fd_gradient_with_ref(T1 ff, double refFit, Eigen::MatrixBase<T2> &point, int numIter,
			  const double eps, Eigen::MatrixBase<T3> &gradOut)
{
	// http://en.wikipedia.org/wiki/Finite_difference
	Eigen::VectorXd p2;
	p2 = point;
	for (int px=0; px < int(point.size()); ++px) {
		double offset = std::max(fabs(point[px] * eps), eps);
		std::vector<double> Gaprox(numIter);
		for(int k = 0; k < numIter; k++) {
			p2[px] = point[px] + offset;
			Gaprox[k] = (ff(p2) - refFit) / offset;
			offset *= .5;
		}
		for(int m = 1; m < numIter; m++) {						// Richardson Step
			for(int k = 0; k < (numIter - m); k++) {
				// NumDeriv Hard-wires 4s for r here. Why?
				Gaprox[k] = (Gaprox[k+1] * pow(4.0, m) - Gaprox[k])/(pow(4.0, m)-1);
			}
		}
		gradOut[px] = Gaprox[0];
		p2[px] = point[px];
	}
}

template <typename T1, typename T2, typename T3>
void fd_gradient(T1 ff, Eigen::MatrixBase<T2> &point, Eigen::MatrixBase<T3> &gradOut)
{
	const double refFit = ff(point);
	fd_gradient_with_ref(ff, refFit, point, 1, 1e-5, gradOut);
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

template <typename T1, typename T2, typename T3>
void central_gradient(T1 ff, Eigen::MatrixBase<T2> &point, int numIter, const double eps,
		 Eigen::MatrixBase<T3> &gradOut)
{
	Eigen::VectorXd p2;
	for (int px=0; px < int(point.size()); ++px) {
		double offset = std::max(fabs(point[px] * eps), eps);
		std::vector<double> Gaprox(numIter);
		for(int k = 0; k < numIter; k++) {
			p2 = point;
			p2[px] += offset;
			double f1 = ff(p2);
			p2 = point;
			p2[px] -= offset;
			double f2 = ff(p2);
			Gaprox[k] = (f1 - f2) / (2.0 * offset);
			offset *= .5;
		}
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
	case GradientAlgorithm_Forward:
		fd_gradient_with_ref(ff, refFit, point, order, eps, gradOut);
		break;
	case GradientAlgorithm_Central:
		central_gradient(ff, point, order, eps, gradOut);
		break;
	default: Rf_error("Unknown gradient algorithm %d", algo);
	}
}

#endif
