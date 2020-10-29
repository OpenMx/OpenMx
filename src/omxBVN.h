#ifndef _OMX_BVN_H
#define _OMX_BVN_H

#ifdef  __cplusplus
extern "C" {
#endif

double F77_CALL(bvnd)(double*, double*, double*);

#ifdef  __cplusplus
}
#endif

inline double pbivnorm1(double pt1, double pt2, double cor)
{
	return F77_CALL(bvnd)(&pt1, &pt2, &cor);
}

inline double pbivnorm(double lower1, double lower2,
		       double upper1, double upper2, double cor)
{
	double p1 = pbivnorm1(upper1, upper2, cor);
	double p2 = -pbivnorm1(lower1, upper2, cor);
	double p3 = -pbivnorm1(upper1, lower2, cor);
	double p4 = pbivnorm1(lower1, lower2, cor);
	return p1 + p2 + p3 + p4;
}

inline double dbivnorm1(double u, double v, double rho)
{
	const double maxRho = .9999;
	if (fabs(rho) > maxRho) rho = rho<0? -maxRho : maxRho;
	double R = 1 - rho * rho;
	return 1.0 / (M_2PI*sqrt(R)) * exp(-0.5*(u*u - 2*rho*u*v + v*v)/R);
}

inline double dbivnorm(double lower1, double lower2,
		       double upper1, double upper2, double cor)
{
	double out = 0;
	if (upper1 < 100 && upper2 < 100)
		out += dbivnorm1(upper1, upper2, cor);
	if (lower1 > -100 && upper2 < 100)
		out -= dbivnorm1(lower1, upper2, cor);
	if (upper1 < 100 && lower2 > -100)
		out -= dbivnorm1(upper1, lower2, cor);
	if (lower1 > -100 && lower2 > -100)
		out += dbivnorm1(lower1, lower2, cor);
	return out;
}

#endif
