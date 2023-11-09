#ifndef STAN__PROB__DISTRIBUTIONS__MULTIVARIATE__CONTINUOUS__MULTI_NORMAL_SUFFICIENT_HPP
#define STAN__PROB__DISTRIBUTIONS__MULTIVARIATE__CONTINUOUS__MULTI_NORMAL_SUFFICIENT_HPP

#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>

#include <stan/math/version.hpp>
#if STAN_MATH_MAJOR >= 4
#include <stan/math/prim/err/check_ldlt_factor.hpp>
#include <stan/math/prim/err/check_symmetric.hpp>
#include <stan/math/prim/err/check_size_match.hpp>
#include <stan/math/prim/err/check_finite.hpp>
#include <stan/math/prim/err/check_not_nan.hpp>
#include <stan/math/prim/err/check_positive.hpp>
#include <stan/math/prim/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/fun/max_size_mvt.hpp>
#include <stan/math/prim/meta/return_type.hpp>
#include <stan/math/prim/meta/include_summand.hpp>
#else
#include <stan/math/prim/mat/err/check_ldlt_factor.hpp>
#include <stan/math/prim/mat/err/check_symmetric.hpp>
#include <stan/math/prim/scal/err/check_size_match.hpp>
#include <stan/math/prim/scal/err/check_finite.hpp>
#include <stan/math/prim/scal/err/check_not_nan.hpp>
#include <stan/math/prim/scal/err/check_positive.hpp>
#include <stan/math/prim/mat/fun/trace_inv_quad_form_ldlt.hpp>
#include <stan/math/prim/mat/fun/log_determinant_ldlt.hpp>
#include <stan/math/prim/scal/meta/return_type.hpp>
#include <stan/math/prim/scal/meta/max_size_mvt.hpp>
#include <stan/math/prim/scal/meta/include_summand.hpp>
#endif

namespace stan {

  namespace prob {

    template <bool propto, typename T_sample, typename T_loc, typename T_covar>
    typename boost::math::tools::promote_args<T_sample, typename scalar_type<T_loc>::type, T_covar>::type
    multi_normal_sufficient_log(const int sampleSize,
				const Eigen::Matrix<T_sample,Eigen::Dynamic,1>& sampleMu,
				const Eigen::Matrix<T_sample,Eigen::Dynamic,Eigen::Dynamic>& sampleSigma,
				const T_loc& mu,
				const Eigen::Matrix<T_covar,Eigen::Dynamic,Eigen::Dynamic>& Sigma) {
      static const char *function("stan::prob::multi_normal_sufficient_log");
      typedef typename boost::math::tools::promote_args<T_sample, typename scalar_type<T_loc>::type, T_covar>::type param_t;
      typedef param_t lp_type;
      lp_type lp(0.0);
	   
      using stan::math::check_size_match;
      using stan::math::check_finite;
      using stan::math::check_not_nan;
      using stan::math::check_positive;
      using stan::math::check_symmetric;
      using stan::math::check_ldlt_factor;
      using stan::math::LOG_TWO_PI;
      using stan::math::include_summand;

      check_size_match(function,
                       "Rows of covariance parameter", sampleSigma.rows(), 
                       "columns of covariance parameter", sampleSigma.cols());
      check_positive(function, "Covariance matrix rows", sampleSigma.rows());
      check_symmetric(function, "Covariance matrix", sampleSigma);

      check_size_match(function,
                       "Rows of covariance parameter", Sigma.rows(), 
                       "columns of covariance parameter", Sigma.cols());
      check_positive(function, "Covariance matrix rows", Sigma.rows());
      check_symmetric(function, "Covariance matrix", Sigma);
      
      check_size_match(function, 
                       "Size of data location", sampleMu.size(),
                       "size of model location", mu.size());
      check_size_match(function, 
                       "Size of data covariance", sampleSigma.rows(), 
                       "size of model covariance", Sigma.rows());
  
#if STAN_MATH_MAJOR >= 4
      auto ldlt_Sigma = stan::math::make_ldlt_factor(Sigma);
#else
      stan::math::LDLT_factor<param_t,Eigen::Dynamic,Eigen::Dynamic> ldlt_Sigma(Sigma);
#endif
      check_ldlt_factor(function, "LDLT_Factor of covariance parameter", ldlt_Sigma);
      
      double crrctn = (double(sampleSize)-1.0)/double(sampleSize);

      Eigen::Matrix<param_t, Eigen::Dynamic, Eigen::Dynamic> ss;
      ss = mdivide_left_ldlt(ldlt_Sigma, sampleSigma);

      if (include_summand<propto>::value)
	      lp += mu.size() * LOG_TWO_PI * sampleSize;

      if (include_summand<propto, T_covar>::value)
	      lp += log_determinant_ldlt(ldlt_Sigma) * sampleSize;

      if (include_summand<propto, T_covar, T_sample>::value)
	      //lp += ss.trace() * (sampleSize-1.0) * crrctn;
        lp += ss.trace() * (sampleSize-1.0);

      if (include_summand<propto, T_covar, T_loc>::value) {
	Eigen::Matrix<param_t, Eigen::Dynamic, 1> y_minus_mu(mu.size());

	for (int j = 0; j < mu.size(); j++)
	  y_minus_mu(j) = mu(j) - sampleMu(j);

	lp += trace_inv_quad_form_ldlt(ldlt_Sigma, y_minus_mu) * sampleSize;
      }
      return lp * -0.5;
    }
  }
}

#endif
