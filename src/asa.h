#ifndef _ASA_H_
#define _ASA_H_
#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************************
* Adaptive Simulated Annealing (ASA)
* Lester Ingber <ingber@ingber.com>
* Copyright (c) 1987-2017 Lester Ingber.  All Rights Reserved.
* ASA-LICENSE file has the license that must be included with ASA code.
***********************************************************************/

  /* $Id: asa.h,v 30.23 2017/11/20 23:17:53 ingber Exp ingber $ */

  /* asa.h for Adaptive Simulated Annealing */

#include "asa_usr_asa.h"

#define ZERO			((double) 0.0)
#define ONE			((double) 1.0)
#define TWO			((double) 2.0)
#define TEN			((double) 10.0)
#define HALF			((double) 0.5)

#define NORMAL_EXIT			((int) 0)
#define P_TEMP_TOO_SMALL		((int) 1)
#define C_TEMP_TOO_SMALL		((int) 2)
#define COST_REPEATING			((int) 3)
#define TOO_MANY_INVALID_STATES		((int) 4)
#define IMMEDIATE_EXIT			((int) 5)
#define INVALID_USER_INPUT		((int) 7)
#define INVALID_COST_FUNCTION		((int) 8)
#define INVALID_COST_FUNCTION_DERIV	((int) 9)
#define CALLOC_FAILED			((int) -1)

#ifndef TIME_STD
#define TIME_STD FALSE
#endif

#ifndef TIME_GETRUSAGE
#define TIME_GETRUSAGE TRUE
#endif

#if TIME_CALC
#if TIME_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#if TIME_STD
#include <sys/syscall.h>
#endif                          /* TIME_STD */
#else                           /* TIME_GETRUSAGE */
#if TRUE                        /* change to FALSE for SunOS 4.1.x */
#include <time.h>
#else
#include </usr/5include/time.h>
#endif
#endif                          /* TIME_GETRUSAGE */
#endif                          /* TIME_CALC */

  /* Set this to TRUE to override the P_TEMP_TOO_SMALL test */
#ifndef NO_PARAM_TEMP_TEST
#define NO_PARAM_TEMP_TEST FALSE
#endif

  /* Set this to TRUE to override the C_TEMP_TOO_SMALL test */
#ifndef NO_COST_TEMP_TEST
#define NO_COST_TEMP_TEST FALSE
#endif

#ifndef SYSTEM_CALL
#define SYSTEM_CALL TRUE
#endif

  /* Printing Options */

#ifndef ASA_PRINT
#define ASA_PRINT TRUE
#endif

#if ASA_PRINT
#else
#if ASA_SAMPLE
#define ASA_PRINT TRUE
#endif
#endif

#ifndef ASA_OUT
#define ASA_OUT "asa_out"
#endif

#ifndef DROPPED_PARAMETERS
#define DROPPED_PARAMETERS FALSE
#endif

  /* You can set ASA_PRINT_INTERMED to TRUE to print out
     intermediate data when SELF_OPTIMIZE is set to TRUE */
#ifndef ASA_PRINT_INTERMED
#if SELF_OPTIMIZE
#define ASA_PRINT_INTERMED FALSE
#else
#define ASA_PRINT_INTERMED TRUE
#endif
#endif

#ifndef ASA_PRINT_MORE
#define ASA_PRINT_MORE FALSE
#endif

  /* The state of the system in terms of parameters and function value */
  typedef struct {
    double cost;
    double *parameter;
#if ASA_PARALLEL
#if USER_ACCEPTANCE_TEST
    int par_user_accept_flag;
    int par_cost_accept_flag;
#endif
#endif
  } STATE;

#if ASA_PARALLEL
  /* parallel generated states */
  STATE *gener_block_state_qsort;
#endif

  /* essential MACROS */

#if USER_REANNEAL_PARAMETERS
#else
  /* FUNCTION_REANNEAL_PARAMS(temperature, tangent, max_tangent)
     determines the reannealed temperature. */
#define FUNCTION_REANNEAL_PARAMS(temperature, tangent, max_tangent) (temperature * (max_tangent / tangent))
#endif

  /* IABS(i)
     absolute value for integers, in stdlib.h on _some_ machines */
#define IABS(i) ((i) < 0? -(i) : (i))

  /*  NO_REANNEAL(x)
     can determine whether to calculate derivatives. */
#define NO_REANNEAL(x)	(IABS(parameter_type[x]) == 2)

  /* VFOR
     is a simple macro to iterate on each parameter index. */

#define VFOR(index_v) for (index_v = 0; index_v < *number_parameters; ++index_v)

#if CHECK_EXPONENT
  /* EXPONENT_CHECK
     checks that an exponent x is within a valid range and,
     if not, adjusts its magnitude to fit in the range. */
#define MIN_EXPONENT (0.9 * F_LOG ((double) MIN_DOUBLE))
#define MAX_EXPONENT (0.9 * F_LOG ((double) MAX_DOUBLE))
#define EXPONENT_CHECK(x) ((x) < MIN_EXPONENT ? MIN_EXPONENT : ((x) > MAX_EXPONENT ? MAX_EXPONENT : (x)))
#else
#define EXPONENT_CHECK(x) (x)
#endif                          /* CHECK_EXPONENT */

  /* PARAMETER_RANGE_TOO_SMALL(x)
     checks if the range of parameter x is too small to work with.
     If user_cost_function changes the parameter ranges,
     this test could be used to adaptively bypass
     some parameters, e.g., depending on constraints. */
#define PARAMETER_RANGE_TOO_SMALL(x) (fabs(parameter_minimum[x] - parameter_maximum[x]) < (double) EPS_DOUBLE)

  /* INTEGER_PARAMETER(x)
     determines if the parameter is an integer type. */
#define INTEGER_PARAMETER(x) (parameter_type[x] > 0)

  /* ROW_COL_INDEX(i, j)
     converts from row i, column j to an index. */
#define ROW_COL_INDEX(i, j) ((i) + *number_parameters * (j))

#if HAVE_ANSI

  /* asa function prototypes */
  void accept_new_state (double (*user_random_generator) (LONG_INT *),
                         LONG_INT * seed,
                         double *parameter_minimum,
                         double *parameter_maximum,
                         double *current_cost_temperature,
#if ASA_SAMPLE
                         double *current_user_parameter_temp,
#endif
                         ALLOC_INT * number_parameters,
                         LONG_INT * recent_number_acceptances,
                         LONG_INT * number_accepted,
                         LONG_INT * index_cost_acceptances,
                         LONG_INT * number_acceptances_saved,
                         LONG_INT * recent_number_generated,
                         LONG_INT * number_generated,
                         LONG_INT * index_parameter_generations,
                         STATE * current_generated_state,
                         STATE * last_saved_state,
#if ASA_SAMPLE
                         FILE * ptr_asa_out,
#endif
                         USER_DEFINES * OPTIONS);

  int generate_new_state (double (*user_random_generator) (LONG_INT *),
                          LONG_INT * seed,
                          double *parameter_minimum,
                          double *parameter_maximum,
                          double *current_parameter_temperature,
#if USER_GENERATING_FUNCTION
                          double *initial_user_parameter_temp,
                          double *temperature_scale_parameters,
#endif
                          ALLOC_INT * number_parameters,
                          int *parameter_type,
                          STATE * current_generated_state,
                          STATE * last_saved_state, USER_DEFINES * OPTIONS);
#if ASA_PARALLEL
  int generate_new_state_par (double (*user_random_generator) (LONG_INT *),
                              LONG_INT * seed,
                              double *parameter_minimum,
                              double *parameter_maximum,
                              double *current_parameter_temperature,
#if USER_GENERATING_FUNCTION
                              double *initial_user_parameter_temp,
                              double *temperature_scale_parameters,
#endif
                              ALLOC_INT * number_parameters,
                              int *parameter_type,
                              LONG_INT i_prll,
                              STATE * gener_block_state,
                              STATE * last_saved_state,
                              USER_DEFINES * OPTIONS);
#endif                          /* ASA_PARALLEL */

  void reanneal (double *parameter_minimum,
                 double *parameter_maximum,
                 double *tangents,
                 double *maximum_tangent,
                 double *current_cost_temperature,
                 double *initial_cost_temperature,
                 double *temperature_scale_cost,
                 double *current_user_parameter_temp,
                 double *initial_user_parameter_temp,
                 double *temperature_scale_parameters,
                 ALLOC_INT * number_parameters,
                 int *parameter_type,
                 LONG_INT * index_cost_acceptances,
                 LONG_INT * index_parameter_generations,
                 STATE * last_saved_state,
                 STATE * best_generated_state, USER_DEFINES * OPTIONS);

  void
    cost_derivatives (double (*user_cost_function)

                       
                      (double *, double *, double *, double *, double *,
                       ALLOC_INT *, int *, int *, int *, USER_DEFINES *),
                      double *parameter_minimum, double *parameter_maximum,
                      double *tangents, double *curvature,
                      double *maximum_tangent, ALLOC_INT * number_parameters,
                      int *parameter_type, int *exit_status,
                      int *curvature_flag, int *valid_state_generated_flag,
                      LONG_INT * number_invalid_generated_states,
                      STATE * current_generated_state,
                      STATE * best_generated_state, FILE * ptr_asa_out,
                      USER_DEFINES * OPTIONS);

  double generate_asa_state (double (*user_random_generator) (LONG_INT *),
                             LONG_INT * seed, double *temp);

  int
    asa_exit (double (*user_cost_function)

               
              (double *, double *, double *, double *, double *, ALLOC_INT *,
               int *, int *, int *, USER_DEFINES *), double *final_cost,
              double *parameter_initial_final, double *parameter_minimum,
              double *parameter_maximum, double *tangents, double *curvature,
              double *maximum_tangent, double *current_cost_temperature,
              double *initial_user_parameter_temp,
              double *current_user_parameter_temp,
              double *accepted_to_generated_ratio,
              ALLOC_INT * number_parameters, int *parameter_type,
              int *valid_state_generated_flag, int *exit_status,
              ALLOC_INT * index_exit_v, ALLOC_INT * start_sequence,
              LONG_INT * number_accepted,
              LONG_INT * best_number_accepted_saved,
              LONG_INT * index_cost_acceptances, LONG_INT * number_generated,
              LONG_INT * number_invalid_generated_states,
              LONG_INT * index_parameter_generations,
              LONG_INT * best_number_generated_saved,
              STATE * current_generated_state, STATE * last_saved_state,
              STATE * best_generated_state, FILE * ptr_asa_out,
              USER_DEFINES * OPTIONS);

  void Exit_ASA (char *statement);

  int asa_test_asa_options (LONG_INT * seed,
                            double *parameter_initial_final,
                            double *parameter_minimum,
                            double *parameter_maximum,
                            double *tangents,
                            double *curvature,
                            ALLOC_INT * number_parameters,
                            int *parameter_type,
                            int *valid_state_generated_flag,
                            int *exit_status,
                            FILE * ptr_asa_out, USER_DEFINES * OPTIONS);

  int cost_function_test (double cost,
                          double *parameter,
                          double *parameter_minimum,
                          double *parameter_maximum,
                          ALLOC_INT * number_parameters,
                          double *xnumber_parameters);

  void print_string (FILE * ptr_asa_out, char *string);
  void print_string_index (FILE * ptr_asa_out, char *string, ALLOC_INT index);

#if ASA_PRINT
  void print_state (double *parameter_minimum,
                    double *parameter_maximum,
                    double *tangents,
                    double *curvature,
                    double *current_cost_temperature,
                    double *current_user_parameter_temp,
                    double *accepted_to_generated_ratio,
                    ALLOC_INT * number_parameters,
                    int *curvature_flag,
                    LONG_INT * number_accepted,
                    LONG_INT * index_cost_acceptances,
                    LONG_INT * number_generated,
                    LONG_INT * number_invalid_generated_states,
                    STATE * last_saved_state,
                    STATE * best_generated_state,
                    FILE * ptr_asa_out, USER_DEFINES * OPTIONS);

  void print_asa_options (FILE * ptr_asa_out, USER_DEFINES * OPTIONS);
#endif                          /* ASA_PRINT */

#if TIME_CALC
#if TIME_GETRUSAGE
  void aux_print_time (struct timeval *time, char *message,
                       FILE * ptr_asa_out);
#if TIME_STD
  int syscall (int sys_option, int who, struct rusage *usage);
#else
  int getrusage (int who, struct rusage *usage);
#endif                          /* TIME_STD */
#else                           /* TIME_GETRUSAGE */
  void aux_print_time (clock_t time, char *message, FILE * ptr_asa_out);
#if FALSE                       /* change to TRUE for SunOS 4.1.x */
  clock_t clock ();
#endif
#endif                          /* TIME_GETRUSAGE */
#endif                          /* TIME_CALC */

#if MULTI_MIN
  static int multi_compare (const void *cost_ii, const void *cost_jj);
  double *multi_cost_qsort;
#endif

#if ASA_PARALLEL
  static int sort_parallel (const void *cost_ii, const void *cost_jj);
#endif

#else                           /* HAVE_ANSI */

  void accept_new_state ();
  int generate_new_state ();
#if ASA_PARALLEL
  int generate_new_state_par ();
#endif                          /* ASA_PARALLEL */
  void reanneal ();
  void cost_derivatives ();
  double generate_asa_state ();
  int asa_exit ();
  void Exit_ASA ();
  int asa_test_asa_options ();
  int cost_function_test ();
  void print_string ();
  void print_string_index ();

#if ASA_PRINT
  void print_state ();
  void print_asa_options ();
#endif                          /* ASA_PRINT */

#if TIME_CALC
  void aux_print_time ();
#if TIME_GETRUSAGE
#if TIME_STD
  int syscall ();
#else
  int getrusage ();
#endif                          /* TIME_STD */
#else                           /* TIME_GETRUSAGE */
#if FALSE                       /* change to TRUE for SunOS 4.1.x */
  clock_t clock ();
#endif
#endif                          /* TIME_GETRUSAGE */
#endif                          /* TIME_CALC */

#if MULTI_MIN
  static int multi_compare ();
  double *multi_cost_qsort;
#endif

#if ASA_PARALLEL
  static int sort_parallel ();
#endif

#endif                          /* HAVE_ANSI */
  static int asa_recursive_max = 0;     /* record of max recursions */

#ifdef __cplusplus
}
#endif
#endif                          /* _ASA_H_ */
