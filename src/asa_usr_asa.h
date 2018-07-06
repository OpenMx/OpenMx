#ifndef _ASA_USER_ASA_H_
#define _ASA_USER_ASA_H_
#ifdef __cplusplus
extern "C" {
#endif

/***********************************************************************
* Adaptive Simulated Annealing (ASA)
* Lester Ingber <ingber@ingber.com>
* Copyright (c) 1987-2017 Lester Ingber.  All Rights Reserved.
* ASA-LICENSE file has the license that must be included with ASA code.
***********************************************************************/

  /* $Id: asa_usr_asa.h,v 30.23 2017/11/20 23:17:55 ingber Exp ingber $ */

  /* asa_usr_asa.h for Adaptive Simulated Annealing */

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>             /* misc defs on most machines */
#include <string.h>

/* required if use machine-defined {DBL_EPSILON DBL_MIN DBL_MAX} */
#include <float.h>

#define	TRUE			1
#define	FALSE			0

#define MIN(x,y)	((x) < (y) ? (x) : (y))
#define MAX(x,y)	((x) > (y) ? (x) : (y))

  /* DEFAULT PARAMETERS SETTINGS */

  /* Pre-Compile Options */

  /* Special ASA_TEMPLATEs */

#ifndef MY_TEMPLATE
#define MY_TEMPLATE TRUE
#endif
#if MY_TEMPLATE                 /* MY_TEMPLATE_asa_user */
#define ASA_TEMPLATE_LIB TRUE
#define ASA_FUZZY TRUE
#define OPTIONAL_DATA_PTR TRUE
#define OPTIONAL_PTR_TYPE struct ComputeGenSA
	//#define ASA_PRINT_MORE TRUE
#define MAX_DOUBLE DBL_MAX
#define MIN_DOUBLE DBL_MIN
#define SMALL_FLOAT MIN_DOUBLE
#define EPS_DOUBLE DBL_EPSILON
#define USER_ASA_OUT TRUE
#define INCL_STDOUT FALSE
	//#define ASA_OUT "STDOUT"
#endif                          /* MY_TEMPLATE */

#ifndef ASA_TEMPLATE_LIB
#define ASA_TEMPLATE_LIB FALSE
#endif
#if ASA_TEMPLATE_LIB
#define ASA_LIB TRUE
#define ASA_TEST TRUE
#endif

#ifndef ASA_TEMPLATE_ASA_OUT_PID
#define ASA_TEMPLATE_ASA_OUT_PID FALSE
#endif
#if ASA_TEMPLATE_ASA_OUT_PID
#define USER_ASA_OUT TRUE
#endif

#ifndef ASA_TEMPLATE_MULTIPLE
#define ASA_TEMPLATE_MULTIPLE FALSE
#endif
#if ASA_TEMPLATE_MULTIPLE
#define COST_FILE FALSE
#define USER_ASA_OUT TRUE
#define ASA_TEST TRUE
#define QUENCH_COST TRUE
#define QUENCH_PARAMETERS TRUE
#define OPTIONS_FILE FALSE
#endif

#ifndef ASA_TEMPLATE_SELFOPT
#define ASA_TEMPLATE_SELFOPT FALSE
#endif
#if ASA_TEMPLATE_SELFOPT
#define COST_FILE FALSE
#define SELF_OPTIMIZE TRUE
#define OPTIONAL_DATA_DBL TRUE
#define USER_ASA_OUT TRUE
#define ASA_TEST TRUE
#define OPTIONS_FILE FALSE
#endif

#ifndef ASA_TEMPLATE_SAMPLE
#define ASA_TEMPLATE_SAMPLE FALSE
#endif
#if ASA_TEMPLATE_SAMPLE
#define COST_FILE FALSE
#define ASA_SAMPLE TRUE
#define USER_ACCEPTANCE_TEST TRUE
#define USER_COST_SCHEDULE TRUE
#define OPTIONS_FILE_DATA FALSE
#define USER_ACCEPT_ASYMP_EXP TRUE
#endif

#ifndef ASA_TEMPLATE_PARALLEL
#define ASA_TEMPLATE_PARALLEL FALSE
#endif
#if ASA_TEMPLATE_PARALLEL
#define COST_FILE FALSE
#define ASA_TEST TRUE
#define ASA_PARALLEL TRUE
#endif

#ifndef ASA_TEMPLATE_SAVE
#define ASA_TEMPLATE_SAVE FALSE
#endif
#if ASA_TEMPLATE_SAVE
#define COST_FILE FALSE
#define ASA_TEST TRUE
#define ASA_SAVE TRUE
#define QUENCH_PARAMETERS TRUE
#define QUENCH_COST TRUE
#endif

#ifndef ASA_TEMPLATE_QUEUE
#define ASA_TEMPLATE_QUEUE FALSE
#endif
#if ASA_TEMPLATE_QUEUE
#define ASA_QUEUE TRUE
#define ASA_RESOLUTION FALSE
#define ASA_TEST TRUE
#define COST_FILE FALSE
#define ASA_PRINT TRUE
#define ASA_PRINT_MORE TRUE
#endif

#ifndef ASA_TEST_POINT
#define ASA_TEST_POINT FALSE
#endif
#if ASA_TEST_POINT
#define ASA_TEST TRUE
#define COST_FILE FALSE
#define SMALL_FLOAT 1.0E-50
#define QUENCH_COST TRUE
#endif

#ifndef ASA_EXIT_ANYTIME
#define ASA_EXIT_ANYTIME FALSE
#endif

#ifndef ADAPTIVE_OPTIONS
#define ADAPTIVE_OPTIONS FALSE
#endif

#ifndef ASA_FUZZY
#define ASA_FUZZY FALSE
#endif

#ifndef ASA_FUZZY_PRINT
#define ASA_FUZZY_PRINT FALSE
#endif

#if ASA_FUZZY
  /* defaults */
#if QUENCH_COST
#else
#define QUENCH_COST TRUE
#endif
#if QUENCH_PARAMETERS
#else
#define QUENCH_PARAMETERS TRUE
#endif
#endif

  /* Standard Pre-Compile Options */

#ifndef USER_COST_FUNCTION
#define USER_COST_FUNCTION cost_function
#endif

#if SELF_OPTIMIZE
#ifndef RECUR_USER_COST_FUNCTION
#define RECUR_USER_COST_FUNCTION recur_cost_function
#endif
#ifndef INCL_STDOUT
#define INCL_STDOUT FALSE
#endif
#endif

#ifndef INCL_STDOUT
#define INCL_STDOUT TRUE
#endif
#if INCL_STDOUT
#ifndef TIME_CALC
#define TIME_CALC FALSE
#endif
#endif

#ifndef OPTIONS_FILE
#define OPTIONS_FILE TRUE
#endif

#if OPTIONS_FILE
#ifndef OPTIONS_FILE_DATA
#define OPTIONS_FILE_DATA TRUE
#endif
#else
#define OPTIONS_FILE_DATA FALSE
#endif

#ifndef RECUR_OPTIONS_FILE
#define RECUR_OPTIONS_FILE FALSE
#endif

#if RECUR_OPTIONS_FILE
#ifndef RECUR_OPTIONS_FILE_DATA
#define RECUR_OPTIONS_FILE_DATA FALSE
#endif
#else
#define RECUR_OPTIONS_FILE_DATA FALSE
#endif

#ifndef COST_FILE
#define COST_FILE TRUE
#endif

#ifndef ASA_LIB
#define ASA_LIB FALSE
#endif

#ifndef HAVE_ANSI
#define HAVE_ANSI TRUE
#endif

#ifndef IO_PROTOTYPES
#define IO_PROTOTYPES FALSE
#endif

#ifndef TIME_CALC
#define TIME_CALC FALSE
#endif

#ifndef INT_LONG
#define INT_LONG TRUE
#endif

#if INT_LONG
#define LONG_INT long int
#else
#define LONG_INT int
#endif

#ifndef INT_ALLOC
#define INT_ALLOC FALSE
#endif

#if INT_ALLOC
#define ALLOC_INT int
#else
#define ALLOC_INT LONG_INT
#endif

  /* You can define SMALL_FLOAT to better correlate to your machine's
     precision, i.e., as used in asa */
#ifndef SMALL_FLOAT
#define SMALL_FLOAT 1.0E-18
/* #define SMALL_FLOAT MINFLOAT */
#endif

  /* You can define your machine's maximum and minimum doubles here */
#ifndef MIN_DOUBLE
#define MIN_DOUBLE ((double) SMALL_FLOAT)
/* #define MIN_DOUBLE MINDOUBLE */
#endif

#ifndef MAX_DOUBLE
#define MAX_DOUBLE ((double) 1.0 / (double) SMALL_FLOAT)
/* #define MIN_DOUBLE MINDOUBLE */
#endif

#ifndef EPS_DOUBLE
#define EPS_DOUBLE ((double) SMALL_FLOAT)
/* #define EPS_DOUBLE DBL_EPSILON */
#endif

#ifndef CHECK_EXPONENT
#define CHECK_EXPONENT FALSE
#endif

#ifndef ASA_TEST
#define ASA_TEST FALSE
#endif

#ifndef ASA_TEMPLATE
#define ASA_TEMPLATE FALSE
#endif

#ifndef USER_INITIAL_COST_TEMP
#define USER_INITIAL_COST_TEMP FALSE
#endif

#ifndef RATIO_TEMPERATURE_SCALES
#define RATIO_TEMPERATURE_SCALES FALSE
#endif

#ifndef USER_INITIAL_PARAMETERS_TEMPS
#define USER_INITIAL_PARAMETERS_TEMPS FALSE
#endif

#ifndef DELTA_PARAMETERS
#define DELTA_PARAMETERS FALSE
#endif

#ifndef QUENCH_PARAMETERS
#define QUENCH_PARAMETERS FALSE
#endif

#ifndef QUENCH_COST
#define QUENCH_COST FALSE
#endif

#ifndef QUENCH_PARAMETERS_SCALE
#define QUENCH_PARAMETERS_SCALE TRUE
#endif

#ifndef QUENCH_COST_SCALE
#define QUENCH_COST_SCALE TRUE
#endif

#ifndef OPTIONAL_DATA_DBL
#define OPTIONAL_DATA_DBL FALSE
#endif

#ifndef OPTIONAL_DATA_INT
#define OPTIONAL_DATA_INT FALSE
#endif

#ifndef OPTIONAL_DATA_PTR
#define OPTIONAL_DATA_PTR FALSE
#endif
#if OPTIONAL_DATA_PTR
/* user must define USER_TYPE; if a struct, it must be declared above */
#ifndef OPTIONAL_PTR_TYPE
#define OPTIONAL_PTR_TYPE USER_TYPE
#endif
#endif                          /* OPTIONAL_DATA_PTR */

#ifndef USER_REANNEAL_COST
#define USER_REANNEAL_COST FALSE
#endif

#ifndef USER_REANNEAL_PARAMETERS
#define USER_REANNEAL_PARAMETERS FALSE
#endif

#ifndef MAXIMUM_REANNEAL_INDEX
#define MAXIMUM_REANNEAL_INDEX 50000
#endif

#ifndef REANNEAL_SCALE
#define REANNEAL_SCALE 10
#endif

#ifndef USER_COST_SCHEDULE
#define USER_COST_SCHEDULE FALSE
#endif

#ifndef USER_ACCEPT_ASYMP_EXP
#define USER_ACCEPT_ASYMP_EXP FALSE
#endif

#ifndef USER_ACCEPT_THRESHOLD
#define USER_ACCEPT_THRESHOLD FALSE
#endif

#ifndef USER_ACCEPTANCE_TEST
#define USER_ACCEPTANCE_TEST FALSE
#endif

#ifndef USER_GENERATING_FUNCTION
#define USER_GENERATING_FUNCTION FALSE
#endif

  /* in asa.c, field-width.precision = G_FIELD.G_PRECISION */
#ifndef G_FIELD
#define G_FIELD 12
#endif
#ifndef G_PRECISION
#define G_PRECISION 7
#endif

#define INTEGER_TYPE		((int) 1)
#define REAL_TYPE		((int) -1)
#define INTEGER_NO_REANNEAL	((int) 2)
#define REAL_NO_REANNEAL	((int) -2)

  /* Set this to TRUE to self-optimize the Program Options */
#ifndef SELF_OPTIMIZE
#define SELF_OPTIMIZE FALSE
#endif

#ifndef USER_ASA_OUT
#define USER_ASA_OUT FALSE
#endif

#ifndef USER_ASA_USR_OUT
#define USER_ASA_USR_OUT FALSE
#endif

#ifndef USER_OUT
#define USER_OUT "asa_usr_out"
#endif

#ifndef ASA_SAMPLE
#define ASA_SAMPLE FALSE
#endif

#ifndef ASA_QUEUE
#define ASA_QUEUE FALSE
#endif

#ifndef ASA_RESOLUTION
#define ASA_RESOLUTION FALSE
#endif

#ifndef ASA_PARALLEL
#define ASA_PARALLEL FALSE
#endif

#ifndef ASA_SAVE_OPT
#define ASA_SAVE_OPT FALSE
#endif
#if ASA_SAVE_OPT
#define ASA_SAVE TRUE
#endif

#ifndef ASA_SAVE_BACKUP
#define ASA_SAVE_BACKUP FALSE
#endif
#if ASA_SAVE_BACKUP
#define ASA_SAVE TRUE
#endif

#ifndef ASA_SAVE
#define ASA_SAVE FALSE
#endif

#ifndef ASA_PIPE
#define ASA_PIPE FALSE
#endif

#ifndef ASA_PIPE_FILE
#define ASA_PIPE_FILE FALSE
#endif

#ifndef FDLIBM_POW
#define FDLIBM_POW FALSE
#endif
#if FDLIBM_POW
#define F_POW s_pow
#else
#define F_POW pow
#endif

#ifndef FDLIBM_LOG
#define FDLIBM_LOG FALSE
#endif
#if FDLIBM_LOG
#define F_LOG s_log
#else
#define F_LOG log
#endif

#ifndef FDLIBM_EXP
#define FDLIBM_EXP FALSE
#endif
#if FDLIBM_EXP
#define F_EXP s_exp
#else
#define F_EXP exp
#endif

#ifndef FITLOC
#define FITLOC FALSE
#endif

#ifndef FITLOC_ROUND
#define FITLOC_ROUND TRUE
#endif

#ifndef FITLOC_PRINT
#define FITLOC_PRINT TRUE
#endif

#ifndef MULTI_MIN
#define MULTI_MIN FALSE
#endif

#if ASA_PARALLEL
#ifdef _OPENMP
/* may need specific path */
#include "omp.h"
#endif                          /* _OPENMP */
#endif                          /* ASA_PARALLEL */

  /* Program Options */

  typedef struct {
    LONG_INT Limit_Acceptances;
    LONG_INT Limit_Generated;
    int Limit_Invalid_Generated_States;
    double Accepted_To_Generated_Ratio;

    double Cost_Precision;
    int Maximum_Cost_Repeat;
    int Number_Cost_Samples;
    double Temperature_Ratio_Scale;
    double Cost_Parameter_Scale_Ratio;
    double Temperature_Anneal_Scale;
#if USER_INITIAL_COST_TEMP
    double *User_Cost_Temperature;
#endif

    int Include_Integer_Parameters;
    int User_Initial_Parameters;
    ALLOC_INT Sequential_Parameters;
    double Initial_Parameter_Temperature;
#if RATIO_TEMPERATURE_SCALES
    double *User_Temperature_Ratio;
#endif
#if USER_INITIAL_PARAMETERS_TEMPS
    double *User_Parameter_Temperature;
#endif

    int Acceptance_Frequency_Modulus;
    int Generated_Frequency_Modulus;
    int Reanneal_Cost;
    int Reanneal_Parameters;

    double Delta_X;
#if DELTA_PARAMETERS
    double *User_Delta_Parameter;
#endif
    int User_Tangents;
    int Curvature_0;

#if QUENCH_PARAMETERS
    double *User_Quench_Param_Scale;
#endif
#if QUENCH_COST
    double *User_Quench_Cost_Scale;
#endif

    LONG_INT N_Accepted;
    LONG_INT N_Generated;
    int Locate_Cost;
    int Immediate_Exit;

    double *Best_Cost;
    double *Best_Parameters;
    double *Last_Cost;
    double *Last_Parameters;

#if OPTIONAL_DATA_DBL
    ALLOC_INT Asa_Data_Dim_Dbl;
    double *Asa_Data_Dbl;
#endif
#if OPTIONAL_DATA_INT
    ALLOC_INT Asa_Data_Dim_Int;
    LONG_INT *Asa_Data_Int;
#endif
#if OPTIONAL_DATA_PTR
    ALLOC_INT Asa_Data_Dim_Ptr;
    OPTIONAL_PTR_TYPE *Asa_Data_Ptr;
#endif
#if USER_ASA_OUT
    const char *Asa_Out_File;
#endif
#if USER_ASA_USR_OUT
    char *Asa_Usr_Out_File;
#endif
    /* Keep OPTIONS_TMP in parameter lists in asa_usr.[ch] as they are
     * needed if using recursively, e.g., with SELF_OPTIMIZE=TRUE.
     * Make (USER_DEFINE *) casts explicit within functions. */
#if USER_COST_SCHEDULE
#if HAVE_ANSI
    double (*Cost_Schedule) (double current_cost_temperature,
                             const void *OPTIONS_TMP);
#else                           /* HAVE_ANSI */
    double (*Cost_Schedule) ();
#endif                          /* HAVE_ANSI */
#endif
#if USER_ACCEPT_ASYMP_EXP
    double Asymp_Exp_Param;
#endif
#if USER_ACCEPTANCE_TEST
#if HAVE_ANSI
    void (*Acceptance_Test) (double cost,
                             double *parameter_minimum,
                             double *parameter_maximum,
                             ALLOC_INT * number_parameters,
                             const void *OPTIONS_TMP);
#else                           /* HAVE_ANSI */
    void (*Acceptance_Test) ();
#endif                          /* HAVE_ANSI */
    int User_Acceptance_Flag;
    int Cost_Acceptance_Flag;
    double Cost_Temp_Curr;
    double Cost_Temp_Init;
    double Cost_Temp_Scale;
    double Prob_Bias;
    LONG_INT *Random_Seed;
#endif
#if USER_GENERATING_FUNCTION
#if HAVE_ANSI
    double (*Generating_Distrib) (LONG_INT * seed,
                                  ALLOC_INT * parameter_dimension,
                                  ALLOC_INT index_v,
                                  double temperature_v,
                                  double init_param_temp_v,
                                  double temp_scale_params_v,
                                  double parameter_v,
                                  double parameter_range_v,
                                  double *last_saved_parameter,
                                  const void *OPTIONS_TMP);
#else                           /* HAVE_ANSI */
    double (*Generating_Distrib) ();
#endif                          /* HAVE_ANSI */
#endif
#if USER_REANNEAL_COST
#if HAVE_ANSI
    int (*Reanneal_Cost_Function) (double *cost_best,
                                   double *cost_last,
                                   double *initial_cost_temperature,
                                   double *current_cost_temperature,
                                   const void *OPTIONS_TMP);
#else                           /* HAVE_ANSI */
    int (*Reanneal_Cost_Function) ();
#endif                          /* HAVE_ANSI */
#endif
#if USER_REANNEAL_PARAMETERS
#if HAVE_ANSI
    double (*Reanneal_Params_Function) (double current_temp,
                                        double tangent,
                                        double max_tangent,
                                        const void *OPTIONS_TMP);
#else                           /* HAVE_ANSI */
    double (*Reanneal_Params_Function) ();
#endif                          /* HAVE_ANSI */
#endif
#if ASA_SAMPLE
    double Bias_Acceptance;
    double *Bias_Generated;
    double Average_Weights;
    double Limit_Weights;
#endif
#if ASA_QUEUE
    ALLOC_INT Queue_Size;
    double *Queue_Resolution;
#endif
#if ASA_RESOLUTION
    double *Coarse_Resolution;
#endif
#if FITLOC
    int Fit_Local;
    int Iter_Max;
    double Penalty;
#endif
#if MULTI_MIN
    int Multi_Number;
    double *Multi_Cost;
    double **Multi_Params;
    double *Multi_Grid;
    int Multi_Specify;
#endif
#if ASA_PARALLEL
    int parallel_id;
    int Gener_Mov_Avr;
    LONG_INT Gener_Block;
    LONG_INT Gener_Block_Max;
#endif
#if ASA_SAVE
    ALLOC_INT Random_Array_Dim;
    double *Random_Array;
#endif
    int Asa_Recursive_Level;
#if ASA_FUZZY
    int NoOfSamples;
    double ThresholdDeviation;
    double Threshold1;
    double Performance_Target;
    double Factor_a;
#endif
  } USER_DEFINES;

  /* system function prototypes */

#if HAVE_ANSI

/* This block gives trouble under some Ultrix */
#if FALSE
  int fprintf (FILE * fp, const char *string, ...);
  int sprintf (char *s, const char *format, ...);
  FILE *popen (const char *command, const char *mode);
  void exit (int code);
#endif

#if IO_PROTOTYPES
  int fprintf ();
  int sprintf ();
  int fflush (FILE * fp);
  int fclose (FILE * fp);
  void exit ();
  int fread ();
  int fwrite ();
  int pclose ();
#endif

  double
    asa (double (*user_cost_function)

          
         (double *, double *, double *, double *, double *, ALLOC_INT *,
          int *, int *, int *, USER_DEFINES *),
         double (*user_random_generator) (LONG_INT *), LONG_INT * rand_seed,
         double *parameter_initial_final, double *parameter_minimum,
         double *parameter_maximum, double *tangents, double *curvature,
         ALLOC_INT * number_parameters, int *parameter_type,
         int *valid_state_generated_flag, int *exit_status,
         USER_DEFINES * OPTIONS);

#if TIME_CALC
  void print_time (char *message, FILE * ptr_out);
#endif

#if FDLIBM_POW
  double s_pow (double x, double y);
#endif
#if FDLIBM_LOG
  double s_log (double x);
#endif
#if FDLIBM_EXP
  double s_exp (double x);
#endif

#else                           /* HAVE_ANSI */

#if IO_PROTOTYPES
  int fprintf ();
  int sprintf ();
  int fflush ();
  int fclose ();
  int fread ();
  int fwrite ();
  FILE *popen ();
  int pclose ();
#endif

  double asa ();

#if TIME_CALC
  void print_time ();
#endif

#if FDLIBM_POW
  double s_pow ();
#endif
#if FDLIBM_LOG
  double s_log ();
#endif
#if FDLIBM_EXP
  double s_exp ();
#endif

#endif                          /* HAVE_ANSI */

#ifdef __cplusplus
}
#endif
#endif                          /* _ASA_USER_ASA_H_ */
