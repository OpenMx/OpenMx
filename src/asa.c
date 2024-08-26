/***********************************************************************
* Adaptive Simulated Annealing (ASA)
* Lester Ingber <ingber@ingber.com>
* Copyright (c) 1987-2017 Lester Ingber.  All Rights Reserved.
* ASA-LICENSE file has the license that must be included with ASA code.
***********************************************************************/

#define ASA_ID "/* $Id: asa.c,v 30.23 2017/11/20 23:17:51 ingber Exp ingber $ */"

#include "asa.h"

char exit_msg[160];             /* temp storage for exit messages */

/***********************************************************************
* asa
*       This procedure implements the full ASA function optimization.
***********************************************************************/
#if HAVE_ANSI
double
asa (double (*user_cost_function)

      
     (double *, double *, double *, double *, double *, ALLOC_INT *, int *,
      int *, int *, USER_DEFINES *),
     double (*user_random_generator) (LONG_INT *), LONG_INT * seed,
     double *parameter_initial_final, double *parameter_minimum,
     double *parameter_maximum, double *tangents, double *curvature,
     ALLOC_INT * number_parameters, int *parameter_type,
     int *valid_state_generated_flag, int *exit_status,
     USER_DEFINES * OPTIONS)
#else
double
asa (user_cost_function,
     user_random_generator,
     seed,
     parameter_initial_final,
     parameter_minimum,
     parameter_maximum,
     tangents,
     curvature,
     number_parameters,
     parameter_type, valid_state_generated_flag, exit_status, OPTIONS)
     double (*user_cost_function) ();
     double (*user_random_generator) ();
     LONG_INT *seed;
     double *parameter_initial_final;
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *curvature;
     ALLOC_INT *number_parameters;
     int *parameter_type;
     int *valid_state_generated_flag;
     int *exit_status;
     USER_DEFINES *OPTIONS;
#endif /* HAVE_ANSI */
{
#if USER_REANNEAL_COST
#else
  int immediate_flag;           /* save Immediate_Exit */
#endif /* USER_REANNEAL_COST */
#if USER_INITIAL_COST_TEMP
#if USER_REANNEAL_COST
#else
  int index_cost_constraint;    /* index cost functions averaged */
#endif /* USER_REANNEAL_COST */
#else /* USER_INITIAL_COST_TEMP */
  int index_cost_constraint;    /* index cost functions averaged */
#endif /* USER_INITIAL_COST_TEMP */

  int index_cost_repeat,        /* test OPTIONS->Cost_Precision when =
                                   OPTIONS->Maximum_Cost_Repeat */
    tmp_var_int, tmp_var_int1, tmp_var_int2;    /* temporary integers */

  int generate_flg;
  ALLOC_INT index_v,            /* iteration index */
   *start_sequence;             /* initial OPTIONS->Sequential_Parameters
                                   used if >= 0 */
  double final_cost,            /* best cost to return to user */
    tmp_var_db, tmp_var_db1, tmp_var_db2;       /* temporary doubles */
  int *curvature_flag;
  FILE *ptr_asa_out;            /* file ptr to output file */
  int ret1_flg;

  /* The 3 states that are kept track of during the annealing process */
  STATE *current_generated_state, *last_saved_state, *best_generated_state;

#if ASA_SAVE
  FILE *ptr_save, *ptr_comm;
  int asa_read;
  char asa_save_comm[100];
#if ASA_SAVE_OPT
  char read_option[80];
  char read_if[4], read_FALSE[6], read_comm1[3], read_ASA_SAVE[9],
    read_comm2[3];
  int read_int;
#if INT_LONG
  LONG_INT read_long;
#endif
  double read_double;
  FILE *ptr_save_opt;
#endif
#endif /* ASA_SAVE */

#if ASA_PIPE_FILE
  FILE *ptr_asa_pipe;
#endif

#if ASA_EXIT_ANYTIME
  FILE *ptr_exit_anytime;
#endif /* ASA_EXIT_ANYTIME */

  int asa_exit_value;
  int best_flag;
  int fscanf_ret;

  double xnumber_parameters[1];

  /* The array of tangents (absolute value of the numerical derivatives),
     and the maximum |tangent| of the array */
  double *maximum_tangent;

  /* ratio of acceptances to generated points - determines when to
     test/reanneal */
  double *accepted_to_generated_ratio;

  /* temperature parameters */
  double temperature_scale, *temperature_scale_parameters;
  /* relative scalings of cost and parameters to temperature_scale */
  double *temperature_scale_cost;
  double *current_user_parameter_temp;
  double *initial_user_parameter_temp;
  double *current_cost_temperature;
  double *initial_cost_temperature;
  double log_new_temperature_ratio;     /* current *temp = initial *temp *
                                           exp(log_new_temperature_ratio) */
  ALLOC_INT *index_exit_v;      /* information for asa_exit */

  /* counts of generated states and acceptances */
  LONG_INT *index_parameter_generations;
  LONG_INT *number_generated, *best_number_generated_saved;
  LONG_INT *recent_number_generated, *number_accepted;
  LONG_INT *recent_number_acceptances, *index_cost_acceptances;
  LONG_INT *number_acceptances_saved, *best_number_accepted_saved;

  /* Flag indicates that the parameters generated were
     invalid according to the cost function validity criteria. */
  LONG_INT *number_invalid_generated_states;
  LONG_INT repeated_invalid_states;

#if ASA_QUEUE
  int queue_new;                /* flag to add new entry */
  int *save_queue_flag;         /* save valid_state_generated_flag */
  LONG_INT queue;               /* index of queue */
  LONG_INT queue_v;             /* index of parameters in queue */
  LONG_INT save_queue_test;     /* test if all parameters are present */
  LONG_INT save_queue;          /* last filled position in queue */
  LONG_INT save_queue_indx;     /* current position in queue */
  double *save_queue_cost, *save_queue_param;   /* saved states */
  ALLOC_INT queue_size_tmp;
#endif /* ASA_QUEUE */

#if MULTI_MIN
  int multi_index;
  int multi_test, multi_test_cmp, multi_test_dim;
  int *multi_sort;
  double *multi_cost;
  double **multi_params;
#endif /* MULTI_MIN */

#if ASA_PARALLEL
  int EXIT_asa_parallel = 0;
  LONG_INT tmp_var_lint;
  LONG_INT *parallel_gen_ratio_block;
  LONG_INT *parallel_sort;
  LONG_INT i_prll, sort_index;  /* count of parallel generated states */
  STATE *gener_block_state;
  int *generate_flg_par;
  LONG_INT *number_invalid_generated_states_par;
  LONG_INT *repeated_invalid_states_par;
  double *tmp_var_db1_par;
  double *tmp_var_db_par;
  int *valid_state_generated_flag_par;
  int valid_state_generated_flag_par_test;
#if ASA_QUEUE
  int *queue_new_par;
  LONG_INT *queue_v_par;
  LONG_INT *save_queue_indx_par;
  LONG_INT *save_queue_test_par;
  LONG_INT *save_queue_par;
  double *queue_par_cost;
  int **save_queue_valid_state_flag_par;
  double **save_queue_cost_par;
  double **save_queue_param_par;
#endif /* ASA_QUEUE */
#endif /* ASA_PARALLEL */

  /* used to index repeated and recursive calls to asa */
  /* This assumes that multiple calls (>= 1) u_or_ recursive
     calls are being made to asa */
  static int asa_open = FALSE;
  static int number_asa_open = 0;
  static int recursive_asa_open = 0;

  /* initializations */

  ret1_flg = 0;
  generate_flg = 0;
  if (generate_flg != 0)
    generate_flg = 0;

  fscanf_ret = 0;               /* stop compiler warning */
  if (fscanf_ret) {
    ;
  }

  if ((curvature_flag = (int *) calloc (1, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): curvature_flag");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((maximum_tangent = (double *) calloc (1, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): maximum_tangent");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((accepted_to_generated_ratio =
       (double *) calloc (1, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): accepted_to_generated_ratio");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((temperature_scale_cost =
       (double *) calloc (1, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): temperature_scale_cost");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((current_cost_temperature =
       (double *) calloc (1, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): current_cost_temperature");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((initial_cost_temperature =
       (double *) calloc (1, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): initial_cost_temperature");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((index_exit_v = (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): index_exit_v");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((start_sequence = (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): start_sequence");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((number_generated =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): number_generated");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((best_number_generated_saved =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): best_number_generated_saved");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((recent_number_generated =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): recent_number_generated");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((number_accepted =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): number_accepted");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((recent_number_acceptances =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): recent_number_acceptances");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((index_cost_acceptances =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): index_cost_acceptances");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((number_acceptances_saved =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): number_acceptances_saved");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((best_number_accepted_saved =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): best_number_accepted_saved");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((number_invalid_generated_states =
       (ALLOC_INT *) calloc (1, sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): number_invalid_generated_states");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  if ((current_generated_state =
       (STATE *) calloc (1, sizeof (STATE))) == NULL) {
    strcpy (exit_msg, "asa(): current_generated_state");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((last_saved_state = (STATE *) calloc (1, sizeof (STATE))) == NULL) {
    strcpy (exit_msg, "asa(): last_saved_state");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((best_generated_state = (STATE *) calloc (1, sizeof (STATE))) == NULL) {
    strcpy (exit_msg, "asa(): best_generated_state");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
#if ASA_PARALLEL
  if ((gener_block_state =
       (STATE *) calloc (OPTIONS->Gener_Block_Max, sizeof (STATE))) == NULL) {
    strcpy (exit_msg, "asa(): gener_block_state");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  gener_block_state_qsort = gener_block_state;
  if ((parallel_sort =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): parallel_sort");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  if ((generate_flg_par =
       (int *) calloc (OPTIONS->Gener_Block_Max, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): generate_flg_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((valid_state_generated_flag_par =
       (int *) calloc (OPTIONS->Gener_Block_Max, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): valid_state_generated_flag_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((number_invalid_generated_states_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): number_invalid_generated_states_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((repeated_invalid_states_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): repeated_invalid_states_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((tmp_var_db1_par =
       (double *) calloc (OPTIONS->Gener_Block_Max,
                          sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): tmp_var_db1_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((tmp_var_db_par =
       (double *) calloc (OPTIONS->Gener_Block_Max,
                          sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): tmp_var_db_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  if ((parallel_gen_ratio_block =
       (LONG_INT *) calloc (OPTIONS->Gener_Mov_Avr,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): parallel_gen_ratio_block");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  for (i_prll = 0; i_prll < OPTIONS->Gener_Mov_Avr; ++i_prll) {
    parallel_gen_ratio_block[i_prll] = OPTIONS->Gener_Block;
  }
#endif /* ASA_PARALLEL */

  fscanf_ret = 0;

  /* set default */
  ptr_asa_out = (FILE *) NULL;

  OPTIONS->Immediate_Exit = FALSE;

  if (asa_open == FALSE) {
    asa_open = TRUE;
    ++number_asa_open;
#if ASA_PRINT
    if (number_asa_open == 1) {
      /* open the output file */
#if USER_ASA_OUT
      if (!strcmp (OPTIONS->Asa_Out_File, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
#if ASA_SAVE
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "a");
#else
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "w");
#endif
      }
#else /* USER_ASA_OUT */
      if (!strcmp (ASA_OUT, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
#if ASA_SAVE
        ptr_asa_out = fopen (ASA_OUT, "a");
#else
        ptr_asa_out = fopen (ASA_OUT, "w");
#endif
      }
#endif /* USER_ASA_OUT */
    } else {
#if USER_ASA_OUT
      if (!strcmp (OPTIONS->Asa_Out_File, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "a");
      }
#else
      if (!strcmp (ASA_OUT, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (ASA_OUT, "a");
      }
#endif
      fprintf (ptr_asa_out, "\n\n\t\t number_asa_open = %d\n",
               number_asa_open);
      fflush (ptr_asa_out);
    }
#endif /* ASA_PRINT */
  } else {
    ++recursive_asa_open;
#if ASA_PRINT
    if (recursive_asa_open == 1) {
      /* open the output file */
#if ASA_SAVE
#if USER_ASA_OUT
      if (!strcmp (OPTIONS->Asa_Out_File, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "a");
      }
#else
      if (!strcmp (ASA_OUT, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (ASA_OUT, "a");
      }
#endif
#else /* ASA_SAVE */
#if USER_ASA_OUT
      if (!strcmp (OPTIONS->Asa_Out_File, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "w");
      }
#else
      if (!strcmp (ASA_OUT, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (ASA_OUT, "w");
      }
#endif
#endif /* ASA_SAVE */
    } else {
#if USER_ASA_OUT
      if (!strcmp (OPTIONS->Asa_Out_File, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (OPTIONS->Asa_Out_File, "a");
      }
#else
      if (!strcmp (ASA_OUT, "STDOUT")) {
#if INCL_STDOUT
        ptr_asa_out = stdout;
#endif /* INCL_STDOUT */
      } else {
        ptr_asa_out = fopen (ASA_OUT, "a");
      }
#endif
      fprintf (ptr_asa_out, "\n\n\t\t recursive_asa_open = %d\n",
               recursive_asa_open);
    }
#endif /* ASA_PRINT */
  }

#if ASA_PIPE_FILE
  ptr_asa_pipe = fopen ("asa_pipe", "a");
  fprintf (ptr_asa_pipe, "%s", "%generate");
  fprintf (ptr_asa_pipe, "\t%s", "accept");
  fprintf (ptr_asa_pipe, "\t%s", "best_cost");
  VFOR (index_v)
#if INT_ALLOC
    fprintf (ptr_asa_pipe, "\t%s-%d", "best_param", index_v);
#else
#if INT_LONG
    fprintf (ptr_asa_pipe, "\t%s-%ld", "best_param", index_v);
#else
    fprintf (ptr_asa_pipe, "\t%s-%d", "best_param", index_v);
#endif
#endif
  fprintf (ptr_asa_pipe, "\t%s", "curr_cost");
  VFOR (index_v)
#if INT_ALLOC
    fprintf (ptr_asa_pipe, "\t%s-%d", "curr_param", index_v);
#else
#if INT_LONG
    fprintf (ptr_asa_pipe, "\t%s-%ld", "curr_param", index_v);
#else
    fprintf (ptr_asa_pipe, "\t%s-%d", "curr_param", index_v);
#endif
#endif
  fprintf (ptr_asa_pipe, "\t%s", "cost_temp");
  VFOR (index_v)
#if INT_ALLOC
    fprintf (ptr_asa_pipe, "\t%s-%d", "param_temp", index_v);
#else
#if INT_LONG
    fprintf (ptr_asa_pipe, "\t%s-%ld", "param_temp", index_v);
#else
    fprintf (ptr_asa_pipe, "\t%s-%d", "param_temp", index_v);
#endif
#endif
  fprintf (ptr_asa_pipe, "\t%s", "last_cost");
  fprintf (ptr_asa_pipe, "\n");
  fflush (ptr_asa_pipe);
#endif /* ASA_PIPE_FILE */

#if ASA_EXIT_ANYTIME
  if ((ptr_exit_anytime = fopen ("asa_exit_anytime", "w")) != NULL) {
    fprintf (ptr_exit_anytime, "%s\n",
             "force IMMEDIATE_EXIT by removing this file if ASA_EXIT_ANYTIME is TRUE");
    fflush (ptr_exit_anytime);
    fclose (ptr_exit_anytime);
  }
#endif /* ASA_EXIT_ANYTIME */

#if ASA_PRINT
  /* print header information as defined by user */
  print_asa_options (ptr_asa_out, OPTIONS);

#if TIME_CALC
  /* print starting time */
  print_time ("start_asa", ptr_asa_out);
#endif
  fflush (ptr_asa_out);
#endif /* ASA_PRINT */

  /* set indices and counts to 0 */
  *best_number_generated_saved =
    *number_generated =
    *recent_number_generated = *recent_number_acceptances = 0;
  *index_cost_acceptances =
    *best_number_accepted_saved =
    *number_accepted = *number_acceptances_saved = 0;
  index_cost_repeat = 0;

  OPTIONS->N_Accepted = *number_accepted;
  OPTIONS->N_Generated = *number_generated;

#if ASA_SAMPLE
  OPTIONS->N_Generated = 0;
  OPTIONS->Average_Weights = 1.0;
#endif

  /* do not calculate curvatures initially */
  *curvature_flag = FALSE;

  /* allocate storage for all parameters */
  if ((current_generated_state->parameter =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): current_generated_state->parameter");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  if ((last_saved_state->parameter =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): last_saved_state->parameter");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  if ((best_generated_state->parameter =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): best_generated_state->parameter");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
#if ASA_PARALLEL
  for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
    if ((gener_block_state[i_prll].parameter =
         (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
      strcpy (exit_msg, "asa(): gener_block_state[i_prll].parameter");
      Exit_ASA (exit_msg);
      *exit_status = CALLOC_FAILED;
      ret1_flg = 1;
      goto RET1_asa;
    } else {
      ;
    }
  }
  OPTIONS->parallel_id = -1;
#endif /* ASA_PARALLEL */

  OPTIONS->Best_Cost = &(best_generated_state->cost);
  OPTIONS->Best_Parameters = best_generated_state->parameter;
  OPTIONS->Last_Cost = &(last_saved_state->cost);
  OPTIONS->Last_Parameters = last_saved_state->parameter;

  if ((initial_user_parameter_temp =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): initial_user_parameter_temp");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  if ((index_parameter_generations =
       (ALLOC_INT *) calloc (*number_parameters,
                             sizeof (ALLOC_INT))) == NULL) {
    strcpy (exit_msg, "asa(): index_parameter_generations");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }

  /* set all temperatures */
  if ((current_user_parameter_temp =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): current_user_parameter_temp");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
#if USER_INITIAL_PARAMETERS_TEMPS
  VFOR (index_v)
    current_user_parameter_temp[index_v] =
    initial_user_parameter_temp[index_v] =
    OPTIONS->User_Parameter_Temperature[index_v];
#else
  VFOR (index_v)
    current_user_parameter_temp[index_v] =
    initial_user_parameter_temp[index_v] =
    OPTIONS->Initial_Parameter_Temperature;
#endif

  if ((temperature_scale_parameters =
       (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): temperature_scale_parameters");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
#if ASA_QUEUE
  if (OPTIONS->Queue_Size > 0) {
    queue_size_tmp = OPTIONS->Queue_Size;
  } else {
    queue_size_tmp = 1;
  }
  if ((save_queue_flag =
       (int *) calloc (queue_size_tmp, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_flag");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_cost =
       (double *) calloc (queue_size_tmp, sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_cost");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_param =
       (double *) calloc ((*number_parameters) * queue_size_tmp,
                          sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_param");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
#if ASA_PARALLEL
  if (OPTIONS->Queue_Size > 0) {
    queue_size_tmp = OPTIONS->Queue_Size;
  } else {
    queue_size_tmp = 1;
  }

  if ((queue_par_cost =
       (double *) calloc (OPTIONS->Gener_Block_Max,
                          sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): queue_par_cost");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((queue_new_par =
       (int *) calloc (OPTIONS->Gener_Block_Max, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): queue_new_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((queue_v_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): queue_v_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_indx_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_indx_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_test_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_test_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_par =
       (LONG_INT *) calloc (OPTIONS->Gener_Block_Max,
                            sizeof (LONG_INT))) == NULL) {
    strcpy (exit_msg, "asa(): save_queue_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  if ((save_queue_valid_state_flag_par =
       (int **) calloc (OPTIONS->Gener_Block_Max, sizeof (int *))) == NULL) {
    strcpy (exit_msg, "asa(): *save_queue_valid_state_flag_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_cost_par =
       (double **) calloc (OPTIONS->Gener_Block_Max,
                           sizeof (double *))) == NULL) {
    strcpy (exit_msg, "asa(): *save_queue_cost_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  if ((save_queue_param_par =
       (double **) calloc (OPTIONS->Gener_Block_Max,
                           sizeof (double *))) == NULL) {
    strcpy (exit_msg, "asa(): *save_queue_param_par");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    return (-1);
  }
  for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
    if ((save_queue_valid_state_flag_par[i_prll] =
         (int *) calloc (queue_size_tmp, sizeof (int))) == NULL) {
      strcpy (exit_msg, "asa(): save_queue_valid_state_flag_par[i_prll]");
      Exit_ASA (exit_msg);
      *exit_status = CALLOC_FAILED;
      return (-1);
    }
    if ((save_queue_cost_par[i_prll] =
         (double *) calloc (queue_size_tmp, sizeof (double))) == NULL) {
      strcpy (exit_msg, "asa(): save_queue_cost_par[i_prll]");
      Exit_ASA (exit_msg);
      *exit_status = CALLOC_FAILED;
      return (-1);
    }
    if ((save_queue_param_par[i_prll] =
         (double *) calloc ((*number_parameters) * queue_size_tmp,
                            sizeof (double))) == NULL) {
      strcpy (exit_msg, "asa(): save_queue_param_par[i_prll]");
      Exit_ASA (exit_msg);
      *exit_status = CALLOC_FAILED;
      return (-1);
    }
  }
#endif /* ASA_PARALLEL */
#endif /* ASA_QUEUE */

#if MULTI_MIN
  if ((multi_cost =
       (double *) calloc (OPTIONS->Multi_Number + 1,
                          sizeof (double))) == NULL) {
    strcpy (exit_msg, "asa(): *multi_cost");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  multi_cost_qsort = multi_cost;
  if ((multi_sort =
       (int *) calloc (OPTIONS->Multi_Number + 1, sizeof (int))) == NULL) {
    strcpy (exit_msg, "asa(): *multi_sort");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  if ((multi_params =
       (double **) calloc (OPTIONS->Multi_Number + 1,
                           sizeof (double *))) == NULL) {
    strcpy (exit_msg, "asa(): *multi_params");
    Exit_ASA (exit_msg);
    *exit_status = CALLOC_FAILED;
    ret1_flg = 1;
    goto RET1_asa;
  }
  for (multi_index = 0; multi_index <= OPTIONS->Multi_Number; ++multi_index) {
    if ((multi_params[multi_index] =
         (double *) calloc (*number_parameters, sizeof (double))) == NULL) {
      strcpy (exit_msg, "asa(): multi_params[multi_index]");
      Exit_ASA (exit_msg);
      *exit_status = CALLOC_FAILED;
      ret1_flg = 1;
      goto RET1_asa;
    }
  }
#endif /* MULTI_MIN */
RET1_asa:
  if (ret1_flg == 1) {
#if ASA_PIPE_FILE
    fclose (ptr_asa_pipe);
#endif
    free (accepted_to_generated_ratio);
    free (best_generated_state);
    free (best_number_accepted_saved);
    free (best_number_generated_saved);
    free (current_cost_temperature);
    free (current_generated_state);
    free (curvature_flag);
    free (index_cost_acceptances);
    free (index_exit_v);
    free (initial_cost_temperature);
    free (last_saved_state);
    free (maximum_tangent);
    free (number_acceptances_saved);
    free (number_accepted);
    free (number_generated);
    free (number_invalid_generated_states);
    free (recent_number_acceptances);
    free (recent_number_generated);
    free (start_sequence);
    free (temperature_scale_cost);

    return (-1);
  }
#if USER_INITIAL_COST_TEMP
#if USER_ACCEPTANCE_TEST
  OPTIONS->Cost_Temp_Curr = OPTIONS->Cost_Temp_Init =
#endif
    *initial_cost_temperature = *current_cost_temperature =
    OPTIONS->User_Cost_Temperature[0];
#endif

  /* set parameters to the initial parameter values */
  VFOR (index_v)
    last_saved_state->parameter[index_v] =
    current_generated_state->parameter[index_v] =
    parameter_initial_final[index_v];
#if USER_ACCEPTANCE_TEST
  OPTIONS->Random_Seed = seed;
  OPTIONS->Random_Seed[0] = *seed;
  OPTIONS->User_Acceptance_Flag = TRUE;
  OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif

#if ASA_PRINT
#if INT_LONG
  fprintf (ptr_asa_out, "Initial Random Seed = %ld\n\n", *seed);
#else
  fprintf (ptr_asa_out, "Initial Random Seed = %d\n\n", *seed);
#endif
#endif /* ASA_PRINT */

  /* save initial user value of OPTIONS->Sequential_Parameters */
  *start_sequence = OPTIONS->Sequential_Parameters;

#if ASA_PRINT
#if INT_ALLOC
  fprintf (ptr_asa_out, "*number_parameters = %d\n\n", *number_parameters);
#else
#if INT_LONG
  fprintf (ptr_asa_out, "*number_parameters = %ld\n\n", *number_parameters);
#else
  fprintf (ptr_asa_out, "*number_parameters = %d\n\n", *number_parameters);
#endif
#endif

  /* print the min, max, current values, and types of parameters */
  fprintf (ptr_asa_out,
           "index_v parameter_minimum parameter_maximum parameter_value parameter_type \n");

#if ASA_PRINT_INTERMED
#if INT_ALLOC
  char *out_format = " %-8d %-*.*g \t\t %-*.*g \t %-*.*g %-7d\n";
#else
#if INT_LONG
  char *out_format = " %-8ld %-*.*g \t\t %-*.*g \t %-*.*g %-7d\n";
#else
  char *out_format = " %-8d %-*.*g \t\t %-*.*g \t %-*.*g %-7d\n";
#endif
#endif
  VFOR (index_v) fprintf (ptr_asa_out,
                          out_format,
                          index_v,
                          G_FIELD, G_PRECISION, parameter_minimum[index_v],
                          G_FIELD, G_PRECISION, parameter_maximum[index_v],
                          G_FIELD, G_PRECISION,
                          current_generated_state->parameter[index_v],
                          parameter_type[index_v]);

  fprintf (ptr_asa_out, "\n\n");
#endif /* ASA_PRINT_INTERMED */
  /* Print out user-defined OPTIONS */

#if DELTA_PARAMETERS
#if INT_ALLOC
  char *delta_format = "OPTIONS->User_Delta_Parameter[%d] = %*.*g\n";
#else
#if INT_LONG
  char *delta_format = "OPTIONS->User_Delta_Parameter[%ld] = %*.*g\n";
#else
  char *delta_format = "OPTIONS->User_Delta_Parameter[%d] = %*.*g\n";
#endif
#endif
  VFOR (index_v) fprintf (ptr_asa_out,
                          delta_format,
                          index_v,
                          G_FIELD, G_PRECISION,
                          OPTIONS->User_Delta_Parameter[index_v]);
  fprintf (ptr_asa_out, "\n");
#endif /* DELTA_PARAMETERS */

#if QUENCH_PARAMETERS
#if INT_ALLOC
  char *quench_format = "OPTIONS->User_Quench_Param_Scale[%d] = %*.*g\n";
#else
#if INT_LONG
  char *quench_format = "OPTIONS->User_Quench_Param_Scale[%ld] = %*.*g\n";
#else
  char *quench_format = "OPTIONS->User_Quench_Param_Scale[%d] = %*.*g\n";
#endif
#endif
  VFOR (index_v) fprintf (ptr_asa_out,
                          quench_format,
                          index_v,
                          G_FIELD, G_PRECISION,
                          OPTIONS->User_Quench_Param_Scale[index_v]);
#endif /* QUENCH_PARAMETERS */

#if QUENCH_COST
  fprintf (ptr_asa_out,
           "\nOPTIONS->User_Quench_Cost_Scale = %*.*g\n\n",
           G_FIELD, G_PRECISION, OPTIONS->User_Quench_Cost_Scale[0]);
#endif /* QUENCH_COST */

#if USER_INITIAL_PARAMETERS_TEMPS
#if INT_ALLOC
  char *temps_format = "OPTIONS->User_Parameter_Temperature[%d] = %*.*g\n";
#else
#if INT_LONG
  char *temps_format = "OPTIONS->User_Parameter_Temperature[%ld] = %*.*g\n";
#else
  char *temps_format = "OPTIONS->User_Parameter_Temperature[%d] = %*.*g\n";
#endif
#endif
  VFOR (index_v) fprintf (ptr_asa_out,
                          temps_format,
                          index_v,
                          G_FIELD, G_PRECISION,
                          initial_user_parameter_temp[index_v]);
#endif /* USER_INITIAL_PARAMETERS_TEMPS */

#if RATIO_TEMPERATURE_SCALES
#if INT_ALLOC
  char *scales_format = "OPTIONS->User_Temperature_Ratio[%d] = %*.*g\n";
#else
#if INT_LONG
  char *scales_format = "OPTIONS->User_Temperature_Ratio[%ld] = %*.*g\n";
#else
  char *scales_format = "OPTIONS->User_Temperature_Ratio[%d] = %*.*g\n";
#endif
#endif
  VFOR (index_v) fprintf (ptr_asa_out,
                          scales_format,
                          index_v,
                          G_FIELD, G_PRECISION,
                          OPTIONS->User_Temperature_Ratio[index_v]);
#endif /* RATIO_TEMPERATURE_SCALES */

#if USER_INITIAL_COST_TEMP
  fprintf (ptr_asa_out,
           "OPTIONS->User_Cost_Temperature[0] = %*.*g\n",
           G_FIELD, G_PRECISION, *initial_cost_temperature);
#endif /* USER_INITIAL_COST_TEMP */

  fflush (ptr_asa_out);
#endif /* ASA_PRINT */

#if MULTI_MIN
#if ASA_PRINT
  fprintf (ptr_asa_out, "\n");
  fprintf (ptr_asa_out, "Multi_Number = %d\n", OPTIONS->Multi_Number);
  fprintf (ptr_asa_out, "Multi_Specify = %d\n", OPTIONS->Multi_Specify);
#if ASA_RESOLUTION
#else
  VFOR (index_v) {
    fprintf (ptr_asa_out,
#if INT_ALLOC
             "Multi_Grid[%d] = %*.*g\n",
#else
#if INT_LONG
             "Multi_Grid[%ld] = %*.*g\n",
#else
             "Multi_Grid[%d] = %*.*g\n",
#endif
#endif
             index_v, G_FIELD, G_PRECISION, OPTIONS->Multi_Grid[index_v]);
  }
#endif /* ASA_RESOLUTION */
  fprintf (ptr_asa_out, "\n");
  fflush (ptr_asa_out);
#endif /* ASA_PRINT */
#endif /* MULTI_MIN */

#if ASA_PARALLEL
#if ASA_PRINT
  fprintf (ptr_asa_out,
#if INT_LONG
           "Initial ASA_PARALLEL OPTIONS->\n\t Gener_Block = %ld\n \t Gener_Block_Max = %ld\n \t Gener_Mov_Avr= %d\n\n",
#else
           "ASA_PARALLEL OPTIONS->\n\t Gener_Block = %d\n \t Gener_Block_Max = %d\n \t Gener_Mov_Avr= %d\n\n",
#endif
           OPTIONS->Gener_Block, OPTIONS->Gener_Block_Max,
           OPTIONS->Gener_Mov_Avr);
#endif
#endif /* ASA_PARALLEL */

#if ASA_SAMPLE
#if ASA_PRINT
  fprintf (ptr_asa_out, "OPTIONS->Limit_Weights = %*.*g\n\n",
           G_FIELD, G_PRECISION, OPTIONS->Limit_Weights);
#endif
#endif
  if (OPTIONS->Asa_Recursive_Level > asa_recursive_max)
    asa_recursive_max = OPTIONS->Asa_Recursive_Level;
#if ASA_SAVE
  if (OPTIONS->Asa_Recursive_Level > 0)
    sprintf (asa_save_comm, "asa_save_%d", OPTIONS->Asa_Recursive_Level);
  else
    sprintf (asa_save_comm, "asa_save");
  if ((ptr_save = fopen (asa_save_comm, "r")) == NULL) {
    asa_read = FALSE;
  } else {
#if ASA_PRINT
    fprintf (ptr_asa_out, "\n\n\trestart after ASA_SAVE\n\n");
#endif
    fclose (ptr_save);
    asa_read = TRUE;

    /* give some value to avoid any problems with other OPTIONS */
#if USER_ACCEPTANCE_TEST
    OPTIONS->Cost_Temp_Curr = OPTIONS->Cost_Temp_Init =
#endif
      current_generated_state->cost
      = *initial_cost_temperature = *current_cost_temperature = 3.1416;
  }
#endif

  tmp_var_int = cost_function_test (current_generated_state->cost,
                                    current_generated_state->parameter,
                                    parameter_minimum,
                                    parameter_maximum, number_parameters,
                                    xnumber_parameters);

  /* compute temperature scales */
  tmp_var_db1 = -F_LOG ((OPTIONS->Temperature_Ratio_Scale));
  tmp_var_db2 = F_LOG (OPTIONS->Temperature_Anneal_Scale);
  temperature_scale =
    tmp_var_db1 * F_EXP (-tmp_var_db2 / *xnumber_parameters);

  /* set here in case not used */
  tmp_var_db = ZERO;

#if QUENCH_PARAMETERS
#if RATIO_TEMPERATURE_SCALES
  VFOR (index_v) temperature_scale_parameters[index_v] = tmp_var_db1 * F_EXP
#if QUENCH_PARAMETERS_SCALE
    (-(tmp_var_db2 * OPTIONS->User_Quench_Param_Scale[index_v])
#else
    (-(tmp_var_db2)
#endif
     / *xnumber_parameters)
    * OPTIONS->User_Temperature_Ratio[index_v];
#else
  VFOR (index_v) temperature_scale_parameters[index_v] = tmp_var_db1 * F_EXP
#if QUENCH_PARAMETERS_SCALE
    (-(tmp_var_db2 * OPTIONS->User_Quench_Param_Scale[index_v])
#else
    (-(tmp_var_db2)
#endif
     / *xnumber_parameters);
#endif /* RATIO_TEMPERATURE_SCALES */
#else /* QUENCH_PARAMETERS */
#if RATIO_TEMPERATURE_SCALES
  VFOR (index_v)
    temperature_scale_parameters[index_v] =
    tmp_var_db1 * F_EXP (-(tmp_var_db2) / *xnumber_parameters)
    * OPTIONS->User_Temperature_Ratio[index_v];
#else
  VFOR (index_v)
    temperature_scale_parameters[index_v] =
    tmp_var_db1 * F_EXP (-(tmp_var_db2) / *xnumber_parameters);
#endif /* RATIO_TEMPERATURE_SCALES */
#endif /* QUENCH_PARAMETERS */

#if USER_ACCEPTANCE_TEST
  OPTIONS->Cost_Temp_Scale =
#endif
    *temperature_scale_cost =
#if QUENCH_COST
#if QUENCH_COST_SCALE
    tmp_var_db1 * F_EXP (-(tmp_var_db2 * OPTIONS->User_Quench_Cost_Scale[0])
#else
    tmp_var_db1 * F_EXP (-(tmp_var_db2)
#endif
                         / *xnumber_parameters) *
    OPTIONS->Cost_Parameter_Scale_Ratio;
#else /* QUENCH_COST */
    tmp_var_db1 * F_EXP (-(tmp_var_db2)
                         / *xnumber_parameters) *
    OPTIONS->Cost_Parameter_Scale_Ratio;
#endif /* QUENCH_COST */

  /* set the initial index of parameter generations to 1 */
  VFOR (index_v) index_parameter_generations[index_v] = 1;

  /* test user-defined options before calling cost function */
  tmp_var_int = asa_test_asa_options (seed,
                                      parameter_initial_final,
                                      parameter_minimum,
                                      parameter_maximum,
                                      tangents,
                                      curvature,
                                      number_parameters,
                                      parameter_type,
                                      valid_state_generated_flag,
                                      exit_status, ptr_asa_out, OPTIONS);
  if (tmp_var_int > 0) {
#if ASA_PRINT
    fprintf (ptr_asa_out, "total number invalid OPTIONS = %d\n", tmp_var_int);
    fflush (ptr_asa_out);
#endif
    *exit_status = INVALID_USER_INPUT;
    goto EXIT_asa;
  }
#if USER_INITIAL_COST_TEMP
#else
#if ASA_SAVE
  if (asa_read == TRUE)
    OPTIONS->Number_Cost_Samples = 1;
#endif
  /* calculate the average cost over samplings of the cost function */
  if (OPTIONS->Number_Cost_Samples < -1) {
    tmp_var_db1 = ZERO;
    tmp_var_db2 = ZERO;
    tmp_var_int = -OPTIONS->Number_Cost_Samples;
  } else {
    tmp_var_db1 = ZERO;
    tmp_var_int = OPTIONS->Number_Cost_Samples;
  }

  OPTIONS->Locate_Cost = 0;     /* initial cost temp */

  for (index_cost_constraint = 0;
       index_cost_constraint < tmp_var_int; ++index_cost_constraint) {
    *number_invalid_generated_states = 0;
    repeated_invalid_states = 0;
    OPTIONS->Sequential_Parameters = *start_sequence - 1;
    do {
#if ASA_EXIT_ANYTIME
      if ((ptr_exit_anytime = fopen ("asa_exit_anytime", "r")) == NULL) {
        *exit_status = IMMEDIATE_EXIT;
        goto EXIT_asa;
      } else {
        fclose (ptr_exit_anytime);
      }
#endif /* ASA_EXIT_ANYTIME */
      ++(*number_invalid_generated_states);
      generate_flg = generate_new_state (user_random_generator,
                                         seed,
                                         parameter_minimum,
                                         parameter_maximum,
                                         current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                                         initial_user_parameter_temp,
                                         temperature_scale_parameters,
#endif
                                         number_parameters,
                                         parameter_type,
                                         current_generated_state,
                                         last_saved_state, OPTIONS);
      *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
      OPTIONS->User_Acceptance_Flag = TRUE;
      OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
      tmp_var_db =
        user_cost_function (current_generated_state->parameter,
                            parameter_minimum,
                            parameter_maximum,
                            tangents,
                            curvature,
                            number_parameters,
                            parameter_type,
                            valid_state_generated_flag, exit_status, OPTIONS);
      if (cost_function_test
          (tmp_var_db, current_generated_state->parameter,
           parameter_minimum, parameter_maximum, number_parameters,
           xnumber_parameters) == 0) {
        *exit_status = INVALID_COST_FUNCTION;
        goto EXIT_asa;
      }

      ++repeated_invalid_states;
      if (repeated_invalid_states > OPTIONS->Limit_Invalid_Generated_States) {
        *exit_status = TOO_MANY_INVALID_STATES;
        goto EXIT_asa;
      }
    }
    while (*valid_state_generated_flag == FALSE);
    --(*number_invalid_generated_states);

    if (OPTIONS->Number_Cost_Samples < -1) {
      tmp_var_db1 += tmp_var_db;
      tmp_var_db2 += (tmp_var_db * tmp_var_db);
    } else {
      tmp_var_db1 += fabs (tmp_var_db);
    }
  }
  if (OPTIONS->Number_Cost_Samples < -1) {
    tmp_var_db1 /= (double) tmp_var_int;
    tmp_var_db2 /= (double) tmp_var_int;
    tmp_var_db = sqrt (fabs ((tmp_var_db2 - tmp_var_db1 * tmp_var_db1)
                             * ((double) tmp_var_int
                                / ((double) tmp_var_int - ONE))))
      + (double) EPS_DOUBLE;
  } else {
    tmp_var_db = tmp_var_db1 / (double) tmp_var_int;
  }

#if USER_ACCEPTANCE_TEST
  OPTIONS->Cost_Temp_Curr = OPTIONS->Cost_Temp_Init =
#endif
    *initial_cost_temperature = *current_cost_temperature = tmp_var_db;
  if (fabs (*initial_cost_temperature) <= SMALL_FLOAT) {
    *initial_cost_temperature = *current_cost_temperature = 2.718;
#if ASA_PRINT
    fprintf (ptr_asa_out,
             "*** invalid too small cost temp = %g, set to = %g ***\n",
             tmp_var_db, *initial_cost_temperature);
    fflush (ptr_asa_out);
#endif
  }
#endif /* USER_INITIAL_COST_TEMP */

  /* set all parameters to the initial parameter values */
  VFOR (index_v)
    best_generated_state->parameter[index_v] =
    last_saved_state->parameter[index_v] =
    current_generated_state->parameter[index_v] =
    parameter_initial_final[index_v];

  OPTIONS->Locate_Cost = 1;     /* initial cost value */

  /* if using user's initial parameters */
  if (OPTIONS->User_Initial_Parameters == TRUE) {
    *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
    OPTIONS->User_Acceptance_Flag = TRUE;
    OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
#if ASA_SAVE
    if (asa_read == FALSE)
#endif
      current_generated_state->cost =
        user_cost_function (current_generated_state->parameter,
                            parameter_minimum,
                            parameter_maximum,
                            tangents,
                            curvature,
                            number_parameters,
                            parameter_type,
                            valid_state_generated_flag, exit_status, OPTIONS);
    if (cost_function_test
        (current_generated_state->cost, current_generated_state->parameter,
         parameter_minimum, parameter_maximum, number_parameters,
         xnumber_parameters) == 0) {
      *exit_status = INVALID_COST_FUNCTION;
      goto EXIT_asa;
    }
#if ASA_PRINT
    if (*valid_state_generated_flag == FALSE)
      fprintf (ptr_asa_out,
               "user's initial parameters generated FALSE *valid_state_generated_flag\n");
#endif
  } else {
    /* let asa generate valid initial parameters */
    repeated_invalid_states = 0;
    OPTIONS->Sequential_Parameters = *start_sequence - 1;
    do {
#if ASA_EXIT_ANYTIME
      if ((ptr_exit_anytime = fopen ("asa_exit_anytime", "r")) == NULL) {
        *exit_status = IMMEDIATE_EXIT;
        goto EXIT_asa;
      } else {
        fclose (ptr_exit_anytime);
      }
#endif /* ASA_EXIT_ANYTIME */
      ++(*number_invalid_generated_states);
      generate_flg = generate_new_state (user_random_generator,
                                         seed,
                                         parameter_minimum,
                                         parameter_maximum,
                                         current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                                         initial_user_parameter_temp,
                                         temperature_scale_parameters,
#endif
                                         number_parameters,
                                         parameter_type,
                                         current_generated_state,
                                         last_saved_state, OPTIONS);
      *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
      OPTIONS->User_Acceptance_Flag = TRUE;
      OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
      current_generated_state->cost =
        user_cost_function (current_generated_state->parameter,
                            parameter_minimum,
                            parameter_maximum,
                            tangents,
                            curvature,
                            number_parameters,
                            parameter_type,
                            valid_state_generated_flag, exit_status, OPTIONS);
      if (cost_function_test
          (current_generated_state->cost,
           current_generated_state->parameter, parameter_minimum,
           parameter_maximum, number_parameters, xnumber_parameters) == 0) {
        *exit_status = INVALID_COST_FUNCTION;
        goto EXIT_asa;
      }
      ++repeated_invalid_states;
      if (repeated_invalid_states > OPTIONS->Limit_Invalid_Generated_States) {
        *exit_status = TOO_MANY_INVALID_STATES;
        goto EXIT_asa;
      }
    }
    while (*valid_state_generated_flag == FALSE);
    --(*number_invalid_generated_states);
  }                             /* OPTIONS->User_Initial_Parameters */

  /* set all states to the last one generated */
  VFOR (index_v) {
#if DROPPED_PARAMETERS
    /* ignore parameters that have too small a range */
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
#endif
    best_generated_state->parameter[index_v] =
      last_saved_state->parameter[index_v] =
      current_generated_state->parameter[index_v];
  }

  /* set all costs to the last one generated */
  best_generated_state->cost = last_saved_state->cost =
    current_generated_state->cost;

  *accepted_to_generated_ratio = ONE;

  /* do not calculate curvatures initially */
  *curvature_flag = FALSE;

#if ASA_PRINT
  fprintf (ptr_asa_out,
           "temperature_scale = %*.*g\n",
           G_FIELD, G_PRECISION, temperature_scale);
#if RATIO_TEMPERATURE_SCALES
#if ASA_PRINT_INTERMED
  VFOR (index_v) {
    fprintf (ptr_asa_out,
#if INT_ALLOC
             "temperature_scale_parameters[%d] = %*.*g\n",
#else
#if INT_LONG
             "temperature_scale_parameters[%ld] = %*.*g\n",
#else
             "temperature_scale_parameters[%d] = %*.*g\n",
#endif
#endif
             index_v,
             G_FIELD, G_PRECISION, temperature_scale_parameters[index_v]);
  }
#endif
#else
  fprintf (ptr_asa_out,
           "temperature_scale_parameters[0] = %*.*g\n",
           G_FIELD, G_PRECISION, temperature_scale_parameters[0]);
#endif /* RATIO_TEMPERATURE_SCALES */
  fprintf (ptr_asa_out,
           "*temperature_scale_cost = %*.*g\n",
           G_FIELD, G_PRECISION, *temperature_scale_cost);
  fprintf (ptr_asa_out, "\n\n");

#if ASA_PRINT_INTERMED
  print_state (parameter_minimum,
               parameter_maximum,
               tangents,
               curvature,
               current_cost_temperature,
               current_user_parameter_temp,
               accepted_to_generated_ratio,
               number_parameters,
               curvature_flag,
               number_accepted,
               index_cost_acceptances,
               number_generated,
               number_invalid_generated_states,
               last_saved_state, best_generated_state, ptr_asa_out, OPTIONS);
#endif
  fprintf (ptr_asa_out, "\n");

  fflush (ptr_asa_out);
#endif

#if ASA_SAMPLE
#if ASA_PRINT
  fprintf (ptr_asa_out,
           ":SAMPLE:   n_accept   cost        cost_temp    bias_accept    aver_weight\n");
  fprintf (ptr_asa_out,
           ":SAMPLE:   index      param[]     temp[]       bias_gener[]   range[]\n");
#endif
#endif

  /* reset the current cost and the number of generations performed */
  *number_invalid_generated_states = 0;
  *best_number_generated_saved =
    *number_generated = *recent_number_generated = 0;
  OPTIONS->N_Generated = *number_generated;
  VFOR (index_v) {
    /* ignore parameters that have too small a range */
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
    index_parameter_generations[index_v] = 1;
  }
#if USER_ACCEPTANCE_TEST
  OPTIONS->User_Acceptance_Flag = FALSE;
  OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif

#if ASA_QUEUE
#if ASA_PRINT
#if INT_ALLOC
  fprintf (ptr_asa_out, "OPTIONS->Queue_Size = %d\n", OPTIONS->Queue_Size);
#else
#if INT_LONG
  fprintf (ptr_asa_out, "OPTIONS->Queue_Size = %ld\n", OPTIONS->Queue_Size);
#else
  fprintf (ptr_asa_out, "OPTIONS->Queue_Size = %d\n", OPTIONS->Queue_Size);
#endif
#endif
  VFOR (index_v) {
    fprintf (ptr_asa_out,
#if INT_ALLOC
             "Queue_Resolution[%d] = %*.*g\n",
#else
#if INT_LONG
             "Queue_Resolution[%ld] = %*.*g\n",
#else
             "Queue_Resolution[%d] = %*.*g\n",
#endif
#endif
             index_v,
             G_FIELD, G_PRECISION, OPTIONS->Queue_Resolution[index_v]);
  }
#endif /* ASA_PRINT */

  /* fill arrays to check allocated memory */
  for (queue = 0; queue < (LONG_INT) queue_size_tmp; ++queue) {
    VFOR (index_v) {
      if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
        continue;
      }
      queue_v = index_v + queue * (LONG_INT) (*number_parameters);
      save_queue_param[queue_v] = current_generated_state->parameter[index_v];
    }
    save_queue_cost[queue] = current_generated_state->cost;
    save_queue_flag[queue] = *valid_state_generated_flag;
  }
  save_queue = save_queue_indx = 0;
#if ASA_PARALLEL
  for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
    for (queue = 0; queue < (LONG_INT) queue_size_tmp; ++queue) {
      VFOR (index_v) {
        if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
          continue;
        }
        queue_v_par[i_prll] =
          index_v + queue * (LONG_INT) (*number_parameters);
        save_queue_param_par[i_prll][queue_v_par[i_prll]] =
          current_generated_state->parameter[index_v];
      }
      save_queue_cost_par[i_prll][queue] = current_generated_state->cost;
      save_queue_valid_state_flag_par[i_prll][queue] =
        *valid_state_generated_flag;
    }
    save_queue_par[i_prll] = save_queue_indx_par[i_prll] = 0;
  }
#endif /* ASA_PARALLEL */
#endif /* ASA_QUEUE */

#if ASA_RESOLUTION
#if ASA_PRINT
  VFOR (index_v) {
    fprintf (ptr_asa_out,
#if INT_ALLOC
             "Coarse_Resolution[%d] = %*.*g\n",
#else
#if INT_LONG
             "Coarse_Resolution[%ld] = %*.*g\n",
#else
             "Coarse_Resolution[%d] = %*.*g\n",
#endif
#endif
             index_v,
             G_FIELD, G_PRECISION, OPTIONS->Coarse_Resolution[index_v]);
  }
#endif /* ASA_PRINT */
#endif /* ASA_RESOLUTION */

#if MULTI_MIN
  multi_sort[OPTIONS->Multi_Number] = OPTIONS->Multi_Number;
  multi_cost[OPTIONS->Multi_Number] = current_generated_state->cost;
  VFOR (index_v) {
    multi_params[OPTIONS->Multi_Number][index_v] =
      current_generated_state->parameter[index_v];
  }
  for (multi_index = 0; multi_index < OPTIONS->Multi_Number; ++multi_index) {
    multi_sort[multi_index] = multi_index;
    multi_cost[multi_index] = OPTIONS->Multi_Cost[multi_index] =
      current_generated_state->cost;
    VFOR (index_v) {
      multi_params[multi_index][index_v] =
        OPTIONS->Multi_Params[multi_index][index_v] =
        current_generated_state->parameter[index_v];
    }
  }
#endif /* MULTI_MIN */

  /* this test is after MULTI_MIN so that params are not all just set to 0 */
  if (*initial_cost_temperature < (double) EPS_DOUBLE) {
#if ASA_PRINT
    fprintf (ptr_asa_out, "*initial_cost_temperature (= %g) < EPS_DOUBLE\n",
             *initial_cost_temperature);
    fflush (ptr_asa_out);
#endif
    *exit_status = INVALID_COST_FUNCTION;
    goto EXIT_asa;
  }

  OPTIONS->Sequential_Parameters = *start_sequence - 1;

  /* MAIN ANNEALING LOOP */
  while (((*number_accepted <= OPTIONS->Limit_Acceptances)
          || (OPTIONS->Limit_Acceptances == 0))
         && ((*number_generated <= OPTIONS->Limit_Generated)
             || (OPTIONS->Limit_Generated == 0))) {

    tmp_var_db1 = -F_LOG ((OPTIONS->Temperature_Ratio_Scale));

    /* compute temperature scales */
    tmp_var_db2 = F_LOG (OPTIONS->Temperature_Anneal_Scale);
    temperature_scale = tmp_var_db1 *
      F_EXP (-tmp_var_db2 / *xnumber_parameters);

#if QUENCH_PARAMETERS
#if RATIO_TEMPERATURE_SCALES
    VFOR (index_v)
      temperature_scale_parameters[index_v] = tmp_var_db1 * F_EXP
#if QUENCH_PARAMETERS_SCALE
      (-(tmp_var_db2 * OPTIONS->User_Quench_Param_Scale[index_v])
#else
      (-(tmp_var_db2)
#endif
       / *xnumber_parameters)
      * OPTIONS->User_Temperature_Ratio[index_v];
#else
    VFOR (index_v)
      temperature_scale_parameters[index_v] = tmp_var_db1 * F_EXP
#if QUENCH_PARAMETERS_SCALE
      (-(tmp_var_db2 * OPTIONS->User_Quench_Param_Scale[index_v])
#else
      (-(tmp_var_db2)
#endif
       / *xnumber_parameters);
#endif /* RATIO_TEMPERATURE_SCALES */
#else /* QUENCH_PARAMETERS */
#if RATIO_TEMPERATURE_SCALES
    VFOR (index_v)
      temperature_scale_parameters[index_v] =
      tmp_var_db1 * F_EXP (-(tmp_var_db2) / *xnumber_parameters)
      * OPTIONS->User_Temperature_Ratio[index_v];
#else
    VFOR (index_v)
      temperature_scale_parameters[index_v] =
      tmp_var_db1 * F_EXP (-(tmp_var_db2) / *xnumber_parameters);
#endif /* RATIO_TEMPERATURE_SCALES */
#endif /* QUENCH_PARAMETERS */

#if USER_ACCEPTANCE_TEST
    OPTIONS->Cost_Temp_Scale =
#endif
      *temperature_scale_cost =
#if QUENCH_COST
#if QUENCH_COST_SCALE
      tmp_var_db1 * F_EXP (-(tmp_var_db2 * OPTIONS->User_Quench_Cost_Scale[0])
#else
      tmp_var_db1 * F_EXP (-(tmp_var_db2)
#endif
                           / *xnumber_parameters) *
      OPTIONS->Cost_Parameter_Scale_Ratio;
#else /* QUENCH_COST */
      tmp_var_db1 * F_EXP (-(tmp_var_db2)
                           / *xnumber_parameters) *
      OPTIONS->Cost_Parameter_Scale_Ratio;
#endif /* QUENCH_COST */

    /* CALCULATE NEW TEMPERATURES */

    /* calculate new parameter temperatures */
    VFOR (index_v) {
      /* skip parameters with too small a range */
      if (PARAMETER_RANGE_TOO_SMALL (index_v))
        continue;

      log_new_temperature_ratio =
        -temperature_scale_parameters[index_v] *
        F_POW ((double) index_parameter_generations[index_v],
#if QUENCH_PARAMETERS
               OPTIONS->User_Quench_Param_Scale[index_v]
#else /* QUENCH_PARAMETERS */
               ONE
#endif /* QUENCH_PARAMETERS */
               / *xnumber_parameters);
      /* check (and correct) for too large an exponent */
      log_new_temperature_ratio = EXPONENT_CHECK (log_new_temperature_ratio);
      current_user_parameter_temp[index_v] =
        initial_user_parameter_temp[index_v]
        * F_EXP (log_new_temperature_ratio);

#if NO_PARAM_TEMP_TEST
      if (current_user_parameter_temp[index_v] < (double) EPS_DOUBLE)
        current_user_parameter_temp[index_v] = (double) EPS_DOUBLE;
#else
      /* check for too small a parameter temperature */
      if (current_user_parameter_temp[index_v] < (double) EPS_DOUBLE) {
        *exit_status = P_TEMP_TOO_SMALL;
        *index_exit_v = index_v;
        goto EXIT_asa;
      }
#endif
    }

    /* calculate new cost temperature */
    log_new_temperature_ratio =
      -*temperature_scale_cost * F_POW ((double) *index_cost_acceptances,
#if QUENCH_COST
                                        OPTIONS->User_Quench_Cost_Scale[0]
#else
                                        ONE
#endif
                                        / *xnumber_parameters);
    log_new_temperature_ratio = EXPONENT_CHECK (log_new_temperature_ratio);
#if USER_ACCEPTANCE_TEST
    OPTIONS->Cost_Temp_Curr = OPTIONS->Cost_Temp_Init =
#endif
      *current_cost_temperature = *initial_cost_temperature
      * F_EXP (log_new_temperature_ratio);

#if NO_COST_TEMP_TEST
    if (*current_cost_temperature < (double) EPS_DOUBLE)
#if USER_ACCEPTANCE_TEST
      OPTIONS->Cost_Temp_Curr =
#endif
        *current_cost_temperature = (double) EPS_DOUBLE;
#else
    /* check for too small a cost temperature */
    if (*current_cost_temperature < (double) EPS_DOUBLE) {
      *exit_status = C_TEMP_TOO_SMALL;
      goto EXIT_asa;
    }
#endif

#if ASA_SAVE
    if (asa_read == TRUE && OPTIONS->Asa_Recursive_Level == asa_recursive_max) {
      if (OPTIONS->Asa_Recursive_Level > 0)
        sprintf (asa_save_comm, "asa_save_%d", OPTIONS->Asa_Recursive_Level);
      else
        sprintf (asa_save_comm, "asa_save");
      ptr_save = fopen (asa_save_comm, "r");

      fread (number_parameters, sizeof (ALLOC_INT), 1, ptr_save);
      fread (xnumber_parameters, sizeof (double), 1, ptr_save);
      fread (parameter_minimum, sizeof (double),
             *number_parameters, ptr_save);
      fread (parameter_maximum, sizeof (double),
             *number_parameters, ptr_save);
      fread (tangents, sizeof (double), *number_parameters, ptr_save);
      fread (current_user_parameter_temp, sizeof (double),
             *number_parameters, ptr_save);
      fread (initial_user_parameter_temp, sizeof (double),
             *number_parameters, ptr_save);
      fread (temperature_scale_parameters, sizeof (double),
             *number_parameters, ptr_save);

      fread (parameter_type, sizeof (int), *number_parameters, ptr_save);
      fread (&index_cost_repeat, sizeof (int), 1, ptr_save);
      fread (&asa_open, sizeof (int), 1, ptr_save);
      fread (&number_asa_open, sizeof (int), 1, ptr_save);
      fread (&recursive_asa_open, sizeof (int), 1, ptr_save);

      fread (current_cost_temperature, sizeof (double), 1, ptr_save);
      fread (initial_cost_temperature, sizeof (double), 1, ptr_save);
      fread (temperature_scale_cost, sizeof (double), 1, ptr_save);
      fread (accepted_to_generated_ratio, sizeof (double), 1, ptr_save);

      fread (curvature_flag, sizeof (int), 1, ptr_save);

      fread (seed, sizeof (LONG_INT), 1, ptr_save);
      fread (number_generated, sizeof (LONG_INT), 1, ptr_save);
      fread (number_accepted, sizeof (LONG_INT), 1, ptr_save);
      fread (number_acceptances_saved, sizeof (LONG_INT), 1, ptr_save);
      fread (recent_number_acceptances, sizeof (LONG_INT), 1, ptr_save);
      fread (recent_number_generated, sizeof (LONG_INT), 1, ptr_save);
      fread (number_invalid_generated_states, sizeof (LONG_INT), 1, ptr_save);
      fread (index_cost_acceptances, sizeof (LONG_INT), 1, ptr_save);
      fread (best_number_generated_saved, sizeof (LONG_INT), 1, ptr_save);
      fread (best_number_accepted_saved, sizeof (LONG_INT), 1, ptr_save);

      fread (index_parameter_generations, sizeof (LONG_INT),
             *number_parameters, ptr_save);

      fread (current_generated_state->parameter,
             sizeof (double), *number_parameters, ptr_save);
      fread (last_saved_state->parameter,
             sizeof (double), *number_parameters, ptr_save);
      fread (best_generated_state->parameter,
             sizeof (double), *number_parameters, ptr_save);
      fread (&(current_generated_state->cost), sizeof (double), 1, ptr_save);
      fread (&(last_saved_state->cost), sizeof (double), 1, ptr_save);
      fread (&(best_generated_state->cost), sizeof (double), 1, ptr_save);

      fread (&(OPTIONS->Limit_Acceptances), sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->Limit_Generated), sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->Limit_Invalid_Generated_States), sizeof (int),
             1, ptr_save);
      fread (&(OPTIONS->Accepted_To_Generated_Ratio), sizeof (double),
             1, ptr_save);
      fread (&(OPTIONS->Cost_Precision), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Maximum_Cost_Repeat), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Number_Cost_Samples), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Temperature_Ratio_Scale), sizeof (double),
             1, ptr_save);
      fread (&(OPTIONS->Cost_Parameter_Scale_Ratio), sizeof (double),
             1, ptr_save);
      fread (&(OPTIONS->Temperature_Anneal_Scale), sizeof (double),
             1, ptr_save);
      fread (&(OPTIONS->Include_Integer_Parameters), sizeof (int),
             1, ptr_save);
      fread (&(OPTIONS->User_Initial_Parameters), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Sequential_Parameters), sizeof (ALLOC_INT), 1,
             ptr_save);
      fread (&(OPTIONS->Initial_Parameter_Temperature), sizeof (double), 1,
             ptr_save);
      fread (&(OPTIONS->Acceptance_Frequency_Modulus), sizeof (int), 1,
             ptr_save);
      fread (&(OPTIONS->Generated_Frequency_Modulus), sizeof (int), 1,
             ptr_save);
      fread (&(OPTIONS->Reanneal_Cost), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Reanneal_Parameters), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Delta_X), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->User_Tangents), sizeof (int), 1, ptr_save);

#if USER_INITIAL_COST_TEMP
      fread (&(OPTIONS->User_Cost_Temperature), sizeof (double), 1, ptr_save);
#endif
#if RATIO_TEMPERATURE_SCALES
      fread (OPTIONS->User_Temperature_Ratio, sizeof (double),
             *number_parameters, ptr_save);
#endif
#if USER_INITIAL_PARAMETERS_TEMPS
      fread (OPTIONS->User_Parameter_Temperature, sizeof (double),
             *number_parameters, ptr_save);
#endif
#if DELTA_PARAMETERS
      fread (OPTIONS->User_Delta_Parameter, sizeof (double),
             *number_parameters, ptr_save);
#endif
#if QUENCH_PARAMETERS
      fread (OPTIONS->User_Quench_Param_Scale, sizeof (double),
             *number_parameters, ptr_save);
#endif
#if QUENCH_COST
      fread (OPTIONS->User_Quench_Cost_Scale, sizeof (double), 1, ptr_save);
#endif
      fread (&(OPTIONS->N_Accepted), sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->N_Generated), sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->Locate_Cost), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Immediate_Exit), sizeof (int), 1, ptr_save);
#if OPTIONAL_DATA_DBL
      fread (&(OPTIONS->Asa_Data_Dim_Dbl), sizeof (ALLOC_INT), 1, ptr_save);
      fread (OPTIONS->Asa_Data_Dbl, sizeof (double),
             OPTIONS->Asa_Data_Dim_Dbl, ptr_save);
#endif
      fread (&(OPTIONS->Random_Array_Dim), sizeof (ALLOC_INT), 1, ptr_save);
      fread (OPTIONS->Random_Array, sizeof (double),
             OPTIONS->Random_Array_Dim, ptr_save);
      fread (&(OPTIONS->Asa_Recursive_Level), sizeof (int), 1, ptr_save);
#if OPTIONAL_DATA_INT
      fread (&(OPTIONS->Asa_Data_Dim_Int), sizeof (ALLOC_INT), 1, ptr_save);
      fread (OPTIONS->Asa_Data_Int, sizeof (LONG_INT),
             OPTIONS->Asa_Data_Dim_Int, ptr_save);
#endif
#if OPTIONAL_DATA_PTR
      fread (&(OPTIONS->Asa_Data_Dim_Ptr), sizeof (ALLOC_INT), 1, ptr_save);
      if (OPTIONS->Asa_Recursive_Level == 0)
        fread (OPTIONS->Asa_Data_Ptr, sizeof (OPTIONAL_PTR_TYPE),
               OPTIONS->Asa_Data_Dim_Ptr, ptr_save);
#if ASA_TEMPLATE_SELFOPT
      if (OPTIONS->Asa_Recursive_Level == 1)
        fread (OPTIONS->Asa_Data_Ptr, sizeof (RECUR_OPTIONAL_PTR_TYPE),
               OPTIONS->Asa_Data_Dim_Ptr, ptr_save);
#endif
#endif
#if USER_ASA_OUT
      fread (OPTIONS->Asa_Out_File, sizeof (char), 80, ptr_save);
#endif
#if USER_ASA_USR_OUT
      fread (OPTIONS->Asa_Usr_Out_File, sizeof (char), 80, ptr_save);
#endif
#if USER_COST_SCHEDULE
      fread (&(OPTIONS->Cost_Schedule), sizeof (char), 1, ptr_save);
#endif
#if USER_ACCEPT_ASYMP_EXP
      fread (&(OPTIONS->Asymp_Exp_Param), sizeof (double), 1, ptr_save);
#endif
#if USER_ACCEPTANCE_TEST
      fread (&(OPTIONS->Acceptance_Test), sizeof (char), 1, ptr_save);
      fread (&(OPTIONS->User_Acceptance_Flag), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Cost_Acceptance_Flag), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Cost_Temp_Curr), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Cost_Temp_Init), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Cost_Temp_Scale), sizeof (double), 1, ptr_save);
#endif
#if USER_GENERATING_FUNCTION
      fread (&(OPTIONS->Generating_Distrib), sizeof (char), 1, ptr_save);
#endif
#if USER_REANNEAL_COST
      fread (&(OPTIONS->Reanneal_Cost_Function), sizeof (char), 1, ptr_save);
#endif
#if USER_REANNEAL_PARAMETERS
      fread (&(OPTIONS->Reanneal_Params_Function), sizeof (char),
             1, ptr_save);
#endif
#if ASA_SAMPLE
      fread (&(OPTIONS->Bias_Acceptance), sizeof (double), 1, ptr_save);
      fread (OPTIONS->Bias_Generated, sizeof (double),
             *number_parameters, ptr_save);
      fread (&(OPTIONS->Average_Weights), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Limit_Weights), sizeof (double), 1, ptr_save);
#endif
#if ASA_QUEUE
      fread (&save_queue, sizeof (LONG_INT), 1, ptr_save);
      fread (&save_queue_indx, sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->Queue_Size), sizeof (ALLOC_INT), 1, ptr_save);
      fread (save_queue_flag, sizeof (int), save_queue, ptr_save);
      fread (save_queue_cost, sizeof (double), save_queue, ptr_save);
      fread (save_queue_param, sizeof (double),
             (*number_parameters) * (OPTIONS->Queue_Size), ptr_save);
#if ASA_RESOLUTION
#else
      fread (OPTIONS->Queue_Resolution, sizeof (double),
             *number_parameters, ptr_save);
#endif
#endif /* ASA_QUEUE */
#if ASA_RESOLUTION
      fread (OPTIONS->Coarse_Resolution, sizeof (double),
             *number_parameters, ptr_save);
#endif
#if FITLOC
      fread (&(OPTIONS->Fit_Local), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Iter_Max), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Penalty), sizeof (double), 1, ptr_save);
#endif
#if ASA_FUZZY
      fread (&(OPTIONS->NoOfSamples), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->ThresholdDeviation), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Performance_Target), sizeof (double), 1, ptr_save);
      fread (&(OPTIONS->Factor_a), sizeof (double), 1, ptr_save);
#endif
#if MULTI_MIN
      fread (OPTIONS->Multi_Number, sizeof (int), 1, ptr_save);
      fread (OPTIONS->Multi_Grid,
             sizeof (double), *number_parameters, ptr_save);
      fread (&(OPTIONS->Multi_Specify), sizeof (int), 1, ptr_save);
      for (multi_index = 0; multi_index < OPTIONS->Multi_Number;
           ++multi_index) {
        fread (&(OPTIONS->Multi_Cost[multi_index]), sizeof (double), 1,
               ptr_save);
        fread (&(OPTIONS->Multi_Params[multi_index]), sizeof (double),
               *number_parameters, ptr_save);
      }
#endif
#if ASA_PARALLEL
      for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
        fread (gener_block_state[i_prll].parameter,
               sizeof (double), *number_parameters, ptr_save);
        fread (&(gener_block_state[i_prll].cost),
               sizeof (double), 1, ptr_save);
#if USER_ACCEPTANCE_TEST
        fread (&
               (gener_block_state[i_prll].par_user_accept_flag),
               sizeof (int), 1, ptr_save);
        fread (&
               (gener_block_state[i_prll].par_cost_accept_flag),
               sizeof (int), 1, ptr_save);
#endif
      }
      fread (&(OPTIONS->Gener_Mov_Avr), sizeof (int), 1, ptr_save);
      fread (&(OPTIONS->Gener_Block), sizeof (LONG_INT), 1, ptr_save);
      fread (&(OPTIONS->Gener_Block_Max), sizeof (LONG_INT), 1, ptr_save);
#endif /* ASA_PARALLEL */

      fclose (ptr_save);

      asa_read = FALSE;
#if ASA_PRINT
      print_state (parameter_minimum,
                   parameter_maximum,
                   tangents,
                   curvature,
                   current_cost_temperature,
                   current_user_parameter_temp,
                   accepted_to_generated_ratio,
                   number_parameters,
                   curvature_flag,
                   number_accepted,
                   index_cost_acceptances,
                   number_generated,
                   number_invalid_generated_states,
                   last_saved_state,
                   best_generated_state, ptr_asa_out, OPTIONS);
#endif /* ASA_PRINT */

#include "asa_opt"
#if ASA_SAVE_OPT
      if ((ptr_save_opt = fopen ("asa_save_opt", "r")) == NULL) {
#if INCL_STDOUT
        printf ("\n\n*** WARNING fopen asa_save_opt failed *** \n\n");
#endif /* INCL_STDOUT */
#if ASA_PRINT
        fprintf (ptr_asa_out,
                 "\n\n*** WARNING fopen asa_save_opt failed *** \n\n");
        fflush (ptr_asa_out);
#endif
      } else {
        fscanf_ret = fscanf (ptr_save_opt, "%s%s%s%s%s",
                             read_if, read_FALSE, read_comm1, read_ASA_SAVE,
                             read_comm2);
        if (strcmp (read_if, "#if") || strcmp (read_FALSE, "FALSE")
            || strcmp (read_comm1, "/*")
            || strcmp (read_ASA_SAVE, "ASA_SAVE")
            || strcmp (read_comm2, "*/")) {
#if INCL_STDOUT
          printf ("\n\n*** EXIT not asa_save_opt for this version *** \n\n");
#endif /* INCL_STDOUT */
#if ASA_PRINT
          fprintf (ptr_asa_out,
                   "\n\n*** not asa_save_opt for this version *** \n\n");
          fflush (ptr_asa_out);
#endif
          *exit_status = INVALID_USER_INPUT;
          goto EXIT_asa;
        }
#if INT_LONG
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%ld", &read_long);
        OPTIONS->Limit_Acceptances = read_long;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%ld", &read_long);
        OPTIONS->Limit_Generated = read_long;
#else
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Limit_Acceptances = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Limit_Generated = read_int;
#endif
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Limit_Invalid_Generated_States = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Accepted_To_Generated_Ratio = read_double;

        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Cost_Precision = read_double;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Maximum_Cost_Repeat = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Number_Cost_Samples = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Temperature_Ratio_Scale = read_double;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Cost_Parameter_Scale_Ratio = read_double;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Temperature_Anneal_Scale = read_double;

        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Include_Integer_Parameters = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->User_Initial_Parameters = read_int;
#if INT_ALLOC
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Sequential_Parameters = read_int;
#else
#if INT_LONG
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%ld", &read_long);
        OPTIONS->Sequential_Parameters = read_long;
#else
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Sequential_Parameters = read_int;
#endif
#endif
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Initial_Parameter_Temperature = read_double;

        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Acceptance_Frequency_Modulus = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Generated_Frequency_Modulus = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Reanneal_Cost = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Reanneal_Parameters = read_int;

        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%lf", &read_double);
        OPTIONS->Delta_X = read_double;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->User_Tangents = read_int;
        fscanf_ret = fscanf (ptr_save_opt, "%s", read_option);
        fscanf_ret = fscanf (ptr_save_opt, "%d", &read_int);
        OPTIONS->Curvature_0 = read_int;

        fclose (ptr_save_opt);
      }
#endif /* ASA_SAVE_OPT */

      goto SAVED_asa;
    }
#endif /* ASA_SAVE */

    if (OPTIONS->Locate_Cost < 0) {
      OPTIONS->Locate_Cost = 12;        /* generate new state from new best */
    } else {
      OPTIONS->Locate_Cost = 2; /* generate new state */
    }

#if USER_ACCEPTANCE_TEST
    OPTIONS->User_Acceptance_Flag = FALSE;
    OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif

#if ASA_EXIT_ANYTIME
    if ((ptr_exit_anytime = fopen ("asa_exit_anytime", "r")) == NULL) {
      goto EXIT_asa;
    } else {
      fclose (ptr_exit_anytime);
    }
#endif /* ASA_EXIT_ANYTIME */

    /* GENERATE NEW PARAMETERS */

    /* generate a new valid set of parameters */
#if ASA_PARALLEL
    /* While this section of code is set to run under OpenMP using the gcc
     * compiler, you may change/add lines of code in this entire ASA_PARALLEL
     * section to correspond to your choice of parallel algorithm and
     * compiler. The entire ASA_PARALLEL section makes assignments to indexed
     * variables to afford flexibility for other such algorithms. */

    /* Note that here the do loop around generated states that tests for
     * invalid states is taken over, not within, blocks of parallel
     * calculated cost functions for these generated states, below. */

    repeated_invalid_states = 0;
    do {
      valid_state_generated_flag_par_test = 0;
      for (i_prll = 0; i_prll < OPTIONS->Gener_Block; ++i_prll) {
        valid_state_generated_flag_par[i_prll] = TRUE;
        number_invalid_generated_states_par[i_prll] =
          *number_invalid_generated_states;
#if USER_ACCEPTANCE_TEST
        gener_block_state[i_prll].par_user_accept_flag =
          OPTIONS->User_Acceptance_Flag;
        gener_block_state[i_prll].par_cost_accept_flag =
          OPTIONS->Cost_Acceptance_Flag;
#endif
      }

      for (i_prll = 0; i_prll < OPTIONS->Gener_Block; ++i_prll) {
        generate_flg_par[i_prll] =
          generate_new_state_par (user_random_generator, seed,
                                  parameter_minimum, parameter_maximum,
                                  current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                                  initial_user_parameter_temp,
                                  temperature_scale_parameters,
#endif
                                  number_parameters,
                                  parameter_type,
                                  i_prll,
                                  gener_block_state, last_saved_state,
                                  OPTIONS);

#if ASA_QUEUE
        /* Binary trees do not seem necessary since we are assuming
           that the cost function calculation is the bottleneck.
           However, see the MISC.DIR/asa_contrib file for
           source code for doubly-linked and hashed lists. */
        save_queue = (LONG_INT) OPTIONS->Queue_Size;
        if (OPTIONS->Queue_Size == 0) {
          queue_new_par[i_prll] = 1;
        } else {
          queue_new_par[i_prll] = 1;
          for (queue = 0; queue < save_queue; ++queue) {
            save_queue_test_par[i_prll] = 0;
            VFOR (index_v) {
              if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
                ++(save_queue_test_par[i_prll]);
              } else {
                queue_v_par[i_prll] =
                  index_v + queue * (LONG_INT) (*number_parameters);
                tmp_var_db_par[i_prll] =
                  fabs (gener_block_state[i_prll].parameter
                        [index_v] -
                        save_queue_param_par[i_prll][queue_v_par[i_prll]]);
                if (
#if ASA_RESOLUTION
                     /* Coarse_Resolution used in gener_block_state */
                     tmp_var_db_par[i_prll]
                     < EPS_DOUBLE
#else
                     tmp_var_db_par[i_prll]
                     < (OPTIONS->Queue_Resolution[index_v] + EPS_DOUBLE)
#endif /* ASA_RESOLUTION */
                  ) {
                  ++(save_queue_test_par[i_prll]);
                }
              }
            }
            if (save_queue_test_par[i_prll] == *number_parameters) {
              queue_par_cost[i_prll] = save_queue_cost_par[i_prll][queue];
              queue_new_par[i_prll] = 0;
              valid_state_generated_flag_par[i_prll] =
                save_queue_valid_state_flag_par[i_prll][queue];
              if (valid_state_generated_flag_par[i_prll] == FALSE) {
#if ASA_PRINT_MORE
#if INT_LONG
                fprintf (ptr_asa_out,
                         "ASA_QUEUE: %ld BlockID: %ld \t previous invalid state",
                         OPTIONS->N_Generated, i_prll);
#else
                fprintf (ptr_asa_out,
                         "ASA_QUEUE: %d BlockID: %d \t previous invalid state",
                         OPTIONS->N_Generated, i_prll);
#endif
#endif /* ASA_PRINT_MORE */
              } else if (valid_state_generated_flag_par[i_prll] == TRUE) {
#if ASA_PRINT_MORE
#if INT_LONG
                fprintf (ptr_asa_out, "ASA_QUEUE: %ld BlockID %ld \t %*.*g\n",
                         OPTIONS->N_Generated, i_prll,
                         G_FIELD, G_PRECISION, tmp_var_db_par[i_prll]);
#else
                fprintf (ptr_asa_out, "ASA_QUEUE: BlockID %d %d \t %*.*g\n",
                         OPTIONS->N_Generated, i_prll,
                         G_FIELD, G_PRECISION, tmp_var_db_par[i_prll]);
#endif
#endif /* ASA_PRINT_MORE */
              }
              break;
            }
          }
        }
#endif /* ASA_QUEUE */
      }

      /* *** ENTER CODE TO SPAWN OFF PARALLEL GENERATED STATES *** */
#ifdef _OPENMP
#pragma omp parallel for
#endif /* _OPENMP */
      for (i_prll = 0; i_prll < OPTIONS->Gener_Block; ++i_prll) {
#if ASA_QUEUE
        if (queue_new_par[i_prll] == 0) {
          gener_block_state[i_prll].cost = queue_par_cost[i_prll];
        } else {
#endif /* ASA_QUEUE */
          OPTIONS->parallel_id = i_prll;
          gener_block_state[i_prll].cost =
            user_cost_function (gener_block_state[i_prll].parameter,
                                parameter_minimum,
                                parameter_maximum,
                                tangents,
                                curvature,
                                number_parameters,
                                parameter_type,
                                &(valid_state_generated_flag_par[i_prll]),
                                exit_status, OPTIONS);
          tmp_var_db1_par[i_prll] =
            cost_function_test (gener_block_state[i_prll].cost,
                                gener_block_state[i_prll].parameter,
                                parameter_minimum, parameter_maximum,
                                number_parameters, xnumber_parameters);
          if (tmp_var_db1_par[i_prll] == 0) {
            EXIT_asa_parallel = 1;
          }
#if ASA_QUEUE
        }
#endif /* ASA_QUEUE */
      }
      /* *** EXIT CODE SPAWNING OFF PARALLEL GENERATED STATES *** */

      if (EXIT_asa_parallel == 1) {
        *exit_status = INVALID_COST_FUNCTION;
        goto EXIT_asa;
      }
#if ASA_QUEUE
      for (i_prll = 0; i_prll < OPTIONS->Gener_Block; ++i_prll) {
        if (valid_state_generated_flag_par[i_prll] == FALSE) {
          ++valid_state_generated_flag_par_test;
        }
        if (queue_new_par[i_prll] == 1) {
          if (OPTIONS->Queue_Size > 0) {        /* in case recursive use */
            VFOR (index_v) {
              if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
                continue;
              }
              queue_v_par[i_prll] = index_v + save_queue_indx_par[i_prll]
                * (LONG_INT) (*number_parameters);
              save_queue_param_par[i_prll][queue_v_par[i_prll]] =
                gener_block_state[i_prll].parameter[index_v];
            }
            save_queue_cost_par[i_prll][save_queue_indx_par[i_prll]] =
              gener_block_state[i_prll].cost;
            save_queue_valid_state_flag_par[i_prll][save_queue_indx_par
                                                    [i_prll]]
              = valid_state_generated_flag_par[i_prll];

            ++(save_queue_par[i_prll]);
            if (save_queue_par[i_prll] == (LONG_INT) OPTIONS->Queue_Size)
              --(save_queue_par[i_prll]);

            ++(save_queue_indx_par[i_prll]);
            if (save_queue_indx_par[i_prll] == (LONG_INT) OPTIONS->Queue_Size)
              save_queue_indx_par[i_prll] = 0;
          }
        }
      }
#endif /* ASA_QUEUE */
      repeated_invalid_states += valid_state_generated_flag_par_test;
    }
    while (valid_state_generated_flag_par_test >= OPTIONS->Gener_Block);

    if (repeated_invalid_states > OPTIONS->Limit_Invalid_Generated_States) {
      *exit_status = TOO_MANY_INVALID_STATES;
      goto EXIT_asa;
    }
#else /* ASA_PARALLEL */
    repeated_invalid_states = 0;
    do {
      ++(*number_invalid_generated_states);
      generate_flg = generate_new_state (user_random_generator,
                                         seed,
                                         parameter_minimum,
                                         parameter_maximum,
                                         current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                                         initial_user_parameter_temp,
                                         temperature_scale_parameters,
#endif
                                         number_parameters,
                                         parameter_type,
                                         current_generated_state,
                                         last_saved_state, OPTIONS);

      *valid_state_generated_flag = TRUE;
#if ASA_QUEUE
      /* Binary trees do not seem necessary since we are assuming
         that the cost function calculation is the bottleneck.
         However, see the MISC.DIR/asa_contrib file for
         source code for doubly-linked and hashed lists. */
      save_queue = (LONG_INT) OPTIONS->Queue_Size;
      if (OPTIONS->Queue_Size == 0) {
        queue_new = 1;
      } else {
        queue_new = 1;
        for (queue = 0; queue < save_queue; ++queue) {
          save_queue_test = 0;
          VFOR (index_v) {
            if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
              ++save_queue_test;
            } else {
              queue_v = index_v + queue * (LONG_INT) (*number_parameters);
              if (
#if ASA_RESOLUTION
                   /* Coarse_Resolution used in current_generated_state */
                   fabs (current_generated_state->parameter[index_v] -
                         save_queue_param[queue_v]) < EPS_DOUBLE
#else
                   fabs (current_generated_state->parameter[index_v] -
                         save_queue_param[queue_v]) <
                   (OPTIONS->Queue_Resolution[index_v] + EPS_DOUBLE)
#endif /* ASA_RESOLUTION */
                ) {
                ++save_queue_test;
              }
            }
          }
          if (save_queue_test == *number_parameters) {
            tmp_var_db = save_queue_cost[queue];
            queue_new = 0;
            *valid_state_generated_flag = save_queue_flag[queue];
            if (*valid_state_generated_flag == FALSE) {
#if ASA_PRINT_MORE
#if INT_LONG
              fprintf (ptr_asa_out,
                       "ASA_QUEUE: %ld \t previous invalid state",
                       OPTIONS->N_Generated);
#else
              fprintf (ptr_asa_out,
                       "ASA_QUEUE: %d \t previous invalid state",
                       OPTIONS->N_Generated);
#endif
#endif /* ASA_PRINT_MORE */
            } else {
#if ASA_PRINT_MORE
#if INT_LONG
              fprintf (ptr_asa_out, "ASA_QUEUE: %ld \t %*.*g\n",
                       OPTIONS->N_Generated,
                       G_FIELD, G_PRECISION, tmp_var_db);
#else
              fprintf (ptr_asa_out, "ASA_QUEUE: %d \t %*.*g\n",
                       OPTIONS->N_Generated,
                       G_FIELD, G_PRECISION, tmp_var_db);
#endif
#endif /* ASA_PRINT_MORE */
            }
            break;
          }
        }
      }
      if (queue_new == 1) {
        tmp_var_db =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if (cost_function_test (tmp_var_db,
                                current_generated_state->parameter,
                                parameter_minimum,
                                parameter_maximum,
                                number_parameters, xnumber_parameters) == 0) {
          *exit_status = INVALID_COST_FUNCTION;
          goto EXIT_asa;
        }
        if (OPTIONS->Queue_Size > 0) {  /* in case recursive use */
          VFOR (index_v) {
            if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
              continue;
            }
            queue_v = index_v + save_queue_indx
              * (LONG_INT) (*number_parameters);
            save_queue_param[queue_v] =
              current_generated_state->parameter[index_v];
          }
          save_queue_cost[save_queue_indx] = tmp_var_db;
          save_queue_flag[save_queue_indx]
            = *valid_state_generated_flag;

          ++save_queue;
          if (save_queue == (LONG_INT) OPTIONS->Queue_Size)
            --save_queue;

          ++save_queue_indx;
          if (save_queue_indx == (LONG_INT) OPTIONS->Queue_Size)
            save_queue_indx = 0;
        }
      }
#else /* ASA_QUEUE */
      tmp_var_db =
        user_cost_function (current_generated_state->parameter,
                            parameter_minimum,
                            parameter_maximum,
                            tangents,
                            curvature,
                            number_parameters,
                            parameter_type,
                            valid_state_generated_flag, exit_status, OPTIONS);
      if (cost_function_test
          (tmp_var_db, current_generated_state->parameter,
           parameter_minimum, parameter_maximum, number_parameters,
           xnumber_parameters) == 0) {
        *exit_status = INVALID_COST_FUNCTION;
        goto EXIT_asa;
      }
#endif /* ASA_QUEUE */
      current_generated_state->cost = tmp_var_db;
      ++repeated_invalid_states;
      if (repeated_invalid_states > OPTIONS->Limit_Invalid_Generated_States) {
        *exit_status = TOO_MANY_INVALID_STATES;
        goto EXIT_asa;
      }
    }
    while (*valid_state_generated_flag == FALSE);
    --(*number_invalid_generated_states);
#endif /* ASA_PARALLEL */

    /* ACCEPT/REJECT NEW PARAMETERS */

#if ASA_PARALLEL
    for (sort_index = 0; sort_index < OPTIONS->Gener_Block; ++sort_index) {
      parallel_sort[sort_index] = sort_index;
    }
    qsort (parallel_sort, OPTIONS->Gener_Block, sizeof (LONG_INT),
           sort_parallel);

    for (sort_index = 0; sort_index < OPTIONS->Gener_Block; ++sort_index) {
      i_prll = parallel_sort[sort_index];
      if (valid_state_generated_flag_par[i_prll] == FALSE) {
        continue;
      }
      current_generated_state->cost = gener_block_state[i_prll].cost;
#if USER_ACCEPTANCE_TEST
      OPTIONS->User_Acceptance_Flag =
        gener_block_state[i_prll].par_user_accept_flag;
      OPTIONS->Cost_Acceptance_Flag =
        gener_block_state[i_prll].par_cost_accept_flag;
#endif
      VFOR (index_v) {
        /* ignore parameters with too small a range */
        if (PARAMETER_RANGE_TOO_SMALL (index_v))
          continue;
        current_generated_state->parameter[index_v] =
          gener_block_state[i_prll].parameter[index_v];
      }
#endif /* ASA_PARALLEL */
      /* decide whether to accept/reject the new state */
      accept_new_state (user_random_generator,
                        seed,
                        parameter_minimum,
                        parameter_maximum, current_cost_temperature,
#if ASA_SAMPLE
                        current_user_parameter_temp,
#endif
                        number_parameters,
                        recent_number_acceptances,
                        number_accepted,
                        index_cost_acceptances,
                        number_acceptances_saved,
                        recent_number_generated,
                        number_generated,
                        index_parameter_generations,
                        current_generated_state, last_saved_state,
#if ASA_SAMPLE
                        ptr_asa_out,
#endif
                        OPTIONS);

#if ASA_PARALLEL
#else
#if ASA_PIPE_FILE
#if INT_ALLOC
      fprintf (ptr_asa_pipe, "%d", *number_generated);
#else
#if INT_LONG
      fprintf (ptr_asa_pipe, "%ld", *number_generated);
#else
      fprintf (ptr_asa_pipe, "%d", *number_generated);
#endif
#endif
#if INT_ALLOC
      fprintf (ptr_asa_pipe, "\t%d", *number_accepted);
#else
#if INT_LONG
      fprintf (ptr_asa_pipe, "\t%ld", *number_accepted);
#else
      fprintf (ptr_asa_pipe, "\t%d", *number_accepted);
#endif
#endif
      fprintf (ptr_asa_pipe, "\t%g", best_generated_state->cost);
      VFOR (index_v)
        fprintf (ptr_asa_pipe, "\t%g",
                 best_generated_state->parameter[index_v]);
      fprintf (ptr_asa_pipe, "\t%g", current_generated_state->cost);
      VFOR (index_v)
        fprintf (ptr_asa_pipe, "\t%g",
                 current_generated_state->parameter[index_v]);
      fprintf (ptr_asa_pipe, "\t%g", *current_cost_temperature);
      VFOR (index_v)
        fprintf (ptr_asa_pipe, "\t%g", current_user_parameter_temp[index_v]);
      fprintf (ptr_asa_pipe, "\t%g", last_saved_state->cost);
      fprintf (ptr_asa_pipe, "\n");
      fflush (ptr_asa_pipe);
#endif /* ASA_PIPE_FILE */
#if INCL_STDOUT
#if ASA_PIPE
#if INT_ALLOC
      printf ("%d", *number_generated);
#else
#if INT_LONG
      printf ("%ld", *number_generated);
#else
      printf ("%d", *number_generated);
#endif
#endif
#if INT_ALLOC
      printf ("\t%d", *number_accepted);
#else
#if INT_LONG
      printf ("\t%ld", *number_accepted);
#else
      printf ("\t%d", *number_accepted);
#endif
#endif
      printf ("\t%g", best_generated_state->cost);
      VFOR (index_v)
        printf ("\t%g", best_generated_state->parameter[index_v]);
      printf ("\t%g", current_generated_state->cost);
      VFOR (index_v)
        printf ("\t%g", current_generated_state->parameter[index_v]);
      printf ("\t%g", *current_cost_temperature);
      VFOR (index_v)
        printf ("\t%g", current_user_parameter_temp[index_v]);
      printf ("\t%g", last_saved_state->cost);
      printf ("\n");
#endif /* ASA_PIPE */
#endif /* INCL_STDOUT */
#endif /* ASA_PARALLEL */

      /* calculate the ratio of acceptances to generated states */
      *accepted_to_generated_ratio =
        (double) (*recent_number_acceptances + 1) /
        (double) (*recent_number_generated + 1);

#if MULTI_MIN
      if (((OPTIONS->Multi_Specify == 0)
           && (current_generated_state->cost <= best_generated_state->cost))
          || ((OPTIONS->Multi_Specify == 1)
              && (current_generated_state->cost <
                  best_generated_state->cost))) {
#if ASA_RESOLUTION
        VFOR (index_v) {
          if (OPTIONS->Multi_Grid[index_v] <
              OPTIONS->Coarse_Resolution[index_v])
            OPTIONS->Multi_Grid[index_v] =
              OPTIONS->Coarse_Resolution[index_v];
        }
#endif /* ASA_RESOLUTION */
        VFOR (index_v) {
          if (OPTIONS->Multi_Grid[index_v] < EPS_DOUBLE)
            OPTIONS->Multi_Grid[index_v] = EPS_DOUBLE;
        }

        multi_test = 0;
        for (multi_index = 0; multi_index < OPTIONS->Multi_Number;
             ++multi_index) {
          multi_test_cmp = 0;
          multi_test_dim = 0;
          VFOR (index_v) {
            if (PARAMETER_RANGE_TOO_SMALL (index_v))
              continue;
            ++multi_test_dim;
            if (fabs (current_generated_state->parameter[index_v]
                      - OPTIONS->Multi_Params[multi_index][index_v])
                < OPTIONS->Multi_Grid[index_v] - EPS_DOUBLE)
              ++multi_test_cmp;
          }
          if (multi_test_cmp == multi_test_dim)
            multi_test = 1;
          if (OPTIONS->Multi_Specify == 1)
            break;
        }

        if (multi_test == 0) {
          multi_cost[OPTIONS->Multi_Number] = current_generated_state->cost;
          VFOR (index_v) {
            multi_params[OPTIONS->Multi_Number][index_v] =
              current_generated_state->parameter[index_v];
          }
          for (multi_index = 0; multi_index < OPTIONS->Multi_Number;
               ++multi_index) {
            multi_cost[multi_index] = OPTIONS->Multi_Cost[multi_index];
            VFOR (index_v) {
              multi_params[multi_index][index_v] =
                OPTIONS->Multi_Params[multi_index][index_v];
            }
          }

          qsort (multi_sort, OPTIONS->Multi_Number + 1, sizeof (int),
                 multi_compare);
          for (multi_index = 0; multi_index < OPTIONS->Multi_Number;
               ++multi_index) {
            OPTIONS->Multi_Cost[multi_index] =
              multi_cost[multi_sort[multi_index]];
            VFOR (index_v) {
              OPTIONS->Multi_Params[multi_index][index_v] =
                multi_params[multi_sort[multi_index]][index_v];
            }
          }
        }
      }
#endif /* MULTI_MIN */

      /* CHECK FOR NEW MINIMUM */
      if (current_generated_state->cost < best_generated_state->cost) {
        best_flag = 1;
      } else {
        best_flag = 0;
      }
#if MULTI_MIN
      if (((OPTIONS->Multi_Specify == 0)
           && (current_generated_state->cost <= best_generated_state->cost))
          || ((OPTIONS->Multi_Specify == 1)
              && (current_generated_state->cost <
                  best_generated_state->cost)))
#else
      if (current_generated_state->cost < best_generated_state->cost)
#endif /* MULTI_MIN */
      {
        /* NEW MINIMUM FOUND */

        OPTIONS->Locate_Cost = -1;

        /* reset the recent acceptances and generated counts */
        *recent_number_acceptances = *recent_number_generated = 0;
        if (best_flag == 1) {
          *best_number_generated_saved = *number_generated;
          *best_number_accepted_saved = *number_accepted;
        }
        index_cost_repeat = 0;

        /* copy the current state into the best_generated state */
        best_generated_state->cost = current_generated_state->cost;
        VFOR (index_v) {
#if DROPPED_PARAMETERS
          /* ignore parameters that have too small a range */
          if (PARAMETER_RANGE_TOO_SMALL (index_v))
            continue;
#endif
          best_generated_state->parameter[index_v] =
            current_generated_state->parameter[index_v];
        }

        /* printout the new minimum state and value */
#if ASA_PRINT
#if INT_LONG
  char *bestc_format = "best...->cost=%-*.*g  *number_accepted=%ld  *number_generated=%ld\n";
#else
  char *bestc_format = "best...->cost=%-*.*g  *number_accepted=%d  *number_generated=%d\n";
#endif /* INT_LONG */
        fprintf (ptr_asa_out,
                 bestc_format,
                 G_FIELD, G_PRECISION, best_generated_state->cost,
                 *number_accepted, *number_generated);
#if ASA_PRINT_MORE
        if (best_flag == 1) {
          fprintf (ptr_asa_out, "\nnew best\n");
        }
#endif /* ASA_PRINT_MORE */
#if ASA_PARALLEL
#if INT_LONG
  char *asapar_format = "OPTIONS->Gener_Block = %ld\n";
#else
  char *asapar_format = "OPTIONS->Gener_Block = %d\n";
#endif /* INT_LONG */
        /* print OPTIONS->Gener_Block just used */
        fprintf (ptr_asa_out,
                 asapar_format,
                 OPTIONS->Gener_Block);
#endif /* ASA_PARALLEL */
        if (best_flag == 1) {
#if ASA_PRINT_MORE
#if INT_ALLOC
          fprintf (ptr_asa_out, "Present Random Seed = %d\n", *seed);
#else
#if INT_LONG
          fprintf (ptr_asa_out, "Present Random Seed = %ld\n", *seed);
#else
          fprintf (ptr_asa_out, "Present Random Seed = %d\n", *seed);
#endif
#endif
          print_state (parameter_minimum,
                       parameter_maximum,
                       tangents,
                       curvature,
                       current_cost_temperature,
                       current_user_parameter_temp,
                       accepted_to_generated_ratio,
                       number_parameters,
                       curvature_flag,
                       number_accepted,
                       index_cost_acceptances,
                       number_generated,
                       number_invalid_generated_states,
                       last_saved_state,
                       best_generated_state, ptr_asa_out, OPTIONS);
#endif /* ASA_PRINT_MORE */
        }
        fflush (ptr_asa_out);
#endif /* ASA_PRINT */
      }
#if ASA_PARALLEL
    }
#endif /* ASA_PARALLEL */

#if ASA_SAVE
    /* These writes are put here with these tests, instead of just
       after a new best state is found, to prevent any confusion with
       any parallel code that might be added by users. */
    if (*recent_number_acceptances == 0
        && *recent_number_generated == 0
        && *best_number_generated_saved == *number_generated
        && *best_number_accepted_saved == *number_accepted
        && OPTIONS->Asa_Recursive_Level == asa_recursive_max
        && index_cost_repeat == 0) {
      if (OPTIONS->Asa_Recursive_Level > 0)
        sprintf (asa_save_comm, "asa_save_%d", OPTIONS->Asa_Recursive_Level);
      else
        sprintf (asa_save_comm, "asa_save");
      ptr_save = fopen (asa_save_comm, "w");

      fwrite (number_parameters, sizeof (ALLOC_INT), 1, ptr_save);
      fwrite (xnumber_parameters, sizeof (double), 1, ptr_save);
      fwrite (parameter_minimum, sizeof (double), *number_parameters,
              ptr_save);
      fwrite (parameter_maximum, sizeof (double), *number_parameters,
              ptr_save);
      fwrite (tangents, sizeof (double), *number_parameters, ptr_save);
      fwrite (current_user_parameter_temp, sizeof (double),
              *number_parameters, ptr_save);
      fwrite (initial_user_parameter_temp, sizeof (double),
              *number_parameters, ptr_save);
      fwrite (temperature_scale_parameters, sizeof (double),
              *number_parameters, ptr_save);

      fwrite (parameter_type, sizeof (int), *number_parameters, ptr_save);
      fwrite (&index_cost_repeat, sizeof (int), 1, ptr_save);
      fwrite (&asa_open, sizeof (int), 1, ptr_save);
      fwrite (&number_asa_open, sizeof (int), 1, ptr_save);
      fwrite (&recursive_asa_open, sizeof (int), 1, ptr_save);

      fwrite (current_cost_temperature, sizeof (double), 1, ptr_save);
      fwrite (initial_cost_temperature, sizeof (double), 1, ptr_save);
      fwrite (temperature_scale_cost, sizeof (double), 1, ptr_save);
      fwrite (accepted_to_generated_ratio, sizeof (double), 1, ptr_save);

      fwrite (curvature_flag, sizeof (int), 1, ptr_save);

      fwrite (seed, sizeof (LONG_INT), 1, ptr_save);
      fwrite (number_generated, sizeof (LONG_INT), 1, ptr_save);
      fwrite (number_accepted, sizeof (LONG_INT), 1, ptr_save);
      fwrite (number_acceptances_saved, sizeof (LONG_INT), 1, ptr_save);
      fwrite (recent_number_acceptances, sizeof (LONG_INT), 1, ptr_save);
      fwrite (recent_number_generated, sizeof (LONG_INT), 1, ptr_save);
      fwrite (number_invalid_generated_states, sizeof (LONG_INT), 1,
              ptr_save);
      fwrite (index_cost_acceptances, sizeof (LONG_INT), 1, ptr_save);
      fwrite (best_number_generated_saved, sizeof (LONG_INT), 1, ptr_save);
      fwrite (best_number_accepted_saved, sizeof (LONG_INT), 1, ptr_save);

      fwrite (index_parameter_generations, sizeof (LONG_INT),
              *number_parameters, ptr_save);

      fwrite (current_generated_state->parameter,
              sizeof (double), *number_parameters, ptr_save);
      fwrite (last_saved_state->parameter,
              sizeof (double), *number_parameters, ptr_save);
      fwrite (best_generated_state->parameter,
              sizeof (double), *number_parameters, ptr_save);
      fwrite (&(current_generated_state->cost), sizeof (double), 1, ptr_save);
      fwrite (&(last_saved_state->cost), sizeof (double), 1, ptr_save);
      fwrite (&(best_generated_state->cost), sizeof (double), 1, ptr_save);

      fwrite (&(OPTIONS->Limit_Acceptances), sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->Limit_Generated), sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->Limit_Invalid_Generated_States), sizeof (int),
              1, ptr_save);
      fwrite (&(OPTIONS->Accepted_To_Generated_Ratio), sizeof (double),
              1, ptr_save);
      fwrite (&(OPTIONS->Cost_Precision), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Maximum_Cost_Repeat), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Number_Cost_Samples), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Temperature_Ratio_Scale), sizeof (double),
              1, ptr_save);
      fwrite (&(OPTIONS->Cost_Parameter_Scale_Ratio), sizeof (double),
              1, ptr_save);
      fwrite (&(OPTIONS->Temperature_Anneal_Scale), sizeof (double),
              1, ptr_save);
      fwrite (&(OPTIONS->Include_Integer_Parameters), sizeof (int),
              1, ptr_save);
      fwrite (&(OPTIONS->User_Initial_Parameters), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Sequential_Parameters), sizeof (ALLOC_INT), 1,
              ptr_save);
      fwrite (&(OPTIONS->Initial_Parameter_Temperature), sizeof (double), 1,
              ptr_save);
      fwrite (&(OPTIONS->Acceptance_Frequency_Modulus), sizeof (int), 1,
              ptr_save);
      fwrite (&(OPTIONS->Generated_Frequency_Modulus), sizeof (int), 1,
              ptr_save);
      fwrite (&(OPTIONS->Reanneal_Cost), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Reanneal_Parameters), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Delta_X), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->User_Tangents), sizeof (int), 1, ptr_save);

#if USER_INITIAL_COST_TEMP
      fwrite (&(OPTIONS->User_Cost_Temperature), sizeof (double), 1,
              ptr_save);
#endif
#if RATIO_TEMPERATURE_SCALES
      fwrite (OPTIONS->User_Temperature_Ratio, sizeof (double),
              *number_parameters, ptr_save);
#endif
#if USER_INITIAL_PARAMETERS_TEMPS
      fwrite (OPTIONS->User_Parameter_Temperature, sizeof (double),
              *number_parameters, ptr_save);
#endif
#if DELTA_PARAMETERS
      fwrite (OPTIONS->User_Delta_Parameter, sizeof (double),
              *number_parameters, ptr_save);
#endif
#if QUENCH_PARAMETERS
      fwrite (OPTIONS->User_Quench_Param_Scale, sizeof (double),
              *number_parameters, ptr_save);
#endif
#if QUENCH_COST
      fwrite (OPTIONS->User_Quench_Cost_Scale, sizeof (double), 1, ptr_save);
#endif
      fwrite (&(OPTIONS->N_Accepted), sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->N_Generated), sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->Locate_Cost), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Immediate_Exit), sizeof (int), 1, ptr_save);
#if OPTIONAL_DATA_DBL
      fwrite (&(OPTIONS->Asa_Data_Dim_Dbl), sizeof (ALLOC_INT), 1, ptr_save);
      fwrite (OPTIONS->Asa_Data_Dbl, sizeof (double),
              OPTIONS->Asa_Data_Dim_Dbl, ptr_save);
#endif
      fwrite (&(OPTIONS->Random_Array_Dim), sizeof (ALLOC_INT), 1, ptr_save);
      fwrite (OPTIONS->Random_Array, sizeof (double),
              OPTIONS->Random_Array_Dim, ptr_save);
      fwrite (&(OPTIONS->Asa_Recursive_Level), sizeof (int), 1, ptr_save);
#if OPTIONAL_DATA_INT
      fwrite (&(OPTIONS->Asa_Data_Dim_Int), sizeof (ALLOC_INT), 1, ptr_save);
      fwrite (OPTIONS->Asa_Data_Int, sizeof (LONG_INT),
              OPTIONS->Asa_Data_Dim_Int, ptr_save);
#endif
#if OPTIONAL_DATA_PTR
      fwrite (&(OPTIONS->Asa_Data_Dim_Ptr), sizeof (ALLOC_INT), 1, ptr_save);
      if (OPTIONS->Asa_Recursive_Level == 0)
        fwrite (OPTIONS->Asa_Data_Ptr, sizeof (OPTIONAL_PTR_TYPE),
                OPTIONS->Asa_Data_Dim_Ptr, ptr_save);
#if ASA_TEMPLATE_SELFOPT
      if (OPTIONS->Asa_Recursive_Level == 1)
        fwrite (OPTIONS->Asa_Data_Ptr, sizeof (RECUR_OPTIONAL_PTR_TYPE),
                OPTIONS->Asa_Data_Dim_Ptr, ptr_save);
#endif
#endif
#if USER_ASA_OUT
      fwrite (OPTIONS->Asa_Out_File, sizeof (char), 80, ptr_save);
#endif
#if USER_ASA_OUT
      fwrite (OPTIONS->Asa_Usr_Out_File, sizeof (char), 80, ptr_save);
#endif
#if USER_COST_SCHEDULE
      fwrite (&(OPTIONS->Cost_Schedule), sizeof (char), 1, ptr_save);
#endif
#if USER_ACCEPT_ASYMP_EXP
      fwrite (&(OPTIONS->Asymp_Exp_Param), sizeof (double), 1, ptr_save);
#endif
#if USER_ACCEPTANCE_TEST
      fwrite (&(OPTIONS->Acceptance_Test), sizeof (char), 1, ptr_save);
      fwrite (&(OPTIONS->User_Acceptance_Flag), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Cost_Acceptance_Flag), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Cost_Temp_Curr), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Cost_Temp_Init), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Cost_Temp_Scale), sizeof (double), 1, ptr_save);
#endif
#if USER_GENERATING_FUNCTION
      fwrite (&(OPTIONS->Generating_Distrib), sizeof (char), 1, ptr_save);
#endif
#if USER_REANNEAL_COST
      fwrite (&(OPTIONS->Reanneal_Cost_Function), sizeof (char), 1, ptr_save);
#endif
#if USER_REANNEAL_PARAMETERS
      fwrite (&(OPTIONS->Reanneal_Params_Function), sizeof (char), 1,
              ptr_save);
#endif
#if ASA_SAMPLE
      fwrite (&(OPTIONS->Bias_Acceptance), sizeof (double), 1, ptr_save);
      fwrite (OPTIONS->Bias_Generated, sizeof (double),
              *number_parameters, ptr_save);
      fwrite (&(OPTIONS->Average_Weights), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Limit_Weights), sizeof (double), 1, ptr_save);
#endif
#if ASA_QUEUE
      fwrite (&save_queue, sizeof (LONG_INT), 1, ptr_save);
      fwrite (&save_queue_indx, sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->Queue_Size), sizeof (ALLOC_INT), 1, ptr_save);
      fwrite (save_queue_flag, sizeof (int), save_queue, ptr_save);
      fwrite (save_queue_cost, sizeof (double), save_queue, ptr_save);
      fwrite (save_queue_param, sizeof (double),
              (*number_parameters) * (OPTIONS->Queue_Size), ptr_save);
#if ASA_RESOLUTION
#else
      fwrite (OPTIONS->Queue_Resolution, sizeof (double),
              *number_parameters, ptr_save);
#endif
#endif /* ASA_QUEUE */
#if ASA_RESOLUTION
      fwrite (OPTIONS->Coarse_Resolution, sizeof (double),
              *number_parameters, ptr_save);
#endif
#if ASA_FUZZY
      fwrite (&(OPTIONS->NoOfSamples), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->ThresholdDeviation), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Performance_Target), sizeof (double), 1, ptr_save);
      fwrite (&(OPTIONS->Factor_a), sizeof (double), 1, ptr_save);
#endif
#if FITLOC
      fwrite (&(OPTIONS->Fit_Local), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Iter_Max), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Penalty), sizeof (double), 1, ptr_save);
#endif
#if MULTI_MIN
      fwrite (OPTIONS->Multi_Number, sizeof (int), 1, ptr_save);
      fwrite (OPTIONS->Multi_Grid,
              sizeof (double), *number_parameters, ptr_save);
      fwrite (&(OPTIONS->Multi_Specify), sizeof (int), 1, ptr_save);
      for (multi_index = 0; multi_index < OPTIONS->Multi_Number;
           ++multi_index) {
        fwrite (&(OPTIONS->Multi_Cost[multi_index]), sizeof (double), 1,
                ptr_save);
        fwrite (&(OPTIONS->Multi_Params[multi_index]), sizeof (double),
                *number_parameters, ptr_save);
      }
#endif
#if ASA_PARALLEL
      for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
        fwrite (gener_block_state[i_prll].parameter,
                sizeof (double), *number_parameters, ptr_save);
        fwrite (&(gener_block_state[i_prll].cost),
                sizeof (double), 1, ptr_save);
#if USER_ACCEPTANCE_TEST
        fwrite (&
                (gener_block_state[i_prll].par_user_accept_flag),
                sizeof (int), 1, ptr_save);
        fwrite (&(gener_block_state[i_prll].par_cost_accept_flag),
                sizeof (int), 1, ptr_save);
#endif
      }
      fwrite (&(OPTIONS->Gener_Mov_Avr), sizeof (int), 1, ptr_save);
      fwrite (&(OPTIONS->Gener_Block), sizeof (LONG_INT), 1, ptr_save);
      fwrite (&(OPTIONS->Gener_Block_Max), sizeof (LONG_INT), 1, ptr_save);
#endif /* ASA_PARALLEL */

      fclose (ptr_save);

    SAVED_asa:
      ;

#if SYSTEM_CALL
#if ASA_SAVE_BACKUP
#if INT_LONG
      if (OPTIONS->Asa_Recursive_Level > 0)
        sprintf (asa_save_comm, "/bin/cp asa_save_%d asa_save_%d.%ld",
                 OPTIONS->Asa_Recursive_Level,
                 OPTIONS->Asa_Recursive_Level, OPTIONS->N_Accepted);
      else
        sprintf (asa_save_comm, "/bin/cp asa_save asa_save.%ld",
                 OPTIONS->N_Accepted);
#else
      if (OPTIONS->Asa_Recursive_Level > 0)
        sprintf (asa_save_comm, "/bin/cp asa_save_%d asa_save_%d.%d",
                 OPTIONS->Asa_Recursive_Level,
                 OPTIONS->Asa_Recursive_Level, OPTIONS->N_Accepted);
      else
        sprintf (asa_save_comm, "/bin/cp asa_save asa_save.%d",
                 OPTIONS->N_Accepted);
#endif
      ptr_comm = popen (asa_save_comm, "r");
      pclose (ptr_comm);
#else /* ASA_SAVE_BACKUP */
      /* extra protection in case run aborts during write */
      if (OPTIONS->Asa_Recursive_Level > 0)
        sprintf (asa_save_comm, "/bin/cp asa_save_%d asa_save_%d.old",
                 OPTIONS->Asa_Recursive_Level, OPTIONS->Asa_Recursive_Level);
      else
        sprintf (asa_save_comm, "/bin/cp asa_save asa_save.old");
      ptr_comm = popen (asa_save_comm, "r");
      pclose (ptr_comm);
#endif /* ASA_SAVE_BACKUP */
#endif /* SYSTEM_CALL */
    }
#endif /* ASA_SAVE */

    if (OPTIONS->Immediate_Exit == TRUE) {
      *exit_status = IMMEDIATE_EXIT;
      goto EXIT_asa;
    }

    /* PERIODIC TESTING/REANNEALING/PRINTING SECTION */

    if (OPTIONS->Acceptance_Frequency_Modulus == 0)
      tmp_var_int1 = FALSE;
    else if ((int) (*number_accepted %
                    ((LONG_INT) OPTIONS->Acceptance_Frequency_Modulus)) == 0
             && *number_acceptances_saved == *number_accepted)
      tmp_var_int1 = TRUE;
    else
      tmp_var_int1 = FALSE;

    if (OPTIONS->Generated_Frequency_Modulus == 0)
      tmp_var_int2 = FALSE;
    else if ((int) (*number_generated %
                    ((LONG_INT) OPTIONS->Generated_Frequency_Modulus)) == 0)
      tmp_var_int2 = TRUE;
    else
      tmp_var_int2 = FALSE;

    if (tmp_var_int1 == TRUE || tmp_var_int2 == TRUE
        || (*accepted_to_generated_ratio
            < OPTIONS->Accepted_To_Generated_Ratio)) {

#if ASA_PARALLEL
      if (OPTIONS->Gener_Mov_Avr > 0) {
        for (i_prll = 1; i_prll < OPTIONS->Gener_Mov_Avr; ++i_prll) {
          parallel_gen_ratio_block[i_prll - 1] =
            parallel_gen_ratio_block[i_prll];
        }
        parallel_gen_ratio_block[OPTIONS->Gener_Mov_Avr - 1] =
          *recent_number_generated;
        tmp_var_lint = 0;
        for (i_prll = 0; i_prll < OPTIONS->Gener_Mov_Avr; ++i_prll) {
          tmp_var_lint += parallel_gen_ratio_block[i_prll];
        }
        OPTIONS->Gener_Block = (LONG_INT)
          ((double) tmp_var_lint / (double) (OPTIONS->Gener_Mov_Avr));
        OPTIONS->Gener_Block =
          MIN (OPTIONS->Gener_Block, OPTIONS->Gener_Block_Max);
      }
#endif /* ASA_PARALLEL */

      if (*accepted_to_generated_ratio <
          (OPTIONS->Accepted_To_Generated_Ratio))
        *recent_number_acceptances = *recent_number_generated = 0;

      /* if best.cost repeats OPTIONS->Maximum_Cost_Repeat then exit */
      if (OPTIONS->Maximum_Cost_Repeat != 0) {
        if (fabs (last_saved_state->cost - best_generated_state->cost)
            < OPTIONS->Cost_Precision) {
          ++index_cost_repeat;
          if (index_cost_repeat == (OPTIONS->Maximum_Cost_Repeat)) {
            *exit_status = COST_REPEATING;
            goto EXIT_asa;
          }
        } else {
          index_cost_repeat = 0;
        }
      }
      if (OPTIONS->Reanneal_Parameters == TRUE) {
        OPTIONS->Locate_Cost = 3;       /* reanneal parameters */

        /* calculate tangents, not curvatures, to reanneal */
        *curvature_flag = FALSE;
        cost_derivatives (user_cost_function,
                          parameter_minimum,
                          parameter_maximum,
                          tangents,
                          curvature,
                          maximum_tangent,
                          number_parameters,
                          parameter_type,
                          exit_status,
                          curvature_flag,
                          valid_state_generated_flag,
                          number_invalid_generated_states,
                          current_generated_state,
                          best_generated_state, ptr_asa_out, OPTIONS);
        if (*exit_status == INVALID_COST_FUNCTION_DERIV) {
          goto EXIT_asa;
        }
      }
#if USER_REANNEAL_COST
#else
      if (OPTIONS->Reanneal_Cost == 0 || OPTIONS->Reanneal_Cost == 1) {
        ;
      } else {
        immediate_flag = OPTIONS->Immediate_Exit;

        if (OPTIONS->Reanneal_Cost < -1) {
          tmp_var_int = -OPTIONS->Reanneal_Cost;
        } else {
          tmp_var_int = OPTIONS->Reanneal_Cost;
        }
        tmp_var_db1 = ZERO;
        tmp_var_db2 = ZERO;

        for (index_cost_constraint = 0;
             index_cost_constraint < tmp_var_int; ++index_cost_constraint) {
          OPTIONS->Locate_Cost = 4;     /* reanneal cost */

          *number_invalid_generated_states = 0;
          repeated_invalid_states = 0;
          OPTIONS->Sequential_Parameters = *start_sequence - 1;
          do {
#if ASA_EXIT_ANYTIME
            if ((ptr_exit_anytime = fopen ("asa_exit_anytime", "r")) == NULL) {
              *exit_status = IMMEDIATE_EXIT;
              goto EXIT_asa;
            } else {
              fclose (ptr_exit_anytime);
            }
#endif /* ASA_EXIT_ANYTIME */
            ++(*number_invalid_generated_states);
            generate_flg = generate_new_state (user_random_generator,
                                               seed,
                                               parameter_minimum,
                                               parameter_maximum,
                                               current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                                               initial_user_parameter_temp,
                                               temperature_scale_parameters,
#endif
                                               number_parameters,
                                               parameter_type,
                                               current_generated_state,
                                               last_saved_state, OPTIONS);
            *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
            OPTIONS->User_Acceptance_Flag = TRUE;
            OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif

#if ASA_QUEUE
            if (OPTIONS->Queue_Size == 0) {
              queue_new = 1;
            } else {
              queue_new = 1;
              for (queue = 0; queue < save_queue; ++queue) {
                save_queue_test = 0;
                VFOR (index_v) {
                  if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
                    ++save_queue_test;
                  } else {
                    queue_v =
                      index_v + queue * (LONG_INT) (*number_parameters);
                    if (
#if ASA_RESOLUTION
                         /* Coarse_Resolution used in current_generated_state */
                         fabs (current_generated_state->parameter[index_v] -
                               save_queue_param[queue_v]) < EPS_DOUBLE
#else
                         fabs (current_generated_state->parameter[index_v] -
                               save_queue_param[queue_v]) <
                         (OPTIONS->Queue_Resolution[index_v] + EPS_DOUBLE)
#endif /* ASA_RESOLUTION */
                      ) {
                      ++save_queue_test;
                    }
                  }
                }
                if (save_queue_test == *number_parameters) {
                  tmp_var_db = save_queue_cost[queue];
                  queue_new = 0;
                  *valid_state_generated_flag = save_queue_flag[queue];
                  if (*valid_state_generated_flag == FALSE) {
#if ASA_PRINT_MORE
#if INT_LONG
                    fprintf (ptr_asa_out,
                             "ASA_QUEUE: %ld \t previous invalid state",
                             OPTIONS->N_Generated);
#else
                    fprintf (ptr_asa_out,
                             "ASA_QUEUE: %d \t previous invalid state",
                             OPTIONS->N_Generated);
#endif
#endif /* ASA_PRINT_MORE */
                  } else {
#if ASA_PRINT_MORE
#if INT_LONG
                    fprintf (ptr_asa_out, "ASA_QUEUE: %ld \t %*.*g\n",
                             OPTIONS->N_Generated,
                             G_FIELD, G_PRECISION, tmp_var_db);
#else
                    fprintf (ptr_asa_out, "ASA_QUEUE: %d \t %*.*g\n",
                             OPTIONS->N_Generated,
                             G_FIELD, G_PRECISION, tmp_var_db);
#endif
#endif /* ASA_PRINT_MORE */
                  }
                  break;
                }
              }
            }
            if (queue_new == 1) {
              tmp_var_db =
                user_cost_function (current_generated_state->parameter,
                                    parameter_minimum, parameter_maximum,
                                    tangents, curvature, number_parameters,
                                    parameter_type,
                                    valid_state_generated_flag, exit_status,
                                    OPTIONS);
              if (cost_function_test
                  (tmp_var_db, current_generated_state->parameter,
                   parameter_minimum, parameter_maximum, number_parameters,
                   xnumber_parameters) == 0) {
                *exit_status = INVALID_COST_FUNCTION;
                goto EXIT_asa;
              }
              if (OPTIONS->Queue_Size > 0) {
                VFOR (index_v) {
                  if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
                    continue;
                  }
                  queue_v = index_v + save_queue
                    * (LONG_INT) (*number_parameters);
                  save_queue_param[queue_v] =
                    current_generated_state->parameter[index_v];
                }
                save_queue_cost[save_queue] = tmp_var_db;
                save_queue_flag[save_queue]
                  = *valid_state_generated_flag;

                ++save_queue;
                if (save_queue == (LONG_INT) OPTIONS->Queue_Size)
                  --save_queue;

                ++save_queue_indx;
                if (save_queue_indx == (LONG_INT) OPTIONS->Queue_Size)
                  save_queue_indx = 0;
              }
            }
#else /* ASA_QUEUE */
            tmp_var_db =
              user_cost_function (current_generated_state->parameter,
                                  parameter_minimum, parameter_maximum,
                                  tangents, curvature, number_parameters,
                                  parameter_type, valid_state_generated_flag,
                                  exit_status, OPTIONS);
            if (cost_function_test
                (tmp_var_db, current_generated_state->parameter,
                 parameter_minimum, parameter_maximum, number_parameters,
                 xnumber_parameters) == 0) {
              *exit_status = INVALID_COST_FUNCTION;
              goto EXIT_asa;
            }
#endif /* ASA_QUEUE */
            ++repeated_invalid_states;
            if (repeated_invalid_states >
                OPTIONS->Limit_Invalid_Generated_States) {
              *exit_status = TOO_MANY_INVALID_STATES;
              goto EXIT_asa;
            }
          }
          while (*valid_state_generated_flag == FALSE);
          --(*number_invalid_generated_states);

          tmp_var_db1 += tmp_var_db;
          tmp_var_db2 += (tmp_var_db * tmp_var_db);
        }
        tmp_var_db1 /= (double) tmp_var_int;
        tmp_var_db2 /= (double) tmp_var_int;
        tmp_var_db =
          sqrt (fabs
                ((tmp_var_db2 -
                  tmp_var_db1 * tmp_var_db1) * ((double) tmp_var_int /
                                                ((double) tmp_var_int -
                                                 ONE))));
        if (OPTIONS->Reanneal_Cost < -1) {
          *current_cost_temperature = *initial_cost_temperature =
            tmp_var_db + (double) EPS_DOUBLE;
        } else {
          *initial_cost_temperature = tmp_var_db + (double) EPS_DOUBLE;
        }
        OPTIONS->Immediate_Exit = immediate_flag;
      }
#endif /* USER_REANNEAL_COST */

      reanneal (parameter_minimum,
                parameter_maximum,
                tangents,
                maximum_tangent,
                current_cost_temperature,
                initial_cost_temperature,
                temperature_scale_cost,
                current_user_parameter_temp,
                initial_user_parameter_temp,
                temperature_scale_parameters,
                number_parameters,
                parameter_type,
                index_cost_acceptances,
                index_parameter_generations,
                last_saved_state, best_generated_state, OPTIONS);
#if ASA_PRINT_INTERMED
#if ASA_PRINT
      print_state (parameter_minimum,
                   parameter_maximum,
                   tangents,
                   curvature,
                   current_cost_temperature,
                   current_user_parameter_temp,
                   accepted_to_generated_ratio,
                   number_parameters,
                   curvature_flag,
                   number_accepted,
                   index_cost_acceptances,
                   number_generated,
                   number_invalid_generated_states,
                   last_saved_state,
                   best_generated_state, ptr_asa_out, OPTIONS);

      fprintf (ptr_asa_out, "\n");
      fflush (ptr_asa_out);
#endif
#endif
    }
  }

  /* FINISHED ANNEALING and MINIMIZATION */

  *exit_status = NORMAL_EXIT;

EXIT_asa:

  asa_exit_value = asa_exit (user_cost_function,
                             &final_cost,
                             parameter_initial_final,
                             parameter_minimum,
                             parameter_maximum,
                             tangents,
                             curvature,
                             maximum_tangent,
                             current_cost_temperature,
                             initial_user_parameter_temp,
                             current_user_parameter_temp,
                             accepted_to_generated_ratio,
                             number_parameters,
                             parameter_type,
                             valid_state_generated_flag,
                             exit_status,
                             index_exit_v,
                             start_sequence,
                             number_accepted,
                             best_number_accepted_saved,
                             index_cost_acceptances,
                             number_generated,
                             number_invalid_generated_states,
                             index_parameter_generations,
                             best_number_generated_saved,
                             current_generated_state,
                             last_saved_state,
                             best_generated_state, ptr_asa_out, OPTIONS);
  if (asa_exit_value == 9) {
    *exit_status = CALLOC_FAILED;
    return (-1);
  }

  free (curvature_flag);
  free (maximum_tangent);
  free (accepted_to_generated_ratio);
  free (temperature_scale_cost);
  free (current_cost_temperature);
  free (initial_cost_temperature);
  free (number_generated);
  free (best_number_generated_saved);
  free (recent_number_generated);
  free (number_accepted);
  free (recent_number_acceptances);
  free (index_cost_acceptances);
  free (number_acceptances_saved);
  free (best_number_accepted_saved);
  free (number_invalid_generated_states);
  free (current_generated_state->parameter);
  free (last_saved_state->parameter);
  free (best_generated_state->parameter);
  free (current_generated_state);
  free (last_saved_state);
  free (best_generated_state);
#if ASA_QUEUE
  free (save_queue_flag);
  free (save_queue_cost);
  free (save_queue_param);
#endif /* ASA_QUEUE */
#if MULTI_MIN
  for (multi_index = 0; multi_index <= OPTIONS->Multi_Number; ++multi_index)
    free (multi_params[multi_index]);
  free (multi_params);
  free (multi_sort);
  free (multi_cost);
#endif
#if ASA_PARALLEL
  for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
    free (gener_block_state[i_prll].parameter);
  }
  free (gener_block_state);
  free (parallel_sort);
  free (parallel_gen_ratio_block);
  free (generate_flg_par);
  free (number_invalid_generated_states_par);
  free (repeated_invalid_states_par);
  free (tmp_var_db1_par);
  free (tmp_var_db_par);
  free (valid_state_generated_flag_par);
#if ASA_QUEUE
  free (queue_par_cost);
  for (i_prll = 0; i_prll < OPTIONS->Gener_Block_Max; ++i_prll) {
    free (save_queue_param_par[i_prll]);
    free (save_queue_valid_state_flag_par[i_prll]);
    free (save_queue_cost_par[i_prll]);
  }
  free (save_queue_valid_state_flag_par);
  free (save_queue_cost_par);
  free (save_queue_param_par);
  free (queue_new_par);
  free (queue_v_par);
  free (save_queue_indx_par);
  free (save_queue_test_par);
  free (save_queue_par);
#endif /* ASA_QUEUE */
#endif /* ASA_PARALLEL */
#if ASA_PIPE_FILE
  fclose (ptr_asa_pipe);
#endif
  free (initial_user_parameter_temp);
  free (index_exit_v);
  free (start_sequence);
  free (index_parameter_generations);
  free (current_user_parameter_temp);
  free (temperature_scale_parameters);
  if (recursive_asa_open == 0)
    asa_open = FALSE;
  return (final_cost);
}

/***********************************************************************
* asa_exit
*	This procedures copies the best parameters and cost into
*       final_cost and parameter_initial_final
***********************************************************************/
#if HAVE_ANSI
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
          USER_DEFINES * OPTIONS)
#else
int

asa_exit (user_cost_function,
          final_cost,
          parameter_initial_final,
          parameter_minimum,
          parameter_maximum,
          tangents,
          curvature,
          maximum_tangent,
          current_cost_temperature,
          initial_user_parameter_temp,
          current_user_parameter_temp,
          accepted_to_generated_ratio,
          number_parameters,
          parameter_type,
          valid_state_generated_flag,
          exit_status,
          index_exit_v,
          start_sequence,
          number_accepted,
          best_number_accepted_saved,
          index_cost_acceptances,
          number_generated,
          number_invalid_generated_states,
          index_parameter_generations,
          best_number_generated_saved,
          current_generated_state,
          last_saved_state, best_generated_state, ptr_asa_out, OPTIONS)
     double (*user_cost_function) ();
     double *final_cost;
     double *parameter_initial_final;
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *curvature;
     double *maximum_tangent;
     double *current_cost_temperature;
     double *initial_user_parameter_temp;
     double *current_user_parameter_temp;
     double *accepted_to_generated_ratio;
     ALLOC_INT *number_parameters;
     int *parameter_type;
     int *valid_state_generated_flag;
     int *exit_status;
     ALLOC_INT *index_exit_v;
     ALLOC_INT *start_sequence;
     LONG_INT *number_accepted;
     LONG_INT *best_number_accepted_saved;
     LONG_INT *index_cost_acceptances;
     LONG_INT *number_generated;
     LONG_INT *number_invalid_generated_states;
     LONG_INT *index_parameter_generations;
     LONG_INT *best_number_generated_saved;
     STATE *current_generated_state;
     STATE *last_saved_state;
     STATE *best_generated_state;
     FILE *ptr_asa_out;
     USER_DEFINES *OPTIONS;
#endif
{
  ALLOC_INT index_v;            /* iteration index */
  int curvatureFlag;
  int exit_exit_status, tmp_locate;
#if MULTI_MIN
  int multi_index;
#endif

  tmp_locate = OPTIONS->Locate_Cost;
  if (tmp_locate) {             /* stop compiler warning */
    ;
  }

  exit_exit_status = 0;

  /* return final function minimum and associated parameters */
  *final_cost = best_generated_state->cost;
  VFOR (index_v) {
    parameter_initial_final[index_v] =
      best_generated_state->parameter[index_v];
  }

  OPTIONS->N_Accepted = *best_number_accepted_saved;
  OPTIONS->N_Generated = *best_number_generated_saved;

#if MULTI_MIN
  for (multi_index = OPTIONS->Multi_Number - 1; multi_index >= 0;
       --multi_index) {
    best_generated_state->cost = OPTIONS->Multi_Cost[multi_index];
    VFOR (index_v) {
      best_generated_state->parameter[index_v] =
        OPTIONS->Multi_Params[multi_index][index_v];
    }
#if ASA_PRINT
    fprintf (ptr_asa_out, "\n\t\t multi_index = %d\n", multi_index);
#endif /* ASA_PRINT */
#endif /* MULTI_MIN */
    if (*exit_status != TOO_MANY_INVALID_STATES
        && *exit_status != IMMEDIATE_EXIT
        && *exit_status != INVALID_USER_INPUT
        && *exit_status != INVALID_COST_FUNCTION
        && *exit_status != INVALID_COST_FUNCTION_DERIV) {
      if (OPTIONS->Curvature_0 != TRUE)
        OPTIONS->Locate_Cost = 5;       /* calc curvatures while exiting asa */

      /* calculate curvatures and tangents at best point */
      curvatureFlag = TRUE;
      cost_derivatives (user_cost_function,
                        parameter_minimum,
                        parameter_maximum,
                        tangents,
                        curvature,
                        maximum_tangent,
                        number_parameters,
                        parameter_type,
                        &exit_exit_status,
                        &curvatureFlag,
                        valid_state_generated_flag,
                        number_invalid_generated_states,
                        current_generated_state,
                        best_generated_state, ptr_asa_out, OPTIONS);
    }
#if ASA_PRINT
    if (exit_exit_status == INVALID_COST_FUNCTION_DERIV)
      fprintf (ptr_asa_out, "\n\n  in asa_exit: INVALID_COST_FUNCTION_DERIV");

    if (*exit_status != INVALID_USER_INPUT
        && *exit_status != INVALID_COST_FUNCTION
        && *exit_status != INVALID_COST_FUNCTION_DERIV)
      print_state (parameter_minimum,
                   parameter_maximum,
                   tangents,
                   curvature,
                   current_cost_temperature,
                   current_user_parameter_temp,
                   accepted_to_generated_ratio,
                   number_parameters,
                   &curvatureFlag,
                   number_accepted,
                   index_cost_acceptances,
                   number_generated,
                   number_invalid_generated_states,
                   last_saved_state,
                   best_generated_state, ptr_asa_out, OPTIONS);
#endif /* ASA_PRINT */

#if MULTI_MIN
  }
  best_generated_state->cost = OPTIONS->Multi_Cost[0];
  VFOR (index_v) {
    best_generated_state->parameter[index_v] =
      OPTIONS->Multi_Params[0][index_v];
  }
#endif /* MULTI_MIN */

#if ASA_PRINT
  switch (*exit_status) {
  case NORMAL_EXIT:
    fprintf (ptr_asa_out,
             "\n\n NORMAL_EXIT exit_status = %d\n", *exit_status);
    break;
  case P_TEMP_TOO_SMALL:
    fprintf (ptr_asa_out,
             "\n\n P_TEMP_TOO_SMALL exit_status = %d\n", *exit_status);
#if INT_ALLOC
  char *small_format = "current_user_parameter_temp[%d] too small = %*.*g\n";
#else
#if INT_LONG
  char *small_format = "current_user_parameter_temp[%ld] too small = %*.*g\n";
#else
  char *small_format = "current_user_parameter_temp[%d] too small = %*.*g\n";
#endif
#endif
    fprintf (ptr_asa_out,
             small_format,
             *index_exit_v,
             G_FIELD, G_PRECISION,
             current_user_parameter_temp[*index_exit_v]);
    break;
  case C_TEMP_TOO_SMALL:
    fprintf (ptr_asa_out,
             "\n\n C_TEMP_TOO_SMALL exit_status = %d\n", *exit_status);
    fprintf (ptr_asa_out,
             "*current_cost_temperature too small = %*.*g\n",
             G_FIELD, G_PRECISION, *current_cost_temperature);
    break;
  case COST_REPEATING:
    fprintf (ptr_asa_out,
             "\n\n COST_REPEATING exit_status = %d\n", *exit_status);
    break;
  case TOO_MANY_INVALID_STATES:
    fprintf (ptr_asa_out,
             "\n\n  TOO_MANY_INVALID_STATES exit_status = %d\n",
             *exit_status);
    break;
  case IMMEDIATE_EXIT:
    fprintf (ptr_asa_out,
             "\n\n  IMMEDIATE_EXIT exit_status = %d\n", *exit_status);
    break;
  case INVALID_USER_INPUT:
    fprintf (ptr_asa_out,
             "\n\n  INVALID_USER_INPUT exit_status = %d\n", *exit_status);
    break;
  case INVALID_COST_FUNCTION:
    fprintf (ptr_asa_out,
             "\n\n  INVALID_COST_FUNCTION exit_status = %d\n", *exit_status);
    break;
  case INVALID_COST_FUNCTION_DERIV:
    fprintf (ptr_asa_out,
             "\n\n  INVALID_COST_FUNCTION_DERIV exit_status = %d\n",
             *exit_status);
    break;
  default:
    fprintf (ptr_asa_out, "\n\n ERR: no exit code available = %d\n",
             *exit_status);
  }

  switch (OPTIONS->Locate_Cost) {
  case 0:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, initial cost temperature\n",
             OPTIONS->Locate_Cost);
    break;
  case 1:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, initial cost value\n", OPTIONS->Locate_Cost);
    break;
  case 2:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, new generated state\n",
             OPTIONS->Locate_Cost);
    break;
  case 12:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, new generated state just after a new best state\n",
             OPTIONS->Locate_Cost);
    break;
  case 3:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, cost derivatives, reannealing parameters\n",
             OPTIONS->Locate_Cost);
    break;
  case 4:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, reannealing cost temperature\n",
             OPTIONS->Locate_Cost);
    break;
  case 5:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, calculating curvatures while exiting asa ()\n",
             OPTIONS->Locate_Cost);
    break;
  case -1:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, exited main asa () loop by user-defined OPTIONS\n",
             OPTIONS->Locate_Cost);
    break;
  default:
    fprintf (ptr_asa_out,
             " Locate_Cost = %d, no index available for Locate_Cost\n",
             OPTIONS->Locate_Cost);
  }

  if (*exit_status != INVALID_USER_INPUT
      && *exit_status != INVALID_COST_FUNCTION
      && *exit_status != INVALID_COST_FUNCTION_DERIV) {
    fprintf (ptr_asa_out,
             "final_cost = best_generated_state->cost = %-*.*g\n",
             G_FIELD, G_PRECISION, *final_cost);
#if INT_LONG
    fprintf (ptr_asa_out,
             "*number_accepted at best_generated_state->cost = %ld\n",
             *best_number_accepted_saved);
    fprintf (ptr_asa_out,
             "*number_generated at best_generated_state->cost = %ld\n",
             *best_number_generated_saved);
#else
    fprintf (ptr_asa_out,
             "*number_accepted at best_generated_state->cost = %d\n",
             *best_number_accepted_saved);
    fprintf (ptr_asa_out,
             "*number_generated at best_generated_state->cost = %d\n",
             *best_number_generated_saved);
#endif
  }
#endif

#if ASA_TEMPLATE_SELFOPT
  if (OPTIONS->Asa_Data_Dbl[0] > (double) MIN_DOUBLE)
    OPTIONS->Asa_Data_Dbl[1] = (double) (*best_number_generated_saved);
#endif

  /* reset OPTIONS->Sequential_Parameters */
  OPTIONS->Sequential_Parameters = *start_sequence;

#if ASA_PRINT
#if TIME_CALC
  /* print ending time */
  print_time ("asa_end", ptr_asa_out);
#endif
  fprintf (ptr_asa_out, "\n\n\n");
  fflush (ptr_asa_out);
  fclose (ptr_asa_out);
#endif

  return (0);
}

/***********************************************************************
* generate_new_state
*       Generates a valid new state from the old state
***********************************************************************/
#if HAVE_ANSI
int

generate_new_state (double (*user_random_generator) (LONG_INT *),
                    LONG_INT * seed,
                    double *parameter_minimum,
                    double *parameter_maximum,
                    double *current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                    double *initial_user_parameter_temp,
                    double *temperature_scale_parameters,
#endif
                    ALLOC_INT * number_parameters,
                    int *parameter_type,
                    STATE * current_generated_state,
                    STATE * last_saved_state, USER_DEFINES * OPTIONS)
#else
int

generate_new_state (user_random_generator,
                    seed,
                    parameter_minimum,
                    parameter_maximum, current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                    initial_user_parameter_temp, temperature_scale_parameters,
#endif
                    number_parameters,
                    parameter_type,
                    current_generated_state, last_saved_state, OPTIONS)
     double (*user_random_generator) ();
     LONG_INT *seed;
     double *parameter_minimum;
     double *parameter_maximum;
     double *current_user_parameter_temp;
#if USER_GENERATING_FUNCTION
     double *initial_user_parameter_temp;
     double *temperature_scale_parameters;
#endif
     ALLOC_INT *number_parameters;
     int *parameter_type;
     STATE *current_generated_state;
     STATE *last_saved_state;
     USER_DEFINES *OPTIONS;
#endif
{
  ALLOC_INT index_v;
  double x;
  double parameter_v, min_parameter_v, max_parameter_v, temperature_v,
    parameter_range_v;
#if USER_GENERATING_FUNCTION
  double init_param_temp_v;
  double temp_scale_params_v;
#endif
#if ASA_RESOLUTION
  double xres, xint, xminus, xplus, dx, dxminus, dxplus;
#endif

  /* generate a new value for each parameter */
  VFOR (index_v) {
    if (OPTIONS->Sequential_Parameters >= -1) {
      ++OPTIONS->Sequential_Parameters;
      if (OPTIONS->Sequential_Parameters == *number_parameters)
        OPTIONS->Sequential_Parameters = 0;
      index_v = OPTIONS->Sequential_Parameters;
    }
    min_parameter_v = parameter_minimum[index_v];
    max_parameter_v = parameter_maximum[index_v];
    parameter_range_v = max_parameter_v - min_parameter_v;

    /* ignore parameters that have too small a range */
    if (fabs (parameter_range_v) < (double) EPS_DOUBLE)
      continue;

    temperature_v = current_user_parameter_temp[index_v];
#if USER_GENERATING_FUNCTION
    init_param_temp_v = initial_user_parameter_temp[index_v];
    temp_scale_params_v = temperature_scale_parameters[index_v];
#endif
    parameter_v = last_saved_state->parameter[index_v];

    /* Handle discrete parameters. */
#if ASA_RESOLUTION
    xres = OPTIONS->Coarse_Resolution[index_v];
    if (xres > EPS_DOUBLE) {
      min_parameter_v -= (xres / TWO);
      max_parameter_v += (xres / TWO);
      parameter_range_v = max_parameter_v - min_parameter_v;
    }
#endif /* ASA_RESOLUTION */
    if (INTEGER_PARAMETER (index_v)) {
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        ;
      } else {
#endif /* ASA_RESOLUTION */
        min_parameter_v -= HALF;
        max_parameter_v += HALF;
        parameter_range_v = max_parameter_v - min_parameter_v;
      }
#if ASA_RESOLUTION
    }
#endif

    /* generate a new state x within the parameter bounds */
    for (;;) {
#if USER_GENERATING_FUNCTION
      x = OPTIONS->Generating_Distrib (seed,
                                       number_parameters,
                                       index_v,
                                       temperature_v,
                                       init_param_temp_v,
                                       temp_scale_params_v,
                                       parameter_v,
                                       parameter_range_v,
                                       last_saved_state->parameter, OPTIONS);
#else
      x = parameter_v
        + generate_asa_state (user_random_generator, seed, &temperature_v)
        * parameter_range_v;
#endif /* USER_GENERATING_FUNCTION */
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        xint = xres * (double) ((LONG_INT) (x / xres));
        xplus = xint + xres;
        xminus = xint - xres;
        dx = fabs (xint - x);
        dxminus = fabs (xminus - x);
        dxplus = fabs (xplus - x);

        if (dx < dxminus && dx < dxplus)
          x = xint;
        else if (dxminus < dxplus)
          x = xminus;
        else
          x = xplus;
      }
#endif /* ASA_RESOLUTION */

      /* exit the loop if within its valid parameter range */
      if (x <= max_parameter_v - (double) EPS_DOUBLE
          && x >= min_parameter_v + (double) EPS_DOUBLE)
        break;
    }

    /* Handle discrete parameters.
       You might have to check rounding on your machine. */
    if (INTEGER_PARAMETER (index_v)) {
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        ;
      } else {
#endif /* ASA_RESOLUTION */
        if (x < min_parameter_v + HALF)
          x = min_parameter_v + HALF + (double) EPS_DOUBLE;
        if (x > max_parameter_v - HALF)
          x = max_parameter_v - HALF + (double) EPS_DOUBLE;

        if (x + HALF > ZERO) {
          x = (double) ((LONG_INT) (x + HALF));
        } else {
          x = (double) ((LONG_INT) (x - HALF));
        }
        if (x > parameter_maximum[index_v])
          x = parameter_maximum[index_v];
        if (x < parameter_minimum[index_v])
          x = parameter_minimum[index_v];
      }
#if ASA_RESOLUTION
    }
    if (xres > EPS_DOUBLE) {
      if (x < min_parameter_v + xres / TWO)
        x = min_parameter_v + xres / TWO + (double) EPS_DOUBLE;
      if (x > max_parameter_v - xres / TWO)
        x = max_parameter_v - xres / TWO + (double) EPS_DOUBLE;

      if (x > parameter_maximum[index_v])
        x = parameter_maximum[index_v];
      if (x < parameter_minimum[index_v])
        x = parameter_minimum[index_v];
    }
#endif /* ASA_RESOLUTION */

    /* save the newly generated value */
    current_generated_state->parameter[index_v] = x;

    if (OPTIONS->Sequential_Parameters >= 0)
      break;
  }
  return (0);
}

#if ASA_PARALLEL
/***********************************************************************
* generate_new_state_par
*       Generates a valid new state from the old state
***********************************************************************/
#if HAVE_ANSI
int

generate_new_state_par (double (*user_random_generator) (LONG_INT *),
                        LONG_INT * seed,
                        double *parameter_minimum,
                        double *parameter_maximum,
                        double *current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                        double *initial_user_parameter_temp,
                        double *temperature_scale_parameters,
#endif
                        ALLOC_INT * number_parameters,
                        int *parameter_type,
                        LONG_INT i_prll,
                        STATE * gener_block_state,
                        STATE * last_saved_state, USER_DEFINES * OPTIONS)
#else
int

generate_new_state_par (user_random_generator,
                        seed,
                        parameter_minimum,
                        parameter_maximum, current_user_parameter_temp,
#if USER_GENERATING_FUNCTION
                        initial_user_parameter_temp,
                        temperature_scale_parameters,
#endif
                        number_parameters,
                        parameter_type,
                        i_prll, gener_block_state, last_saved_state, OPTIONS)
     double (*user_random_generator) ();
     LONG_INT *seed;
     double *parameter_minimum;
     double *parameter_maximum;
     double *current_user_parameter_temp;
#if USER_GENERATING_FUNCTION
     double *initial_user_parameter_temp;
     double *temperature_scale_parameters;
#endif
     ALLOC_INT *number_parameters;
     int *parameter_type;
     LONG_INT i_prll;
     STATE *gener_block_state;
     STATE *last_saved_state;
     USER_DEFINES *OPTIONS;
#endif
{
  ALLOC_INT index_v;
  double x;
  double parameter_v, min_parameter_v, max_parameter_v, temperature_v,
    parameter_range_v;
#if USER_GENERATING_FUNCTION
  double init_param_temp_v;
  double temp_scale_params_v;
#endif
#if ASA_RESOLUTION
  double xres, xint, xminus, xplus, dx, dxminus, dxplus;
#endif

  /* generate a new value for each parameter */
  VFOR (index_v) {
    if (OPTIONS->Sequential_Parameters >= -1) {
      ++OPTIONS->Sequential_Parameters;
      if (OPTIONS->Sequential_Parameters == *number_parameters)
        OPTIONS->Sequential_Parameters = 0;
      index_v = OPTIONS->Sequential_Parameters;
    }
    min_parameter_v = parameter_minimum[index_v];
    max_parameter_v = parameter_maximum[index_v];
    parameter_range_v = max_parameter_v - min_parameter_v;

    /* ignore parameters that have too small a range */
    if (fabs (parameter_range_v) < (double) EPS_DOUBLE)
      continue;

    temperature_v = current_user_parameter_temp[index_v];
#if USER_GENERATING_FUNCTION
    init_param_temp_v = initial_user_parameter_temp[index_v];
    temp_scale_params_v = temperature_scale_parameters[index_v];
#endif
    parameter_v = last_saved_state->parameter[index_v];

    /* Handle discrete parameters. */
#if ASA_RESOLUTION
    xres = OPTIONS->Coarse_Resolution[index_v];
    if (xres > EPS_DOUBLE) {
      min_parameter_v -= (xres / TWO);
      max_parameter_v += (xres / TWO);
      parameter_range_v = max_parameter_v - min_parameter_v;
    }
#endif /* ASA_RESOLUTION */
    if (INTEGER_PARAMETER (index_v)) {
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        ;
      } else {
#endif /* ASA_RESOLUTION */
        min_parameter_v -= HALF;
        max_parameter_v += HALF;
        parameter_range_v = max_parameter_v - min_parameter_v;
      }
#if ASA_RESOLUTION
    }
#endif

    /* generate a new state x within the parameter bounds */
    for (;;) {
#if USER_GENERATING_FUNCTION
      x = OPTIONS->Generating_Distrib (seed,
                                       number_parameters,
                                       index_v,
                                       temperature_v,
                                       init_param_temp_v,
                                       temp_scale_params_v,
                                       parameter_v,
                                       parameter_range_v,
                                       last_saved_state->parameter, OPTIONS);
#else
      x = parameter_v
        + generate_asa_state (user_random_generator, seed, &temperature_v)
        * parameter_range_v;
#endif /* USER_GENERATING_FUNCTION */
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        xint = xres * (double) ((LONG_INT) (x / xres));
        xplus = xint + xres;
        xminus = xint - xres;
        dx = fabs (xint - x);
        dxminus = fabs (xminus - x);
        dxplus = fabs (xplus - x);

        if (dx < dxminus && dx < dxplus)
          x = xint;
        else if (dxminus < dxplus)
          x = xminus;
        else
          x = xplus;
      }
#endif /* ASA_RESOLUTION */

      /* exit the loop if within its valid parameter range */
      if (x <= max_parameter_v - (double) EPS_DOUBLE
          && x >= min_parameter_v + (double) EPS_DOUBLE)
        break;
    }

    /* Handle discrete parameters.
       You might have to check rounding on your machine. */
    if (INTEGER_PARAMETER (index_v)) {
#if ASA_RESOLUTION
      if (xres > EPS_DOUBLE) {
        ;
      } else {
#endif /* ASA_RESOLUTION */
        if (x < min_parameter_v + HALF)
          x = min_parameter_v + HALF + (double) EPS_DOUBLE;
        if (x > max_parameter_v - HALF)
          x = max_parameter_v - HALF + (double) EPS_DOUBLE;

        if (x + HALF > ZERO) {
          x = (double) ((LONG_INT) (x + HALF));
        } else {
          x = (double) ((LONG_INT) (x - HALF));
        }
        if (x > parameter_maximum[index_v])
          x = parameter_maximum[index_v];
        if (x < parameter_minimum[index_v])
          x = parameter_minimum[index_v];
      }
#if ASA_RESOLUTION
    }
    if (xres > EPS_DOUBLE) {
      if (x < min_parameter_v + xres / TWO)
        x = min_parameter_v + xres / TWO + (double) EPS_DOUBLE;
      if (x > max_parameter_v - xres / TWO)
        x = max_parameter_v - xres / TWO + (double) EPS_DOUBLE;

      if (x > parameter_maximum[index_v])
        x = parameter_maximum[index_v];
      if (x < parameter_minimum[index_v])
        x = parameter_minimum[index_v];
    }
#endif /* ASA_RESOLUTION */

    /* save the newly generated value */
    gener_block_state[i_prll].parameter[index_v] = x;

    if (OPTIONS->Sequential_Parameters >= 0)
      break;
  }
  return (0);
}

#endif /* ASA_PARALLEL */

/***********************************************************************
* generate_asa_state
*       This function generates a single value according to the
*       ASA generating function and the passed temperature
***********************************************************************/
#if HAVE_ANSI
double

generate_asa_state (double (*user_random_generator) (LONG_INT *),
                    LONG_INT * seed, double *temp)
#else
double
generate_asa_state (user_random_generator, seed, temp)
     double (*user_random_generator) ();
     LONG_INT *seed;
     double *temp;
#endif
{
  double x, y, z;

  x = (*user_random_generator) (seed);
  y = x < HALF ? -ONE : ONE;
  z = y * *temp * (F_POW ((ONE + ONE / *temp), fabs (TWO * x - ONE)) - ONE);

  return (z);

}

/***********************************************************************
* accept_new_state
*	This procedure accepts or rejects a newly generated state,
*	depending on whether the difference between new and old
*	cost functions passes a statistical test. If accepted,
*	the current state is updated.
***********************************************************************/
#if HAVE_ANSI
void

accept_new_state (double (*user_random_generator) (LONG_INT *),
                  LONG_INT * seed,
                  double *parameter_minimum,
                  double *parameter_maximum, double *current_cost_temperature,
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
                  STATE * current_generated_state, STATE * last_saved_state,
#if ASA_SAMPLE
                  FILE * ptr_asa_out,
#endif
                  USER_DEFINES * OPTIONS)
#else
void

accept_new_state (user_random_generator,
                  seed,
                  parameter_minimum,
                  parameter_maximum, current_cost_temperature,
#if ASA_SAMPLE
                  current_user_parameter_temp,
#endif
                  number_parameters,
                  recent_number_acceptances,
                  number_accepted,
                  index_cost_acceptances,
                  number_acceptances_saved,
                  recent_number_generated,
                  number_generated,
                  index_parameter_generations,
                  current_generated_state, last_saved_state,
#if ASA_SAMPLE
                  ptr_asa_out,
#endif
                  OPTIONS)
     double (*user_random_generator) ();
     LONG_INT *seed;
     double *parameter_minimum;
     double *parameter_maximum;
     double *current_cost_temperature;
#if ASA_SAMPLE
     double *current_user_parameter_temp;
#endif
     ALLOC_INT *number_parameters;
     LONG_INT *recent_number_acceptances;
     LONG_INT *number_accepted;
     LONG_INT *index_cost_acceptances;
     LONG_INT *number_acceptances_saved;
     LONG_INT *recent_number_generated;
     LONG_INT *number_generated;
     LONG_INT *index_parameter_generations;
     STATE *current_generated_state;
     STATE *last_saved_state;
#if ASA_SAMPLE
     FILE *ptr_asa_out;
#endif
     USER_DEFINES *OPTIONS;

#endif
{
#if USER_ACCEPTANCE_TEST
#else
  double delta_cost;
#if USER_ACCEPT_ASYMP_EXP
  double q;
#endif
#endif
  double prob_test, unif_test;
  double curr_cost_temp;
  ALLOC_INT index_v;
#if ASA_SAMPLE
  LONG_INT active_params;
  double weight_param_ind, weight_aver, range;
#endif

  /* update accepted and generated count */
  ++*number_acceptances_saved;
  ++*recent_number_generated;
  ++*number_generated;
  OPTIONS->N_Generated = *number_generated;

  /* increment the parameter index generation for each parameter */
  if (OPTIONS->Sequential_Parameters >= 0) {
    /* ignore parameters with too small a range */
    if (!PARAMETER_RANGE_TOO_SMALL (OPTIONS->Sequential_Parameters))
      ++(index_parameter_generations[OPTIONS->Sequential_Parameters]);
  } else {
    VFOR (index_v) {
      if (!PARAMETER_RANGE_TOO_SMALL (index_v))
        ++(index_parameter_generations[index_v]);
    }
  }

  /* effective cost function for testing acceptance criteria,
     calculate the cost difference and divide by the temperature */
  curr_cost_temp = *current_cost_temperature;
#if USER_ACCEPTANCE_TEST
  if (OPTIONS->Cost_Acceptance_Flag == TRUE) {
    if (OPTIONS->User_Acceptance_Flag == TRUE) {
      unif_test = ZERO;
      OPTIONS->User_Acceptance_Flag = FALSE;
      OPTIONS->Cost_Acceptance_Flag = FALSE;
    } else {
      unif_test = ONE;
      OPTIONS->Cost_Acceptance_Flag = FALSE;
    }
  } else {
    OPTIONS->Acceptance_Test (current_generated_state->cost,
                              parameter_minimum,
                              parameter_maximum, number_parameters, OPTIONS);
    if (OPTIONS->User_Acceptance_Flag == TRUE) {
      unif_test = ZERO;
      OPTIONS->User_Acceptance_Flag = FALSE;
    } else {
      unif_test = ONE;
    }
  }
  prob_test = OPTIONS->Prob_Bias;
#else /* USER_ACCEPTANCE_TEST */

#if USER_COST_SCHEDULE
  curr_cost_temp =
    (OPTIONS->Cost_Schedule (*current_cost_temperature, OPTIONS)
     + (double) EPS_DOUBLE);
#endif
  delta_cost = (current_generated_state->cost - last_saved_state->cost)
    / (curr_cost_temp + (double) EPS_DOUBLE);

#if USER_ACCEPT_ASYMP_EXP
  q = OPTIONS->Asymp_Exp_Param;
  if (fabs (ONE - q) < (double) EPS_DOUBLE)
    prob_test = MIN (ONE, (F_EXP (EXPONENT_CHECK (-delta_cost))));
  else if ((ONE - (ONE - q) * delta_cost) < (double) EPS_DOUBLE)
    prob_test = MIN (ONE, (F_EXP (EXPONENT_CHECK (-delta_cost))));
  else
    prob_test = MIN (ONE, F_POW ((ONE - (ONE - q) * delta_cost),
                                 (ONE / (ONE - q))));
#else /* USER_ACCEPT_ASYMP_EXP */

#if USER_ACCEPT_THRESHOLD       /* USER_ACCEPT_THRESHOLD */
  prob_test = delta_cost <= 1.0 ? 1.0 : 0.0;
#else /* Metropolis */
  prob_test = MIN (ONE, (F_EXP (EXPONENT_CHECK (-delta_cost))));
#endif /* USER_ACCEPT_THRESHOLD */

#endif /* USER_ACCEPT_ASYMP_EXP */

  unif_test = (*user_random_generator) (seed);
#endif /* USER_ACCEPTANCE_TEST */

#if ASA_SAMPLE
  active_params = 0;
  weight_aver = ZERO;
  VFOR (index_v) {
    /* ignore parameters with too small a range */
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
    ++active_params;
    range = parameter_maximum[index_v] - parameter_minimum[index_v];
    weight_param_ind = TWO * (fabs ((last_saved_state->parameter[index_v]
                                     -
                                     current_generated_state->parameter
                                     [index_v]) / range)
                              + current_user_parameter_temp[index_v])
      * F_LOG (ONE + ONE / current_user_parameter_temp[index_v]);
    weight_aver += weight_param_ind;
    OPTIONS->Bias_Generated[index_v] = ONE / weight_param_ind;
  }
  weight_aver /= (double) active_params;
  OPTIONS->Average_Weights = weight_aver;
  if (prob_test >= unif_test) {
    OPTIONS->Bias_Acceptance = prob_test;
  } else {
    OPTIONS->Bias_Acceptance = ONE - prob_test;
  }

#if ASA_PRINT
  if (OPTIONS->Limit_Weights < OPTIONS->Average_Weights) {
    fprintf (ptr_asa_out, ":SAMPLE#\n");
    if (prob_test >= unif_test) {
      fprintf (ptr_asa_out,
#if INT_LONG
               ":SAMPLE+ %10ld %*.*g %*.*g %*.*g %*.*g\n",
#else
               ":SAMPLE+ %10d %*.*g %*.*g %*.*g\n",
#endif
               OPTIONS->N_Accepted,
               G_FIELD, G_PRECISION, current_generated_state->cost,
               G_FIELD, G_PRECISION, *current_cost_temperature,
               G_FIELD, G_PRECISION, OPTIONS->Bias_Acceptance,
               G_FIELD, G_PRECISION, OPTIONS->Average_Weights);
      VFOR (index_v) {
        /* ignore parameters with too small a range */
        if (PARAMETER_RANGE_TOO_SMALL (index_v))
          continue;
        range = parameter_maximum[index_v] - parameter_minimum[index_v];
        fprintf (ptr_asa_out,
#if INT_ALLOC
                 ":SAMPLE %11d %*.*g %*.*g %*.*g %*.*g\n",
#else
#if INT_LONG
                 ":SAMPLE %11ld %*.*g %*.*g %*.*g %*.*g\n",
#else
                 ":SAMPLE %11d %*.*g %*.*g %*.*g %*.*g\n",
#endif
#endif
                 index_v,
                 G_FIELD, G_PRECISION,
                 current_generated_state->parameter[index_v], G_FIELD,
                 G_PRECISION, current_user_parameter_temp[index_v],
                 G_FIELD, G_PRECISION, OPTIONS->Bias_Generated[index_v],
                 G_FIELD, G_PRECISION, range);
      }
    } else {
      fprintf (ptr_asa_out,
#if INT_LONG
               ":SAMPLE %11ld %*.*g %*.*g %*.*g %*.*g\n",
#else
               ":SAMPLE %11d %*.*g %*.*g %*.*g\n",
#endif
               OPTIONS->N_Accepted,
               G_FIELD, G_PRECISION, last_saved_state->cost,
               G_FIELD, G_PRECISION, *current_cost_temperature,
               G_FIELD, G_PRECISION, OPTIONS->Bias_Acceptance,
               G_FIELD, G_PRECISION, OPTIONS->Average_Weights);
      VFOR (index_v) {
        /* ignore parameters with too small a range */
        if (PARAMETER_RANGE_TOO_SMALL (index_v))
          continue;
        range = parameter_maximum[index_v] - parameter_minimum[index_v];
        fprintf (ptr_asa_out,
#if INT_ALLOC
                 ":SAMPLE %11d %*.*g %*.*g %*.*g %*.*g\n",
#else
#if INT_LONG
                 ":SAMPLE %11ld %*.*g %*.*g %*.*g %*.*g\n",
#else
                 ":SAMPLE %11d %*.*g %*.*g %*.*g %*.*g\n",
#endif
#endif
                 index_v,
                 G_FIELD, G_PRECISION,
                 last_saved_state->parameter[index_v], G_FIELD,
                 G_PRECISION, current_user_parameter_temp[index_v],
                 G_FIELD, G_PRECISION, OPTIONS->Bias_Generated[index_v],
                 G_FIELD, G_PRECISION, range);
      }
    }
  }
#endif
#endif /* ASA_SAMPLE */

  /* accept/reject the new state */
  if (prob_test >= unif_test) {
    /* copy current state to the last saved state */

    last_saved_state->cost = current_generated_state->cost;
    VFOR (index_v) {
      /* ignore parameters with too small a range */
      if (PARAMETER_RANGE_TOO_SMALL (index_v))
        continue;
      last_saved_state->parameter[index_v] =
        current_generated_state->parameter[index_v];
    }

    /* update acceptance counts */
    ++*recent_number_acceptances;
    ++*number_accepted;
    ++*index_cost_acceptances;
    *number_acceptances_saved = *number_accepted;
    OPTIONS->N_Accepted = *number_accepted;
  }
}

/***********************************************************************
* reanneal
*	Readjust temperatures of generating and acceptance functions
***********************************************************************/
#if HAVE_ANSI
void

reanneal (double *parameter_minimum,
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
          STATE * best_generated_state, USER_DEFINES * OPTIONS)
#else
void

reanneal (parameter_minimum,
          parameter_maximum,
          tangents,
          maximum_tangent,
          current_cost_temperature,
          initial_cost_temperature,
          temperature_scale_cost,
          current_user_parameter_temp,
          initial_user_parameter_temp,
          temperature_scale_parameters,
          number_parameters,
          parameter_type,
          index_cost_acceptances,
          index_parameter_generations,
          last_saved_state, best_generated_state, OPTIONS)
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *maximum_tangent;
     double *current_cost_temperature;
     double *initial_cost_temperature;
     double *temperature_scale_cost;
     double *current_user_parameter_temp;
     double *initial_user_parameter_temp;
     double *temperature_scale_parameters;
     ALLOC_INT *number_parameters;
     int *parameter_type;
     LONG_INT *index_cost_acceptances;
     LONG_INT *index_parameter_generations;
     STATE *last_saved_state;
     STATE *best_generated_state;
     USER_DEFINES *OPTIONS;
#endif
{
  ALLOC_INT index_v;
  int cost_test;
  double tmp_var_db3;
  double new_temperature;
  double log_new_temperature_ratio;
  double log_init_cur_temp_ratio;
  double temperature_rescale_power;
  double cost_best, cost_last;
  double tmp_dbl, tmp_dbl1;

  double xnumber_parameters[1];

  cost_test = cost_function_test (last_saved_state->cost,
                                  last_saved_state->parameter,
                                  parameter_minimum,
                                  parameter_maximum, number_parameters,
                                  xnumber_parameters);

  if (OPTIONS->Reanneal_Parameters == TRUE) {
    VFOR (index_v) {
      if (NO_REANNEAL (index_v))
        continue;

      /* use the temp double to prevent overflow */
      tmp_dbl = (double) index_parameter_generations[index_v];

      /* skip parameters with too small range or integer parameters */
      if (OPTIONS->Include_Integer_Parameters == TRUE) {
        if (PARAMETER_RANGE_TOO_SMALL (index_v))
          continue;
      } else {
        if (PARAMETER_RANGE_TOO_SMALL (index_v) ||
            INTEGER_PARAMETER (index_v))
          continue;
      }

      /* ignore parameters with too small tangents */
      if (fabs (tangents[index_v]) < (double) EPS_DOUBLE)
        continue;

      /* reset the index of parameter generations appropriately */
#if USER_REANNEAL_PARAMETERS
      new_temperature =
        fabs (OPTIONS->Reanneal_Params_Function (current_user_parameter_temp
                                                 [index_v], tangents[index_v],
                                                 *maximum_tangent, OPTIONS));
#else
      new_temperature =
        fabs (FUNCTION_REANNEAL_PARAMS
              (current_user_parameter_temp[index_v], tangents[index_v],
               *maximum_tangent));
#endif
      if (new_temperature < initial_user_parameter_temp[index_v]) {
        log_init_cur_temp_ratio =
          fabs (F_LOG (((double) EPS_DOUBLE
                        + initial_user_parameter_temp[index_v])
                       / ((double) EPS_DOUBLE + new_temperature)));
        tmp_dbl = (double) EPS_DOUBLE
          + F_POW (log_init_cur_temp_ratio
                   / temperature_scale_parameters[index_v],
                   *xnumber_parameters
#if QUENCH_PARAMETERS
                   / OPTIONS->User_Quench_Param_Scale[index_v]);
#else
          );
#endif
      } else {
        tmp_dbl = ONE;
      }

      /* Reset index_parameter_generations if index reset too large,
         and also reset the initial_user_parameter_temp, to achieve
         the same new temperature. */
      while (tmp_dbl > ((double) MAXIMUM_REANNEAL_INDEX)) {
        log_new_temperature_ratio =
          -temperature_scale_parameters[index_v] * F_POW (tmp_dbl,
#if QUENCH_PARAMETERS
                                                          OPTIONS->User_Quench_Param_Scale
                                                          [index_v]
#else
                                                          ONE
#endif
                                                          /
                                                          *xnumber_parameters);
        log_new_temperature_ratio =
          EXPONENT_CHECK (log_new_temperature_ratio);
        new_temperature =
          initial_user_parameter_temp[index_v] *
          F_EXP (log_new_temperature_ratio);
        tmp_dbl /= (double) REANNEAL_SCALE;
        temperature_rescale_power = ONE / F_POW ((double) REANNEAL_SCALE,
#if QUENCH_PARAMETERS
                                                 OPTIONS->User_Quench_Param_Scale
                                                 [index_v]
#else
                                                 ONE
#endif
                                                 / *xnumber_parameters);
        initial_user_parameter_temp[index_v] =
          new_temperature * F_POW (initial_user_parameter_temp[index_v] /
                                   new_temperature,
                                   temperature_rescale_power);
      }
      /* restore from temporary double */
      index_parameter_generations[index_v] = (LONG_INT) tmp_dbl;
    }
  }

  if (OPTIONS->Reanneal_Cost == 0) {
    ;
  } else if (OPTIONS->Reanneal_Cost < -1) {
    *index_cost_acceptances = 1;
  } else {
    /* reanneal : Reset the current cost temp and rescale the
       index of cost acceptances. */

    cost_best = best_generated_state->cost;
    cost_last = last_saved_state->cost;
#if USER_REANNEAL_COST
    cost_test = OPTIONS->Reanneal_Cost_Function (&cost_best,
                                                 &cost_last,
                                                 initial_cost_temperature,
                                                 current_cost_temperature,
                                                 OPTIONS);
    tmp_dbl1 = *current_cost_temperature;
#else
    cost_test = TRUE;
    if (OPTIONS->Reanneal_Cost == 1) {
      /* (re)set the initial cost_temperature */
      tmp_dbl = MAX (fabs (cost_last), fabs (cost_best));
      tmp_dbl = MAX (tmp_dbl, fabs (cost_best - cost_last));
      tmp_dbl = MAX ((double) EPS_DOUBLE, tmp_dbl);
      *initial_cost_temperature = MIN (*initial_cost_temperature, tmp_dbl);
    }

    tmp_dbl = (double) *index_cost_acceptances;

    tmp_dbl1 = MAX (fabs (cost_last - cost_best), *current_cost_temperature);
    tmp_dbl1 = MAX ((double) EPS_DOUBLE, tmp_dbl1);
    tmp_dbl1 = MIN (tmp_dbl1, *initial_cost_temperature);
#endif /* USER_REANNEAL_COST */
    if (cost_test == TRUE && (*current_cost_temperature > tmp_dbl1)) {
      tmp_var_db3 =
        fabs (F_LOG (((double) EPS_DOUBLE + *initial_cost_temperature) /
                     (tmp_dbl1)));
      tmp_dbl = (double) EPS_DOUBLE + F_POW (tmp_var_db3
                                             / *temperature_scale_cost,
                                             *xnumber_parameters
#if QUENCH_COST
                                             /
                                             OPTIONS->User_Quench_Cost_Scale
                                             [0]);
#else
        );
#endif
    } else {
      log_init_cur_temp_ratio =
        fabs (F_LOG (((double) EPS_DOUBLE + *initial_cost_temperature) /
                     ((double) EPS_DOUBLE + *current_cost_temperature)));
      tmp_dbl = (double) EPS_DOUBLE
        + F_POW (log_init_cur_temp_ratio
                 / *temperature_scale_cost, *xnumber_parameters
#if QUENCH_COST
                 / OPTIONS->User_Quench_Cost_Scale[0]
#else
#endif
        );
    }

    /* reset index_cost_temperature if index reset too large */
    while (tmp_dbl > ((double) MAXIMUM_REANNEAL_INDEX)) {
      log_new_temperature_ratio = -*temperature_scale_cost * F_POW (tmp_dbl,
#if QUENCH_COST
                                                                    OPTIONS->User_Quench_Cost_Scale
                                                                    [0]
#else
                                                                    ONE
#endif
                                                                    /
                                                                    *xnumber_parameters);
      log_new_temperature_ratio = EXPONENT_CHECK (log_new_temperature_ratio);
      new_temperature =
        *initial_cost_temperature * F_EXP (log_new_temperature_ratio);
      tmp_dbl /= (double) REANNEAL_SCALE;
      temperature_rescale_power = ONE / F_POW ((double) REANNEAL_SCALE,
#if QUENCH_COST
                                               OPTIONS->User_Quench_Cost_Scale
                                               [0]
#else
                                               ONE
#endif
                                               / *xnumber_parameters);
      *initial_cost_temperature =
        new_temperature * F_POW (*initial_cost_temperature /
                                 new_temperature, temperature_rescale_power);
    }
    *index_cost_acceptances = (LONG_INT) tmp_dbl;
#if USER_ACCEPTANCE_TEST
    OPTIONS->Cost_Temp_Init = *initial_cost_temperature;
#endif
  }
}

/***********************************************************************
* cost_derivatives
*	This procedure calculates the derivatives of the cost function
*	with respect to its parameters.  The first derivatives are
*	used as a sensitivity measure for reannealing.  The second
*	derivatives are calculated only if *curvature_flag=TRUE;
*	these are a measure of the covariance of the fit when a
*	minimum is found.
***********************************************************************/
  /* Calculate the numerical derivatives of the best
     generated state found so far */

  /* Assuming no information is given about the metric of the parameter
     space, use simple Cartesian space to calculate curvatures. */

#if HAVE_ANSI
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
                  USER_DEFINES * OPTIONS)
#else
void

cost_derivatives (user_cost_function,
                  parameter_minimum,
                  parameter_maximum,
                  tangents,
                  curvature,
                  maximum_tangent,
                  number_parameters,
                  parameter_type,
                  exit_status,
                  curvature_flag,
                  valid_state_generated_flag,
                  number_invalid_generated_states,
                  current_generated_state,
                  best_generated_state, ptr_asa_out, OPTIONS)
     double (*user_cost_function) ();
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *curvature;
     double *maximum_tangent;
     ALLOC_INT *number_parameters;
     int *parameter_type;
     int *exit_status;
     int *curvature_flag;
     int *valid_state_generated_flag;
     LONG_INT *number_invalid_generated_states;
     STATE *current_generated_state;
     STATE *best_generated_state;
     FILE *ptr_asa_out;
     USER_DEFINES *OPTIONS;
#endif
{
  ALLOC_INT index_v, index_vv, index_v_vv, index_vv_v;
  LONG_INT saved_num_invalid_gen_states;
#if ASA_PRINT
  LONG_INT tmp_saved;
#endif
  double parameter_v, parameter_vv, parameter_v_offset, parameter_vv_offset;
  double recent_best_cost;
  double new_cost_state_1, new_cost_state_2, new_cost_state_3;
  double delta_parameter_v, delta_parameter_vv;
  int immediate_flag;

  if (OPTIONS->Curvature_0 == TRUE)
    *curvature_flag = FALSE;
  if (OPTIONS->Curvature_0 == -1)
    *curvature_flag = TRUE;

  /* save Immediate_Exit flag */
  immediate_flag = OPTIONS->Immediate_Exit;

  /* save the best cost */
  recent_best_cost = best_generated_state->cost;

  /* copy the best state into the current state */
  VFOR (index_v) {
    /* ignore parameters with too small ranges */
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
    current_generated_state->parameter[index_v] =
      best_generated_state->parameter[index_v];
  }

  saved_num_invalid_gen_states = (*number_invalid_generated_states);

  /* set parameters (& possibly constraints) to best state */
  *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
  OPTIONS->User_Acceptance_Flag = TRUE;
  OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
  current_generated_state->cost =
    user_cost_function (current_generated_state->parameter,
                        parameter_minimum,
                        parameter_maximum,
                        tangents,
                        curvature,
                        number_parameters,
                        parameter_type,
                        valid_state_generated_flag, exit_status, OPTIONS);
  if ((*valid_state_generated_flag == FALSE)
      || ((current_generated_state->cost) != (current_generated_state->cost))
      || current_generated_state->cost < -MAX_DOUBLE
      || current_generated_state->cost > MAX_DOUBLE) {
    *exit_status = INVALID_COST_FUNCTION_DERIV;
    return;
  }
  if (*valid_state_generated_flag == FALSE)
    ++(*number_invalid_generated_states);

  if (OPTIONS->User_Tangents == TRUE) {
    *valid_state_generated_flag = -1;
#if USER_ACCEPTANCE_TEST
    OPTIONS->User_Acceptance_Flag = TRUE;
    OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
    current_generated_state->cost =
      user_cost_function (current_generated_state->parameter,
                          parameter_minimum,
                          parameter_maximum,
                          tangents,
                          curvature,
                          number_parameters,
                          parameter_type,
                          valid_state_generated_flag, exit_status, OPTIONS);
    if ((*valid_state_generated_flag == FALSE)
        || ((current_generated_state->cost) !=
            (current_generated_state->cost))
        || current_generated_state->cost < -MAX_DOUBLE
        || current_generated_state->cost > MAX_DOUBLE) {
      *exit_status = INVALID_COST_FUNCTION_DERIV;
      return;
    }
    if (*valid_state_generated_flag == FALSE)
      ++(*number_invalid_generated_states);
  } else {
    /* calculate tangents */
    VFOR (index_v) {
      if (NO_REANNEAL (index_v)) {
        tangents[index_v] = ZERO;
        continue;
      }
      /* skip parameters with too small range or integer parameters */
      if (OPTIONS->Include_Integer_Parameters == TRUE) {
        if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
          tangents[index_v] = ZERO;
          continue;
        }
      } else {
        if (PARAMETER_RANGE_TOO_SMALL (index_v) ||
            INTEGER_PARAMETER (index_v)) {
          tangents[index_v] = ZERO;
          continue;
        }
      }
#if DELTA_PARAMETERS
      delta_parameter_v = OPTIONS->User_Delta_Parameter[index_v];
#else
      delta_parameter_v = OPTIONS->Delta_X;
#endif
      if (delta_parameter_v < SMALL_FLOAT) {
        tangents[index_v] = 0;
        continue;
      }

      /* save the v_th parameter and delta_parameter */
      parameter_v = best_generated_state->parameter[index_v];

      parameter_v_offset = (ONE + delta_parameter_v) * parameter_v;
      if (parameter_v_offset > parameter_maximum[index_v] ||
          parameter_v_offset < parameter_minimum[index_v]) {
        delta_parameter_v = -delta_parameter_v;
        parameter_v_offset = (ONE + delta_parameter_v) * parameter_v;
      }

      /* generate the first sample point */
      current_generated_state->parameter[index_v] = parameter_v_offset;
      *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
      OPTIONS->User_Acceptance_Flag = TRUE;
      OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
      current_generated_state->cost =
        user_cost_function (current_generated_state->parameter,
                            parameter_minimum,
                            parameter_maximum,
                            tangents,
                            curvature,
                            number_parameters,
                            parameter_type,
                            valid_state_generated_flag, exit_status, OPTIONS);
      if ((*valid_state_generated_flag == FALSE)
          || ((current_generated_state->cost) !=
              (current_generated_state->cost))
          || current_generated_state->cost < -MAX_DOUBLE
          || current_generated_state->cost > MAX_DOUBLE) {
        *exit_status = INVALID_COST_FUNCTION_DERIV;
        return;
      }
      if (*valid_state_generated_flag == FALSE)
        ++(*number_invalid_generated_states);
      new_cost_state_1 = current_generated_state->cost;

      /* restore the parameter state */
      current_generated_state->parameter[index_v] = parameter_v;

      /* calculate the numerical derivative */
      tangents[index_v] = (new_cost_state_1 - recent_best_cost)
        / (delta_parameter_v * parameter_v + (double) EPS_DOUBLE);

    }
  }

  /* find the maximum |tangent| from all tangents */
  *maximum_tangent = 0;
  VFOR (index_v) {
    if (NO_REANNEAL (index_v))
      continue;

    /* ignore too small ranges and integer parameters types */
    if (OPTIONS->Include_Integer_Parameters == TRUE) {
      if (PARAMETER_RANGE_TOO_SMALL (index_v))
        continue;
    } else {
      if (PARAMETER_RANGE_TOO_SMALL (index_v)
          || INTEGER_PARAMETER (index_v))
        continue;
    }

    /* find the maximum |tangent| (from all tangents) */
    if (fabs (tangents[index_v]) > *maximum_tangent) {
      *maximum_tangent = fabs (tangents[index_v]);
    }
  }

  if (*curvature_flag == TRUE || *curvature_flag == -1) {
    /* calculate diagonal curvatures */
    VFOR (index_v) {
      /* index_v_vv: row index_v, column index_v */
      index_v_vv = ROW_COL_INDEX (index_v, index_v);

      if (NO_REANNEAL (index_v)) {
        curvature[index_v_vv] = ZERO;
        continue;
      }
      /* skip parameters with too small range or integer parameters */
      if (OPTIONS->Include_Integer_Parameters == TRUE) {
        if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
          curvature[index_v_vv] = ZERO;
          continue;
        }
      } else {
        if (PARAMETER_RANGE_TOO_SMALL (index_v) ||
            INTEGER_PARAMETER (index_v)) {
          curvature[index_v_vv] = ZERO;
          continue;
        }
      }
#if DELTA_PARAMETERS
      delta_parameter_v = OPTIONS->User_Delta_Parameter[index_v];
#else
      delta_parameter_v = OPTIONS->Delta_X;
#endif
      if (delta_parameter_v < SMALL_FLOAT) {
        curvature[index_v_vv] = ZERO;
        continue;
      }

      /* save the v_th parameter and delta_parameter */
      parameter_v = best_generated_state->parameter[index_v];

      if (parameter_v + delta_parameter_v * fabs (parameter_v)
          > parameter_maximum[index_v]) {
        /* generate the first sample point */
        current_generated_state->parameter[index_v] =
          parameter_v - TWO * delta_parameter_v * fabs (parameter_v);
        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_1 = current_generated_state->cost;

        /* generate the second sample point */
        current_generated_state->parameter[index_v] =
          parameter_v - delta_parameter_v * fabs (parameter_v);

        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_2 = current_generated_state->cost;

        /* restore the parameter state */
        current_generated_state->parameter[index_v] = parameter_v;

        /* calculate and store the curvature */
        curvature[index_v_vv] =
          (recent_best_cost - TWO * new_cost_state_2
           + new_cost_state_1) / (delta_parameter_v * delta_parameter_v
                                  * parameter_v * parameter_v +
                                  (double) EPS_DOUBLE);
      } else if (parameter_v - delta_parameter_v * fabs (parameter_v)
                 < parameter_minimum[index_v]) {
        /* generate the first sample point */
        current_generated_state->parameter[index_v] =
          parameter_v + TWO * delta_parameter_v * fabs (parameter_v);
        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_1 = current_generated_state->cost;

        /* generate the second sample point */
        current_generated_state->parameter[index_v] =
          parameter_v + delta_parameter_v * fabs (parameter_v);

        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_2 = current_generated_state->cost;

        /* restore the parameter state */
        current_generated_state->parameter[index_v] = parameter_v;

        /* index_v_vv: row index_v, column index_v */
        index_v_vv = ROW_COL_INDEX (index_v, index_v);

        /* calculate and store the curvature */
        curvature[index_v_vv] =
          (recent_best_cost - TWO * new_cost_state_2
           + new_cost_state_1) / (delta_parameter_v * delta_parameter_v
                                  * parameter_v * parameter_v +
                                  (double) EPS_DOUBLE);
      } else {
        /* generate the first sample point */
        parameter_v_offset = (ONE + delta_parameter_v) * parameter_v;
        current_generated_state->parameter[index_v] = parameter_v_offset;
        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_1 = current_generated_state->cost;

        /* generate the second sample point */
        current_generated_state->parameter[index_v] =
          (ONE - delta_parameter_v) * parameter_v;

        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_2 = current_generated_state->cost;

        /* restore the parameter state */
        current_generated_state->parameter[index_v] = parameter_v;

        /* calculate and store the curvature */
        curvature[index_v_vv] =
          (new_cost_state_2 - TWO * recent_best_cost
           + new_cost_state_1) / (delta_parameter_v * delta_parameter_v
                                  * parameter_v * parameter_v +
                                  (double) EPS_DOUBLE);
      }
    }

    /* calculate off-diagonal curvatures */
    VFOR (index_v) {
#if DELTA_PARAMETERS
      delta_parameter_v = OPTIONS->User_Delta_Parameter[index_v];
#else
      delta_parameter_v = OPTIONS->Delta_X;
#endif
      if (delta_parameter_v < SMALL_FLOAT) {
        VFOR (index_vv) {
          /* index_v_vv: row index_v, column index_vv */
          index_v_vv = ROW_COL_INDEX (index_v, index_vv);
          index_vv_v = ROW_COL_INDEX (index_vv, index_v);
          curvature[index_vv_v] = curvature[index_v_vv] = ZERO;
        }
        continue;
      }

      /* save the v_th parameter and delta_x */
      parameter_v = current_generated_state->parameter[index_v];

      VFOR (index_vv) {
        /* index_v_vv: row index_v, column index_vv */
        index_v_vv = ROW_COL_INDEX (index_v, index_vv);
        index_vv_v = ROW_COL_INDEX (index_vv, index_v);

        if (NO_REANNEAL (index_vv) || NO_REANNEAL (index_v)) {
          curvature[index_vv_v] = curvature[index_v_vv] = ZERO;
          continue;
        }
        /* calculate only the upper diagonal */
        if (index_v <= index_vv) {
          continue;
        }
        /* skip parms with too small range or integer parameters */
        if (OPTIONS->Include_Integer_Parameters == TRUE) {
          if (PARAMETER_RANGE_TOO_SMALL (index_v) ||
              PARAMETER_RANGE_TOO_SMALL (index_vv)) {
            curvature[index_vv_v] = curvature[index_v_vv] = ZERO;
            continue;
          }
        } else {
          if (INTEGER_PARAMETER (index_v) ||
              INTEGER_PARAMETER (index_vv) ||
              PARAMETER_RANGE_TOO_SMALL (index_v) ||
              PARAMETER_RANGE_TOO_SMALL (index_vv)) {
            curvature[index_vv_v] = curvature[index_v_vv] = ZERO;
            continue;
          }
        }
#if DELTA_PARAMETERS
        delta_parameter_vv = OPTIONS->User_Delta_Parameter[index_vv];
#else
        delta_parameter_vv = OPTIONS->Delta_X;
#endif
        if (delta_parameter_vv < SMALL_FLOAT) {
          curvature[index_vv_v] = curvature[index_v_vv] = ZERO;
          continue;
        }

        /* save the vv_th parameter and delta_parameter */
        parameter_vv = current_generated_state->parameter[index_vv];

        /* generate first sample point */
        parameter_v_offset = current_generated_state->parameter[index_v] =
          (ONE + delta_parameter_v) * parameter_v;
        parameter_vv_offset = current_generated_state->parameter[index_vv] =
          (ONE + delta_parameter_vv) * parameter_vv;
        if (parameter_v_offset > parameter_maximum[index_v] ||
            parameter_v_offset < parameter_minimum[index_v]) {
          delta_parameter_v = -delta_parameter_v;
          current_generated_state->parameter[index_v] =
            (ONE + delta_parameter_v) * parameter_v;
        }
        if (parameter_vv_offset > parameter_maximum[index_vv] ||
            parameter_vv_offset < parameter_minimum[index_vv]) {
          delta_parameter_vv = -delta_parameter_vv;
          current_generated_state->parameter[index_vv] =
            (ONE + delta_parameter_vv) * parameter_vv;
        }

        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_1 = current_generated_state->cost;

        /* restore the v_th parameter */
        current_generated_state->parameter[index_v] = parameter_v;

        /* generate second sample point */
        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_2 = current_generated_state->cost;

        /* restore the vv_th parameter */
        current_generated_state->parameter[index_vv] = parameter_vv;

        /* generate third sample point */
        current_generated_state->parameter[index_v] =
          (ONE + delta_parameter_v) * parameter_v;
        *valid_state_generated_flag = TRUE;
#if USER_ACCEPTANCE_TEST
        OPTIONS->User_Acceptance_Flag = TRUE;
        OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
        current_generated_state->cost =
          user_cost_function (current_generated_state->parameter,
                              parameter_minimum,
                              parameter_maximum,
                              tangents,
                              curvature,
                              number_parameters,
                              parameter_type,
                              valid_state_generated_flag,
                              exit_status, OPTIONS);
        if ((*valid_state_generated_flag == FALSE)
            || ((current_generated_state->cost) !=
                (current_generated_state->cost))
            || current_generated_state->cost < -MAX_DOUBLE
            || current_generated_state->cost > MAX_DOUBLE) {
          *exit_status = INVALID_COST_FUNCTION_DERIV;
          return;
        }
        if (*valid_state_generated_flag == FALSE)
          ++(*number_invalid_generated_states);
        new_cost_state_3 = current_generated_state->cost;

        /* restore the v_th parameter */
        current_generated_state->parameter[index_v] = parameter_v;

        /* calculate and store the curvature */
        curvature[index_vv_v] = curvature[index_v_vv] =
          (new_cost_state_1 - new_cost_state_2
           - new_cost_state_3 + recent_best_cost)
          / (delta_parameter_v * delta_parameter_vv
             * parameter_v * parameter_vv + (double) EPS_DOUBLE);
      }
    }
  }

  /* restore Immediate_Exit flag */
  OPTIONS->Immediate_Exit = immediate_flag;

  /* restore the best cost function value */
  current_generated_state->cost = recent_best_cost;
#if ASA_PRINT
  tmp_saved = *number_invalid_generated_states - saved_num_invalid_gen_states;
  if (tmp_saved > 0)
#if INT_LONG
    fprintf (ptr_asa_out,
             "Generated %ld invalid states when calculating the derivatives\n",
             tmp_saved);
#else
    fprintf (ptr_asa_out,
             "Generated %d invalid states when calculating the derivatives\n",
             tmp_saved);
#endif
#endif /* ASA_PRINT */
  *number_invalid_generated_states = saved_num_invalid_gen_states;
#if USER_ACCEPTANCE_TEST
  OPTIONS->User_Acceptance_Flag = TRUE;
  OPTIONS->Cost_Acceptance_Flag = FALSE;
#endif
}

/***********************************************************************
* asa_test_asa_options
*       Tests user's selected options
***********************************************************************/
#if HAVE_ANSI
int

asa_test_asa_options (LONG_INT * seed,
                      double *parameter_initial_final,
                      double *parameter_minimum,
                      double *parameter_maximum,
                      double *tangents,
                      double *curvature,
                      ALLOC_INT * number_parameters,
                      int *parameter_type,
                      int *valid_state_generated_flag,
                      int *exit_status,
                      FILE * ptr_asa_out, USER_DEFINES * OPTIONS)
#else
int

asa_test_asa_options (seed,
                      parameter_initial_final,
                      parameter_minimum,
                      parameter_maximum,
                      tangents,
                      curvature,
                      number_parameters,
                      parameter_type,
                      valid_state_generated_flag,
                      exit_status, ptr_asa_out, OPTIONS)
     LONG_INT *seed;
     double *parameter_initial_final;
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *curvature;
     ALLOC_INT *number_parameters;
     int *parameter_type;
     int *valid_state_generated_flag;
     int *exit_status;
     FILE *ptr_asa_out;
     USER_DEFINES *OPTIONS;
#endif /* HAVE_ANSI */
{
  int invalid, index_v;

  invalid = 0;

  if (seed == NULL) {
    strcpy (exit_msg, "*** seed == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (parameter_initial_final == NULL) {
    strcpy (exit_msg, "*** parameter_initial_final == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (parameter_minimum == NULL) {
    strcpy (exit_msg, "*** parameter_minimum == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (parameter_maximum == NULL) {
    strcpy (exit_msg, "*** parameter_maximum == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (tangents == NULL) {
    strcpy (exit_msg, "*** tangents == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if ((OPTIONS->Curvature_0 == FALSE) || (OPTIONS->Curvature_0 == -1)) {
    if (curvature == NULL) {
      strcpy (exit_msg, "*** curvature == NULL ***");
      print_string (ptr_asa_out, exit_msg);
      ++invalid;
    }
  }
  if (number_parameters == NULL) {
    strcpy (exit_msg, "*** number_parameters == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (parameter_type == NULL) {
    strcpy (exit_msg, "*** parameter_type == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (valid_state_generated_flag == NULL) {
    strcpy (exit_msg, "*** valid_state_generated_flag == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (exit_status == NULL) {
    strcpy (exit_msg, "*** exit_status == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS == NULL) {
    strcpy (exit_msg, "*** OPTIONS == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }

  VFOR (index_v) if (parameter_minimum[index_v] > parameter_maximum[index_v]) {
    strcpy (exit_msg, "*** parameter_minimum[] > parameter_maximum[] ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
  VFOR (index_v)
    if (parameter_initial_final[index_v] < parameter_minimum[index_v]) {
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
    strcpy (exit_msg, "*** parameter_initial[] < parameter_minimum[] ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
  VFOR (index_v)
    if (parameter_initial_final[index_v] > parameter_maximum[index_v]) {
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
    strcpy (exit_msg, "*** parameter_initial[] > parameter_maximum[] ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
  if (*number_parameters < 1) {
    strcpy (exit_msg, "*** *number_parameters < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  VFOR (index_v)
    if (parameter_type[index_v] != -2 && parameter_type[index_v] != 2
        && parameter_type[index_v] != -1 && parameter_type[index_v] != 1) {
    strcpy (exit_msg,
            "*** parameter_type[] != -2 && parameter_type[] != 2 && parameter_type[] != -1 && parameter_type[] != 1 ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }

  if (OPTIONS_FILE != FALSE && OPTIONS_FILE != TRUE) {
    strcpy (exit_msg,
            "*** OPTIONS_FILE != FALSE && OPTIONS_FILE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS_FILE_DATA != FALSE && OPTIONS_FILE_DATA != TRUE) {
    strcpy (exit_msg,
            "*** OPTIONS_FILE_DATA != FALSE && OPTIONS_FILE_DATA != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (RECUR_OPTIONS_FILE != FALSE && RECUR_OPTIONS_FILE != TRUE) {
    strcpy (exit_msg,
            "*** RECUR_OPTIONS_FILE != FALSE && RECUR_OPTIONS_FILE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (RECUR_OPTIONS_FILE_DATA != FALSE && RECUR_OPTIONS_FILE_DATA != TRUE) {
    strcpy (exit_msg,
            "*** RECUR_OPTIONS_FILE_DATA != FALSE && RECUR_OPTIONS_FILE_DATA != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (COST_FILE != FALSE && COST_FILE != TRUE) {
    strcpy (exit_msg, "*** COST_FILE != FALSE && COST_FILE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_LIB != FALSE && ASA_LIB != TRUE) {
    strcpy (exit_msg, "*** ASA_LIB != FALSE && ASA_LIB != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (MY_TEMPLATE != FALSE && MY_TEMPLATE != TRUE) {
    strcpy (exit_msg, "*** MY_TEMPLATE != FALSE && MY_TEMPLATE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_LIB != FALSE && ASA_TEMPLATE_LIB != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_LIB != FALSE && ASA_TEMPLATE_LIB != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (HAVE_ANSI != FALSE && HAVE_ANSI != TRUE) {
    strcpy (exit_msg, "*** HAVE_ANSI != FALSE && HAVE_ANSI != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (IO_PROTOTYPES != FALSE && IO_PROTOTYPES != TRUE) {
    strcpy (exit_msg,
            "*** IO_PROTOTYPES != FALSE && IO_PROTOTYPES != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (TIME_CALC != FALSE && TIME_CALC != TRUE) {
    strcpy (exit_msg, "*** TIME_CALC != FALSE && TIME_CALC != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (TIME_STD != FALSE && TIME_STD != TRUE) {
    strcpy (exit_msg, "*** TIME_STD != FALSE && TIME_STD != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (TIME_GETRUSAGE != FALSE && TIME_GETRUSAGE != TRUE) {
    strcpy (exit_msg,
            "*** TIME_GETRUSAGE != FALSE && TIME_GETRUSAGE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (INT_LONG != FALSE && INT_LONG != TRUE) {
    strcpy (exit_msg, "*** INT_LONG != FALSE && INT_LONG != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (INT_ALLOC != FALSE && INT_ALLOC != TRUE) {
    strcpy (exit_msg, "*** INT_ALLOC != FALSE && INT_ALLOC != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (SMALL_FLOAT < ZERO) {
    strcpy (exit_msg, "*** SMALL_FLOAT < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (MIN_DOUBLE < ZERO) {
    strcpy (exit_msg, "*** MIN_DOUBLE < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (MAX_DOUBLE < ZERO) {
    strcpy (exit_msg, "*** MAX_DOUBLE < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (EPS_DOUBLE < ZERO) {
    strcpy (exit_msg, "*** EPS_DOUBLE < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (CHECK_EXPONENT != FALSE && CHECK_EXPONENT != TRUE) {
    strcpy (exit_msg,
            "*** CHECK_EXPONENT != FALSE && CHECK_EXPONENT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (NO_PARAM_TEMP_TEST != FALSE && NO_PARAM_TEMP_TEST != TRUE) {
    strcpy (exit_msg,
            "*** NO_PARAM_TEMP_TEST != FALSE && NO_PARAM_TEMP_TEST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (NO_COST_TEMP_TEST != FALSE && NO_COST_TEMP_TEST != TRUE) {
    strcpy (exit_msg,
            "*** NO_COST_TEMP_TEST != FALSE && NO_COST_TEMP_TEST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (SELF_OPTIMIZE != FALSE && SELF_OPTIMIZE != TRUE) {
    strcpy (exit_msg,
            "*** SELF_OPTIMIZE != FALSE && SELF_OPTIMIZE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEST != FALSE && ASA_TEST != TRUE) {
    strcpy (exit_msg, "*** ASA_TEST != FALSE && ASA_TEST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEST_POINT != FALSE && ASA_TEST_POINT != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEST_POINT != FALSE && ASA_TEST_POINT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_EXIT_ANYTIME != FALSE && ASA_EXIT_ANYTIME != TRUE) {
    strcpy (exit_msg,
            "*** ASA_EXIT_ANYTIME != FALSE && ASA_EXIT_ANYTIME != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE != FALSE) {
    strcpy (exit_msg, "*** ASA_TEMPLATE != FALSE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_ASA_OUT_PID != FALSE && ASA_TEMPLATE_ASA_OUT_PID != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_ASA_OUT_PID != FALSE && ASA_TEMPLATE_ASA_OUT_PID != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_MULTIPLE != FALSE && ASA_TEMPLATE_MULTIPLE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_MULTIPLE != FALSE && ASA_TEMPLATE_MULTIPLE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_SELFOPT != FALSE && ASA_TEMPLATE_SELFOPT != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_SELFOPT != FALSE && ASA_TEMPLATE_SELFOPT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_SAMPLE != FALSE && ASA_TEMPLATE_SAMPLE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_SAMPLE != FALSE && ASA_TEMPLATE_SAMPLE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_QUEUE != FALSE && ASA_TEMPLATE_QUEUE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_QUEUE != FALSE && ASA_TEMPLATE_QUEUE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_PARALLEL != FALSE && ASA_TEMPLATE_PARALLEL != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_PARALLEL != FALSE && ASA_TEMPLATE_PARALLEL != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_TEMPLATE_SAVE != FALSE && ASA_TEMPLATE_SAVE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_TEMPLATE_SAVE != FALSE && ASA_TEMPLATE_SAVE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_INITIAL_COST_TEMP != FALSE && USER_INITIAL_COST_TEMP != TRUE) {
    strcpy (exit_msg,
            "*** USER_INITIAL_COST_TEMP != FALSE && USER_INITIAL_COST_TEMP != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (RATIO_TEMPERATURE_SCALES != FALSE && RATIO_TEMPERATURE_SCALES != TRUE) {
    strcpy (exit_msg,
            "*** RATIO_TEMPERATURE_SCALES != FALSE && RATIO_TEMPERATURE_SCALES != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_INITIAL_PARAMETERS_TEMPS != FALSE
      && USER_INITIAL_PARAMETERS_TEMPS != TRUE) {
    strcpy (exit_msg,
            "*** USER_INITIAL_PARAMETERS_TEMPS != FALSE && USER_INITIAL_PARAMETERS_TEMPS != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (DELTA_PARAMETERS != FALSE && DELTA_PARAMETERS != TRUE) {
    strcpy (exit_msg,
            "*** DELTA_PARAMETERS != FALSE && DELTA_PARAMETERS != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (QUENCH_PARAMETERS != FALSE && QUENCH_PARAMETERS != TRUE) {
    strcpy (exit_msg,
            "*** QUENCH_PARAMETERS != FALSE && QUENCH_PARAMETERS != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (QUENCH_COST != FALSE && QUENCH_COST != TRUE) {
    strcpy (exit_msg, "*** QUENCH_COST != FALSE && QUENCH_COST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (QUENCH_PARAMETERS_SCALE != FALSE && QUENCH_PARAMETERS_SCALE != TRUE) {
    strcpy (exit_msg,
            "*** QUENCH_PARAMETERS_SCALE != FALSE && QUENCH_PARAMETERS_SCALE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (QUENCH_COST_SCALE != FALSE && QUENCH_COST_SCALE != TRUE) {
    strcpy (exit_msg,
            "*** QUENCH_COST_SCALE != FALSE && QUENCH_COST_SCALE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONAL_DATA_DBL != FALSE && OPTIONAL_DATA_DBL != TRUE) {
    strcpy (exit_msg,
            "*** OPTIONAL_DATA_DBL != FALSE && OPTIONAL_DATA_DBL != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONAL_DATA_INT != FALSE && OPTIONAL_DATA_INT != TRUE) {
    strcpy (exit_msg,
            "*** OPTIONAL_DATA_INT != FALSE && OPTIONAL_DATA_INT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONAL_DATA_PTR != FALSE && OPTIONAL_DATA_PTR != TRUE) {
    strcpy (exit_msg,
            "*** OPTIONAL_DATA_PTR != FALSE && OPTIONAL_DATA_PTR != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_COST_SCHEDULE != FALSE && USER_COST_SCHEDULE != TRUE) {
    strcpy (exit_msg,
            "*** USER_COST_SCHEDULE != FALSE && USER_COST_SCHEDULE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_ACCEPT_ASYMP_EXP != FALSE && USER_ACCEPT_ASYMP_EXP != TRUE) {
    strcpy (exit_msg,
            "*** USER_ACCEPT_ASYMP_EXP != FALSE && USER_ACCEPT_ASYMP_EXP != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_ACCEPT_THRESHOLD != FALSE && USER_ACCEPT_THRESHOLD != TRUE) {
    strcpy (exit_msg,
            "*** USER_ACCEPT_THRESHOLD != FALSE && USER_ACCEPT_THRESHOLD != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_ACCEPTANCE_TEST != FALSE && USER_ACCEPTANCE_TEST != TRUE) {
    strcpy (exit_msg,
            "*** USER_ACCEPTANCE_TEST != FALSE && USER_ACCEPTANCE_TEST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_GENERATING_FUNCTION != FALSE && USER_GENERATING_FUNCTION != TRUE) {
    strcpy (exit_msg,
            "*** USER_GENERATING_FUNCTION != FALSE && USER_GENERATING_FUNCTION != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_REANNEAL_COST != FALSE && USER_REANNEAL_COST != TRUE) {
    strcpy (exit_msg,
            "*** USER_REANNEAL_COST != FALSE && USER_REANNEAL_COST != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_REANNEAL_PARAMETERS != FALSE && USER_REANNEAL_PARAMETERS != TRUE) {
    strcpy (exit_msg,
            "*** USER_REANNEAL_PARAMETERS != FALSE && USER_REANNEAL_PARAMETERS != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (MAXIMUM_REANNEAL_INDEX < 1) {
    strcpy (exit_msg, "*** MAXIMUM_REANNEAL_INDEX < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (REANNEAL_SCALE < ZERO) {
    strcpy (exit_msg, "*** REANNEAL_SCALE < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_SAMPLE != FALSE && ASA_SAMPLE != TRUE) {
    strcpy (exit_msg, "*** ASA_SAMPLE != FALSE && ASA_SAMPLE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ADAPTIVE_OPTIONS != FALSE && ADAPTIVE_OPTIONS != TRUE) {
    strcpy (exit_msg,
            "*** ADAPTIVE_OPTIONS != FALSE && ADAPTIVE_OPTIONS != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_QUEUE != FALSE && ASA_QUEUE != TRUE) {
    strcpy (exit_msg, "*** ASA_QUEUE != FALSE && ASA_QUEUE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_RESOLUTION != FALSE && ASA_RESOLUTION != TRUE) {
    strcpy (exit_msg,
            "*** ASA_RESOLUTION != FALSE && ASA_RESOLUTION != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_FUZZY != FALSE && ASA_FUZZY != TRUE) {
    strcpy (exit_msg, "*** ASA_FUZZY != FALSE && ASA_FUZZY != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_FUZZY_PRINT != FALSE && ASA_FUZZY_PRINT != TRUE) {
    strcpy (exit_msg,
            "*** ASA_FUZZY_PRINT != FALSE && ASA_FUZZY_PRINT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FITLOC != FALSE && FITLOC != TRUE) {
    strcpy (exit_msg, "*** FITLOC != FALSE && FITLOC != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FITLOC_ROUND != FALSE && FITLOC_ROUND != TRUE) {
    strcpy (exit_msg,
            "*** FITLOC_ROUND != FALSE && FITLOC_ROUND != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FITLOC_PRINT != FALSE && FITLOC_PRINT != TRUE) {
    strcpy (exit_msg,
            "*** FITLOC_PRINT != FALSE && FITLOC_PRINT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (MULTI_MIN != FALSE && MULTI_MIN != TRUE) {
    strcpy (exit_msg, "*** MULTI_MIN != FALSE && MULTI_MIN != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#if MULTI_MIN
  if (OPTIONS->Multi_Number <= 0) {
    strcpy (exit_msg, "*** OPTIONS->Multi_Number <= 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  VFOR (index_v) {
    if (((OPTIONS->Multi_Grid[index_v]) != (OPTIONS->Multi_Grid[index_v]))
        || OPTIONS->Multi_Grid[index_v] < 0) {
      strcpy (exit_msg,
              "*** (OPTIONS->Multi_Grid[]) != (OPTIONS->Multi_Grid[]) || OPTIONS->Multi_Grid[] < 0 ***");
      print_string_index (ptr_asa_out, exit_msg, index_v);
      ++invalid;
    }
  }
  if (OPTIONS->Multi_Specify != 0 && OPTIONS->Multi_Specify != 1) {
    strcpy (exit_msg,
            "*** OPTIONS->Multi_Specify != 0 && OPTIONS->Multi_Specify != 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
  if (ASA_PARALLEL != FALSE && ASA_PARALLEL != TRUE) {
    strcpy (exit_msg,
            "*** ASA_PARALLEL != FALSE && ASA_PARALLEL != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_SAVE != FALSE && ASA_SAVE != TRUE) {
    strcpy (exit_msg, "*** ASA_SAVE != FALSE && ASA_SAVE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_SAVE_OPT != FALSE && ASA_SAVE_OPT != TRUE) {
    strcpy (exit_msg,
            "*** ASA_SAVE_OPT != FALSE && ASA_SAVE_OPT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_SAVE_BACKUP != FALSE && ASA_SAVE_BACKUP != TRUE) {
    strcpy (exit_msg,
            "*** ASA_SAVE_BACKUP != FALSE && ASA_SAVE_BACKUP != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_PIPE != FALSE && ASA_PIPE != TRUE) {
    strcpy (exit_msg, "*** ASA_PIPE != FALSE && ASA_PIPE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_PIPE_FILE != FALSE && ASA_PIPE_FILE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_PIPE_FILE != FALSE && ASA_PIPE_FILE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (SYSTEM_CALL != FALSE && SYSTEM_CALL != TRUE) {
    strcpy (exit_msg, "*** SYSTEM_CALL != FALSE && SYSTEM_CALL != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FDLIBM_POW != FALSE && FDLIBM_POW != TRUE) {
    strcpy (exit_msg, "*** FDLIBM_POW != FALSE && FDLIBM_POW != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FDLIBM_LOG != FALSE && FDLIBM_LOG != TRUE) {
    strcpy (exit_msg, "*** FDLIBM_LOG != FALSE && FDLIBM_LOG != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (FDLIBM_EXP != FALSE && FDLIBM_EXP != TRUE) {
    strcpy (exit_msg, "*** FDLIBM_EXP != FALSE && FDLIBM_EXP != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_PRINT != FALSE && ASA_PRINT != TRUE) {
    strcpy (exit_msg, "*** ASA_PRINT != FALSE && ASA_PRINT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_ASA_OUT != FALSE && USER_ASA_OUT != TRUE) {
    strcpy (exit_msg,
            "*** USER_ASA_OUT != FALSE && USER_ASA_OUT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (USER_ASA_USR_OUT != FALSE && USER_ASA_USR_OUT != TRUE) {
    strcpy (exit_msg,
            "*** USER_ASA_USR_OUT != FALSE && USER_ASA_USR_OUT != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_PRINT_INTERMED != FALSE && ASA_PRINT_INTERMED != TRUE) {
    strcpy (exit_msg,
            "*** ASA_PRINT_INTERMED != FALSE && ASA_PRINT_INTERMED != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (ASA_PRINT_MORE != FALSE && ASA_PRINT_MORE != TRUE) {
    strcpy (exit_msg,
            "*** ASA_PRINT_MORE != FALSE && ASA_PRINT_MORE != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (G_FIELD < 0) {
    strcpy (exit_msg, "*** G_FIELD < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (G_PRECISION < 0) {
    strcpy (exit_msg, "*** G_PRECISION < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }

  if (OPTIONS->Limit_Acceptances < 0) {
    strcpy (exit_msg, "*** Limit_Acceptances < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Limit_Generated < 0) {
    strcpy (exit_msg, "*** Limit_Generated < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Limit_Invalid_Generated_States < 0) {
    strcpy (exit_msg, "*** Limit_Invalid_Generated_States < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Accepted_To_Generated_Ratio <= ZERO) {
    strcpy (exit_msg, "*** Accepted_To_Generated_Ratio <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Cost_Precision <= ZERO) {
    strcpy (exit_msg, "*** Cost_Precision <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Maximum_Cost_Repeat < 0) {
    strcpy (exit_msg, "*** Maximum_Cost_Repeat < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Number_Cost_Samples == 0 || OPTIONS->Number_Cost_Samples == -1) {
    strcpy (exit_msg,
            "*** Number_Cost_Samples == 0 || Number_Cost_Samples == -1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Temperature_Ratio_Scale <= ZERO) {
    strcpy (exit_msg, "*** Temperature_Ratio_Scale <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Cost_Parameter_Scale_Ratio <= ZERO) {
    strcpy (exit_msg, "*** Cost_Parameter_Scale_Ratio <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Temperature_Anneal_Scale <= ZERO) {
    strcpy (exit_msg, "*** Temperature_Anneal_Scale <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#if USER_INITIAL_COST_TEMP
  if (OPTIONS->User_Cost_Temperature[0] <= ZERO) {
    strcpy (exit_msg, "*** User_Cost_Temperature[0] <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
  if (OPTIONS->Include_Integer_Parameters != FALSE
      && OPTIONS->Include_Integer_Parameters != TRUE) {
    strcpy (exit_msg, "");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->User_Initial_Parameters != FALSE
      && OPTIONS->User_Initial_Parameters != TRUE) {
    strcpy (exit_msg,
            "*** User_Initial_Parameters != FALSE && User_Initial_Parameters != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Sequential_Parameters >= *number_parameters) {
    strcpy (exit_msg, "*** Sequential_Parameters >= *number_parameters ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Initial_Parameter_Temperature <= ZERO) {
    strcpy (exit_msg, "*** Initial_Parameter_Temperature <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#if RATIO_TEMPERATURE_SCALES
  VFOR (index_v) if (OPTIONS->User_Temperature_Ratio[index_v] <= ZERO) {
    strcpy (exit_msg, "*** User_Temperature_Ratio[] <= ZERO ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
#endif
#if USER_INITIAL_PARAMETERS_TEMPS
  VFOR (index_v) if (OPTIONS->User_Parameter_Temperature[index_v] <= ZERO) {
    strcpy (exit_msg, "*** User_Parameter_Temperature[] <= ZERO ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
#endif
  if (OPTIONS->Acceptance_Frequency_Modulus < 0) {
    strcpy (exit_msg, "*** Acceptance_Frequency_Modulus < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Generated_Frequency_Modulus < 0) {
    strcpy (exit_msg, "*** Generated_Frequency_Modulus < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Reanneal_Cost == -1) {
    strcpy (exit_msg, "*** Reanneal_Cost == -1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Reanneal_Parameters != FALSE
      && OPTIONS->Reanneal_Parameters != TRUE) {
    strcpy (exit_msg,
            "*** Reanneal_Parameters != FALSE && Reanneal_Parameters != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Delta_X < ZERO) {
    strcpy (exit_msg, "*** Delta_X < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#if DELTA_PARAMETERS
  VFOR (index_v) if (OPTIONS->User_Delta_Parameter[index_v] < ZERO) {
    strcpy (exit_msg, "*** User_Delta_Parameter[] < ZERO ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
#endif
  if (OPTIONS->User_Tangents != FALSE && OPTIONS->User_Tangents != TRUE) {
    strcpy (exit_msg,
            "*** User_Tangents != FALSE && User_Tangents != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Curvature_0 != -1 && OPTIONS->Curvature_0 != FALSE
      && OPTIONS->Curvature_0 != TRUE) {
    strcpy (exit_msg,
            "*** Curvature_0 -1 && Curvature_0 != FALSE && Curvature_0 != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#if QUENCH_PARAMETERS
  VFOR (index_v) if (OPTIONS->User_Quench_Param_Scale[index_v] <= ZERO) {
    strcpy (exit_msg, "*** User_Quench_Param_Scale[] <= ZERO ***");
    print_string_index (ptr_asa_out, exit_msg, index_v);
    ++invalid;
  }
#endif
#if QUENCH_COST
  if (OPTIONS->User_Quench_Cost_Scale[0] <= ZERO) {
    strcpy (exit_msg, "*** User_Quench_Cost_Scale[0] <= ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if OPTIONAL_DATA_DBL
  if (OPTIONS->Asa_Data_Dim_Dbl < 1) {
    strcpy (exit_msg, "*** Asa_Data_Dim_Dbl < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Asa_Data_Dbl == NULL) {
    strcpy (exit_msg, "*** Asa_Data_Dbl == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if ASA_SAVE
  if (OPTIONS->Random_Array_Dim < 1) {
    strcpy (exit_msg, "*** Random_Array_Dim < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Random_Array == NULL) {
    strcpy (exit_msg, "*** Random_Array == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if OPTIONAL_DATA_INT
  if (OPTIONS->Asa_Data_Dim_Int < 1) {
    strcpy (exit_msg, "*** Asa_Data_Dim_Int < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Asa_Data_Int == NULL) {
    strcpy (exit_msg, "*** Asa_Data_Int == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if OPTIONAL_DATA_PTR
  if (OPTIONS->Asa_Data_Dim_Ptr < 1) {
    strcpy (exit_msg, "*** Asa_Data_Dim_Ptr < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Asa_Data_Ptr == NULL) {
    strcpy (exit_msg, "*** Asa_Data_Ptr == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_ASA_OUT
  if (OPTIONS->Asa_Out_File == NULL) {
    strcpy (exit_msg, "*** Asa_Out_File == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_COST_SCHEDULE
  if (OPTIONS->Cost_Schedule == NULL) {
    strcpy (exit_msg, "*** Cost_Schedule == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_ACCEPTANCE_TEST
  if (OPTIONS->Acceptance_Test == NULL) {
    strcpy (exit_msg, "*** Acceptance_Test == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->User_Acceptance_Flag != FALSE
      && OPTIONS->User_Acceptance_Flag != TRUE) {
    strcpy (exit_msg,
            "*** User_Acceptance_Flag != FALSE && User_Acceptance_Flag != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Cost_Acceptance_Flag != FALSE
      && OPTIONS->Cost_Acceptance_Flag != TRUE) {
    strcpy (exit_msg,
            "*** Cost_Acceptance_Flag != FALSE && Cost_Acceptance_Flag != TRUE ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_GENERATING_FUNCTION
  if (OPTIONS->Generating_Distrib == NULL) {
    strcpy (exit_msg, "*** Generating_Distrib == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_REANNEAL_COST
  if (OPTIONS->Reanneal_Cost_Function == NULL) {
    strcpy (exit_msg, "*** Reanneal_Cost_Function == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if USER_REANNEAL_PARAMETERS
  if (OPTIONS->Reanneal_Params_Function == NULL) {
    strcpy (exit_msg, "*** Reanneal_Params_Function == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if ASA_SAMPLE
  if (OPTIONS->Bias_Generated == NULL) {
    strcpy (exit_msg, "*** Bias_Generated == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Limit_Weights < ZERO) {
    strcpy (exit_msg, "*** Limit_Weights < ZERO ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if ASA_QUEUE
  if (OPTIONS->Queue_Size < 0) {
    strcpy (exit_msg, "*** Queue_Size < 0 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Queue_Size > 0) {
    if (OPTIONS->Queue_Resolution == NULL) {
      strcpy (exit_msg, "*** Queue_Resolution == NULL ***");
      print_string (ptr_asa_out, exit_msg);
      ++invalid;
    }
  }
#endif /* ASA_QUEUE */
#if ASA_RESOLUTION
  if (OPTIONS->Coarse_Resolution == NULL) {
    strcpy (exit_msg, "*** Coarse_Resolution == NULL ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif
#if ASA_PARALLEL
  if (OPTIONS->Gener_Block < 1) {
    strcpy (exit_msg, "*** Gener_Block < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Gener_Block_Max < 1) {
    strcpy (exit_msg, "*** Gener_Block_Max < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
  if (OPTIONS->Gener_Mov_Avr < 1) {
    strcpy (exit_msg, "*** Gener_Mov_Avr < 1 ***");
    print_string (ptr_asa_out, exit_msg);
    ++invalid;
  }
#endif /* ASA_PARALLEL */

  return (invalid);
}

/***********************************************************************
* cost_function_test
*       Tests user's returned cost function values and parameters
***********************************************************************/
#if HAVE_ANSI
int

cost_function_test (double cost,
                    double *parameter,
                    double *parameter_minimum,
                    double *parameter_maximum,
                    ALLOC_INT * number_parameters, double *xnumber_parameters)
#else
int

cost_function_test (cost,
                    parameter,
                    parameter_minimum, parameter_maximum,
                    number_parameters, xnumber_parameters)
     double cost;
     double *parameter;
     double *parameter_minimum;
     double *parameter_maximum;
     ALLOC_INT *number_parameters;
     double *xnumber_parameters;
#endif /* HAVE_ANSI */
{
  ALLOC_INT index_v;
  int test_flag;

  test_flag = 1;

  if (((cost) != (cost)) || (cost < -MAX_DOUBLE || cost > MAX_DOUBLE))
    test_flag = 0;

  *xnumber_parameters = (double) *number_parameters;
  VFOR (index_v) {
    if (PARAMETER_RANGE_TOO_SMALL (index_v)) {
      *xnumber_parameters -= 1.0;
      continue;
    }
    if (parameter[index_v] < parameter_minimum[index_v] ||
        parameter[index_v] > parameter_maximum[index_v]) {
      test_flag = 0;
    }
  }

  return (test_flag);
}

/***********************************************************************
* print_string
*	This prints the designated string
***********************************************************************/
#if HAVE_ANSI
void
print_string (FILE * ptr_asa_out, char *string)
#else
void
print_string (ptr_asa_out, string)
     FILE *ptr_asa_out;
     char *string;
#endif /* HAVE_ANSI */
{
#if INCL_STDOUT
  printf ("\n\n%s\n\n", string);
  fflush (stdout);
#endif /* INCL_STDOUT */
#if ASA_PRINT
  fprintf (ptr_asa_out, "\n\n%s\n\n", string);
  fflush (ptr_asa_out);
#else
#endif
}

/***********************************************************************
* print_string_index
*	This prints the designated string and index
***********************************************************************/
#if HAVE_ANSI
void
print_string_index (FILE * ptr_asa_out, char *string, ALLOC_INT index)
#else
void
print_string_index (ptr_asa_out, string, index)
     FILE *ptr_asa_out;
     char *string;
     ALLOC_INT index;
#endif /* HAVE_ANSI */
{
#if INCL_STDOUT
#if INT_ALLOC
  printf ("\n\n%s index = %d\n\n", string, index);
#else /* INT_ALLOC */
#if INT_LONG
  printf ("\n\n%s index = %ld\n\n", string, index);
#else /* INT_LONG */
  printf ("\n\n%s index = %ld\n\n", string, index);
#endif /* INT_LONG */
#endif /* INT_ALLOC */
  fflush (stdout);
#endif /* INCL_STDOUT */

#if ASA_PRINT
#if INT_ALLOC
  fprintf (ptr_asa_out, "\n\n%s index = %d\n\n", string, index);
#else /* INT_ALLOC */
#if INT_LONG
  fprintf (ptr_asa_out, "\n\n%s index = %ld\n\n", string, index);
#else /* INT_LONG */
  fprintf (ptr_asa_out, "\n\n%s index = %d\n\n", string, index);
#endif /* INT_LONG */
#endif /* INT_ALLOC */
  fflush (ptr_asa_out);
#else /* ASA_PRINT */
  ;
#endif /* ASA_PRINT */
}

#if ASA_PRINT
/***********************************************************************
* print_state
*	Prints a description of the current state of the system
***********************************************************************/
#if HAVE_ANSI
void

print_state (double *parameter_minimum,
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
             FILE * ptr_asa_out, USER_DEFINES * OPTIONS)
#else
void

print_state (parameter_minimum,
             parameter_maximum,
             tangents,
             curvature,
             current_cost_temperature,
             current_user_parameter_temp,
             accepted_to_generated_ratio,
             number_parameters,
             curvature_flag,
             number_accepted,
             index_cost_acceptances,
             number_generated,
             number_invalid_generated_states,
             last_saved_state, best_generated_state, ptr_asa_out, OPTIONS)
     double *parameter_minimum;
     double *parameter_maximum;
     double *tangents;
     double *curvature;
     double *current_cost_temperature;
     double *current_user_parameter_temp;
     double *accepted_to_generated_ratio;
     ALLOC_INT *number_parameters;
     int *curvature_flag;
     LONG_INT *number_accepted;
     LONG_INT *index_cost_acceptances;
     LONG_INT *number_generated;
     LONG_INT *number_invalid_generated_states;
     STATE *last_saved_state;
     STATE *best_generated_state;
     FILE *ptr_asa_out;
     USER_DEFINES *OPTIONS;
#endif /* HAVE_ANSI */
{
  ALLOC_INT index_v;
  ALLOC_INT index_vv, index_v_vv;

  fprintf (ptr_asa_out, "\n");
#if TIME_CALC
  print_time ("", ptr_asa_out);
#endif

  if (OPTIONS->Curvature_0 == TRUE)
    *curvature_flag = FALSE;
  if (OPTIONS->Curvature_0 == -1)
    *curvature_flag = TRUE;

#if INT_LONG
  fprintf (ptr_asa_out,
           "*index_cost_acceptances = %ld, *current_cost_temperature = %*.*g\n",
           *index_cost_acceptances,
           G_FIELD, G_PRECISION, *current_cost_temperature);
  fprintf (ptr_asa_out,
           "*accepted_to_generated_ratio = %*.*g, *number_invalid... = %ld\n",
           G_FIELD, G_PRECISION, *accepted_to_generated_ratio,
           (*number_invalid_generated_states));
  fprintf (ptr_asa_out, "*number_generated = %ld, *number_accepted = %ld\n",
           *number_generated, *number_accepted);
#else
  fprintf (ptr_asa_out,
           "*index_cost_acceptances = %d, *current_cost_temperature = %*.*g\n",
           *index_cost_acceptances,
           G_FIELD, G_PRECISION, *current_cost_temperature);
  fprintf (ptr_asa_out,
           "*accepted_to_generated_ratio = %*.*g, *number_invalid... = %d\n",
           G_FIELD, G_PRECISION, *accepted_to_generated_ratio,
           *number_invalid_generated_states);
  fprintf (ptr_asa_out, "*number_generated = %d, *number_accepted = %d\n",
           *number_generated, *number_accepted);
#endif

  fprintf (ptr_asa_out, "best...->cost = %*.*g, last...->cost = %*.*g\n",
           G_FIELD, G_PRECISION, best_generated_state->cost, G_FIELD,
           G_PRECISION, last_saved_state->cost);

  /* Note that tangents will not be calculated until reanneal
     is called, and therefore their listing in the printout only
     is relevant then */

  fprintf (ptr_asa_out,
           "index_v  best...->parameter current_parameter_temp\ttangent\n");
  VFOR (index_v) {
    /* ignore too small ranges */
#if DROPPED_PARAMETERS
    if (PARAMETER_RANGE_TOO_SMALL (index_v))
      continue;
#endif
#if INT_ALLOC
  char *dropped_format = "%d\t%*.*g\t\t%*.*g\t%*.*g\n";
#else
#if INT_LONG
  char *dropped_format = "%ld\t%*.*g\t\t%*.*g\t%*.*g\n";
#else
  char *dropped_format = "%d\t%*.*g\t\t%*.*g\t%*.*g\n";
#endif
#endif
    fprintf (ptr_asa_out,
             dropped_format,
             index_v,
             G_FIELD, G_PRECISION, best_generated_state->parameter[index_v],
             G_FIELD, G_PRECISION, current_user_parameter_temp[index_v],
             G_FIELD, G_PRECISION, tangents[index_v]);
  }

  if (*curvature_flag == TRUE) {
    /* print curvatures */
    VFOR (index_v) {
      /* ignore too small ranges */
      if (PARAMETER_RANGE_TOO_SMALL (index_v))
        continue;
      fprintf (ptr_asa_out, "\n");
      VFOR (index_vv) {
        /* only print upper diagonal of matrix */
        if (index_v < index_vv)
          continue;
        /* ignore too small ranges (index_vv) */
        if (PARAMETER_RANGE_TOO_SMALL (index_vv))
          continue;

        /* index_v_vv: row index_v, column index_vv */
        index_v_vv = ROW_COL_INDEX (index_v, index_vv);

        if (index_v == index_vv) {
#if INT_ALLOC
  char *curvature1_format = "curvature[%d][%d] = %*.*g\n";
#else
#if INT_LONG
  char *curvature1_format = "curvature[%ld][%ld] = %*.*g\n";
#else
  char *curvature1_format = "curvature[%d][%d] = %*.*g\n";
#endif
#endif
          fprintf (ptr_asa_out,
                   curvature1_format,
                   index_v, index_vv,
                   G_FIELD, G_PRECISION, curvature[index_v_vv]);
        } else {
#if INT_ALLOC
  char *curvature2_format = "curvature[%d][%d] = %*.*g \t = curvature[%d][%d]\n";
#else
#if INT_LONG
  char *curvature2_format = "curvature[%ld][%ld] = %*.*g \t = curvature[%ld][%ld]\n";
#else
  char *curvature2_format = "curvature[%d][%d] = %*.*g \t = curvature[%d][%d]\n";
#endif
#endif
          fprintf (ptr_asa_out,
                   curvature2_format,
                   index_v, index_vv,
                   G_FIELD, G_PRECISION, curvature[index_v_vv],
                   index_vv, index_v);
        }
      }
    }
  }
  fprintf (ptr_asa_out, "\n");
  fflush (ptr_asa_out);

}

/***********************************************************************
* print_asa_options
*	Prints user's selected options
***********************************************************************/
#if HAVE_ANSI
void
print_asa_options (FILE * ptr_asa_out, USER_DEFINES * OPTIONS)
#else
void
print_asa_options (ptr_asa_out, OPTIONS)
     FILE *ptr_asa_out;
     USER_DEFINES *OPTIONS;
#endif /* HAVE_ANSI */
{
  fprintf (ptr_asa_out, "\t\tADAPTIVE SIMULATED ANNEALING\n\n");

  fprintf (ptr_asa_out, "%s\n\n", ASA_ID);

  fprintf (ptr_asa_out, "OPTIONS_FILE = %d\n", (int) OPTIONS_FILE);
  fprintf (ptr_asa_out, "OPTIONS_FILE_DATA = %d\n", (int) OPTIONS_FILE_DATA);
  fprintf (ptr_asa_out, "RECUR_OPTIONS_FILE = %d\n",
           (int) RECUR_OPTIONS_FILE);
  fprintf (ptr_asa_out, "RECUR_OPTIONS_FILE_DATA = %d\n",
           (int) RECUR_OPTIONS_FILE_DATA);
  fprintf (ptr_asa_out, "COST_FILE = %d\n", (int) COST_FILE);
  fprintf (ptr_asa_out, "ASA_LIB = %d\n", (int) ASA_LIB);
  fprintf (ptr_asa_out, "HAVE_ANSI = %d\n", (int) HAVE_ANSI);
  fprintf (ptr_asa_out, "IO_PROTOTYPES = %d\n", (int) IO_PROTOTYPES);
  fprintf (ptr_asa_out, "TIME_CALC = %d\n", (int) TIME_CALC);
  fprintf (ptr_asa_out, "TIME_STD = %d\n", (int) TIME_STD);
  fprintf (ptr_asa_out, "TIME_GETRUSAGE = %d\n", (int) TIME_GETRUSAGE);
  fprintf (ptr_asa_out, "INT_LONG = %d\n", (int) INT_LONG);
  fprintf (ptr_asa_out, "INT_ALLOC = %d\n", (int) INT_ALLOC);
  fprintf (ptr_asa_out, "SMALL_FLOAT = %*.*g\n",
           G_FIELD, G_PRECISION, (double) SMALL_FLOAT);
  fprintf (ptr_asa_out, "MIN_DOUBLE = %*.*g\n",
           G_FIELD, G_PRECISION, (double) MIN_DOUBLE);
  fprintf (ptr_asa_out, "MAX_DOUBLE = %*.*g\n",
           G_FIELD, G_PRECISION, (double) MAX_DOUBLE);
  fprintf (ptr_asa_out, "EPS_DOUBLE = %*.*g\n",
           G_FIELD, G_PRECISION, (double) EPS_DOUBLE);
  fprintf (ptr_asa_out, "CHECK_EXPONENT = %d\n", (int) CHECK_EXPONENT);
  fprintf (ptr_asa_out, "NO_PARAM_TEMP_TEST = %d\n",
           (int) NO_PARAM_TEMP_TEST);
  fprintf (ptr_asa_out, "NO_COST_TEMP_TEST = %d\n", (int) NO_COST_TEMP_TEST);
  fprintf (ptr_asa_out, "SELF_OPTIMIZE = %d\n", (int) SELF_OPTIMIZE);
  fprintf (ptr_asa_out, "ASA_TEST = %d\n", (int) ASA_TEST);
  fprintf (ptr_asa_out, "ASA_TEST_POINT = %d\n", (int) ASA_TEST_POINT);
  fprintf (ptr_asa_out, "ASA_EXIT_ANYTIME = %d\n", (int) ASA_EXIT_ANYTIME);
  fprintf (ptr_asa_out, "ASA_TEMPLATE = %d\n", (int) ASA_TEMPLATE);
  fprintf (ptr_asa_out, "MY_TEMPLATE = %d\n", (int) MY_TEMPLATE);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_LIB = %d\n", (int) ASA_TEMPLATE_LIB);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_ASA_OUT_PID = %d\n",
           (int) ASA_TEMPLATE_ASA_OUT_PID);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_MULTIPLE = %d\n",
           (int) ASA_TEMPLATE_MULTIPLE);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_SELFOPT = %d\n",
           (int) ASA_TEMPLATE_SELFOPT);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_SAMPLE = %d\n",
           (int) ASA_TEMPLATE_SAMPLE);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_QUEUE = %d\n",
           (int) ASA_TEMPLATE_QUEUE);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_PARALLEL = %d\n",
           (int) ASA_TEMPLATE_PARALLEL);
  fprintf (ptr_asa_out, "ASA_TEMPLATE_SAVE = %d\n", (int) ASA_TEMPLATE_SAVE);
  fprintf (ptr_asa_out, "USER_INITIAL_COST_TEMP = %d\n",
           (int) USER_INITIAL_COST_TEMP);
  fprintf (ptr_asa_out, "RATIO_TEMPERATURE_SCALES = %d\n",
           (int) RATIO_TEMPERATURE_SCALES);
  fprintf (ptr_asa_out, "USER_INITIAL_PARAMETERS_TEMPS = %d\n",
           (int) USER_INITIAL_PARAMETERS_TEMPS);
  fprintf (ptr_asa_out, "DELTA_PARAMETERS = %d\n", (int) DELTA_PARAMETERS);
  fprintf (ptr_asa_out, "QUENCH_PARAMETERS = %d\n", (int) QUENCH_PARAMETERS);
  fprintf (ptr_asa_out, "QUENCH_COST = %d\n", (int) QUENCH_COST);
  fprintf (ptr_asa_out, "QUENCH_PARAMETERS_SCALE = %d\n",
           (int) QUENCH_PARAMETERS_SCALE);
  fprintf (ptr_asa_out, "QUENCH_COST_SCALE = %d\n", (int) QUENCH_COST_SCALE);
  fprintf (ptr_asa_out, "OPTIONAL_DATA_DBL = %d\n", (int) OPTIONAL_DATA_DBL);
  fprintf (ptr_asa_out, "OPTIONAL_DATA_INT = %d\n", (int) OPTIONAL_DATA_INT);
  fprintf (ptr_asa_out, "OPTIONAL_DATA_PTR = %d\n", (int) OPTIONAL_DATA_PTR);
  fprintf (ptr_asa_out, "USER_COST_SCHEDULE = %d\n",
           (int) USER_COST_SCHEDULE);
  fprintf (ptr_asa_out, "USER_ACCEPT_ASYMP_EXP = %d\n",
           (int) USER_ACCEPT_ASYMP_EXP);
  fprintf (ptr_asa_out, "USER_ACCEPT_THRESHOLD = %d\n",
           (int) USER_ACCEPT_THRESHOLD);
  fprintf (ptr_asa_out, "USER_ACCEPTANCE_TEST = %d\n",
           (int) USER_ACCEPTANCE_TEST);
  fprintf (ptr_asa_out, "USER_GENERATING_FUNCTION = %d\n",
           (int) USER_GENERATING_FUNCTION);
  fprintf (ptr_asa_out, "USER_REANNEAL_COST = %d\n",
           (int) USER_REANNEAL_COST);
  fprintf (ptr_asa_out, "USER_REANNEAL_PARAMETERS = %d\n",
           (int) USER_REANNEAL_PARAMETERS);
#if INT_LONG
  fprintf (ptr_asa_out, "MAXIMUM_REANNEAL_INDEX = %ld\n",
           (LONG_INT) MAXIMUM_REANNEAL_INDEX);
#else
  fprintf (ptr_asa_out, "MAXIMUM_REANNEAL_INDEX = %d\n",
           (LONG_INT) MAXIMUM_REANNEAL_INDEX);
#endif
  fprintf (ptr_asa_out, "REANNEAL_SCALE = %*.*g\n",
           G_FIELD, G_PRECISION, (double) REANNEAL_SCALE);
  fprintf (ptr_asa_out, "ASA_SAMPLE = %d\n", (int) ASA_SAMPLE);
  fprintf (ptr_asa_out, "ADAPTIVE_OPTIONS = %d\n", (int) ADAPTIVE_OPTIONS);
  fprintf (ptr_asa_out, "ASA_QUEUE = %d\n", (int) ASA_QUEUE);
  fprintf (ptr_asa_out, "ASA_RESOLUTION = %d\n", (int) ASA_RESOLUTION);
  fprintf (ptr_asa_out, "ASA_FUZZY = %d\n", (int) ASA_FUZZY);
  fprintf (ptr_asa_out, "ASA_FUZZY_PRINT = %d\n", (int) ASA_FUZZY_PRINT);
  fprintf (ptr_asa_out, "FITLOC = %d\n", (int) FITLOC);
  fprintf (ptr_asa_out, "FITLOC_ROUND = %d\n", (int) FITLOC_ROUND);
  fprintf (ptr_asa_out, "FITLOC_PRINT = %d\n", (int) FITLOC_PRINT);
  fprintf (ptr_asa_out, "MULTI_MIN = %d\n", (int) MULTI_MIN);
  fprintf (ptr_asa_out, "ASA_PARALLEL = %d\n", (int) ASA_PARALLEL);
  fprintf (ptr_asa_out, "FDLIBM_POW = %d\n", (int) FDLIBM_POW);
  fprintf (ptr_asa_out, "FDLIBM_LOG = %d\n", (int) FDLIBM_LOG);
  fprintf (ptr_asa_out, "FDLIBM_EXP = %d\n\n", (int) FDLIBM_EXP);

  fprintf (ptr_asa_out, "ASA_PRINT = %d\n", (int) ASA_PRINT);
  fprintf (ptr_asa_out, "USER_OUT = %s\n", USER_OUT);
#if USER_ASA_OUT
  fprintf (ptr_asa_out, "ASA_OUT = %s\n", OPTIONS->Asa_Out_File);
#else
  fprintf (ptr_asa_out, "ASA_OUT = %s\n", ASA_OUT);
#endif
  fprintf (ptr_asa_out, "USER_ASA_OUT = %d\n", (int) USER_ASA_OUT);
  fprintf (ptr_asa_out, "USER_ASA_USR_OUT = %d\n", (int) USER_ASA_USR_OUT);
  fprintf (ptr_asa_out, "ASA_PRINT_INTERMED = %d\n",
           (int) ASA_PRINT_INTERMED);
  fprintf (ptr_asa_out, "ASA_PRINT_MORE = %d\n", (int) ASA_PRINT_MORE);
  fprintf (ptr_asa_out, "INCL_STDOUT = %d\n", (int) INCL_STDOUT);
  fprintf (ptr_asa_out, "G_FIELD = %d\n", (int) G_FIELD);
  fprintf (ptr_asa_out, "G_PRECISION = %d\n", (int) G_PRECISION);
  fprintf (ptr_asa_out, "ASA_SAVE = %d\n", (int) ASA_SAVE);
  fprintf (ptr_asa_out, "ASA_SAVE_OPT = %d\n", (int) ASA_SAVE_OPT);
  fprintf (ptr_asa_out, "ASA_SAVE_BACKUP = %d\n", (int) ASA_SAVE_BACKUP);
  fprintf (ptr_asa_out, "ASA_PIPE = %d\n", (int) ASA_PIPE);
  fprintf (ptr_asa_out, "ASA_PIPE_FILE = %d\n", (int) ASA_PIPE_FILE);
  fprintf (ptr_asa_out, "SYSTEM_CALL = %d\n\n", (int) SYSTEM_CALL);

#if INT_LONG
  fprintf (ptr_asa_out, "OPTIONS->Limit_Acceptances = %ld\n",
           (LONG_INT) OPTIONS->Limit_Acceptances);
  fprintf (ptr_asa_out, "OPTIONS->Limit_Generated = %ld\n",
           (LONG_INT) OPTIONS->Limit_Generated);
#else
  fprintf (ptr_asa_out, "OPTIONS->Limit_Acceptances = %d\n",
           (LONG_INT) OPTIONS->Limit_Acceptances);
  fprintf (ptr_asa_out, "OPTIONS->Limit_Generated = %d\n",
           (LONG_INT) OPTIONS->Limit_Generated);
#endif
  fprintf (ptr_asa_out, "OPTIONS->Limit_Invalid_Generated_States = %d\n",
           OPTIONS->Limit_Invalid_Generated_States);
  fprintf (ptr_asa_out, "OPTIONS->Accepted_To_Generated_Ratio = %*.*g\n\n",
           G_FIELD, G_PRECISION, OPTIONS->Accepted_To_Generated_Ratio);

  fprintf (ptr_asa_out, "OPTIONS->Cost_Precision = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Cost_Precision);
  fprintf (ptr_asa_out, "OPTIONS->Maximum_Cost_Repeat = %d\n",
           OPTIONS->Maximum_Cost_Repeat);
  fprintf (ptr_asa_out, "OPTIONS->Number_Cost_Samples = %d\n",
           OPTIONS->Number_Cost_Samples);
  fprintf (ptr_asa_out, "OPTIONS->Temperature_Ratio_Scale = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Temperature_Ratio_Scale);
  fprintf (ptr_asa_out, "OPTIONS->Cost_Parameter_Scale_Ratio = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Cost_Parameter_Scale_Ratio);
  fprintf (ptr_asa_out, "OPTIONS->Temperature_Anneal_Scale = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Temperature_Anneal_Scale);

  fprintf (ptr_asa_out, "OPTIONS->Include_Integer_Parameters = %d\n",
           OPTIONS->Include_Integer_Parameters);
  fprintf (ptr_asa_out, "OPTIONS->User_Initial_Parameters = %d\n",
           OPTIONS->User_Initial_Parameters);
#if INT_ALLOC
  fprintf (ptr_asa_out, "OPTIONS->Sequential_Parameters = %d\n",
           (int) OPTIONS->Sequential_Parameters);
#else
#if INT_LONG
  fprintf (ptr_asa_out, "OPTIONS->Sequential_Parameters = %ld\n",
           (LONG_INT) OPTIONS->Sequential_Parameters);
#else
  fprintf (ptr_asa_out, "OPTIONS->Sequential_Parameters = %d\n",
           (LONG_INT) OPTIONS->Sequential_Parameters);
#endif
#endif
  fprintf (ptr_asa_out, "OPTIONS->Initial_Parameter_Temperature = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Initial_Parameter_Temperature);

  fprintf (ptr_asa_out, "OPTIONS->Acceptance_Frequency_Modulus = %d\n",
           OPTIONS->Acceptance_Frequency_Modulus);
  fprintf (ptr_asa_out, "OPTIONS->Generated_Frequency_Modulus = %d\n",
           OPTIONS->Generated_Frequency_Modulus);
  fprintf (ptr_asa_out, "OPTIONS->Reanneal_Cost = %d\n",
           OPTIONS->Reanneal_Cost);
  fprintf (ptr_asa_out, "OPTIONS->Reanneal_Parameters = %d\n\n",
           OPTIONS->Reanneal_Parameters);

  fprintf (ptr_asa_out, "OPTIONS->Delta_X = %*.*g\n",
           G_FIELD, G_PRECISION, OPTIONS->Delta_X);
  fprintf (ptr_asa_out, "OPTIONS->User_Tangents = %d\n",
           OPTIONS->User_Tangents);
  fprintf (ptr_asa_out, "OPTIONS->Curvature_0 = %d\n", OPTIONS->Curvature_0);
  fprintf (ptr_asa_out, "OPTIONS->Asa_Recursive_Level = %d\n\n",
           OPTIONS->Asa_Recursive_Level);

  fprintf (ptr_asa_out, "\n");
}
#endif /* ASA_PRINT */

#if TIME_CALC
#if TIME_GETRUSAGE
/***********************************************************************
* print_time
*	This calculates the time and runtime and prints it.
***********************************************************************/
#if HAVE_ANSI
void
print_time (char *message, FILE * ptr_asa_out)
#else
void
print_time (message, ptr_asa_out)
     char *message;
     FILE *ptr_asa_out;
#endif /* HAVE_ANSI */
{
  int who = RUSAGE_SELF;        /* Check our own time */
  struct rusage usage;

  /* get the resource usage information */
#if TIME_STD
  syscall (SYS_GETRUSAGE, who, &usage);
#else
  getrusage (who, &usage);
#endif

  /* print the usage time in reasonable form */
  aux_print_time (&usage.ru_utime, message, ptr_asa_out);
}

/***********************************************************************
* aux_print_time
*      auxiliary print the time routine
***********************************************************************/
#if HAVE_ANSI
void
aux_print_time (struct timeval *time, char *message, FILE * ptr_asa_out)
#else
void
aux_print_time (time, message, ptr_asa_out)
     struct timeval *time;
     char *message;
     FILE *ptr_asa_out;
#endif /* HAVE_ANSI */
{
  static double sx;
  double us, s, m, h;
  double ds, dm, dh;

  /* calculate the new microseconds, seconds, minutes, hours
     and the differences since the last call */
  us = (double) ((int) ((double) EPS_DOUBLE + time->tv_usec)) / 1.E6;
  s = (double) ((int) ((double) EPS_DOUBLE + time->tv_sec)) + us;
  ds = s - sx;
  sx = s;

  h = (int) ((double) EPS_DOUBLE + s / 3600.);
  m = (int) ((double) EPS_DOUBLE + s / 60.) - 60. * h;
  s -= (3600. * h + 60. * m);
  dh = (int) ((double) EPS_DOUBLE + ds / 3600.);
  dm = (int) ((double) EPS_DOUBLE + ds / 60.) - 60. * dh;
  ds -= (3600. * dh + 60. * dm);

  /* print the statistics */
  fprintf (ptr_asa_out,
           "%s:time: %gh %gm %gs; incr: %gh %gm %gs\n",
           message, h, m, s, dh, dm, ds);
}
#else /* TIME_GETRUSAGE */
  /* Note that on many machines the time resolution of this algorithm
   * may be less than the other alternatives, e.g., rounding off the
   * number of ticks to the nearest tens of thousands.  Also, because
   * time here is typically indexed by a long integer, there typically
   * is a cycle of time in periods of fractions of an hour.  For
   * example, under Solaris 2.5.1:  The value returned by clock() is
   * defined in microseconds, since the first call to clock(), for
   * compatibility with  systems that have CPU clocks with much higher
   * resolution.  Because of this, the value returned will wrap around
   * after accumulating only 2147 seconds of CPU time (about 36 minutes).
   *
   * Set TIME_GETRUSAGE to FALSE and TIME_STD to TRUE under
   * Cygwin with -mno-cygwin
   *
   * See asa.h for two places where some additional modifications should
   * be made under SunOS 4.1.x. */

#if HAVE_ANSI
void
print_time (char *message, FILE * ptr_asa_out)
#else
void
print_time (message, ptr_asa_out)
     char *message;
     FILE *ptr_asa_out;
#endif /* HAVE_ANSI */
{
  aux_print_time (clock (), message, ptr_asa_out);

}                               /*print_time */

/***********************************************************************
* aux_print_time
*      auxiliary print the time routine
***********************************************************************/
#if HAVE_ANSI
void
aux_print_time (clock_t time, char *message, FILE * ptr_asa_out)
#else
void
aux_print_time (time, message, ptr_asa_out)
     clock_t time;
     char *message;
     FILE *ptr_asa_out;
#endif /* HAVE_ANSI */
{
  static clock_t previousTime = -1;
  clock_t diffTime;
  double clocksPerSecF = CLOCKS_PER_SEC;
  double timeF, diffF;
  double s, m, h;
  double ds, dm, dh;

  if (previousTime != -1) {
    diffTime = time - previousTime;
    timeF = time;
    diffF = diffTime;
    previousTime = time;

    s = timeF / clocksPerSecF;
    ds = diffF / clocksPerSecF;

    h = (int) ((double) EPS_DOUBLE + s / 3600.);
    m = (int) ((double) EPS_DOUBLE + s / 60.) - 60. * h;
    s -= (3600. * h + 60. * m);
    dh = (int) ((double) EPS_DOUBLE + ds / 3600.);
    dm = (int) ((double) EPS_DOUBLE + ds / 60.) - 60. * dh;
    ds -= (3600. * dh + 60. * dm);

    fprintf (ptr_asa_out,
             "%s:time: %gh %gm %gs; incr: %gh %gm %gs\n",
             message, h, m, s, dh, dm, ds);
  } else {
    /* The first call will be invalid - don't output anything. */
    fprintf (ptr_asa_out, "TIMING PARAMETERS: ticks/sec: %lu\n",
             CLOCKS_PER_SEC);
    previousTime = time;
    /* Output initial message. */
    print_time (message, ptr_asa_out);
  }
}                               /* aux_print_time */

#endif /* TIME_GETRUSAGE */

#endif /* TIME_CALC */

#if MULTI_MIN
#if HAVE_ANSI
static int
multi_compare (const void *ii, const void *jj)
#else /* HAVE_ANSI */
static int
multi_compare (ii, jj)
     char *ii;
     char *jj;
#endif /* HAVE_ANSI */
{
  int i;
  int j;

  i = *(int *) ii;
  j = *(int *) jj;

  if (multi_cost_qsort[i] > multi_cost_qsort[j] + (double) EPS_DOUBLE)
    return (1);
  else if (multi_cost_qsort[i] < multi_cost_qsort[j] - (double) EPS_DOUBLE)
    return (-1);
  else
    return (0);
}
#endif /* MULTI_MIN */

#if ASA_PARALLEL
#if HAVE_ANSI
static int
sort_parallel (const void *ii, const void *jj)
#else /* HAVE_ANSI */
static int
sort_parallel (ii, jj)
     void *ii;
     void *jj;
#endif /* HAVE_ANSI */
{
  LONG_INT i;
  LONG_INT j;

  i = *(LONG_INT *) ii;
  j = *(LONG_INT *) jj;

  if (gener_block_state_qsort[i].cost > gener_block_state_qsort[j].cost)
    return (1);
  else if (gener_block_state_qsort[i].cost < gener_block_state_qsort[j].cost)
    return (-1);
  else
    return (0);
}
#endif /* ASA_PARALLEL */
#if HAVE_ANSI
void
Exit_ASA (char *statement)
#else /* HAVE_ANSI */
void
Exit_ASA (statement)
     char *statement;
#endif /* HAVE_ANSI */
{
#if INCL_STDOUT
  printf ("\n\n*** EXIT calloc failed in ASA *** %s\n\n", statement);
#else
  ;
#endif /* INCL_STDOUT */
}
