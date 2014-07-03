#ifndef _FORTRAN_INTERFACE_
#define _FORTRAN_INTERFACE_

#include "class.h"

/* This library is meant to be called by a third party code, possibly written in
another language, which does not have (or does not bother) access the complicated
internal structure of CLASS. The parameters for a CLASS run are set either
directly via the 'class_interface_set' function or indirectly by filling
a 'file_content' structure. The latter method is preferred. The below number
is the maximum number of parameters that can be stored in such structure. */
#define _MAX_NUM_PARAMETERS_ 4096

struct class_run {
  
  // =========================================================================
  // =                          CLASS structures                             =
  // =========================================================================
  
  struct precision * ppr;        /* for precision parameters */
  struct background * pba;       /* for cosmological background */
  struct thermo * pth;           /* for thermodynamics */
  struct perturbs * ppt;         /* for source functions */
  struct transfers * ptr;        /* for transfer functions */
  struct primordial * ppm;       /* for primordial spectra */
  struct spectra * psp;          /* for output spectra */
  struct nonlinear * pnl;        /* for non-linear spectra */
  struct lensing * ple;          /* for lensed spectra */
  struct output * pop;           /* for output files */


  // ========================================================================
  // =                        Input of parameters                           =
  // ========================================================================
  
  struct file_content * params; /* input parameters for the run */


  // ========================================================================
  // =                            Status flags                              =
  // ========================================================================

  int structs_are_allocated; /* have CLASS structures been initialised correctly? */
  int params_are_set; /* have the parameters been set to their default values? */


  // ========================================================================
  // =                       Technical parameters                           =
  // ========================================================================

  ErrorMsg error_message;        /* for error messages */
  int class_verbose;             /* verbosity level */ 

};

int class_interface_init (
       struct class_run ** pcr
       );

int class_interface_free (
       struct class_run * pcr
       );

int class_interface_compute (
       struct class_run * pcr
       );

int class_interface_reset_cosmology_params (
       struct class_run * pcr
       );

int class_interface_reset_precision_params (
       struct class_run * pcr
       );

int class_interface_reset (
       struct class_run * pcr
       );

int class_interface_set_param (
       struct class_run * pcr,
       char * name,
       char * value,
       enum entry_operation what_to_do
       );

int class_interface_load_params (
       struct class_run * pcr
       );

int class_interface_print_error (
       struct class_run * pcr
       );

int class_interface_set_verbose (
       struct class_run * pcr,
       int class_verbose,
       int module_verbose
       );

#endif