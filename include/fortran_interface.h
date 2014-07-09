#ifndef _FORTRAN_INTERFACE_
#define _FORTRAN_INTERFACE_

#define _QUADRATURE_RELTOL_BASE_ (1e-5)

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
  int cls_are_ready; /* have the C_l's been computed? */


  // ========================================================================
  // =                       Precision parameters                           =
  // ========================================================================

  double accuracy_level; /* overall precision of CLASS */
  double quadrature_reltol; /* relative tolerance for numerical quadrature */

  // ========================================================================
  // =                       Technical parameters                           =
  // ========================================================================

  ErrorMsg error_message;        /* for error messages */
  int class_verbose;             /* verbosity level */

};


/**************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

// ==========================================================================================
// =                                      Main functions                                    =
// ==========================================================================================


int class_interface_init (
       struct class_run ** pcr
       );

int class_interface_free (
       struct class_run * pcr
       );


// =======================================================================================
// =                                    Input/output                                     =
// =======================================================================================

int class_interface_set_precision (
       struct class_run * pcr
       );

int class_interface_compute_cls (
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

int class_interface_get_cls (
      struct class_run * pcr,
      char * cl_name,
      int l_min,
      int l_max,
      double * out
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
       int module_verbose,
       int output_verbose
       );


// =========================================================================================
// =                              Compute specific stuff                                   =
// =========================================================================================

int class_interface_derived_at_z (
       struct class_run * pcr,
       double z,
       double * result
       );

int class_interface_derived (
       struct class_run * pcr,
       double * result
       );
       
int class_interface_damping_scale (
       struct class_run * pcr,
       double * k_d
       );

double damping_scale_integrand (
         double a,
         void * pcr_void
         );


// =========================================================================================
// =                                  Numerical functions                                  =
// =========================================================================================

int qromb (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      double EPS,
      double * result,
      ErrorMsg errmsg
      );
      
int trapzd (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      int n,
      double * result,
      ErrorMsg errmsg
      );
      
int polint (
      double xa[],
      double ya[],
      int n,
      double x,
      double *y,
      double *dy,
      ErrorMsg errmsg
      );

int qtrap (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      double EPS,
      double * result,
      ErrorMsg errmsg
      );
          

#ifdef __cplusplus
}
#endif

#endif
