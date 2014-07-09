#ifndef _INTERFACE_
#define _INTERFACE_


#include "class.h"

/* Maximum number of parameters that can be stored in the parameter structure. */
#define _MAX_NUM_PARAMETERS_ (4096)

/* Baseline values for setting precision parameters in CLASS */
#define _QUADRATURE_RELTOL_BASE_ (1e-5)

/* Flags that represent a certain stage in the execution workflow of CLASS */
enum execution_stages {
  CLASS_ALLOCATED = 0,
  DATA_ALLOCATED = 1,
  PARAMS_SET = 2,
  BACKGROUND_COMPUTED = 3,
  THERMODYNAMICS_COMPUTED = 4,
  TRANSFERS_COMPUTED = 5,
  CLS_COMPUTED = 6
};

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

  enum execution_stages execution_stage; /* keep track of the status of execution */


  // ========================================================================
  // =                       Precision parameters                           =
  // ========================================================================

  double accuracy_level; /* overall precision of CLASS */
  double quadrature_reltol; /* relative tolerance for numerical quadrature */

  // ========================================================================
  // =                       Technical parameters                           =
  // ========================================================================

  ErrorMsg error_message;        /* for error messages */
  int class_verbose;             /* CLASS super-structure verbosity level */
  int modules_verbose;           /* Internal modules verbosity level */ /* TODO */
  int output_verbose;            /* File production verbosity level */  /* TODO */

};




/**************************************************************/

/*
 * Boilerplate for C++
 */
#ifdef __cplusplus
extern "C" {
#endif

// =================================================================================
// =                                Main functions                                 =
// =================================================================================


int class_interface_init (
       struct class_run ** pcr
       );


int class_interface_compute_cls (
       struct class_run * pcr
       );


// =====================================================================================
// =                              Memory management                                    =
// =====================================================================================


int class_interface_allocate_data (
       struct class_run * pcr
       );

int class_interface_free_data (
       struct class_run * pcr
       );

int class_interface_free (
       struct class_run ** ppcr
       );



// =======================================================================================
// =                                    Input/output                                     =
// =======================================================================================

int class_interface_set_precision (
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
