#ifndef _FORTRAN_INTERFACE_
#define _FORTRAN_INTERFACE_

#include "class.h"

struct class_run {
  
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

int class_interface_reset_params (
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