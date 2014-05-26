#include "fortran_interface.h"


/**
 * Initialise the CLASS superstructure. First allocate memory, then fill the substructures
 * with default values for the cosmological and precision parameters.
 */
int class_interface_init (
       struct class_run ** ppcr
       )
{

  /* Allocate memory for the substructures in cr */
  class_alloc ((*ppcr), sizeof(struct class_run), (*ppcr)->error_message);        /* precision parameters */
  class_alloc ((*ppcr)->ppr, sizeof(struct precision), (*ppcr)->error_message);   /* precision parameters */
  class_alloc ((*ppcr)->pba, sizeof(struct background), (*ppcr)->error_message);  /* cosmological background */
  class_alloc ((*ppcr)->pth, sizeof(struct thermo), (*ppcr)->error_message);      /* thermodynamics */
  class_alloc ((*ppcr)->ppt, sizeof(struct perturbs), (*ppcr)->error_message);    /* source functions */
  class_alloc ((*ppcr)->ptr, sizeof(struct transfers), (*ppcr)->error_message);   /* transfer functions */
  class_alloc ((*ppcr)->ppm, sizeof(struct primordial), (*ppcr)->error_message);  /* primordial spectra */
  class_alloc ((*ppcr)->psp, sizeof(struct spectra), (*ppcr)->error_message);     /* output spectra */
  class_alloc ((*ppcr)->pnl, sizeof(struct nonlinear), (*ppcr)->error_message);   /* non-linear spectra */
  class_alloc ((*ppcr)->ple, sizeof(struct lensing), (*ppcr)->error_message);     /* lensed spectra */
  class_alloc ((*ppcr)->pop, sizeof(struct output), (*ppcr)->error_message);      /* output files */

  /* Set CLASS cosmological & precision parameters to their default values */
  class_interface_reset_params (*ppcr);

  return _SUCCESS_;
    
}


/**
 * Empty and deallocate the CLASS superstructure by sequentially calling the xxx_free
 * functions in each module.
 */
int class_interface_free (
       struct class_run * pcr
       )
{
  
  if (pcr->class_verbose > 1)
    printf (" @ Freeing CLASS structures\n");
  
  /* Empty substructures in pcr */
  class_call (lensing_free(pcr->ple), pcr->ple->error_message, pcr->error_message);
  class_call (spectra_free(pcr->psp), pcr->psp->error_message, pcr->error_message);
  class_call (transfer_free(pcr->ptr), pcr->ptr->error_message, pcr->error_message);
  class_call (nonlinear_free(pcr->pnl), pcr->pnl->error_message, pcr->error_message);
  class_call (primordial_free(pcr->ppm), pcr->ppm->error_message, pcr->error_message);
  class_call (perturb_free(pcr->ppt), pcr->ppt->error_message, pcr->error_message);
  class_call (thermodynamics_free(pcr->pth), pcr->pth->error_message, pcr->error_message);
  class_call (background_free(pcr->pba), pcr->pba->error_message, pcr->error_message);
  
  /* Deallocate substructures in pcr */
  free (pcr->ppr);
  free (pcr->pba);
  free (pcr->pth);
  free (pcr->ppt);
  free (pcr->ptr);
  free (pcr->ppm);
  free (pcr->psp);
  free (pcr->pnl);
  free (pcr->ple);
  free (pcr->pop);

  return _SUCCESS_;
    
}


/**
 * Run CLASS by calling its modules sequentially. The results will be stored in the CLASS
 * superstructure pointed by 'pcr'.
 */
int class_interface_compute (
       struct class_run * pcr
       )
{
  
  if (pcr->class_verbose > 0)
    printf (" @ Running CLASS\n");
  
  class_call (background_init(pcr->ppr,pcr->pba),pcr->pba->error_message, pcr->error_message);
  class_call (thermodynamics_init(pcr->ppr,pcr->pba,pcr->pth), pcr->pth->error_message, pcr->error_message);
  class_call (perturb_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt), pcr->ppt->error_message, pcr->error_message);
  class_call (primordial_init(pcr->ppr,pcr->ppt,pcr->ppm), pcr->pnl->error_message, pcr->error_message);
  class_call (nonlinear_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt,pcr->ppm,pcr->pnl), pcr->pnl->error_message, pcr->error_message);
  class_call (transfer_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt,pcr->pnl,pcr->ptr), pcr->ptr->error_message, pcr->error_message);
  class_call (spectra_init(pcr->ppr,pcr->pba,pcr->ppt,pcr->ppm,pcr->pnl,pcr->ptr,pcr->psp), pcr->psp->error_message, pcr->error_message);
  class_call (lensing_init(pcr->ppr,pcr->ppt,pcr->psp,pcr->pnl,pcr->ple), pcr->ple->error_message, pcr->error_message);
  
  return _SUCCESS_;
    
}




/**
 * Reset the cosmological and precision parameters in the CLASS superstructure to their
 * default values.
 */
int class_interface_reset_params (
       struct class_run * pcr
       )
{

  /* Set CLASS cosmological parameters to their default values */
  input_default_params(pcr->pba,
                       pcr->pth,
                       pcr->ppt,
                       pcr->ptr,
                       pcr->ppm,
                       pcr->psp,
                       pcr->pnl,
                       pcr->ple,
                       pcr->pop);
  
  /* Set CLASS precision parameters to their default values */
  input_default_precision(pcr->ppr);

  if (pcr->class_verbose > 1)
    printf (" @ Reset CLASS parameters\n");

  return _SUCCESS_;
    
}



/**
 * Print to screen the error message stored in the CLASS superstructure, inside
 * pcr->error_message.
 */
int class_interface_print_error (
       struct class_run * pcr
       )
{
  
  printf ("%s\n", pcr->error_message);
  
  return _SUCCESS_;
    
}


/**
 * Set the level of information displayed on screen by CLASS.
 */
int class_interface_set_verbose (
       struct class_run * pcr,
       int class_verbose,
       int modules_verbose
       )
{

  pcr->class_verbose = class_verbose;
  pcr->pba->background_verbose = modules_verbose;
  pcr->pth->thermodynamics_verbose = modules_verbose;
  pcr->ppt->perturbations_verbose = modules_verbose;
  pcr->ptr->transfer_verbose = modules_verbose;
  pcr->ppm->primordial_verbose = modules_verbose;
  pcr->psp->spectra_verbose = modules_verbose;
  pcr->pnl->nonlinear_verbose = modules_verbose;
  pcr->ple->lensing_verbose = modules_verbose;
  pcr->pop->output_verbose = modules_verbose;
  
  return _SUCCESS_;
    
}






