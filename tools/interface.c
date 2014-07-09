/** @file fortra_inteface.c Documented interface module
 *
 * Guido Pettinari, 09.07.2014
 *
 * This library is meant to be used by a third party code, possibly written in
 * another language, which does not have (or does not bother) access the complicated
 * internal structure of CLASS.
 * 
 * A CLASS session is started with the function 'class_interface_init'.
 * The parameters for a CLASS run are set via the 'class_interface_set_param' function.
 * They are then loaded into CLASS using 'class_interface_load_params'.
 * CLASS is run with functiosn like 'class_interface_compute_cls'.
 * Data is accessed via 'class_interface_get_cls', which just require simple information
 * to extract the required quantity.
 *
 */

#include "interface.h"


// =================================================================================
// =                                Main functions                                 =
// =================================================================================


/**
 * Initialise the CLASS superstructure by:
 * - Allocating memory for the structure itself.
 * - Allocating memory for its params substructure.
 * - Allocating memory for CLASS internal substructures (background, thermo, pertubs ...)
 * - Setting the overall precision parameters.
 */
int class_interface_init (
       struct class_run ** ppcr
       )
{

  /* Check that the structure has not alreaby been allocated */
  if (*ppcr != NULL) {
    printf ("ERROR in %s(%d): input CLASS super-structure is already allocated", __func__, __LINE__);
    return _FAILURE_;
  }
  
  /* Initialise the superstructure for this CLASS run */
  *ppcr = malloc (sizeof (struct class_run));

  /* Check that the memory allocation went fine */
  if (*ppcr == NULL) {
    printf ("ERROR in %s(%d): could not allocate CLASS super-structure", __func__, __LINE__);
    return _FAILURE_;
  }

  /* The structure is ready to be filled */
  (*ppcr)->execution_stage = CLASS_ALLOCATED;

  /* Initialise the parameter structure */
  class_alloc ((*ppcr)->params, sizeof(struct file_content), (*ppcr)->error_message);
  class_call (parser_init((*ppcr)->params, _MAX_NUM_PARAMETERS_, "params", (*ppcr)->error_message),
    (*ppcr)->error_message, (*ppcr)->error_message);
  (*ppcr)->params->size = 0;

  /* Allocate memory for CLASS structures */
  class_call (class_interface_allocate_data (*ppcr),
    (*ppcr)->error_message, (*ppcr)->error_message);

  /* Set the overall numerical accuracy of CLASS */
  class_call (class_interface_set_precision (*ppcr), (*ppcr)->error_message, (*ppcr)->error_message);

  return _SUCCESS_;

}



/**
 * Run CLASS by calling its modules sequentially. The results will be stored in the CLASS
 * superstructure pointed by 'pcr'.
 */
int class_interface_compute_cls (
       struct class_run * pcr
       )
{
  
  if (pcr->class_verbose > 0)
    printf (" @ Running CLASS\n");
    
  class_test (pcr->execution_stage < PARAMS_SET,
    pcr->error_message,
    "cannot run CLASS: parameters have not been set");
    
  class_call (background_init(pcr->ppr,pcr->pba),
    pcr->pba->error_message, pcr->error_message);
  class_call (thermodynamics_init(pcr->ppr,pcr->pba,pcr->pth),
    pcr->pth->error_message, pcr->error_message);
  class_call (perturb_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt),
    pcr->ppt->error_message, pcr->error_message);
  class_call (primordial_init(pcr->ppr,pcr->ppt,pcr->ppm),
    pcr->pnl->error_message, pcr->error_message);
  class_call (nonlinear_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt,pcr->ppm,pcr->pnl),
    pcr->pnl->error_message, pcr->error_message);
  class_call (transfer_init(pcr->ppr,pcr->pba,pcr->pth,pcr->ppt,pcr->pnl,pcr->ptr),
    pcr->ptr->error_message, pcr->error_message);
  class_call (spectra_init(pcr->ppr,pcr->pba,pcr->ppt,pcr->ppm,pcr->pnl,pcr->ptr,pcr->psp),
    pcr->psp->error_message, pcr->error_message);
  class_call (lensing_init(pcr->ppr,pcr->ppt,pcr->psp,pcr->pnl,pcr->ple),
    pcr->ple->error_message, pcr->error_message);
  class_call (output_init(pcr->pba,pcr->pth,pcr->ppt,pcr->ppm,pcr->ptr,pcr->psp,pcr->pnl,pcr->ple,pcr->pop),
    pcr->ple->error_message, pcr->error_message);
    
  /* C_l's are ready to be extracted */
  pcr->execution_stage = CLS_COMPUTED;

  return _SUCCESS_;
    
}



// =======================================================================================
// =                                       Input                                         =
// =======================================================================================

/**
 * Set the overall numerical accuracy level of CLASS. This function is not properly
 * implemented yet. It could be useful as now there is no way to globally set
 * precision in CLASS. The question is whether such function should be here (in pcr)
 * or should rather belong to the input module.
 */
int class_interface_set_precision (
       struct class_run * pcr
       )
{

  /* Main accuracy parameters */
  pcr->accuracy_level = 1;

  /* Derived accuracy parameters */
  pcr->quadrature_reltol = _QUADRATURE_RELTOL_BASE_ / pow(10,-1+pcr->accuracy_level);

  return _SUCCESS_;
    
}


/**
 * Reset the cosmological parameters in the CLASS superstructure to their
 * default values.
 */
int class_interface_reset_cosmology_params (
       struct class_run * pcr
       )
{

  if ((pcr->class_verbose > 2) && (pcr->execution_stage >= PARAMS_SET))
    printf (" @ Reset CLASS cosmology parameters\n");

  input_default_params(pcr->pba,
                       pcr->pth,
                       pcr->ppt,
                       pcr->ptr,
                       pcr->ppm,
                       pcr->psp,
                       pcr->pnl,
                       pcr->ple,
                       pcr->pop);

  return _SUCCESS_;
    
}


/**
 * Reset the precision parameters in the CLASS superstructure to their
 * default values.
 */
int class_interface_reset_precision_params (
       struct class_run * pcr
       )
{

  if ((pcr->class_verbose > 2) && (pcr->execution_stage >= PARAMS_SET))
    printf (" @ Reset CLASS precision parameters\n");

  input_default_precision(pcr->ppr);

  return _SUCCESS_;
    
}


/**
 * Reset all parameters in CLASS to their default values, while mantaining the
 * general settings in the superstructure untouched.
 */
int class_interface_reset_params (
       struct class_run * pcr
       )
{

  if (pcr->class_verbose > 1)
    printf (" @ Reset CLASS parameters\n");

  /* Reset all parameters */
  class_call (class_interface_reset_precision_params (pcr),
    pcr->error_message, pcr->error_message);

  class_call (class_interface_reset_cosmology_params (pcr),
    pcr->error_message, pcr->error_message);

  /* CLASS is ready to be run */
  pcr->execution_stage = PARAMS_SET;

  return _SUCCESS_;
    
}


/**
 * Set the CLASS parameters 'name' to the value 'value'. The parameter names
 * should match those given in ini files (see explanatory.ini).
 */
int class_interface_set_param (
       struct class_run * pcr,
       char * name,
       char * value,
       enum entry_operation what_to_do
       )
{
  
  class_test (pcr->execution_stage < CLASS_ALLOCATED,
    pcr->error_message,
    "cannot write CLASS parameters, memory not allocated");
  
  int index = 0;
  
  class_call (parser_add_entry (
                pcr->params,
                name,
                value,
                what_to_do,
                &index,
                pcr->error_message),
    pcr->error_message,
    pcr->error_message);
    
  if (pcr->class_verbose > 9) {
    printf (" @  Entry #%d is now '%s = %s'\n",
      index, pcr->params->name[index], pcr->params->value[index]);
  }
  
  return _SUCCESS_;
  
}


/**
 * Read the parameters stored in the structure 'pcr->params' and load them
 * into CLASS substructures. The parameters are stored in 'pcr->params'
 * via the 'class_interface_set_param' function.
 */
int class_interface_load_params (
       struct class_run * pcr
       )
{
  
  class_test (pcr->execution_stage < DATA_ALLOCATED,
    pcr->error_message,
    "cannot load CLASS parameters, memory not allocated");
  
  if (pcr->class_verbose > 4) {
    printf (" @ Setting up CLASS with %d parameters:\n", pcr->params->size);
    parser_print (pcr->params);
  }
  
  /* Remove any trailing commas from the params structure */
  class_call (parser_remove_extra_commas (pcr->params), pcr->error_message, pcr->error_message);
  
  /* Launch CLASS parameter reading routine using our 'pcr->params' structure as input */
  class_call (input_read_parameters (
                pcr->params,
                pcr->ppr,
                pcr->pba,
                pcr->pth,
                pcr->ppt,
                pcr->ptr,
                pcr->ppm,
                pcr->psp,
                pcr->pnl,
                pcr->ple,
                pcr->pop,
                pcr->error_message),
    pcr->error_message,
    pcr->error_message);
    
  /* CLASS is now ready to run */
  pcr->execution_stage = PARAMS_SET;
    
  return _SUCCESS_;
  
}


/**
 * Set the level of information displayed on screen by CLASS.
 * - 'class_verbose' controls the level of information relative to the
 *   operations performed in this module on the superstructure.
 * - 'modules_verbose' controls the info printed by CLASS internal
 *   modules.
 * - 'output_verbose', if larger than zero, prints to file CLASS
 *   outputs (C_l, P(K)...).
 */
int class_interface_set_verbose (
       struct class_run * pcr,
       int class_verbose,
       int modules_verbose,
       int output_verbose
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
  pcr->pop->output_verbose = output_verbose;
  
  return _SUCCESS_;
    
}



// =======================================================================================
// =                                      Output                                         =
// =======================================================================================

/**
 * Extract the angular power spectrum C_l corresponding to 'cl_name' from the CLASS run.
 * The string 'cl_name' can include the following:
 *
 * - "TT", "TE", "EE", "BB", "PP", "TP", "EP" to choose the fields, where T=temperature,
 *   E,B=polarisation, P=lensing potential.
 * - "scalar", "vector", "tensor" to choose the mode; if none of the three labels is given,
 *   the function will output the total C_l's.
 * - "lens" for lensed C_l, nothing for unlensed C_l.
 *
 * The output power spectra DO include a l*(l+1)/(2*_PI_) factor.
 */
int class_interface_get_cls (
      struct class_run * pcr,
      char * cl_name,
      int l_min,
      int l_max,
      double * out
      )
{
  
  /* Shortcuts */
  struct perturbs * ppt = pcr->ppt;
  struct spectra * psp = pcr->psp;
  struct lensing * ple = pcr->ple;
  int want_lensing = (strstr (cl_name, "lens") != NULL);
  int want_scalar = (strstr (cl_name, "scalar") != NULL);
  int want_vector = (strstr (cl_name, "vector") != NULL);
  int want_tensor = (strstr (cl_name, "tensor") != NULL);

  /* Checks */
  class_test (pcr->execution_stage < CLS_COMPUTED,
    pcr->error_message,
    "cannot load CLASS C_l, they have not been computed yet!");

  class_test ((l_min < 2) || (l_max > psp->l_max_tot) || (l_min >= l_max),
    pcr->error_message,
    "your requested power spectrum (%s) was not computed for all the values in l=[%d,%d]",
    cl_name, l_min, l_max);
    
  class_test (want_lensing && (ple->has_lensed_cls == _FALSE_),
    pcr->error_message,
    "you requested lensed C_l from CLASS, but lensing was not computed");

  class_test (want_scalar && (ppt->has_scalars == _FALSE_),
    pcr->error_message,
    "you requested scalar C_l from CLASS, but scalars were not computed");

  class_test (want_vector && (ppt->has_vectors == _FALSE_),
    pcr->error_message,
    "you requested vector C_l from CLASS, but vectors were not computed");

  class_test (want_tensor && (ppt->has_tensors == _FALSE_),
    pcr->error_message,
    "you requested tensor C_l from CLASS, but tensors were not computed");

  class_test (want_lensing && (want_scalar || want_vector || want_tensor),
    pcr->error_message,
    "you requested C_l(%s), but with lensing we cannot single out scalar, vector or tensor contributions",
    cl_name);

  /* Determine which mode to output. When no mode is specified, we set index_mode=-1 and
  return the total C_l. */
  int index_mode = -1;

  if (want_scalar)
    index_mode = ppt->index_md_scalars;
  else if (want_vector)
    index_mode = ppt->index_md_vectors;
  else if (want_tensor)
    index_mode = ppt->index_md_tensors;
  
  /* Determine which spectrum type to output. This index will point to the column of the
  'cl' table corresponding to the requested power spectrum */
  int index_cl = -1;

       if (strstr (cl_name, "TT") != NULL) index_cl = want_lensing?ple->index_lt_tt:psp->index_ct_tt;
  else if (strstr (cl_name, "TE") != NULL) index_cl = want_lensing?ple->index_lt_te:psp->index_ct_te;
  else if (strstr (cl_name, "EE") != NULL) index_cl = want_lensing?ple->index_lt_ee:psp->index_ct_ee;
  else if (strstr (cl_name, "BB") != NULL) index_cl = want_lensing?ple->index_lt_bb:psp->index_ct_bb;
  else if (strstr (cl_name, "PP") != NULL) index_cl = want_lensing?ple->index_lt_pp:psp->index_ct_pp;
  else if (strstr (cl_name, "TP") != NULL) index_cl = want_lensing?ple->index_lt_tp:psp->index_ct_tp;
  else if (strstr (cl_name, "EP") != NULL) index_cl = want_lensing?ple->index_lt_ep:psp->index_ct_ep;
  else class_stop (pcr->error_message, "requested C_l(%s) does not exist", cl_name);
  
  /* Check that the required l_max is within bounds */
  int l_max_computed = 0;

  if (!want_lensing) {
    if (index_mode !=-1) 
      l_max_computed = psp->l_max_ct[index_mode][index_cl];
    else 
      for (int index_md=0; index_md < ppt->md_size; ++index_md)
        l_max_computed = MAX(l_max_computed, psp->l_max_ct[index_md][index_cl]);
  } 
  else {
    l_max_computed = ple->l_max_lt[index_cl];
  }

  if (pcr->class_verbose > 0)
  class_test_permissive (l_max > l_max_computed,
    pcr->error_message,
    "you asked for C_l(%s) up to l_max=%d, but it was computed only up to l=%d; will set to zero.\n",
    cl_name, l_max, l_max_computed);
  


  // =====================================================================================
  // =                            Interpolate unlensed C_l                               =
  // =====================================================================================

  if (pcr->class_verbose > 6)
    printf (" @ Extracting C_l(%s) from CLASS\n", cl_name);

  /* Array that will contain the interpolated C_l as a function of type and multipole */
  double * cl;
  class_calloc (cl, (l_max-l_min+1), sizeof(double), pcr->error_message);

  if (!want_lensing) {

    /* Allocate memory for interpolated results */
    double * cl_unlensed, ** cl_md_ic, ** cl_md;
    class_alloc(cl_unlensed, psp->ct_size*sizeof(double), pcr->error_message);
    class_alloc(cl_md_ic, psp->md_size*sizeof(double *), pcr->error_message);
    class_alloc(cl_md, psp->md_size*sizeof(double *), pcr->error_message);
    for (int index_md = 0; index_md < psp->md_size; index_md++) {
      class_alloc(cl_md[index_md], psp->ct_size*sizeof(double), pcr->error_message);
      if (psp->ic_size[index_md] > 1)
        class_alloc(cl_md_ic[index_md], psp->ic_ic_size[index_md]*psp->ct_size*sizeof(double),
          pcr->error_message);
    }

    /* Interpolate all C_l in requested l-values. If some of the C_l are not available up
    to 'l_max', they will be set to zero. */    
    for (int l=l_min; l<=l_max; l++) {

      class_call(spectra_cl_at_l(psp,l,cl_unlensed,cl_md,cl_md_ic),
                 psp->error_message,
                 pcr->error_message);
                 
      /* Pick the C_l's that we need based on the requested mode (scalar, vector, tensor) */
      cl[l-l_min] = cl_unlensed[index_cl];
      if ((psp->md_size > 1) && (index_mode >= 0)) /* important to enter here only if psp->md_size>1 */
        cl[l-l_min] = cl_md[index_mode][index_cl];
    }
    
    /* Free memory */
    for (int index_md = 0; index_md < psp->md_size; index_md++) {
      free(cl_md[index_md]);
      if (psp->ic_size[index_md] > 1)
        free(cl_md_ic[index_md]);
    }
    free(cl_md_ic);
    free(cl_md);
    free(cl_unlensed);

  } // end of unlensed C_l

  // =====================================================================================
  // =                             Interpolate lensed C_l                                =
  // =====================================================================================

  else {  /* if requesting lensed C_l */

    /* Allocate memory for interpolated results */
    double * cl_lensed;
    class_alloc(cl_lensed, ple->lt_size*sizeof(double), pcr->error_message);

    /* Interpolate all C_l in requested l-values. If some of the C_l are not available up
    to 'l_max', they will be set to zero. */    
    for (int l=l_min; l<=l_max; l++) {

      class_call(lensing_cl_at_l(ple,l,cl_lensed),
                 ple->error_message,
                 pcr->error_message);
                 
      /* When lensing is included, we cannot distinguish between scalar, vector and tensor */
      cl[l-l_min] = cl_lensed[index_cl];
    }

    free (cl_lensed);

  } // end of lensed C_l

  // ======================================================================================
  // =                              Extract requested C_l                                 =
  // ======================================================================================

  /* Write the requested power spectrum in the output array */
  for (int l=l_min; l <= l_max; ++l)
    out [l-l_min] = l*(l+1)/(2*_PI_) * cl[l-l_min];

  /* Debug - print result */
  // printf ("~~~ cl_name=%s, l_min=%d, l_max=%d\n", cl_name, l_min, l_max);
  // for (int l=l_min; l <= l_max; ++l)
  //   printf ("%4d %24.12g\n", l, out[l-l_min]);
  

  /* Free memory */
  free(cl);

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



// =====================================================================================
// =                              Memory management                                    =
// =====================================================================================

/**
 * Allocate memory for the internal structures of CLASS.
 */
int class_interface_allocate_data (
       struct class_run * pcr
       )
{
  
  /* Check that the structures are not already allocated */
  class_test (pcr->execution_stage >= DATA_ALLOCATED,
    pcr->error_message,
    "stopping to avoid memory leakage");

  /* Allocate memory for the substructures in cr */
  class_alloc (pcr->ppr, sizeof(struct precision), pcr->error_message);   /* precision parameters */
  class_alloc (pcr->pba, sizeof(struct background), pcr->error_message);  /* cosmological background */
  class_alloc (pcr->pth, sizeof(struct thermo), pcr->error_message);      /* thermodynamics */
  class_alloc (pcr->ppt, sizeof(struct perturbs), pcr->error_message);    /* source functions */
  class_alloc (pcr->ptr, sizeof(struct transfers), pcr->error_message);   /* transfer functions */
  class_alloc (pcr->ppm, sizeof(struct primordial), pcr->error_message);  /* primordial spectra */
  class_alloc (pcr->psp, sizeof(struct spectra), pcr->error_message);     /* output spectra */
  class_alloc (pcr->pnl, sizeof(struct nonlinear), pcr->error_message);   /* non-linear spectra */
  class_alloc (pcr->ple, sizeof(struct lensing), pcr->error_message);     /* lensed spectra */
  class_alloc (pcr->pop, sizeof(struct output), pcr->error_message);      /* output files */

  /* CLASS structures are ready to be filled */
  pcr->execution_stage = DATA_ALLOCATED;

  return _SUCCESS_;

}


/**
 * Deallocate CLASS internal structures by sequentially calling the xxx_free
 * functions in each module, while leaving unchanged the information in the
 * top level of the CLASS superstructure, such as verbosity levels and parameters.
 * This function must be called between successive calls to functions that compute
 * stuff, like 'class_interface_compute_cls', in order to avoid memory leakage.
 */
int class_interface_free_data (
       struct class_run * pcr
       )
{
  
  /* Do not free data if not needed */
  if (pcr->execution_stage < CLS_COMPUTED)
    goto update_and_return;
  
  if (pcr->class_verbose > 2)
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
  
  update_and_return:
  pcr->execution_stage = DATA_ALLOCATED;

  return _SUCCESS_;
    
}


/**
 * Free the CLASS superstructure and all the data therein.
 */
int class_interface_free (
       struct class_run ** ppcr
       )
{
  
  if ((*ppcr)->class_verbose > 4)
    printf (" @ Freeing CLASS\n");

  /* Empty CLASS structure */
  if ((*ppcr)->execution_stage >= CLS_COMPUTED)
    class_call (class_interface_free_data (*ppcr),
      (*ppcr)->error_message, (*ppcr)->error_message);
      
  /* Free CLASS actual structures */
  free ((*ppcr)->ppr);
  free ((*ppcr)->pba);
  free ((*ppcr)->pth);
  free ((*ppcr)->ppt);
  free ((*ppcr)->ptr);
  free ((*ppcr)->ppm);
  free ((*ppcr)->psp);
  free ((*ppcr)->pnl);
  free ((*ppcr)->ple);
  free ((*ppcr)->pop);

  /* Free parameter structure */
  (*ppcr)->params->size = _MAX_NUM_PARAMETERS_;
  parser_free((*ppcr)->params);
  free ((*ppcr)->params);

  /* Free super structure */
  free (*ppcr);

  return _SUCCESS_;
    
}



// =========================================================================================
// =                              Compute specific stuff                                   =
// =========================================================================================

/**
 * Compute the quantities needed for the BAO likelihood at the given redshift
 * 'z' for the given cosmology 'pcr'.
 *
 * Return the result as an array of three elements (must be preallocated):
 * result[0] = ratio between the sound horizon at baryon drag and the volume distance (R_s/D_v).
 * result[1] = expansion rate (H=100*h) 
 * result[2] = the angular diameter distance (D_a)
 */
int class_interface_derived_at_z (
       struct class_run * pcr,
       double z,
       double * result
       )
{
  
  class_test (pcr->execution_stage < THERMODYNAMICS_COMPUTED,
    pcr->error_message,
    "needs thermodynamics");
  
  struct background * pba = pcr->pba;
  struct thermo * pth = pcr->pth;
  
  /* Comoving sound horizon at baryon drag */
  double R_s = pth->rs_d;
  
  /* Find tau(z) */
  double tau;
  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pcr->error_message);
  
  /* Interpolate D_a(z) and H(z) */
  double * pvecback;
  class_alloc (pvecback, pba->bg_size*sizeof(double), pcr->error_message);
  int junk;
  class_call(background_at_tau(pba,
                               tau,
                               pba->long_info,
                               pba->inter_normal,
                               &junk,
                               pvecback),
             pba->error_message,
             pcr->error_message);
             
  double D_a = pvecback[pba->index_bg_ang_distance];
  double H = pvecback[pba->index_bg_H];
  
  /* Compute volume distance D_v(z), taken from CosmoMC */
  double ADD = D_a*(1+z);
  double D_v = pow((ADD*ADD)*z/H,1./3);

  /* Build result array */
  result[0] = R_s/D_v;
  result[1] = H * _c_ * 1e-3;
  result[2] = D_a;
  
  free (pvecback);
  
  /* Debug - print results */
  // printf ("CLASS: z=%g, Rs/Dv=%g\n", z, result[0]);
  // printf ("CLASS: z=%g, H=%g\n", z, result[1]);
  // printf ("CLASS: z=%g, Da=%g\n", z, result[2]);
  // printf ("CLASS: z=%g, rs=%g\n", z, R_s);
  
  return _SUCCESS_;
  
}

/**
 * Compute the ten derived parameters needed by CosmoMC. These are hard-coded in
 * CAMB (inithermo function, modules.f90):
 * Age, zstar, rstar, thetastar, zdrag, rdrag, kD, thetaD , zEQ, thetaEQ,
 * where 'star' stands for recombination, 'drag' for end of baryon drag, 'D'
 * for damping and 'EQ' for matter-radiation equality (including massive)
 * neutrinos.
 *
 * Return the result as an array of ten elements, in the above order (must be
 * preallocated).
 */
int class_interface_derived (
       struct class_run * pcr,
       double * result
       )
{
  
  class_test (pcr->execution_stage < THERMODYNAMICS_COMPUTED,
    pcr->error_message,
    "needs thermodynamics");

  struct background * pba = pcr->pba;
  struct thermo * pth = pcr->pth;

  /* Compute damping scale k_d */
  double k_d;
  class_call (class_interface_damping_scale (pcr, &k_d), pcr->error_message, pcr->error_message);

  /* Build result array */
  result[0] = pba->age; /* Age of the Universe in Gyr */
  result[1] = pth->z_rec; /* Redshift at the peak of recombination */
  result[2] = pth->rs_rec; /* Comoving sound horizon at recombination */
  result[3] = 100*pth->rs_rec/pth->ra_rec; /* Recombination angular scale today */
  result[4] = pth->z_d; /* Redshift at the end of baryon drag */
  result[5] = pth->rs_d; /* Comoving sound horizon at baryon drag */
  result[6] = k_d; /* Damping scale */
  result[7] = 100*_PI_/k_d/pth->ra_rec; /* Damping angular scale today */
  result[8] = 1/pba->a_eq-1; /* Redshift at matter radiation equality */
  result[9] = 100*(1/pba->k_eq)/pth->ra_rec; /* Equality angular scale today (why da_rec and not da_eq?) */
  
  return _SUCCESS_;
  
}


/**
 * Compute the damping scale k_d at recombination, needed by CosmoMC as a parameter.
 */
int class_interface_damping_scale (
       struct class_run * pcr,
       double * k_d
       )
{
 
  class_test (pcr->execution_stage < THERMODYNAMICS_COMPUTED,
    pcr->error_message,
    "needs thermodynamics");
 
   /* Here we use 'qromb' on the interpolated integrand function. I think this is not
  ideal: if we use interpolation, we could just sum over the node points rather than
  treating the interpolated function like an analytical one! */
  
  double a_ini = pcr->ppr->a_ini_over_a_today_default * pcr->pba->a_today;
  double a_end = 1/(pcr->pth->z_rec+1);
  
  class_call (qromb (
                &damping_scale_integrand,
                (void *)pcr,
                a_ini,
                a_end,
                pcr->quadrature_reltol,
                k_d, /* result */
                pcr->error_message),
    pcr->error_message,
    pcr->error_message);

  /* In the tau domain */
  // class_call (qromb (
  //               &damping_scale_integrand_of_tau,
  //               (void *)pcr,
  //               1,
  //               tau_today,
  //               pcr->quadrature_reltol,
  //               k_d, /* result */
  //               pcr->error_message),
  //   pcr->error_message,
  //   pcr->error_message);

  *k_d = sqrt(1/(*k_d));
    
  class_test (fabs(*k_d) < (1/_HUGE_),
    pcr->error_message,
    "damping scale k_d too small! k_d=%g", *k_d);
    
  class_test (*k_d < 0,
    pcr->error_message,
    "damping scale k_d negative! k_d=%g", *k_d);

  
  return _SUCCESS_;
  
}

double damping_scale_integrand_of_tau (
         double tau,
         void * pcr_void
         )
{

  /* Shortcuts */
  struct class_run * pcr = (struct class_run *)pcr_void;
  struct background * pba = pcr->pba;
  struct thermo * pth = pcr->pth;

  /* Find R = 3/4 rho_b/rho_g by interpolation. Could do it analytically but interpolation
  is needed anyway to obtain the thermodynamical quantities. */
  double * pvecback;
  class_alloc (pvecback, pba->bg_size*sizeof(double), pcr->error_message);
  int junk;
  class_call(background_at_tau(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_normal,
                               &junk,
                               pvecback),
             pba->error_message,
             pcr->error_message);

  double a = pvecback[pba->index_bg_a];
  double R = 3./4. * pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];
  double dtau_da = pvecback[pba->index_bg_H_prime]/a;

  /* Find the interaction rate dkappa_dtau per unit of conformal time */
  double * pvecthermo;
  class_alloc (pvecthermo, pth->th_size*sizeof(double), pcr->error_message);
  double z = 1/a-1;
  class_call(thermodynamics_at_z(pba,
                                 pth,
                                 z,
                                 pth->inter_normal,
                                 &junk,
                                 pvecback,
                                 pvecthermo),
             pth->error_message,
             pcr->error_message);

  double kappa_dot = pvecthermo[pth->index_th_dkappa];

  /* Build the integrand as it is done in CAMB */
  double integrand = (R*R + 16*(1+R)/15.)/pow(1+R,2) / kappa_dot / 6.;

  /* Debug */
  printf ("%17.7g %17.7g %17.7g %17.7g %17.7g\n", tau, integrand, R, dtau_da, kappa_dot);

  // ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(a)*a**2/(Recombination_xe(a)*akthom)

  return integrand;

}

/**
 * Support function for 'class_interface_damping_scale'. Return the integrand
 * of the damping scale integral for a given value of the scale factor.
 * Here we use the same approach of CAMB (function inithermo, modules.f90); the
 * integrand differs from Eq. 8.40 page 232 of Dodelson's Modern Cosmology as it
 * has a different coefficient in the square brackets (16/15 rather than 8/9).
 *
 */
double damping_scale_integrand (
         double a,
         void * pcr_void
         )
{

  /* Shortcuts */
  struct class_run * pcr = (struct class_run *)pcr_void;
  struct background * pba = pcr->pba;
  struct thermo * pth = pcr->pth;
  double z = 1/a-1;

  /* Find tau(a) */
  double tau;
  class_call(background_tau_of_z(pba,z,&tau),
             pba->error_message,
             pcr->error_message);

  /* Find R = 3/4 rho_b/rho_g by interpolation. Could do it analytically but interpolation
  is needed anyway to obtain the thermodynamical quantities. */
  double * pvecback;
  class_alloc (pvecback, pba->bg_size*sizeof(double), pcr->error_message);
  int junk;
  class_call(background_at_tau(pba,
                               tau,
                               pba->normal_info,
                               pba->inter_normal,
                               &junk,
                               pvecback),
             pba->error_message,
             pcr->error_message);

  double R = 3./4. * pvecback[pba->index_bg_rho_b]/pvecback[pba->index_bg_rho_g];
  double dtau_da = 1/(a*a*pvecback[pba->index_bg_H]); /* da_dtau = a*Hc = a*a*H */

  /* Find the interaction rate dkappa_dtau per unit of conformal time */
  double * pvecthermo;
  class_alloc (pvecthermo, pth->th_size*sizeof(double), pcr->error_message);
  class_call(thermodynamics_at_z(pba,
                                 pth,
                                 z,
                                 pth->inter_normal,
                                 &junk,
                                 pvecback,
                                 pvecthermo),
             pth->error_message,
             pcr->error_message);

  double dkappa_dtau = pvecthermo[pth->index_th_dkappa];

  /* Build the integrand as it is done in CAMB. */
  double integrand = (R*R + 16*(1+R)/15.)/pow(1+R,2) * dtau_da / dkappa_dtau / 6.;

  /* Debug */
  // printf ("%17.7g %17.7g %17.7g %17.7g %17.7g\n", a, integrand, R, dtau_da, dkappa_dtau);

  /* Camb version: */
  // ddamping_da = (R**2 + 16*(1+R)/15)/(1+R)**2*dtauda(a)*a**2/(Recombination_xe(a)*akthom)

  return integrand;

}
  
  
// =========================================================================================
// =                                  Numerical functions                                  =
// =========================================================================================

/**
 * FROM NUMERICAL RECIPES, 2ND EDITION, PAGE 140:
 * Returns the integral of the function func from a to b. Integration is performed by
 * Romberg’s method of order 2K, where, e.g., K=2 is Simpson’s rule. 
 * Here EPS is the fractional accuracy desired, as determined by the extrapolation error
 * estimate; JMAX limits the total number of steps; K is the number of points used in
 * the extrapolation.
 *
 * UPDATE BY GWP:
 * This version of 'qromb' is improved with respect to the NR original, because the
 * integrand function can accept parameters. This is achieved via the extra argument
 * 'func_params'. (This is a much simplified version of what is done in GSL with the
 * gsl_function structure.) The extra argument is a 'double *' that can contain any number
 * of arguments, or NULL for no arguments. It is the caller's responsibility to make sure
 * that this array is properly allocated, as the integrating routine does not know how
 * many parameters the integrand function has.
 *
 */
int qromb (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      double EPS,
      double * result,
      ErrorMsg errmsg
      )
{

  int JMAX = 20;
  int JMAXP = JMAX+1;
  int K = 5;

  double ss, dss;
  double s[JMAXP], h[JMAXP+1];   /* Successive trapezoidal approximations and relative stepsizes. */
  int j;

  h[1]=1.0;

  for (j=1;j<=JMAX;j++) {

    class_call (trapzd (func,func_params,a,b,j,&s[j],errmsg), errmsg, errmsg);

    if (j >= K) {

      class_call ( polint (&h[j-K], &s[j-K], K, 0.0, &ss, &dss, errmsg),
         errmsg, errmsg);

      if (fabs(dss) <= EPS*fabs(ss)) {
        *result = ss;
        return _SUCCESS_;
      }
    }

    /* This is a key step: The factor is 0.25 even though the stepsize is decreased by
    only 0.5. This makes the extrapolation a polynomial in h2 as allowed by equation
    (4.2.1), not just a polynomial in h. */
    h[j+1]=0.25*h[j];
  }

  /* Should not be here */
  class_stop(errmsg, "Too many steps in routine qromb");

}
  
  
/**
 * FROM NUMERICAL RECIPES, 2ND EDITION, PAGE 137.
 * REFINED USING GSL LIBRARY, FILE INTEGRATION.C (no more use of static variable).
 * This routine computes the n-th stage of refinement of an extended trapezoidal rule.
 * The function pointer 'func'is the function to be integrated between limits 'a' and 'b',
 * also input. When called with n=1, the routine returns the crudest estimate of int_a^b(f(x)dx).
 *
 * Subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy by
 * adding 2^(n-2) additional interior points.  The code-fragment for doing so looks like:
 *                                                                      
 *       double answer;                                                 
 *       for(j=1; j<=M+1; j++)                                          
 *         trapezoid_rule(func, a, b, j, &answer);
 * 
 * It is however better to use 'qtrap' because it does the same but has a convergence criterion.
 *
 * UPDATE BY GWP:
 * This version of 'trapzd' is improved with respect to the NR original, because the
 * integrand function can accept parameters. This is achieved via the extra argument
 * 'func_params'. (This is a much simplified version of what is done in GSL with the
 * gsl_function structure.) The extra argument is a 'void *' that can be anything,
 * or NULL for no arguments. It is the caller's responsibility to make sure
 * that this array is properly working.
 *
 */
int trapzd (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      int n,
      double * result, /* result */
      ErrorMsg errmsg
      )
{
  double x, tnm, sum, del;
  int it, j;
  double s = 0.0;

  if(n==1){
    s = 0.5 * (b-a) * ((*func)(b,func_params) + (*func)(a,func_params));
  }
  else {

    for(it=1, j=1; j < n-1; j++)
      it <<= 1;

    tnm = (double) it;
    del = (b-a) / tnm; /* This is the spacing of the points to be added.  */
    x = a + 0.5 * del;
    
    for(sum=0., j=1; j<=it; j++, x+=del)
      sum += (*func)(x,func_params);

    s = 0.5 * (s + del * sum); /* This replaces s by its refined value. */
    
  }
  
  *result = s;
  
  return _SUCCESS_;
}


/**
 * FROM NUMERICAL RECIPES, 2ND EDITION, PAGE 137.
 * Returns the integral of the function func from a to b. The parameters EPS can
 * be set to the desired fractional accuracy and JMAX so that 2 to the power JMAX-1
 * is the maximum allowed number of steps. Integration is performed by the trapezoidal
 * rule.
 *
 * UPDATE BY GWP:
 * This version of 'qtrap' is improved with respect to the NR original, because the
 * integrand function can accept parameters. This is achieved via the extra argument
 * 'func_params'. (This is a much simplified version of what is done in GSL with the
 * gsl_function structure.) The extra argument is a 'double *' that can contain any number
 * of arguments, or NULL for no arguments. It is the caller's responsibility to make sure
 * that this array is properly allocated, as the integrating routine does not know how
 * many parameters the integrand function has.
 *
 * ATTENTION: NEVER TESTED BEFORE, JUST COPIED DOWN FROM THE BOOK!
 * 
 */
int qtrap (double (*func)(double, void *),
      double * func_params,
      double a,
      double b,
      double EPS,
      double * result, /* result */
      ErrorMsg errmsg
      )
{
  
  int JMAX = 20;
  int j;
  double s, olds;

  olds = -1.0e30; /* Any number unlikely to be the function average at its endpoints will do. */

  for (j=1; j<=JMAX; j++) {

    class_call (trapzd(func,func_params,a,b,j,&s,errmsg), errmsg, errmsg);

    if (j > 5) { /* Avoid spurious early convergence */
      if (fabs(s-olds) < EPS*fabs(olds) || (s == 0.0 && olds == 0.0)) {
        *result = s;
        return _SUCCESS_;
      }
    }

    olds = s;
  }

  /* Should not be here */
  class_stop (errmsg, "Too many steps in routine qtrap");

}

/**
 * FROM NUMERICAL RECIPES, 2ND EDITION, PAGE 109:
 * Routine for polynomial interpolation or extrapolation from N input points.
 * Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns
 * a value y, and an error estimate dy. If P (x) is the polynomial of degree N − 1
 * such that P (xai) = yai, i = 1,...,n, then the returned value y = P(x).
 * 
 * Note that the input arrays are assumed to be unit-offset. If you have zero-offset
 * arrays, remember to subtract 1.
 */
int polint (
      double xa[],
      double ya[],
      int n,
      double x,
      double *y,
      double *dy,
      ErrorMsg errmsg)
{
  
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  double *c,*d;
  
  dif=fabs(x-xa[1]);
  class_alloc (c, sizeof(double)*(n+1), errmsg);
  class_alloc (d, sizeof(double)*(n+1), errmsg);
  
  for (i=1;i<=n;i++) { /* Here we find the index ns of the closest table entry, */
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i]; /* and initialize the tableau of c’s and d’s. */
    d[i]=ya[i];
  }
  *y=ya[ns--]; /* This is the initial approximation to y. */
  for (m=1;m<n;m++) { /* For each column of the tableau, */
    for (i=1;i<=n-m;i++) { /* we loop over the current c’s and d’s and update them. */
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      class_test ((den=ho-hp) == 0.0, errmsg,
        "This error can occur only if two input xa’s are (to within roundoff) identical.");
      den=w/den;
      d[i]=hp*den; /* Here the c’s and d’s are updated. */
      c[i]=ho*den;
  }
  *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  /* After each column in the tableau is completed, we decide which correction, c or d, we
  want to add to our accumulating value of y, i.e., which path to take through the tableau—forking
  up or down. We do this in such a way as to take the most “straight line” route through the
  tableau to its apex, updating ns accordingly to keep track of where we are. This route keeps
  the partial approximations centered (insofar as possible) on the target x. The last dy added
  is thus the error indication. */
  }
  free(c);
  free(d);
  
  return _SUCCESS_;
  
}
