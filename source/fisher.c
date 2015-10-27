/** @file fisher.c
 *
 * Compute the observability of the CMB bispectra.
 *
 * Created by Guido W. Pettinari on 19.07.2012.
 * Last modified by Guido W. Pettinari on 23.10.2015
 */

#include "fisher.h"


int fisher_init (
     struct precision * ppr,
     struct background * pba,
     struct thermo * pth,
     struct perturbs * ppt,
     struct bessels * pbs,
     struct transfers * ptr,
     struct primordial * ppm,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     struct fisher * pfi
     )
{


  /* Check whether we need to compute spectra at all */  
  if ((pbi->has_bispectra == _FALSE_) || (pfi->has_fisher == _FALSE_)) {
  
    pfi->has_fisher = _FALSE_;
  
    printf_log_if (pfi->fisher_verbose, 0, 
      "No forecasts requested. Fisher module skipped.\n");
  
    return _SUCCESS_;
  }
  else {
    printf_log_if (pfi->fisher_verbose, 0, 
      "Computing Fisher matrix\n");
  }


  // ====================================================================================
  // =                                  Preparations                                    =
  // ====================================================================================
  
  /* Initialize indices & arrays in the Fisher structure */
  
  class_call (fisher_indices(
                ppr,pba,ppt,pbs,ptr,
                ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
    

  /* Load intrinsic bispectra from disk if they were not already computed. Note that
  we are only loading the non-separable and intrinsic bispectra, since the primary ones
  are very quick to recompute. */

  if (ppr->load_bispectra_from_disk == _TRUE_) {
    
    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt)   
      if ((pbi->bispectrum_type[index_bt] == non_separable_bispectrum)
       || (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum))
        class_call (bispectra_load_from_disk (
                      pbi,
                      index_bt),
          pbi->error_message,
          pfi->error_message);
  }


  // ====================================================================================
  // =                          Compute experimental noise                              =
  // ====================================================================================

  /* Compute the noise associated to the considered experiment */
  
  class_call (fisher_noise(
                ppr,pba,ppt,pbs,ptr,
                ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
  
  

  // ====================================================================================
  // =                              Compute cross C_l                                   =
  // ====================================================================================
    
  /* Compute C_l needed to build the covariance matrix */
    
  class_call (fisher_cross_cls(
                ppr,pba,ppt,pbs,ptr,
                ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
  


  // ====================================================================================
  // =                               Compute Fisher matrix                              =
  // ====================================================================================
  
  /* Compute the Fisher matrix, add the effect of lensing variance, store the
  results and print them */
  
  class_call (fisher_compute(
                ppr,pba,ppt,pbs,ptr,
                ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);


  return _SUCCESS_;

}











/**
 * This routine frees all the memory space allocated by fisher_init().
 *
 * @param pfi Input: pointer to fisher structure (which fields must be freed)
 * @return the error status
 */

int fisher_free(
     struct bispectra * pbi,
     struct fisher * pfi
     )
{

  if (pfi->has_fisher == _TRUE_) {

    for (int X = 0; X < pfi->ff_size; ++X) {
      for (int Y = 0; Y < pfi->ff_size; ++Y) {
        for (int Z = 0; Z < pfi->ff_size; ++Z) {        
          /* lmax */
          for (int index_l3=0; index_l3 < pfi->l3_size; ++index_l3) {          
            for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
              free (pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft]);
              free (pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l3][index_ft]);          
            }
            free (pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3]);
            free (pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l3]);
          }
          /* lmin */
          for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {          
            for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
              free (pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft]);
              free (pfi->fisher_matrix_XYZ_lmin[X][Y][Z][index_l1][index_ft]);          
            }
            free (pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1]);
            free (pfi->fisher_matrix_XYZ_lmin[X][Y][Z][index_l1]);
          }
          /* common */
          free (pfi->fisher_matrix_XYZ_largest[X][Y][Z]);
          free (pfi->fisher_matrix_XYZ_lmax[X][Y][Z]);
          free (pfi->fisher_matrix_XYZ_smallest[X][Y][Z]);
          free (pfi->fisher_matrix_XYZ_lmin[X][Y][Z]);
        }
        free (pfi->fisher_matrix_XYZ_largest[X][Y]);
        free (pfi->fisher_matrix_XYZ_lmax[X][Y]);
        free (pfi->fisher_matrix_XYZ_smallest[X][Y]);
        free (pfi->fisher_matrix_XYZ_lmin[X][Y]);
      }
      free (pfi->fisher_matrix_XYZ_largest[X]);
      free (pfi->fisher_matrix_XYZ_lmax[X]);
      free (pfi->fisher_matrix_XYZ_smallest[X]);
      free (pfi->fisher_matrix_XYZ_lmin[X]);
    }
    free (pfi->fisher_matrix_XYZ_largest);
    free (pfi->fisher_matrix_XYZ_lmax);
    free (pfi->fisher_matrix_XYZ_smallest);
    free (pfi->fisher_matrix_XYZ_lmin);

    for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3) {
    
      for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
        free (pfi->fisher_matrix_lmax[index_l3][index_ft]);
        free (pfi->fisher_matrix_largest[index_l3][index_ft]);
        free (pfi->inverse_fisher_matrix_lmax[index_l3][index_ft]);
      }
    
      free (pfi->fisher_matrix_lmax[index_l3]);
      free (pfi->fisher_matrix_largest[index_l3]);
      free (pfi->inverse_fisher_matrix_lmax[index_l3]);
      free (pfi->sigma_fnl_lmax[index_l3]);
    
    } // end of for(index_l3)
    
    free (pfi->fisher_matrix_lmax);
    free (pfi->fisher_matrix_largest);
    free (pfi->inverse_fisher_matrix_lmax);
    free (pfi->sigma_fnl_lmax);

    for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
      for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
        free (pfi->fisher_matrix_lmin[index_l1][index_ft]);
        free (pfi->fisher_matrix_smallest[index_l1][index_ft]);
        free (pfi->inverse_fisher_matrix_lmin[index_l1][index_ft]);
      }
    
      free (pfi->fisher_matrix_lmin[index_l1]);
      free (pfi->fisher_matrix_smallest[index_l1]);
      free (pfi->inverse_fisher_matrix_lmin[index_l1]);
      free (pfi->sigma_fnl_lmin[index_l1]);
    
    } // end of for(index_l1)

    free (pfi->fisher_matrix_lmin);
    free (pfi->fisher_matrix_smallest);
    free (pfi->inverse_fisher_matrix_lmin);
    free (pfi->sigma_fnl_lmin);    
    
    if (pfi->include_lensing_effects == _TRUE_) {

      int N = pfi->ff_size*pfi->ft_size;
    
      for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {  
        for (int i=0; i < N; ++i) {
          free (pfi->fisher_matrix_CZ_smallest[index_l1][i]);
        }
        free (pfi->fisher_matrix_CZ_smallest[index_l1]);
      }
      free (pfi->fisher_matrix_CZ_smallest);

      for(int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
        for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
          free (pfi->fisher_matrix_lensvar_lmin[index_l1][index_ft]);
          free (pfi->fisher_matrix_lensvar_smallest[index_l1][index_ft]);
          free (pfi->inverse_fisher_matrix_lensvar_lmin[index_l1][index_ft]);
        }
    
        free (pfi->fisher_matrix_lensvar_lmin[index_l1]);
        free (pfi->fisher_matrix_lensvar_smallest[index_l1]);
        free (pfi->inverse_fisher_matrix_lensvar_lmin[index_l1]);
        free (pfi->sigma_fnl_lensvar_lmin[index_l1]);
    
      } // end of for(index_l1)

      free (pfi->fisher_matrix_lensvar_lmin);
      free (pfi->fisher_matrix_lensvar_smallest);
      free (pfi->inverse_fisher_matrix_lensvar_lmin);
      free (pfi->sigma_fnl_lensvar_lmin);
      
      if (pfi->compute_lensing_variance_lmax == _TRUE_) {

        for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {
          for (int index_l3=0; index_l3 < pfi->l3_size; ++index_l3) {
            for (int i=0; i < N; ++i) {
              free (pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][i]);
            }
            free (pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3]);
          }
          free (pfi->fisher_matrix_CZ_smallest_largest[index_l1]);
        } 
        free (pfi->fisher_matrix_CZ_smallest_largest);

        for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
          for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
            free (pfi->fisher_matrix_lensvar_lmax[index_l1][index_ft]);
            free (pfi->inverse_fisher_matrix_lensvar_lmax[index_l1][index_ft]);
          }
    
          free (pfi->fisher_matrix_lensvar_lmax[index_l1]);
          free (pfi->inverse_fisher_matrix_lensvar_lmax[index_l1]);
          free (pfi->sigma_fnl_lensvar_lmax[index_l1]);
    
        } // end of for(index_l1)
    
        free (pfi->fisher_matrix_lensvar_lmax);
        free (pfi->inverse_fisher_matrix_lensvar_lmax);
        free (pfi->sigma_fnl_lensvar_lmax);
      } // end of(compute_lensing_variance_lmax)

    } // end of(include_lensing_effects)

    free (pfi->l1);
    free (pfi->l2);
    free (pfi->l3);
    
    /* Free 3j symbols */
    if ((pbi->interpolation_method != mesh_interpolation_3D) &&
        (pbi->interpolation_method != mesh_interpolation_2D) &&
        (pbi->interpolation_method != bilinear_interpolation)) 
      free (pfi->I_l1_l2_l3);
    
    /* Free cross-power spectrum */
    for (int l=pbi->l_min; l <= pbi->l_max; ++l) {

      for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff) {
        free (pfi->C[l-2][index_ff]);
        free (pfi->inverse_C[l-2][index_ff]);
      }
      free (pfi->C[l-2]);
      free (pfi->inverse_C[l-2]);
    }
    free (pfi->C);
    free (pfi->inverse_C);
      
    /* Free noise spectrum */
    for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff)
      free(pfi->N_l[index_ff]);
    free(pfi->N_l);
    
  } // end of if(has_fisher)
  
  return _SUCCESS_;
 
}




/**
 * Define indices and allocates tables in the Fisher structure
 */

int fisher_indices (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        )
{


  // ====================================================================================
  // =                             Which fields to use?                                 =
  // ====================================================================================

  /* Find out which fields (index_ff=T,E,B...) to include in the Fisher matrix analysis. By default use
  all the fields that were computed in the bispectrum module (denoted by the indices index_bf), unless
  they are explicitly ignored by the pfi->ignore_X flags. */
  
  int index_ff = 0;
  
  pfi->has_fisher_t = _FALSE_;
  pfi->has_fisher_e = _FALSE_;
  pfi->has_fisher_b = _FALSE_;
    
  if ((pbi->has_bispectra_t == _TRUE_) && (pfi->ignore_t == _FALSE_)) {
    pfi->has_fisher_t = _TRUE_;
    pfi->index_bf_of_ff[index_ff] = pbi->index_bf_t;
    pfi->index_ff_t = index_ff++;
  }
  
  if ((pbi->has_bispectra_e == _TRUE_) && (pfi->ignore_e == _FALSE_)) {
    pfi->has_fisher_e = _TRUE_;
    pfi->index_bf_of_ff[index_ff] = pbi->index_bf_e;
    pfi->index_ff_e = index_ff++;
  }
  
  if ((pbi->has_bispectra_b == _TRUE_) && (pfi->ignore_b == _FALSE_)) {
    pfi->has_fisher_b = _TRUE_;
    pfi->index_bf_of_ff[index_ff] = pbi->index_bf_b;
    pfi->index_ff_b = index_ff++;
  }
    
  pfi->ff_size = index_ff;
  pfi->n_probes = pow(pfi->ff_size, 3);

  class_test (pfi->ff_size < 1,
    pfi->error_message,
    "no fields (T,E,B...) requested");


  /* Fill the arrays that relate the fields to quantities computed throughout the code */

  for (int X = 0; X < pfi->ff_size; ++X) {
    
    int index_bf_X = pfi->index_bf_of_ff[X];
    pfi->index_ct_of_phi_ff[X] = pbi->index_ct_of_phi_bf[index_bf_X];
    strcpy (pfi->ff_labels[X], pbi->bf_labels[index_bf_X]);
    
    for (int Y = 0; Y < pfi->ff_size; ++Y) {

      int index_bf_Y = pfi->index_bf_of_ff[Y];
      pfi->index_ct_of_ff_ff[X][Y] = pbi->index_ct_of_bf_bf[index_bf_X][index_bf_Y];

      if (pfi->include_lensing_effects == _TRUE_)
        pfi->index_lt_of_ff_ff[X][Y] = pbi->index_lt_of_bf_bf[index_bf_X][index_bf_Y];
      
      for (int Z = 0; Z < pfi->ff_size; ++Z) {
        
        int index_bf_Z = pfi->index_bf_of_ff[Z];
        strcpy (pfi->ffff_labels[X][Y][Z], pbi->bfff_labels[index_bf_X][index_bf_Y][index_bf_Z]);

      } // Z
    } // Y
  } // X
  

  /* Print information on the Fisher matrix to be computed */

  if (pfi->fisher_verbose > 0) {
    printf_log (" -> Fisher matrix will include %d field%s (",
      pfi->ff_size, (pfi->ff_size!=1)?"s":"");
    for (int index_ff=0; index_ff < (pfi->ff_size-1); ++index_ff)
      printf_log ("%s,", pfi->ff_labels[index_ff]);
    printf_log ("%s)", pfi->ff_labels[pfi->ff_size-1]);
    int n_ignored = pbi->bf_size - pfi->ff_size;
    if (n_ignored > 0) {
      printf_log (" and ignore %d\n", n_ignored);
    }
    else {
      printf_log ("\n");
    }
    
    if (pfi->squeezed_ratio > 1)
      printf_log (" -> only squeezed triangles with L2/L1 > %g\n",
        pfi->squeezed_ratio);

    if (pfi->squeezed_ratio < -1)
      printf_log (" -> only equilateral triangles with L3/L1 < %g\n",
        fabs(pfi->squeezed_ratio));
  }
    
  // ====================================================================================
  // =                              Which bispectra to use?                             =
  // ====================================================================================

  /* We won't need all the bispectra computed in the bispectra.c module in the Fisher
  matrix. Here we select those we are interested into, and create a correspondence
  between rows of the Fisher matrix and bispectra position in pbi->bispectra[index_bt] */

  int index_ft = 0;
  
  if (pbi->has_local_model) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_local;
    pfi->index_ft_local = index_ft++;
  }
  if (pbi->has_equilateral_model) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_equilateral;
    pfi->index_ft_equilateral = index_ft++;
  }  
  if (pbi->has_orthogonal_model) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_orthogonal;
    pfi->index_ft_orthogonal = index_ft++;
  }
  if (pbi->has_galileon_model) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_galileon_gradient;
    pfi->index_ft_galileon_gradient = index_ft++;
    
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_galileon_time;
    pfi->index_ft_galileon_time = index_ft++;
  }
  if (pbi->has_local_squeezed == _TRUE_) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_local_squeezed;
    pfi->index_ft_local_squeezed = index_ft++;
  }
  if (pbi->has_intrinsic_squeezed == _TRUE_) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_intrinsic_squeezed;
    pfi->index_ft_intrinsic_squeezed = index_ft++;
  }
  if (pbi->has_cosine_shape == _TRUE_) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_cosine;
    pfi->index_ft_cosine = index_ft++;
  }
  /* If asking for lensing effects, the CMB-lensing bispectrum is automatically included
  in the Fisher matrix. It is important to avoid duplicates as otherwise the Fisher
  matrix would be singular, and therefore impossible to be inverted to compute the effect
  of lensing variance. */
  if (pfi->include_lensing_effects == _FALSE_) {
    if (pbi->has_cmb_lensing == _TRUE_) {
      pfi->index_bt_of_ft[index_ft] = pbi->index_bt_cmb_lensing;
      pfi->index_ft_cmb_lensing = index_ft++;
    }
    if (pbi->has_cmb_lensing_squeezed == _TRUE_) {
      pfi->index_bt_of_ft[index_ft] = pbi->index_bt_cmb_lensing_squeezed;
      pfi->index_ft_cmb_lensing_squeezed = index_ft++;
    }
  }
  if (pbi->has_cmb_lensing_kernel == _TRUE_) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_cmb_lensing_kernel;
    pfi->index_ft_cmb_lensing_kernel = index_ft++;
  }
  if (pbi->has_intrinsic == _TRUE_) {
    pfi->index_bt_of_ft[index_ft] = pbi->index_bt_intrinsic;
    pfi->index_ft_intrinsic = index_ft++;
  }
    
  pfi->ft_size = index_ft;  

  /* Fill the arrays that relate the fields to quantities computed throughout the code */
  for (int index_ft = 0; index_ft < pfi->ft_size; ++index_ft) {
    int index_bt = pfi->index_bt_of_ft[index_ft];
    strcpy (pfi->ft_labels[index_ft], pbi->bt_labels[index_bt]);
  }

  /* Print information on the Fisher matrix to be computed */
  if (pfi->fisher_verbose > 0) {
    printf_log (" -> Fisher matrix will have %d row%s: ",
      pfi->ft_size, (pfi->ft_size!=1)?"s":"");
    for (int index_ft=0; index_ft < (pfi->ft_size-1); ++index_ft)
      printf_log ("%s, ", pfi->ft_labels[index_ft]);
    printf_log ("%s\n", pfi->ft_labels[pfi->ft_size-1]);
  }


  // ====================================================================================
  // =                                Choose l-grid                                     =
  // ====================================================================================

  /* The Fisher matrix is obtained as a sum over (l1,l2,l3). Having computed the
  bispectrum on a grid in (l1,l2,l3), we need to resort to some interpolation or
  integration scheme in order to evaluate the sum. Depending on the adopted
  method, we consider different grids for the three l-directions, and store
  them in pfi->l1, pfi->l2 and pfi->l3. */


  /* In the simplest approach, the trilinear integration, we evaluate the Fisher
  matrix as a weighted sum over the nodes, without thus needing to interpolate
  the bispectrum. In this case, we take the three l samplings to coincide with
  the nodes sampling, pbi->l. */

  if ((pbi->interpolation_method == trilinear_interpolation) ||
      (pbi->interpolation_method == smart_interpolation)) {
    
    pfi->l1_is_full = _FALSE_;
    pfi->l2_is_full = _FALSE_;
    pfi->l3_is_full = _FALSE_;
    
  }

  /* In the bilinear interpolation approach, we cycle over all the (l2,l3) pairs,
  from l_min=2 to l_max, while for l1 we only consider those multipoles in pbi->l.
  The bispectrum is obtained in (l2,l3) via bilinear interpolation, while the sum
  over l1 is weighted with trapezoidal weights. The mesh_2D approach differs from
  the bilinear approach only in how the interpolation is performed; the sums are
  identical. In particular, l1 is drawn from pbi->l, while l2 and l3 go from
  l_min=2 to l_max. */

  else if ((pbi->interpolation_method == bilinear_interpolation) ||
           (pbi->interpolation_method == mesh_interpolation_2D)) {
  
    pfi->l1_is_full = _FALSE_;
    pfi->l2_is_full = _TRUE_;
    pfi->l3_is_full = _TRUE_;

  }
  
  /* The mesh_3D scheme is the only scheme that allows us to interpolate
  the bispectrum in any (l1,l2,l3) configuration, and not only in (l2,l3).
  We use this property to choose a full grid in (l1,l2,l3). This approach
  is very time consuming but it gives a nice countercheck to the other
  methods. */

  else if (pbi->interpolation_method == mesh_interpolation_3D) {
    
    pfi->l1_is_full = _FALSE_; /* set to _FALSE_ to use bilinear interpolation */
    pfi->l2_is_full = _TRUE_;
    pfi->l3_is_full = _TRUE_;
    
  }

  /* Build the actual grids */

  pfi->l1_size = pfi->l1_is_full ? pbi->full_l_size : pbi->l_size;
  pfi->l2_size = pfi->l2_is_full ? pbi->full_l_size : pbi->l_size;
  pfi->l3_size = pfi->l3_is_full ? pbi->full_l_size : pbi->l_size;

  class_alloc (pfi->l1, pfi->l1_size*sizeof(int), pfi->error_message);
  class_alloc (pfi->l2, pfi->l2_size*sizeof(int), pfi->error_message);
  class_alloc (pfi->l3, pfi->l3_size*sizeof(int), pfi->error_message);

  for (int index_l=0; index_l < pfi->l1_size; ++index_l)
    pfi->l1[index_l] = pfi->l1_is_full ? (pbi->l_min + index_l) : pbi->l[index_l];
  for (int index_l=0; index_l < pfi->l2_size; ++index_l)
    pfi->l2[index_l] = pfi->l2_is_full ? (pbi->l_min + index_l) : pbi->l[index_l];
  for (int index_l=0; index_l < pfi->l3_size; ++index_l)
    pfi->l3[index_l] = pfi->l3_is_full ? (pbi->l_min + index_l) : pbi->l[index_l];


  /* Which l-configurations should be included in the forecast? We set global limits
  on l based on the input form the user (pfi->l_min_estimator and pfi->l_max_estimator)
  and on what has been computed so far (pbi->l_min, pbi->l_max). Modify by hand if you
  want to include specific configurations only, eg. only squeezed ones. Remember that
  in the Fisher sum we consider only l1<=l2<=l3. */

  pfi->l1_min_global = MAX (pfi->l_min_estimator, pbi->l_min);
  pfi->l2_min_global = MAX (pfi->l_min_estimator, pbi->l_min);
  pfi->l3_min_global = MAX (pfi->l_min_estimator, pbi->l_min); 

  pfi->l1_max_global = MIN (pfi->l_max_estimator, pbi->l_max);
  pfi->l2_max_global = MIN (pfi->l_max_estimator, pbi->l_max);
  pfi->l3_max_global = MIN (pfi->l_max_estimator, pbi->l_max);



  // ====================================================================================
  // =                            Allocate Fisher matrix                                =
  // ====================================================================================
  
  /* Allocate the arrays where we shall store the various contributions to the Fisher
  matrix. It is crucial that these arrays are initialised with calloc, as they will
  be incremented in the Fisher matrix estimator function.  */
  
  class_alloc (pfi->fisher_matrix_XYZ_largest, pfi->ff_size*sizeof(double *****), pfi->error_message);
  class_alloc (pfi->fisher_matrix_XYZ_lmax, pfi->ff_size*sizeof(double *****), pfi->error_message);
  class_alloc (pfi->fisher_matrix_XYZ_smallest, pfi->ff_size*sizeof(double *****), pfi->error_message);
  class_alloc (pfi->fisher_matrix_XYZ_lmin, pfi->ff_size*sizeof(double *****), pfi->error_message);

  for (int X = 0; X < pfi->ff_size; ++X) {

    class_alloc (pfi->fisher_matrix_XYZ_largest[X], pfi->ff_size*sizeof(double ****), pfi->error_message);
    class_alloc (pfi->fisher_matrix_XYZ_lmax[X], pfi->ff_size*sizeof(double ****), pfi->error_message);
    class_alloc (pfi->fisher_matrix_XYZ_smallest[X], pfi->ff_size*sizeof(double ****), pfi->error_message);
    class_alloc (pfi->fisher_matrix_XYZ_lmin[X], pfi->ff_size*sizeof(double ****), pfi->error_message);
        
    for (int Y = 0; Y < pfi->ff_size; ++Y) {

      class_alloc (pfi->fisher_matrix_XYZ_largest[X][Y], pfi->ff_size*sizeof(double ***), pfi->error_message);
      class_alloc (pfi->fisher_matrix_XYZ_lmax[X][Y], pfi->ff_size*sizeof(double ***), pfi->error_message);
      class_alloc (pfi->fisher_matrix_XYZ_smallest[X][Y], pfi->ff_size*sizeof(double ***), pfi->error_message);
      class_alloc (pfi->fisher_matrix_XYZ_lmin[X][Y], pfi->ff_size*sizeof(double ***), pfi->error_message);

      for (int Z = 0; Z < pfi->ff_size; ++Z) {

        /* Arrays as a function of lmax (involves l3) */
        class_alloc (pfi->fisher_matrix_XYZ_largest[X][Y][Z], pfi->l3_size*sizeof(double **), pfi->error_message);
        class_alloc (pfi->fisher_matrix_XYZ_lmax[X][Y][Z], pfi->l3_size*sizeof(double **), pfi->error_message);

        for (int index_l3=0; index_l3 < pfi->l3_size; ++index_l3) {
          
          class_alloc (pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3], pfi->ft_size*sizeof(double *), pfi->error_message);
          class_alloc (pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l3], pfi->ft_size*sizeof(double *), pfi->error_message);
          
          for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
            class_calloc (pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft],
              pfi->ft_size, sizeof(double), pfi->error_message);
            class_calloc (pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l3][index_ft],
              pfi->ft_size, sizeof(double), pfi->error_message);            
          }
        }

        /* Arrays as a function of lmin (involves l1) */
        class_alloc (pfi->fisher_matrix_XYZ_smallest[X][Y][Z], pfi->l1_size*sizeof(double **), pfi->error_message);
        class_alloc (pfi->fisher_matrix_XYZ_lmin[X][Y][Z], pfi->l1_size*sizeof(double **), pfi->error_message);

        for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {
          
          class_alloc (pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
          class_alloc (pfi->fisher_matrix_XYZ_lmin[X][Y][Z][index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
          
          for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
            class_calloc (pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft],
              pfi->ft_size, sizeof(double), pfi->error_message);
            class_calloc (pfi->fisher_matrix_XYZ_lmin[X][Y][Z][index_l1][index_ft],
              pfi->ft_size, sizeof(double), pfi->error_message);
          }
        } 
        
      } // Z
    } // Y 
  } // X
  
  /* More arrays as a function of lmax (involving l3) */
  class_alloc (pfi->fisher_matrix_lmax, pfi->l3_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->fisher_matrix_largest, pfi->l3_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->inverse_fisher_matrix_lmax, pfi->l3_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->sigma_fnl_lmax, pfi->l3_size*sizeof(double *), pfi->error_message);
  
  for (int index_l3=0; index_l3 < pfi->l3_size; ++index_l3) {

    class_alloc (pfi->fisher_matrix_lmax[index_l3], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->fisher_matrix_largest[index_l3], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->inverse_fisher_matrix_lmax[index_l3], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_calloc (pfi->sigma_fnl_lmax[index_l3], pfi->ft_size, sizeof(double), pfi->error_message);
    
    for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
      class_calloc (pfi->fisher_matrix_lmax[index_l3][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->fisher_matrix_largest[index_l3][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->inverse_fisher_matrix_lmax[index_l3][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
    }    
  }

  /* More arrays as a function of lmin (involving l1) */
  class_alloc (pfi->fisher_matrix_lmin, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->fisher_matrix_smallest, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->inverse_fisher_matrix_lmin, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->sigma_fnl_lmin, pfi->l1_size*sizeof(double *), pfi->error_message);
  
  for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {

    class_alloc (pfi->fisher_matrix_lmin[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->fisher_matrix_smallest[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->inverse_fisher_matrix_lmin[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
    class_calloc (pfi->sigma_fnl_lmin[index_l1], pfi->ft_size, sizeof(double), pfi->error_message);
    
    for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
      class_calloc (pfi->fisher_matrix_lmin[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->fisher_matrix_smallest[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->inverse_fisher_matrix_lmin[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
    }    
  }

  if (pfi->include_lensing_effects == _TRUE_) {
    
    int N = pfi->ff_size*pfi->ft_size;

    /* Allocate fisher_matrix_CZ_smallest, needed to compute the lensing variance. It is equivalent to \bar{F}_{l_1}^{ip}
    in Eq. 5.25 of Lewis et al. 2011 (see header file). */
    class_alloc (pfi->fisher_matrix_CZ_smallest, pfi->l1_size*sizeof(double **), pfi->error_message);
    for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {
      class_alloc (pfi->fisher_matrix_CZ_smallest[index_l1], N*sizeof(double *), pfi->error_message);
      for (int i=0; i < N; ++i)
        class_calloc (pfi->fisher_matrix_CZ_smallest[index_l1][i], N, sizeof(double), pfi->error_message);
    }
    
    /* Arrays that will contain the Fisher matrix including lensing variance */
    class_alloc (pfi->fisher_matrix_lensvar_lmin, pfi->l1_size*sizeof(double **), pfi->error_message);
    class_alloc (pfi->fisher_matrix_lensvar_smallest, pfi->l1_size*sizeof(double **), pfi->error_message);
    class_alloc (pfi->inverse_fisher_matrix_lensvar_lmin, pfi->l1_size*sizeof(double **), pfi->error_message);
    class_alloc (pfi->sigma_fnl_lensvar_lmin, pfi->l1_size*sizeof(double *), pfi->error_message);
  
    for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {

      class_alloc (pfi->fisher_matrix_lensvar_lmin[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
      class_alloc (pfi->fisher_matrix_lensvar_smallest[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
      class_alloc (pfi->inverse_fisher_matrix_lensvar_lmin[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
      class_calloc (pfi->sigma_fnl_lensvar_lmin[index_l1], pfi->ft_size, sizeof(double), pfi->error_message);
    
      for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
        class_calloc (pfi->fisher_matrix_lensvar_lmin[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
        class_calloc (pfi->fisher_matrix_lensvar_smallest[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
        class_calloc (pfi->inverse_fisher_matrix_lensvar_lmin[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
      }    
    }
    
    /* Same as above, but for l_max. Note that we allocate these arrays with pfi->l1_size
    rather than pfi->l3_size, even though so far we associated l3 to the largest multipole
    and l1 to the smallest one. The reason is that, for mesh interpolation, l3 has a fine
    grid (e.g. pfi->l3_size=1999 for for l_max=2000), while l1 is always sampled (e.g.
    pfi->l1_size=O(100) for l_max=2000). Therefore, using l1 intead of l3 allows us to
    compute the lensing variance correction a much smaller number of times */
    if (pfi->compute_lensing_variance_lmax == _TRUE_) {

      class_alloc (pfi->fisher_matrix_lensvar_lmax, pfi->l1_size*sizeof(double **), pfi->error_message);
      class_alloc (pfi->inverse_fisher_matrix_lensvar_lmax, pfi->l1_size*sizeof(double **), pfi->error_message);
      class_alloc (pfi->sigma_fnl_lensvar_lmax, pfi->l1_size*sizeof(double *), pfi->error_message);
  
      for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {

        class_alloc (pfi->fisher_matrix_lensvar_lmax[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
        class_alloc (pfi->inverse_fisher_matrix_lensvar_lmax[index_l1], pfi->ft_size*sizeof(double *), pfi->error_message);
        class_calloc (pfi->sigma_fnl_lensvar_lmax[index_l1], pfi->ft_size, sizeof(double), pfi->error_message);
    
        for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
          class_calloc (pfi->fisher_matrix_lensvar_lmax[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
          class_calloc (pfi->inverse_fisher_matrix_lensvar_lmax[index_l1][index_ft], pfi->ft_size, sizeof(double), pfi->error_message);
        }    
      }
      
      /* Allocate fisher_matrix_CZ_smallest_largest, needed to compute the lensing variance as a function of l_max. */
      class_alloc (pfi->fisher_matrix_CZ_smallest_largest, pfi->l1_size*sizeof(double ***), pfi->error_message);
      for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {
        class_alloc (pfi->fisher_matrix_CZ_smallest_largest[index_l1], pfi->l3_size*sizeof(double **), pfi->error_message);
        for (int index_l3=0; index_l3 < pfi->l3_size; ++index_l3) {
          class_alloc (pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3], N*sizeof(double *), pfi->error_message);
          for (int i=0; i < N; ++i)
            class_calloc (pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][i], N, sizeof(double), pfi->error_message);
        }
      }
      
    } // if(compute_lensing_variance_lmax)
    
  } // if(include_lensing_effects)


  /* Allocate array for the noise in the C_l */
  class_alloc (pfi->N_l, pfi->ff_size*sizeof(double *), pfi->error_message);
  for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff)
    class_calloc (pfi->N_l[index_ff], pbi->full_l_size, sizeof(double), pfi->error_message);



  // ====================================================================================
  // =                                   Compute 3j                                     =
  // ====================================================================================
  
  /* Allocate and fill the 3j array containing I_l1_l2_l3, defined in eq. 13 of
  Komatsu, Spergel & Wandelt (2005):
  
     I_l1_l2_l3 = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/4*pi) *  3J(l1 l2 l3)(0  0  0)
  
  which is needed to compute the angle-averaged bispectra, which in turn are
  used by the f_NL estimator. */
  
  /* For mesh interpolation, we shall compute the 3j's directly in the estimator */
  if ((pbi->interpolation_method != mesh_interpolation_3D) &&
      (pbi->interpolation_method != mesh_interpolation_2D) &&
      (pbi->interpolation_method != bilinear_interpolation)) {
  
    printf_log_if (pfi->fisher_verbose, 1, 
      " -> allocating and computing ~ %.3g MB (%ld doubles) of 3j-symbols\n",
      pbi->n_independent_configurations*sizeof(double)/1e6, pbi->n_independent_configurations);
  
    /* Allocate memory for pfi->I_l1_l2_l3 */
    class_alloc (pfi->I_l1_l2_l3, pbi->n_independent_configurations*sizeof(double), pfi->error_message);
  
    /* Temporary values needed for the computation of the 3j symbol. The temporary array will contain at most 2*l_max+1
      values when l1=l2=l_max, with the +1 accounting for the l=0 case. */
    double * temporary_three_j;
    class_alloc (temporary_three_j, (2*pbi->l_max+1)*sizeof(double), pfi->error_message);
    double l3_min_D, l3_max_D;
             
    /* Fill the I_l1_l2_l3 array */
    /* TODO:   parallelize the loop by creating a 3j array per thread */
    // #pragma omp parallel \
    // private (index_l1, index_l2, index_l3, l3_min, l3_max, temporary_three_j_size, index_l1_l2_l3)
    
    for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {
      
      int l1 = pbi->l[index_l1];
    
      for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
        
        int l2 = pbi->l[index_l2];
        
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (l2 < l1/2)
          continue;
        
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);

        /* Compute the actual 3j symbol for all allowed values of l3 */
        class_call (drc3jj (
                      l1, l2, 0, 0,
                      &l3_min_D, &l3_max_D,
                      temporary_three_j,
                      (2*pbi->l_max+1),
                      pfi->error_message       
                      ),
          pfi->error_message,
          pfi->error_message);

        int l3_min = (int)(l3_min_D+_EPS_);
    
        /* Fill the 3j array */
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
                
          int l3 = pbi->l[index_l3];
          
          /* Index of the current (l1,l2,l3) configuration */
          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
    
          /* This 3j symbol must vanish when l1+l2+l3 is even. */
          class_test (((l1+l2+l3)%2!=0) && (temporary_three_j[l3-l3_min] != 0.),
            pfi->error_message,
            "error in computation of 3j, odd is not zero!");
        
          /* The following is what it takes to convert a reduced bispectrum to an angle-averaged one */
          pfi->I_l1_l2_l3[index_l1_l2_l3] =
            sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.)/(4.*_PI_)) * temporary_three_j[l3-l3_min];
  
          /* Some debug */
          // fprintf (stderr, "%5d %5d %5d %5d/%5d %17.7g %17.7g %17.7g\n",
          //   l1, l2, l3, index_l1_l2_l3, pbi->n_independent_configurations-1,
          //     pfi->I_l1_l2_l3[index_l1_l2_l3], sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.)/(4.*_PI_)), temporary_three_j[l3-l3_min]);
  
        } // for(index_l3)
      } // for(index_l2)
    } // for(index_l1)
  
    free (temporary_three_j);
  
  } // if not mesh interpolation


  /* Debug: check the 3j routine */
  // int l1=1200;
  // int l2=1500;
  // 
  // double * temporary_three_j;
  // class_alloc (temporary_three_j, (l1+l2+1)*sizeof(double), pfi->error_message);
  // double l3_min_D, l3_max_D, temporary_three_j_size;
  // 
  // class_call (drc3jj (
  //               l1, l2, 0, 0,
  //               &l3_min_D, &l3_max_D,
  //               temporary_three_j,
  //               (2*pbi->l_max+1),
  //               pfi->error_message       
  //               ),
  //   pfi->error_message,
  //   pfi->error_message);
  //   
  // int l3_min = (int)(l3_min_D+_EPS_);
  // int l3_max = (int)(l3_max_D+_EPS_);
  //  
  // for (int l3=l3_min; l3<=l3_max; ++l3)
  //   printf ("I(%d,%d,%d)=%.15g\n", l1,l2,l3,temporary_three_j[l3-l3_min]);



  // ====================================================================================
  // =                          Prepare interpolation mesh                              =
  // ====================================================================================

  /* Compute the 3D interpolation mesh if needed. Given a bispectrum b_l1_l2_l3,
  the mesh is basically its interpolation table conveniently binned in l space.
  Note that we cannot compute the 3D mesh at the last moment, like we do for the
  2D mesh, because it would create a race condition in fisher_compute_matrix(). 
  See comment in bispectra_mesh_interpolate() for more details. */

  if (pbi->interpolation_method == mesh_interpolation_3D) {

    for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {

      if (pbi->interpolate_me[pfi->index_bt_of_ft[index_ft]]) {

        /* The array pbi->meshes contains the table required to interpolate the bispectra
        with the 2D and 3D mesh methods. Its first level, pbi->meshes[index_l1], is used
        only for the 2D mesh interpolation. For 3D interpolation, this level is redundant
        because a single 3D mesh is enough to interpolate the whole (l1,l2,l3) space.
        Therefore, for 3D interpolation we have only one entry: pbi->meshes[0]. */
        class_call (bispectra_mesh_create_for_all_probes(
                      ppr, psp, ple, pbi,
                      pfi->index_bt_of_ft[index_ft],
                      -1,
                      pbi->meshes[0]),
          pbi->error_message,
          pfi->error_message);

      }
    }
  }


  // ====================================================================================
  // =                              Compute wasted points                               =
  // ====================================================================================

  /* Count the number of bispectra configurations that won't be used directly in
  the Fisher sum, ie. those with l1+l2+l3 odd. These points will be still crucial
  to interpolate the bispectrum in the configurations with l1+l2+l3 even. */

  long int count_wasted = 0;

  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {
    int l1 = pbi->l[index_l1];
    for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
      int l2 = pbi->l[index_l2];
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
        if ((l1+l2+pbi->l[index_l3])%2!=0)
          count_wasted++;
    }
  }

  printf_log_if (pfi->fisher_verbose, 1, 
    "     * %.3g%% of the computed bispectra configurations will be (directly) used to compute the Fisher matrix\n",
    100-count_wasted/(double)(pbi->n_independent_configurations)*100);


  /* Initialise info strings */
  sprintf (pfi->info, "");
  if (pfi->include_lensing_effects == _TRUE_)
    sprintf (pfi->info_lensvar, "");
  

  return _SUCCESS_;

}



/**
 * Fill the pfi->N_l arrays with the observational noise for the C_l's.
 *
 * For the time being, we have implemented the noise model in astro-ph/0506396v2, a co-added
 * Gaussian noise where the beam for each frequency channel of the experiment is considered
 * to be Gaussian, and the noise homogeneous.
 *
 */
int fisher_noise (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        )
{
  
  // ===================================================================================
  // =                       Read in the noise values for the fields                   =
  // ===================================================================================

  double noise[pfi->ff_size][pfi->n_channels];
  
  for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff) {
    for (int index_channel=0; index_channel < pfi->n_channels; ++index_channel) {

      if ((pfi->has_fisher_t == _TRUE_) && (index_ff == pfi->index_ff_t))
        noise[index_ff][index_channel] = pfi->noise_t[index_channel];

      if ((pfi->has_fisher_e == _TRUE_) && (index_ff == pfi->index_ff_e))
        noise[index_ff][index_channel] = pfi->noise_e[index_channel];

    }
  }


  // ===================================================================================
  // =                                Print information                                =
  // ===================================================================================

  if (pfi->fisher_verbose > 0) {

    printf_log (" -> noise model: beam_fwhm=");
    for (int index_channel=0; index_channel < pfi->n_channels-1; ++index_channel)
      printf ("%g,", pfi->beam[index_channel] * 60 * 180./_PI_);
    printf_log ("%g", pfi->beam[pfi->n_channels-1] * 60 * 180./_PI_);
    printf_log (" (arcmin)");
    
    for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff) {

      printf_log (", sigma_%s=", pfi->ff_labels[index_ff]);
      for (int index_channel=0; index_channel < pfi->n_channels-1; ++index_channel)
        printf_log ("%g,", sqrt(noise[index_ff][index_channel]) * 1e6*pba->T_cmb / pfi->beam[index_channel]);
      printf_log ("%g", sqrt(noise[index_ff][pfi->n_channels-1]) * 1e6*pba->T_cmb / pfi->beam[pfi->n_channels-1]);
      printf_log (" (uK)");
    }
    printf_log ("\n");
  }



  // ===================================================================================
  // =                             Co-added Gaussian noise                             =
  // ===================================================================================
  
  /* Eeach frequency channel contributes to the inverse noise with a 1/N contribution */
  for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff) {
    for (int l=pbi->l_min; l <= pbi->l_max; ++l) {
      for (int index_channel=0; index_channel < pfi->n_channels; ++index_channel) {

        double beam_exponent = (l*(l+1.)*pfi->beam[index_channel]*pfi->beam[index_channel])/(8.*log(2.));

        class_test (beam_exponent > DBL_MAX_EXP,
          pfi->error_message,
          "stop to prevent overflow. Reduce experiment_beam_fwhm to something more reasonable, or set it to zero.")
        
        /* Inverse contribution to l from the considered channel */
        pfi->N_l[index_ff][l-2] += 1. / (noise[index_ff][index_channel] * exp(beam_exponent));
      }

      /* Invert back */
      pfi->N_l[index_ff][l-2] = 1. / pfi->N_l[index_ff][l-2];

      /* Is this really needed? If the C_l's are infinite, it means that the signal is
      swamped by the noise, which is not an error */
      class_test (fabs(pfi->N_l[index_ff][l-2]) > _HUGE_,
        pfi->error_message,
        "stopping to avoid infinities, reduce noise parameters: pfi->N_%d(%s) = %g",
        l, pfi->ff_labels[index_ff], pfi->N_l[index_ff][l-2]);
    
      /* Some debug */
      // double factor = 1e12*pba->T_cmb*pba->T_cmb*l*(l+1)/(2*_PI_);
      // fprintf (stderr, "%d %g %g\n", l, pfi->N_l[pfi->index_ff_t][l-2], factor*pfi->N_l[pfi->index_ff_t][l-2]);
    
    } // end of for(l=2)
  } // end of for(T,E,...)
  

  return _SUCCESS_;
  
}



/**
 * Compute the Fisher matrix for all pairs of bispectra.
 *
 * In detail, this function does:
 *
 * -# Initialise a workspace with the interpolation weights.
 *
 * -# Compute the Fisher matrix via fisher_compute_matrix().
 *
 * -# Add the effect of lensing variance to the Fisher matrix,
 *    via fisher_lensing_variance().
 * 
 * Called by bispectra_init().
 */

int fisher_compute (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        )
{
 

  // ====================================================================================
  // =                            Initialise working space                              =
  // ====================================================================================

  /* Allocate workspace needed to estimate f_NL */
  struct fisher_workspace * pwf;
  class_alloc (pwf, sizeof(struct fisher_workspace), pfi->error_message);


  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #ifdef _OPENMP
  #pragma omp parallel private (thread)
  number_of_threads = omp_get_num_threads();
  #endif
    
  /* Interpolation weights for the estimator in the rectangular directions
  (l1 and l2). These are the usual trapezoidal integration weights for
  in discrete space */
  class_alloc (pwf->delta_l, pbi->l_size*sizeof(double), pfi->error_message);

  pwf->delta_l[0] = (pbi->l[1] - pbi->l[0] + 1.)/2.;
  
  for (int index_l=1; index_l < pbi->l_size-1; ++index_l)
    pwf->delta_l[index_l] = (pbi->l[index_l+1] - pbi->l[index_l-1])/2.;
      
  pwf->delta_l[pbi->l_size-1] = (pbi->l[pbi->l_size-1] - pbi->l[pbi->l_size-2] + 1.)/2.;


  /* Allocate memory for the interpolation weights in the triangular direction
  (l3). Being dependent on (l1,l2), they will be computed later inside a loop
  on (l1,l2). */

  if (pbi->interpolation_method == smart_interpolation) {

    class_alloc (pwf->delta_l3, number_of_threads*sizeof(double*), pfi->error_message);

    #pragma omp parallel private (thread)
    {

      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
    
      class_calloc_parallel (pwf->delta_l3[thread], pbi->l_size, sizeof(double), pfi->error_message);

    } if (abort == _TRUE_) return _FAILURE_; // end of parallel region
  
  }
  
  
  if (pfi->fisher_verbose > 0)  {

    char buffer[128];

    if (pbi->interpolation_method == mesh_interpolation_2D)
      strcpy (buffer, "mesh (2D)");
    else if (pbi->interpolation_method == mesh_interpolation_3D)
      strcpy (buffer, "mesh (3D)");
    else if (pbi->interpolation_method == bilinear_interpolation)
      strcpy (buffer, "bilinear");
    else if (pbi->interpolation_method == trilinear_interpolation)
      strcpy (buffer, "trilinear");
    else if (pbi->interpolation_method == smart_interpolation)
      strcpy (buffer, "smart");

    printf_log (" -> computing Fisher matrix for l_max = %d with %s interpolation\n",
      MIN (pfi->l_max_estimator, pbi->l_max), buffer);
    
  }



  // ====================================================================================
  // =                              Compute Fisher matrix                               =
  // ====================================================================================
    

  /* Fill the "pfi->fisher_matrix_XYZ_largest" and "fisher_matrix_XYZ_smallest" Fisher matrices
  using mesh interpolation */
  if ((pbi->interpolation_method == mesh_interpolation_2D) ||
      (pbi->interpolation_method == mesh_interpolation_3D) ||
      (pbi->interpolation_method == bilinear_interpolation)) {
  
    class_call (fisher_compute_matrix(
                  ppr,
                  ptr,
                  psp,
                  ple,
                  pbi,
                  pfi,
                  pwf),
      pfi->error_message,
      pfi->error_message);
      
  }
  
  /* Fill the "pfi->fisher_matrix_XYZ_largest" and "fisher_matrix_XYZ_smallest" Fisher matrices
  using linear interpolation */
  else if ((pbi->interpolation_method == smart_interpolation) ||
           (pbi->interpolation_method == trilinear_interpolation)) {
            
    class_call (fisher_compute_matrix_nodes(
                  ppr,
                  psp,
                  ple,
                  pbi,
                  pfi,
                  pwf),
      pfi->error_message,
      pfi->error_message);
    
  } 
  
  else {

    class_stop (pfi->error_message, "could not understand bispectrum interpolation method!");

  }
  
  
  
  // ====================================================================================
  // =                            Compute lensing variance                              =
  // ====================================================================================

  /* Include the extra variance induced by lensing in the Fisher matrix
  estimator, as explained in Section 5 of Lewis, Challinor & Hanson
  (2011, http://uk.arxiv.org/abs/1101.2234). See the documentation of
  fisher_lensing_variance() for further details */ 

  if (pfi->include_lensing_effects == _TRUE_) {

    printf_log_if (pfi->fisher_verbose, 0, 
      " -> adding lensing-induced noise according to arxiv:1101.2234\n");


    // -------------------------------------------------------------------------------
    // -                          As a function of l_min                             -
    // -------------------------------------------------------------------------------

    if (pfi->compute_lensing_variance_lmax == _FALSE_) {
      
      class_call (fisher_lensing_variance (
                    ppr, pba, ppt, pbs, ptr, ppm,
                    psp, ple, pbi, pfi, pwf),
        pfi->error_message,
        pfi->error_message);
       
    }


    // -------------------------------------------------------------------------------
    // -                           As a function of l_max                            -
    // -------------------------------------------------------------------------------

    /* Compute lensing variance as a function of l_max, the maximum resolution
    of the considered survey. In absence of lensing variance, this can be easily
    extracted from the Fisher matrix summation in fisher_compute_matrix(). To
    include lensing variance, however, one needs to be compute the effect separately
    for each value of l_max. The reason is that the lensing variance formula we adopt
    (Sec. 5 of http://uk.arxiv.org/abs/1101.2234) includes a sum where l1 is the
    smallest multipole in the Fisher summation, while SONG computes the Fisher matrix
    assuming that l1 is the largest multipole in the Fisher summation. */
       
    else {

      for (int index_lmax=0; index_lmax<pfi->l1_size; ++index_lmax) {

        int l_max = pfi->l1[index_lmax];
          
        /* Uncomment to stop at a certain l_max */
        // if (index_lmax < (pfi->l1_size-1))
        //   continue;
        
        printf_log_if (pfi->fisher_verbose, 1, 
          " -> adding lensing noise for l_max=%d\n", l_max);

        /* The function fisher_lensing_variance() computes the lensing variance correction
        based on the content of pfi->fisher_matrix_CZ_smallest, which represents F_bar as
        defined in Section 5 of http://uk.arxiv.org/abs/1101.2234. This array was filled
        above in fisher_compute_matrix() function and is a sum of all the l3 contributions,
        l3 being the largest multipole in the Fisher matrix sum. To break down the effect
        of lensing variance as a function of l_max, the maximum angular resolution of a
        CMB survey, we loop over l_max and overwrite pfi->fisher_matrix_CZ_smallest with a
        sum over l3 which is truncated at l3=l_max. That is, we build F_bar(l3<=l_max). */
        
        /* Initialise pfi->fisher_matrix_CZ_smallest */
        for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1)
          for (int i=0; i < pfi->ff_size*pfi->ft_size; ++i)
            for (int j=0; j < pfi->ff_size*pfi->ft_size; ++j)
              pfi->fisher_matrix_CZ_smallest[index_l1][i][j] = 0;

        /* Overwrite it with F_bar(l3<=l_max) */
        for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1)
          for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3)
            for (int i=0; i < pfi->ff_size*pfi->ft_size; ++i)
              for (int j=0; j < pfi->ff_size*pfi->ft_size; ++j)
                if (pfi->l3[index_l3] <= l_max)
                  pfi->fisher_matrix_CZ_smallest[index_l1][i][j]
                    += pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][i][j];

        /* Compute lensing variance for this value of l_max */
        class_call (fisher_lensing_variance (
                      ppr, pba, ppt, pbs, ptr, ppm,
                      psp, ple, pbi, pfi, pwf),
          pfi->error_message,
          pfi->error_message);

        /* Debug - print the Fisher matrix for the considered value of l_max */
        // printf (" -> lmax=%d\n", l_max);
        // PrintMatrix(pfi->fisher_matrix_lensvar_lmax[index_lmax], pfi->ft_size);
        
        /* Fill the l_max, inverse_lmax and sigma_fnl arrays */
        for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
          for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {

            pfi->fisher_matrix_lensvar_lmax[index_lmax][index_ft_1][index_ft_2]
              = pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_2];

            if (index_ft_1 == index_ft_2)
              pfi->sigma_fnl_lensvar_lmax[index_lmax][index_ft_1]
                = 1/sqrt(pfi->fisher_matrix_lensvar_lmax[index_lmax][index_ft_1][index_ft_1]);
          }
        }

        InverseMatrix (pfi->fisher_matrix_lensvar_lmax[index_lmax], pfi->ft_size,
          pfi->inverse_fisher_matrix_lensvar_lmax[index_lmax]);
                
      } // end of for(index_lmax)

    } // end of if(compute_lensing_variance_lmax)
         
  } // end of lensing variance
  


  // ====================================================================================
  // =                              Compute S/N and f_nl                                =
  // ====================================================================================
  
  /* So far we have computed the Fisher matrix for each value of l1; now we sum
  those values to obtain the Fisher matrix for each value of l_max */
  for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
    for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {

      /* We need to build the cumulative Fisher distribution, so we need
      accumulators */
      double accumulator;
      double accumulator_XYZ[pfi->ff_size][pfi->ff_size][pfi->ff_size];

      /* Initialize the accumulators */
      accumulator = 0;
      for (int X = 0; X < pfi->ff_size; ++X)
        for (int Y = 0; Y < pfi->ff_size; ++Y)
          for (int Z = 0; Z < pfi->ff_size; ++Z)
            accumulator_XYZ[X][Y][Z] = 0.;


      // -------------------------------------------------------------------------------
      // -                            As a function of lmin                            -
      // -------------------------------------------------------------------------------
      
      /* Build the Fisher matrix as a function of the minimum l in the survey:
      sum_{lmin<=l1<=lmax} fisher_matrix_XYZ_smallest(l1) */

      for (int index_l1=(pfi->l1_size-1); index_l1>=0; --index_l1) {

        /* When including lensing variance, the squeezed kernel is used in the Fisher
        matrix. This has a C_l^{Z\phi} factored out and therefore is much larger than
        any other bispectrum. In fisher_lensing_variance() we have rescaled
        fisher_matrix_CZ_smallest so that it includes this C_l, and therefore it
        contains the actual squeezed CMB-lensing bispectrum rather than its kernel.
        The rescaling cannot be performed for the other Fisher arrays because they
        don't have information on the field C in the Fisher matrix sum. Therefore, we
        use fisher_matrix_CZ_smallest rather than fisher_matrix_XYZ_smallest to build
        the Fisher matrix. This prevents us to fill the array fisher_matrix_XYZ_lmin
        and we thus leave it zero-valued. */

        if (pfi->include_lensing_effects == _TRUE_) {

          for (int C=0; C < pfi->ff_size; ++C) {
            for (int Z=0; Z < pfi->ff_size; ++Z) {
                  
              double l1_contribution = pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z];
              pfi->fisher_matrix_smallest[index_l1][index_ft_1][index_ft_2] += l1_contribution;
              if (!pfi->l1_is_full)
                l1_contribution *= pwf->delta_l[index_l1];
              accumulator += l1_contribution;
              pfi->fisher_matrix_lmin[index_l1][index_ft_1][index_ft_2] = accumulator;
  
            }
          }
        }
        
        else {

          for (int X = 0; X < pfi->ff_size; ++X) {
            for (int Y = 0; Y < pfi->ff_size; ++Y) {
              for (int Z = 0; Z < pfi->ff_size; ++Z) {

                /* Contribution from this l1 to the total (X+Y+Z) Fisher matrix */
                double l1_contribution = pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft_1][index_ft_2];
                    
                /* Sum over  all possible XYZ */
                pfi->fisher_matrix_smallest[index_l1][index_ft_1][index_ft_2] += l1_contribution;

                /* Include the interpolation weight for the l1 direction, if needed, in order
                to correcly sum over l1 (not needed for the sum over XYZ, above) */
                if (pfi->l1_size != pbi->full_l_size)
                  l1_contribution *= pwf->delta_l[index_l1];

                /* Contribution of all the l's larger than l1 to the total (X+Y+Z) Fisher matrix */
                accumulator += l1_contribution;
                pfi->fisher_matrix_lmin[index_l1][index_ft_1][index_ft_2] = accumulator;

                /* Contribution of all the l's larger than l1 to the Fisher matrix of XYZ */
                accumulator_XYZ[X][Y][Z] += l1_contribution;
                pfi->fisher_matrix_XYZ_lmin[X][Y][Z][index_l1][index_ft_1][index_ft_2] = accumulator_XYZ[X][Y][Z];

              } // Z
            } // Y
          } // X
        } // if not lensing effects
            
        /* Detectability of the bispectrum, assuming that the bispectrum is a
        template with f_nl = 1 */
        if (index_ft_1 == index_ft_2)
          pfi->sigma_fnl_lmin[index_l1][index_ft_1] = 1/sqrt(pfi->fisher_matrix_lmin[index_l1][index_ft_1][index_ft_1]);

        /* Print fnl(l_min) */
        // fprintf (stderr, "%12d %12g\n", pfi->l1[index_l1], pfi->sigma_fnl_lmin[index_l1][0]);

      } // end of for(index_l1)
            
      /* Compute the inverse Fisher matrix */
      for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1)
        InverseMatrix (pfi->fisher_matrix_lmin[index_l1], pfi->ft_size, pfi->inverse_fisher_matrix_lmin[index_l1]);

            
      // -------------------------------------------------------------------------------
      // -                            As a function of lmax                            -
      // -------------------------------------------------------------------------------
        
      /* Build the Fisher matrix as a function of the maximum l in the survey lmax, as:
      sum_{lmin<=l3<=lmax} fisher_matrix_XYZ_smallest */
      accumulator = 0;
      for (int X = 0; X < pfi->ff_size; ++X)
        for (int Y = 0; Y < pfi->ff_size; ++Y)
          for (int Z = 0; Z < pfi->ff_size; ++Z)
            accumulator_XYZ[X][Y][Z] = 0.;

      for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3) {
        for (int X = 0; X < pfi->ff_size; ++X) {
          for (int Y = 0; Y < pfi->ff_size; ++Y) {
            for (int Z = 0; Z < pfi->ff_size; ++Z) {
      
              double l3_contribution = pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft_1][index_ft_2];
              pfi->fisher_matrix_largest[index_l3][index_ft_1][index_ft_2] += l3_contribution;
              accumulator += l3_contribution;
              pfi->fisher_matrix_lmax[index_l3][index_ft_1][index_ft_2] = accumulator;
              accumulator_XYZ[X][Y][Z] += l3_contribution;
              pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l3][index_ft_1][index_ft_2] = accumulator_XYZ[X][Y][Z];
      
            } // Z
          } // Y
        } // X
        
        if (index_ft_1 == index_ft_2)
          pfi->sigma_fnl_lmax[index_l3][index_ft_1] = 1/sqrt(pfi->fisher_matrix_lmax[index_l3][index_ft_1][index_ft_1]);

        /* Print fnl(l_max) */
        // fprintf (stderr, "%12d %12g\n", pfi->l3[index_l3], pfi->sigma_fnl_lmax[index_l3][0]);

      } // end of for(index_l3)
      
    } // end of for(ft_2)
  } // end of for (ft_1)
    
  for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3)
    InverseMatrix (pfi->fisher_matrix_lmax[index_l3], pfi->ft_size, pfi->inverse_fisher_matrix_lmax[index_l3]);

  /* Simple consistency check: the sum from above (l3) should be equal to the sum from below (l1).
  The test cannot be performed for the lensing variance row/column of the Fisher matrix because,
  in that case, the l3 arrays are not rescaled while the l1 ones are. */
  for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
    for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {
  
      if ((pfi->include_lensing_effects==_TRUE_)
      && ((index_ft_1==pfi->index_ft_cmb_lensing_kernel)
      || (index_ft_2==pfi->index_ft_cmb_lensing_kernel)))
        continue;
       
      double diff = 1-pfi->fisher_matrix_lmax[pfi->l3_size-1][index_ft_1][index_ft_2]
        /pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_2];
      
      class_test (fabs (diff) > _SMALL_,
        pfi->error_message,
        "error in the Fisher matrix sum, diff=%g.", diff);
    }
  }

  /* Print the signal to noise as a function of l_max for all the analised bispectra */
  // fprintf (stderr, "%4s ", "l");
  // for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft)    
  //   fprintf (stderr, "%20s ", pfi->ft_labels[index_ft]);
  // fprintf (stderr, "\n");
  // 
  // for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
  //   fprintf (stderr, "%4d ", pfi->l1[index_l1]);
  //   for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft)
  //     fprintf (stderr, "%20.7g ", sqrt(pfi->fisher_matrix_lmax[index_l1][index_ft][index_ft]));
  //   fprintf (stderr, "\n");
  // }
  
  

  // ====================================================================================
  // =                                  Print results                                   =
  // ====================================================================================
    
  /* Print the Fisher matrix */
  if (pfi->fisher_verbose > 0) {

    sprintf (pfi->info, "%sFisher matrix for l_max = %d:\n",
      pfi->info, MIN (pfi->l_max_estimator, pbi->l_max));
    
    for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
      sprintf (pfi->info, "%s\t%20s\t", pfi->info, pfi->ft_labels[index_ft_1]);
      
      sprintf (pfi->info, "%s(", pfi->info);
      
      for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
        sprintf (pfi->info, "%s %13.6g ", pfi->info, pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_2]);

      sprintf (pfi->info, "%s)\n", pfi->info);
    }

    if (pfi->include_lensing_effects == _TRUE_) {
      
      sprintf (pfi->info_lensvar, "%sFisher matrix for l_max = %d:\n",
        pfi->info_lensvar, MIN (pfi->l_max_estimator, pbi->l_max));
    
      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
        sprintf (pfi->info_lensvar, "%s\t%20s\t", pfi->info_lensvar, pfi->ft_labels[index_ft_1]);
      
        sprintf (pfi->info_lensvar, "%s(", pfi->info_lensvar);
      
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
          sprintf (pfi->info_lensvar, "%s %13.6g ", pfi->info_lensvar, pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_2]);

        sprintf (pfi->info_lensvar, "%s)\n", pfi->info_lensvar);
      }
    }
  }
  
  /* Print the correlation matrix */
  if ((pfi->fisher_verbose > 0) && (pfi->ft_size>1)) {

    sprintf (pfi->info, "%sCorrelation matrix for l_max = %d:\n",
      pfi->info, MIN (pfi->l_max_estimator, pbi->l_max));
    
    for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
      sprintf (pfi->info, "%s\t%20s\t", pfi->info, pfi->ft_labels[index_ft_1]);
      
      sprintf (pfi->info, "%s(", pfi->info);
      
      for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
        sprintf (pfi->info, "%s %13.6g ", pfi->info,
        pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_2]
        /sqrt(pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_1]
        *pfi->fisher_matrix_lmin[0][index_ft_2][index_ft_2]));

      sprintf (pfi->info, "%s)\n", pfi->info);
    }

    if (pfi->include_lensing_effects == _TRUE_) {
      
      sprintf (pfi->info_lensvar, "%sCorrelation matrix for l_max = %d:\n",
        pfi->info_lensvar, MIN (pfi->l_max_estimator, pbi->l_max));
    
      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
        sprintf (pfi->info_lensvar, "%s\t%20s\t", pfi->info_lensvar, pfi->ft_labels[index_ft_1]);
      
        sprintf (pfi->info_lensvar, "%s(", pfi->info_lensvar);
      
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
          sprintf (pfi->info_lensvar, "%s %13.6g ", pfi->info_lensvar,
          pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_2]
          /sqrt(pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_1]
          *pfi->fisher_matrix_lensvar_lmin[0][index_ft_2][index_ft_2]));

        sprintf (pfi->info_lensvar, "%s)\n", pfi->info_lensvar);
      }
    }
  }
  
  /* Print the inverse matrix */
  if (pfi->fisher_verbose > 0) {
  
    sprintf (pfi->info, "%sReciprocal of the inverse Fisher matrix:\n", pfi->info);
    
    for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
      sprintf (pfi->info, "%s\t%20s\t", pfi->info, pfi->ft_labels[index_ft_1]);
      
      sprintf (pfi->info, "%s(", pfi->info);
      
      for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
        sprintf (pfi->info, "%s %13.6g ", pfi->info, 1/pfi->inverse_fisher_matrix_lmin[0][index_ft_1][index_ft_2]);
  
      sprintf (pfi->info, "%s)\n", pfi->info);
    }
    
    if (pfi->include_lensing_effects == _TRUE_) {
    
      sprintf (pfi->info_lensvar, "%sReciprocal of the inverse Fisher matrix:\n", pfi->info_lensvar);
    
      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
        sprintf (pfi->info_lensvar, "%s\t%20s\t", pfi->info_lensvar, pfi->ft_labels[index_ft_1]);
      
        sprintf (pfi->info_lensvar, "%s(", pfi->info_lensvar);
      
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
          sprintf (pfi->info_lensvar, "%s %13.6g ", pfi->info_lensvar,
          1/pfi->inverse_fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_2]);
  
        sprintf (pfi->info_lensvar, "%s)\n", pfi->info_lensvar);
      }
    }
  }
  
  /* Print the fnl matrix */
  if (pfi->fisher_verbose > 0)  {
  
    sprintf (pfi->info, "%sfnl matrix (diagonal: 1/sqrt(F_ii), upper: F_12/F_11, lower: F_12/F_22):\n", pfi->info);
    
    for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
      sprintf (pfi->info, "%s\t%20s\t", pfi->info, pfi->ft_labels[index_ft_1]);
      
      sprintf (pfi->info, "%s(", pfi->info);
      
      for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {
  
        /* Diagonal elements = 1/sqrt(F_ii) */ 
        if (index_ft_1==index_ft_2)
          sprintf (pfi->info, "%s %13.6g ", pfi->info, 1/sqrt(pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_1]));
        /* Upper triangle = F_12/F_11, lower triangle = F_12/F_22. */
        else
          sprintf (pfi->info, "%s %13.6g ", pfi->info, pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_2]
            /pfi->fisher_matrix_lmin[0][index_ft_1][index_ft_1]);
      }
      sprintf (pfi->info, "%s)\n", pfi->info);
    }
    
    if (pfi->include_lensing_effects == _TRUE_) {
    
      sprintf (pfi->info_lensvar, "%sfnl matrix:\n", pfi->info_lensvar);
    
      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
        sprintf (pfi->info_lensvar, "%s\t%20s\t", pfi->info_lensvar, pfi->ft_labels[index_ft_1]);
      
        sprintf (pfi->info_lensvar, "%s(", pfi->info_lensvar);
      
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {
  
          /* Diagonal elements = 1/sqrt(F_ii) */ 
          if (index_ft_1==index_ft_2)
            sprintf (pfi->info_lensvar, "%s %13.6g ", pfi->info_lensvar,
            1/sqrt(pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_1]));
          /* Upper triangle = F_12/F_11, lower triangle = F_12/F_22. */
          else
            sprintf (pfi->info_lensvar, "%s %13.6g ", pfi->info_lensvar, pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_2]
              /pfi->fisher_matrix_lensvar_lmin[0][index_ft_1][index_ft_1]);
        }
        sprintf (pfi->info_lensvar, "%s)\n", pfi->info_lensvar);
      }
    }
  }
  
  /* Print the Fisher matrix for all the bispectra (TTT,TTE,TET...). This tells us which bispectrum
  contributes the most to the total signal-to-noise. */
  if ((pfi->fisher_verbose > 2) && (pfi->ff_size > 1)) {
  
    for (int X = 0; X < pfi->ff_size; ++X) {
      for (int Y = 0; Y < pfi->ff_size; ++Y) {
        for (int Z = 0; Z < pfi->ff_size; ++Z) {
  
          sprintf (pfi->info, "%s     * %s contribution (percent)\n", pfi->info, pfi->ffff_labels[X][Y][Z]);
  
          for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      
            sprintf (pfi->info, "%s\t%20s\t(", pfi->info, "");
      
            for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
              sprintf (pfi->info, "%s %8.3g ", pfi->info,
                pfi->fisher_matrix_XYZ_lmax[X][Y][Z][pfi->l3_size-1][index_ft_1][index_ft_2]/
                pfi->fisher_matrix_lmax[pfi->l3_size-1][index_ft_1][index_ft_2]*100);
  
            sprintf (pfi->info, "%s)\n", pfi->info);
          }
        }
      }
    }
  }
  
  /* Print to screen the Fisher matrix */
  printf_log ("\n");
  printf_log ("%s", pfi->info);
  printf_log ("\n");
  
  if (pfi->include_lensing_effects == _TRUE_) {
    printf_log (" -> optimal estimator including lensing variance:\n");
    printf_log ("\n");
    printf_log ("%s", pfi->info_lensvar);
    printf_log ("\n");
  }
    
  /* Free memory */
  free (pwf->delta_l);
  free (pwf);
  
  return _SUCCESS_;
  
}



/**
 * Computes the Fisher matrix for all bispectra types using the trapezoidal
 * integration along the l1 direction.
 *
 * The bispectra in the l2 and l3 directions are computed using interpolation.
 * Two methods are supported:
 * 
 * - Bilinear interpolation via the bispectra_at_l2l3_bilinear_bilinear() function.
 * - Mesh interpolation via the bispectra_at_l2l3_bilinear_mesh() function.
 * 
 * The following arrays are filled:
 *
 * - pfi->fisher_matrix_XYZ_largest
 * - pfi->fisher_matrix_XYZ_smallest
 * - pfi->fisher_matrix_CZ_smallest
 */

int fisher_compute_matrix (
       struct precision * ppr,
       struct transfers * ptr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi,
       struct fisher_workspace * pwf
       )
{  

  /* We shall count the (l1,l2,l3) configurations over which we compute the estimator,
  and the configurations that we shall skip */
  long int counter = 0;
  long int counter_skipped = 0;

  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel private (thread)
  #ifdef _OPENMP
  number_of_threads = omp_get_num_threads();
  #endif

  /* Temporary array to store the interpolated bispectrum */
  double (*interpolated_bispectra)[pfi->ft_size][pfi->ff_size][pfi->ff_size][pfi->ff_size]
    = malloc (number_of_threads * pfi->ft_size * pfi->ff_size * pfi->ff_size * pfi->ff_size * sizeof(double));
    
  #pragma omp parallel shared (abort,counter) private (thread)
  {
    
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    // -------------------------------------------------------------------------------
    // -                                  Sum over l1                                -
    // -------------------------------------------------------------------------------

    #pragma omp for schedule (dynamic)
    for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
      int l1 = pfi->l1[index_l1];
    
      if ((l1 < pfi->l1_min_global) || (l1 > pfi->l1_max_global))
        continue;

      /* Arrays that will contain the 3j symbols for a given (l1,l2) */
      double threej_000[2*pbi->l_max+1];

      /* First l-multipole stored in the above arrays */
      int l3_min_000=0;
    
      printf_log_if (pfi->fisher_verbose, 1, 
        "     * computing Fisher matrix for l1=%d\n", l1);

      // -------------------------------------------------------------------------------
      // -                                  Sum over l2                                -
      // -------------------------------------------------------------------------------
    
      for (int l2=l1; l2 <= pbi->l_max; ++l2) {
  
        if ((l2 < pfi->l2_min_global) || (l2 > pfi->l2_max_global))
          continue;

        // -------------------------------------------------------------------------------
        // -                                Compute 3J                                   -
        // -------------------------------------------------------------------------------
      
        /* Compute the 3j-symbol that enters the estimator, (l1 l2 l3)(0 0 0) */
        double min_D, max_D;
  
        class_call_parallel (drc3jj (
                               l1, l2, 0, 0,
                               &min_D, &max_D,
                               threej_000,
                               (2*pbi->l_max+1),
                               pfi->error_message),
          pfi->error_message,
          pfi->error_message);          

        l3_min_000 = (int)(min_D + _EPS_);


        // -------------------------------------------------------------------------------
        // -                                 Sum over l3                                 -
        // -------------------------------------------------------------------------------

        /* The estimator only includes even contributions of (l1+l2+l3), as parity invariance kills the others */
        for (int l3 = l2 + l1%2; l3 <= MIN(pbi->l_max,l1+l2); l3=l3+2) {
              
          if ((l3 < pfi->l3_min_global) || (l3 > pfi->l3_max_global))
            continue;

          /* Consider only squeezed or equilateral configurations according to pfi->squeezed_ratio */
          if (((pfi->squeezed_ratio>1)&&(l2<pfi->squeezed_ratio*l1))
          || ((pfi->squeezed_ratio<-1)&&(l3>fabs(pfi->squeezed_ratio)*l1))) {
            #pragma omp atomic
            counter++;
            #pragma omp atomic
            counter_skipped++;
            continue;
          }

          /* Uncomment to consider only squeezed configurations, with l1 superhorizon at recombination.
          This matches the "Maldacena" limit used in Creminelli et al. 2004. */
          // if (!((l1<=50) && (l2>=20*l1))) {
          //   #pragma omp atomic
          //   counter++;
          //   continue;
          // }        


          // ---------------------------------------------------------------------------------
          // -                                 Obtain bispectra                              -
          // ---------------------------------------------------------------------------------

          /* Compute 3J ratios, needed only for polarised bispectra such as CMB-lensing
          and the quadratic correction. Note that with respect to the bispectrum module,
          here we switch the m's below in the l1 and l3 columns, as now the smallest multipole
          is l1 and not l3. */
          double threej_ratio_20m2, threej_ratio_m220, threej_ratio_0m22;

          if (pbi->need_3j_symbols == _TRUE_) {          

            class_call_parallel (threej_ratio_M (l2, l3, l1, 2, &threej_ratio_20m2, pbi->error_message),
              pbi->error_message, pbi->error_message);

            class_call_parallel (threej_ratio_M (l1, l3, l2, 2, &threej_ratio_m220, pbi->error_message),
              pbi->error_message, pbi->error_message);

            class_call_parallel (threej_ratio_M (l3, l2, l1, 2, &threej_ratio_0m22, pbi->error_message),
              pbi->error_message, pbi->error_message);

          }

          /* Factor that relates the reduced bispectrum to the angle-averaged one */
          double I_l1_l2_l3 = sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.)/(4.*_PI_)) * threej_000[l3-l3_min_000];

          /* Obtain all the bispectra in this (l1,l2,l3) configuration, either via interpolation or
          by direct computation (if the bispectrum can be expressed in a simple analytical form).
          These loops go over all the possible bispectra, eg. local_ttt, equilateral_ete, intrinsic_ttt. */
          for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
            
            int index_bt = pfi->index_bt_of_ft[index_ft];
            
            for (int X = 0; X < pfi->ff_size; ++X) {
              for (int Y = 0; Y < pfi->ff_size; ++Y) {
                for (int Z = 0; Z < pfi->ff_size; ++Z) {

                  /* Uncomment to ignore high-l configurations for temperature */
                  // int l_max_T = 0;
                  // if (pbi->has_bispectra_t) {
                  //   if (((X==pbi->index_bf_t) && (l3>l_max_T))
                  //   || ((Y==pbi->index_bf_t) && (l2>l_max_T))
                  //   || ((Z==pbi->index_bf_t) && (l1>l_max_T))) {
                  //     interpolated_bispectra[thread][index_ft][X][Y][Z] = 0;
                  //     continue;
                  //   }
                  // }

                  /* Uncomment to ignore high-l configurations for polarisation */
                  // int l_max_E = 0;
                  // if (pbi->has_bispectra_e) {
                  //   if (((X==pbi->index_bf_e) && (l3>l_max_E))
                  //   || ((Y==pbi->index_bf_e) && (l2>l_max_E))
                  //   || ((Z==pbi->index_bf_e) && (l1>l_max_E))) {
                  //     interpolated_bispectra[thread][index_ft][X][Y][Z] = 0;
                  //     continue;
                  //   }
                  // }

                  /* Corresponding field indices in the bispectrum module */
                  int X_ = pfi->index_bf_of_ff[X];
                  int Y_ = pfi->index_bf_of_ff[Y];
                  int Z_ = pfi->index_bf_of_ff[Z];


                  /* Compute analytical bispectra */

                  if (!pbi->interpolate_me[index_bt]) {

                    /* Check that the current bispectrum has a function associated to it */
                    class_test_parallel (pbi->bispectrum_function[index_bt]==NULL,
                      pbi->error_message,
                      "no function associated for the bispectrum '%s'. Maybe it's not analytical?",
                      pbi->bt_labels[index_bt]);

                    class_call_parallel ((*pbi->bispectrum_function[index_bt]) (
                                             ppr, psp, ple, pbi,
                                             l3, l2, l1, /* smallest one goes in third position  */
                                             X_, Y_, Z_, /* need underscore because they are bispectra field indices */
                                             pbi->lens_me[index_bt],
                                             threej_ratio_20m2,
                                             threej_ratio_m220,
                                             threej_ratio_0m22,
                                             &interpolated_bispectra[thread][index_ft][X][Y][Z]),
                      pbi->error_message,
                      pfi->error_message);

                  } 

                  /* Interpolate all other bispectra */

                  else {
                    
                    if (pbi->interpolation_method == bilinear_interpolation) {

                      class_call_parallel (bispectra_at_l2l3_bilinear (
                                             ptr,
                                             pbi,
                                             index_bt,
                                             index_l1, l2, l3,
                                             Z, Y, X,
                                             _TRUE_,
                                             &interpolated_bispectra[thread][index_ft][X][Y][Z],
                                             NULL),
                        pbi->error_message,
                        pfi->error_message);
                    }

                    else if (pbi->interpolation_method == mesh_interpolation_2D) {

                      class_call_parallel (bispectra_mesh_interpolate(
                                             ppr, psp, ple, pbi,
                                             index_bt,
                                             index_l1,
                                             l3, l2, l1, /* smallest one goes in third position  */
                                             X_, Y_, Z_,
                                             pbi->meshes[index_l1][index_bt][X_][Y_][Z_][0],
                                             pbi->meshes[index_l1][index_bt][X_][Y_][Z_][1],
                                             &interpolated_bispectra[thread][index_bt][X_][Y_][Z_]),
                        pbi->error_message,
                        pfi->error_message);
                    }

                    else if (pbi->interpolation_method == mesh_interpolation_3D) {

                      class_call_parallel (bispectra_mesh_interpolate(
                                             ppr, psp, ple, pbi,
                                             index_bt,
                                             -1,
                                             l3, l2, l1, /* smallest one goes in third position  */
                                             X_, Y_, Z_,
                                             pbi->meshes[0][index_bt][X_][Y_][Z_][0], /* There is one 3D mesh for all l1 values */
                                             pbi->meshes[0][index_bt][X_][Y_][Z_][1], /* There is one 3D mesh for all l1 values */
                                             &interpolated_bispectra[thread][index_bt][X_][Y_][Z_]),
                        pbi->error_message,
                        pfi->error_message);
                    }

                    /* Debug - print the interpolated values */
                    // if ((l1==1000)) {
                    // // if ((l1==63)&&(l2==38)&&(l3==27)) {
                    //   printf ("INTERPOLATED: %s_%s(%d,%d,%d) = %g\n",
                    //     pfi->ft_labels[index_ft], pfi->ffff_labels[X][Y][Z], l1, l2, l3,
                    //     interpolated_bispectra[thread][index_ft][X][Y][Z]);
                    // }

                    /* Compensate the effect of the window function */
                    double inverse_window = 1;
                
                    if (pbi->window_function[index_bt] != NULL)
                      class_call_parallel ((*pbi->window_function[index_bt]) (
                                              ppr, psp, ple, pbi,
                                              l3, l2, l1, /* smallest one goes in third position  */
                                              X_, Y_, Z_,
                                              &inverse_window),
                        pbi->error_message,
                        pfi->error_message);

                    interpolated_bispectra[thread][index_ft][X][Y][Z] *= inverse_window;
                    
                  } 
              
                  /* All bispectra have to be multiplied by the parity factor */
                  interpolated_bispectra[thread][index_ft][X][Y][Z] *= I_l1_l2_l3;

                } // end of for(Z)
              } // end of for(Y)
            } // end of for(X)
          } // end of for(index_ft)
          

          // -------------------------------------------------------------------------------
          // -                             Build the Fisher matrix                         -
          // -------------------------------------------------------------------------------

          /* To compute the covariance matrix, we use the formula first given in eq. 17 of Yadav,
          Komatsu & Wandelt 2007 (http://uk.arxiv.org/abs/astro-ph/0701921), which is a numerical
          improvement over the treatment in Babich & Zaldarriaga 2004. (Note that explicit form
          given in appendix B of the same paper is valid only when l1=l2=l3, as noted in appendix E
          of Lewis, Challinor & Hanson 2011.)
          For the case where only one probe is requested, the below loops reduce to one iteration
          over a single bispectrum (eg. T -> TTT). When two probes are requested, each loop corresponds
          to one of 2^3 bispectra (eg. TTT,TTE,TET,ETT,TEE,ETE,EET,EEE), for a total of 64 elements in
          the covariance matrix. */

          /* Shortcuts to the inverse of the cross-power spectrum */
          double ** C1 = pfi->inverse_C[l1-2];
          double ** C2 = pfi->inverse_C[l2-2];
          double ** C3 = pfi->inverse_C[l3-2];
          
          /* Correction factor for the variance */
          double one_over_delta = 1;
          if ((l1==l2) && (l1==l3))
            one_over_delta = 1/6.;
          else if ((l1==l2) || (l1==l3) || (l2==l3))
            one_over_delta = 0.5;
                    
          /* Compute the Fisher matrix. In the simple case of only one probe (eg. temperature -> TTT),
          each element of the Fisher matrix is simply given by the square of the bispectrum divided
          by three C_l's. When dealing with two probes (eg. T and E -> TTT,TTE,TET,ETT,TEE,ETE,EET,EEE),
          each element of the Fisher matrix, contains a quadratic sum over the 8 available bispectra,
          for a total of 64 terms. */
          for (int X = 0; X < pfi->ff_size; ++X) { for (int A = 0; A < pfi->ff_size; ++A) {
          for (int Y = 0; Y < pfi->ff_size; ++Y) { for (int B = 0; B < pfi->ff_size; ++B) {
          for (int Z = 0; Z < pfi->ff_size; ++Z) { for (int C = 0; C < pfi->ff_size; ++C) {

            double inverse_covariance = C1[C][Z] * C2[B][Y] * C3[A][X];

            for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
              for (int index_ft_2=index_ft_1; index_ft_2 < pfi->ft_size; ++index_ft_2) {              

                double fisher = one_over_delta *
                                inverse_covariance *
                                interpolated_bispectra[thread][index_ft_1][A][B][C] *
                                interpolated_bispectra[thread][index_ft_2][X][Y][Z];

                /* Fisher matrix as a function of the largest multipole. Since this quantity is summed
                over the incomplete l1 direction, we have to include the interpolation weight. */
                if (!pfi->l1_is_full) {
                  #pragma omp atomic
                  pfi->fisher_matrix_XYZ_largest[X][Y][Z][l3-2][index_ft_1][index_ft_2] += fisher * pwf->delta_l[index_l1];
                }
                else {
                  #pragma omp atomic
                  pfi->fisher_matrix_XYZ_largest[X][Y][Z][l3-2][index_ft_1][index_ft_2] += fisher;
                }

                /* Fisher matrix as a function of the smallest multipole */
                pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft_1][index_ft_2] += fisher;
                
                /* Increment the array needed to compute the lensing variance. This is equivalent to
                \bar{F} in Eq. 5.25 of http://uk.arxiv.org/abs/1101.2234. We do not include the
                interpolation weight because we have to process the array before summing over it.  */
                if (pfi->include_lensing_effects == _TRUE_) {
                  
                  pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z]
                    += fisher;
                  
                  /* Store information on l3 in order to compute the lensing variance as a function of l_max */
                  if (pfi->compute_lensing_variance_lmax == _TRUE_)
                    pfi->fisher_matrix_CZ_smallest_largest[index_l1][l3-2][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z]
                      += fisher;
 
                } // end of if(include_lensing_effects)
                
              } // ft_2
            } // ft_1
          }}} // XYZ
          }}} // ABC
                
          #pragma omp atomic
          counter++;          

        } // end of for(l3)
      } // end of for(l2)
      
      /* Free the memory associated to the 2D mesh for this l1 value. If you will
      be needing to interpolate the bispectra again, comment to gain a speed up. */
      if (pbi->interpolation_method == mesh_interpolation_2D)
        class_call_parallel (bispectra_mesh_empty (pbi, &pbi->meshes[index_l1]),
          pbi->error_message,
          pfi->error_message);

    } // for(index_l1)
    
  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region
  
  free (interpolated_bispectra);
  
  if ((pfi->fisher_verbose > 1) && (counter_skipped > 0)) {
    printf (" -> skipped %.3g percent of %ld configurations whereby ",
      100*(double)counter_skipped/counter, counter);
    if (pfi->squeezed_ratio > 1)
      printf ("L2/L1 < %g", pfi->squeezed_ratio);
    else if (pfi->squeezed_ratio < -1)
      printf ("L3/L1 > %g", fabs(pfi->squeezed_ratio));
    printf ("\n");
  }
  

  /* Include the effect of partial sky coverage and symmetrise the Fisher matrix */
  class_call (fisher_sky_coverage(pfi),
    pfi->error_message,
    pfi->error_message);


  return _SUCCESS_;
  
}





/**
 * 
 * This function computes the Fisher matrix and fills pfi->fisher_matrix_XYZ_largest and 
 * pfi->fisher_matrix_XYZ_smallest using a linear interpolation of the bispectrum.
 */
int fisher_compute_matrix_nodes (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi,
       struct fisher_workspace * pwf
       )
{  

  /* We shall count the (l1,l2,l3) configurations over which we compute the estimator */
  long int counter = 0;

  // =======================================================================================
  // =                            Compute Fisher matrix elements                           =
  // =======================================================================================

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel private (thread)
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    // ------------------------------------------------
    // -                  Sum over l1                 -
    // ------------------------------------------------

    #pragma omp for schedule (dynamic)
    for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
      int l1 = pfi->l1[index_l1];
    
      if ((l1 < pfi->l1_min_global) || (l1 > pfi->l1_max_global))
        continue;
    
      printf_log_if (pfi->fisher_verbose, 2, 
        "     * computing Fisher matrix for l1=%d\n", l1);

      // ------------------------------------------------
      // -                  Sum over l2                 -
      // ------------------------------------------------
    
      for (int index_l2=0; index_l2<=index_l1; ++index_l2) {
        
        int l2 = pfi->l2[index_l2];
  
        if ((l2 < pfi->l2_min_global) || (l2 > pfi->l2_max_global))
          continue;
        
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (l2 < l1/2)
          continue;
        
        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
        
        if (pbi->interpolation_method == smart_interpolation) {

          class_call_parallel (fisher_interpolation_weights(
                                 ppr,
                                 psp,
                                 ple,
                                 pbi,
                                 pfi,
                                 index_l1,
                                 index_l2,
                                 pwf->delta_l3[thread],
                                 pwf),
            pfi->error_message,
            pfi->error_message);
        }

        // ------------------------------------------------
        // -                  Sum over l3                 -
        // ------------------------------------------------

        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
        
          int l3 = pfi->l3[index_l3];
                                
          if ((l3 < pfi->l3_min_global) || (l3 > pfi->l3_max_global))
            continue;

          /* Parity invariance enforces that the contribution of the modes with odd l1+l2+l3 is zero. Even if you
          comment this out, the result won't change as the I_l1_l2_l3 coefficient would vanish in that case. */
          if ((l1+l2+l3)%2!=0) {
            #pragma omp atomic
            counter++;
            continue;
          }

          /* Index of the current (l1,l2,l3) configuration */
          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

          // ------------------------------------------------------------------------
          // -                        Build the Fisher matrix                       -
          // ------------------------------------------------------------------------

          /* To compute the covariance matrix, we use the formula first given in eq. 17 of Yadav,
          Komatsu & Wandelt 2007 (http://uk.arxiv.org/abs/astro-ph/0701921), which is a numerical
          improvement over the treatment in Babich & Zaldarriaga 2004. (Note that explicit form
          given in appendix B of the same paper is valid only when l1=l2=l3, as noted in appendix E
          of Lewis, Challinor & Hanson 2011.)
          For the case where only one probe is requested, the below loops reduce to one iteration
          over a single bispectrum (eg. T -> TTT). When two probes are requested, each loop corresponds
          to one of 2^3 bispectra (eg. TTT,TTE,TET,ETT,TEE,ETE,EET,EEE), for a total of 64 elements in
          the covariance matrix. */
          
          /* Shortcuts to the inverse of the cross-power spectrum */
          double ** C1 = pfi->inverse_C[l1-2];
          double ** C2 = pfi->inverse_C[l2-2];
          double ** C3 = pfi->inverse_C[l3-2];
          
          /* Weights for the interpolation. Its default value is unity, in case interpolation is not
          needed. */
          double interpolation_weight = 1;
          
          if (pbi->interpolation_method == smart_interpolation) {
            interpolation_weight = pwf->delta_l[index_l1] * pwf->delta_l[index_l2] * pwf->delta_l3[thread][index_l3];
          }
          else if (pbi->interpolation_method == trilinear_interpolation) {
            /* Weight for trilinear interpolation. When the trilinear interpolation is turned on, we always consider
            configurations where all l's are even. The factor 1/2 comes from the fact that when interpolating over even
            configurations only, we are implicitly assuming non-zero values for the odd points, which instead are
            forced to vanish by the 3J. (Note that this factor, in the case of the smart interpolation, is already
            included in pwf->delta_l3) */
            interpolation_weight = 0.5 * pwf->delta_l[index_l1] * pwf->delta_l[index_l2] * pwf->delta_l[index_l3];
          }
          
          /* If we are taking all the l-points, then any interpolation scheme reduces to a simple sum over the
          multipoles, and then we can naturally implement the delta factor in KSW2005. */
          double one_over_delta = 1;
          if (ppr->l_linstep == 1) {

            if ((l1==l2) && (l1==l3))
              one_over_delta = 1/6.;

            else if ((l1==l2) || (l1==l3) || (l2==l3))
              one_over_delta = 0.5;
          }
          
          /* Include the double 3j symbol */
          double threej_000_squared = pfi->I_l1_l2_l3[index_l1_l2_l3] * pfi->I_l1_l2_l3[index_l1_l2_l3];
          
          /* Compute the Fisher matrix. In the simple case of only one probe (eg. temperature -> TTT),
          each element of the Fisher matrix is simply given by the square of the bispectrum divided
          by three C_l's. When dealing with two probes (eg. T and E -> TTT,TTE,TET,ETT,TEE,ETE,EET,EEE),
          each element of the Fisher matrix, contains a quadratic sum over the 8 available bispectra,
          for a total of 64 terms. */
          for (int X = 0; X < pfi->ff_size; ++X) { for (int A = 0; A < pfi->ff_size; ++A) {
          for (int Y = 0; Y < pfi->ff_size; ++Y) { for (int B = 0; B < pfi->ff_size; ++B) {
          for (int Z = 0; Z < pfi->ff_size; ++Z) { for (int C = 0; C < pfi->ff_size; ++C) {

            double inverse_covariance = C1[A][X] * C2[B][Y] * C3[C][Z];

            for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
              for (int index_ft_2=index_ft_1; index_ft_2 < pfi->ft_size; ++index_ft_2) {
                
                double b_1 = pbi->bispectra[pfi->index_bt_of_ft[index_ft_1]]
                             [pfi->index_bf_of_ff[A]][pfi->index_bf_of_ff[B]][pfi->index_bf_of_ff[C]][index_l1_l2_l3];

                double b_2 = pbi->bispectra[pfi->index_bt_of_ft[index_ft_2]]
                             [pfi->index_bf_of_ff[X]][pfi->index_bf_of_ff[Y]][pfi->index_bf_of_ff[Z]][index_l1_l2_l3];

                double fisher = one_over_delta *
                                interpolation_weight *
                                threej_000_squared *
                                inverse_covariance *
                                b_1 * b_2;
                                
                #pragma omp atomic
                pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l1][index_ft_1][index_ft_2] += fisher;
                #pragma omp atomic
                pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l3][index_ft_1][index_ft_2] += fisher;
                  
              } // bt_2
            } // bt_1
          }}} // XYZ
          }}} // ABC
                
          #pragma omp atomic
          counter++;          

        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_l1)
  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region
      
  /* Include the effect of partial sky coverage and symmetrise the Fisher matrix */
  class_call (fisher_sky_coverage(pfi), pfi->error_message, pfi->error_message);
      
  return _SUCCESS_;
  
}




int fisher_interpolation_weights (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi,
       int index_l1,
       int index_l2,
       double * delta_l3,
       struct fisher_workspace * pwf
       )
{  
  
  int l1 = pbi->l[index_l1];
  int l2 = pbi->l[index_l2];

  /* Limits of our l3 sampling */
  int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
  int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      
  // -----------------------------------------------------------------------------------------
  // -                            Exclude points with l1+l2+l3 odd                           -
  // -----------------------------------------------------------------------------------------
      
  // if ((l1==40) && (l2==40)) {
  //   printf("(l1=%d,l2=%d) BEFORE:", l1, l2);
  //   for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
  //     printf("%d,",pbi->l[index_l3]);
  //   printf("\n");
  // }
      
  /* Restrict the l3 range so that the first and last points are both even (if l1+l2 is even) or odd
  (if l1+l2 is odd). If this cannot be achieved, it means that our l-sampling fails to capture a valid
  point for the considered (l1,l2) pair. Then, the Fisher matrix computation will neglect this pair. */
  if ((l1+l2)%2==0) {
    while ((pbi->l[index_l3_min]%2!=0) && (index_l3_min<=index_l3_max))
      index_l3_min++;
    while ((pbi->l[index_l3_max]%2!=0) && (index_l3_max>=index_l3_min))
      index_l3_max--;
  }
  else {
    while ((pbi->l[index_l3_min]%2==0) && (index_l3_min<=index_l3_max))
      index_l3_min++;
    while ((pbi->l[index_l3_max]%2==0) && (index_l3_max>=index_l3_min))
      index_l3_max--;
  }

  // if ((l1==40) && (l2==40)) {
  //   printf("(l1=%d,l2=%d) AFTER:", l1, l2);
  //   for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
  //     printf("%d,",pbi->l[index_l3]);
  //   printf("\n");
  // }

  // ------------------------------------------------------------------------------------------
  // -                                Generate linear weights                                 -
  // ------------------------------------------------------------------------------------------
  
  /* First, we generate "normal" interpolation weights, as if the l3-range was rectangular.
  We also include an extra 1/2 factor to account for the fact that the Fisher estimator
  only includes those configuration where l1+l2+l3 is even. */
  /* TODO: This formula is valid only if l1+l2 is even. When l1+l2 is odd, then l3 must be
  odd to and the formula has to be adjusted */
  if (index_l3_min < index_l3_max) {

    delta_l3[index_l3_min] = (pbi->l[index_l3_min+1] - pbi->l[index_l3_min] + 2.)/4.;
      
    for (int index_l=index_l3_min+1; index_l<=(index_l3_max-1); ++index_l)
      delta_l3[index_l] = (pbi->l[index_l+1] - pbi->l[index_l-1])/4.;
      
    delta_l3[index_l3_max] = (pbi->l[index_l3_max] - pbi->l[index_l3_max-1] + 2.)/4.;

   if (index_l3_min == index_l3_max)
    delta_l3[index_l3_min] = 1;
  }
  /* If only one point is valid, its weight must be unity */
  else if (index_l3_min == index_l3_max) {
    delta_l3[index_l3_min] = 1;
  }
  /* If there are no usable points at all, it really does not matter what we do here, because the (l1,l2)
  configuration will be skipped anyway during the Fisher matrix l-sum. */
  else {
    return _SUCCESS_;
  }

  /* Some debug */
  // for (int index_l=index_l3_min; index_l <= (index_l3_max); ++index_l)
  //   printf ("delta_l3[index_l] = %g\n", delta_l3[index_l]);

  
  // ----------------------------------------------------------------------------------
  // -                          Account for out-of-bonds points                       -
  // ----------------------------------------------------------------------------------
  
  /* We assign to those points outside our l3-sampling (but inside the triangular condition) the
  bispectrum-value in l3. To do so, however, we need to properly count which of these
  out-of-bounds points are even and odd */
  for (int l3 = MAX(pbi->l_min, abs(l1-l2)); l3 < pbi->l[index_l3_min]; ++l3)
    if((l1+l2+l3)%2==0)
      delta_l3[index_l3_min] += 1.;
  
  for (int l3 = pbi->l[index_l3_max]+1; l3 <= MIN(l2,l1+l2); ++l3)
    if((l1+l2+l3)%2==0)
      delta_l3[index_l3_max] += 1;
  
  /* Check that, for a function that always evaluates to one, the weights sum up to the number
  of even/odd points in the integration range. For a given pair of (l1,l2), the range goes
  from MAX(2,|l1-l2|) to MIN(l2,l1+l2). 
  Remember that now 'index_l3_min' and 'index_l3_max' have been modified with respect to those in
  pbi->index_l_triangular_min and pbi->index_l_triangular_max to be even/odd. */
  double sum = 0;
  for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
    sum += delta_l3[index_l3];
  
  int n_even=0, n_odd=0;
  for (int l3 = MAX(pbi->l_min,abs(l1-l2)); l3 <= MIN(l2,l1+l2); ++l3) {
    if(l3%2==0) n_even++;
    else n_odd++;
  }
  
  if ((l1+l2)%2==0) {
    class_test (sum != n_even,
      pfi->error_message,
      "EVEN: l1=%d, l2=%d, sum=%g, expected=%d\n", l1, l2, sum, n_even);
  }
  else {
    class_test (sum != n_odd,
      pfi->error_message,
      "ODD:  l1=%d, l2=%d, sum=%g, expected=%d\n", l1, l2, sum, n_odd);
  }

  return _SUCCESS_;

}













/**
 * Account for partial sky coverage by multiplying the Fisher matrix by the sky fraction.
 * Also, symmetrise the Fisher matrix.
 */
int fisher_sky_coverage (
       struct fisher * pfi
       )
{  

  for (int X = 0; X < pfi->ff_size; ++X) {
    for (int Y = 0; Y < pfi->ff_size; ++Y) {
      for (int Z = 0; Z < pfi->ff_size; ++Z) {
    
        for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
          for (int index_ft_2=index_ft_1; index_ft_2 < pfi->ft_size; ++index_ft_2) {

            for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3) {
              pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft_1][index_ft_2] *= pfi->f_sky;
              pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft_2][index_ft_1]
                = pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l3][index_ft_1][index_ft_2];   
            }
            
            /* Do the same for the Fisher matrix as a function of the smallest multipole */
            for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
              pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft_1][index_ft_2] *= pfi->f_sky;
              pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft_2][index_ft_1]
                = pfi->fisher_matrix_XYZ_smallest[X][Y][Z][index_l1][index_ft_1][index_ft_2]; 
  }}}}}}

  /* Do the same for F_bar(l3,C,Z), which is needed to compute the lensing variance correction */
  if (pfi->include_lensing_effects == _TRUE_) {
  
    for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
      for (int Z = 0; Z < pfi->ff_size; ++Z) {
        for (int C = 0; C < pfi->ff_size; ++C) {
          for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
            for (int index_ft_2=index_ft_1; index_ft_2 < pfi->ft_size; ++index_ft_2) {              
              pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z] *= pfi->f_sky;
              pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_2*pfi->ff_size+C][index_ft_1*pfi->ff_size+Z]
                = pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_1*pfi->ff_size+Z][index_ft_2*pfi->ff_size+C];
    }}}}}

    /* Previously (wrong because you apply f_sky twice) */
    // for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    //   for (int i=0; i < pfi->ff_size*pfi->ft_size; ++i) {
    //     for (int j=0; j < pfi->ff_size*pfi->ft_size; ++j) {
    //         pfi->fisher_matrix_CZ_smallest[index_l1][i][j] *= pfi->f_sky;
    //         pfi->fisher_matrix_CZ_smallest[index_l1][j][i] = pfi->fisher_matrix_CZ_smallest[index_l1][i][j];
    // }}}

    if (pfi->compute_lensing_variance_lmax == _TRUE_) {
      for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
        for (int index_l3=0; index_l3<pfi->l3_size; ++index_l3) {
          for (int Z = 0; Z < pfi->ff_size; ++Z) {
            for (int C = 0; C < pfi->ff_size; ++C) {
              for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
                for (int index_ft_2=index_ft_1; index_ft_2 < pfi->ft_size; ++index_ft_2) {              
                  pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z] *= pfi->f_sky;
                  pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][index_ft_2*pfi->ff_size+C][index_ft_1*pfi->ff_size+Z]
                    = pfi->fisher_matrix_CZ_smallest_largest[index_l1][index_l3][index_ft_1*pfi->ff_size+Z][index_ft_2*pfi->ff_size+C];
    }}}}}}}
    
  } // end of if(include_lensing_effects)


  return _SUCCESS_;
  
}


/** 
 * 
 */
int fisher_cross_cls (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi
        )
{

  /* The cross-power spectrum is defined as the symmetrix matrix
  C_l^ij = ( C_l^XX C_l^XY )
           ( C_l^YX C_l^YY )
  Its inverse, which we compute below, is needed to compute the covariance matrix
  (see eq. 17 of http://uk.arxiv.org/abs/astro-ph/0701921).
  
  When only one probe is considered (eg. T or E), the cross-power spectrum collapses
  to one number. */
  
  // ====================================================================================
  // =                                Allocate memory                                   =
  // ====================================================================================

  /* Allocate memory for the cross-power spectrum and its inverse. We tabulate them
  for all values of l between l_min and l_max, using interpolation. */
  class_alloc (pfi->C, pbi->full_l_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->inverse_C, pbi->full_l_size*sizeof(double **), pfi->error_message);

  for (int l=pbi->l_min; l <= pbi->l_max; ++l) {

    class_alloc (pfi->C[l-2], pfi->ff_size*sizeof(double*), pfi->error_message);
    class_alloc (pfi->inverse_C[l-2], pfi->ff_size*sizeof(double*), pfi->error_message);

    for (int index_ff=0; index_ff < pfi->ff_size; ++index_ff) {
      class_alloc (pfi->C[l-2][index_ff], pfi->ff_size*sizeof(double), pfi->error_message);
      class_alloc (pfi->inverse_C[l-2][index_ff], pfi->ff_size*sizeof(double), pfi->error_message);
    }
  }
  
  

  // ====================================================================================
  // =                              Lensed or not lensed?                               =
  // ====================================================================================
  
  /* The C_l's should be the observed ones, i.e. the LENSED ones (see for example
  Eq. 5.2 of Lewis and Challinor, 2011). */

  double ** cls = pbi->cls;
  int index_cls[pfi->ff_size][pfi->ff_size];
  for (int X=0; X < pfi->ff_size; ++X)
    for (int Y=0; Y < pfi->ff_size; ++Y)
      index_cls[X][Y] = pfi->index_ct_of_ff_ff[X][Y];

  if (pfi->include_lensing_effects == _TRUE_) {
  
    cls = pbi->lensed_cls;
    for (int X=0; X < pfi->ff_size; ++X)
      for (int Y=0; Y < pfi->ff_size; ++Y)
        index_cls[X][Y] = pfi->index_lt_of_ff_ff[X][Y];
  }



  // ====================================================================================
  // =                            Compute cross-power spectrum                          =
  // ====================================================================================

  /* The cross-power spectrum is used to build the covariance matrix and, therefore, needs to be
  include the instrumental noise. We assume that this affects only the direct C_l's,
  eg. C_l^TT and C_l^EE have noise but C_l^TE does not. */
  
  for (int l=pbi->l_min; l <= pbi->l_max; ++l) {  
    for (int X = 0; X < pfi->ff_size; ++X) {
      for (int Y = 0; Y < pfi->ff_size; ++Y) {

        /* Compute ij element */
        pfi->C[l-2][X][Y] = cls[index_cls[X][Y]][l-2];

        /* Add noise only for diagonal combinations (TT, EE). We assume a vanishing noise
        cross-correlation between T and E, as the temperature and polarization detectors
        in a given experiment are usually built independently. */
        if (X==Y) pfi->C[l-2][X][Y] += pfi->N_l[X][l-2];

        /* If a C_l vanishes, then we risk to have infinities in the inverse */
        class_test (fabs(pfi->C[l-2][X][Y]) < _MINUSCULE_, pfi->error_message,
          "C_%d was found ~ zero. Stopping to prevent seg fault. k-max too small?",l);
      }
    }
    
    /* Compute the inverse of C_l^ij for the considered l */
    InverseMatrix (pfi->C[l-2], pfi->ff_size, pfi->inverse_C[l-2]);
    
    /* Debug - print the Cl's */
    // double factor = 1e12*pba->T_cmb*pba->T_cmb*l*(l+1)/(2*_PI_);
    // fprintf (stderr, "%8d %10g %10g %10g %10g\n", l,
    //   factor*pfi->C[l-2][0][0], factor*pfi->C[l-2][1][1],
    //   factor*pfi->C[l-2][0][1], factor*pfi->C[l-2][1][0]);

    /* Debug - print the inverse Cl's */
    // double factor = 1e12*pba->T_cmb*pba->T_cmb*l*(l+1)/(2*_PI_);
    // fprintf (stderr, "%8d %10g %10g %10g %10g\n", l,
    //   pfi->inverse_C[l-2][0][0]/factor, pfi->inverse_C[l-2][1][1]/factor,
    //   pfi->inverse_C[l-2][0][1]/factor, pfi->inverse_C[l-2][1][0]/factor);

  } // end of for(l)
  
  
  return _SUCCESS_;
  
}




/**
 * Add the lensing variance to the Fisher matrix estimator, as explained in Section 5 of
 * Lewis, Challinor & Hanson (2011, http://uk.arxiv.org/abs/1101.2234). This added variance
 * comes from the fact that, in presence of lensing, the assumption of almost Gaussian CMB
 * assumed to obtain the standard form of the covariance matrix fails.
 *
 * In the function we code Eq. 5.35 of Lewis et al., which reduces to Eq. 5.25 in the presence
 * of only one bispectrum type. With respect to that reference we have l1->l3, i->C, j->Z.
 * We also factor out C_l^{X\phi} to avoid inverting singular matrices, as suggested in the
 * paper. Extra credits to Antony also because this function is inspired by the
 * SeparableBispectrum.F90 module of CAMB.
 *
 * This function fills the array 'pfi->fisher_matrix_lensvar_smallest'.
 */
int fisher_lensing_variance (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        struct fisher_workspace * pwf
        )
{

  /* Dimension of the matrices that we will need to invert. All the arrays will have
  the field (T,E...) and bispectrum types (local template, intrinsic, lensing...) 
  dimensions flattened in a single one. This mixing is explained in the last page
  of Sec. 5 in arxiv:1101.2234  */
  int N = pfi->ft_size * pfi->ff_size;

  /* Initialise the output arrays */
  for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {
    for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
      for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {
        pfi->fisher_matrix_lensvar_smallest[index_l1][index_ft_1][index_ft_2] = 0;
        pfi->fisher_matrix_lensvar_lmin[index_l1][index_ft_1][index_ft_2] = 0;
  }}}

  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #ifdef _OPENMP
  #pragma omp parallel private (thread)
  number_of_threads = omp_get_num_threads();
  #endif

  /* Temporary arrays */
  double *** inverse_f, *** f;
  class_alloc (inverse_f, number_of_threads*sizeof(double**), pfi->error_message);
  class_alloc (f, number_of_threads*sizeof(double**), pfi->error_message);

  #pragma omp parallel shared (abort) private (thread)
  { 
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    class_alloc_parallel (inverse_f[thread], N*sizeof(double *), pfi->error_message);
    class_alloc_parallel (f[thread], N*sizeof(double *), pfi->error_message);
    for (int i=0; i < N; ++i) {
      class_calloc_parallel (inverse_f[thread][i], N, sizeof(double), pfi->error_message);
      class_calloc_parallel (f[thread][i], N, sizeof(double), pfi->error_message);
    }

    /* Start sum over the smallest multipole: l1<=l2<=l3 */
    #pragma omp for schedule (dynamic)
    for (int index_l1=0; index_l1 < pfi->l1_size; ++index_l1) {

      /* Shortcuts to F_bar(l1,C,Z) and F_bar^-1(l1,C,Z) */
      double ** f_bar = f[thread];
      double ** inverse_f_bar = inverse_f[thread];

      /* By default, assume there is no lensing variance. This instruction is needed
      in order for the 'goto' statements to work */
      for (int i=0; i < N; ++i)
        for (int j=0; j < N; ++j)
          f_bar[i][j] = pfi->fisher_matrix_CZ_smallest[index_l1][i][j];
    
      /* Determine current l1 value */
      int l1 = pfi->l1[index_l1];
  
      // ----------------------------------------------------------------------------------
      // -                                Perform checks                                  -
      // ----------------------------------------------------------------------------------
  
      /* Skip l1 contributions where the Fisher matrix has not been computed */
      if ((l1 < pfi->l1_min_global) || (l1 > pfi->l1_max_global))
        goto lensing_variance_final_step;

      /* Limit the computation of the lensing variance only to those l's where we can trust
      the linear C_l^{\phi\phi}. Not doing so results in ~1% differences in the Fisher
      matrix elements of the local and intrinsic bispectra. */
      int l1_max = 0;
    
      if (pfi->has_fisher_t == _TRUE_)
        l1_max = MAX (l1_max, pbi->lmax_lensing_corrT);
      if (pfi->has_fisher_e == _TRUE_)
        l1_max = MAX (l1_max, pbi->lmax_lensing_corrE);
    
      if (l1>l1_max)
        goto lensing_variance_final_step;

      /* Skip l1 contributions where the Fisher matrix vanishes. This might happen, for example,
      when l1 is very close to the maximum l, since l1 is the smallest multipole in the sum. */
      int skip = _TRUE_;
      for (int i=0; i < N; ++i)
        if (fabs(f_bar[i][i]) > _MINUSCULE_)
          skip = _FALSE_;

      if (skip == _TRUE_) {
        printf_log_if (pfi->fisher_verbose, 2, 
          "     * skipping the contribution from l1=%d to the lensing variance (Fisher matrix too small)\n", l1);        
        goto lensing_variance_final_step;
      }
        
      /* Check that the matrix is not singular. If this is the case, pretend there is no lensing
      variance for this value of l1, and jump to the end of the l1 loop, where we rescale the
      Fisher matrix. This can happen if, for example, for this value of (l1,C,Z):
      1) a bispectrum has underflown for certain values of l1 (row of zeros)
      2) a bispectrum has been artificially set to zero for this value of (l1,C,Z) (row of zeros)
      3) two of the bispectra included in the Fisher matrix analysis are too similar, e.g.,
      because you have asked for lensing=yes and for cmb-lensing at the same time (two equal rows) */
      double det = Determinant(f_bar, N);
      if (fabs(det) < _MINUSCULE_) {

        if (pfi->fisher_verbose > 2) {
          printf_log ("     * skipping the contribution from l1=%d to the lensing variance (det=%g too small)\n", l1, det);
          printf_log ("     * F_bar matrix:\n");
          PrintMatrix (f_bar, N);
        }

        goto lensing_variance_final_step;
      }
    
        
      // ------------------------------------------------------------------------------------
      // -                             Compute l1,C,Z contribution                          -
      // ------------------------------------------------------------------------------------

      printf_log_if (pfi->fisher_verbose, 2, 
        "     * processing l1=%d...\n", l1);

      /* Contribution to the considered l1 (smallest multipole in the sum) to the Fisher
      matrix, including lensing variance. This quantity will be summed over C and Z. */
      double l1_contribution = 0;
  
      /* Invert the Fisher matrix with indices Z and C (i and p in Eq. 5.25 of 1101.2234) */
      InverseMatrix (f_bar, N, inverse_f_bar);
    
      /* Debug - print the matrix and its inverse */
      // if ((l1==1000) || (l1==986)) {
      //   printf ("F_bar(l1) for l1=%d, determinant = %g\n", l1, det);
      //   PrintMatrix (pfi->fisher_matrix_CZ_smallest[index_l1], N);
      //   printf ("F_bar^{-1}(l1) for l1=%d:\n", l1);
      //   PrintMatrix (inverse_f_bar, N);
      // }
  
      /* Find the lensing-induced noise contribution to the Fisher matrix coming from the considered
      multipole. With respect to Eq. 5.35 of Lewis et al. 2011, this is the term in square brackets.
      It should be noted that there are regions in l-space where the CMB-lensing bispectrum vanishes,
      but the effect of lensing on the Fisher matrix doesn't. In Eq. 5.35, this corresponds to the case
      when the first term in square brackets (C_ZP*C_CP) vanishes while the second one (C_tot_CZ*C_PP)
      doesn't. This might happen for l>100 where C_l^{T\phi} is very small. See also comment to that
      equation in the paper.  */
      for (int C=0; C < pfi->ff_size; ++C) {
        for (int Z=0; Z < pfi->ff_size; ++Z) {
            
          /* Compute the extra noise due to lensing. The cross-correlation between the two fields
          (C_tot_CZ) should include the instrument noise, as specified below Eq. 5.1 of Lewis
          et al. 2011. The second term (C_tot_CZ*C_PP) dominates over the first one (C_ZP*C_CP).
          Note that the noise_correction with respect to Eq. 5.35 (ibidem) is rescaled by a factor
          C_ZP*C_CP in order to avoid numerical issues (we cannot trust neither C_ZP nor C_CP
          for l>~100) */
          double C_tot_CZ = pfi->C[l1-2][C][Z];
          double C_CP = pbi->cls[pfi->index_ct_of_phi_ff[ C ]][l1-2];
          double C_ZP = pbi->cls[pfi->index_ct_of_phi_ff[ Z ]][l1-2];
          double C_PP = pbi->cls[psp->index_ct_pp][l1-2];
          double noise_correction = (C_ZP*C_CP + C_tot_CZ*C_PP)/(2*l1+1);
  
          /* Contribution to the inverse Fisher matrix from this (l1,Z,C) */
          int index_C = pfi->index_ft_cmb_lensing_kernel*pfi->ff_size+C;
          int index_Z = pfi->index_ft_cmb_lensing_kernel*pfi->ff_size+Z;
          inverse_f_bar[index_C][index_Z] += noise_correction;

          /* Debug - print out the correction to F_bar(l1,C,Z)^-1. To plot the file in gnuplot, do:
          set term wxt 1; set log; plot [2:100] "noise_correction.dat" u 1:(abs($5)) every\
          4::0 w li t "noise TT", "" u 1:(abs($5)) every 4::1 w li t "noise TE", "" u 1:(abs($5))\
          every 4::3 w li t "noise EE” */
          // if (l1<200)
          //   fprintf (stderr, "%4d %4s %4s %17g %17g %17g %17g %17g %17g %17g\n",
          //     l1, pfi->ff_labels[C], pfi->ff_labels[Z],
          //     noise_correction,
          //     C_ZP*C_CP/(2*l1+1),
          //     C_tot_CZ*C_PP/(2*l1+1),
          //     inverse_f_bar[index_C][index_Z],
          //     (pbi->has_intrinsic==_TRUE_)?
          //     C_ZP*C_CP*inverse_f_bar[pfi->index_ft_intrinsic*pfi->ff_size+C]
          //     [pfi->index_ft_intrinsic*pfi->ff_size+Z]:0,
          //     (pbi->has_intrinsic_squeezed==_TRUE_)?
          //     C_ZP*C_CP*inverse_f_bar[pfi->index_ft_intrinsic_squeezed*pfi->ff_size+C]
          //     [pfi->index_ft_intrinsic_squeezed*pfi->ff_size+Z]:0,
          //     (pbi->has_local_model==_TRUE_)?
          //     C_ZP*C_CP*inverse_f_bar[pfi->index_ft_local*pfi->ff_size+C]
          //     [pfi->index_ft_local*pfi->ff_size+Z]:0
          //   );

          /* Debug - print r_l. This should match Fig. 3 of Lewis et al. 2011 when considering only
          temperature or only polarisation. Note that the r-plot thus produced will look as the
          absolute value of the curves in Fig. 3, as we compute it as 1/sqrt(r^-2). */
          // double r_minus_2 = C_tot_CZ*C_PP/(C_ZP*C_CP);
          // double r = 1/sqrt(r_minus_2);
          // double l_factor = l1*(l1+1)/(2*_PI_);
          // double t_factor = 2.7255*1e6;
          // double factor = pow(t_factor,2) * l_factor;
          // fprintf (stderr, "%4d %17g %17g %17g %17g %17g\n",
          //   l1, r, factor*C_tot_CZ, l_factor*l1*(l1+1)*C_PP,
          //   t_factor*l_factor*sqrt(l1*(l1+1))*C_ZP, t_factor*l_factor*sqrt(l1*(l1+1))*C_CP);
  
        } // Z
      } // C
    
      /* Debug - print the inverse matrix after adding the noise correction */
      // if ((l1==1000) || (l1==986)) {
      //   printf ("F^{-1}(l1) for l1=%d:\n", l1);
      //   PrintMatrix (inverse_f_bar, N);
      // }
    
      /* Check that the matrix is not singular */
      det = Determinant(inverse_f_bar, N);
      if (fabs(det) < _MINUSCULE_) {

        if (pfi->fisher_verbose > 2) {
          printf_log ("     * skipping the contribution from l1=%d to the lensing variance (det=%g too small)\n", l1, det);
          printf_log ("     * F_bar^{-1} + noise matrix:\n");
          PrintMatrix (inverse_f_bar, N);
        }

        goto lensing_variance_final_step;
      }

      /* Invert the contribution from (l1,Z,C) and overwrite f_bar */
      InverseMatrix (inverse_f_bar, N, f_bar);

      /* Check that the inversion didn't produce nans */
      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1)
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2)
          for (int C=0; C < pfi->ff_size; ++C)
            for (int Z=0; Z < pfi->ff_size; ++Z)
              class_test_parallel (isnan (f_bar[index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z]),
                pfi->error_message,
                "stopping to prevent nans");

    
      // ---------------------------------------------------------------------------------------------
      // -                                   Add up the contributions                                -
      // ---------------------------------------------------------------------------------------------

      /* Sum over (Z,C) to obtain the optimal Fisher matrix */
      lensing_variance_final_step:

      for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
        for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {
      
          for (int C=0; C < pfi->ff_size; ++C) {
            for (int Z=0; Z < pfi->ff_size; ++Z) {

              /* So far we have been using the kernel of the CMB-lensing bispectrum (Eq. 5.20 of
              arXiv:1101.2234). The reason for doing so is to avoid numerical noise in the term
              C_l1^{i\psi}*C_l1^{j\psi} at the denominator of Eq. 5.35. Here we rescale the Fisher
              matrix so that it contains the actual CMB-lensing bispectrum by multiplying it
              by C_l1^{i\psi}*C_l1^{j\psi}. */
              double correction = 1;
              if (index_ft_1 == pfi->index_ft_cmb_lensing_kernel)
                correction *= pbi->cls[pfi->index_ct_of_phi_ff[ C ]][l1-2];
              if (index_ft_2 == pfi->index_ft_cmb_lensing_kernel)
                correction *= pbi->cls[pfi->index_ct_of_phi_ff[ Z ]][l1-2];

              #pragma omp atomic
              pfi->fisher_matrix_lensvar_smallest[index_l1][index_ft_1][index_ft_2]
                += correction * f_bar[index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z];
    
              /* Rescale also the zero-signal matrix */
              #pragma omp atomic
              pfi->fisher_matrix_CZ_smallest[index_l1][index_ft_1*pfi->ff_size+C][index_ft_2*pfi->ff_size+Z]
                *= correction;
            
            } // Z
          } // C
        } // ft_2
      } // ft_1    
    } // index_l1
    
    /* Free memory */
    for (int X=0; X < N; ++X) {
      free (inverse_f[thread][X]);
      free (f[thread][X]);
    }
    free (inverse_f[thread]);
    free (f[thread]);
    
  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region

  free (inverse_f);
  free (f);


  // -----------------------------------------------------------------------------------------
  // -                      Lensing variance as a function of l_min                          -
  // -----------------------------------------------------------------------------------------  

  /* Find the total Fisher matrix by summing 'fisher_matrix_lensvar_smallest' over l1 */
  for (int index_ft_1=0; index_ft_1 < pfi->ft_size; ++index_ft_1) {
    for (int index_ft_2=0; index_ft_2 < pfi->ft_size; ++index_ft_2) {

      double accumulator = 0;

      for (int index_l1=(pfi->l1_size-1); index_l1>=0; --index_l1) {

        double l1_contribution = pfi->fisher_matrix_lensvar_smallest[index_l1][index_ft_1][index_ft_2];
        if (!pfi->l1_is_full)
          l1_contribution *= pwf->delta_l[index_l1];
        accumulator += l1_contribution;
        pfi->fisher_matrix_lensvar_lmin[index_l1][index_ft_1][index_ft_2] = accumulator;
        if (index_ft_1 == index_ft_2)
          pfi->sigma_fnl_lensvar_lmin[index_l1][index_ft_1]
            = 1/sqrt(pfi->fisher_matrix_lensvar_lmin[index_l1][index_ft_1][index_ft_1]);

      } // end of for(index_l1)  
    } // end of for(ft_2)
  } // end of for (ft_1)

 for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1)
   InverseMatrix (pfi->fisher_matrix_lensvar_lmin[index_l1],
     pfi->ft_size, pfi->inverse_fisher_matrix_lensvar_lmin[index_l1]);

  return _SUCCESS_;
  
}



