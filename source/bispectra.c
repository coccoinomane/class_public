/** @file bispectra.c documented spectra module for second-order perturbations
 *
 * Guido W Pettinari, 19.07.2012.
 */

#include "bispectra.h"


/**
 * This routine initializes the spectra structure (in particular, 
 * computes table of anisotropy and Fourier spectra \f$ C_l^{X}, P(k), ... \f$)
 * 
 * @param ppr Input : pointer to precision structure
 * @param pba Input : pointer to background structure (will provide H, Omega_m at redshift of interest)
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfer structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Output: pointer to initialized spectra structure
 * @return the error status
 */

int bispectra_init (
     struct precision * ppr,
     struct background * pba,
     struct thermo * pth,
     struct perturbs * ppt,
     struct bessels * pbs,
     struct transfers * ptr,
     struct primordial * ppm,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi
     )
{

  /* Check whether we need to compute bispectra at all */  
  if (pbi->has_bispectra == _FALSE_) {

    if (pbi->bispectra_verbose > 0)
      printf("No bispectra requested. Bispectra module skipped.\n");

    pbi->bt_size=0;

    return _SUCCESS_;
  }
  else {
    if (pbi->bispectra_verbose > 0)
      printf("Computing bispectra\n");
  }


  // =========================================================
  // =                    Preparations                       =
  // =========================================================

  /* Initialize indices & arrays in the bispectra structure */

  class_call (bispectra_indices (ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);





  // =======================================================
  // =                  Compute bispectra                  =
  // =======================================================
  
  class_call (bispectra_harmonic (ppr,pba,ppt,pbs,ptr,ppm,psp,pbi),
    pbi->error_message,
    pbi->error_message);



  // ============================================================================
  // =                      Check bispectra against nan's                       =
  // ============================================================================

  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    for (int i = 0; i < pbi->bf_size; ++i) {
    for (int j = 0; j < pbi->bf_size; ++j) {
    for (int k = 0; k < pbi->bf_size; ++k) {
    
      for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
        for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
       
          /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
          int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
          int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
       
          for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

            /* Index of the current (l1,l2,l3) configuration */
            long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
            
            double bispectrum = pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3];
     
            if (isnan(bispectrum))
              printf ("@@@ WARNING: b(%d,%d,%d) = %g for bispectrum '%s_%s'.\n",
              pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], bispectrum,
              pbi->bt_labels[index_bt], pbi->bfff_labels[i][j][k]);
  
          } // end of for(index_l3)
        } // end of for(index_l2)
      } // end of for(index_l1)
    }}} // end of for(ijk)
  } // end of for(index_bt)


	/* Debug - Print some bispectra configuration */
  // for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {
  //   for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
  //     for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  //       int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
  //       int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  //       for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
  //         long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  //         int l1 = pbi->l[index_l1]; int l2 = pbi->l[index_l2]; int l3 = pbi->l[index_l3];
  //         double cl_1 = pbi->cls[psp->index_ct_tt][l1-2];;
  //         double cl_2 = pbi->cls[psp->index_ct_tt][l2-2];;
  //         double cl_3 = pbi->cls[psp->index_ct_tt][l3-2];;
  //         double squeezed_normalisation = 1/(12*cl_1*cl_3);
  //         double equilateral_normalisation = 1e16 * l1*l1 * (l1+1)*(l1+1) / ((2*_PI_)*(2*_PI_));
  //         if (strstr(pbi->bt_labels[index_bt], "equilateral")!=NULL) {
  //           /* Squeezed configuration */
  //           if ((l3==6) && (l1==l2)) {
  //             fprintf (stderr, "%12d %12g %12g %12g %12g %12g %12g %12g %12g\n",
  //               l1,
  //               pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][1][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //               pbi->bispectra[index_bt][1][1][0][index_l1_l2_l3] * equilateral_normalisation
  //             );
  //           }
  //           /* Equilateral configuration */
  //           // if ((l1==l2) && (l2==l3)) {
  //           //   fprintf (stderr, "%12d %12g %12g %12g\n",
  //           //   l1, pbi->bispectra[index_bt][index_l1_l2_l3], equilateral_normalisation,
  //           //   pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3]*equilateral_normalisation);
  //           // }
  //         }
  //       } // end of for(index_l3)
  //     } // end of for(index_l2)
  //   } // end of for(index_l1)
  // } // end of for(index_bt)

  return _SUCCESS_;

}






/**
 * This routine frees all the memory space allocated by bispectra_init().
 *
 * @param pbi Input: pointer to bispectra structure (which fields must be freed)
 * @return the error status
 */

int bispectra_free(
     struct perturbs * ppt,
     struct spectra * psp,
     struct bispectra * pbi
     )
{

  if (pbi->has_bispectra == _TRUE_) {

    free(pbi->l);
    free(pbi->pk);

    for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {

      free(pbi->l_triangular_size[index_l1]);
      free(pbi->index_l_triangular_min[index_l1]);
      free(pbi->index_l_triangular_max[index_l1]);
    
    } // end of for(index_l1)

    free(pbi->l_triangular_size);
    free(pbi->index_l_triangular_min);
    free(pbi->index_l_triangular_max);


    /* Free pbi->bispectra */
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      for (int i = 0; i < pbi->bf_size; ++i) {
        for (int j = 0; j < pbi->bf_size; ++j) {
          for (int k = 0; k < pbi->bf_size; ++k)
            free (pbi->bispectra[index_bt][i][j][k]);
          free (pbi->bispectra[index_bt][i][j]);
        }
        free (pbi->bispectra[index_bt][i]);
      }
      free (pbi->bispectra[index_bt]);
    }
    free (pbi->bispectra);


    /* Arrays specific to the primordial models */
  
    if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
        
      for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
        for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
          free (pbi->alpha[index_bf][index_l]);
          free (pbi->beta[index_bf][index_l]);
        }
        free (pbi->alpha[index_bf]);
        free (pbi->beta[index_bf]);
      }
      free (pbi->alpha);
      free (pbi->beta);
    } // end of local model

    if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
  
      for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
        for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
          free (pbi->gamma[index_bf][index_l]);
          free (pbi->delta[index_bf][index_l]);
        }
        free (pbi->gamma[index_bf]);
        free (pbi->delta[index_bf]);
      }
      free (pbi->gamma);
      free (pbi->delta);
    } // end of equilateral and orthogonal models  
    
    free (pbi->delta_k);
    
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->cls[index_ct]);
    free (pbi->cls);
    // if (pbi->has_intrinsic_squeezed == _TRUE_) {
    //   for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    //     free (pbi->dlog_cls[index_ct]);
    //   free (pbi->dlog_cls);
    // }
      
    
  } // end of if(has_bispectra)

  
  return _SUCCESS_;
 
}








/**
 * This routine defines indices and allocates tables in the bispectra structure 
 *
 */

int bispectra_indices (
        struct precision * ppr,
        struct background * pba,
        struct perturbs * ppt,
        struct bessels * pbs,
        struct transfers * ptr,
        struct primordial * ppm,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi
        )
{ 

  // ================================================================================================
  // =                                        Count bispectra                                       =
  // ================================================================================================

  /* Find out which kind of bispectra to compute and assign them indices and labels */


  // ------------------------------------------------------------------
  // -                        Bispectra fields                        -
  // ------------------------------------------------------------------  

  /* Generate indices for the probes (T for temperature, E for E-mode polarisation,
  R for Rayleigh...) that we will use to build the bispectra and the Fisher matrix
  elements. */
  
  int index_bf = 0;
  
  if (ppt->has_cl_cmb_temperature == _TRUE_) {
    pbi->has_bispectra_t = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "t");
    pbi->index_bf_t = index_bf++;
  }
  
  if (ppt->has_cl_cmb_polarization == _TRUE_) {
    pbi->has_bispectra_e = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "e");
    pbi->index_bf_e = index_bf++;
  }

  pbi->bf_size = index_bf;
  pbi->n_probes = pow(pbi->bf_size, 3);
      
  class_test (pbi->bf_size > _MAX_NUM_FIELDS_,
    pbi->error_message,
    "we cannot compute the bispectrum for more than two probes (e.g. T and E), reduce your expectations :-)");

  class_test (pbi->bf_size < 1,
    pbi->error_message,
    "no probes requested");

  /* Create labels for the full bispectra */
  for (int i = 0; i < pbi->bf_size; ++i)
  for (int j = 0; j < pbi->bf_size; ++j)
  for (int k = 0; k < pbi->bf_size; ++k)
    sprintf (pbi->bfff_labels[i][j][k], "%s%s%s",
    pbi->bf_labels[i], pbi->bf_labels[j], pbi->bf_labels[k]);


  /* Associate to each field (T,E,...) its transfer function, which was computed in the transfer.c module,
  and to each possible pair of fields (TT,EE,TE,...) their power spectra, which were computed
  in the spectra.c module. */
  for (int i = 0; i < pbi->bf_size; ++i) {
    
    if ((pbi->has_bispectra_t == _TRUE_) && (i == pbi->index_bf_t)) {
      pbi->index_tt_of_bf[i] = ptr->index_tt_t;
      pbi->index_ct_of_bf[i][i] = psp->index_ct_tt;
    }

    if ((pbi->has_bispectra_e == _TRUE_) && (i == pbi->index_bf_e)) {
      pbi->index_tt_of_bf[i] = ptr->index_tt_e;
      pbi->index_ct_of_bf[i][i] = psp->index_ct_ee; 
    }

    for (int j = 0; j < pbi->bf_size; ++j) {

      if (((pbi->has_bispectra_t == _TRUE_) && (i == pbi->index_bf_t))
       && ((pbi->has_bispectra_e == _TRUE_) && (j == pbi->index_bf_e)))
        pbi->index_ct_of_bf[i][j] = pbi->index_ct_of_bf[j][i] = psp->index_ct_te;
    }
  }
  

  // /* The bte bispectrum is just implemented as a test, ignore */
  // if (ppt->has_cmb_polarization2 == _TRUE_) {
  // 
  //   pbi->has_bte = _TRUE_;
  //   strcpy (pbi->bp_labels[index_bp], "bte");
  //   pbi->bispectrum_multiplicity[index_bp] = 1;
  //   pbi->index_bp_bte = index_bp++;
  // 
  //   /* DIRTYNESS!!! Reset the way the probes (ttt,bte...) interact with the bispectra (primordial, intrinsic) */
  //   pbi->has_ttt = _FALSE_;  
  // }
  // 
  // pbi->bp_size = index_bp;


  // ------------------------------------------------------------------
  // -                        Bispectra types                         -
  // ------------------------------------------------------------------  
  
  int index_bt = 0;
  pbi->bt_size = 0;
  
  pbi->n[separable_bispectrum] = 0;
  pbi->n[non_separable_bispectrum] = 0;
  pbi->n[analytical_bispectrum] = 0;
  pbi->n[intrinsic_bispectrum] = 0;

  // *** Separable bispectra
  
  if (pbi->has_local_model) {
    pbi->index_bt_local = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_equilateral_model) {
    pbi->index_bt_equilateral = index_bt;
    strcpy (pbi->bt_labels[index_bt], "equilateral");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_orthogonal_model) {
    pbi->index_bt_orthogonal = index_bt;
    strcpy (pbi->bt_labels[index_bt], "orthogonal");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    index_bt++;
  }

  // *** Non-separable bispectra

  if (pbi->has_galileon_model) {

    /* Bispectrum induced by pi_dot*pi_grad^2 */
    pbi->index_bt_galileon_gradient = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_grad");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    index_bt++;

    /* Bispectrum induced by pi_dot^3 */
    pbi->index_bt_galileon_time = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_time");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    index_bt++;
  }

  // *** Analytical bispectra

  if (pbi->has_local_squeezed == _TRUE_) {
    pbi->index_bt_local_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "l_squeezed");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  /* TODO: implemented intrinsic squeezed approximation */
  // if (pbi->has_intrinsic_squeezed == _TRUE_) {
  //   pbi->index_bt_intrinsic_squeezed = index_bt;
  //   strcpy (pbi->bt_labels[index_bt], "i_squeezed");
  //   pbi->bispectrum_type[index_bt] = analytical_bispectrum;
  //   pbi->n[analytical_bispectrum]++;
  //   index_bt++;
  // }

  if (pbi->has_cosine_shape == _TRUE_) {
    pbi->index_bt_cosine = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cosine");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  if (pbi->has_isw_lensing == _TRUE_) {
    pbi->index_bt_isw_lensing = index_bt;
    strcpy (pbi->bt_labels[index_bt], "isw-lensing");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    index_bt++;
  }

  // *** Intrinsic (i.e. second-order) bispectra

  if ((pbi->has_intrinsic) && (ppt->has_cmb_temperature_2==_TRUE_)) {
    pbi->index_bt_intrinsic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic");
    pbi->bispectrum_type[index_bt] = intrinsic_bispectrum;
    pbi->n[intrinsic_bispectrum]++;
    index_bt++;
  }
    
  pbi->bt_size = index_bt;


  /* Perform some checks */
  // class_test (
  //   (pbi->has_intrinsic_squeezed==_TRUE_) && ((pbi->bf_size>1) || ((pbi->bf_size==1 && pbi->has_bispectra_t==_FALSE_))),
  //   pbi->error_message,
  //   "the squeezed bispectrum is implemented only for temperature.");

  class_test (
    (pbi->has_isw_lensing==_TRUE_) && ((pbi->bf_size>1) || ((pbi->bf_size==1 && pbi->has_bispectra_t==_FALSE_))),
    pbi->error_message,
    "the isw_lensing bispectrum is implemented only for temperature.");

  class_test (
    (pbi->has_intrinsic==_TRUE_) && ((pbi->bf_size>1) || ((pbi->bf_size==1 && pbi->has_bispectra_t==_FALSE_))),
    pbi->error_message,
    "the intrinsic bispectrum is implemented only for temperature.");



  // =================================================================================================
  // =                                      Determine l-sampling                                     =
  // =================================================================================================
  
  /* We compute the bispectrum on a mesh where l1>=l2>=l3, with l3 determined by the triangular
  condition. All the multipoles are drawn from pbi->l, which is a copy of ptr->l.  */

  pbi->l_size = ptr->l_size_max;
  class_alloc (pbi->l, pbi->l_size*sizeof(int), pbi->error_message);
  for(int index_l=0; index_l<pbi->l_size; ++index_l)
    pbi->l[index_l] = ptr->l[index_l];

  /* Some debug - Print the full l-list */
  // printf ("     * ");
  // for (int index_l=0; index_l < (pbi->l_size-1); ++index_l)
  //   printf ("%d,", pbi->l[index_l]);
  // printf ("%d\n", pbi->l[pbi->l_size-1]);


  /* Maximum value in pbi->l */
  pbi->l_max = pbi->l[pbi->l_size-1];
  pbi->full_l_size = pbi->l_max - 2 + 1;

  // *** Allocate & fill pb1->l_triangular_size and pb1->index_l_triangular_min
  
  /* The three multipole indexes l1, l2, l3 satisfy the triangular condition |l1-l2| <= l3 <= l1+l2.
  We choose to impose the condition on l3, for a fixed couple (l1,l2).  As a consequence, the
  allowed values for l3 will only be a subset of those contained in the vector pbs->l. The subset
  is determined by pbi->index_l_triangular_min[index_l1][index_l2] and pbi->l_triangular_size[index_l1][index_l2],
  which we fill below. */

  class_alloc (pbi->l_triangular_size, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l_triangular_min, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l_triangular_max, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->index_l1_l2_l3, pbi->l_size*sizeof(long int **), pbi->error_message);

  /* Number of multipole configurations (l1,l2,l3) that satisfy the triangular condition and l1>=l2>=l3 */
  pbi->n_independent_configurations = 0;

  /* Number of multipole configurations (l1,l2,l3) that satisfy the triangular condition */
  pbi->n_total_configurations = 0;

  /* Initialise the index for a given (l1,l2,l3) triplet */
  long int index_l1_l2_l3 = 0;

  for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {

    int l1 = pbi->l[index_l1];

    /* Allocate l1 level */
    class_alloc (pbi->l_triangular_size[index_l1], pbi->l_size*sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_min[index_l1], pbi->l_size*sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_max[index_l1], pbi->l_size*sizeof(int), pbi->error_message);

    /* We consider only configurations whereby l1>=l2 */
    class_alloc (pbi->index_l1_l2_l3[index_l1], (index_l1+1)*sizeof(long int *), pbi->error_message);

    /* Fill pbi->l_triangular_size and pbi->index_l_triangular_min */
    for(int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
      
      int l2 = pbi->l[index_l2];

      /* Limits imposed on l3 by the triangular condition */
      int l_triangular_min = abs(l1-l2);
      int l_triangular_max = l1+l2;

      /* Find the index corresponding to l_triangular_min inside pbi->l */
      class_test (l_triangular_min >= pbi->l[pbs->l_size-1],
        pbi->error_message,
        "could not fulfill triangular condition for l3 using the multipoles in pbi->l");
      
      int index_l_triangular_min = 0;
      while (pbi->l[index_l_triangular_min] < l_triangular_min) ++index_l_triangular_min;


      /* Find the index corresponding to l_triangular_max inside pbi->l */
      class_test (l_triangular_max <= pbi->l[0],
        pbi->error_message,
        "could not fulfill triangular condition for l3 using the multipoles in pbi->l");
      
      int index_l_triangular_max = pbi->l_size-1;
      while (pbi->l[index_l_triangular_max] > l_triangular_max) --index_l_triangular_max;
    

      /* Fill pbi->index_l_triangular_min and pbi->l_triangular_size */
      pbi->index_l_triangular_min[index_l1][index_l2] = index_l_triangular_min;
      pbi->index_l_triangular_max[index_l1][index_l2] = index_l_triangular_max;
      pbi->l_triangular_size[index_l1][index_l2] = index_l_triangular_max - index_l_triangular_min + 1;

      /* Update counter of triangular configurations */
      pbi->n_total_configurations += pbi->l_triangular_size[index_l1][index_l2];

      /* We shall store the bispectra only for those configurations that contemporaneously satisfy 
      the triangular condition and the l1>=l2>=l3 condition. We use the pbi->index_l1_l2_l3 array
      to keep track of the index assigned to a given allowed configuration. */
      if (index_l2<=index_l1) {

        /* When the triangular condition is not compatible with index_l3<=index_l2, then index_l3_max < index_l3_min+1 will
        be either zero or negative. */ 
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
        int l3_size = MAX (0, index_l3_max-index_l3_min+1);
        class_alloc (pbi->index_l1_l2_l3[index_l1][index_l1-index_l2], l3_size*sizeof(long int), pbi->error_message);
        
        /* The indexing of pbi->index_l1_l2_l3 reflects the l1>=l2>=l3 constraint */
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
          pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3] = index_l1_l2_l3++;

      }
      
      /* Print some debug */
      // printf("l1=%d, l2=%d, l_triangular_min[%d]=%d, l_triangular_max[%d]=%d, l[index_l_triangular_min]=%d, l=[index_l_triangular_max]=%d, l_triangular_size=%d\n", 
      //   l1, l2, index_l_triangular_min, l_triangular_min, index_l_triangular_max, l_triangular_max, pbs->l[index_l_triangular_min], pbs->l[index_l_triangular_max], pbi->l_triangular_size[index_l1][index_l2]);
      
    } // end of for(index_l2)
  } // end of for(index_l1)

  /* Each bispectrum will be computed for the following number of configurations */
  pbi->n_independent_configurations = index_l1_l2_l3;
  

  /* Inform the user on how much her machine will have to suffer */
  if (pbi->bispectra_verbose > 1) {
    printf(" -> we shall compute %dx%d=%d bispectr%s for %ld configurations of (l1,l2,l3)\n",
      pbi->bt_size, pbi->n_probes, pbi->bt_size*pbi->n_probes,
      ((pbi->bt_size*pbi->n_probes)!=1?"a":"um"), pbi->n_independent_configurations);
    // printf("    with (L1,L2,L3) ranging from l=%d to %d (l_size=%d)\n", pbi->l[0], pbi->l[pbi->l_size-1], pbi->l_size);
  }
  
  
  
  
  
  // ============================================================================================
  // =                               Integration grid in k1 and k2                              =
  // ============================================================================================
  
  
  /* Sampling of the first-order transfer functions  */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;


  /* Allocate & fill delta_k, which is needed for the trapezoidal integration of the bispectrum.
  This array is has ptr->q_size elements, and is defined as k(i+1)-k(i-1) except for the first
  k(1)-k(0) and last k(N)-k(N-1) elements.  Note that when dealing with non-separable shapes,
  we shall not use this array for the integration over k3, as in that case the grid in k3 is
  not fixed but it depends on k1 and k2. */
  class_alloc (pbi->delta_k, k_tr_size * sizeof(double), pbi->error_message);
  
  /* Fill pbi->delta_k */
  pbi->delta_k[0] = k_tr[1] - k_tr[0];
      
  for (int index_k=1; index_k < k_tr_size-1; ++index_k)
    pbi->delta_k[index_k] = k_tr[index_k+1] - k_tr[index_k-1];
      
  pbi->delta_k[k_tr_size-1] = k_tr[k_tr_size-1] - k_tr[k_tr_size-2];
  
  
  
  
  
  
  // ==============================================================================================
  // =                             Allocate memory for bispectra                                  =
  // ==============================================================================================
  
  class_alloc (pbi->bispectra, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
  for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {

    class_alloc (pbi->bispectra[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

    for (int i = 0; i < pbi->bf_size; ++i) {
      
      class_alloc (pbi->bispectra[index_bt][i], pbi->bf_size*sizeof(double **), pbi->error_message);

      for (int j = 0; j < pbi->bf_size; ++j) {
        
        class_alloc (pbi->bispectra[index_bt][i][j], pbi->bf_size*sizeof(double *), pbi->error_message);
        
        for (int k = 0; k < pbi->bf_size; ++k)
          class_alloc (pbi->bispectra[index_bt][i][j][k], pbi->n_independent_configurations*sizeof(double), pbi->error_message);
        
      }
    }
  }
  
  pbi->count_allocated_for_bispectra = pbi->bt_size*pbi->n_probes*pbi->n_independent_configurations;
  
  if (pbi->bispectra_verbose > 1)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the bispectra array\n",
      pbi->count_allocated_for_bispectra*sizeof(double)/1e6, pbi->count_allocated_for_bispectra);



  // ===========================================================================================
  // =                                     Interpolate P(k)                                    =
  // ===========================================================================================
  
  /* The primordial power spectrum for the comoving curvature perturbatio R was already
  computed when we initialized the ppm structure, but we store it locally for the k-values
  of ptr-k to access it faster */
  class_call (bispectra_primordial_power_spectrum(
                pba,
                ppt,
                ptr,
                ppm,
                pbi),
    pbi->error_message,
    pbi->error_message);





  // =============================================================================================
  // =                                   Interpolate the C_l's                                   =
  // =============================================================================================

  /* Interpolate the Cl's in all l-values */
  class_call (bispectra_cls(
                ppt,
                psp,
                ple,
                pbi),
    pbi->error_message,
    pbi->error_message);





  return _SUCCESS_;

}







/**
 * Compute and store in pbi->pk[index_k] the primordial power spectrum of the Newtonian
 * time potential psi. The power spectrum is computed in the points contained in ptr->q, as
 * these are the points where the transfer functions are computed.
 *
 * The purpose of this function is twofold. First, it stores the primordial power spectrum
 * into memory for faster access by the bispectra module, as calling primordial_spectrum_at_k
 * is fairly expensive. Secondly, we convert the dimensionless spectrum of the curvature
 * perturbation R outputted by the primordial module of CLASS, Delta_R(k), into the power
 * spectrum for the Newtonian curvature potential, P_Phi(k). The
 * two spectra are related by:
 *  
 *  P_Phi(k) = 2*Pi^2/k^3 * Delta_R(k)
 *
 * where
 * 
 * Delta_R(k) = A_s * (k/k_pivot)^(n_s-1)
 *
 */  
int bispectra_primordial_power_spectrum (
    struct background * pba,
    struct perturbs * ppt,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi
    )
{

  /* Allocate the pbi->pk vector so that it contains ptr->q_size values */
  int k_size = ptr->q_size;
  class_alloc (pbi->pk, k_size*sizeof(double), pbi->error_message);
  
  /* Fill pk with the values of the primordial power spectrum, as obtained in the ppm module */

  for (int index_k=0; index_k<k_size; ++index_k) {
    
    double k = ptr->q[index_k];

    class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k,
                  &(pbi->pk[index_k])),
      ppm->error_message,
      pbi->error_message);

    /* Convert CLASS dimensionless power spectrum for the curvature perturbation into the dimensional one. */
    pbi->pk[index_k] = 2*_PI_*_PI_/(k*k*k) * pbi->pk[index_k];
    
  } // end of for(index_k)
  
  

  /* Do the same, but with ppt->k */
  int k_pt_size = ppt->k_size;
  class_alloc (pbi->pk_pt, k_pt_size*sizeof(double), pbi->error_message);
  
  /* Fill pk with the values of the primordial power spectrum, as obtained in the ppm module */

  for (int index_k_pt=0; index_k_pt<k_pt_size; ++index_k_pt) {
    
    double k_pt = ppt->k[index_k_pt];

    class_call (primordial_spectrum_at_k (
                  ppm,
                  ppt->index_md_scalars,
                  linear,
                  k_pt,
                  &(pbi->pk_pt[index_k_pt])),
      ppm->error_message,
      pbi->error_message);

    pbi->pk_pt[index_k_pt] = 2*_PI_*_PI_/(k_pt*k_pt*k_pt) * pbi->pk_pt[index_k_pt];
    
  } // end of for(index_k_pt)
  
  
  return _SUCCESS_;
  
}



int bispectra_cls (
    struct perturbs * ppt,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi
    )
{
  

  // ==============================================================================================
  // =                                          Allocate arrays                                   =
  // ==============================================================================================

  /* Allocate the array that will contain the C_l's for all types and all l's. */
  class_alloc (pbi->cls, psp->ct_size*sizeof(double*), pbi->error_message);
  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    class_alloc (pbi->cls[index_ct], pbi->full_l_size*sizeof(double), pbi->error_message);
  
  /* Do the same for the logarithmic derivative of the C_l's */
  // if (pbi->has_intrinsic_squeezed == _TRUE_) {
  //   class_alloc (pbi->dlog_cls, psp->ct_size*sizeof(double*), pbi->error_message);
  //   for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
  //     class_alloc (pbi->dlog_cls[index_ct], pbi->full_l_size*sizeof(double), pbi->error_message);
  // }
  
  
  /* We shall now call the CLASS function 'spectra_cl_at_l'. This gives three outputs:
  the total Cl's for each probe (T, E, B...); the Cl's divided by probe and mode
  (scalar, vector, tensor); the Cl's divided by probe, mode, and initial condition
  (adiabatic, isocurvature...). We have to allocate these three arrays before being
  able to call 'spectra_cl_at_l'. We do so copying what is done in the function
  'output_total_cl_at_l' in the output module */
  double * cl;        /* cl_md_ic[index_ct] */
  double ** cl_md;    /* cl_md[index_mode][index_ct] */
  double ** cl_md_ic; /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */
  class_alloc(cl, psp->ct_size*sizeof(double), pbi->error_message);	
  class_alloc(cl_md_ic, psp->md_size*sizeof(double *), pbi->error_message);
  class_alloc(cl_md, psp->md_size*sizeof(double *), pbi->error_message);
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {
    if (psp->md_size > 1)
      class_alloc(cl_md[index_mode], psp->ct_size*sizeof(double), pbi->error_message);	
    if (psp->ic_size[index_mode] > 1)
      class_alloc(cl_md_ic[index_mode], psp->ic_ic_size[index_mode]*psp->ct_size*sizeof(double), pbi->error_message);
  }
  
  
  
  // ==========================================================================================================
  // =                                                Store C_l's                                             =
  // ==========================================================================================================
  
  for (int l=2; l<=pbi->l_max; ++l) {
    
    class_call(spectra_cl_at_l(
                 psp,
                 (double)l,
                 cl,
                 cl_md,
                 cl_md_ic),
      psp->error_message,
      pbi->error_message);
      
    /* Store the total Cl's into an array as a function of l and probe. By 'total' we mean the power
    spectrum summed over all the modes and initial conditions */
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      pbi->cls[index_ct][l-2] = cl[index_ct];

    /* To compute the squeezed limit approximation case, we need the logarithimc derivative of the
    C_l's. We interpolate exactly in the same way as we did for the normal C_l's. */
    // if (pbi->has_intrinsic_squeezed == _TRUE_) {
    // 
    //   class_call(spectra_dlogcl_at_l(
    //                psp,
    //                (double)l,
    //                cl,
    //                cl_md,
    //                cl_md_ic),
    //     psp->error_message,
    //     pbi->error_message);
    //   
    //   /* Store the total Cl's into an array as a function of l and probe. By 'total' we mean the power
    //   spectrum summed over all the modes and initial conditions */
    //   for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    //     pbi->dlog_cls[index_ct][l-2] = cl[index_ct];
    // 
    //   /* Some debug - print the cls */
    //   // printf ("%d %g %g %g\n", l,
    //   //   pbi->cls[psp->index_ct_tt][l-2],
    //   //   pbi->cls[psp->index_ct_tz][l-2],
    //   //   pbi->dlog_cls[psp->index_ct_tt][l-2]);
    // 
    // }
  } // end of for loop on the l's


  /* Free memory */
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {    
    if (psp->md_size > 1) free(cl_md[index_mode]);  
    if (psp->ic_size[index_mode] > 1) free(cl_md_ic[index_mode]);
  }  
  free(cl_md_ic);
  free(cl_md);
      
  return _SUCCESS_;
  
}



/**
 * This routine computes a table of values for all harmonic spectra C_l's,
 * given the transfer functions and primordial spectra.
 * 
 * @param ppt Input : pointer to perturbation structure
 * @param ptr Input : pointer to transfers structure
 * @param ppm Input : pointer to primordial structure
 * @param psp Input/Output: pointer to spectra structure 
 * @return the error status
 */

int bispectra_harmonic (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
    struct bispectra * pbi
    )
{


  /* Initialize counter for the values that will be memorised in pbi->bispectra */
  pbi->count_memorised_for_bispectra = 0;


  // =============================================================================================
  // =                             Compute separable bispectra                                   =
  // =============================================================================================

  /* We first obtain the bispectra whose primordial shape function is separable in (k1,k2,k3) because they
  are quick to compute. */
      
  if (pbi->n[separable_bispectrum] > 0) {

    struct bispectra_workspace_separable * pwb_sep;
    class_alloc (pwb_sep, sizeof(struct bispectra_workspace_separable), pbi->error_message);

    /* Compute the separable integrals, and integrate them together in r according to the chosen models */
    class_call (bispectra_separable_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  pwb_sep),
      pbi->error_message,
      pbi->error_message);
  

    /* Free the 'pwb_sep' workspace */
    class_call (bispectra_separable_workspace_free(
                  pbi,
                  pwb_sep),
      pbi->error_message,
      pbi->error_message);

    /* IMPORTANT: No symmetrization is needed, as the primary bispectrum is automatically symmetrised if we feed a
    primordial bispectrum B(k1,k2,k3) symmetrised in k, as we do. */

  }






  // ======================================================================================================
  // =                                    Compute analytic bispectra                                      =
  // ======================================================================================================
  
  
  /* Compute the bispectra obtained from simple analytical formulas, such as the lensing one */
  
  if (pbi->n[analytical_bispectrum] > 0) {
  
    class_call (bispectra_analytical_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  psp,
                  pbi),
      pbi->error_message,
      pbi->error_message);
      
  }
  
  
  
  
  // ===========================================================================================================
  // =                                   Compute non-separable bispectra                                       =
  // ===========================================================================================================
  
  
  if (pbi->n[non_separable_bispectrum] > 0) {
  
    struct bispectra_workspace_non_separable * pwb_nonsep;
    class_alloc (pwb_nonsep, sizeof(struct bispectra_workspace_non_separable), pbi->error_message);
  
    /* Compute the non-separable bispectra */
    class_call (bispectra_non_separable_init(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  pwb_nonsep),
      pbi->error_message,
      pbi->error_message);
  
  
    /* Free the 'pwb_nonsep' workspace */
    /* TODO: MEMORY CHECKSUM ERROR, FIX */
    // class_call (bispectra_non_separable_workspace_free(
    //               pbi,
    //               pwb_nonsep),
    //   pbi->error_message,
    //   pbi->error_message);
  
  }
  
  /* Print information on memory usage */
  if ((pbi->bispectra_verbose > 1) && (pbi->n[intrinsic_bispectrum]) < 1)
    printf(" -> memorised ~ %.3g MB (%ld doubles) in the bispectra array\n",
      pbi->count_memorised_for_bispectra*sizeof(double)/1e6, pbi->count_memorised_for_bispectra);
  
  
  /* Check that we correctly filled the bispectra array (but only if there are no
  intrinsic bispectra left to be computed)*/
  if (pbi->n[intrinsic_bispectrum] < 1)
    class_test_permissive (pbi->count_allocated_for_bispectra != pbi->count_memorised_for_bispectra,
      pbi->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!",
      pbi->count_allocated_for_bispectra, pbi->count_memorised_for_bispectra);
  
  return _SUCCESS_;

}








int bispectra_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{

  // ===================================================================================================
  // =                                    Prepare integration grid                                     =
  // ===================================================================================================
  
  /* We set the r-sampling as if it were a time sampling.  We do so because 'r' has the right dimensions, and
  it always appears in the argument of a Bessel function multiplying a wavemode, just as it was for conformal
  time in the line-of-sight integral.  This is the only place in the module where the background structure
  is accessed. */
  pwb->r_min = pbi->r_min;
  pwb->r_max = pbi->r_max;
  pwb->r_size = pbi->r_size;
    
  /* We decide to sample r linearly */
  class_alloc (pwb->r, pwb->r_size*sizeof(double), pbi->error_message);
  lin_space (pwb->r, pwb->r_min, pwb->r_max, pwb->r_size);
  
  
  /* Allocate & fill delta_r, the measure for the trapezoidal integration over r */
  class_alloc (pwb->delta_r, pwb->r_size * sizeof(double), pbi->error_message);

  /* Fill pwb->delta_r */
  pwb->delta_r[0] = pwb->r[1] - pwb->r[0];
      
  for (int index_r=1; index_r < pwb->r_size-1; ++index_r)
    pwb->delta_r[index_r] = pwb->r[index_r+1] - pwb->r[index_r-1];
      
  pwb->delta_r[pwb->r_size-1] = pwb->r[pwb->r_size-1] - pwb->r[pwb->r_size-2];






  // =====================================================================================================
  // =                                  Allocate memory for filter functions                             =
  // =====================================================================================================
  
  /* Cycle variables */
  int index_l;
  

  /* Allocate the separable integrals needed for each requested model of primordial non-Gaussianity.
  Refer to the header file for details on the models. The integrand arrays need to be written
  by different threads at the same time, hence we allocate one of them for each thread. */
      
  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel
  {
    #ifdef _OPENMP
    number_of_threads = omp_get_num_threads();
    #endif
  }


  // -------------------------------------------------------
  // -            Needed by all bispectra models           -
  // -------------------------------------------------------
  
  /* The alpha and beta functions need to be computed for all models */      
  if ((pbi->has_local_model == _TRUE_) ||
      (pbi->has_equilateral_model == _TRUE_) ||
      (pbi->has_orthogonal_model == _TRUE_)) {

    class_alloc (pbi->alpha, pbi->bf_size*sizeof(double **), pbi->error_message);
    class_alloc (pbi->beta, pbi->bf_size*sizeof(double **), pbi->error_message);

    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

      class_alloc (pbi->alpha[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
      class_alloc (pbi->beta[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
  
      /* Allocate r-level of the integral arrays */
      for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
        class_calloc (pbi->alpha[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
        class_calloc (pbi->beta[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
      }
    }
    
    /* Allocate memory for the integrand functions */
    class_alloc (pwb->alpha_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    class_alloc (pwb->beta_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    
    #pragma omp parallel private (thread)
    {
      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
      
      class_calloc_parallel (pwb->alpha_integrand[thread], ptr->q_size, sizeof(double), pbi->error_message);
      class_calloc_parallel (pwb->beta_integrand[thread], ptr->q_size, sizeof(double), pbi->error_message);      
    }
    
  } // end of if all models



  // ---------------------------------------------------------
  // -          Specific to equilateral & orthogonal         -
  // ---------------------------------------------------------

  if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {

    class_alloc (pbi->gamma, pbi->bf_size*sizeof(double **), pbi->error_message);
    class_alloc (pbi->delta, pbi->bf_size*sizeof(double **), pbi->error_message);

    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

      class_alloc (pbi->gamma[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
      class_alloc (pbi->delta[index_bf], pbi->l_size*sizeof(double *), pbi->error_message);
  
      /* Allocate r-level of the integral arrays */
      for (int index_l=0; index_l<pbi->l_size; ++index_l) {    
        class_calloc (pbi->gamma[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
        class_calloc (pbi->delta[index_bf][index_l], pwb->r_size, sizeof(double), pbi->error_message);
      }
    }
    
    /* Allocate memory for the integrand functions */
    class_alloc (pwb->gamma_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    class_alloc (pwb->delta_integrand, number_of_threads*sizeof(double *), pbi->error_message);
    
    #pragma omp parallel private (thread)
    {
      #ifdef _OPENMP
      thread = omp_get_thread_num();
      #endif
      
      class_calloc_parallel (pwb->gamma_integrand[thread], ptr->q_size, sizeof(double), pbi->error_message);
      class_calloc_parallel (pwb->delta_integrand[thread], ptr->q_size, sizeof(double), pbi->error_message);      
    }
    
  } // end of if equilateral || orthogonal

  return _SUCCESS_;
  
}




int bispectra_separable_filter_functions (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bf,
    struct bispectra_workspace_separable * pwb
    )
{
  
  
  /* We shall integrate over the same k-points where we computed the first-order transfer functions */
  int k_size = ptr->q_size;

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel                                             \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)                         \
    private (thread)
  {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    /* Cycle variables, allocate inside parallel loop */
    int index_l, index_k, index_r;
      
    #pragma omp for schedule (dynamic)
    for (int index_l = 0; index_l < pbi->l_size; ++index_l) {

      if (pbi->bispectra_verbose > 2)
        printf("      \\ computing filter functions for l=%d, index_l=%d\n", pbi->l[index_l], index_l);

      /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
      int tt_size = ptr->tt_size[ppt->index_md_scalars];
      int l_size = ptr->l_size[ppt->index_md_scalars];

      double * transfer = &(ptr->transfer
        [ppt->index_md_scalars]
        [((ppt->index_ic_ad * tt_size + pbi->index_tt_of_bf[index_bf]) * l_size + index_l) * k_size]);


      // ==============================================
      // =          Build integrand functions         =
      // ==============================================

      /* Build the integrand functions, leaving out the k^2 * j_l(k*r) bit and any other numerical factor. */
    
      // *** Local and equilateral model ***
      if ((pbi->has_local_model == _TRUE_) ||
          (pbi->has_orthogonal_model == _TRUE_) ||
          (pbi->has_equilateral_model == _TRUE_)) {

        for (int index_k=0; index_k < ptr->q_size; ++index_k) {
        
          double pk = pbi->pk[index_k];
          double tr = transfer[index_k];

          pwb->alpha_integrand[thread][index_k] = tr;
          pwb->beta_integrand[thread][index_k] = tr * pk;
          
          /* Print out the transfer function for a given field */
          // if (ptr->l[index_l]==284)
          //   if ((ppt->has_cl_cmb_polarization == _TRUE_) && (pbi->index_tt_of_bf[index_bf] == ptr->index_tt_e))
          //     fprintf (stderr, "%12g %12g\n", ptr->q[index_k], tr);
          
        }

      } // end of all models
    
      // *** Equilateral and orthogonal models ***
      if ((pbi->has_equilateral_model == _TRUE_) ||
          (pbi->has_orthogonal_model == _TRUE_)) {
     
        /* Here we basically copy eqs. 15-18 in Creminelli et al. 2006, keeping out the 2/pi factor and the
        Bessel function, and multiplying the rest by k (we use the dimensional power spectrum and we already
        factored out a k^2 factor) */
        for (int index_k=0; index_k < ptr->q_size; ++index_k) {      

          double pk_one_third = pow(pbi->pk[index_k], _ONE_THIRD_);
          double pk_two_thirds = pk_one_third*pk_one_third;
          double tr = transfer[index_k];

          pwb->gamma_integrand[thread][index_k] = tr * pk_one_third;
          pwb->delta_integrand[thread][index_k] = tr * pk_two_thirds;

        }
     
      } // end of equilateral and orthogonal model





      // =========================================================
      // =         Convolve with the Bessel function             =
      // =========================================================
      
      for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
        double r = pwb->r[index_r];
      
        if (pbi->bispectra_verbose > 3)
          printf("       \\ r=%g, index_r=%d\n", r, index_r);
        
        
        // *** All models ***
        if ((pbi->has_local_model == _TRUE_) ||
            (pbi->has_equilateral_model == _TRUE_) ||
            (pbi->has_orthogonal_model == _TRUE_)) {
                      
          /* alpha integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->q,
                        pbi->delta_k,
                        ptr->q_size,
                        pwb->alpha_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->alpha[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);

          /* beta integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->q,
                        pbi->delta_k,
                        ptr->q_size,
                        pwb->beta_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->beta[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);
                
        } // end of local model
        
        // /* Some debug */
        // if (index_r == 70)
        //   fprintf (stderr, "%10d %17.7g %17.7g\n",
        //   pbi->l[index_l], pbi->alpha[index_bf][index_l][index_r], pbi->beta[index_bf][index_l][index_r]);

        
        // *** Equilateral and orthogonal models ***
        if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {

          /* gamma integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->q,
                        pbi->delta_k,
                        ptr->q_size,
                        pwb->gamma_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->gamma[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);

          /* delta integral */
          class_call_parallel (bessel_convolution (
                        ppr,
                        pbs,
                        ptr->q,
                        pbi->delta_k,
                        ptr->q_size,
                        pwb->delta_integrand[thread],
                        NULL,
                        index_l,
                        r,
                        &(pbi->delta[index_bf][index_l][index_r]),
                        pbi->error_message
                        ),
            pbi->error_message,
            pbi->error_message);
      
        } // end of equilateral and orthogonal models
        
        #pragma omp flush(abort)
      
      } // end of for(index_r)
  
    } // end of for(index_l)
    
  } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */

    
  /* Output the filter functions */
  // int index_l = 60;
  // 
  // if ((pbi->has_bispectra_t == _TRUE_) && (index_bf == pbi->index_bf_t) ) {
  // 
  //   /* Output beta = b_l */ 
  //   fprintf (stderr, "# Printing BETA filter function for l=%d\n", pbi->l[index_l]);
  //   fprintf (stderr, "%17s %17s\n", "r", "beta(l,r)");
  //   for (int index_r = 0; index_r < pwb->r_size; ++index_r)
  //     fprintf (stderr, "%17.7g %17.7g\n", pwb->r[index_r], pbi->beta[index_bf][index_l][index_r]);
  // 
  //   /* Output alpha = b_nl */  
  //   fprintf (stderr, "# Printing ALPHA filter function for l=%d\n", pbi->l[index_l]);
  //   fprintf (stderr, "%17s %17s\n", "r", "alpha(l,r)");
  //   for (int index_r = 0; index_r < pwb->r_size; ++index_r)
  //     fprintf (stderr, "%17.7g %17.7g\n", pwb->r[index_r], pbi->alpha[index_bf][index_l][index_r]);
  //   }
  // }

  return _SUCCESS_;
  
}









int bispectra_separable_integrate_over_r (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int i, // index for the first field (T,B,E)
    int j, // index for the second field (T,B,E)
    int k, // index for the third field (T,B,E)
    struct bispectra_workspace_separable * pwb
    )

{

  if (pbi->bispectra_verbose > 2)
    printf("     * computing the r-integral for the bispectrum %s_%s%s%s\n",
    pbi->bt_labels[index_bt], pbi->bf_labels[i], pbi->bf_labels[j], pbi->bf_labels[k]);

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  /* Shortcuts */
  double *** alpha = pbi->alpha;
  double *** beta = pbi->beta;
  double *** gamma = pbi->gamma;
  double *** delta = pbi->delta;
  
  // =======================================================================================
  // =                              Compute bispectrum r-integral                          =
  // =======================================================================================

  /* We now proceed to the final step of the bispectrum computation with the smooth_separable technique. This time
  we shall loop over all the l1-l2-l3 configurations, compute the simple integral over r, and store the bispectrum
  in pbi->bispectra.

    When computing the primordial bispectra, we shall assume a primordial f_NL of unity for phi. So far, we
  derived transfer functions for R. The curvature perturbation R and the potential during matter domination
  phi are related by R = -5/3 phi. The bispectrum integral has two power spectra (equivalent to 4 R's)
  and 3 transfer functions.  When converting to phi, each R brings a (-5/3) while each transfer function
  brings a (-3/5). We are then left with an overall factor of (-5/3), which means that a bispectrum
  with fnl_phi=1 is equivalent to a bispectrum with fnl_R=-3/5.
  
     We write down the r-integral using the same notation as in Yadav & Wandelt 2010 (eq. 32,
  35 and 38 of http://arxiv.org/abs/1006.0275v3), so that the i,j,k indices there correspond 
  to those in this function (T, E, ...). See also eq. 12 of Yadav & Wandelt 2007. */

  double fnl_R = -3/5.;

  /* We parallelize the outer loop over 'l1'. */
  #pragma omp parallel                                     \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)       \
    private (thread)
  {
  
    /* Cycle variables, declared inside parallel loop */ 
    int index_l1, index_l2, index_l3, index_r;
    
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      if (pbi->bispectra_verbose > 2)
        printf("      \\ computing r-integral for l1=%d, index_l1=%d\n", pbi->l[index_l1], index_l1);
    
      for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (pbi->l[index_l2] < pbi->l[index_l1]/2)
          continue;
  
        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

          /* Index of the current (l1,l2,l3) configuration */
          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  
          double integral = 0;
  
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
            /* Value of the considered r */
            double r = pwb->r[index_r];
            double integrand;
    
            // ------------------------------------------------------------------------
            // -                               Local model                            -
            // ------------------------------------------------------------------------

            if ((pbi->has_local_model == _TRUE_) && (index_bt == pbi->index_bt_local)) {
  
              /* The primordial bispectrum for the local model has only two extra contributions due to the symmetrization
              P(k1)*P(k2) + P(k1)*P(k3) + P(k2)*P(k3). It is easy to see that each contribution has the same value, but
              we sum them nonetheless to be consistent. */
              integrand = 2 * fnl_R * r*r * (
                  alpha[i][index_l1][index_r] * beta[j][index_l2][index_r]  * beta[k][index_l3][index_r]
                + beta[i][index_l1][index_r]  * beta[j][index_l2][index_r]  * alpha[k][index_l3][index_r]
                + beta[i][index_l1][index_r]  * alpha[j][index_l2][index_r] * beta[k][index_l3][index_r]
                );
            
            } // end of local model

            // -------------------------------------------------------------------------
            // -                           Equilateral model                           -
            // -------------------------------------------------------------------------
            
            else if ((pbi->has_equilateral_model == _TRUE_) && (index_bt == pbi->index_bt_equilateral)) {
  
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - alpha[i][index_l1][index_r] * beta[j][index_l2][index_r]  * beta[k][index_l3][index_r]
                - beta[i][index_l1][index_r]  * beta[j][index_l2][index_r]  * alpha[k][index_l3][index_r]
                - beta[i][index_l1][index_r]  * alpha[j][index_l2][index_r] * beta[k][index_l3][index_r]
                    
                /* Symmetrical part */
                - 2 * delta[i][index_l1][index_r] * delta[j][index_l2][index_r] * delta[k][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + beta[i][index_l1][index_r]  * gamma[j][index_l2][index_r] * delta[k][index_l3][index_r]
                + beta[i][index_l1][index_r]  * delta[j][index_l2][index_r] * gamma[k][index_l3][index_r]
                + gamma[i][index_l1][index_r] * beta[j][index_l2][index_r]  * delta[k][index_l3][index_r]
                + gamma[i][index_l1][index_r] * delta[j][index_l2][index_r] * beta[k][index_l3][index_r]
                + delta[i][index_l1][index_r] * beta[j][index_l2][index_r]  * gamma[k][index_l3][index_r]
                + delta[i][index_l1][index_r] * gamma[j][index_l2][index_r] * beta[k][index_l3][index_r]
              
              );

            } // end of equilateral model

            // -------------------------------------------------------------------------------
            // -                              Orthogonal model                               -
            // -------------------------------------------------------------------------------
            
            else if ((pbi->has_orthogonal_model == _TRUE_) && (index_bt == pbi->index_bt_orthogonal)) {
  
              /* We take the formula from Senatore et al. 2010, also shown in Komatsu et al. 2011 (WMAP7 paper, eq. 64) */
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - 3 * alpha[i][index_l1][index_r] * beta[j][index_l2][index_r]  * beta[k][index_l3][index_r]
                - 3 * beta[i][index_l1][index_r]  * beta[j][index_l2][index_r]  * alpha[k][index_l3][index_r]
                - 3 * beta[i][index_l1][index_r]  * alpha[j][index_l2][index_r] * beta[k][index_l3][index_r]
                    
                /* Symmetrical part. We found what we think is a typo in eq. 38 of http://arxiv.org/abs/1006.0275v3,
                where the coefficient is -2/3*3 = -2 instead of -8. We think -8 is the correct coefficient, as it can
                be verified from eq. 3.2 of Senatore, Smith & Zaldarriaga 2010.  */
                - 8 * delta[i][index_l1][index_r] * delta[j][index_l2][index_r] * delta[k][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + 3 * beta[i][index_l1][index_r]  * gamma[j][index_l2][index_r] * delta[k][index_l3][index_r]
                + 3 * beta[i][index_l1][index_r]  * delta[j][index_l2][index_r] * gamma[k][index_l3][index_r]
                + 3 * gamma[i][index_l1][index_r] * beta[j][index_l2][index_r]  * delta[k][index_l3][index_r]
                + 3 * gamma[i][index_l1][index_r] * delta[j][index_l2][index_r] * beta[k][index_l3][index_r]             
                + 3 * delta[i][index_l1][index_r] * beta[j][index_l2][index_r]  * gamma[k][index_l3][index_r]
                + 3 * delta[i][index_l1][index_r] * gamma[j][index_l2][index_r] * beta[k][index_l3][index_r]             
              
              );

            } // end of orthogonal model

            
            /* Increment the estimate of the integral */
            integral += integrand * pwb->delta_r[index_r];

            /* Print integrand as a function of r */
            // if ((index_l1 == pbi->l_size-1) && (index_l2 == pbi->l_size-1) && (index_l3 == pbi->l_size-1))
            // if ((index_l1 == 0) && (index_l2 == 0) && (index_l3 == 0))
            // if ((index_l1 == 0) && (index_l2 == 0) && (index_l3 == 0))
            //   fprintf (stderr, "%15.7g %15.7g\n", r, integrand);

  
          } // end of for(index_r)
  

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3), including the factor 1/2
          from trapezoidal rule */
          pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] = 0.5 * integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. In KSW2005, this factor
          was split between the alpha and beta integrals, but from the numerical point of view it is preferable
          to include it at the end of the computation (as in eq. 17 of Fergusson & Shellard 2007). */
          pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] *= pow(2./_PI_,3);
  
          /* Some debug */
          // if ((index_l1==index_l2) && (index_l2==index_l3))
          //   fprintf (stderr, "%10d %17.7g\n", pbi->l[index_l1], pbi->bispectra[index_bt][index_l1][index_l2][index_l3-index_l3_min]);
    
          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_for_bispectra++;

        } // end of for(index_l3)
      } // end of for(index_l2)
      
      #pragma omp flush(abort)
  
    } // end of for(index_l1)
  
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  return _SUCCESS_; 
  
}
  


/**
 * Compute the angular bispectrum for the separable primordial bispectra.
 *
*/
int bispectra_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{


  // ===================================================================================
  // =                                Initialise workspace                             =
  // ===================================================================================

  /* Allocate arrays inside the integration workspace */
  class_call (bispectra_separable_workspace_init(
                ppr,
                pba,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);


  if (pbi->bispectra_verbose > 1)
    printf (" -> computing separable bispectra; r sampled %d times in [%g,%g]\n", pwb->r_size, pwb->r_min, pwb->r_max);
    




  // ===================================================================================
  // =                             Compute filter functions                            =
  // ===================================================================================

  /* Compute the filter functions: alpha(l,r), beta(l,r), gamma(l,r), delta(l,r) for each field 
  (T,E,B). The filter functions are convolution of P(k)^A*T(k)^B with a Bessel function. */

  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

    if (pbi->bispectra_verbose > 1)
        printf("     * computing %s filter functions for the separable bispectra ...\n", pbi->bf_labels[index_bf]);

    class_call (bispectra_separable_filter_functions(
                  ppr,
                  pba,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  index_bf,
                  pwb),
      pbi->error_message,
      pbi->error_message);

  }


  // ==========================================================================================
  // =                               Perform final integration                                =
  // ==========================================================================================

  /* Compute the bispectra by integrating the filter functions in r. This is where the primordial shape
  functions are used. */
  
  /* Loop on the type of bispectrum (local, equilateral, orthogonal...) */
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {
  
    if (pbi->bispectrum_type[index_bt] != separable_bispectrum)
      continue;

    /* Loop on the considered probe (TTT, TTE, TEE, EEE...) */
    for (int index_bf_1 = 0; index_bf_1 < pbi->bf_size; ++index_bf_1) {
    for (int index_bf_2 = 0; index_bf_2 < pbi->bf_size; ++index_bf_2) {
    for (int index_bf_3 = 0; index_bf_3 < pbi->bf_size; ++index_bf_3) {
  
      class_call (bispectra_separable_integrate_over_r(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    index_bt,
                    index_bf_1,
                    index_bf_2,
                    index_bf_3,
                    pwb),
        pbi->error_message,
        pbi->error_message);
  
    }}} // end of for(ijk)
    
  } // end of for(index_bf)
  
  
  return _SUCCESS_; 
  
}






int bispectra_separable_workspace_free (
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{

  free (pwb->r);
  free (pwb->delta_r);
 
  /* Arrays specific to the primordial models */

  int thread = 0;
  
  #pragma omp parallel private (thread)
  {
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
      free (pwb->alpha_integrand[thread]);
      free (pwb->beta_integrand[thread]);
    }

    if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
      free (pwb->gamma_integrand[thread]);
      free (pwb->delta_integrand[thread]);
    }
  }

  if ((pbi->has_local_model == _TRUE_) || (pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
    free (pwb->alpha_integrand);
    free (pwb->beta_integrand);
  }

  if ((pbi->has_equilateral_model == _TRUE_) || (pbi->has_orthogonal_model == _TRUE_)) {
    free (pwb->gamma_integrand);
    free (pwb->delta_integrand);
  }
 
  return _SUCCESS_; 
  
}






/**
 * Compute the angular bispectrum for the analytical primordial bispectra.
 *
 */
int bispectra_analytical_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
    struct bispectra * pbi
    )
{  

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
    
    if (pbi->bispectrum_type[index_bt] != analytical_bispectrum)
      continue;

      for (int i = 0; i < pbi->bf_size; ++i) {
      for (int j = 0; j < pbi->bf_size; ++j) {
      for (int k = 0; k < pbi->bf_size; ++k) {

      if (pbi->bispectra_verbose > 1)
        printf(" -> computing the bispectrum (%s_%s)\n",
          pbi->bt_labels[index_bt], pbi->bfff_labels[i][j][k]);

      /* We parallelize the outer loop over 'l1'. */
      #pragma omp parallel                                     \
        shared (ppt,pbs,ptr,ppm,pbi,abort)       \
        private (thread)
      {
  
        /* Cycle variables, declared inside parallel loop */ 
        int index_l1, index_l2, index_l3, index_r;
    
        #pragma omp for schedule (dynamic)
        for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
    
          int l1 = pbi->l[index_l1];
          double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    
          for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
            /* Skip those configurations that are forbidden by the triangular condition (optimization) */
            if (pbi->l[index_l2] < pbi->l[index_l1]/2)
              continue;

            int l2 = pbi->l[index_l2];
            double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
  
            /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
            int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
            int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
            for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

              int l3 = pbi->l[index_l3];
              double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];

              /* Index of the current (l1,l2,l3) configuration */
              long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

              /* Initialize bispectrum */
              pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] = 0;
          
              // ---------------------------------------------------------------------------------------------
              // -                          Squeezed limit of the intrinsic bispectrum                       -
              // ---------------------------------------------------------------------------------------------
          
              /* We use the Eq. 1 from Huang & Vernizzi 2013; in the language of Lewis et al. 2011,
              it consists of a redshift modulation term plus Ricci focussing. */
              // if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed)) {
              // 
              //   double C_l1_tz = pbi->cls[psp->index_ct_tz][l1-2];
              //   double dlog_lsq_cl_l1 = pbi->dlog_cls[psp->index_ct_tt][l1-2];
              // 
              //   double C_l2_tz = pbi->cls[psp->index_ct_tz][l2-2];
              //   double dlog_lsq_cl_l2 = pbi->dlog_cls[psp->index_ct_tt][l2-2];
              // 
              //   double C_l3_tz = pbi->cls[psp->index_ct_tz][l3-2];
              //   double dlog_lsq_cl_l3 = pbi->dlog_cls[psp->index_ct_tt][l3-2];
              // 
              //   /* Redshift (local modulation) */
              //   pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] += (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
              //   
              //   /* Ricci focussing (scale redefinition) */
              //   pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] += _ONE_THIRD_ *
              //     (- 0.5 * C_l1_tz * (C_l2*dlog_lsq_cl_l2) + C_l2_tz * (C_l1*dlog_lsq_cl_l1) +
              //      - 0.5 * C_l1_tz * (C_l3*dlog_lsq_cl_l3) + C_l3_tz * (C_l1*dlog_lsq_cl_l1) +
              //      - 0.5 * C_l2_tz * (C_l3*dlog_lsq_cl_l3) + C_l3_tz * (C_l2*dlog_lsq_cl_l2));
              // 
              //   /* Uncomment to use the unsymmetrised version */
              //   // pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] = 0;
              //   // pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] += C_l3 * (C_l1 + C_l2);
              //   // pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] += - 0.5 * C_l3_tz * (C_l1 * dlog_lsq_cl_l1 + C_l2 * dlog_lsq_cl_l2);
              // 
              // }

              // ---------------------------------------------------------------------------------------------
              // -                          Squeezed limit of the local bispectrum                           -
              // ---------------------------------------------------------------------------------------------

              if ((pbi->has_local_squeezed == _TRUE_) && (index_bt == pbi->index_bt_local_squeezed)) {
              
                pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] = 6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);

              }

              // ---------------------------------------------------------------------------------------------
              // -                          Squeezed limit of the local bispectrum                           -
              // ---------------------------------------------------------------------------------------------

              if ((pbi->has_cosine_shape == _TRUE_) && (index_bt == pbi->index_bt_cosine)) {
              
                double cosine = 0;
                pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] = 6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) * cosine;

              }

              // ----------------------------------------------------------------------------------------------
              // -                                   Lensing-ISW bispectrum                                   -
              // ----------------------------------------------------------------------------------------------
            
              /* TODO: This is not the correct formula, it is just a test. For the correct formula, refer
              to eq. 4.5 of Lewis, Challinor & Hanson 2011. */
              if ((pbi->has_isw_lensing == _TRUE_) && (index_bt == pbi->index_bt_isw_lensing)) {
              
                pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3] +=
                  + (l1*(l1+1) + l2*(l2+1) - l3*(l3+1))
                  + l1*(l1+1) - l2*(l2+1) + l3*(l3+1)
                  - l1*(l1+1) + l2*(l2+1) + l3*(l3+1);
                
              }

          
              /* Update the counter */
              #pragma omp atomic
              pbi->count_memorised_for_bispectra++;

            } // end of for(index_l3)
          } // end of for(index_l2)
      
          #pragma omp flush(abort)
  
        } // end of for(index_l1)
  
      } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
    
    }}} // end of for(ijk)

  } // end of for(index_bt)

  
  return _SUCCESS_;

}







int bispectra_non_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Allocate arrays inside the integration workspace */
  class_call (bispectra_non_separable_workspace_init(
                ppr,
                pba,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);
  

  if (pbi->bispectra_verbose > 1)
    printf (" -> computing non-separable bispectra; r sampled %d times in [%g,%g]\n", pwb->r_size, pwb->r_min, pwb->r_max);
  
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* Skip the bispectrum if it not of the non-separable type */
    if (pbi->bispectrum_type[index_bt] != non_separable_bispectrum)
      continue;

    for (int i = 0; i < pbi->bf_size; ++i) {
    for (int j = 0; j < pbi->bf_size; ++j) {
    for (int k = 0; k < pbi->bf_size; ++k) {

      if (pbi->bispectra_verbose > 1)
        printf("     * computing bispectrum (%s_%s)\n",
        pbi->bt_labels[index_bt], pbi->bfff_labels[i][j][k]);
    
      // ==================================================================================================
      // =                              Determine the bispectrum to compute                               =
      // ==================================================================================================
        
      /* Which primordial shape is needed? */
      if ((pbi->has_galileon_model==_TRUE_) && (index_bt==pbi->index_bt_galileon_gradient))
        pwb->shape_function = bispectra_galileon_gradient;

      else if ((pbi->has_galileon_model==_TRUE_) && (index_bt==pbi->index_bt_galileon_time))
        pwb->shape_function = bispectra_galileon_time;

      /* Debug - Uncomment to compute a custom bispectrum; useful to test against the separable shapes */
      // pwb->shape_function = bispectra_orthogonal_model;

      /* Which probe is needed? We define three indices that point to the transfer functions 
      computed in the transfer module. For example, for a TTE bispectrum we could have
      index_tt_k1=T, index_tt_k2=T, index_tt_k3=E. */
      int index_tt_k1 = pbi->index_tt_of_bf[i];
      int index_tt_k2 = pbi->index_tt_of_bf[j];
      int index_tt_k3 = pbi->index_tt_of_bf[k];

      /* Compute fist integral over k3 */
      class_call (bispectra_non_separable_integrate_over_k3(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    index_bt,
                    index_tt_k3,
                    pwb),
        pbi->error_message,
        pbi->error_message);

      /* Compute second integral over k2 */
      class_call (bispectra_non_separable_integrate_over_k2(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    index_bt,
                    index_tt_k2,
                    pwb),
        pbi->error_message,
        pbi->error_message);
      
      /* Compute the third integral over k1 */
      class_call (bispectra_non_separable_integrate_over_k1(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    index_bt,
                    index_tt_k1,
                    pwb),
        pbi->error_message,
        pbi->error_message);
      
      /* Compute the fourth and last integral over r */
      class_call (bispectra_non_separable_integrate_over_r(
                    ppr,
                    pba,
                    ppt,
                    pbs,
                    ptr,
                    ppm,
                    pbi,
                    pbi->bispectra[index_bt][i][j][k],
                    pwb),
        pbi->error_message,
        pbi->error_message);

    }}} // end of for(ijk)

  } // end of for(index_bt)
  
  return _SUCCESS_;
  
}






int bispectra_non_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{

  // ===================================================================================================
  // =                                    Prepare integration grid                                     =
  // ===================================================================================================


  // ------------------------------------------------------------
  // -                         Grid in r                        -
  // ------------------------------------------------------------
  
  /* We set the r-sampling as if it were a time sampling.  We do so because 'r' has the right dimensions, and
  it always appears in the argument of a Bessel function multiplying a wavemode, just as it was for conformal
  time in the line-of-sight integral.  This is the only place in the module where the background structure
  is accessed. */
  pwb->r_min = pbi->r_min;
  pwb->r_max = pbi->r_max;
  pwb->r_size = pbi->r_size;
    

  /* We decide to sample r linearly */
  class_alloc (pwb->r, pwb->r_size*sizeof(double), pbi->error_message);
  lin_space (pwb->r, pwb->r_min, pwb->r_max, pwb->r_size);
    
  /* Allocate & fill delta_r, the measure for the trapezoidal integration over r */
  class_alloc (pwb->delta_r, pwb->r_size * sizeof(double), pbi->error_message);

  /* Fill pwb->delta_r */
  pwb->delta_r[0] = pwb->r[1] - pwb->r[0];
      
  for (int index_r=1; index_r < pwb->r_size-1; ++index_r)
    pwb->delta_r[index_r] = pwb->r[index_r+1] - pwb->r[index_r-1];
      
  pwb->delta_r[pwb->r_size-1] = pwb->r[pwb->r_size-1] - pwb->r[pwb->r_size-2];



  // -----------------------------------------------------------------------
  // -                       Grid for the shape function                   -
  // -----------------------------------------------------------------------
  
  
  /* We shall sample the primordial shape function only in the k-points determined in the perturbation
  module, which are much sparsely sampled than those of the transfer functions. Hence, we are assuming
  that the shape function is a smooth function of (k1,k2,k3) */
    
  pwb->k_smooth_grid = ppt->k;
  pwb->k_smooth_size = ppt->k_size;
  
    
  
  
  // -----------------------------------------------------------------------
  // -                       Grid for the shape function                   -
  // -----------------------------------------------------------------------

  /* Here we set the integration limits on k3. These can be made equal to those of k1 and k2 (that is equal to
  the range where we computed the transfer functions) but we can do better than that. In fact, 
  the r-integral enforces the triangular condition |k1-k2| <= k3 <= k1+k2 which means that we can restrict
  our range to those configurations. However, we should not get too close to the triangular limits otherwise
  the integral becomes numerically instable. */

  /* Sampling of the first-order transfer functions.  */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  pwb->k3_size_max = k_tr_size-1;
  
  class_alloc (pwb->index_k3_lower, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  class_alloc (pwb->index_k3_upper, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  class_alloc (pwb->k3_grid_size, pwb->k_smooth_size*sizeof(int *), pbi->error_message);
  
  for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
  
    double k1 = pwb->k_smooth_grid[index_k1];
    class_alloc (pwb->index_k3_lower[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
    class_alloc (pwb->index_k3_upper[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
    class_alloc (pwb->k3_grid_size[index_k1], (index_k1+1)*sizeof(int), pbi->error_message);
        
    for (int index_k2=0; index_k2 <= index_k1; ++index_k2) {

      double k2 = pwb->k_smooth_grid[index_k2];
      
      /* Uncomment to use the triangular condition. This will give an imprecise result, one
      needs to extend the range a bit, same as we do for the second-order bispectrum. */
      // double k_tr_min = fabs(k1-k2);
      // double k_tr_max = k1+k2;

      /* Uncomment to take the whole k3 range for the integration (safer) */
      double k_tr_min = k_tr[0];
      double k_tr_max = k_tr[k_tr_size-1];
      
      /* Find the index corresponding to k3_lower inside k_tr */
      int index_k3_lower = 0;
      while (k_tr[index_k3_lower] < k_tr_min) ++index_k3_lower;
      pwb->index_k3_lower[index_k1][index_k2] = index_k3_lower;

      /* Find the index corresponding to k3_upper inside ptr->k */
      int index_k3_upper = k_tr_size - 1;
      while (k_tr[index_k3_upper] > k_tr_max) --index_k3_upper;
      pwb->index_k3_upper[index_k1][index_k2] = index_k3_upper;
      
      /* Number of points in the k_tr grid between 'k3_lower' and 'k3_upper' */
      pwb->k3_grid_size[index_k1][index_k2] = index_k3_upper - index_k3_lower + 1;

      /* Some debug - print out the k3_grid list for a special configuration */      
      // if ((index_k1==100) && (index_k2>-1)) {
      //   fprintf (stderr, "k1[%d]=%.5e, k2[%d]=%.5e, k3_grid_size=%d, k3_min=%.5e, k3_max=%.5e\n",
      //     index_k1, k_tr[index_k1], index_k2, k_tr[index_k2], pwb->k3_grid_size[index_k1][index_k2], k_tr_min, k_tr_max);
      //   for (int index_k3=0; index_k3 < pwb->k3_grid_size[index_k1][index_k2]; ++index_k3)
      //     fprintf(stderr, "%d %.17f /\\ ", index_k3, k_tr[index_k3]);
      // 
      //   fprintf (stderr, "\n\n");
      // }

    
    } // end of for(k2)
  } // end of for(k1)



  // ========================================================================================
  // =                               Determine window function                              =
  // ========================================================================================
  
  
  /* Determine the window function for the interpolation of the k-space bispectrum in k1 and k2. We need
  a window function because the primordial bispectrum will usually contain products of power spectra 
  that diverge as k^-3 for k -> 0. */
    
  /* Window function will have pbi->k_smooth_size elements */
  class_alloc (pwb->k_window, pwb->k_smooth_size*sizeof(double), pbi->error_message);

  for (int index_k=0; index_k < pwb->k_smooth_size; ++index_k) {
    
    double k = pwb->k_smooth_grid[index_k];
    double pk = pbi->pk_pt[index_k];
    
    pwb->k_window[index_k] = pow(k,2);

  }

  /* Inverse window function will have ptr->q_size elements */
  class_alloc (pwb->k_window_inverse, ptr->q_size*sizeof(double), pbi->error_message);
  for (int index_k=0; index_k < ptr->q_size; ++index_k) {

    double k = ptr->q[index_k];
    double pk = pbi->pk[index_k];

    pwb->k_window_inverse[index_k] = 1./pow(k,2);

  }



  // =========================================================================================
  // =                                    Allocate arrays                                    =
  // =========================================================================================
  
  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel private (thread)
  {
    #ifdef _OPENMP
    number_of_threads = omp_get_num_threads();
    #endif
  }
  
  /* We need a k3_grid per thread because it varies with k1 and k2 due to the triangular condition */
  class_alloc (pwb->k3_grid, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->delta_k3, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->integral_splines, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->interpolated_integral, number_of_threads*sizeof(double*), pbi->error_message);
  class_alloc (pwb->f, number_of_threads*sizeof(double*), pbi->error_message);
    
  /* Beginning of parallel region */
  abort = _FALSE_;
  #pragma omp parallel shared(pwb,pbi) private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    /* Allocate the integration grid in k3 with the maximum possible number of k3-values. This is given by the
    number of k-values in ptr->k */
    class_alloc_parallel(pwb->k3_grid[thread], pwb->k3_size_max*sizeof(double), pbi->error_message);
    class_alloc_parallel(pwb->delta_k3[thread], pwb->k3_size_max*sizeof(double), pbi->error_message);
  
    /* Allocate memory for the interpolation arrays (used only for the k2 and k3 integrations) */
    class_alloc_parallel (pwb->integral_splines[thread], ptr->q_size*sizeof(double), pbi->error_message);
    class_alloc_parallel (pwb->interpolated_integral[thread], ptr->q_size*sizeof(double), pbi->error_message);
    class_alloc_parallel (pwb->f[thread], pwb->k_smooth_size*sizeof(double), pbi->error_message);
  
  } // end of parallel region
  
  if (abort == _TRUE_) return _FAILURE_;

  return _SUCCESS_;
  
} // end of bispectra_non_separable_workspace_init







int bispectra_non_separable_integrate_over_k3 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k3,
    struct bispectra_workspace_non_separable * pwb
    )
{
  

  /* Indices */
  int index_r;
  int index_k1, index_k2, index_k3;
  int index_l1, index_l2, index_l3;
  
  /* Integration grid */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;


  // =========================================================================================
  // =                             Allocate memory for I_l3(r,k1,k2)                         =
  // =========================================================================================
  
    
  /* Initialize counter */
  pwb->count_allocated_for_integral_over_k3 = 0;
  
  
  /* Allocate l3-level.  Note that, even if l3 must satisfy the triangular
  inequality, we allocate this level for all the possible l3 values.  We do so because this
  array is going to be used by all (l1,l2) computations that follow, which means that l3 will
  eventually cover all the allowed range. */
  class_alloc (pwb->integral_over_k3, pbi->l_size*sizeof(double ***), pbi->error_message);
  
  
  for (index_l3=0; index_l3<pbi->l_size; ++index_l3) {
    
    /* Allocate r-level */
    class_alloc (pwb->integral_over_k3[index_l3], pwb->r_size*sizeof(double **), pbi->error_message);
  
    /* Allocate 'k1' level */
    for (index_r=0; index_r < pwb->r_size; ++index_r) {
  
      int k1_size = pwb->k_smooth_size;
      class_alloc (pwb->integral_over_k3[index_l3][index_r], k1_size*sizeof(double *), pbi->error_message);

      /* Allocate 'k2' level */
      for (index_k1=0; index_k1<k1_size; ++index_k1) {
  
        int k2_size = index_k1 + 1;
        class_calloc (pwb->integral_over_k3[index_l3][index_r][index_k1], k2_size, sizeof(double), pbi->error_message);
        
        /* Increase memory counter */
        pwb->count_allocated_for_integral_over_k3 += k2_size;
  
      } // end of for(index_k1)
  
    } // end of for(index_r)
    
  } // end of for(index_l3)
    
  if (pbi->bispectra_verbose > 1)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the k3-integral array (k_size=%d)\n",
      pwb->count_allocated_for_integral_over_k3*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k3, pwb->k_smooth_size);
  



  // ===================================================================================================
  // =                               Compute the INT_l3(r, k1, k2)  integral                           =
  // ===================================================================================================
  
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k3 = 0;

  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)       \
    private (index_l3,index_k1,index_k2,index_k3,index_r,thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
  
    #pragma omp for schedule (dynamic)
    for (index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {

      double k1 = pwb->k_smooth_grid[index_k1];
      double pk_1 = pbi->pk_pt[index_k1];
      
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k3 integral for k1=%g, index_k1=%d\n", pwb->k_smooth_grid[index_k1], index_k1);
  
      /* We only need to consider those k2's that are equal to or larger than k1,
        as the shape function is assumed to be symmetric woth respect to k1<->k2<->k3 */      
      for (index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
    
        double k2 = pwb->k_smooth_grid[index_k2];
        double pk_2 = pbi->pk_pt[index_k2];
    
        /* Get the size of the integration grid. Note that when extrapolation is turned on, the k3-grid will also
          include values that do not satisfty the triangular condition k1 + k2 = k3. */
        int k3_size = pwb->k3_grid_size[index_k1][index_k2];
        int index_k3_lower = pwb->index_k3_lower[index_k1][index_k2];
        int index_k3_upper = pwb->index_k3_upper[index_k1][index_k2];
  
        /* Pointer to the lower limit of integration. */
        pwb->k3_grid[thread] = k_tr + index_k3_lower;
                    
        /* If there are no points in k_tr that satisfy the triangular condition for the current (k1,k2), then
          there is no contribution to the integral */
        class_test_parallel (k3_size<=0, pbi->error_message, "include when triangular condition does not fit with ptr->k");
          
        /* Determine the measure for the trapezoidal rule for k3 */  
        pwb->delta_k3[thread][0] = pwb->k3_grid[thread][1] - pwb->k3_grid[thread][0];
  
        for (index_k3=1; index_k3<(k3_size-1); ++index_k3)
          pwb->delta_k3[thread][index_k3] = pwb->k3_grid[thread][index_k3 + 1] - pwb->k3_grid[thread][index_k3 - 1];
  
        pwb->delta_k3[thread][k3_size-1] = pwb->k3_grid[thread][k3_size - 1] - pwb->k3_grid[thread][k3_size - 2];
        
        /* Shape function for this (k1,k2) slice */
        for (index_k3=index_k3_lower; index_k3 <= index_k3_upper; ++index_k3) {
          
          double k3 = k_tr[index_k3];
          double pk_3 = pbi->pk[index_k3];

          class_call_parallel (pwb->shape_function (
                        ppm,
                        pbi,
                        k1, k2, k3,
                        pk_1, pk_2, pk_3,
                        &pwb->interpolated_integral[thread][index_k3]),
            pbi->error_message,
            pbi->error_message);
          
        }
        
  
        /* We compute the integral over k3 for all possible l-values */
        for (index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

          /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
          int tt_size = ptr->tt_size[ppt->index_md_scalars];
          int l_size = ptr->l_size[ppt->index_md_scalars];

          double * transfer = &(ptr->transfer
            [ppt->index_md_scalars]
            [((ppt->index_ic_ad * tt_size + index_tt_k3) * l_size + index_l3) * k_tr_size]);

          /* Some debug - print transfer function */
          // for (index_k3=0; index_k3 < k3_size; ++index_k3)
          //   if ((index_k1==3) && (index_k2==2) && (index_l3==0))
          //     fprintf(stderr, "%g %g\n", pwb->k3_grid[thread][index_k3], transfer[index_k3]);
    
          for (index_r = 0; index_r < pwb->r_size; ++index_r) {
  
            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr + index_k3_lower,
                          pwb->delta_k3[thread],
                          k3_size,
                          transfer + index_k3_lower,
                          pwb->interpolated_integral[thread] + index_k3_lower,
                          index_l3,
                          pwb->r[index_r],
                          &(pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
                

            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k3;
  
            #pragma omp flush(abort)
  
          } // end of for(index_r)          
        } // end of for(index_l3)
      } // end of for(index_k2)
    } // end of for(index_k1)
  } if (abort == _TRUE_) return _FAILURE_; /* end of parallel region */  
  
  
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k3-integral array\n",
      pwb->count_memorised_for_integral_over_k3*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k3);
  
  /* Check that we correctly filled the array */
  class_test_permissive (pwb->count_memorised_for_integral_over_k3 != pwb->count_allocated_for_integral_over_k3,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k3, pwb->count_memorised_for_integral_over_k3);

  return _SUCCESS_;
  
}








int bispectra_non_separable_integrate_over_k2 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k2,
    struct bispectra_workspace_non_separable * pwb
    )
{
  

  /* Indices */
  int index_r;
  int index_k1, index_k2, index_k3;
  int index_l1, index_l2, index_l3;
  
  /* Integration grid */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  
  
  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  

  // ======================================================================================================
  // =                                   Allocate memory for INT_l2_l3(r,k1)                              =
  // ======================================================================================================
  
  /* Initialize counter */
  pwb->count_allocated_for_integral_over_k2 = 0;
  
  /* Allocate l2-level */
  class_alloc (pwb->integral_over_k2, pbi->l_size*sizeof(double ***), pbi->error_message);
  
  for (index_l2=0; index_l2<pbi->l_size; ++index_l2) {
  
    /* Allocate l3-level. We only need l3<=l2 because of the k2<->k3 symmetry of the shape function */
    class_alloc (pwb->integral_over_k2[index_l2], (index_l2+1)*sizeof(double **), pbi->error_message);
  
    for (index_l3=0; index_l3<=index_l2; ++index_l3) {
    
      /* Allocate r-level */
      class_alloc (pwb->integral_over_k2[index_l2][index_l3], pwb->r_size*sizeof(double *), pbi->error_message);
  
      /* Allocate 'k1' level */
      for (index_r=0; index_r < pwb->r_size; ++index_r) {
  
        int k1_size = pwb->k_smooth_size;
        class_alloc (pwb->integral_over_k2[index_l2][index_l3][index_r], k1_size*sizeof(double), pbi->error_message);
  
        /* Increase memory counter */
        pwb->count_allocated_for_integral_over_k2 += k1_size;
  
      } // end of for(index_k1)
  
    } // end of for(index_r)
    
  } // end of for(index_l3)
    
  if (pbi->bispectra_verbose > 1)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the k2-integral array\n",
      pwb->count_allocated_for_integral_over_k2*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k2);
  



  // ==============================================================================================================
  // =                                   Compute the INT_l2_l3(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k2 = 0;
    
  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)       \
    private (index_r,index_l2,index_l3,index_k1,thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k2 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);
  
      for (index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {
          
        for (index_l2 = 0; index_l2 < pbi->l_size; ++index_l2) {
      
          /* Interpolate the integral I_l3(k1,k2,r) that we computed above in the integration grid of k2.
          Note that we pass the integral_splines, interpolated_integral and pwb->f arrays separately rather
          than accessing them from pwb, because they are thread dependent (while pwb isn't). */
          class_call_parallel (bispectra_non_separable_interpolate_over_k2(
                      ppr,
                      ppt,
                      pbs,
                      ptr,
                      ppm,
                      pbi,
                      index_r,
                      index_k1,
                      index_l2,
                      pwb->integral_splines[thread],
                      pwb->interpolated_integral[thread],
                      pwb->f[thread],
                      pwb),
            pbi->error_message,
            pbi->error_message);
      
          for (index_l3 = 0; index_l3 <= index_l2; ++index_l3) {  
  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k2) * l_size + index_l3) * k_tr_size]);

            /* Some debug - print transfer function */
            // if (index_r==0)
            //   for (index_k3=0; index_k3 < k_tr_size; ++index_k3)
            //     if ((index_k1==3) && (index_l2==2) && (index_l3==0))
            //       fprintf(stderr, "%g %g\n", k_tr[index_k3], transfer[index_k3]);

  
            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l3,
                          pwb->r[index_r],
                          &(pwb->integral_over_k2[index_l2][index_l3][index_r][index_k1]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);

  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k2;
  
            #pragma omp flush(abort)
  
          } // end of for(index_l3)
        } // end of for(index_l2)
      } // end of for(index_k1)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k2-integral array\n",
      pwb->count_memorised_for_integral_over_k2*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k2);
  
  /* Check that we correctly filled the array */
  class_test_permissive (pwb->count_memorised_for_integral_over_k2 != pwb->count_allocated_for_integral_over_k2,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k2, pwb->count_memorised_for_integral_over_k2);
  
  
  
  /* We can free the memory that was allocated for the I_l3 integral, as it is no longer needed */
  for (index_l3=0; index_l3 < pbi->l_size; ++index_l3) {
    for (index_r=0; index_r < pwb->r_size; ++index_r) {
      for (index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
        free (pwb->integral_over_k3[index_l3][index_r][index_k1]);
      } // end of for(index_k1)
      free (pwb->integral_over_k3[index_l3][index_r]);
    } // end of for(index_r)
    free (pwb->integral_over_k3[index_l3]);
  } // end of for(index_l3)
    
  free (pwb->integral_over_k3);
  
  
  
  return _SUCCESS_;
  
}








int bispectra_non_separable_interpolate_over_k2 (
    struct precision * ppr,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_k1,
    int index_l3,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  int index_k, index_k_tr, index_k2;

  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;


  /* So far, we always assumed that k1>=k2 because the shape function is symmetric wrt k1<->k2.
  Interpolating an array with such property is complicated, hence we build a temporary array
  where f(index_k1, index_k2) = f(index_k2, index_k1) when index_k1 < index_k2. */
  for (index_k2=0; index_k2 < k_pt_size; ++index_k2) {

    f[index_k2] = (index_k1 > index_k2 ? pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]:
                                         pwb->integral_over_k3[index_l3][index_r][index_k2][index_k1]);

    /* Multiply by window function */
    f[index_k2] *= pwb->k_window[index_k2];

  }
  
  
  if (ppr->transfers_k2_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider (desired size of the slow index) */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;

    /* Interpolate for each value of l3, r, k1 */
    if (ppr->transfers_k2_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k2_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1] + ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Revert the effect of the window function */
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  } // end of for (index_k_tr)


  /* Some debug - print the original array and the interpolation */
  // if ((index_k1==1) && (index_l3==0) && (index_r==0)) {
  // 
  //   fprintf (stderr, "\n\n");
  // 
  //   /* Node points before window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]/pwb->k_window[index_k]);
  // 
  //   fprintf (stderr, "\n");
  //   
  //   /* Node points after window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]);
  // 
  //   fprintf (stderr, "\n");
  // 
  //   /* Interpolation after inverse window */  
  //   for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  //  }

  return _SUCCESS_;

}






int bispectra_non_separable_integrate_over_k1 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_bt,
    int index_tt_k1,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Indices */
  int index_r;
  int index_k1, index_k2, index_k3;
  int index_l1, index_l2, index_l3, index_l1_l2_l3;
  
  /* Integration grid */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;


  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  
  // ==========================================================================================================
  // =                                    Allocate memory for INT_l1_l2_l3(r)                                 =
  // ==========================================================================================================
  
  /* The integral over k1 yields a function of (l1,l2,l3) that is symmetric under permutations of (l1,l2,l3).
  Hence, we compute and store it only for l1>=l2>=l3 configurations (that satisfy the triangular condition) */
  pwb->count_allocated_for_integral_over_k1 = pwb->r_size * pbi->n_independent_configurations;
  
  /* Allocate (l1,l2,l3)-level */
  class_alloc (pwb->integral_over_k1, pbi->n_independent_configurations*sizeof(double *), pbi->error_message);
  
  for (index_l1_l2_l3=0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
    class_alloc (pwb->integral_over_k1[index_l1_l2_l3], pwb->r_size*sizeof(double), pbi->error_message);
    
  if (pbi->bispectra_verbose > 1)
    printf("     * allocated ~ %.3g MB (%ld doubles) for the k1-integral array\n",
      pwb->count_allocated_for_integral_over_k1*sizeof(double)/1e6, pwb->count_allocated_for_integral_over_k1);
  
  
  
  
  // ==========================================================================================================
  // =                                    Compute  INT_l1_l2_l3(r)  integral                                  =
  // ==========================================================================================================
  
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k1 = 0;

  
  /* As for the other integrals, we parallelize the loop over 'r'. */
  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)       \
    private (index_r,index_l1_l2_l3,index_l1,index_l2,index_l3, thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      if (pbi->bispectra_verbose > 2)
        printf("     * computing the k1 integral for r=%g, index_r=%d\n", pwb->r[index_r], index_r);

      for (index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
    
        for (index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
          /* Skip those configurations that are forbidden by the triangular condition (optimization) */
          if (pbi->l[index_l2] < pbi->l[index_l1]/2)
            continue;

          /* Interpolate the integral I_l1_l2(k1,r) that we computed above in the integration grid of k1 */
          class_call_parallel (bispectra_non_separable_interpolate_over_k1 (
                        ppr,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_r,
                        index_l1,
                        index_l2,
                        pwb->integral_splines[thread],
                        pwb->interpolated_integral[thread],
                        pwb->f[thread],
                        pwb),
            pbi->error_message,
            pbi->error_message);      

  
          /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
          int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
          int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
          for (index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

            /* Index of the current (l1,l2,l3) configuration */
            index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k1) * l_size + index_l3) * k_tr_size]);
  
  
            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l3,
                          pwb->r[index_r],
                          &(pwb->integral_over_k1[index_l1_l2_l3][index_r]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
  
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k1;
              
            #pragma omp flush(abort)
  
          } // end of for(index_l3)
        } // end of for(index_l2)
      } // end of for(index_l1)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  if (pbi->bispectra_verbose > 2)
    printf(" -> memorised ~ %.3g MB (%ld doubles) for the k1-integral array\n",
      pwb->count_memorised_for_integral_over_k1*sizeof(double)/1e6, pwb->count_memorised_for_integral_over_k1);
  
  /* Check that we correctly filled the array */
  class_test_permissive (pwb->count_memorised_for_integral_over_k1 != pwb->count_allocated_for_integral_over_k1,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k1, pwb->count_memorised_for_integral_over_k1);
  
  /* We can free the memory that was allocated for the I_l2_l3 integral, as it is no longer needed */
  for (index_l2=0; index_l2<pbi->l_size; ++index_l2) {
    for (index_l3=0; index_l3<=index_l2; ++index_l3) {
      for (index_r=0; index_r < pwb->r_size; ++index_r) {      
        free (pwb->integral_over_k2[index_l2][index_l3][index_r]);
      } // end of for(index_r)
      free (pwb->integral_over_k2[index_l2][index_l3]);
    } // end of for(index_l3)
    free (pwb->integral_over_k2[index_l2]);
  } // end of for(index_l2)
    
  free (pwb->integral_over_k2);

  
  return _SUCCESS_;
  
}








int bispectra_non_separable_interpolate_over_k1 (
    struct precision * ppr,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    int index_r,
    int index_l1,
    int index_l2,
    double * integral_splines,
    double * interpolated_integral,
    double * f,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  int index_k, index_k_tr, index_k1;
  
  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  
  /* Define the function to be interpolated, and multiply it by a window function */
  for (index_k1=0; index_k1 < k_pt_size; ++index_k1) {
    
    f[index_k1] = pwb->integral_over_k2[index_l1][index_l2][index_r][index_k1];
    f[index_k1] *= pwb->k_window[index_k1];
  
  }

  
  if (ppr->transfers_k1_interpolation == cubic_interpolation) {
    
    class_call (array_spline_table_columns (
                  k_pt,
                  k_pt_size,
                  f,
                  1,  /* How many columns to consider (desired size of the slow index) */
                  integral_splines,
                  _SPLINE_EST_DERIV_,
                  pbi->error_message),
         pbi->error_message,
         pbi->error_message);
  }


  /* Interpolate at each k value using the usual spline interpolation algorithm */
  index_k = 0;
  double h = k_pt[index_k+1] - k_pt[index_k];
    
  for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr) {
    
    while (((index_k+1) < k_pt_size) && (k_pt[index_k+1] < k_tr[index_k_tr])) {
      index_k++;
      h = k_pt[index_k+1] - k_pt[index_k];
    }
    
    class_test(h==0., pbi->error_message, "stop to avoid division by zero");
    
    double b = (k_tr[index_k_tr] - k_pt[index_k])/h;
    double a = 1.-b;
      
    /* Interpolate for each value of l3, r, k1 */
    if (ppr->transfers_k1_interpolation == linear_interpolation) {
      interpolated_integral[index_k_tr] = a * f[index_k] + b * f[index_k+1];
    }
    else if (ppr->transfers_k1_interpolation == cubic_interpolation) {
      interpolated_integral[index_k_tr] =  
        a * f[index_k] + b * f[index_k+1] + ((a*a*a-a) * integral_splines[index_k] +(b*b*b-b) * integral_splines[index_k+1])*h*h/6.0;
    }

    /* Revert the effect of the window function */
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  } // end of for (index_k_tr)


  /* Some debug - print the original array and the interpolation */
  // if ((index_l1==48) && (index_l2==35) && (index_r==50)) {
  // 
  //   fprintf (stderr, "\n\n");
  // 
  //   /* Node points before window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]/pwb->k_window[index_k]);
  // 
  //   fprintf (stderr, "\n");
  //   
  //   /* Node points after window */
  //   for (index_k=0; index_k < k_pt_size; ++index_k)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_pt[index_k], f[index_k]);
  // 
  //   fprintf (stderr, "\n");
  // 
  //   /* Interpolation after inverse window */  
  //   for (index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
  //     fprintf (stderr, "%17.7g %17.7g\n", k_tr[index_k_tr], interpolated_integral[index_k_tr]);
  // 
  //   fprintf (stderr, "\n\n");
  //   
  //  }

  return _SUCCESS_;
  
}









int bispectra_non_separable_integrate_over_r (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    double * bispectrum,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* We now proceed to the final step of the bispectrum computation where we integrate over
  the r-dependence and obtain the bispectrum. We also multiply the result by the appropriate
  coefficients.  */  
    
  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  /* We parallelize the outer loop over 'l1'. */
  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)
  {
  
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      if (pbi->bispectra_verbose > 1)
        printf("     * computing the r-integral for l1=%d, index_l1=%d\n", pbi->l[index_l1], index_l1);
    
      for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
  
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (pbi->l[index_l2] < pbi->l[index_l1]/2)
          continue;

        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

          long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];    
          double * I = pwb->integral_over_k1[index_l1_l2_l3];
          double integral = 0;

          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
            double r = pwb->r[index_r];
            double integrand = r*r * I[index_r];
      
            if (integrand == 0.)
              continue;
      
            integral += integrand * pwb->delta_r[index_r];

            /* Some debug - output intermediate results on stderr for a custom (l1,l2,l3) configuration */
            // if ( (l1==100) && (l2==200) && (l3==300) ) {
            //   if (index_r==0) {
            //     fprintf(stderr, "##########    l1 = %d, l2 = %d, l3 = %d, n_rows = %d    ##########\n\n",
            //       l1, l2, l3, r, pwb->r_size);
            //     fprintf(stderr, "%12s %17s %17s %17s\n", "r", "r^2*I_l1_l2_l3(r)", "integral","delta_r");
            //   }
            //   else {
            //     fprintf(stderr, "%12.7g %17.7g %17.7g %17.7g\n", r, integrand, integral, pwb->delta_r[index_r]);
            //   }
            // }
   
          } // end of for(index_r)

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3) with 1/2 from trapezoidal rule */
          bispectrum[index_l1_l2_l3] = 0.5 * integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. This factor is seen
            (see, for instance, eq. 17 of Fergusson & Shellard 2007). */
          bispectrum[index_l1_l2_l3] *= pow(2./_PI_,3);

          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_for_bispectra++;

          /* Some debug - output the integral as a function of r on stderr for a custom (l1,l2,l3) */
          // if ( (l1==l2) && (l2==l3) ) {
          //   fprintf(stderr, "%12d %17.7g\n", l1, pwb->integral_over_r[index_l1][index_l2][index_l3-index_l3_min]);
          // }

        } // end of for(index_l3)
      } // end of for(index_l2)
      
      #pragma omp flush(abort)
      
    } // end of for(index_l1)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  /* We can free the memory that was allocated for the I_l1_l2_l3(r) integral, as it is no longer needed */
  for (long int index_l1_l2_l3 = 0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
    free (pwb->integral_over_k1[index_l1_l2_l3]);

  free (pwb->integral_over_k1);
    
  return _SUCCESS_;
  
}




int bispectra_non_separable_workspace_free (
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  free (pwb->r);
  free (pwb->delta_r);  
  
  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  // int i;
  // for (i=0; i < 200 ; ++i) {
  //   printf ("pwb->k3_grid[0][i] = %g\n", pwb->k3_grid[0][i]);
  // }
  
  #pragma omp parallel private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    free(pwb->k3_grid[thread]);
    free(pwb->delta_k3[thread]);
    free(pwb->integral_splines[thread]);
    free(pwb->interpolated_integral[thread]);
    free(pwb->f[thread]);
    
  }  if (abort == _TRUE_) return _FAILURE_;
  
  free(pwb->k3_grid);
  free(pwb->delta_k3);
  free(pwb->integral_splines);
  free(pwb->interpolated_integral);
  free(pwb->f);
  
  return _SUCCESS_;
  
}





/**
 * Bispectrum produced by pi_dot * grad_pi^2 in the Galileon Lagrangian, from eq. 22 of arXiv:0905.3746v3
 *
 */
int bispectra_galileon_gradient (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3,
  double pk_1, double pk_2, double pk_3,
  double * out
  )
{ 
  
  double K1 = k1 + k2 + k3;
  double K1_squared = K1*K1;
  double K1_cubed = K1_squared*K1;
  double K1_fourth = K1_squared*K1_squared;
  double K1_sixth = K1_fourth*K1_squared;
  
  double K2_squared = k1*k2 + k2*k3 + k3*k1;
  double K2_fourth = K2_squared*K2_squared;
  
  double K3_cubed = k1*k2*k3;
  double K3_sixth = K3_cubed*K3_cubed;
  double K3_ninth = K3_sixth*K3_cubed;
  
  *out = -27/17. *

      (

        + 24*K3_sixth - 8*K2_squared*K3_cubed*K1 - 8*K2_fourth*K1_squared
        + 22*K3_cubed*K1_cubed - 6*K2_squared*K1_fourth + 2*K1_sixth

      ) / (K3_ninth*K1_cubed);
  
  /* Multiply by the amplitude of primordial fluctuations. The factor -3/5 derives by the
  fact that we defined the transfer functions wrt R, and that a bispectrum with fnl_phi=1
  (which is what we want) is equivalent to a bispectrum with fnl_R=-3/5.  */
  double fnl_R = -3/5.; 
  *out *= fnl_R * (2*_PI_*_PI_*ppm->A_s) * (2*_PI_*_PI_*ppm->A_s);
  
  return _SUCCESS_; 
  
}



/**
 * Bispectrum produced by pi_dot^3 in the Galileon Lagrangian, from eq. 22 of arXiv:0905.3746v3
 *
 */
int bispectra_galileon_time (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double K1 = k1+k2+k3;
  double K1_cubed = K1*K1*K1;
  double K3_cubed = k1*k2*k3;
  
  *out = 162 / (K3_cubed * K1_cubed);
  
  /* Multiply by the amplitude of primordial fluctuations */
  double fnl_R = -3/5.; 
  *out *= fnl_R * (2*_PI_*_PI_*ppm->A_s) * (2*_PI_*_PI_*ppm->A_s);
  
  return _SUCCESS_; 
  
}



/**
 * Local template bispectrum. Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_local_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  *out = 2 * (-3/5.) * ( pk_1*pk_2 + pk_1*pk_3 + pk_2*pk_3 );
  
  return _SUCCESS_; 
  
}


/**
 * Equilateral template bispectrum. Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_equilateral_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double pk_one_third_1 = pow(pk_1, _ONE_THIRD_);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, _ONE_THIRD_);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, _ONE_THIRD_);
  double pk_two_thirds_3 = pk_one_third_3*pk_one_third_3;

  
  *out = 6 * (-3/5.) * (
      
      - pk_1 * pk_2
      - pk_1 * pk_3
      - pk_2 * pk_3
        
      - 2 * pk_two_thirds_1 * pk_two_thirds_2 * pk_two_thirds_3
        
      + pk_1 * pk_one_third_2 * pk_two_thirds_3
      + pk_1 * pk_one_third_3 * pk_two_thirds_2
      + pk_2 * pk_one_third_1 * pk_two_thirds_3
      + pk_2 * pk_one_third_3 * pk_two_thirds_1
      + pk_3 * pk_one_third_1 * pk_two_thirds_2
      + pk_3 * pk_one_third_2 * pk_two_thirds_1
        
  );
  

  return _SUCCESS_; 
  
}
 
 
   
/**
 * Orthogonal template bispectrum (arXiv:0905.3746). Note that this is debug only, as all separable bispectra
 * are computed more efficiently in separable_init.
 *
 */
int bispectra_orthogonal_model (
  struct primordial * ppm,
  struct bispectra * pbi,
  double k1, double k2, double k3, 
  double pk_1, double pk_2, double pk_3, 
  double * out
  )
{ 
  
  double pk_one_third_1 = pow(pk_1, _ONE_THIRD_);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, _ONE_THIRD_);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, _ONE_THIRD_);
  double pk_two_thirds_3 = pk_one_third_3*pk_one_third_3;

  
  *out = 6 * (-3/5.) * (
      
      - 3 * pk_1 * pk_2
      - 3 * pk_1 * pk_3
      - 3 * pk_2 * pk_3
        
      - 8 * pk_two_thirds_1 * pk_two_thirds_2 * pk_two_thirds_3
        
      + 3 * pk_1 * pk_one_third_2 * pk_two_thirds_3
      + 3 * pk_1 * pk_one_third_3 * pk_two_thirds_2
      + 3 * pk_2 * pk_one_third_1 * pk_two_thirds_3
      + 3 * pk_2 * pk_one_third_3 * pk_two_thirds_1
      + 3 * pk_3 * pk_one_third_1 * pk_two_thirds_2
      + 3 * pk_3 * pk_one_third_2 * pk_two_thirds_1
        
  );
  

  return _SUCCESS_; 
  
}








