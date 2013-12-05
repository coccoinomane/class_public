/** @file fisher.c documented Fisher matrix module for second-order perturbations
 *
 * Guido W Pettinari, 19.07.2012
 *
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
  if (pfi->has_fisher == _FALSE_) {
  
    if (pfi->fisher_verbose > 0)
      printf("No forecasts requested. Fisher module skipped.\n");
  
    return _SUCCESS_;
  }
  else {
    if (pfi->fisher_verbose > 0)
      printf("Computing fisher matrix for %d bispectra\n", pbi->bt_size);
  }


  // =========================================================
  // =                    Preparations                       =
  // =========================================================
  
  /* Initialize indices & arrays in the Fisher structure */
  
  class_call (fisher_indices(ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
    

  /* Load intrinsic bispectra from disk if they were not already computed. Note that we are only loading the
  intrinsic bispectra, since the primary ones are very quick to recompute. */

  if (pbi->load_bispectra_from_disk == _TRUE_) {
    
    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt)   
      if ((pbi->bispectrum_type[index_bt]==non_separable_bispectrum)
       || (pbi->bispectrum_type[index_bt]==intrinsic_bispectrum))
        class_call (bispectra_load_from_disk (
                      pbi,
                      index_bt),
          pfi->error_message,
          pfi->error_message);
  }


  // ======================================================================
  // =                   Compute noise power spectrum                     =
  // ======================================================================
  
  class_call (fisher_noise(ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
  
  

  // ====================================================================
  // =                      Compute cross C_l's                         =
  // ====================================================================
    
  class_call (fisher_cross_cls(ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi,pfi),
    pfi->error_message,
    pfi->error_message);
  
  
  
  // =====================================================================
  // =                      Prepare interpolation mesh                   =
  // =====================================================================
  
  
  if ((pfi->bispectra_interpolation == mesh_interpolation) || (pfi->bispectra_interpolation == mesh_interpolation_2d))
    
    class_call (fisher_create_interpolation_mesh(ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi,pfi),
      pfi->error_message,
      pfi->error_message);
  
  
  
  // ====================================================================
  // =                       Compute Fisher matrix                      =
  // ====================================================================
  
  class_call (fisher_compute (ppr,pba,ppt,pbs,ptr,ppm,psp,ple,pbi,pfi),
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
     struct perturbs * ppt,
     struct bispectra * pbi,
     struct fisher * pfi
     )
{

  if (pfi->has_fisher == _TRUE_) {

    for(int index_l=0; index_l<pfi->l1_size; ++index_l) {

      for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
        free (pfi->fisher_matrix[index_l][index_bt]);
        free (pfi->fisher_matrix_of_l1[index_l][index_bt]);
        free (pfi->inverse_fisher_matrix[index_l][index_bt]);
      }

      free (pfi->fisher_matrix[index_l]);
      free (pfi->fisher_matrix_of_l1[index_l]);
      free (pfi->inverse_fisher_matrix[index_l]);
      free (pfi->sigma_fnl[index_l]);

    } // end of for(index_l)
    
    free (pfi->fisher_matrix);
    free (pfi->fisher_matrix_of_l1);
    free (pfi->inverse_fisher_matrix);
    free (pfi->sigma_fnl);

    free (pfi->l1);
    free (pfi->l2);
    free (pfi->l3);
    

    /* Free 3j symbols */
    if ((pfi->bispectra_interpolation != mesh_interpolation)
       && (pfi->bispectra_interpolation != mesh_interpolation_2d)) 
      free (pfi->I_l1_l2_l3);


    /* Free meshes. Note that we start from the last bispectrum as the grid, which is shared between
    the various meshes, belongs to the index_bt=0 bispectrum, and hence should be freed last.  */
    if ((pfi->bispectra_interpolation == mesh_interpolation)
     || (pfi->bispectra_interpolation == mesh_interpolation_2d)) {
         
      for (int index_bt=(pbi->bt_size-1); index_bt >= 0; --index_bt) {
        for (int i = (pbi->bf_size-1); i >= 0; --i) {
          for (int j = (pbi->bf_size-1); j >= 0; --j) {
            for (int k = (pbi->bf_size-1); k >= 0; --k) {
              for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
     
                if ((index_mesh == 0) && (pfi->l_turnover[0] <= pbi->l[0]))
                  continue;
     
                if ((index_mesh == 1) && (pfi->l_turnover[0] > pbi->l[pbi->l_size-1]))
                  continue;
     
                class_call (mesh_free (pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]),
                  pfi->error_message,
                  pfi->error_message);      
              }
              free (pfi->mesh_workspaces[index_bt][i][j][k]);
            }
            free (pfi->mesh_workspaces[index_bt][i][j]);
          }
          free (pfi->mesh_workspaces[index_bt][i]);
        }
        free (pfi->mesh_workspaces[index_bt]);
      }
      free (pfi->mesh_workspaces);
    }  // end of if(mesh)
    
    /* Free cross-power spectrum */
    for (int index_l=0; index_l < pbi->l_size; ++index_l) {
      for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
        free (pfi->C[index_l][index_bf]);
        free (pfi->inverse_C[index_l][index_bf]);
      }
      free (pfi->C[index_l]);
      free (pfi->inverse_C[index_l]);
    }
    free (pfi->C);
    free (pfi->inverse_C);

    /* Free noise spectrum */
    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf)
      free(pfi->N_l[index_bf]);
    free(pfi->N_l);

    
  } // end of if(has_fisher)
  
  return _SUCCESS_;
 
}




/**
 * This routine defines indices and allocates tables in the Fisher structure 
 *
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
  // =                                Determine l-sampling                              =
  // ====================================================================================
  
  /* Here we choose the 3D grid in (l1,l2,l3) over which to sum the estimator. For the mesh
  interpolation, where we interpolate the bispectra for all configurations, we need the full
  grid. For trilinear interpolation, where we sum over the support points, we only need 
  the l's where we computed the bispectra.  */

  pfi->l_min = 2;
  pfi->l_max = pbi->l[pbi->l_size-1];
  pfi->full_l_size = pfi->l_max - pfi->l_min + 1;

  if (pfi->bispectra_interpolation == mesh_interpolation) {

    pfi->l1_size = pfi->l2_size = pfi->l3_size = pfi->full_l_size;

    class_alloc (pfi->l1, pfi->full_l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l2, pfi->full_l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l3, pfi->full_l_size*sizeof(int), pfi->error_message);
    
    for (int index_l=0; index_l < pfi->full_l_size; ++index_l) {
      pfi->l1[index_l] = pfi->l_min + index_l;
      pfi->l2[index_l] = pfi->l_min + index_l;
      pfi->l3[index_l] = pfi->l_min + index_l;
    }    
  }
  /* mesh_2d means that we fully interpolate in the l2 and l3 level, but we only sum over
  the support points for the l1 level. */
  else if (pfi->bispectra_interpolation == mesh_interpolation_2d) {

    pfi->l1_size = pbi->l_size;
    pfi->l2_size = pfi->l3_size = pfi->full_l_size;

    class_alloc (pfi->l1, pbi->l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l2, pfi->full_l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l3, pfi->full_l_size*sizeof(int), pfi->error_message);
    
    for (int index_l=0; index_l < pbi->l_size; ++index_l) {
      pfi->l1[index_l] = pbi->l[index_l];
    }
    
    for (int index_l=0; index_l < pfi->full_l_size; ++index_l) {
      pfi->l2[index_l] = pfi->l_min + index_l;
      pfi->l3[index_l] = pfi->l_min + index_l;
    }
  }
  /* Any other type of interpolation sums over the node points only */
  else {

    pfi->l1_size = pfi->l2_size = pfi->l3_size = pbi->l_size;

    class_alloc (pfi->l1, pbi->l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l2, pbi->l_size*sizeof(int), pfi->error_message);
    class_alloc (pfi->l3, pbi->l_size*sizeof(int), pfi->error_message);
    
    for (int index_l=0; index_l < pbi->l_size; ++index_l) {
      pfi->l1[index_l] = pbi->l[index_l];
      pfi->l2[index_l] = pbi->l[index_l];
      pfi->l3[index_l] = pbi->l[index_l];
    }
  }


  /* Allocate arrays for the noise in the C_l's */
  class_alloc (pfi->N_l, pbi->bf_size*sizeof(double *), pfi->error_message);
  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf)
    class_calloc (pfi->N_l[index_bf], pfi->full_l_size, sizeof(double), pfi->error_message);


  /* Allocate arrays needed to perform the mesh interpolation */
  class_alloc (pfi->fisher_matrix, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->fisher_matrix_of_l1, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->inverse_fisher_matrix, pfi->l1_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->sigma_fnl, pfi->l1_size*sizeof(double *), pfi->error_message);
  
  for (int index_l=0; index_l < pfi->l1_size; ++index_l) {

    class_alloc (pfi->fisher_matrix[index_l], pbi->bt_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->fisher_matrix_of_l1[index_l], pbi->bt_size*sizeof(double *), pfi->error_message);
    class_alloc (pfi->inverse_fisher_matrix[index_l], pbi->bt_size*sizeof(double *), pfi->error_message);
    class_calloc (pfi->sigma_fnl[index_l], pbi->bt_size, sizeof(double), pfi->error_message);
    
    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
      class_calloc (pfi->fisher_matrix[index_l][index_bt], pbi->bt_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->fisher_matrix_of_l1[index_l][index_bt], pbi->bt_size, sizeof(double), pfi->error_message);
      class_calloc (pfi->inverse_fisher_matrix[index_l][index_bt], pbi->bt_size, sizeof(double), pfi->error_message);
    }
    
  }


  // =========================================================================================
  // =                                 Interpolation arrays                                  = 
  // =========================================================================================
  
  if ((pfi->bispectra_interpolation == mesh_interpolation)
   || (pfi->bispectra_interpolation == mesh_interpolation_2d)) {
  
    /* TODO: generalize this, or modify the mesh_sort/mesh_int functions */
    pfi->n_meshes = 2;
    
    /* Allocate arrays that will contain the parameters for the mesh interpolation */
    class_alloc (pfi->link_lengths, pfi->n_meshes*sizeof(double), pfi->error_message);
    class_alloc (pfi->group_lengths, pfi->n_meshes*sizeof(double), pfi->error_message);
    class_alloc (pfi->soft_coeffs, pfi->n_meshes*sizeof(double), pfi->error_message);
  
    /* Turnover point between the two meshes */
    class_alloc (pfi->l_turnover, (pfi->n_meshes-1)*sizeof(int), pfi->error_message);

    /* Allocate one mesh interpolation workspace per type of bispectrum */
    class_alloc (pfi->mesh_workspaces, pbi->bt_size*sizeof(struct mesh_interpolation_workspace *****), pfi->error_message);
  
    /* Allocate worspaces and intialize counters */
    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
      
      class_alloc (pfi->mesh_workspaces[index_bt], pbi->bf_size*sizeof(struct mesh_interpolation_workspace ****), pfi->error_message);

        for (int i = 0; i < pbi->bf_size; ++i) {

          class_alloc (pfi->mesh_workspaces[index_bt][i], pbi->bf_size*sizeof(struct mesh_interpolation_workspace ***), pfi->error_message);

          for (int j = 0; j < pbi->bf_size; ++j) {
            
            class_alloc (pfi->mesh_workspaces[index_bt][i][j], pbi->bf_size*sizeof(struct mesh_interpolation_workspace **), pfi->error_message);

            for (int k = 0; k < pbi->bf_size; ++k) {

              class_alloc (pfi->mesh_workspaces[index_bt][i][j][k], pfi->n_meshes*sizeof(struct mesh_interpolation_workspace *), pfi->error_message);
      
              for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
        
                class_alloc (pfi->mesh_workspaces[index_bt][i][j][k][index_mesh], sizeof(struct mesh_interpolation_workspace), pfi->error_message);
                pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->n_allocated_in_grid = 0;
                pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->n_allocated_in_mesh = 0; 
                pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->count_interpolations = 0;
                pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->count_range_extensions = 0;
        
            } // end of for(index_mesh)
          } // end of for(i)
        } // end of for(j)
      } // end of for(k)
    } // end of for(index_bt)
  } // end of if(mesh_interpolation)
  


  // =================================================================================================
  // =                                    Compute three_j symbol                                     =
  // =================================================================================================
  
  /* In the following lines, we shall allocate and fill the (l1,l2,l3) array containing the quantity
    I_l1_l2_l3 in eq. 13 of Komatsu, Spergel & Wandelt (2005):
  
     I_l1_l2_l3 = sqrt( (2L1+1)(2L2+1)(2L3+1)/4*pi ) *  (L1 L2 L3)
                                                        (0  0  0 )
  
    which is needed to compute the angle-averaged bispectra, which in turn are used by the f_NL
    estimator. */
  
  /* For mesh interpolation, we shall compute the 3j's directly in the estimator */
  if ((pfi->bispectra_interpolation != mesh_interpolation) 
   && (pfi->bispectra_interpolation != mesh_interpolation_2d)) {
  
    if (pfi->fisher_verbose > 1)
      printf(" -> allocating and computing ~ %.3g MB (%ld doubles) of 3j-symbols\n",
        pbi->n_independent_configurations*sizeof(double)/1e6, pbi->n_independent_configurations);
  
    /* Allocate memory for pfi->I_l1_l2_l3 */
    class_alloc (pfi->I_l1_l2_l3, pbi->n_independent_configurations*sizeof(double), pfi->error_message);
  
    /* Temporary values needed for the computation of the 3j symbol. The temporary array will contain at most 2*l_max+1
      values when l1=l2=l_max, with the +1 accounting for the l=0 case. */
    double * temporary_three_j;
    class_alloc (temporary_three_j, (2*pfi->l_max+1)*sizeof(double), pfi->error_message);
    int l3_min, l3_max, temporary_three_j_size;
             
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
        int index_l3_max = min (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
    
        /* Compute the actual 3j symbol for all allowed values of l3 */
        class_call (threej_l1 (
                      l1, l2, 0, 0,
                      &l3_min, &l3_max,
                      &temporary_three_j,
                      &temporary_three_j_size,
                      pfi->error_message       
                      ),
          pfi->error_message,
          pfi->error_message);
    
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
  
        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_l1)
  
    free (temporary_three_j);
  
  } // end of if not mesh interpolation


  // ==================================================================================================
  // =                                      Compute wasted points                                     =
  // ==================================================================================================

  /* Count the number of bispectra configurations that won't be used for the interpolation */
  long int count_wasted = 0;

  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {
    int l1 = pbi->l[index_l1];
    for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
      int l2 = pbi->l[index_l2];
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = min (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
        if ((l1+l2+pbi->l[index_l3])%2!=0)
          count_wasted++;
    }
  }

  if (pfi->fisher_verbose > 1)
    printf ("     * %.3g%% of the computed bispectra configurations will be (directly) used to compute the Fisher matrix (%d wasted)\n",
    100-count_wasted/(double)(pbi->n_independent_configurations)*100, count_wasted);
  
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

  double noise[pbi->bf_size][pfi->n_channels];
  
  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
    for (int index_channel=0; index_channel < pfi->n_channels; ++index_channel) {
      if ((pbi->has_bispectra_t == _TRUE_) && (index_bf == pbi->index_bf_t))
        noise[index_bf][index_channel] = pfi->noise_t[index_channel];

      if ((pbi->has_bispectra_e == _TRUE_) && (index_bf == pbi->index_bf_e))
        noise[index_bf][index_channel] = pfi->noise_e[index_channel];
    }
  }


  // ===================================================================================
  // =                                Print information                                =
  // ===================================================================================

  if (pfi->fisher_verbose > 0) {

    printf (" -> noise model: beam_fwhm=");
    for (int index_channel=0; index_channel < pfi->n_channels-1; ++index_channel)
      printf ("%g,", pfi->beam[index_channel] * 60 * 180./_PI_);
    printf ("%g", pfi->beam[pfi->n_channels-1] * 60 * 180./_PI_);
    printf (" (arcmin)");
    
    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

      printf (", sigma_%s=", pbi->bf_labels[index_bf]);
      for (int index_channel=0; index_channel < pfi->n_channels-1; ++index_channel)
        printf ("%g,", sqrt(noise[index_bf][index_channel]) * 1e6*pba->T_cmb / pfi->beam[index_channel]);
      printf ("%g", sqrt(noise[index_bf][pfi->n_channels-1]) * 1e6*pba->T_cmb / pfi->beam[pfi->n_channels-1]);
      printf (" (uK)");
    }
    printf ("\n");
  }



  // ===================================================================================
  // =                             Co-added Gaussian noise                             =
  // ===================================================================================
  
  /* Eeach frequency channel contributes to the inverse noise with a 1/N contribution */
  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
    for (int l=pfi->l_min; l <= pfi->l_max; ++l) {
      for (int index_channel=0; index_channel < pfi->n_channels; ++index_channel) {

        double beam_exponent = (l*(l+1.)*pfi->beam[index_channel]*pfi->beam[index_channel])/(8.*log(2.));

        class_test (beam_exponent > DBL_MAX_EXP,
          pfi->error_message,
          "stop to prevent overflow. Reduce experiment_beam_fwhm to something more reasonable, or set it to zero.")
        
        /* Inverse contribution to l from the considered channel */
        pfi->N_l[index_bf][l-2] += 1. / (noise[index_bf][index_channel] * exp(beam_exponent));
      }
    
      /* Invert back */
      pfi->N_l[index_bf][l-2] = 1. / pfi->N_l[index_bf][l-2];
    
      /* Some debug */
      // double factor = 1e12*pba->T_cmb*pba->T_cmb*l*(l+1)/(2*_PI_);
      // fprintf (stderr, "%d %g %g\n", l, factor*pfi->N_l[pbi->index_bf_t][l-2]);
    
    } // end of for(l=2)
  } // end of for(T,E,...)
  

  return _SUCCESS_;
  
}



/**
 * Build the mesh for the interpolation of the bispectra.
 *
 */
int fisher_create_interpolation_mesh(
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
  

  if (pfi->fisher_verbose > 1)
    printf (" -> preparing the interpolation mesh for %d bispectr%s\n",
    pbi->bt_size*pbi->n_probes, ((pbi->bt_size*pbi->n_probes)!=1?"a":"um"));


  // =====================================================
  // =              Interpolation parameters             =
  // =====================================================
  
  // ***  Fine mesh parameters  ***
        
  pfi->link_lengths[0] = 2 * (pbi->l[1] - pbi->l[0]);
  pfi->group_lengths[0] = 0.1 * (pbi->l[1] - pbi->l[0]);
  pfi->soft_coeffs[0] = 0.5;  

  // ***  Coarse mesh parameters  ***

  pfi->link_lengths[1] = 0.5/sqrt(2) * ppr->l_linstep;
  pfi->group_lengths[1] = 0.1 * ppr->l_linstep;
  pfi->soft_coeffs[1] = 0.5;
  
  for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh)
    class_test (pfi->link_lengths[index_mesh] <= pfi->group_lengths[index_mesh],
      pfi->error_message,
      "the linking length must be larger than the grouping length.")
  
  
  

  // ==============================================
  // =       Determine the turnover points        =
  // ==============================================
  
  /* Never use the fine grid if the large step is used since the beginning */
  if ((pbi->l[1]-pbi->l[0]) >= ppr->l_linstep) {
    pfi->l_turnover[0] = pbi->l[0];
  }
  /* Never use the coarse grid if the linstep (i.e. the largest possible step) is never used.
  The max() is needed because the last point is fixed and does not depend on the grid spacing. */
  else if (max(pbi->l[pbi->l_size-1]-pbi->l[pbi->l_size-2], pbi->l[pbi->l_size-2]-pbi->l[pbi->l_size-3]) < ppr->l_linstep) {
    pfi->l_turnover[0] = pbi->l[pbi->l_size-1] + 1;
  }
  else {
    int index_l = 0;
    while ((index_l < pbi->l_size-1) && ((pbi->l[index_l+1] - pbi->l[index_l]) < ppr->l_linstep))
      index_l++;
    pfi->l_turnover[0] = pbi->l[index_l-1];
  }


  if (pfi->fisher_verbose > 1)
    printf ("     * l_turnover=%d, n_boxes=[%d,%d], linking lengths=[%g,%g], grouping lengths=[%g,%g]\n",
      pfi->l_turnover[0],
      (long int)ceil(pfi->l_turnover[0]/ (pfi->link_lengths[0]*(1.+pfi->soft_coeffs[0]))),
      (long int)ceil(pbi->l[pbi->l_size-1] / (pfi->link_lengths[1]*(1.+pfi->soft_coeffs[1]))),
      pfi->link_lengths[0], pfi->link_lengths[1],
      pfi->group_lengths[0], pfi->group_lengths[1]);
    




  // ===============================================
  // =                 Create meshes               =
  // ===============================================

  /* TODO: parallelize the filling of values, it should improve speed when using a very
  fine l-grid */

  if (pfi->fisher_verbose > 1)
    printf ("     * allocating memory for values...\n");

  /* Allocate array that will contain the values of the bispectra arranged in the right order for mesh sorting */
  double ** values;
  class_alloc (values, pbi->n_independent_configurations*sizeof(double *), pfi->error_message);

  for (long int index_l1_l2_l3=0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
    class_alloc (values[index_l1_l2_l3], 4*sizeof(double), pfi->error_message);

  if (pfi->fisher_verbose > 1)
    printf ("     * rearranging values of the bispectra...\n");  
  
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    for (int i = 0; i < pbi->bf_size; ++i) {
    for (int j = 0; j < pbi->bf_size; ++j) {
    for (int k = 0; k < pbi->bf_size; ++k) {

      // -------------------------------------------------------------
      // -                  Rearrange the bispectra                  -
      // -------------------------------------------------------------

      for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

        int l1 = pbi->l[index_l1];
        /* TODO: choose Cl's according to the considered probe (ie not only TTT) */
        double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
        
        for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
          
          int l2 = pbi->l[index_l2];
          double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];

          /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
          int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
          int index_l3_max = min (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
  
          // ***  Rearrange the bispectra into 'values' ***
          for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

            int l3 = pbi->l[index_l3];
            /* TODO: choose Cl's according to the considered probe (ie not only TTT) */
            double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];

            /* Index of the current (l1,l2,l3) configuration */
            long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
  
            /* The 0th element of 'values' is the function to be interpolated. The natural scaling for the bispectrum
            (both the templates and the second-order one) is given by Cl1*Cl2 + Cl2*Cl3 + Cl1*Cl3, which we adopt
            as a window function. */
            double window = 1./(C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3);
            values[index_l1_l2_l3][0] = window * pbi->bispectra[index_bt][i][j][k][index_l1_l2_l3];

            /* 1st to 3rd arguments are the coordinates */
            values[index_l1_l2_l3][1] = (double)(l1);
            values[index_l1_l2_l3][2] = (double)(l2);
            values[index_l1_l2_l3][3] = (double)(l3);
            
          } // end of for(index_l3)
        } // end of for(index_l2)
      } // end of for(index_l1)

      // -------------------------------------------------------------
      // -                    Generate the meshes                    -
      // -------------------------------------------------------------      

      for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
  
        /* Set the maximum l. The fine grid does not need to go up to l_max */
        if (index_mesh == 0)
          pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->l_max = (double)pfi->l_turnover[0];
        else 
          pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->l_max = (double)pbi->l[pbi->l_size-1];

        pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->n_points = pbi->n_independent_configurations;
        pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->link_length = pfi->link_lengths[index_mesh];
        pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->group_length = pfi->group_lengths[index_mesh];
        pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->soft_coeff = pfi->soft_coeffs[index_mesh];

        /* Since the grid is shared between different bispectra, it needs to be computed only for the first one */
        if ((index_bt == 0) && (i == 0) && (j == 0) && (k == 0)) {
          pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->compute_grid = _TRUE_;
        }
        else {
          pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->compute_grid = _FALSE_;
          pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->grid = pfi->mesh_workspaces[0][0][0][0][index_mesh]->grid;
        }
  
        if (pfi->fisher_verbose > 2) 
          printf ("     * computing mesh for bispectrum %s_%s(%d)\n",
          pbi->bt_labels[index_bt], pbi->bfff_labels[i][j][k], index_mesh);
              
        /* Skip the fine grid if the turnover point is smaller than than the smallest multipole */
        if ((index_mesh == 0) && (pfi->l_turnover[0] <= pbi->l[0])) {
          if (pfi->fisher_verbose > 2)
            printf ("      \\ fine grid not needed because l_turnover <= l_min (%d <= %d)\n", pfi->l_turnover[0], pbi->l[0]);
          continue;
        }
  
        /* Skip the coarse grid if the turnover point is larger than than the largest multipole */
        if ((index_mesh == 1) && (pfi->l_turnover[0] > pbi->l[pbi->l_size-1])) {
          if (pfi->fisher_verbose > 2)
            printf ("      \\ coarse grid not needed because l_turnover > l_max (%d > %d)\n", pfi->l_turnover[0], pbi->l[pbi->l_size-1]);
          continue;
        }
            
        /* Generate the mesh */
        class_call (mesh_sort (
                      pfi->mesh_workspaces[index_bt][i][j][k][index_mesh],
                      values),
          pfi->error_message,
          pfi->error_message);
  
        if (pfi->fisher_verbose > 2)
          printf ("      \\ allocated (grid,mesh)=(%g,%g) MBs\n",
            pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->n_allocated_in_grid*8/1e6,
            pfi->mesh_workspaces[index_bt][i][j][k][index_mesh]->n_allocated_in_mesh*8/1e6);
  
      } // end of for(index_mesh)
            
    }}} // end of for(ijk)
  }  // end of for(index_bt)
  
  /* We do not need the node values anymore */
  for (long int index_l1_l2_l3=0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
    free (values[index_l1_l2_l3]);
  
  free (values);
  
  /* Check whether we recover the node point b(2,2,2) */
  // double interpolated_bispectrum;
  // 
  // class_call (mesh_int(pbi, 0, 2., 2., 2., &interpolated_bispectrum), pfi->error_message, pfi->error_message);
  // 
  // 
  // 
  // printf ("%g %g\n", pbi->bispectrum_mesh[0][0][0][0][0][0], interpolated_bispectrum);
  
  
  /* Print a slice of interpolated bispectrum */
  // int index_l1_debug = 50;
  // int index_l2_debug = 50;
  // int index_l3;
  // 
  // int index_l_triangular_min = pbi->index_l_triangular_min[index_l1_debug][index_l2_debug];
  // int l_triangular_size = pbi->l_triangular_size[index_l1_debug][index_l2_debug];
  // 
  // double l1 = pbi->l[index_l1_debug];
  // double l2 = pbi->l[index_l2_debug];
  // double l3;
  // 
  // fprintf (stderr, "### l1 = %g, l2 = %g\n", l1, l2);
  // printf ("### l1 = %g, l2 = %g\n", l1, l2);
      
  /* Node values */  
  // for (index_l3=index_l_triangular_min; index_l3<(index_l_triangular_min + l_triangular_size); ++index_l3) {  
  //   
  //   l3 = pbi->l[index_l3];
  //   double interpolated_bispectrum;
  //   class_call (bispectra_interpolate (pbi, 0, l1, l2, l3, &interpolated_bispectrum), pfi->error_message, pfi->error_message);
  // 
  //   double C_l1 = psp->cl[pwf->index_mode][(index_l1_debug * psp->ic_ic_size[pwf->index_mode] + index_ic1_ic2) * psp->ct_size + pwf->index_ct];
  //   double C_l2 = psp->cl[pwf->index_mode][(index_l2_debug * psp->ic_ic_size[pwf->index_mode] + index_ic1_ic2) * psp->ct_size + pwf->index_ct];
  //   double C_l3 = psp->cl[pwf->index_mode][(index_l3 * psp->ic_ic_size[pwf->index_mode] + index_ic1_ic2) * psp->ct_size + pwf->index_ct];
  //   double bispectrum = pbi->bispectra[0][index_l1_debug][index_l2_debug][index_l3-index_l_triangular_min]/sqrt(C_l1*C_l2*C_l3);
  //   
  //   fprintf (stderr, "%g %g %g\n", l3, bispectrum, interpolated_bispectrum);
  // }
    
  /* Interpolation */
  // int L;
  //   
  // for (L=max(2, abs(l1-l2)); L<=min(pbi->l[pbi->l_size-1], l1+l2); L+=2) {  
  //   
  //   double interpolated_bispectrum;
  //   
  //   class_call (bispectra_interpolate (pbi, 0, l1, l2, (double)L, &interpolated_bispectrum), pfi->error_message, pfi->error_message);
  //   
  //   fprintf (stderr, "%d %g\n", L, interpolated_bispectrum);
  // }  
  
  /* Check how much time it takes to interpolate the full bispectrum */
  // int delta = 5;
  // int counter = 0;
  // double interpolated_bispectrum;
  // 
  // int L1, L2, L3;
  // int L_MAX = pbi->l[pbi->l_size-1];
  // 
  // for (L1=2; L1 <= L_MAX; L1+=delta) {
  //   for (L2=2; L2 <= L_MAX; L2+=delta) {
  //     for (L3=max(2, abs(L1-L2)); L3 <= min(L_MAX, L1+L2); ++L3) {
  //       
  //       counter++;
  //       
  //       if ((L1+L2+L3)%2 !=0)
  //         continue;
  //       
  //       class_call (bispectra_interpolate (pbi, 0, L1, L2, L3, &interpolated_bispectrum), pfi->error_message, pfi->error_message)
  // 
  //     }
  //   }
  // }
  // 
  // printf("counter = %d\n", counter);


  return _SUCCESS_;

}



/**
 * Interpolate the bispectrum corresponding to index_bt in the (l1,l2,l3) configuration.
 *
 * Note that the returned value needs to be divided by the window function. We do not do it
 * inside the interpolate function for optimization purposes.
 *
 */
 
int fisher_interpolate_bispectrum (
    struct bispectra * pbi,
    struct fisher * pfi,
    int index_bt,
    int i,
    int j,
    int k,
    double l1,
    double l2,
    double l3,
    double * interpolated_value
    )
{
  
  /* Use the fine mesh when all of the multipoles are small. */
  if ((l1<pfi->l_turnover[0]) && (l2<pfi->l_turnover[0]) && (l3<pfi->l_turnover[0])) {
    class_call (mesh_int (pfi->mesh_workspaces[index_bt][i][j][k][0], l1, l2, l3, interpolated_value),
    pfi->error_message, pfi->error_message);
  }
  /* Use the coarse mesh when any of the multipoles is large. */
  else {
    class_call (mesh_int (pfi->mesh_workspaces[index_bt][i][j][k][1], l1, l2, l3, interpolated_value),
    pfi->error_message, pfi->error_message);
  }
  
  /* Check for nan's */
  if (isnan(*interpolated_value))
    printf ("@@@ WARNING: Interpolated b(%g,%g,%g) = %g for bispectrum %s_%s!!!\n",
    l1, l2, l3, *interpolated_value, pbi->bt_labels[index_bt], pbi->bfff_labels[i][j][k]);
  
  return _SUCCESS_;
  
}



/**
 * Compute the Fisher matrix for all pairs of bispectra.
 *
 * Called by bispectra_init.
 *
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
 
  /* Allocate workspace needed to estimate f_NL */
  struct fisher_workspace * pwf;
  class_alloc (pwf, sizeof(struct fisher_workspace), pfi->error_message);

  // ====================================================================================
  // =                         Compute interpolation weights                            =
  // ====================================================================================

  // ----------------------------------------------------
  // -               Rectangular directions             -
  // ----------------------------------------------------
    
  /* Interpolation weights for the estimator in the rectangular directions (l1 and l2) */
  class_alloc (pwf->delta_l, pbi->l_size*sizeof(double), pfi->error_message);

  /* For the "normal" directions l1 and l2, the weights are the usual ones of the trapezoidal rule
  for a discrete sum. */  
  pwf->delta_l[0] = (pbi->l[1] - pbi->l[0] + 1.)/2.;
  
  for (int index_l=1; index_l < pbi->l_size-1; ++index_l)
    pwf->delta_l[index_l] = (pbi->l[index_l+1] - pbi->l[index_l-1])/2.;
      
  pwf->delta_l[pbi->l_size-1] = (pbi->l[pbi->l_size-1] - pbi->l[pbi->l_size-2] + 1.)/2.;

  // ----------------------------------------------------
  // -                Triangular direction              -
  // ----------------------------------------------------

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

  /* Interpolation weights for the estimator in the triangular direction (l3) */
  class_alloc (pwf->delta_l3, number_of_threads*sizeof(double*), pfi->error_message);

  #pragma omp parallel private (thread)
  {
    
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
    
    class_calloc_parallel (pwf->delta_l3[thread], pbi->l_size, sizeof(double), pfi->error_message);

  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region

  /* Print some info */
  if (pfi->fisher_verbose > 0) {

    char buffer[128];

    if (pfi->bispectra_interpolation == mesh_interpolation)
      strcpy (buffer, "mesh (3D)");
    else if (pfi->bispectra_interpolation == mesh_interpolation_2d)
      strcpy (buffer, "mesh (2D)");
    else if (pfi->bispectra_interpolation == trilinear_interpolation)
      strcpy (buffer, "trilinear");
    else if (pfi->bispectra_interpolation == sum_over_all_multipoles)
      strcpy (buffer, "no");
    else if (pfi->bispectra_interpolation == smart_interpolation)
      strcpy (buffer, "smart");

    printf (" -> computing Fisher matrix with %s interpolation...\n", buffer);
  }



  // ===============================================================
  // =                    Compute Fisher matrix                    =
  // ===============================================================
    

  /* Fill the "pfi->fisher_matrix_of_l1" matrix using mesh interpolation */
  if ((pfi->bispectra_interpolation == mesh_interpolation)
   || (pfi->bispectra_interpolation == mesh_interpolation_2d)) {
  
    class_call (fisher_cross_correlate_mesh(
                  ppr,
                  psp,
                  ple,
                  pbi,
                  pfi,
                  pwf),
      pfi->error_message,
      pfi->error_message);
      
  }
  
  /* Fill the "pfi->fisher_matrix_of_l1" matrix using linear interpolation */
  else if ((pfi->bispectra_interpolation == smart_interpolation)
        || (pfi->bispectra_interpolation == trilinear_interpolation)
        || (pfi->bispectra_interpolation == sum_over_all_multipoles)) {
            
    class_call (fisher_cross_correlate_nodes(
                  ppr,
                  psp,
                  ple,
                  pbi,
                  pfi,
                  pwf),
      pfi->error_message,
      pfi->error_message);
    
  } 
  
  /* So far we have computed the Fisher matrix for each value of l_1; now we sum those values
  to obtain the Fisher matrix for each value of l_max */

  for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
    for (int index_bt_2=index_bt_1; index_bt_2 < pbi->bt_size; ++index_bt_2) {

      double accumulator = 0;
      
      for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {

        /* Perform the sum over l1 */
        accumulator += pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2];
        pfi->fisher_matrix[index_l1][index_bt_1][index_bt_2] = accumulator;
      
        /* The Fisher matrix is symmetric */
        pfi->fisher_matrix[index_l1][index_bt_2][index_bt_1] = pfi->fisher_matrix[index_l1][index_bt_1][index_bt_2];
       
        /* Detectability of the bispectrum, assuming that the bispectrum is a template with f_nl = 1. */
        if (index_bt_1 == index_bt_2)
          pfi->sigma_fnl[index_l1][index_bt_1] = 1./sqrt(pfi->fisher_matrix[index_l1][index_bt_1][index_bt_1]);

        /* Print fnl(l_max) */
        // fprintf (stderr, "%12d %12g\n", pfi->l1[index_l1], pfi->sigma_fnl[index_l1][0]);

      }
    }    
  }
  
  /* Compute inverse Fisher matrix. If there is only one bispectrum type, the function just
  computes 1/fisher. */
  for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1)
    InverseMatrix (pfi->fisher_matrix[index_l1], pbi->bt_size, pfi->inverse_fisher_matrix[index_l1]);
  
  
  /* Print the signal to noise as a function of l_max for all the analised bispectra */
  // fprintf (stderr, "%4s ", "l");
  // for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)    
  //   fprintf (stderr, "%20s ", pbi->bt_labels[index_bt]);
  // fprintf (stderr, "\n");
  // 
  // for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
  //   fprintf (stderr, "%4d ", pfi->l1[index_l1]);
  //   for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
  //     fprintf (stderr, "%20.7g ", sqrt(pfi->fisher_matrix[index_l1][index_bt][index_bt]));
  //   fprintf (stderr, "\n");
  // }

  
  // ===============================================================
  // =                       Print results                         =
  // ===============================================================
  
  /* Print the Fisher matrix */
  if (pfi->fisher_verbose > 0) {

    printf (" -> Fisher matrix:\n");
    
    for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
      
      printf ("\t%20s\t", pbi->bt_labels[index_bt_1]);
      
      printf ("(");
      
      for (int index_bt_2=0; index_bt_2 < pbi->bt_size; ++index_bt_2)
        printf (" %+.5e ", pfi->fisher_matrix[pfi->l1_size-1][index_bt_1][index_bt_2]);

      printf (")\n");
    }
  }

  /* Print the inverse matrix */
  if (pfi->fisher_verbose > 0) {

    printf (" -> Inverse Fisher matrix:\n");
    
    for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
      
      printf ("\t%20s\t", pbi->bt_labels[index_bt_1]);
      
      printf ("(");
      
      for (int index_bt_2=0; index_bt_2 < pbi->bt_size; ++index_bt_2)
        printf (" %+.5e ", pfi->inverse_fisher_matrix[pfi->l1_size-1][index_bt_1][index_bt_2]);

      printf (")\n");
    }
  }

  /* Print the fnl matrix */
  if (pfi->fisher_verbose > 0) {

    printf (" -> fnl matrix (diagonal: 1/sqrt(F_ii), upper: F_12/F_11, lower: F_12/F_22):\n");
    
    for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
      
      printf ("\t%20s\t", pbi->bt_labels[index_bt_1]);
      
      printf ("(");
      
      for (int index_bt_2=0; index_bt_2 < pbi->bt_size; ++index_bt_2) {

        /* Diagonal elements = 1/sqrt(F_ii) */ 
        if (index_bt_1==index_bt_2)
          printf (" %+.5e ", 1./sqrt(pfi->fisher_matrix[pfi->l1_size-1][index_bt_1][index_bt_1]));
        /* Upper triangle = F_12/F_11, lower triangle = F_12/F_22. */
        else
          printf (" %+.5e ", pfi->fisher_matrix[pfi->l1_size-1][index_bt_1][index_bt_2]/pfi->fisher_matrix[pfi->l1_size-1][index_bt_1][index_bt_1]);

      }

      printf (")\n");
    }
  }
 
  return _SUCCESS_;
  
}






/**
 *
 * This function fills pfi->fisher_matrix_of_l1.
 *
 */
int fisher_cross_correlate_nodes (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi,
       struct fisher_workspace * pwf
       )
{  

  /* Physical limits in l */
  int l_min = pfi->l_min;
  int l_max = pfi->l_max;
  
  /* Edit here in order to consider only a subset of the sum */
  int l1_min_global=l_min, l1_max_global=l_max,
      l2_min_global=l_min, l2_max_global=l_max,
      l3_min_global=l_min, l3_max_global=l_max;

  /* We shall count the (l1,l2,l3) configurations over which we compute the estimator */
  long int counter = 0;

  // =======================================================================================
  // =                            Compute Fisher matrix elements                           =
  // =======================================================================================

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
    
      if ((l1 < l1_min_global) || (l1 > l1_max_global))
        continue;
    
      if (pfi->fisher_verbose > 2)
        printf ("     * computing Fisher matrix for l1=%d\n", l1);

      // ------------------------------------------------
      // -                  Sum over l2                 -
      // ------------------------------------------------
    
      for (int index_l2=0; index_l2<=index_l1; ++index_l2) {
        
        int l2 = pfi->l2[index_l2];
  
        if ((l2 < l2_min_global) || (l2 > l2_max_global))
          continue;
        
        /* Skip those configurations that are forbidden by the triangular condition (optimization) */
        if (l2 < l1/2)
          continue;
        
        /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = min (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
        
        if (pfi->bispectra_interpolation == smart_interpolation) {

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
                                
          if ((l3 < l3_min_global) || (l3 > l3_max_global))
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
          
          /* Weights for the interpolation. We use this variable to also store some other coefficients,
          for optimisation purposes. */
          double interpolation_weight;
          
          if (pfi->bispectra_interpolation == smart_interpolation) {
            interpolation_weight = pwf->delta_l[index_l1] * pwf->delta_l[index_l2] * pwf->delta_l3[thread][index_l3];
          }
          else if (pfi->bispectra_interpolation == trilinear_interpolation) {
            /* Weight for trilinear interpolation. The factor 1/2 comes from the fact that when interpolating over even
            configurations only, we are implicitly assuming non-zero values for the odd points, which instead are
            forced to vanish by the 3J. (Note that this factor, in the case of the smart interpolation, is already
            included in pwf->delta_l3) */
            interpolation_weight = 0.5 * pwf->delta_l[index_l1] * pwf->delta_l[index_l2] * pwf->delta_l[index_l3];
          }
          /* All other kind of interpolation are sums over all l's and require no interpolation weight */
          else {            
            interpolation_weight = 1;
          }
          
          /* If we are taking all the l-points, then any interpolation scheme reduces to a simple sum over the
          multipoles, and then we can naturally implement the delta factor in KSW2005. */
          if (ppr->l_linstep == 1) {
            if ((l1==l2) && (l1==l3)) 
              interpolation_weight /= 6.;
            else if ((l1==l2) || (l1==l3) || (l2==l3))
              interpolation_weight /= 2.;
          }
          
          /* Include the double 3j symbol */
          interpolation_weight *= pfi->I_l1_l2_l3[index_l1_l2_l3] * pfi->I_l1_l2_l3[index_l1_l2_l3];
          
          /* Compute the Fisher matrix. In the simple case of only one probe (eg. temperature -> TTT),
          each element of the Fisher matrix is simply given by the square of the bispectrum divided
          by three C_l's. When dealing with two probes (eg. T and E -> TTT,TTE,TET,ETT,TEE,ETE,EET,EEE),
          each element of the Fisher matrix, contains a quadratic sum over the 8 available bispectra,
          for a total of 64 terms. */
          double inverse_covariance = 1;
          for (int i = 0; i < pbi->bf_size; ++i) { for (int p = 0; p < pbi->bf_size; ++p) {
          for (int j = 0; j < pbi->bf_size; ++j) { for (int q = 0; q < pbi->bf_size; ++q) {
          for (int k = 0; k < pbi->bf_size; ++k) { for (int r = 0; r < pbi->bf_size; ++r) {

            double inverse_covariance = C1[i][p] * C2[j][q] * C3[k][r];

            for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
              for (int index_bt_2=index_bt_1; index_bt_2 < pbi->bt_size; ++index_bt_2) {              

                  pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2] += 
                    inverse_covariance *
                    pbi->bispectra[index_bt_1][i][j][k][index_l1_l2_l3] *
                    pbi->bispectra[index_bt_2][p][q][r][index_l1_l2_l3] *
                    interpolation_weight;

                }}} // ijk
              }}} // pqr
            } // bt_2
          } // bt_1
                
          #pragma omp atomic
          counter++;          

        } // end of for(index_l3)

      } // end of for(index_l2)

      /* Now that we have completed the computation for a given value of l1, we can perform some tasks on the result. */
      for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
        for (int index_bt_2=index_bt_1; index_bt_2 < pbi->bt_size; ++index_bt_2) {
      
          /* Include sky coverage */
          pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2] *= pfi->f_sky;
                              
          /* The Fisher matrix is symmetric */
          pfi->fisher_matrix_of_l1[index_l1][index_bt_2][index_bt_1] = pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2];
        }
      }
      
    } // end of for(index_l1)


    // ------------------------------------------------------------------------
    // -                            End of all loops                          -
    // ------------------------------------------------------------------------
    
  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region
  

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
  int index_l3_max = min (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      
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
  for (int l3 = max(pfi->l_min, abs(l1-l2)); l3 < pbi->l[index_l3_min]; ++l3)
    if((l1+l2+l3)%2==0)
      delta_l3[index_l3_min] += 1.;
  
  for (int l3 = pbi->l[index_l3_max]+1; l3 <= min(l2,l1+l2); ++l3)
    if((l1+l2+l3)%2==0)
      delta_l3[index_l3_max] += 1;
  
  /* Check that, for a function that always evaluates to one, the weights sum up to the number
  of even/odd points in the integration range. For a given pair of (l1,l2), the range goes
  from max(2,|l1-l2|) to min(l2,l1+l2). 
  Remember that now 'index_l3_min' and 'index_l3_max' have been modified with respect to those in
  pbi->index_l_triangular_min and pbi->index_l_triangular_max to be even/odd. */
  double sum = 0;
  for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
    sum += delta_l3[index_l3];
  
  int n_even=0, n_odd=0;
  for (int l3 = max(pfi->l_min,abs(l1-l2)); l3 <= min(l2,l1+l2); ++l3) {
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
 *
 * This function fills pfi->fisher_matrix_of_l1.
 *
 */
int fisher_cross_correlate_mesh (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi,
       struct fisher_workspace * pwf
       )
{  

  /* Physical limits in l */
  int l_min = pfi->l_min;
  int l_max = pfi->l_max;
  
  /* Edit here in order to consider only a subset of the sum */
  int l1_min_global=l_min, l1_max_global=l_max,
      l2_min_global=l_min, l2_max_global=l_max,
      l3_min_global=l_min, l3_max_global=l_max;

  /* We shall count the (l1,l2,l3) configurations over which we compute the estimator */
  long int counter = 0;

  // =======================================================================================
  // =                            Compute Fisher matrix elements                           =
  // =======================================================================================

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

  /* Temporary values needed for the computation of the 3j symbol */
  double ** temporary_three_j;
  class_alloc (temporary_three_j, number_of_threads*sizeof(double *), pfi->error_message);

  /* Temporary values needed to store the bispectra interpolated in a given (l1,l2,l3) configuration 
  and for a certain type of bispectrum (e.g. local_tte). The array is indexed as
  interpolated_bispectra[thread][index_bt][i][j][k] */
  double ***** interpolated_bispectra;
  class_alloc (interpolated_bispectra, number_of_threads*sizeof(double *****), pfi->error_message);  

  #pragma omp parallel shared (abort,counter) private (thread)
  {
    
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    /* The temporary array will contain at most 2*l_max+1 values when l1=l2=l_max, with the +1
    accounting for the l=0 case. */
    int l3_min_3j, l3_max_3j, temporary_three_j_size;
    class_alloc_parallel (temporary_three_j[thread], (2*l_max+1)*sizeof(double ****), pfi->error_message);

    /* As many interpolations as the total number of bispectra (pbi->bt_size * pbi->bf_size^3) */
    class_alloc_parallel (interpolated_bispectra[thread], pbi->bt_size*sizeof(double ***), pfi->error_message);
    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
      class_alloc_parallel (interpolated_bispectra[thread][index_bt], pbi->bf_size*sizeof(double **), pfi->error_message);
      for (int i = 0; i < pbi->bf_size; ++i) {
        class_alloc_parallel (interpolated_bispectra[thread][index_bt][i], pbi->bf_size*sizeof(double *), pfi->error_message);
        for (int j = 0; j < pbi->bf_size; ++j)
          class_alloc_parallel (interpolated_bispectra[thread][index_bt][i][j], pbi->bf_size*sizeof(double), pfi->error_message);
      }
    }


    // ------------------------------------------------
    // -                  Sum over l1                 -
    // ------------------------------------------------

    #pragma omp for schedule (dynamic)
    for (int index_l1=0; index_l1<pfi->l1_size; ++index_l1) {
    
      int l1 = pfi->l1[index_l1];
    
      if ((l1 < l1_min_global) || (l1 > l1_max_global))
        continue;
    
      if (pfi->fisher_verbose > 2)
        printf ("     * computing Fisher matrix for l1=%d\n", l1);

      /* The triangular condition imposes l3 >= l1-l2, while the bound on the summation imposes l3 <= l2. Hence,
      the minimum allowed value for l2 is l1/2. */
      int l2_min = max(l1/2,2);

      // ------------------------------------------------
      // -                  Sum over l2                 -
      // ------------------------------------------------
    
      for (int l2=l2_min; l2<=l1; ++l2) {
  
        if ((l2 < l2_min_global) || (l2 > l2_max_global))
          continue;
      
        /* Compute the actual 3j symbol for all allowed values of l3 */
        class_call_parallel (threej_l1 (
                      l1, l2, 0, 0,
                      &l3_min_3j, &l3_max_3j,
                      &temporary_three_j[thread],
                      &temporary_three_j_size,
                      pfi->error_message       
                      ),
          pfi->error_message,
          pfi->error_message);
  
        int l3_min = max(l1-l2,2);

        // ------------------------------------------------
        // -                  Sum over l3                 -
        // ------------------------------------------------

        for (int l3=l3_min; l3<=l2; ++l3) {
                                
          if ((l3 < l3_min_global) || (l3 > l3_max_global))
            continue;
        
          /* Parity invariance enforces that the contribution of the modes with odd l1+l2+l3 is zero. Even if you
          comment this out, the result won't change as the I_l1_l2_l3 coefficient would vanish in that case. */
          if ((l1+l2+l3)%2!=0) {
            #pragma omp atomic
            counter++;
            continue;
          }

          if ((l3 < l3_min_global) || (l3 > l3_max_global))
            continue;
  
          // ---------------------------------------------------------------------------
          // -                            Interpolate bispectra                        -
          // ---------------------------------------------------------------------------

          /* Compensate the effect of the window function. Note that the C_l's appearing here should not
          be corrected for instrumental noise. */  

          /* TODO: choose Cl's according to the considered probe (ie not only TTT) */
          double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
          double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
          double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
          double inverse_window = C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3;

          /* Factor that relates the reduced bispectrum to the angle-averaged one */
          double I_l1_l2_l3 = sqrt((2.*l1+1.)*(2.*l2+1.)*(2.*l3+1.)/(4.*_PI_)) * temporary_three_j[thread][l3-l3_min_3j];

          /* Interpolate all bispectra in (l1,l2,l3). These two loops go over all the possible bispectra,
          eg. local_ttt, equilateral_ete, intrinsic_ttt. */
          for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
            for (int i = 0; i < pbi->bf_size; ++i) {
            for (int j = 0; j < pbi->bf_size; ++j) {
            for (int k = 0; k < pbi->bf_size; ++k) {

              class_call_parallel (fisher_interpolate_bispectrum(
                                     pbi,
                                     pfi,
                                     index_bt,
                                     i,j,k,
                                     l1,
                                     l2,
                                     l3,
                                     &interpolated_bispectra[thread][index_bt][i][j][k]),
                pfi->error_message,
                pfi->error_message);

              interpolated_bispectra[thread][index_bt][i][j][k] *= I_l1_l2_l3 * inverse_window;

            }}} // end of for(ijk)
          } // end of for(index_bt)
          

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
          
          /* Correction factor for the variance */
          double one_over_delta = 1;
          if ((l1==l2) && (l1==l3))
            one_over_delta = 1/6.;
          else if ((l1==l2) || (l1==l3) || (l2==l3))
            one_over_delta = 0.5;
          
          /* Shortcuts to the inverse of the cross-power spectrum */
          double ** C1 = pfi->inverse_C[l1-2];
          double ** C2 = pfi->inverse_C[l2-2];
          double ** C3 = pfi->inverse_C[l3-2];
          
          /* Compute the Fisher matrix. In the simple case of only one probe (eg. temperature -> TTT),
          each element of the Fisher matrix is simply given by the square of the bispectrum divided
          by three C_l's. When dealing with two probes (eg. T and E -> TTT,TTE,TET,ETT,TEE,ETE,EET,EEE),
          each element of the Fisher matrix, contains a quadratic sum over the 8 available bispectra,
          for a total of 64 terms. */
          double inverse_covariance = 1;
          for (int i = 0; i < pbi->bf_size; ++i) { for (int p = 0; p < pbi->bf_size; ++p) {
          for (int j = 0; j < pbi->bf_size; ++j) { for (int q = 0; q < pbi->bf_size; ++q) {
          for (int k = 0; k < pbi->bf_size; ++k) { for (int r = 0; r < pbi->bf_size; ++r) {

            double inverse_covariance = C1[i][p] * C2[j][q] * C3[k][r];

            for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
              for (int index_bt_2=index_bt_1; index_bt_2 < pbi->bt_size; ++index_bt_2) {              

                  pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2] += 
                    one_over_delta *
                    inverse_covariance *
                    interpolated_bispectra[thread][index_bt_1][i][j][k] *
                    interpolated_bispectra[thread][index_bt_2][p][q][r];

                }}} // ijk
              }}} // pqr
            } // bt_2
          } // bt_1
                
          #pragma omp atomic
          counter++;          

        } // end of for(index_l3)
      } // end of for(index_l2)

      /* Now that we have completed the computation for a given value of l1, we can perform some tasks on the result. */
      for (int index_bt_1=0; index_bt_1 < pbi->bt_size; ++index_bt_1) {
        for (int index_bt_2=index_bt_1; index_bt_2 < pbi->bt_size; ++index_bt_2) {
      
          /* Include sky coverage */
          pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2] *= pfi->f_sky;
                    
          /* Include the interpolation weight for the l1 direction, if needed */
          if (pfi->bispectra_interpolation == mesh_interpolation_2d)
            pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2] *= pwf->delta_l[index_l1];
          
          /* The Fisher matrix is symmetric */
          pfi->fisher_matrix_of_l1[index_l1][index_bt_2][index_bt_1] = pfi->fisher_matrix_of_l1[index_l1][index_bt_1][index_bt_2];   
        }
      }
      
    } // end of for(index_l1)


    // ------------------------------------------------------------------------
    // -                            End of all loops                          -
    // ------------------------------------------------------------------------
    
    free (temporary_three_j[thread]);
    
    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
      for (int i = 0; i < pbi->bf_size; ++i) {
        for (int j = 0; j < pbi->bf_size; ++j) {
          free (interpolated_bispectra[thread][index_bt][i][j]); }
        free (interpolated_bispectra[thread][index_bt][i]); }
      free (interpolated_bispectra[thread][index_bt]); }
    free (interpolated_bispectra[thread]);
    
  } if (abort == _TRUE_) return _FAILURE_; // end of parallel region
  
  free (temporary_three_j);
  free (interpolated_bispectra);

  /* LOCAL */
  // double expected = 691351; // l_max = 200, constant function
  // double expected = 26.2101; // l_max = 200, vanishing function
  // double expected = -5.7540e-18; // l_max = 100, bispectrum only with 50r14to15
  // double expected = -0.9982; // l_max = 100, bispectrum/sqrt(cl1cl2cl3) only with 50r14to15
  // double expected = -8.38648e-13; // l_max = 100, bispectrum*l1^4/3*l2^4/3*l3^4/3 only with 50r14to15
  // double expected = -2.46446e-12; // l_max = 200, bispectrum*l1^4/3*l2^4/3*l3^4/3 only with 50r14to15
  // double expected = 1.57974e6; // l_max = 200, bispectrum*l1^4*l2^4*l3^4 only with 50r14to15
  // double expected = -5.03355e-18; // l_max = 200, bispectrum only with 50r14to15
  // double expected = 0.058806; // l_max = 2000, F^ii with 1000r13to16 with isw
  // double expected = 0.000479757; // l_max = 200, F^ii with 1000r13to16 with iws
  // double expected = 0.000513081; // l_max = 200, F^ii with 1000r13to16, NO DELTA_l1_l2_l3 with iws
  
  /* EQUILATERAL */
  // double expected = 4.26226e-18; // l_max = 200, bispectrum only with 50r14to15
  // double expected = 1.42016e-05; // l_max = 200, equilateral F^ii with 1000r13to16 with isw
  
  /* ORTHOGONAL */
  // double expected = 1.52214e-17; // l_max = 200, bispectrum only with 50r14to15
  // double expected = 3.54087e-05; // l_max = 200, orthogonal F^ii with 1000r13to16 with isw
  // double expected = 3.53817e-05; // l_max = 200, bispectrum only with 50r14to15

  /* RESULT */
  // double found = F_22[l1_iteration_max];
  // printf ("#####   found = %g, expected = %g, diff = %g   #####\n", found, expected, 100*(found/expected - 1.)); 

    
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

  
  // ==========================================================================
  // =                           Allocate memory                              =
  // ==========================================================================

  /* Allocate memory for the cross-power spectrum and its inverse */
  class_alloc (pfi->C, pfi->full_l_size*sizeof(double **), pfi->error_message);
  class_alloc (pfi->inverse_C, pfi->full_l_size*sizeof(double **), pfi->error_message);

  for (int l=pfi->l_min; l <= pfi->l_max; ++l) {

    class_alloc (pfi->C[l-2], pbi->bf_size*sizeof(double*), pfi->error_message);
    class_alloc (pfi->inverse_C[l-2], pbi->bf_size*sizeof(double*), pfi->error_message);

    for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {
      class_alloc (pfi->C[l-2][index_bf], pbi->bf_size*sizeof(double), pfi->error_message);
      class_alloc (pfi->inverse_C[l-2][index_bf], pbi->bf_size*sizeof(double), pfi->error_message);
    }
  }


  // ==========================================================================
  // =                      Compute cross-power spectrum                      =
  // ==========================================================================

  /* The cross-power spectrum is used to build the covariance matrix and, therefore, needs to be
  include the instrumental noise. We assume that this affects only the direct C_l's,
  eg. C_l^TT and C_l^EE have noise but C_l^TE does not. */
  
  for (int l=pfi->l_min; l <= pfi->l_max; ++l) {  
    for (int i = 0; i < pbi->bf_size; ++i) {
      for (int j = 0; j < pbi->bf_size; ++j) {

        /* Compute ij element */
        pfi->C[l-2][i][j] = pbi->cls[pbi->index_ct_of_bf[i][j]][l-2];

        /* Add noise only for diagonal combinations (TT, EE) */
        if (i==j) pfi->C[l-2][i][j] += pfi->N_l[i][l-2];

        /* If a C_l vanishes, then we risk to have infinities in the inverse */
        class_test (pfi->C[l-2][i][j] == 0., pfi->error_message,
          "a C_l (l=%d) was found to be zero. Stopping to prevent seg fault. Maybe your k-max is too small so that C_l=0 for large l's?",l);
      }
    }
    /* Compute the inverse of C_l^ij for the considered l */
    InverseMatrix (pfi->C[l-2], pbi->bf_size, pfi->inverse_C[l-2]);
    
    /* Some debug */
    // double factor = 1e12*pba->T_cmb*pba->T_cmb*l*(l+1)/(2*_PI_);
    // fprintf (stderr, "%d %g %g %g %g\n", l,
    // factor*pfi->C[l-2][0][0], factor*pfi->C[l-2][1][1], factor*pfi->C[l-2][0][1], factor*pfi->C[l-2][1][0]);

  } // end of for(l)
  
  
  return _SUCCESS_;
  
}








