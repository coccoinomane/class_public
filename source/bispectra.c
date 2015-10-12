/** @file bispectra.c
 * 
 * Compute bispectra in harmonic (l1,l2,l3) space.
 *
 * The bispectrum is practically the three-point function of a random field. In cosmology,
 * it is important because it encodes the non-linearities in our Universe. For example,
 * the bispectrum of cold dark matter contains informations on the non-linear processes
 * of gravity, while the presence of a primordial CMB bispectrum might be indicative of
 * non-Gaussian initial conditions.
 *
 * In this module we compute the CMB bispectrum in temperature and polarisation. Details
 * on the theory and computation can be found in chapter 6 of my thesis
 * (http://arxiv.org/abs/1405.2280).
 *
 * We deal with three types of bispectra:
 *
 * -# The analytical bispectra, which are built using the C_l power spectra computed
 *    in the spectra.c module. These are the simplest bispectra to compute as they do
 *    not require to solve numerical integrals. Examples of such bispectra are the
 *    CMB-lensing bispectrum (http://uk.arxiv.org/abs/1101.2234) and the squeezed limit
 *    of the intrinsic bispectrum (http://arxiv.org/abs/1204.5018 and
 *    http://arxiv.org/abs/1109.1822).
 *
 * -# The separable bispectra, which can be computed by solving three 1D integrals
 *    involving the transfer functions computed in the transfer.c module,
 *    as first shown in Komatsu et al. 2001 (http://arxiv.org/abs/astro-ph/0005036).
 *    Examples are the local, equilateral and orthogonal templates of non-Gaussianity
 *    (see Planck paper, http://arxiv.org/abs/1303.5084).
 *
 * -# The non-separable bispectra, which require solving a non-separable 3D integral in
 *    (k1,k2,k3), involving the transfer functions and a primordial shape function in
 *    (k1,k2,k3). The user can specify the arbitrary shape shape function in k1,k2,k3.
 *    Examples are the Galileon models in http://arxiv.org/abs/0905.3746.
 *
 * The intrinsic bispectrum of the CMB is also a non-separable bispectrum, but we leave it
 * for the bispectra2.c module, as it requires the computation of second-order perturbations
 * (see chapter 6 of http://arxiv.org/abs/1405.2280).
 *
 * The main functions in this module that can be called externally are:
 * -# bispectra_init() to run the module, requires the background, thermodynamics,
 *    perturbations, bessel, transfer, primordial, spectra and lensing modules.
 * -# bispectra_free() to free all the memory associated to the module.
 * 
 * If the user specified 'store_bispectra_to_disk=yes', the module will save all the 
 * bispectra but the intrinsic one in pbi->bispectra to disk after computation.
 *
 *
 * REDUCED BISPECTRA
 *
 * This module computes the REDUCED bispectrum b_l1_l2_l3, that is, the full bispectrum
 * <a_l1m1 a_l2m2 a_l3m3> divided by the Gaunt symbol. The full bispectrum vanishes
 * for even values of l1+l2+l3, due to the presence of the 3j symbol (l1,l2,l3)(0,0,0)
 * in the Gaunt symbol. The reduced bispectrum, on the other hand, can be defined to
 * be continuous for all l.
 *
 * For certain bispectra, the 3j-symbol (l1,l2,l3)(0,0,0) cannot be pulled out
 * analytically because it does not appear explicitly. This is the case for
 * the CMB-lensing bispectrum in presence of polarisation and for the non-scalar
 * (m != 0) intrinsic bispectrum. In these cases, the reduced bispectrum can be
 * obtained by explicitly dividing the full bispectrum by 3J(l1,l2,l3)(0,0,0).
 * The result will be then undefined for odd values of l1+l2+l3, where the 3J
 * vanishes. However, this is not much of an issue because these configurations
 * are not physical. For example, in the Fisher matrix only the even l1+l2+l3
 * are considered.
 * 
 * The explicit approach is problematic for another reason. It highlights one of
 * SONG's own limitations: the three multipoles (l1,l2,l3) are drawn from the same
 * 1D grid. Including an odd value of l1 means that we will have plenty of
 * combinations of (l1,l2,l3) with odd l1+l2+l3, where the 3J symbol
 * in the denominator is not defined.
 *
 * To circumvent this issue, an option is to have an l-grid where all the l are
 * even. This is not completely satisfactory because half of the configurations
 * (those with even l1+l2+l3 but two odd components, like 2,3,3 or 2,3,7) will be
 * skipped, even if they are perfectly valid. These gaps in the (l1,l2,l3)
 * sampling make the interpolation of the bispectrum harder.
 *
 * A better option is to use the recursive relation in Schulten & Gordon, 1961
 * to express the ratio of 3J as a sum of analytic functions. In this way,
 * we are able to pull out the 3J symbol from the full bispectrum even for
 * configurations with odd l1+l2+l3. This allows us to use an arbitrary 
 * (l1,l2,l3) grid and to avoid gaps in the bispectrum, thus allowing for
 * easier interpolation.
 * 
 * As of now, SONG implements the recursive relation for the analytical
 * bispectra such as the CMB-lensing bispectrum and the quadratic bispectrum.
 * It is not used yet to compute the intrinsic bispectrum, for which the
 * relation needs to be adjusted. This means that if you want the non-scalar
 * contributions to the intrinsic bispectrum, SONG will automatically switch
 * to an all-even l-grid.
 *
 * NOTE: for odd parity bispectra such as those involving odd combinations
 * of B-modes, there is no obvious definition for the reduced bispectrum. 
 * Since SONG does not compute odd parity bispectra, it is not an issue (yet).
 * See last of paragraph of Sec. 4 of Lewis, Challinor & Hanson 2011
 * (http://arxiv.org/abs/1101.2234) for more detail.
 * 
 * TODO: Implement an independent l3 array, so that we can tune it to have
 * only the even l1+l2+l3 configurations without sacrificing the (2,3,3)
 * configurations.
 *
 * Created by Guido W. Pettinari on 19.07.2012.
 * Last modified by Guido W. Pettinari on 06.10.2015
 */

#include "bispectra.h"


/**
 * Fill the fields in the bispectra structure, especially the pbi->bispectra
 * array.
 * 
 * This function calls the other bispectra_XXX functions in the order needed to 
 * compute the requested bispectra. It requires the following modules to be
 * filled and available: background, thermodynamics, perturbations, bessel,
 * transfer, primordial, spectra, lensing.
 *
 * Details on the physics and on the adopted method can be found in my thesis
 * (http://arxiv.org/abs/1405.2280), especially chapters 6. The code itself is
 * extensively documented and hopefully will give you further insight; in doubt,
 * feel free to drop an email to guido.pettinari@gmail.com.
 *
 * In detail, this function does:
 *
 * -# Determine which bispectra need to be computed and their l-sampling via
 *    bispectra_indices(), according to the content of the previous modules.
 *
 * -# Compute the three types of bispectra (analytic, separable, non-separable)
 *    via bispectra_harmonic().
 *
 * -# If requested, store to disk the content of pbi->bispectra via
 *    bispectra_store_to_disk().
 *
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
  if (ppt->has_cmb_bispectra == _FALSE_) {

    pbi->has_bispectra = _FALSE_;

    printf_log_if (pbi->bispectra_verbose, 0,
      "No bispectra requested. Bispectra module skipped.\n");

    pbi->bt_size=0;

    return _SUCCESS_;
  }
  else {
    
    pbi->has_bispectra = _TRUE_;
    
    printf_log_if (pbi->bispectra_verbose, 0,
      "Computing bispectra.\n");
  }



  // =====================================================================================
  // =                                    Preparations                                   =
  // =====================================================================================

  /* Initialize indices & arrays in the bispectra structure */

  class_call (bispectra_indices (ppr,pba,pth,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);



  // =====================================================================================
  // =                                 Compute bispectra                                 =
  // =====================================================================================
  
  /* Compute the three types of CMB bispectra: analytic, separable and non-separable */
  
  class_call (bispectra_harmonic (ppr,pba,pth,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);



  // =====================================================================================
  // =                                  Compute lensing                                  =
  // =====================================================================================

  /* Compute the three types of CMB bispectra: analytic, separable and non-separable */
  
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    if (pbi->lens_me[index_bt] == _TRUE_) {

      class_call (bispectra_lensing (ppr,pba,pth,ppt,pbs,ptr,ppm,psp,ple,pbi,index_bt),
        pbi->error_message,
        pbi->error_message);
        
    }
  }


  // ====================================================================================
  // =                              Prepare interpolation                               =
  // ====================================================================================



  /* Apart from pbi->bispectra, at this point all the arrays in the module have been
  filled.  If the user requested to load the bispectra from disk, we can stop the
  execution of this module now without regrets. */
    
  if (ppr->load_bispectra_from_disk == _TRUE_) {

    printf_log_if (pbi->bispectra_verbose, 0, 
      " -> the intrinsic and non-separable bispectra will be read from disk\n");

    /* Uncomment to produce the bispectra output files again */
    // class_call (bispectra_output (ppr,pba,pth,ppt,pbs,ptr,ppm,psp,ple,pbi),
    //   pbi->error_message,
    //   pbi->error_message);

    return _SUCCESS_;

  }


  // =====================================================================================
  // =                                  Produce output                                   =
  // =====================================================================================

  /* Create output files containing the bispectra */

  class_call (bispectra_output (ppr,pba,pth,ppt,pbs,ptr,ppm,psp,ple,pbi),
    pbi->error_message,
    pbi->error_message);



  // =====================================================================================
  // =                              Store bispectra to disk                              =
  // =====================================================================================
  
  if (ppr->store_bispectra_to_disk==_TRUE_) {

    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

      /* Save the non-separable bispectra to disk. We do not save the other ones because
      they take no time to recompute. */
      if (pbi->bispectrum_type[index_bt] == non_separable_bispectrum)
        class_call (bispectra_store_to_disk (
                      pbi,
                      index_bt),
          pbi->error_message,
          pbi->error_message);
    }  
  }
  
  /* Check that we correctly filled the bispectra array (but only if there are no
  intrinsic bispectra left to be computed)*/
  if (pbi->n[intrinsic_bispectrum] < 1)
    class_test_permissive (pbi->count_allocated_for_bispectra != pbi->count_memorised_for_bispectra,
      pbi->error_message,
      "there is a mismatch between allocated (%ld) and used (%ld) space!",
      pbi->count_allocated_for_bispectra, pbi->count_memorised_for_bispectra);


	/* Debug - Print some bispectra configurations */
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
  //         if (strstr(pbi->bt_labels[index_bt], "cmb_lensing")!=NULL) {
  //           /* Triangular configuration */
  //           // if ((l1==2000) && (l2==1600)) {
  //           //   fprintf (stderr, "%12d %12g %12g %12g %12g %12g %12g %12g %12g\n",
  //           //     l3,
  //           //     pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][1][1][1][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][0][0][1][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][0][1][0][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][1][0][0][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][0][1][1][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][1][0][1][index_l1_l2_l3],
  //           //     pbi->bispectra[index_bt][1][1][0][index_l1_l2_l3]
  //           //   );
  //           // }
  //           /* Squeezed configuration */
  //           // if ((l3==6) && (l1==l2)) {
  //           //   fprintf (stderr, "%12d %12g %12g %12g %12g %12g %12g %12g %12g\n",
  //           //     l1,
  //           //     pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][1][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][0][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][0][1][0][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][1][0][0][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][0][1][1][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][1][0][1][index_l1_l2_l3] * equilateral_normalisation,
  //           //     pbi->bispectra[index_bt][1][1][0][index_l1_l2_l3] * equilateral_normalisation
  //           //   );
  //           // }
  //           /* Equilateral configuration */
  //           if ((l1==l2) && (l2==l3)) {
  //             fprintf (stderr, "%12d %12g %12g %12g\n",
  //               l1,
  //               pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3],
  //               equilateral_normalisation,
  //               pbi->bispectra[index_bt][0][0][0][index_l1_l2_l3]*equilateral_normalisation
  //             );
  //           }
  //         }
  //       } // end of for(index_l3)
  //     } // end of for(index_l2)
  //   } // end of for(index_l1)
  // } // end of for(index_bt)

  return _SUCCESS_;

}



/**
 * Evaluate the reduced bispectrum at a given (l1,l2,l3) configuration inside pbi->l.
 * 
 * This function returns b(l1,l2,l3) for any configuration, as long all l belong to
 * pbi->l and satisfy the triangular condition. If the bispectrum has been lensed,
 * it also returns the unlensed bispectrum.
 *
 * We need this function rather than accessing directly pbi->bispectra because we
 * computed the bispectrum only for l1>=l2>=l3. We obtain the values outside that
 * range using the symmetry of the bispectrum with respect to permutations of (1,2,3).
 *
 * Note that this approach won't work for any squeezed-limit approximations, such
 * as the intrinsic_squeezed  and local_squeezed bispectra, because they are asymmetric
 * by construction. For these bispectra, we set the output to zero outside the l1>=l2>=l3
 * regime.
 */

int bispectra_at_l1l2l3 (
    struct bispectra * pbi,
    int index_bt,
    int index_l1, int index_l2, int index_l3,
    int X, int Y, int Z,
    double * bispectrum,
    double * bispectrum_unlensed
    )
{

  int l1 = pbi->l[index_l1];
  int l2 = pbi->l[index_l2];
  int l3 = pbi->l[index_l3];

  /* Debug - print arguments */
  // printf ("(%d^%s,%d^%s,%d^%s)\n", l1, pbi->bf_labels[X],
  //   l2, pbi->bf_labels[Y], l3, pbi->bf_labels[Z]);
  
#ifdef DEBUG
  /* Test that l1, l2 and l3 satisfy the triangular condition */
  class_test (!is_triangular_int(l1,l2,l3),
    pbi->error_message,
    "(l1=%d, l2=%d, l3=%d) is not a triangular configuration", l1, l2, l3);
#endif // DEBUG
  
  /* Find the ordering of (l1,l2,l3) */
  
  int index_l[4] = {0, index_l1, index_l2, index_l3};

  int order[4];

  class_call (ordering_int (index_l, order, pbi->error_message),
    pbi->error_message, pbi->error_message);


  /* Index of the current (l1,l2,l3) configuration */

  int index_1 = index_l[order[1]]; /* Smallest l */
  int index_2 = index_l[order[2]]; /* Mid l */
  int index_3 = index_l[order[3]]; /* Largest l */

  int index_1_max = MIN (index_2, pbi->index_l_triangular_max[index_3][index_2]);
  long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_3][index_3-index_2][index_1_max-index_1];


  /* Extract the bispectrum in (l1,l2,l3), using the fact that a permutation of
  (l1,l2,l3) is cancelled by the same permutation of (X,Y,Z). For example:
  b^TTE(l1,l2,l3) = b^TET(l1,l3,l2). */

  int XYZ[4] = {0, X, Y, Z};

  *bispectrum = pbi->bispectra[index_bt][XYZ[order[3]]][XYZ[order[2]]][XYZ[order[1]]][index_l1_l2_l3];
  
  if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt] == _TRUE_))
    *bispectrum_unlensed = *bispectrum - pbi->lensing_correction[index_bt][XYZ[order[3]]][XYZ[order[2]]][XYZ[order[1]]][index_l1_l2_l3];
    

  /* The squeezed approximation in SONG are computed assuming that l3>=l2>=l1 */ 

  if (((((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed)))
     ||((pbi->has_local_squeezed == _TRUE_) && (index_bt == pbi->index_bt_local_squeezed))
     ||((pbi->has_cmb_lensing_squeezed == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing_squeezed)))
     && ((l2>l3) || (l1>l3) || (l1>l2))) {
      *bispectrum = 0;
      *bispectrum_unlensed = 0;
  }
  
  return _SUCCESS_;
  
}




/**
 * Save the bispectra in pbi->bispectra to disk for a given type.
 * 
 * The bispectra will be saved to the file in pbi->bispectra_files[index_bt].
 */
int bispectra_store_to_disk (
    struct bispectra * pbi,
    int index_bt
    )
{

  /* Open file for writing */
  class_open (pbi->bispectra_files[index_bt], pbi->bispectra_paths[index_bt], "a+b", pbi->error_message);

  /* Print some debug */
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * writing bispectra to disk for index_bt=%d on '%s'\n",
    index_bt, pbi->bispectra_paths[index_bt]);

  /* Write all the independent (l1,l2,l3) triplets for this bispectrum */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        fwrite(
              pbi->bispectra[index_bt][X][Y][Z],
              sizeof(double),
              pbi->n_independent_configurations,
              pbi->bispectra_files[index_bt]
              );

  /* Close file */
  fclose(pbi->bispectra_files[index_bt]);
  
  return _SUCCESS_;
  
}


/**
 * Free all the memory space allocated by bispectra_init() in the
 * bispectrum structure.
 */
int bispectra_free(
     struct precision * ppr,
     struct perturbs * ppt,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi
     )
{

  if (pbi->has_bispectra == _TRUE_) {

    free(pbi->l);
    free(pbi->pk);
    free(pbi->pk_pt);

    for(int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {

      for (int index_l2=0; index_l2 <= index_l1; ++index_l2)
        free (pbi->l3[index_l1][index_l2]);

      free(pbi->l_triangular_size[index_l1]);
      free(pbi->index_l_triangular_min[index_l1]);
      free(pbi->index_l_triangular_max[index_l1]);
      free(pbi->l3_size[index_l1]);
      free(pbi->l3[index_l1]);
          
    } // end of for(index_l1)

    free(pbi->l_triangular_size);
    free(pbi->index_l_triangular_min);
    free(pbi->index_l_triangular_max);
    free(pbi->l3_size);
    free(pbi->l3);    
    

    for(int index_l1=0; index_l1<pbi->l_size; ++index_l1) {
      for(int index_l2=0; index_l2<=index_l1; ++index_l2)
          free (pbi->index_l1_l2_l3[index_l1][index_l1-index_l2]);
      free (pbi->index_l1_l2_l3[index_l1]);
    } free (pbi->index_l1_l2_l3);

    /* Free pbi->bispectra */
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt)
      class_call (bispectra_free_type_level (pbi, index_bt),
        pbi->error_message,
        pbi->error_message);

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
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->d_lsq_cls[index_ct]);
    free (pbi->d_lsq_cls);
    if (ppr->extend_lensed_cls == _TRUE_) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_d_lsq_cls[index_lt]);
      free (pbi->lensed_d_lsq_cls);
    }
    
    if (ppr->extend_lensed_cls == _TRUE_) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_cls[index_lt]);
      free (pbi->lensed_cls);
    }
      
    /* Free file arrays */
    if ((ppr->store_bispectra_to_disk == _TRUE_) || (ppr->load_bispectra_from_disk == _TRUE_)) {
    
      // fclose(pbi->bispectra_status_file);
    
      for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt)
        free (pbi->bispectra_paths[index_bt]);
    
      free (pbi->bispectra_files);
      free (pbi->bispectra_paths);
    
    }
    
  } // end of if(has_bispectra)

  
  return _SUCCESS_;
 
}



/**
 * Free the index_bt level of the bispectra array (pbi->bispectra).
 */

int bispectra_free_type_level(
     struct bispectra * pbi,
     int index_bt
     )
{


  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        free (pbi->bispectra[index_bt][X][Y][Z]);
      free (pbi->bispectra[index_bt][X][Y]);
    }
    free (pbi->bispectra[index_bt][X]);
  }
  free (pbi->bispectra[index_bt]);


  if (pbi->lens_me[index_bt] == _TRUE_) {
    
    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z)
          free (pbi->lensing_correction[index_bt][X][Y][Z]);
        free (pbi->lensing_correction[index_bt][X][Y]);
      }
      free (pbi->lensing_correction[index_bt][X]);
    }
    free (pbi->lensing_correction[index_bt]);
    
  }


  return _SUCCESS_;
 
}








/**
 * Initialize indices and arrays in the second-order transfer functions structure.
 *
 * In detail, this function does:
 *
 *  -# Determine which bispectra to compute and assign them the pbi->index_bt_XXX indices
 *     (for the type: local, cmb-lensing, intrinsic...) and the pbi->index_bf_XXX indices
 *     (for the field: temperature, polarisation...).
 *
 *  -# Determine the common (l1,l2,l3) sampling of all bispectra; for each (l1,l2) couple,
 *     define the extent of the l3 direction based on the triangular condition, with the 
 *     constraint that l1, l2 and l3 all belong to ptr->l (the multipole list where we have
 *     compute the transfer functions).
 *
 *  -# Determine window function for the bispectrum interpolation in harmonic space.
 *
 *  -# Define the integration grid in the k1 and k2 directions for the integration of the
 *     separable and non-separable bispectra.
 *
 *  -# Allocate all levels of the pbi->bispectra array, which will contain the bispectra.
 *
 *  -# Open the files where we will store the bispectra at the end of the computation.
 *
 *  -# Store P(k) (power spectrum of phi) in pbi->pk for all the k-values in the transfer
 *     module (ptr->q) using interpolation, to speed up the integration of the separable
 *     and non-separable bispectra.
 *
 *  -# Store C_l (angular power spectrum of the CMB) in pbi->cls for all the multipoles
 *     between l=2 and l=ptr->l[size-1], to speed up the computation of the analytic
 *     bispectra.
 *
 */

int bispectra_indices (
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

  // =========================================================================================
  // =                                   Count bispectra fields                                     =
  // =========================================================================================

  /* Find out which kind of bispectra to compute and assign them indices and labels
  Generate indices for the probes (T for temperature, E for E-mode polarisation,
  R for Rayleigh...) that we will use to build the bispectra and the Fisher matrix
  elements. */
  
  int index_bf = 0;
  
  pbi->has_bispectra_t = _FALSE_;
  pbi->has_bispectra_e = _FALSE_;
  pbi->has_bispectra_b = _FALSE_;
  for (int i=0; i < _MAX_NUM_FIELDS_; ++i)
    for (int j=0; j < _MAX_LENGTH_LABEL_; ++j)
      pbi->bf_labels[i][j] = '\0';

    
  if (ppt->has_bi_cmb_temperature == _TRUE_) {
    pbi->has_bispectra_t = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "t");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 0;
    pbi->index_bf_t = index_bf++;
  }
  
  if (ppt->has_bi_cmb_polarization == _TRUE_) {
    pbi->has_bispectra_e = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "e");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 2;
    pbi->index_bf_e = index_bf++;
  }

  pbi->bf_size = index_bf;
  pbi->n_probes = pow(pbi->bf_size, 3);

  class_test (pbi->bf_size > 2,
    pbi->error_message,
    "we cannot compute the bispectrum for more than %d fields (e.g. T and E), reduce your expectations :-)", pbi->bf_size);

  class_test (pbi->bf_size < 1,
    pbi->error_message,
    "no probes requested");

  /* Create labels for the full bispectra */
  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {
        for (int i=0; i < _MAX_LENGTH_LABEL_; ++i)
          pbi->bfff_labels[X][Y][Z][i] = '\0';
        sprintf (pbi->bfff_labels[X][Y][Z], "%s%s%s", pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);
      }
    }
  }


  /* Associate to each field X=T,E,... its transfer function, which was computed in the transfer.c
  module, the C_l correlation <X phi> with the lensing potential, the C_l correlation <X zeta> with
  the curvature perturbation and, to each possible pair of fields (TT,EE,TE,...), their power spectra,
  which were computed in the spectra.c module. */
  for (int X = 0; X < pbi->bf_size; ++X) {
    
    if ((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_t;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_tp;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_tt;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_tt;
      if (ppt->has_cl_cmb_zeta == _TRUE_)
        pbi->index_ct_of_zeta_bf[X] = psp->index_ct_tz;
      if (ppr->extend_lensed_cls == _TRUE_)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_tt;
    }

    if ((pbi->has_bispectra_e == _TRUE_) && (X == pbi->index_bf_e)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_e;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_ep;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_te;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_ee;
      if (ppt->has_cl_cmb_zeta == _TRUE_)
        pbi->index_ct_of_zeta_bf[X] = psp->index_ct_ez;
      if (ppr->extend_lensed_cls == _TRUE_)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_ee;
    }

    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      if (((pbi->has_bispectra_t == _TRUE_) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_e == _TRUE_) && (Y == pbi->index_bf_e))) {

        pbi->index_ct_of_bf_bf[X][Y] = pbi->index_ct_of_bf_bf[Y][X] = psp->index_ct_te;
        if (ppr->extend_lensed_cls == _TRUE_)
          pbi->index_lt_of_bf_bf[X][Y] = pbi->index_lt_of_bf_bf[Y][X] = ple->index_lt_te;
      }
    }
  }

  // ================================================================================================
  // =                                    Count bispectra types                                     =
  // ================================================================================================
  
  int index_bt = 0;
  
  pbi->n[separable_bispectrum] = 0;
  pbi->n[non_separable_bispectrum] = 0;
  pbi->n[analytical_bispectrum] = 0;
  pbi->n[intrinsic_bispectrum] = 0;

  for (int i=0; i < _MAX_NUM_BISPECTRA_; ++i)
    for (int j=0; j < _MAX_LENGTH_LABEL_; ++j)
      pbi->bt_labels[i][j] = '\0';


  // *** Separable bispectra
  
  if (pbi->has_local_model) {
    pbi->index_bt_local = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_equilateral_model) {
    pbi->index_bt_equilateral = index_bt;
    strcpy (pbi->bt_labels[index_bt], "equilateral");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_orthogonal_model) {
    pbi->index_bt_orthogonal = index_bt;
    strcpy (pbi->bt_labels[index_bt], "orthogonal");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;
  }

  // *** Non-separable bispectra

  if (pbi->has_galileon_model) {

    /* Bispectrum induced by pi_dot*pi_grad^2 */
    pbi->index_bt_galileon_gradient = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_grad");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;

    /* Bispectrum induced by pi_dot^3 */
    pbi->index_bt_galileon_time = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_time");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;
  }

  // *** Analytical bispectra

  if (pbi->has_local_squeezed == _TRUE_) {
    pbi->index_bt_local_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_intrinsic_squeezed == _TRUE_) {
    pbi->index_bt_intrinsic_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_cosine_shape == _TRUE_) {
    pbi->index_bt_cosine = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cosine");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_cmb_lensing == _TRUE_) {
    pbi->index_bt_cmb_lensing = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  if (pbi->has_cmb_lensing_squeezed == _TRUE_) {
    pbi->index_bt_cmb_lensing_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  /* The kernel for the squeezed CMB-lensing bispectrum is needed to compute
  the lensing contribution to the variance. In the final Fisher matrix, the kernel
  will be multiplied by C_l^{X\phi} to give the actual squeezed bispectrum (see
  eq. 5.20 of http://uk.arxiv.org/abs/1101.2234); it will therefore show up as
  CMB-lensing_sqz rather than kernel_sqz. */
  if (pbi->has_cmb_lensing_kernel == _TRUE_) {
    pbi->index_bt_cmb_lensing_kernel = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  if (pbi->has_quadratic_correction == _TRUE_) {
    pbi->index_bt_quadratic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "quadratic");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    index_bt++;
  }

  // *** Intrinsic (i.e. second-order) bispectra

  if (pbi->has_intrinsic == _TRUE_) {
    pbi->index_bt_intrinsic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic");
    pbi->bispectrum_type[index_bt] = intrinsic_bispectrum;
    pbi->n[intrinsic_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _TRUE_;
    index_bt++;
  }

  pbi->bt_size = index_bt;

  class_test (pbi->bt_size > _MAX_NUM_BISPECTRA_,
   "exceeded maximum number of allowed bispectra, increase _MAX_NUM_BISPECTRA_ in common.h",
   pbi->error_message);

#ifndef WITH_SONG_SUPPORT
   class_test (pbi->has_intrinsic == _TRUE_,
     pbi->error_message,
     "cannot compute the intrinsic bispectrum without SONG support\n");
#endif // WITH_SONG_SUPPORT

  /* Are the Wigner 3j-symbols needed to compute the requested bispectra? */
  pbi->need_3j_symbols = ((pbi->has_bispectra_e) &&
    ((pbi->has_quadratic_correction == _TRUE_) ||
     (pbi->has_cmb_lensing == _TRUE_) ||
     (pbi->has_cmb_lensing_squeezed == _TRUE_) ||
     (pbi->has_cmb_lensing_kernel == _TRUE_)));

  /* No lensing unless the user asked for it explititly */
  if (pbi->has_lensed_bispectra == _FALSE_)
    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
      pbi->lens_me[index_bt] = _FALSE_;
  
  
  
  // =================================================================================================
  // =                                      Determine l-sampling                                     =
  // =================================================================================================
  
  /* We compute the bispectrum on a mesh where l1>=l2>=l3, with l3 determined by the triangular
  condition. All the multipoles are drawn from pbi->l, which is a copy of pbs->l.  */

  pbi->l_size = pbs->l_size;
  class_alloc (pbi->l, pbi->l_size*sizeof(int), pbi->error_message);
  for(int index_l=0; index_l<pbi->l_size; ++index_l)
    pbi->l[index_l] = pbs->l[index_l];

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
  class_alloc (pbi->l3_size, pbi->l_size*sizeof(int *), pbi->error_message);
  class_alloc (pbi->l3, pbi->l_size*sizeof(int **), pbi->error_message);
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
    class_calloc (pbi->l_triangular_size[index_l1], pbi->l_size, sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_min[index_l1], pbi->l_size*sizeof(int), pbi->error_message);
    class_alloc (pbi->index_l_triangular_max[index_l1], pbi->l_size*sizeof(int), pbi->error_message);

    /* We consider only configurations whereby l1>=l2 */
    class_calloc (pbi->l3_size[index_l1], (index_l1+1), sizeof(int), pbi->error_message);
    class_alloc (pbi->l3[index_l1], (index_l1+1)*sizeof(int *), pbi->error_message);
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
      while (pbi->l[index_l_triangular_min] < l_triangular_min)
        ++index_l_triangular_min;

      /* Find the index corresponding to l_triangular_max inside pbi->l */
      class_test (l_triangular_max <= pbi->l[0],
        pbi->error_message,
        "could not fulfill triangular condition for l3 using the multipoles in pbi->l");
      
      int index_l_triangular_max = pbi->l_size-1;
      while (pbi->l[index_l_triangular_max] > l_triangular_max)
        --index_l_triangular_max;    

      /* Fill pbi->index_l_triangular_min and pbi->l_triangular_size */
      pbi->index_l_triangular_min[index_l1][index_l2] = index_l_triangular_min;
      pbi->index_l_triangular_max[index_l1][index_l2] = index_l_triangular_max;
      pbi->l_triangular_size[index_l1][index_l2] = index_l_triangular_max - index_l_triangular_min + 1;

      /* Update counter of triangular configurations */
      pbi->n_total_configurations += pbi->l_triangular_size[index_l1][index_l2];

      /* We shall store the bispectra only for those configurations that simultaneously satisfy 
      the triangular condition and the l1>=l2>=l3 condition. We use the pbi->index_l1_l2_l3 array
      to keep track of the index assigned to a given allowed configuration. */
      if (index_l2 <= index_l1) {

        /* When the triangular condition is not compatible with index_l3<=index_l2, then index_l3_max < index_l3_min+1 will
        be either zero or negative. */ 
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
        int l3_size = MAX (0, index_l3_max-index_l3_min+1);
        pbi->l3_size[index_l1][index_l2] = l3_size;
        class_alloc (pbi->l3[index_l1][index_l2], l3_size*sizeof(int), pbi->error_message);
        class_alloc (pbi->index_l1_l2_l3[index_l1][index_l1-index_l2], l3_size*sizeof(long int), pbi->error_message);
        
        /* The indexing of pbi->index_l1_l2_l3 reflects the l1>=l2>=l3 constraint */
        for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
          pbi->l3[index_l1][index_l2][index_l3-index_l3_min] = pbi->l[index_l3];
          pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3] = index_l1_l2_l3++;
        }

      }
      
      /* Debug - Print out the l3 list for a special configuration */
      // if ((index_l1 == 45) && (index_l2 == 38)) {
      //   if (index_l2 <= index_l1) {
      //
      //     int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      //     int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      //     int l3_size = MAX (0, index_l3_max-index_l3_min+1);
      //     printf ("index_l3_min = %d\n", index_l3_min);
      //     printf ("index_l3_max = %d\n", index_l3_max);
      //
      //     fprintf (stderr, "l1[%d]=%4d, l2[%d]=%4d, l3_size=%4d, l3_min=%4d, l3_max=%4d\n",
      //       index_l1, l1, index_l2, l2, l3_size, pbi->l3[index_l1][index_l2][0], pbi->l3[index_l1][index_l2][l3_size-1]);
      //
      //     for (int index_l3=0; index_l3 < l3_size; ++index_l3)
      //       fprintf(stderr, "%8d %12d\n", index_l3, pbi->l3[index_l1][index_l2][index_l3]);
      //
      //     fprintf (stderr, "\n\n");
      //   }
      // }

      
    } // end of for(index_l2)
  } // end of for(index_l1)

  /* Each bispectrum will be computed for the following number of configurations */
  pbi->n_independent_configurations = index_l1_l2_l3;
  

  /* Inform the user on how much her machine will have to suffer */
  if (pbi->bispectra_verbose > 0) {
    printf_log (" -> we shall compute %dx%d=%d bispectr%s for %ld configurations of (l1,l2,l3)\n",
      pbi->bt_size, pbi->n_probes, pbi->bt_size*pbi->n_probes,
      ((pbi->bt_size*pbi->n_probes)!=1?"a":"um"), pbi->n_independent_configurations);
    // printf("    with (L1,L2,L3) ranging from l=%d to %d (l_size=%d)\n", pbi->l[0], pbi->l[pbi->l_size-1], pbi->l_size);
    if (pbi->bispectra_verbose > 1) {
      printf_log ("     * bispectra types: ");
      for (int index_bt=0; index_bt < (pbi->bt_size-1); ++index_bt) {
        printf_log ("%s, ", pbi->bt_labels[index_bt]);
      }
      printf_log ("%s\n", pbi->bt_labels[pbi->bt_size-1]);
    }
  }
  
  
  

  // ====================================================================================
  // =                           Integration grid in k1 and k2                          =
  // ====================================================================================
  
  /* Sampling of the first-order transfer functions  */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;

  /* Allocate & fill delta_k, which is needed for the trapezoidal integration of both separable and
  non-separable bispectra. This array is has k_size elements, and is defined as k(i+1)-k(i-1)
  except for the first k(1)-k(0) and last k(N)-k(N-1) elements.  Note that when dealing with
  non-separable shapes, we shall not use this array for the integration over k3, as in that
  case the grid in k3 is not fixed but it depends on k1 and k2. */
  class_alloc (pbi->delta_k, k_tr_size * sizeof(double), pbi->error_message);
  
  /* Fill pbi->delta_k */
  pbi->delta_k[0] = k_tr[1] - k_tr[0];
      
  for (int index_k=1; index_k < k_tr_size-1; ++index_k)
    pbi->delta_k[index_k] = k_tr[index_k+1] - k_tr[index_k-1];
      
  pbi->delta_k[k_tr_size-1] = k_tr[k_tr_size-1] - k_tr[k_tr_size-2];
  
  
  
  // ====================================================================================
  // =                            Allocate memory for bispectra                         =
  // ====================================================================================
  
  class_alloc (pbi->bispectra, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
  for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {

    class_alloc (pbi->bispectra[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

    for (int X = 0; X < pbi->bf_size; ++X) {
      
      class_alloc (pbi->bispectra[index_bt][X], pbi->bf_size*sizeof(double **), pbi->error_message);

      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        
        class_alloc (pbi->bispectra[index_bt][X][Y], pbi->bf_size*sizeof(double *), pbi->error_message);
        
        for (int Z = 0; Z < pbi->bf_size; ++Z) 
          class_calloc (pbi->bispectra[index_bt][X][Y][Z], pbi->n_independent_configurations, sizeof(double), pbi->error_message);

      }
    }
  }
  
  pbi->count_allocated_for_bispectra = pbi->bt_size*pbi->n_probes*pbi->n_independent_configurations;


  /* Do the same for the lensing correction */
  
  if (pbi->has_lensed_bispectra == _TRUE_) {

    class_alloc (pbi->lensing_correction, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {

      if (pbi->lens_me[index_bt] == _FALSE_)
        continue;

      class_alloc (pbi->lensing_correction[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

      for (int X = 0; X < pbi->bf_size; ++X) {
      
        class_alloc (pbi->lensing_correction[index_bt][X], pbi->bf_size*sizeof(double **), pbi->error_message);

        for (int Y = 0; Y < pbi->bf_size; ++Y) {
        
          class_alloc (pbi->lensing_correction[index_bt][X][Y], pbi->bf_size*sizeof(double *), pbi->error_message);
        
          for (int Z = 0; Z < pbi->bf_size; ++Z) {
            class_calloc (pbi->lensing_correction[index_bt][X][Y][Z], pbi->n_independent_configurations, sizeof(double), pbi->error_message);
            pbi->count_allocated_for_bispectra += pbi->n_independent_configurations;

          }
        }
      }
    }
  }

  
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * allocated ~ %.3g MB (%ld doubles) for the bispectra array\n",
    pbi->count_allocated_for_bispectra*sizeof(double)/1e6, pbi->count_allocated_for_bispectra);


  
  // ====================================================================================
  // =                               Create storage files                               =
  // ====================================================================================
  
  /* Create the files to store the bispectra in */
  if ((ppr->store_bispectra_to_disk == _TRUE_) || (ppr->load_bispectra_from_disk == _TRUE_)) {

    /* We are going to store the bispectra in n=bt_size files, one for each requested type of bispectrum */
    class_alloc (pbi->bispectra_files, pbi->bt_size*sizeof(FILE *), pbi->error_message);
    class_alloc (pbi->bispectra_paths, pbi->bt_size*sizeof(char *), pbi->error_message);
  
    for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      
      /* Include the name of the bispectrum in its file */
      class_alloc (pbi->bispectra_paths[index_bt], _FILENAMESIZE_*sizeof(char), pbi->error_message);
      sprintf (pbi->bispectra_paths[index_bt], "%s/bispectra_%s.dat", pbi->bispectra_dir, pbi->bt_labels[index_bt]);
      
    } // end of for(index_bt)

    if (ppr->store_bispectra_to_disk == _TRUE_)
      printf_log_if (pbi->bispectra_verbose, 1, 
        "     * will create %d files for the bispectra\n", pbi->bt_size);
      
  } // end of if(ppr->store_bispectra_to_disk)
  
  

  // ====================================================================================
  // =                                 Interpolate P(k)                                 =
  // ====================================================================================
  
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


  // ====================================================================================
  // =                               Interpolate the C_l                                =
  // ====================================================================================

  /* Interpolate the Cl's in all l-values */
  class_call (bispectra_cls(
                ppr,
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
 * time potential psi. The power spectrum is computed in the points contained in ptr->q,
 * as these are the points where the transfer functions are computed.
 *
 * The purpose of this function is twofold. First, it stores the primordial power spectrum
 * into memory for faster access by the bispectra module, as calling primordial_spectrum_at_k()
 * is fairly expensive. Secondly, we convert the dimensionless spectrum of the curvature
 * perturbation R outputted by the primordial module of CLASS, Delta_R(k), into the power
 * spectrum for the Newtonian curvature potential, P_Phi(k). The two spectra are related by:
 *  
 * P_Phi(k) = 2*Pi^2/k^3 * Delta_R(k)
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
  int k_pt_size = ppt->k_size[ppt->index_md_scalars];
  class_alloc (pbi->pk_pt, k_pt_size*sizeof(double), pbi->error_message);
  
  /* Fill pk with the values of the primordial power spectrum, as obtained in the ppm module */

  for (int index_k_pt=0; index_k_pt<k_pt_size; ++index_k_pt) {
    
    double k_pt = ppt->k[ppt->index_md_scalars][index_k_pt];

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


/**
 * Interpolate the C_l and their derivatives from the spectra.c module, and store them
 * in the bispectra structure for easy access.
 *
 * This function allocates and fills the following fields:
 * -# pbi->cls
 * -# pbi->d_lsq_cls
 * -# pbi->lensed_cls
 * -# pbi->lensed_d_lsq_cls
 */
int bispectra_cls (
    struct precision * ppr,
    struct perturbs * ppt,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi
    )
{

  // ====================================================================================
  // =                                    Allocate arrays                               =
  // ====================================================================================

  /* Allocate the array that will contain the C_l's for all types and all l's. */
  class_alloc (pbi->cls, psp->ct_size*sizeof(double*), pbi->error_message);
  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    class_calloc (pbi->cls[index_ct], pbi->full_l_size, sizeof(double), pbi->error_message);
  
  /* Do the same for the logarithmic derivative of the C_l's */
  class_alloc (pbi->d_lsq_cls, psp->ct_size*sizeof(double*), pbi->error_message);
  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
    class_calloc (pbi->d_lsq_cls[index_ct], pbi->full_l_size, sizeof(double), pbi->error_message);

  /* If they have been computed, also store the lensed C_l up to l_max */
  if (ppr->extend_lensed_cls == _TRUE_) {

    class_alloc (pbi->lensed_cls, ple->lt_size*sizeof(double*), pbi->error_message);
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
      class_calloc (pbi->lensed_cls[index_lt], pbi->full_l_size, sizeof(double), pbi->error_message);
    
    /* The squeezed approximation of the intrinsic bispectrum requires the computation of the
    C_l derivatives */
    class_alloc (pbi->lensed_d_lsq_cls, ple->lt_size*sizeof(double*), pbi->error_message);
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
      class_calloc (pbi->lensed_d_lsq_cls[index_lt], pbi->full_l_size, sizeof(double), pbi->error_message);
  }
  
  
  /* We shall call the CLASS function 'spectra_cl_at_l'. This gives three outputs:
  the total Cl's for each probe (T, E, B...); the Cl's divided by probe and mode
  (scalar, vector, tensor); the Cl's divided by probe, mode, and initial condition
  (adiabatic, isocurvature...). We have to allocate these three arrays before being
  able to call 'spectra_cl_at_l'. We do so copying what is done in the function
  'output_total_cl_at_l' in the output module */

  double * cl;        /* cl_md_ic[index_ct] */
  double ** cl_md;    /* cl_md[index_mode][index_ct] */
  double ** cl_md_ic; /* cl_md_ic[index_mode][index_ic1_ic2*psp->ct_size+index_ct] */
  class_alloc(cl, MAX(psp->ct_size,ple->lt_size)*sizeof(double), pbi->error_message);	
  class_alloc(cl_md_ic, psp->md_size*sizeof(double *), pbi->error_message);
  class_alloc(cl_md, psp->md_size*sizeof(double *), pbi->error_message);
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {
    if (psp->md_size > 1)
      class_alloc(cl_md[index_mode], psp->ct_size*sizeof(double), pbi->error_message);	
    if (psp->ic_size[index_mode] > 1)
      class_alloc(cl_md_ic[index_mode], psp->ic_ic_size[index_mode]*psp->ct_size*sizeof(double), pbi->error_message);
  }
  
  // ====================================================================================
  // =                                     Store C_l                                    =
  // ====================================================================================

  for (int l=2; l<=pbi->l_max; ++l) {
    
    class_call(spectra_cl_at_l(
                 psp,
                 (double)l,
                 cl,
                 cl_md,
                 cl_md_ic),
      psp->error_message,
      pbi->error_message);
      
    /* Store the total Cl's into an array as a function of l and probe. By 'total' we mean
    the power spectrum summed over all the modes and initial conditions */
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      pbi->cls[index_ct][l-2] = cl[index_ct];
    

    /* Store the lensed C_l up to l_max, if they have been computed. Note that for the
    lensing potential C_l we shall always use the unlensed array pbi->cls. */
    if (ppr->extend_lensed_cls == _TRUE_) {

        /* The lensed C_l's must have been computed up to the maximum required multipole,
        that is, up to pbi->l_max. This check is made inside lensing_cl_at_l() */
        class_call(lensing_cl_at_l(
                     ple,
                     l,
                     cl),
          ple->error_message,
          pbi->error_message);
    
        for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
          pbi->lensed_cls[index_lt][l-2] = cl[index_lt];
        
        /* Debug - print temperature-lensing potential C_l */
        // double factor = l*(l+1.)/(2*_PI_);
        // fprintf (stderr, "%4d %16g %16g %16g %16g\n",
        //   l, factor*pbi->cls[psp->index_ct_tt][l-2], factor*pbi->lensed_cls[ple->index_lt_tt][l-2],
        //   factor*sqrt(l*(l+1))*pbi->cls[psp->index_ct_tp][l-2], factor*sqrt(l*(l+1))*pbi->lensed_cls[ple->index_lt_tp][l-2]);
        
    }
    
    /* Uncomment to turn the CMB-lensing C_l to zero on small scales, where we cannot
    trust them. This won't change the result because these C_l are very small for large
    l, but it will considerably speed up the lensing variance computations in the Fisher
    module. In CAMB, Antony sets C_l^TP=0 for l>300 and C_l^EP=0 for l>40. */
    if ((pbi->has_cmb_lensing == _TRUE_)
    || (pbi->has_cmb_lensing_squeezed == _TRUE_)
    || (pbi->has_cmb_lensing_kernel == _TRUE_)) {

      pbi->lmax_lensing_corrT = 300;
      pbi->lmax_lensing_corrE = 300;

      if ((l > pbi->lmax_lensing_corrT) && (pbi->has_bispectra_t == _TRUE_))
        pbi->cls[psp->index_ct_tp][l-2] = 0;

      if ((l > pbi->lmax_lensing_corrE) && (pbi->has_bispectra_e == _TRUE_))
        pbi->cls[psp->index_ct_ep][l-2] = 0;

    }

    /* To compute the squeezed limit approximation case, we need the derivative of 
    l*l*C_l. The best way to obtain these is to take the derivative of the
    spline-interpolated C_l's that we have just stored in pbi->cls. The reason is
    that the C_l's are smooth and the cubic interpolation does a very good job at
    filling the gaps in psp->cl, which contains only the C_l's computed by the spectra.c
    module. On the other hand, directly taking the derivative of the sparsely sampled C_l's
    in psp->cl will give a bad numerical convergence. We have tested that the first
    approach (derivative of spline-interpolated C_l's) converges much better than
    the second one (derivative of computed C_l's), with respect to the number of points
    taken in the l-grid. Uncomment the following lines to use the first approach. */

    //   class_call(spectra_dcl_at_l(
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
    //     pbi->d_lsq_cls[index_ct][l-2] = cl[index_ct];
    // 
    //   /* Some debug - print the cls */
    //   // printf ("%d %g %g %g\n", l,
    //   //   pbi->cls[psp->index_ct_tt][l-2],
    //   //   pbi->cls[psp->index_ct_tz][l-2],
    //   //   pbi->d_lsq_cls[psp->index_ct_tt][l-2]);
    
  } // end of for loop on the l's

  /* Compute the derivative of l*l*C_l from the spline-interpolated C_l's. This corresponds
  to the second method mentioned in the long comment above. */
  double * l_array;
  class_alloc (l_array, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);
  for (int l=2; l<=pbi->l_max; ++l)
    l_array[l-2] = l;

  double * lsq_cl, * dd_lsq_cl;
  class_alloc (lsq_cl, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);
  class_alloc (dd_lsq_cl, (pbi->l_max-2+1)*sizeof(double), pbi->error_message);

  for (int index_ct=0; index_ct < psp->ct_size; ++index_ct) {

    for (int l=2; l<=pbi->l_max; ++l)
      lsq_cl[l-2] = l*l*pbi->cls[index_ct][l-2];
  
    /* Compute the second derivatives of l^2*C_l */
    class_call(array_spline_table_lines(
                 l_array,
                 pbi->l_max-2+1,
                 lsq_cl,
                 1,
                 dd_lsq_cl,
                 _SPLINE_EST_DERIV_,
                 pbi->error_message),
      pbi->error_message,
      pbi->error_message);
      
    /* Compute the first derivative using the above information */
    class_call (array_spline_derive_table_lines(
                  l_array,
                  pbi->l_max-2+1,
                  lsq_cl,
                  dd_lsq_cl,
                  1,
                  pbi->d_lsq_cls[index_ct],
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);

  } // end of loop on index_ct
  
  /* Do the same for the lensed C_l's */
  if (ppr->extend_lensed_cls == _TRUE_) {
  
    for (int index_lt=0; index_lt < ple->lt_size; ++index_lt) {

      for (int l=2; l<=pbi->l_max; ++l)
        lsq_cl[l-2] = l*l*pbi->lensed_cls[index_lt][l-2];
  
      /* Compute the second derivatives of l^2*C_l */
      class_call(array_spline_table_lines(
                   l_array,
                   pbi->l_max-2+1,
                   lsq_cl,
                   1,
                   dd_lsq_cl,
                   _SPLINE_EST_DERIV_,
                   pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      
      /* Compute the first derivative using the above information */
      class_call (array_spline_derive_table_lines(
                    l_array,
                    pbi->l_max-2+1,
                    lsq_cl,
                    dd_lsq_cl,
                    1,
                    pbi->lensed_d_lsq_cls[index_lt],
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);

    } // end of loop on index_lt
  } // end of if lensing
  
  free (l_array);
  free(lsq_cl);
  free(dd_lsq_cl);

  /* Free memory */
  for (int index_mode = 0; index_mode < psp->md_size; index_mode++) {    
    if (psp->md_size > 1) free(cl_md[index_mode]);  
    if (psp->ic_size[index_mode] > 1) free(cl_md_ic[index_mode]);
  }  
  free(cl_md_ic);
  free(cl_md);
  free(cl);
  
  return _SUCCESS_;
  
}



/**
 * Compute the CMB bispectra as a function of (l1,l2,l3) and store them
 * in pbi->bispectra.
 *
 * In detail, this function does:
 *
 * -# Compute and store the separable bispectra (local, equilateral, orthogonal, etc.) in
 * bispectra_separable_init(). 
 *
 * -# Compute and store the analytical bispectra (cmb-lensing, squeezed-limit approximations,
 * etc.) in bispectra_analytical_init().
 *
 * -# Compute and store the non-separable bispectra (galileon or any arbitrary shape) in
 * bispectra_non_separable_init().
 *
 * -# Chec that the bispectra thus computed and stored in pbi->bispectra do not contain
 * invalid entries (nans).
 */

int bispectra_harmonic (
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
                  pth,
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

    /* IMPORTANT: No symmetrization is needed, as the primary bispectrum is automatically
    symmetrised if we feed a primordial bispectrum B(k1,k2,k3) symmetrised in k, as we do.
    For the second-order bispectrum, we will need to perform the symmetrisation manually. */

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
                  ple,
                  pbi),
      pbi->error_message,
      pbi->error_message);
      
  }
  
  /* If we are loading the bispectra from disk, nothing else needs to be done */

  if (ppr->load_bispectra_from_disk == _TRUE_) {
        
    return _SUCCESS_;

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
                  pth,
                  ppt,
                  pbs,
                  ptr,
                  ppm,
                  pbi,
                  pwb_nonsep),
      pbi->error_message,
      pbi->error_message);
  
  
    /* Free the 'pwb_nonsep' workspace */
    class_call (bispectra_non_separable_workspace_free(
                  pbi,
                  pwb_nonsep),
      pbi->error_message,
      pbi->error_message);
  }
  
  /* Print information on memory usage */
  if ((pbi->bispectra_verbose > 1) && (pbi->n[intrinsic_bispectrum]) < 1)
    printf(" -> memorised ~ %.3g MB (%ld doubles) in the bispectra array\n",
      pbi->count_memorised_for_bispectra*sizeof(double)/1e6, pbi->count_memorised_for_bispectra);
  
  
  
  // ============================================================================
  // =                      Check bispectra against nan's                       =
  // ============================================================================

  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* We have not computed the intrinsic bispectra yet, so we skip them */
    if (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {
    
          #pragma omp parallel for
          for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {
            for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {
     
              /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
              int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
              int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
     
              for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

                /* Index of the current (l1,l2,l3) configuration */
                long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
          
                double bispectrum = pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
   
                if (isnan(bispectrum))
                  printf ("@@@ WARNING: b(%d,%d,%d) = %g for bispectrum '%s_%s'.\n",
                  pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], bispectrum,
                  pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

              } // end of for(index_l3)
            } // end of for(index_l2)
          } // end of for(index_l1)
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)
  } // end of for(index_bt)
  
  return _SUCCESS_;

}



/**
 * Produce output files for the bispectrum.
 *
 * Three types of files will be created:
 *
 * - A text file with the bispectra configurations for l1=l1_out and
 *   l2=l2_out, tabulated as a function of l3, named bispectra_1D_LXXX.txt.
 *   If l2_out<0, then the configurations with l1=l1_out and l2=l3 will
 *   be printed to file.
 *
 * - A larger text file with the bispectra configurations for l1=l1_out, tabulated
 *   as a function of l2 and l3, named bispectra_2D_LXXX.txt.
 *
 * - A binary file with all bispectra configurations for all bispectra types,
 *   named bispectra.dat.
 *
 */

int bispectra_output (
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
  
  /* Load the bispectra from disk if they were not already computed */

  if (ppr->load_bispectra_from_disk == _TRUE_) {
  
    for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt)

      if ((pbi->bispectrum_type[index_bt] == non_separable_bispectrum) ||
          (pbi->bispectrum_type[index_bt] == intrinsic_bispectrum))

        class_call (bispectra_load_from_disk (
                      pbi,
                      index_bt),
          pbi->error_message,
          pbi->error_message);    
  
  }


  // ====================================================================================
  // =                                     2D output                                    =
  // ====================================================================================

  /* The 2D output consists in a text file with the bispectrum tabulated as a function
  of (l2,l3) for a given l1. We generate one file for each probe (TTT, EEE, TTE...). */

  for (int index_l_out=0; index_l_out < ppr->l_out_size; ++index_l_out) {

    for (int X = 0; X < pbi->bf_size; ++X) {

      for (int Y = 0; Y < pbi->bf_size; ++Y) {

        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          /* Multipole value and index of this output file */
          int l1 = ppr->l1_out[index_l_out];
          int index_l1 = ppr->index_l1_out[index_l_out];

          /* Build filenames */ 
          int index_probe = X*pbi->bf_size*pbi->bf_size + Y*pbi->bf_size + Z;
          sprintf (ppr->l_out_paths_1D[index_l_out][index_probe], "%sbispectra_1D_L%03d_%s.txt",
            ppr->l_out_paths_1D[index_l_out][index_probe], index_l_out, pbi->bfff_labels[X][Y][Z]);
          sprintf (ppr->l_out_paths_2D[index_l_out][index_probe], "%sbispectra_2D_L%03d_%s.txt",
            ppr->l_out_paths_2D[index_l_out][index_probe], index_l_out, pbi->bfff_labels[X][Y][Z]);

          /* Open files */
          FILE * file_1D = ppr->l_out_files_1D[index_l_out][index_probe];
          FILE * file_2D = ppr->l_out_files_2D[index_l_out][index_probe];
          class_open(file_1D, ppr->l_out_paths_1D[index_l_out][index_probe], "w", pbi->error_message);
          class_open(file_2D, ppr->l_out_paths_2D[index_l_out][index_probe], "w", pbi->error_message);
          

          // -------------------------------------------------------------------------------
          // -                               Print information                             -
          // -------------------------------------------------------------------------------
          
          char line[1024];
          
          /* Write the information header of the 1D and 2D files */
          if (ppr->l2_out[index_l_out] > 0)
            sprintf (line, "CMB reduced bispectra b_l1_l2_l3 tabulated as a function of l3 and bispectrum type for a fixed (l1,l2) pair.");
          else
            sprintf (line, "CMB reduced bispectra b_l1_l_l tabulated as a function of l for a fixed l1 value.");
          fprintf (file_1D, "%s%s\n", _COMMENT_, line);
          sprintf (line, "CMB reduced bispectra b_l1_l2_l3 tabulated as a function of (l2,l3) and bispectrum type for a fixed l1 value.");
          fprintf (file_2D, "%s%s\n", _COMMENT_, line);
          sprintf (line, "This file was generated by SONG %s (%s) on %s.", _SONG_VERSION_, _SONG_URL_, ppr->date);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s\n", _COMMENT_);
          if (pbi->has_lensed_bispectra == _TRUE_) {
            sprintf (line, "The suffix _u denotes unlensed bispectra.");
            fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
            fprintf_2way (1, file_1D, 0, file_2D, 0, "%s\n", _COMMENT_);
          }
          sprintf (line, "The column 'norm' is 6/5 * [cl1_%sZ*(cl2_%s%s + cl3_%s%s) + cl2_%sZ*(cl3_%s%s + cl1_%s%s) + cl3_%sZ*(cl1_%s%s + cl2_%s%s)],",
            pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Y], pbi->bf_labels[Z], pbi->bf_labels[Z],
            pbi->bf_labels[Y], pbi->bf_labels[Z], pbi->bf_labels[Z], pbi->bf_labels[X], pbi->bf_labels[X],
            pbi->bf_labels[Z], pbi->bf_labels[X], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Y]);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
          sprintf (line, "where Z is the curvature perturbation zeta at recombination.");
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
          sprintf (line, "The column 'norm_positive' is -24/5 * (cl1_%s%s*cl2_%s%s + cl2_%s%s*cl3_%s%s + cl3_%s%s*cl1_%s%s).",
            pbi->bf_labels[X], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Y],
            pbi->bf_labels[Y], pbi->bf_labels[Y], pbi->bf_labels[Z], pbi->bf_labels[Z],
            pbi->bf_labels[Z], pbi->bf_labels[Z], pbi->bf_labels[X], pbi->bf_labels[X]);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s\n", _COMMENT_);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s", pba->info);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s\n", _COMMENT_);
          fprintf (file_1D, "%sInformation on the output (l1,l2):\n", _COMMENT_);
          fprintf (file_1D, "%sl1 = %d, index_l1 = %d/%d\n", _COMMENT_, l1, index_l1, pbi->l_size-1);
          fprintf (file_2D, "%sInformation on the output l1:\n", _COMMENT_);
          fprintf (file_2D, "%sl1 = %d, index_l1 = %d/%d\n", _COMMENT_, l1, index_l1, pbi->l_size-1);

          /* For the 2D file, we shall print a row for each (l2,l3) configuration; for the
          1D file, a row for each l3 configuration */

          int n_rows_1D = 0;
          int n_rows_2D = 0;

          for (int index_l2=0; index_l2 < pbi->l_size; ++index_l2) {

            int l2 = pbi->l[index_l2];
            int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
            int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
            int l3_size = index_l3_max - index_l3_min + 1;

            for (int index_l3=index_l3_min; index_l3 <= index_l3_max; ++index_l3) {

              int l3 = pbi->l[index_l3];

              // -------------------------------------------------------------------------------
              // -                              Extract bispectra                              -
              // -------------------------------------------------------------------------------
              
              /* We are going to print all (l2,l3) configurations even though we computed the
              bispectrum only for l1>=l2>=l3. We obtain the values outside that range using
              the symmetry of the bispectrum with respect to permutations of (l1,l2,l3). This
              approach won't work for any squeezed-limit approximations, such as the i_squeezed
              and l_squeezed bispectra, because they are asymmetric by construction. For these
              bispectra, we set the output to zero outside the l1>=l2>=l3 regime. */

              double bispectrum[pbi->bt_size];
              double bispectrum_unlensed[pbi->bt_size];
  
              for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
  
                class_call (bispectra_at_l1l2l3 (
                              pbi,
                              index_bt,
                              index_l1, index_l2, index_l3,
                              X, Y, Z,
                              &bispectrum[index_bt],
                              &bispectrum_unlensed[index_bt]),
                  pbi->error_message,
                  pbi->error_message);

              }

              
              // -------------------------------------------------------------------------------
              // -                                 Normalisation                               -
              // -------------------------------------------------------------------------------
              
              /* Include normalisation factors in the output files. The bispectra divided by
              this normalisation will be of order unity. */

              double normalisation = 0;

              class_call (bispectra_normalisation (
                            ppr, psp, ple, pbi,
                            l1, l2, l3,
                            X, Y, Z,
                            0, 0, 0,
                            &normalisation),
                pbi->error_message,
                pbi->error_message);

              if (fabs(normalisation) < _MINUSCULE_)
                printf ("WARNING: normalisation=%g is small; beware of inf\n", normalisation);

              double normalisation_positive = 0;

              class_call (bispectra_normalisation_positive (
                            ppr, psp, ple, pbi,
                            l1, l2, l3,
                            X, Y, Z,
                            0, 0, 0,
                            &normalisation_positive),
                pbi->error_message,
                pbi->error_message);

              if (fabs(normalisation_positive) < _MINUSCULE_)
                printf ("WARNING: normalisation_positive=%g is small; beware of inf\n", normalisation_positive);


              // -------------------------------------------------------------------------------
              // -                                   Build row                                 -
              // -------------------------------------------------------------------------------

              /* Arrays containing all the information on the columns to be printed, labels included */
              char label[_MAX_NUM_COLUMNS_][_MAX_LENGTH_LABEL_];
              double value[_MAX_NUM_COLUMNS_];
              short condition[_MAX_NUM_COLUMNS_];
  
              /* Initialise column arrays */
              for (int i=0; i < _MAX_NUM_COLUMNS_; ++i)
                condition[i] = _TRUE_;

              /* Shortcut for file verbosity */
              int v = 1;
  
              /* Initialise column counter  */
              int i = -1;

              /* Multipole l2 (won't be printed on the 1D file) */
              strcpy (label[++i], "l2");
              value[i] = pbi->l[index_l2];
              
              /* Multipole l3 */
              if (ppr->l2_out[index_l_out] > 0)
                strcpy (label[++i], "l3");
              else
                strcpy (label[++i], "l");
              value[i] = pbi->l[index_l3];
              
              /* Bispectra */
              for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

                sprintf (label[++i], "%s", pbi->bt_labels[index_bt]);
                value[i] = bispectrum[index_bt];

                if (pbi->lens_me[index_bt] == _TRUE_) {
                  sprintf (label[++i], "%s_u", pbi->bt_labels[index_bt]);
                  value[i] = bispectrum_unlensed[index_bt];
                }

              }

              /* Normalisation */
              sprintf (label[++i], "norm");
              value[i] = normalisation;

              /* Strictly positive normalisation */
              sprintf (label[++i], "norm_positive");
              value[i] = normalisation_positive;


              // -------------------------------------------------------------------------------
              // -                             Print row to 2D file                             -
              // -------------------------------------------------------------------------------

              /* Maximum number of columns that will be written */
              int n_max_columns = i+1;  
              class_test (n_max_columns > _MAX_NUM_COLUMNS_,
                pbi->error_message,
                "too many columns; raise _MAX_NUM_COLUMNS_ to at least %d",
                _MAX_NUM_COLUMNS_);

              /* Choose how label & values should be formatted */
              char format_label[64] = "%18s(%02d) ";
              char format_value[64] = "%22.11g ";

              /* Write row with labels */
              int n_columns_2D = 1;
              if (n_rows_2D++ == 0) {
                for (int i=0; i < n_max_columns; ++i)
                  if (condition[i])
                    fprintf (file_2D, format_label, label[i], n_columns_2D++);
                fprintf (file_2D, "\n");
              }

              /* Write row with data to file */
              for (int i=0; i < n_max_columns; ++i)
                if (condition[i])
                  fprintf (file_2D, format_value, value[i]);
              fprintf (file_2D, "\n");


              // -------------------------------------------------------------------------------
              // -                             Print row to 1D file                             -
              // -------------------------------------------------------------------------------

              if (index_l2 == ppr->index_l2_out[index_l_out]) {

                /* Write row with labels and append information on l2 to the header */
                int n_columns_1D = 1;
                if (n_rows_1D++ == 0) {
                  fprintf (file_1D, "%sl2 = %d, index_l2 = %d/%d\n", _COMMENT_, l2, index_l2, pbi->l_size-1);
                  fprintf (file_1D, "%sl3 spans %d values from %d to %d\n", _COMMENT_, l3_size, pbi->l[index_l3_min], pbi->l[index_l3_max]);
                  for (int i=1; i < n_max_columns; ++i) /* skip the first (l2) column */
                    if (condition[i])
                      fprintf (file_1D, format_label, label[i], n_columns_1D++);
                  fprintf (file_1D, "\n");
                }

                /* Write row with data to file */
                for (int i=1; i < n_max_columns; ++i)  /* skip the first (l2) column */
                  if (condition[i])
                    fprintf (file_1D, format_value, value[i]);
                fprintf (file_1D, "\n");

              }
              
              /* If the user gave -1 as a value for l2_out, then print the configuration with
              l1=l1_out and l2=l3 */
              
              else if ((ppr->l2_out[index_l_out] < 0) && (l2 == l3)) {
                
                /* Write row with labels and append information on l2 to the header */
                int n_columns_1D = 1;
                if (n_rows_1D++ == 0) {
                  fprintf (file_1D, "%swill print b(%d,l,l) as a function of l\n", _COMMENT_, l1, pbi->l_size);
                  for (int i=1; i < n_max_columns; ++i) /* skip the first (l2) column */
                    if (condition[i])
                      fprintf (file_1D, format_label, label[i], n_columns_1D++);
                  fprintf (file_1D, "\n");
                }

                /* Write row with data to file */
                for (int i=1; i < n_max_columns; ++i)  /* skip the first (l2) column */
                  if (condition[i])
                    fprintf (file_1D, format_value, value[i]);
                fprintf (file_1D, "\n");                
                
              }

              
            } // for l3
          } // for l2


          /* Close the file */
          fclose (file_2D);

        } // Z
      } // Y
    } // X

  } // for l_out
  
  

  // ====================================================================================
  // =                                     3D output                                    =
  // ====================================================================================

  if (pbi->output_binary_bispectra == _TRUE_) {

    /* Binary files produced by this function will have a human-readable ASCII
    header. It will include information on the cosmological parameters and on the
    bispectrum, plus a binary map useful to understand how to access the data
    in the binary file. */
    int header_size = 
      _MAX_INFO_SIZE_ + /* For background information */
      _MAX_INFO_SIZE_ + /* For other information */
      _MAX_INFO_SIZE_ + _MAX_HEADER_LINE_LENGTH_*pbi->n_probes; /* For binary map */

    /* Define a new binary file structure */
    struct binary_file * file;
    class_alloc (file, sizeof(struct binary_file), pbi->error_message);
  
    /* Open the binary file */
    class_call (binary_init (
                  file,
                  &(ppr->l_out_file_3D),
                  ppr->l_out_path_3D,
                  "w",
                  header_size,
                  ppr->output_single_precision),
      file->error_message,
      pbi->error_message);

  
    // ---------------------------------------------------------------------------
    // -                            Print information                            -
    // ---------------------------------------------------------------------------

    /* Add information to the file header */
    binary_sprintf (file, "CMB reduced bispectra b_l1_l2_l3 tabulated as a function of (l1,l2,l3) and bispectrum type.");
    binary_sprintf (file, "This binary file was generated by SONG (%s) on %s.", _SONG_URL_, ppr->date);
    binary_sprintf (file, "");
    sprintf (file->header, "%s%s", file->header, pba->info);
    file->header_size += strlen (pba->info) + 1;
    if (ppt->gauge == newtonian) binary_sprintf(file, "gauge = newtonian");
    if (ppt->gauge == synchronous) binary_sprintf(file, "gauge = synchronous");


    // --------------------------------------------------------------------------
    // -                                Build blocks                             -
    // ---------------------------------------------------------------------------

    char desc[1024];
    char name[1024];

    int label_size = _MAX_LENGTH_LABEL_;
    sprintf (desc, "length of a label (=%d)", _MAX_LENGTH_LABEL_);
    sprintf (name, "_MAX_LENGTH_LABEL_");
    class_call (binary_append_int (file, &label_size, 1, desc, name),
      file->error_message,
      pbi->error_message);

    /* l sampling */

    sprintf (desc, "size of l array (=%d)", pbi->l_size);
    sprintf (name, "pbi->l_size");
    class_call (binary_append_int (file, &pbi->l_size, 1, desc,name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "l array");
    sprintf (name, "pbi->l");
    class_call (binary_append_int (file, pbi->l, pbi->l_size, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "size of l3 grid: l3_size[index_l1][index_l2] with index_l1 < pbi->l_size, index_l2 <= index_l1");
    sprintf (name, "pbi->l3_size");
    int index_l3_size_block = file->n_blocks;
  
    for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {
      for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
        class_call (binary_add_block (
                      file,
                      &pbi->l3_size[index_l1][index_l2],
                      1,
                      sizeof (int),
                      desc,
                      "int",
                      name,
                      index_l3_size_block),
          file->error_message,
          pbi->error_message);
      }
    }

    sprintf (desc, "l3 array: l3[index_l1][index_l2] with index_l1 < pbi->l_size, index_l2 <= index_l1");
    sprintf (name, "pbi->l3");
    int index_l3_block = file->n_blocks;
  
    for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1)
      for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
                
        class_call (binary_add_block (
                      file,
                      pbi->l3[index_l1][index_l2],
                      pbi->l3_size[index_l1][index_l2],
                      sizeof (int),
                      desc,
                      "int",
                      name,
                      index_l3_block),
          file->error_message,
          pbi->error_message);
    }

    sprintf (desc, "number of (l1,l2,l3) configurations in each bispectrum (=%d)", pbi->n_independent_configurations);
    sprintf (name, "pbi->n_independent_configurations");
    class_call (binary_append_long_int (file, &pbi->n_independent_configurations, 1, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "array of indices associated to the (l1,l2,l3) configurations");
    sprintf (name, "pbi->index_l1_l2_l3");
    class_call (binary_append_int (file, pbi->index_l1_l2_l3, pbi->n_independent_configurations, desc, name),
      file->error_message,
      pbi->error_message);

    /* First-order C_l */

    sprintf (desc, "number of l-values for which the C_l are stored (=%d)", pbi->full_l_size);
    sprintf (name, "pbi->full_l_size");
    class_call (binary_append_int (file, &pbi->full_l_size, 1, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "number of C_l types (=%d)", psp->ct_size);
    sprintf (name, "psp->ct_size");
    class_call (binary_append_int (file, &psp->ct_size, 1, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "array of names of C_l types (each has %d char)", _MAX_LENGTH_LABEL_);
    sprintf (name, "psp->ct_labels");
    class_call (binary_append_string (file, psp->ct_labels, psp->ct_size*_MAX_LENGTH_LABEL_, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "unlensed C_l: cls[index_ct] with index_ct < psp->ct_size");
    sprintf (name, "pbi->cls");
    int index_unlensed_cl = file->n_blocks;
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      class_call (binary_add_block (
                    file,
                    pbi->cls[index_ct],
                    pbi->full_l_size,
                    sizeof (double),
                    desc,
                    "double",
                    name,
                    index_unlensed_cl),
        file->error_message,
        pbi->error_message);

    sprintf (desc, "unlensed C_l logarithmic derivative: d_lsq_cls[index_ct] with index_ct < psp->ct_size");
    sprintf (name, "pbi->d_lsq_cls");
    index_unlensed_cl = file->n_blocks;
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      class_call (binary_add_block (
                    file,
                    pbi->d_lsq_cls[index_ct],
                    pbi->full_l_size,
                    sizeof (double),
                    desc,
                    "double",
                    name,
                    index_unlensed_cl),
        file->error_message,
        pbi->error_message);

    if (ppr->extend_lensed_cls == _TRUE_) {

      sprintf (desc, "lensed C_l: cls[index_ct] with index_ct < psp->ct_size");
      sprintf (name, "pbi->lensed_cls");
      int index_lensed_cl = file->n_blocks;
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        class_call (binary_add_block (
                      file,
                      pbi->lensed_cls[index_lt],
                      pbi->full_l_size,
                      sizeof (double),
                      desc,
                      "double",
                      name,
                      index_lensed_cl),
          file->error_message,
          pbi->error_message);

      sprintf (desc, "lensed C_l logarithmic derivative: d_lsq_cls[index_ct] with index_ct < psp->ct_size");
      sprintf (name, "pbi->lensed_d_lsq_cls");
      index_lensed_cl = file->n_blocks;
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        class_call (binary_add_block (
                      file,
                      pbi->lensed_d_lsq_cls[index_lt],
                      pbi->full_l_size,
                      sizeof (double),
                      desc,
                      "double",
                      name,
                      index_lensed_cl),
          file->error_message,
          pbi->error_message);
    }


    /* Bispectrum */

    sprintf (desc, "number of bispectra types (local, equilateral, intrinsic...) (=%d)", pbi->bt_size);
    sprintf (name, "pbi->bt_size");
    class_call (binary_append_int (file, &pbi->bt_size, 1, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "array of names of bispectra types (each has %d char)", _MAX_LENGTH_LABEL_);
    sprintf (name, "pbi->bt_labels");
    class_call (binary_append_string (file, pbi->bt_labels, pbi->bt_size*_MAX_LENGTH_LABEL_, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "number of fields (T,E,B...) (=%d)", pbi->bf_size);
    sprintf (name, "pbi->bf_size");
    class_call (binary_append_int (file, &pbi->bf_size, 1, desc, name),
      file->error_message,
      pbi->error_message);

    sprintf (desc, "array of names of fields (each has %d char)", _MAX_LENGTH_LABEL_);
    sprintf (name, "pbi->bf_labels");
    class_call (binary_append_string (file, pbi->bf_labels, pbi->bf_size*_MAX_LENGTH_LABEL_, desc, name),
      file->error_message,
      pbi->error_message);

    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

      sprintf (desc, "%s CMB bispectrum for all values of (l1,l2,l3) and for all fields (X,Y,Z)", pbi->bt_labels[index_bt]);
      sprintf (name, "pbi->bispectra[index_bt=%d]", index_bt);
      int index_bispectrum = file->n_blocks;      

      for (int X = 0; X < pbi->bf_size; ++X) {
        for (int Y = 0; Y < pbi->bf_size; ++Y) {
          for (int Z = 0; Z < pbi->bf_size; ++Z) {

            class_call (binary_add_block (
                          file,
                          pbi->bispectra[index_bt][X][Y][Z],
                          pbi->n_independent_configurations,
                          sizeof (double),
                          desc,
                          "double",
                          name,
                          index_bispectrum),
              file->error_message,
              pbi->error_message);

          } // Z
        } // Y
      } // X
    } // bt
    
  
    // ---------------------------------------------------------------------------
    // -                              Write to file                              -
    // ---------------------------------------------------------------------------

    class_call (binary_write (
                  file),
      file->error_message,
      pbi->error_message);


    // ---------------------------------------------------------------------------
    // -                                Clean up                                 -
    // ---------------------------------------------------------------------------

    class_call (binary_free (
                  file),
      file->error_message,
      pbi->error_message);
      
  } // if output_binary_bispectra


  return _SUCCESS_;
  
}





/**
 * Determine the integration grid in r.
 *
 * In the bispectrum integral, the integration variable r appears both as a
 * prefactor and in the argument of the three Bessel functions, which in turn
 * are multiplied by three transfer functions:
 * 
 *  b_l1_l2_l3(r) ~ (k1*k2*k3*r)^2 * j_l1(r*k1) * j_l2(r*k2) * j_l3(r*k3)
 *                                 * T_l1(k1)   * T_l2(k2)   * T_l3(k3)
 *                                 * F(k1,k2,k3)
 *
 * Please refer to Fergusson and Shellard 2007 for the full formula and its
 * derivation.
 *
 * The transfer function T_l(k) is the convolution of the source function
 * S(k,tau) with the Bessel function j_l(k*(tau0-tau)). The product of this
 * Bessel with j_l(k*r) is maximised when they both have the same frequency,
 * that is, for r~tau0-tau, because both are oscillatory functions.
 *
 * This means that any recombination physics will show up in the bispectrum
 * integral close to r~tau0-tau_rec, where the visibility function is
 * strongly peaked.
 *
 * For the same reason we expect late-time effects such as ISW and reionisation
 * to contribute especially to small r. What we see, however, is that the
 * contribution from r~tau0-tau_rec dominates also for the late time effects.
 * We find two reasons for this behaviour:
 *
 * - The r^2 factor penalises the contributions from small r.
 *
 * - Each Bessel function exponentially suppresses contributions from r<l/k.
 *
 * - Contributions from r>l/k are penalised by a factor 1/r only, rather than
 *   exponentially; this allows the late-time effects to overlap in r with the
 *   recombination effects and thus generate a correlation.
 *
 * Let us consider the late ISW effect as an example. We expect its effect to
 * be maximised when l1, l2 and l3 are all smaller than 10, because that's
 * where it shows up in the transfer functions. Indeed, we see that b_2_2_2(r)
 * for r<4000 is completely dominated by the late ISW effect (local template).
 * However, this contribution is still four orders of magnitude smaller than 
 * b_2_2_2(r) around r~tau0-tau_rec, which means that values of r smaller
 * than tau0-tau_rec can be safely excluded from the bispectrum integral, even
 * where the late ISW is strongest.
 *
 * On the other hand, the late ISW does contribute significantly to b_2_2_2(r) 
 * around r=tau0-tau_rec, as it dampens its amplitude by roughly 40%. This 
 * "early time" contribution by a late-time effect is explained by the r^2
 * boost in the bispectrum integral and by the correlation between the late
 * ISW effect and the recombination effects. The correlation, which manifests
 * itself in the triple product of Bessel and transfer functions, is made
 * possible by the fact that the Bessel functions j_l(k*r) penalise only mildly
 * the configurations with r > l/k. (Note that while the late time effects
 * can show up at high r, the recombination effects cannot show up at small r,
 * because the Bessel functions suppress r < l/k exponentially.)
 *
 * In absence of recombination effects, the bispectrum will be made up entirely
 * of late time effects. This is never the case for temperature, which has
 * recombination sources even on superhorizon scales (SW effect), but it does
 * happen for polarisation, which has no sources on scales that were
 * superhorizon at recombination (l<200). On these scales, the only polarised
 * signal is generated by reionisation, in a similar way to what happens for
 * the C_l^EE. We see this directly in the EEE bispectrum for l1=l2=l3 = 2,
 * which has a smooth peak at r~tau0-tau_reio (~11000 Mpc for a standard LCDM
 * model) with a tiny protuberance at r~tau0-tau_rec. For l1=l2=l3 = 10, the
 * two peaks are almost merged, while for higher l the recombination peak
 * dominates over the reionisation one.
 * 
 * In summary:
 *
 * - For temperature, the late time effects (ISW and reionisation) will
 *   contribute to the bispectrum mostly via their correlation with the
 *   recombination effects (SW, doppler, early ISW). Therefore, r needs
 *   to be sampled only around r~tau0-tau_rec.
 * 
 * - For polarisation, reionisation generates a signal on large scales (l<40)
 *   that does not correlate with the recombination sources (l>200), because
 *   there are none. Therefore, the dominant effect from reionisation does not
 *   sit at r~tau0-tau_rec, and r should be sampled also for r<tau0-tau_rec.
 *
 * I have tested that the above considerations hold for the local template, the
 * equilateral template and the intrinsic bispectrum. For the equilateral
 * template, it is advised to extend r_max to at least tau0+10*tau_rec. For the
 * intrinsic TTT bispectrum with l1,l2,l3 < O(10), the small-r contribution
 * from the ISW effect is comparable to the r~tau0-tau_rec contribution. This
 * is probably an effect of mode coupling that mixes small and large scales.
 * However, the net effect on the SN of the intrinsic bispectrum is negligible.
 *
 * One last thing to take into account at second-order. To correctly compute
 * the bispectrum in the r<tau0-tau_rec regime, and thus capture the effect
 * of reionisation on polarisation and the tiny ISW effect, make sure to use
 * a fine enough time sampling for the line of sight integral at late times.
 * This can be achieved by either increasing the sampling of the source functions
 * via the parameter perturb_sampling_late_time_boost or by increasing the line
 * of sight integration grid with the parameter tau_linstep_song. The first way
 * turns out to be less computationally expensive.
 */

int bispectra_get_r_grid (
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
    struct perturbs * ppt,
    struct bispectra * pbi,
    double * tau_sampling, /**< input, time sampling for the line-of-sight sources. 
                           Used only if ppr->bispectra_r_sampling==sources_r_sampling */
    int tau_size, /**< input, size of tau_sampling. Used only if
                  ppr->bispectra_r_sampling==sources_r_sampling */
    double ** r_grid, /**< output, integration grid of size r_size */
    int * r_size, /**< output, size of the integration grid */
    double * r_min, /**< output, minimum value in the r-grid */
    double * r_max, /**< output, maximum value in the r-grid */
    double ** delta_r /**< output, trapezoidal measure of the r-grid */
    )
{

  /* - Sample r linearly between the given limits */
  
  if (ppr->bispectra_r_sampling == custom_r_sampling) {

    *r_size = ppr->r_size;
    *r_min = ppr->r_min;
    *r_max = ppr->r_max;

    /* Sample r linearly */
    class_alloc (*r_grid, *r_size*sizeof(double), pbi->error_message);
    lin_space (*r_grid, *r_min, *r_max, *r_size);
  }

  /* - Centre the r-grid on tau0-tau_rec */

  else if (ppr->bispectra_r_sampling == centred_r_sampling) {

    *r_size = ppr->r_size;
    double centre = pba->conformal_age - pth->tau_rec;
    *r_min = MAX (0, centre - ppr->r_left*pth->tau_rec);
    *r_max = centre + ppr->r_right*pth->tau_rec;
    
    /* Sample r linearly */
    class_alloc (*r_grid, *r_size*sizeof(double), pbi->error_message);
    lin_space (*r_grid, *r_min, *r_max, *r_size);
    
  }

  /* - Determine the r array as tau0-tau */

  else if (ppr->bispectra_r_sampling == sources_r_sampling) {

    class_alloc (*r_grid,
      TAU_SIZE_MAX*sizeof(double),
      pbi->error_message);

    /* Find range in which we shall sample r densely */

    double centre = pba->conformal_age - pth->tau_rec;
    double r_left = MAX (0, centre - ppr->r_left*pth->tau_rec);
    double r_right = centre + ppr->r_right*pth->tau_rec;
    
    /* Determine the time step between r_left and r_right based on user input */

    double tau_step = (r_right - r_left)/(ppr->r_size-1);

    /* Determine the time step automatically */

    if (ppr->r_size < 0) {

      int index_tau_rec = 0;
      while ((index_tau_rec < (tau_size-1)) && (tau_sampling[index_tau_rec] < pth->tau_rec))
        ++index_tau_rec;

      class_test ((index_tau_rec==0) || (index_tau_rec==(tau_size-1)),
        pbi->error_message,
        "something wrong in the tau sampling");

      tau_step = tau_sampling[index_tau_rec] - tau_sampling[index_tau_rec-1];
    }

    /* First, add points to the left drawing from the sources sampling */

    int index_tau = tau_size-2;
    int index_r = 0;
    (*r_grid)[0] = MIN (r_left, pba->conformal_age - tau_sampling[index_tau]);

    while ((*r_grid)[index_r] < r_left) {
      --index_tau;
      ++index_r;
      (*r_grid)[index_r] = pba->conformal_age - tau_sampling[index_tau];
    }
    
    /* Then, add linearly-spaced points until r_right */
    
    while ((*r_grid)[index_r] < r_right) {
      ++index_r;
      (*r_grid)[index_r] = (*r_grid)[index_r-1] + tau_step;
    }
    
    *r_size = index_r+1;
    
    class_realloc (*r_grid,
      *r_grid,
      *r_size*sizeof(double),
      pbi->error_message);
    
    *r_min = (*r_grid)[0];
    *r_max = (*r_grid)[*r_size-1];
    
    /* Debug - Check the parameters of the grid */
    // printf ("centre = %g\n", centre);
    // printf ("r_left = %g\n", r_left);
    // printf ("r_right = %g\n", r_right);
    // printf ("pth->tau_rec = %g\n", pth->tau_rec);
    // printf ("tau_step = %g\n", tau_step);
    // printf ("*r_size = %d\n", *r_size);

  }
  
  /* Debug - Print r-integration grid */
  // printf ("# ~~~ r-sampling for the bispectrum integration ~~~\n");
  // for (int index_r=0; index_r < *r_size; ++index_r) {
  //   printf ("%12d %16g\n", index_r, (*r_grid)[index_r]);
  // }

  /* Issue a warning if the integration grid starts too late to sample the effect of
  reionisation on polarisation */
  double tau_reio;
  class_call (background_tau_of_z (pba,pth->z_reio,&tau_reio),
    pba->error_message,
    pbi->error_message);

  if ((pth->reio_parametrization != reio_none)
  && ((pbi->has_bispectra_e == _TRUE_) || (pbi->has_bispectra_b == _TRUE_)))
    class_test_permissive (*r_min > (pba->conformal_age-pth->tau_reio),
      pbi->error_message,
      "the r-sampling might be inadequate to compute reionisation");

  /* Check that the r-grid is strictly ascending */
  for(int index_r=0; index_r<(*r_size-1); ++index_r) {
    class_test ((*r_grid)[index_r] >= (*r_grid)[index_r+1],
      pbi->error_message,
      "the r grid should be stricty ascending");
  }


  /* Allocate & fill delta_r, the measure for the trapezoidal integration over r */
  class_alloc ((*delta_r), *r_size * sizeof(double), pbi->error_message);
  (*delta_r)[0] = (*r_grid)[1] - (*r_grid)[0];
  for (int index_r=1; index_r < *r_size-1; ++index_r)
    (*delta_r)[index_r] = (*r_grid)[index_r+1] - (*r_grid)[index_r-1];
  (*delta_r)[*r_size-1] = (*r_grid)[*r_size-1] - (*r_grid)[*r_size-2];
  
  
  return _SUCCESS_;
  
}



/**
 * Initialise the workspace that will be required to integrate the separable bispectra.
 */

int bispectra_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_separable * pwb
    )
{

  // ====================================================================================
  // =                             Prepare integration grid                             =
  // ====================================================================================

  class_call (bispectra_get_r_grid (
                ppr,
                pba,
                pth,
                ppt,
                pbi,
                ppt->tau_sampling,
                ppt->tau_size,
                &(pwb->r),
                &(pwb->r_size),
                &(pwb->r_min),
                &(pwb->r_max),
                &(pwb->delta_r)),
    pbi->error_message,
    pbi->error_message);


  // =====================================================================================================
  // =                                  Allocate memory for filter functions                             =
  // =====================================================================================================
  
  /* Allocate the separable integrals needed for each requested model of primordial non-Gaussianity.
  Refer to the header file for details on the models. The integrand arrays need to be written
  by different threads at the same time, hence we allocate one of them for each thread. */
      
  /* Parallelization variables */
  int number_of_threads = 1;
  int thread = 0;
  int abort = _FALSE_;
  
  #pragma omp parallel
  #ifdef _OPENMP
  number_of_threads = omp_get_num_threads();
  #endif

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

      printf_log_if (pbi->bispectra_verbose, 2, 
        "      \\ computing filter functions for l=%d, index_l=%d\n", pbi->l[index_l], index_l);

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

        for (int index_k=0; index_k < k_size; ++index_k) {
        
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
        for (int index_k=0; index_k < k_size; ++index_k) {      

          double pk_one_third = pow(pbi->pk[index_k], 1/3.);
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
      
        printf_log_if (pbi->bispectra_verbose, 3, 
          "       \\ r=%g, index_r=%d\n", r, index_r);
        
        
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
                        k_size,
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
                        k_size,
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
                        k_size,
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
                        k_size,
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
    int X, // index for the first field (T,B,E)
    int Y, // index for the second field (T,B,E)
    int Z, // index for the third field (T,B,E)
    struct bispectra_workspace_separable * pwb
    )

{

  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * computing the r-integral for the bispectrum %s_%s%s%s\n",
    pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);

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
  
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      printf_log_if (pbi->bispectra_verbose, 2, 
        "      \\ computing r-integral for l1=%d, index_l1=%d\n", pbi->l[index_l1], index_l1);
    
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
                  alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                );
            
            } // end of local model

            // -------------------------------------------------------------------------
            // -                           Equilateral model                           -
            // -------------------------------------------------------------------------
            
            else if ((pbi->has_equilateral_model == _TRUE_) && (index_bt == pbi->index_bt_equilateral)) {
  
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                - beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                - beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                    
                /* Symmetrical part */
                - 2 * delta[X][index_l1][index_r] * delta[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + beta[X][index_l1][index_r]  * gamma[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * delta[Y][index_l2][index_r] * gamma[Z][index_l3][index_r]
                + gamma[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * delta[Z][index_l3][index_r]
                + gamma[X][index_l1][index_r] * delta[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                + delta[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * gamma[Z][index_l3][index_r]
                + delta[X][index_l1][index_r] * gamma[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
              
              );

            } // end of equilateral model

            // -------------------------------------------------------------------------------
            // -                              Orthogonal model                               -
            // -------------------------------------------------------------------------------
            
            else if ((pbi->has_orthogonal_model == _TRUE_) && (index_bt == pbi->index_bt_orthogonal)) {
  
              /* We take the formula from Senatore et al. 2010, also shown in Komatsu et al. 2011 (WMAP7 paper, eq. 64) */
              integrand = 6 * fnl_R * r*r * (

                /* Local part */
                - 3 * alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                - 3 * beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                - 3 * beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                    
                /* Symmetrical part. We found what we think is a typo in eq. 38 of http://arxiv.org/abs/1006.0275v3,
                where the coefficient is -2/3*3 = -2 instead of -8. We think -8 is the correct coefficient, as it can
                be verified from eq. 3.2 of Senatore, Smith & Zaldarriaga 2010.  */
                - 8 * delta[X][index_l1][index_r] * delta[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                    
                /* Completely asymmetrical part */
                + 3 * beta[X][index_l1][index_r]  * gamma[Y][index_l2][index_r] * delta[Z][index_l3][index_r]
                + 3 * beta[X][index_l1][index_r]  * delta[Y][index_l2][index_r] * gamma[Z][index_l3][index_r]
                + 3 * gamma[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * delta[Z][index_l3][index_r]
                + 3 * gamma[X][index_l1][index_r] * delta[Y][index_l2][index_r] * beta[Z][index_l3][index_r]             
                + 3 * delta[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * gamma[Z][index_l3][index_r]
                + 3 * delta[X][index_l1][index_r] * gamma[Y][index_l2][index_r] * beta[Z][index_l3][index_r]             
              
              );

            } // end of orthogonal model

            
            /* Increment the estimate of the integral */
            integral += integrand * pwb->delta_r[index_r];

            /* Debug - Print integrand as a function of r */
            // if (index_bt == pbi->index_bt_local)
            //   if ((X==pbi->index_bf_t) && (Y==pbi->index_bf_t) && (Z==pbi->index_bf_t))
            //     if ((pbi->l[index_l1] == 2) && (pbi->l[index_l2] == 2) && (pbi->l[index_l3] == 2))
            //       fprintf (stderr, "%15.7g %15.7g\n", r, integrand);

  
          } // end of for(index_r)
  

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3), including the factor 1/2
          from trapezoidal rule */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = 0.5 * integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. In KSW2005, this factor
          was split between the alpha and beta integrals, but from the numerical point of view it is preferable
          to include it at the end of the computation (as in eq. 17 of Fergusson & Shellard 2007). */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] *= pow(2./_PI_,3);
  
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
 * Compute the CMB bispectra that are separable in k-space.
 *
 * Separable bispectra are characterised by a primordial shape function S(k1,k2,k3)
 * that can be expressed as the product of three functions:
 *     S(k1,k2,k3) = A(k1)*B(k2)*C(k3).
 * As a result, the integration in Fourier space of separable bispectra can be reduced
 * to three simple one-dimensional integrals, which can be solved quickly. Examples of
 * separable bispectra are the local, equilateral and orthogonal primordial templates
 * (see Planck paper, http://arxiv.org/abs/1303.5084).
 */
int bispectra_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
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
                pth,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);


  printf_log_if (pbi->bispectra_verbose, 1, 
    " -> computing separable bispectra; r sampled %d times in [%g,%g]\n",
    pwb->r_size, pwb->r_min, pwb->r_max);
    




  // ===================================================================================
  // =                             Compute filter functions                            =
  // ===================================================================================

  /* Compute the filter functions: alpha(l,r), beta(l,r), gamma(l,r), delta(l,r) for each field 
  (T,E,B). The filter functions are convolution of P(k)^A*T(k)^B with a Bessel function. */

  for (int index_bf=0; index_bf < pbi->bf_size; ++index_bf) {

    printf_log_if (pbi->bispectra_verbose, 1, 
      "     * computing %s filter functions for the separable bispectra ...\n",
      pbi->bf_labels[index_bf]);

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

    printf_log_if (pbi->bispectra_verbose, 1, 
      "     * integrating the filter functions over r ...\n");

    /* Loop on the considered probe (TTT, TTE, TEE, EEE...) */
    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {
  
          class_call (bispectra_separable_integrate_over_r(
                        ppr,
                        pba,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_bt,
                        X,
                        Y,
                        Z,
                        pwb),
            pbi->error_message,
            pbi->error_message);
  
        } // end of for(X)
      } // end of for(Y)
    } // end of for(Z)
    
  } // end of for(index_bt)
  
  
  return _SUCCESS_; 
  
}



/**
 * Free the memory associated to the workspace used for the computation of the
 * separable bispectra.
 */
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
 
  free (pwb);
 
  return _SUCCESS_; 
  
}






/**
 * Compute the analytic CMB bispectra.
 *
 * Analytic bispectra are obtained using a closed formula usually involving
 * three-j symbols and angular power spectra C_l. They do not require a numerical
 * integrations and therefore they are quicker to compute. Examples are the
 * cmb-lensing bispectrum and any squeezed-limit approximation.
 *
 * The user has to provide the formula for the analytic bispectra by writing actual
 * functions that take l1, l2 and l3 as arguments.
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
    struct lensing * ple,
    struct bispectra * pbi
    )
{  

  printf_log_if (pbi->bispectra_verbose, 0, 
    " -> computing analytic bispectra\n");

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

  // ===================================================================================
  // =                           Define analytic bispectra                             =
  // ===================================================================================

  /* Associate to each analytic bispectrum a function. In order to add your custom
  function, just add a line here, */
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    pbi->bispectrum_function[index_bt] = NULL;

    if ((pbi->has_cmb_lensing == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_bispectrum;

    else if ((pbi->has_cmb_lensing_squeezed == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_bispectrum;

    else if ((pbi->has_cmb_lensing_kernel == _TRUE_) && (index_bt == pbi->index_bt_cmb_lensing_kernel))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_kernel;

    else if ((pbi->has_local_squeezed == _TRUE_) && (index_bt == pbi->index_bt_local_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_local_squeezed_bispectrum;
    
    else if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_intrinsic_squeezed_bispectrum;
     
    else if ((pbi->has_intrinsic_squeezed_unlensed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed_unlensed))
      pbi->bispectrum_function[index_bt] = bispectra_intrinsic_squeezed_unlensed_bispectrum;
     
    else if ((pbi->has_quadratic_correction == _TRUE_) && (index_bt == pbi->index_bt_quadratic))
      pbi->bispectrum_function[index_bt] = bispectra_quadratic_correction;
    
    else if ((pbi->has_cosine_shape == _TRUE_) && (index_bt == pbi->index_bt_cosine))
      pbi->bispectrum_function[index_bt] = bispectra_cosine_bispectrum;

  }
  

  // ===================================================================================
  // =                                   Main loop                                     =
  // ===================================================================================
    
  /* We parallelize the outer loop over 'l1'. */
  #pragma omp parallel for shared (ppt,pbs,ptr,ppm,pbi,abort) private (thread)
  for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif

    int l1 = pbi->l[index_l1];

    for (int index_l2 = 0; index_l2 <= index_l1; ++index_l2) {

      int l2 = pbi->l[index_l2];

      /* Skip those configurations that are forbidden by the triangular condition (optimization) */
      if (l2 < l1/2)
        continue;

      /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);

      /* The reduced bispectrum is defined as the angular averaged bispectrum divided
      by 3J(l1,l2,l3)(0,0,0). In some cases, it is possible to perform this ratio 
      analytically, and obtain a reduced bispectrum that is defined also for odd
      values if l1+l2+l3. When polarisation is involved, the ratio cannot always be
      performed analytically in an obvious way, eg. for the cmb lensing bispectrum.
      Doing it numerically would imply throwing away all odd configurations of
      l1+l2+l3. Rather, we use the recursive relation in threej_ratio_M_recursive()
      to express the ratio of 3J as a sum of analytic functions. Here we test that
      the result does not change whether we use the recursive relation or the 
      numerical ratio between 3J symbols. */        
      // int M=4;
      // double threej_num[2*pbi->l_max+1], threej_den[2*pbi->l_max+1];
      // int l3_min_num, l3_min_den;
      // double min_D, max_D;
      // class_call_parallel (drc3jj (
      //                        MAX(l1,M), MAX(l2,M), 0, -M,
      //                        &min_D, &max_D,
      //                        &(threej_num[0]),
      //                        (2*pbi->l_max+1),
      //                        pbi->error_message),
      //   pbi->error_message,
      //   pbi->error_message);
      // l3_min_num = (int)(min_D + _EPS_);
      // class_call_parallel (drc3jj (
      //                        MAX(l1,M), MAX(l2,M), 0, 0,
      //                        &min_D, &max_D,
      //                        &(threej_den[0]),
      //                        (2*pbi->l_max+1),
      //                        pbi->error_message),
      //   pbi->error_message,
      //   pbi->error_message);
      // l3_min_den = (int)(min_D + _EPS_);
      //
      // for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
      //   int l3 = pbi->l[index_l3];
      //   if ((l1+l2+l3)%2==0) {
      //     if ((l1<M) || (l2<M) || (l3<M)) continue;
      //     double * ratio = malloc(sizeof(double)*(M+1));
      //     class_call_parallel (threej_ratio_M_recursive(l1, l2, l3, M, ratio, pbi->error_message),
      //       pbi->error_message, pbi->error_message);
      //     double result_1 = threej_num[l3-l3_min_num]/threej_den[l3-l3_min_den];
      //     double result_2 = ratio[M];
      //     double frac = 1-result_1/result_2;
      //     class_test_parallel (fabs(frac) > _SMALL_,
      //       pbi->error_message,
      //       "(%3d,%3d,%3d,M=%d), res_1=%14.6g, res_2=%14.6g, diff=%14.6g\n",
      //       l1, l2, l3, M, result_1, result_2, frac);
      //   }
      // }
          

      // -----------------------------------------------------------------------------------
      // -                               Compute bispectra                                 -
      // -----------------------------------------------------------------------------------

      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {  

        /* Index of the current (l1,l2,l3) configuration */
        int l3 = pbi->l[index_l3];
        long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
        
        /* Compute 3J ratios, needed only for polarised bispectra such as CMB-lensing
        and the quadratic correction. */
        double threej_ratio_20m2, threej_ratio_m220, threej_ratio_0m22;

        if (pbi->need_3j_symbols == _TRUE_) {          

          class_call_parallel (threej_ratio_M (l2, l1, l3, 2, &threej_ratio_20m2, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l3, l1, l2, 2, &threej_ratio_m220, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l1, l2, l3, 2, &threej_ratio_0m22, pbi->error_message),
            pbi->error_message, pbi->error_message);

        } // end of 3j computation

        for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

          if (pbi->bispectrum_type[index_bt] != analytical_bispectrum)
            continue;
      
          /* Check that the current bispectrum has a function associated to it */
          class_test_parallel (pbi->bispectrum_function[index_bt]==NULL,
            pbi->error_message,
            "no function associated for the bispectrum '%s'. Maybe it's not analytical?",
            pbi->bt_labels[index_bt]);

          for (int X = 0; X < pbi->bf_size; ++X) {
            for (int Y = 0; Y < pbi->bf_size; ++Y) {
              for (int Z = 0; Z < pbi->bf_size; ++Z) {

                  /* Compute the bispectrum using the function associated to index_bt */
                  class_call_parallel ((*pbi->bispectrum_function[index_bt]) (
                                ppr, psp, ple, pbi,
                                l1, l2, l3,
                                X, Y, Z,
                                threej_ratio_20m2,
                                threej_ratio_m220,
                                threej_ratio_0m22,
                                &(pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3])),
                    pbi->error_message,
                    pbi->error_message);
      
                  /* Update the counter */
                  increase_counter:;
                  #pragma omp atomic
                  pbi->count_memorised_for_bispectra++;

              } // end of for(Z)
            } // end of for(Y)
          } // end of for(X)
        } // end of for(index_l3)
      } // end of for(index_l2)
    } // end of for(index_bt)
    #pragma omp flush(abort)
  } // end of for(index_l1)
  if (abort == _TRUE_) return _FAILURE_;  // end of parallel region

  return _SUCCESS_;

}






/**
 * Compute the CMB bispectra that have an arbitrary form in k-space.
 *
 * The non-separable bispectra can be obtained only by brute force, that is, convolving
 * the arbitrary shape function S(k1,k2,k3) with the Bessel functions and the transfer
 * functions inside a 3D integral, according to the formula in Eq. 12 of Fergusson
 * and Shellard 2007.
 *
 * The user provides the shape functions by writing actual functions that take
 * k1, k2 and k3 as arguments. The integration of the shape functions is then
 * broken down in the following steps:
 *
 * -# Determine the k-sampling for the shape functions and the r-sampling for
 * the bispectrum integral in bispectra_non_separable_workspace_init().
 * 
 * -# For each non-separable bispectrum, determine its shape function and let the
 * pointer pbi->shape_function point to it.
 * 
 * -# Integrate the shape function S(k1,k2,k3,r) in the k3 direction via
 * bispectra_non_separable_integrate_over_k3(); in the future, this will
 * be done in an optimised way that takes into account the triangular 
 * condition.
 *
 * -# Integrate the resulting array I_l3(k1,k2,r) in the k2 direction via
 * bispectra_non_separable_integrate_over_k2().
 *
 * -# Integrate the resulting array I_l2_l3(k1,r) in the k1 direction via
 * bispectra_non_separable_integrate_over_k1().
 *
 * -# Finally, in bispectra_non_separable_integrate_over_r() integrate the
 * resulting array I_l1_l2_l3(r) in the r direction and store the resulting
 * bispectrum b_l1_l2_l3 in pbi->bispectra.
 *
 * When integrating the k-directions, we always assume that the primordial shape functions
 * are smoother than the transfer functions in k1, k2 and k3; in particular, we assume that
 * their features are well captured by the k-sampling in ppt->k. There are certain models
 * of cosmic inflation for which this might not be true; in these cases, make sure to modify 
 * the sampling accordingly in bispectra_non_separable_workspace_init().
 *
 * Note that the hard number grinding is done by the bessel_convolution() function in the
 * bessel.c module.
 */
int bispectra_non_separable_init (
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  /* Allocate arrays inside the integration workspace; in particular, 
  determine the k-grid for */
  class_call (bispectra_non_separable_workspace_init(
                ppr,
                pba,
                pth,
                ppt,
                pbs,
                ptr,
                ppm,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);
  

  printf_log_if (pbi->bispectra_verbose, 1, 
    " -> computing non-separable bispectra; r sampled %d times in [%g,%g]\n",
    pwb->r_size, pwb->r_min, pwb->r_max);
  
  for (int index_bt = 0; index_bt < pbi->bt_size; ++index_bt) {

    /* Skip the bispectrum if it not of the non-separable type */
    if (pbi->bispectrum_type[index_bt] != non_separable_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {

      pwb->X = X;

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

      printf_log_if (pbi->bispectra_verbose, 0, 
        " -> computing %s bispectrum involving %s^(2)\n",
        pbi->bt_labels[index_bt], pbi->bf_labels[X]);

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
                    pbi->index_tt_of_bf[X],
                    pwb),
        pbi->error_message,
        pbi->error_message);
        
      for (int Y = 0; Y < pbi->bf_size; ++Y) {

        pwb->Y = Y;

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
                      pbi->index_tt_of_bf[Y],
                      pwb),
          pbi->error_message,
          pbi->error_message);

        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          pwb->Z = Z;

          printf_log_if (pbi->bispectra_verbose, 1, 
            "   \\ computing bispectrum %s_%s%s%s\n",
            pbi->bt_labels[index_bt], pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);

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
                        pbi->index_tt_of_bf[Z],
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
                        pbi->bispectra[index_bt][X][Y][Z],
                        pwb),
            pbi->error_message,
            pbi->error_message);

        } // end of for(Z)
      } // end of for(Y)
    } // end of for(X)
  } // end of for(index_bt)
  
  return _SUCCESS_;
  
}


/**
 * Initialise the workspace that will be required to integrate the non-separable
 * bispectra.
 *
 * In particular:
 *
 * -# Set up the integration grid in the time direction r.
 *
 * -# Determine the sampling in k of the primordial shape function S(k1,k2,k3); in
 * doing so, we assume that S is smoother than the transfer functions.
 *
 * -# For each (k1,k2) pair, enforce the limits imposed by the triangular condition
 * on the integration grid in the k3 direction.
 *
 * -# Determine the window function to apply to the bispectrum integral during the
 * integration in k2 and k1; the window function is supposed to make the function
 * smoother and therefore easier to interpolate.
 *
 */
int bispectra_non_separable_workspace_init (
    struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{

  // ====================================================================================
  // =                             Prepare integration grid                             =
  // ====================================================================================


  // ------------------------------------------------------------
  // -                         Grid in r                        -
  // ------------------------------------------------------------
  
  class_call (bispectra_get_r_grid (
                ppr,
                pba,
                pth,
                ppt,
                pbi,
                ppt->tau_sampling,
                ppt->tau_size,
                &(pwb->r),
                &(pwb->r_size),
                &(pwb->r_min),
                &(pwb->r_max),
                &(pwb->delta_r)),
    pbi->error_message,
    pbi->error_message);


  // -----------------------------------------------------------------------
  // -                       Grid for the shape function                   -
  // -----------------------------------------------------------------------
  
  /* We shall sample the primordial shape function only in the k-points determined
  in the perturbation module, which are much sparsely sampled than those of the
  transfer functions. This means that we are assuming that the shape function is a
  smooth function of (k1,k2,k3) */
  
  pwb->k_smooth_size = ppt->k_size[ppt->index_md_scalars];
  
  class_alloc (pwb->k_smooth_grid, pwb->k_smooth_size*sizeof(double), pbi->error_message);
  
  for (int index_k=0; index_k < pwb->k_smooth_size; ++index_k)
    pwb->k_smooth_grid[index_k] = ppt->k[ppt->index_md_scalars][index_k];
    
  
  
  // -----------------------------------------------------------------------
  // -                     Enforce triangular condition                    -
  // -----------------------------------------------------------------------

  /* Here we set the integration limits on k3. These can be made equal to those of k1
  and k2 (that is equal to the range where we computed the transfer functions) but we
  can do better than that. In fact, the r-integral enforces the triangular condition
  |k1-k2| <= k3 <= k1+k2 which means that we can restrict our range to those configurations.
  However, we should not get too close to the triangular limits otherwise the integral
  becomes numerically unstable. */

  /* Sampling of the first-order transfer functions.  */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  pwb->k3_size_max = k_tr_size;
  
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
      
      /* Uncomment to take the whole k3 range for the integration (safe approach) */
      double k_tr_min = k_tr[0];
      double k_tr_max = k_tr[k_tr_size-1];
      
      /* Uncomment to restrict the integration to the k3 values that satisfy the triangular
      condition |k1-k2| <= k3 <= k1+k2. This will give faster but imprecise result. The reason
      is that the Dirac delta function that enforces the triangular condition, \delta(k1+k2-k3),
      is not explicit in the integral. In fact, we expanded \delta in Bessel functions using the 
      Rayleigh expansion of a plane wave (see eq. 6.30 of my thesis http://arxiv.org/abs/1405.2280
      or eqs. 9 and 10 of Fergusson and Shellard 2007). Analytically, this does not make any
      difference, but numerically it renders the integral unstable unless we take extra oscillations
      in the k3 direction. This means that to make the integral stable, we should extend the range
      a bit, similarly to what we do for the second-order bispectrum in bispectra2.c. */
      // double k_tr_min = fabs(k1-k2);
      // double k_tr_max = k1+k2;
      
      /* Find the index corresponding to k3_lower inside k_tr */
      int index_k3_lower = 0;
      while (k_tr[index_k3_lower] < k_tr_min) ++index_k3_lower;
      pwb->index_k3_lower[index_k1][index_k2] = index_k3_lower;

      /* Find the index corresponding to k3_upper inside ptr->q */
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
  
  /* Determine the window function for the interpolation of the k-space bispectrum in k1
  and k2. We need a window function because the primordial bispectrum will usually contain
  products of power spectra that diverge as k^-3 for k->0. */
    
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
  #ifdef _OPENMP
  number_of_threads = omp_get_num_threads();
  #endif
  
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
  
    /* Allocate the integration grid in k3 with the maximum possible number of k3-values. This is
    given by the number of k-values in ptr->q */
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
  
  for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
    
    /* Allocate r-level */
    class_alloc (pwb->integral_over_k3[index_l3], pwb->r_size*sizeof(double **), pbi->error_message);
  
    /* Allocate 'k1' level */
    for (int index_r=0; index_r < pwb->r_size; ++index_r) {
  
      int k1_size = pwb->k_smooth_size;
      class_alloc (pwb->integral_over_k3[index_l3][index_r], k1_size*sizeof(double *), pbi->error_message);

      /* Allocate 'k2' level */
      for (int index_k1=0; index_k1<k1_size; ++index_k1) {
  
        int k2_size = index_k1 + 1;
        class_calloc (pwb->integral_over_k3[index_l3][index_r][index_k1], k2_size, sizeof(double), pbi->error_message);
        
        /* Increase memory counter */
        pwb->count_allocated_for_integral_over_k3 += k2_size;
      } // end of for(index_k1)
    } // end of for(index_r)
  } // end of for(index_l3)
    
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * allocated ~ %.3g MB (%ld doubles) for the k3-integral array (k_size=%d)\n",
    pwb->count_allocated_for_integral_over_k3*sizeof(double)/1e6,
    pwb->count_allocated_for_integral_over_k3, pwb->k_smooth_size);
  



  // ===================================================================================================
  // =                               Compute the INT_l3(r, k1, k2)  integral                           =
  // ===================================================================================================
  
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k3 = 0;

  abort = _FALSE_;
  #pragma omp parallel              \
    shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
  
    #pragma omp for schedule (dynamic)
    for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {

      double k1 = pwb->k_smooth_grid[index_k1];
      double pk_1 = pbi->pk_pt[index_k1];
      
      printf_log_if (pbi->bispectra_verbose, 2, 
        "     * computing the k3 integral for k1=%g, index_k1=%d\n",
        pwb->k_smooth_grid[index_k1], index_k1);
  
      /* We only need to consider those k2's that are equal to or larger than k1,
      as the shape function is assumed to be symmetric woth respect to k1<->k2<->k3 */      
      for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
    
        double k2 = pwb->k_smooth_grid[index_k2];
        double pk_2 = pbi->pk_pt[index_k2];
    
        /* Get the size of the integration grid. Note that when extrapolation is turned on, the k3-grid will also
        include values that do not satisfty the triangular condition k1 + k2 = k3. */
        int k3_size = pwb->k3_grid_size[index_k1][index_k2];
        int index_k3_lower = pwb->index_k3_lower[index_k1][index_k2];
        int index_k3_upper = pwb->index_k3_upper[index_k1][index_k2];
  
        /* Determine the integration grid. This is given by the portion of the transfer function grid
        starting at 'index_k3_lower' and ending at 'index_k3_upper' */
        for (int index_k3=index_k3_lower; index_k3<=index_k3_upper; ++index_k3)
          pwb->k3_grid[thread][index_k3-index_k3_lower] = k_tr[index_k3];

        /* If there are no points in k_tr that satisfy the triangular condition for the current
        (k1,k2), then there is no contribution to the integral */
        class_test_parallel (k3_size<=0, pbi->error_message,
          "include when triangular condition does not fit with ptr->q");
          
        /* Determine the measure for the trapezoidal rule for k3 */  
        pwb->delta_k3[thread][0] = pwb->k3_grid[thread][1] - pwb->k3_grid[thread][0];
  
        for (int index_k3=1; index_k3<(k3_size-1); ++index_k3)
          pwb->delta_k3[thread][index_k3] = pwb->k3_grid[thread][index_k3 + 1] - pwb->k3_grid[thread][index_k3 - 1];
  
        pwb->delta_k3[thread][k3_size-1] = pwb->k3_grid[thread][k3_size - 1] - pwb->k3_grid[thread][k3_size - 2];
        
        /* Shape function for this (k1,k2) slice */
        for (int index_k3=index_k3_lower; index_k3 <= index_k3_upper; ++index_k3) {
          
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
        for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

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
    
          for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
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
  
  printf_log_if (pbi->bispectra_verbose, 2, 
    " -> memorised ~ %.3g MB (%ld doubles) for the k3-integral array\n",
    pwb->count_memorised_for_integral_over_k3*sizeof(double)/1e6,
    pwb->count_memorised_for_integral_over_k3);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k3 != pwb->count_allocated_for_integral_over_k3,
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

  /* Integration grid */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
    
  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;
  
  // ======================================================================================================
  // =                                   Allocate memory for INT_l3_l2(r,k1)                              =
  // ======================================================================================================
  
  /* Because we shall recycle the array for more than one fields (T,E...) we make sure we allocate it
  only once. This is achieved by performing the allocation only at the beginning of the loop over
  the field (that is, when pwb->Y==0) */

  if (pwb->Y == 0) {
  
    /* Initialize counter */
    pwb->count_allocated_for_integral_over_k2 = 0;
  
    /* Allocate l3-level */
    class_alloc (pwb->integral_over_k2, pbi->l_size*sizeof(double ***), pbi->error_message);
  
    for (int index_l3=0; index_l3<pbi->l_size; ++index_l3) {
  
      /* Allocate l2-level. We only need l2<=l3 because of the k2<->k3 symmetry of the shape function */
      class_alloc (pwb->integral_over_k2[index_l3], (index_l3+1)*sizeof(double **), pbi->error_message);
  
      for (int index_l2=0; index_l2<=index_l3; ++index_l2) {
    
        /* Allocate r-level */
        class_alloc (pwb->integral_over_k2[index_l3][index_l2], pwb->r_size*sizeof(double *), pbi->error_message);
  
        /* Allocate 'k1' level */
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {
  
          int k1_size = pwb->k_smooth_size;
          class_alloc (pwb->integral_over_k2[index_l3][index_l2][index_r], k1_size*sizeof(double), pbi->error_message);
  
          /* Increase memory counter */
          pwb->count_allocated_for_integral_over_k2 += k1_size;
  
        } // end of for(index_k1)
      } // end of for(index_r)
    } // end of for(index_l2)
    
    printf_log_if (pbi->bispectra_verbose, 2, 
      "     * allocated ~ %.3g MB (%ld doubles) for the k2-integral array\n",
      pwb->count_allocated_for_integral_over_k2*sizeof(double)/1e6,
      pwb->count_allocated_for_integral_over_k2);
  
  } // end of if (pwb->Y==0)


  // ==============================================================================================================
  // =                                   Compute the INT_l3_l2(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k2 = 0;
    
  abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      printf_log_if (pbi->bispectra_verbose, 2, 
        "     * computing the k2 integral for r=%g, index_r=%d\n",
        pwb->r[index_r], index_r);
  
      for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {
          
        for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
      
          /* Interpolate the integral I_l2(k1,k2,r) that we computed above in the integration grid of k2.
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
                      index_l3,
                      pwb->integral_splines[thread],
                      pwb->interpolated_integral[thread],
                      pwb->f[thread],
                      pwb),
            pbi->error_message,
            pbi->error_message);
      
          for (int index_l2 = 0; index_l2 <= index_l3; ++index_l2) {  
  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k2) * l_size + index_l2) * k_tr_size]);

            /* Some debug - print transfer function */
            // if (index_r==0)
            //   for (index_k3=0; index_k3 < k_tr_size; ++index_k3)
            //     if ((index_k1==3) && (index_l3==2) && (index_l2==0))
            //       fprintf(stderr, "%g %g\n", k_tr[index_k3], transfer[index_k3]);

            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l2,
                          pwb->r[index_r],
                          &(pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k2;
  
            #pragma omp flush(abort)
  
          } // end of for(index_l2)
        } // end of for(index_l3)
      } // end of for(index_k1)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  
  printf_log_if (pbi->bispectra_verbose, 2, 
    " -> memorised ~ %.3g MB (%ld doubles) for the k2-integral array\n",
    pwb->count_memorised_for_integral_over_k2*sizeof(double)/1e6,
    pwb->count_memorised_for_integral_over_k2);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k2 != pwb->count_allocated_for_integral_over_k2,
    pbi->error_message,
    "there is a mismatch between allocated (%ld) and used (%ld) space!",
    pwb->count_allocated_for_integral_over_k2, pwb->count_memorised_for_integral_over_k2);
  
  /* Free the memory that was allocated for the I_l_k3 integral, but only if we have already computed it
  for all the required probes */
  if (pwb->Y == (pbi->bf_size-1)) {
    for (int index_l2=0; index_l2 < pbi->l_size; ++index_l2) {
      for (int index_r=0; index_r < pwb->r_size; ++index_r) {
        for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
          free (pwb->integral_over_k3[index_l2][index_r][index_k1]);
        } // end of for(index_k1)
        free (pwb->integral_over_k3[index_l2][index_r]);
      } // end of for(index_r)
      free (pwb->integral_over_k3[index_l2]);
    } // end of for(index_l2)    
    free (pwb->integral_over_k3);
  }
    
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

    f[index_k2] = (index_k1 > index_k2 ?
      pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]:
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
  
  /* Integration grid */
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

  // ==========================================================================================================
  // =                                    Allocate memory for INT_l3_l2_l1(r)                                 =
  // ==========================================================================================================
  
  /* We shall allocate the array pwb->integral_over_k1[index_l_k2][index_l_k3][index_l_k1-index_l_k1_min][index_r]
  so that the l_k1 level is the one satisfying the triangular inequality (|l_k2-l_k3| <= l_k1 <= l_k2+l_k3). 
  pwb->integral_over_k1 is recycled by the different iterations in Z-field loop, hence we allocate
  it only at the first iteration (pwb->Z==0) */
  
  if (pwb->Z==0) {

    /* The integral over k1 yields a function of (l3,l2,l1) that is symmetric under permutations of (l3,l2,l1).
    Hence, we compute and store it only for l3>=l2>=l1 configurations (that satisfy the triangular condition) */
    pwb->count_allocated_for_integral_over_k1 = pwb->r_size * pbi->n_independent_configurations;
  
    /* Allocate (l3,l2,l1)-level */
    class_alloc (pwb->integral_over_k1, pbi->n_independent_configurations*sizeof(double *), pbi->error_message);
  
    for (long int index_l3_l2_l1=0; index_l3_l2_l1 < pbi->n_independent_configurations; ++index_l3_l2_l1)
      class_alloc (pwb->integral_over_k1[index_l3_l2_l1], pwb->r_size*sizeof(double), pbi->error_message);
    
    printf_log_if (pbi->bispectra_verbose, 2, 
      "     * allocated ~ %.3g MB (%ld doubles) for the k1-integral array\n",
      pwb->count_allocated_for_integral_over_k1*sizeof(double)/1e6,
      pwb->count_allocated_for_integral_over_k1);
  }
  
  // ==========================================================================================================
  // =                                    Compute  INT_l3_l2_l1(r)  integral                                  =
  // ==========================================================================================================

  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k1 = 0;
  
  /* As for the other integrals, we parallelize the loop over 'r'. */
  abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort) private (thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
    #pragma omp for schedule (dynamic)
    for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
  
      printf_log_if (pbi->bispectra_verbose, 2, 
        "     * computing the k1 integral for r=%g, index_r=%d\n",
        pwb->r[index_r], index_r);

      for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {
    
        for (int index_l2 = 0; index_l2 <= index_l3; ++index_l2) {
  
          /* Skip those configurations that are forbidden by the triangular condition (optimization) */
          if (pbi->l[index_l2] < pbi->l[index_l3]/2)
            continue;

          /* Interpolate the integral I_l3_l2(k1,r) that we computed above in the integration grid of k1 */
          class_call_parallel (bispectra_non_separable_interpolate_over_k1 (
                        ppr,
                        ppt,
                        pbs,
                        ptr,
                        ppm,
                        pbi,
                        index_r,
                        index_l3,
                        index_l2,
                        pwb->integral_splines[thread],
                        pwb->interpolated_integral[thread],
                        pwb->f[thread],
                        pwb),
            pbi->error_message,
            pbi->error_message);      
  
          /* Determine the limits for l1, which come from the triangular inequality |l3-l2| <= l1 <= l3+l2 */
          int index_l1_min = pbi->index_l_triangular_min[index_l3][index_l2];
          int index_l1_max = MIN (index_l2, pbi->index_l_triangular_max[index_l3][index_l2]);
  
          for (int index_l1=index_l1_min; index_l1<=index_l1_max; ++index_l1) {  

            /* Index of the current (l3,l2,l1) configuration */
            long int index_l3_l2_l1 = pbi->index_l1_l2_l3[index_l3][index_l3-index_l2][index_l1_max-index_l1];
  
            /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
            int tt_size = ptr->tt_size[ppt->index_md_scalars];
            int l_size = ptr->l_size[ppt->index_md_scalars];
  
            double * transfer = &(ptr->transfer
              [ppt->index_md_scalars]
              [((ppt->index_ic_ad * tt_size + index_tt_k1) * l_size + index_l1) * k_tr_size]);

            class_call_parallel (bessel_convolution (
                          ppr,
                          pbs,
                          k_tr,
                          pbi->delta_k,
                          k_tr_size,
                          transfer,
                          pwb->interpolated_integral[thread],
                          index_l1,
                          pwb->r[index_r],
                          &(pwb->integral_over_k1[index_l3_l2_l1][index_r]),
                          pbi->error_message
                          ),
              pbi->error_message,
              pbi->error_message);
  
            /* Update the counter */
            #pragma omp atomic
            ++pwb->count_memorised_for_integral_over_k1;
              
            #pragma omp flush(abort)
  
          } // end of for(index_l1)
        } // end of for(index_l2)
      } // end of for(index_l3)
    } // end of for(index_r)
  } if (abort == _TRUE_) return _FAILURE_;  // end of parallel region
  
  printf_log_if (pbi->bispectra_verbose, 2, 
    " -> memorised ~ %.3g MB (%ld doubles) for the k1-integral array\n",
    pwb->count_memorised_for_integral_over_k1*sizeof(double)/1e6,
    pwb->count_memorised_for_integral_over_k1);
  
  /* Check that we correctly filled the array */
  class_test (pwb->count_memorised_for_integral_over_k1 != pwb->count_allocated_for_integral_over_k1,
              pbi->error_message,
              "there is a mismatch between allocated (%ld) and used (%ld) space!",
                pwb->count_allocated_for_integral_over_k1, pwb->count_memorised_for_integral_over_k1);
  
  /* Free the memory that was allocated for the integral over k2, but do that only when we are at
  the last iteration of the Z and Y loops. */
  if ((pwb->Y == (pbi->bf_size-1)) && (pwb->Z == (pbi->bf_size-1))) {
    for (int index_l2=0; index_l2<pbi->l_size; ++index_l2) {
      for (int index_l1=0; index_l1<=index_l2; ++index_l1) {
        for (int index_r=0; index_r < pwb->r_size; ++index_r) {      
          free (pwb->integral_over_k2[index_l2][index_l1][index_r]);
        } // end of for(index_r)
        free (pwb->integral_over_k2[index_l2][index_l1]);
      } // end of for(index_l1)
      free (pwb->integral_over_k2[index_l2]);
    } // end of for(index_l2)
    free (pwb->integral_over_k2);
  }
  
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
    int index_l3,
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
    
    f[index_k1] = pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1];
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
      
    /* Interpolate for each value of l3, l2, r */
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
  // if ((index_l3==48) && (index_l2==35) && (index_r==50)) {
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
  
  /* We parallelize the outer loop over 'l1'. */
  int abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,pbi,pwb,abort)
  {
  
    #pragma omp for schedule (dynamic)
    for (int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

      printf_log_if (pbi->bispectra_verbose, 2, 
        "     * computing the r-integral for l1=%d, index_l1=%d\n",
        pbi->l[index_l1], index_l1);
    
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
  
  /* Free the memory that was allocated for the integral over k1, but do that only when we are at
  the last iteration of the Z and Y . */
  if (pwb->Z == (pbi->bf_size-1)) {
    for (long int index_l1_l2_l3 = 0; index_l1_l2_l3 < pbi->n_independent_configurations; ++index_l1_l2_l3)
      free (pwb->integral_over_k1[index_l1_l2_l3]);
    free (pwb->integral_over_k1);
  }
    
  return _SUCCESS_;
  
}



/**
 * Free the memory associated to the workspace used for the computation of the
 * non-separable bispectra.
 */
int bispectra_non_separable_workspace_free (
    struct bispectra * pbi,
    struct bispectra_workspace_non_separable * pwb
    )
{
  
  free (pwb->r);
  free (pwb->delta_r);  
  free (pwb->k_smooth_grid);

  for (int index_k1=0; index_k1 < pwb->k_smooth_size; ++index_k1) {
    free (pwb->index_k3_lower[index_k1]);
    free (pwb->index_k3_upper[index_k1]);
    free (pwb->k3_grid_size[index_k1]);    
  }
  free (pwb->index_k3_lower);
  free (pwb->index_k3_upper);
  free (pwb->k3_grid_size);

  /* Parallelization variables */
  int thread = 0;
  int abort = _FALSE_;

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
 * Compute the lensing correction to the CMB bispectrum.
 *
 * We use a generalised version of eq. 12 in Hanson, Smith, Challinor & Liguori 2009
 * (http://arxiv.org/abs/0905.4732) which includes polarisation (credits to Christian
 * Fidler).
 */

int bispectra_lensing (
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
    int index_bt
    )
{

  printf_log_if (pbi->bispectra_verbose, 0,
    " -> computing lensing of the %s bispectrum\n", pbi->bt_labels[index_bt]);

  /* Compute the RMS of the lensing deflection angle */

  double R = 0;
  
  for (int l=2; l<=pbi->l_max; ++l)
    R += l*(l+1)*(2*l+1)/(4*_PI_) * pbi->cls[psp->index_ct_pp][l-2];
  
  printf_log_if (pbi->bispectra_verbose, 0,
    " -> RMS of the lensing deflection angle = %g arcmin\n", sqrt(R) * 180/_PI_ * 60);
  

  /* Loop over all bispectra to compute the lensing correction */

  int abort = _FALSE_;

  #pragma omp parallel for schedule (dynamic)
  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {

    int l1 = pbi->l[index_l1];
    
    printf_log_if (pbi->bispectra_verbose, 1,
     "     * computing lensing for l1=%d\n", l1);
    
    for (int X1=0; X1 < pbi->bf_size; ++X1) {

      int F_X1 = pbi->field_spin[X1];

      for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {

        int l2 = pbi->l[index_l2];      
        int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
        int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);

        for (int X2=0; X2 < pbi->bf_size; ++X2) {

          int F_X2 = pbi->field_spin[X2];

          for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

            int l3 = pbi->l[index_l3];
            long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

            for (int X3=0; X3 < pbi->bf_size; ++X3) {

              int F_X3 = pbi->field_spin[X3];

              /* Overall R factor */

              double R_factor = - 0.5 * R * (
                (l1+F_X1)*(l1-F_X1+1) + (l1-F_X1)*(l1+F_X1+1) +
                (l2+F_X2)*(l2-F_X2+1) + (l2-F_X2)*(l2+F_X2+1) +
                (l3+F_X3)*(l3-F_X3+1) + (l3-F_X3)*(l3+F_X3+1)
              );
                            
              /* Convolution factor */
              
              double convolution_factor_123 = 0;
              double convolution_factor_231 = 0;
              double convolution_factor_312 = 0;
              
              class_call_parallel (bispectra_lensing_convolution (ppr, ptr, psp, ple, pbi,
                                     index_l1, index_l2, index_l3,
                                     X1, X2, X3,
                                     index_bt,
                                     &convolution_factor_123),
                pbi->error_message,
                pbi->error_message);

              /* In presence of multiple fields (ie. T and E) permute them */

              if (pbi->bf_size > 1) {

                class_call_parallel (bispectra_lensing_convolution (ppr, ptr, psp, ple, pbi,
                                       index_l2, index_l3, index_l1,
                                       X2, X3, X1,
                                       index_bt,
                                       &convolution_factor_231),
                  pbi->error_message,
                  pbi->error_message);

                class_call_parallel (bispectra_lensing_convolution (ppr, ptr, psp, ple, pbi,
                                       index_l3, index_l1, index_l2,
                                       X3, X1, X2,
                                       index_bt,
                                       &convolution_factor_312),
                  pbi->error_message,
                  pbi->error_message);
                  
              }

              /* With one field, the permutation would yield the same result */

              else {
                
                double convolution_factor_231 = convolution_factor_123;
                double convolution_factor_312 = convolution_factor_123;
                
              }


              /* Final lensing correction */

              double b = pbi->bispectra[index_bt][X1][X2][X3][index_l1_l2_l3];

              pbi->lensing_correction[index_bt][X1][X2][X3][index_l1_l2_l3] =
                b * R_factor/4
                + convolution_factor_123
                + convolution_factor_231
                + convolution_factor_312;
              
            }
          }
        }
      }
    }
  }
  
  if (abort == _TRUE_)
    return _FAILURE_;
  

  /* Add the lensing correction to the bispectrum */
  
  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {
    for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
        long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];
        for (int X1=0; X1 < pbi->bf_size; ++X1)
          for (int X2=0; X2 < pbi->bf_size; ++X2)
            for (int X3=0; X3 < pbi->bf_size; ++X3)
              pbi->bispectra[index_bt][X1][X2][X3][index_l1_l2_l3] += pbi->lensing_correction[index_bt][X1][X2][X3][index_l1_l2_l3];
      }
    }
  }



  return _SUCCESS_;
  
}



/**
 * Compute the convolution part of the lensing correction to the bispectrum
 * configuration b_l1_l2_l3.
 *
 * The convolution term is computed by solving a sum over the three dummy
 * multipoles (l,p,q). The lensing potential C_l^PP is computed only in l,
 * while the bispectrum is computed in (l1,p,q). The (p,q) multipoles are
 * symmetric.
 */

int bispectra_lensing_convolution (
     struct precision * ppr,
     struct transfers * ptr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int index_l1, int index_l2, int index_l3,
     int X1, int X2, int X3,
     int index_bt,
     double * result
     )
{
  
  int l1 = pbi->l[index_l1], F_X1 = pbi->field_spin[X1];
  int l2 = pbi->l[index_l2], F_X2 = pbi->field_spin[X2];
  int l3 = pbi->l[index_l3], F_X3 = pbi->field_spin[X3];

  class_test (!is_triangular_int(l1,l2,l3),
    pbi->error_message,
    "l1=%d, l2=%d and l3=%d do not form a triangle", l1, l2, l3);

  /* The convolution sum goes over (l,p,q). We take all values in l between 2 and
  l_max, while we only take the values of p and q in pbi->l and later interpolate */
  // double correction[pbi->l_size][pbi->l_size];
  // for (int index_p=0; index_p < pbi->l_size; ++index_p)
  //   for (int index_q=0; index_q < pbi->l_size; ++index_q)
  //     correction[index_p][index_q] = 0;

  /* Since the bispectrum does not depend on l, we store it in a (p,q) array
  and reuse it when l changes */
  double bispectrum[pbi->l_size][pbi->l_size];
  short filled[pbi->l_size][pbi->l_size];
  for (int i=0; i < pbi->l_size; ++i)
    for (int j=0; j < pbi->l_size; ++j)
      filled[i][j] = _FALSE_;

  /* Count how many (l,p,q) configurations we shall consider */
  long int count = 0;

  /* Initialise output */
  *result = 0;

  // -------------------------------------------------------------------------------
  // -                                  Loop on l                                  -
  // -------------------------------------------------------------------------------

  for (int l=2; l <= pbi->l_max; ++l) {

    double C_PP = pbi->cls[psp->index_ct_pp][l-2];

    /* The lensing potential enters as an overall factor; if it vanishes, the whole
    term vanishes. This happens whenever pbi->lmax_lensing_corrT for temperature or
    pbi->lmax_lensing_corrE for polarisation are smaller than pbi->l_max. */
    if (C_PP == 0)
      continue;

    /* Range of p dictated by the triangular inequality on l,l2,p */
    int p_min = MAX (abs(l-l2), 2);
    int p_max = MIN (l+l2, pbi->l_max);

    int index_p_min = 0;
    while (pbi->l[index_p_min] < p_min)
      ++index_p_min;

    int index_p_max = pbi->l_size-1;
    while (pbi->l[index_p_max] > p_max)
      --index_p_max;    

    int p_size = index_p_max - index_p_min + 1;

    /* Weights for the trapezoidal sum along the p direction */

    double * delta_p;

    if (p_size > 0) {

      class_alloc (delta_p, p_size*sizeof(double), pbi->error_message);
      
      if (p_size == 1) {

        delta_p[0] = p_max - p_min;

      }

      if (p_size > 1) {
        
        /* Build integration grid */
        
        double p[p_size];
        
        p[0] = p_min;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          p[index_p] = pbi->l[index_p_min+index_p];

        p[p_size-1] = p_max;

        /* Build integration measure */

        delta_p[0] = (p[1] - p[0] + 1)/2.0;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          delta_p[index_p] = (p[index_p+1] - p[index_p-1])/2.0;

        delta_p[p_size-1] = (p[p_size-1] - p[p_size-2] + 1)/2.0;

      }
      
    } // if(p_size>0)


    /* Compute 3j symbol with p */

    double threej_p[2*pbi->l_max+1];
    double min_D, max_D;
    class_call (drc3jj (
                  l, l2, 0, -F_X2,
                  &min_D, &max_D,
                  threej_p,
                  (2*pbi->l_max+1),
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);
    int p_min_3j = (int)(min_D + _EPS_);
    int p_max_3j = (int)(max_D + _EPS_);
    
    for (int p=p_min; p <= p_max; ++p)
      threej_p[p-p_min_3j] *= sqrt((2*l+1)*(2*p+1)*(2*l2+1)/(4*_PI_)) * 0.5 * (l*(l+1)+p*(p+1)-l2*(l2+1));


    /* Range of q dictated by the triangular inequality on l,q,l3 */
    int q_min = MAX (abs(l-l3), 2);
    int q_max = MIN (l+l3, pbi->l_max);


    // -------------------------------------------------------------------------------
    // -                                  Loop on p                                  -
    // -------------------------------------------------------------------------------

    for (int index_p=index_p_min; index_p <= index_p_max; ++index_p) {

      int p = pbi->l[index_p];
      
      /* Range of q dictated by the triangular inequality on l1,q,p */
      q_min = MAX (q_min, abs(l1-p));
      q_max = MIN (q_max, l1+p);
      
      int index_q_min = pbi->index_l_triangular_min[index_l1][index_p];
      int index_q_max = pbi->index_l_triangular_max[index_l1][index_p];
      int q_size = index_q_max - index_q_min + 1;

      /* Weights for the trapezoidal sum along the q direction */

      double * delta_q;

      if (q_size > 0) {

        class_alloc (delta_q, q_size*sizeof(double), pbi->error_message);
      
        if (q_size == 1) {

          delta_q[0] = q_max - q_min;

        }

        if (q_size > 1) {
        
          /* Build integration grid */
        
          double q[q_size];
        
          q[0] = q_min;

          for (int index_q=1; index_q < q_size-1; ++index_q)
            q[index_q] = pbi->l[index_q_min+index_q];

          q[q_size-1] = q_max;

          /* Build integration measure */

          delta_q[0] = (q[1] - q[0] + 1)/2.0;

          for (int index_q=1; index_q < q_size-1; ++index_q)
            delta_q[index_q] = (q[index_q+1] - q[index_q-1])/2.0;

          delta_q[q_size-1] = (q[q_size-1] - q[q_size-2] + 1)/2.0;

        }

      } // if(q_size>0)


      /* Comqute 3j symbol with q */

      double threej_q[2*pbi->l_max+1];
      double min_D, max_D;
      class_call (drc3jj (
                    l, l3, 0, -F_X3,
                    &min_D, &max_D,
                    threej_q,
                    (2*pbi->l_max+1),
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      int q_min_3j = (int)(min_D + _EPS_);
      int q_max_3j = (int)(max_D + _EPS_);
    
      for (int q=q_min; q <= q_max; ++q)
        threej_q[q-q_min_3j] *= sqrt((2*l+1)*(2*q+1)*(2*l3+1)/(4*_PI_)) * 0.5 * (l*(l+1)+q*(q+1)-l3*(l3+1));


      /* Compute the 6j symbol for all values of q */
      double sixj_q[2*pbi->l_max+1];
      class_call (drc6j (
                    /*q,*/ p, l1, l2, l3, l,
                    &min_D, &max_D,
                    sixj_q,
                    (2*pbi->l_max+1),
                    pbi->error_message       
                    ),
        pbi->error_message,
        pbi->error_message);
      int q_min_6j = (int)(min_D + _EPS_);
      int q_max_6j = (int)(max_D + _EPS_);


      /* Computing alternatign sign */
      int alternating_sign = ALTERNATING_SIGN (l1+l2+p);


      // -------------------------------------------------------------------------------
      // -                                  Loop on q                                  -
      // -------------------------------------------------------------------------------

      for (int index_q=index_q_min; index_q <= index_q_max; ++index_q) {

        int q = pbi->l[index_q];

        if ((q < q_min) || (q > q_max))
          continue;

        int X2_ = X2;
        int X3_ = X3;

        /* Interpolate the bispectrum in (l1,p,q). We do so only at the first
        iteration of l, because the bispectrum does not depend on it. */

        // if (count_l == 0) {
        //
        //   int X2_ = X2;
        //   int X3_ = X3;

          /* Extract closest node to the left of q  */

          // int index_q_left = ptr->index_l_left[q];
          // int q_left = pbi->l[index_q_left];
          // double b_left = 0;
          //
          // printf ("q = %d\n", q);
          // printf ("q_left = %d\n", q_left);

        if (filled[index_p][index_q] == _FALSE_) {

          class_call (bispectra_at_l1l2l3 (
                        pbi,
                        index_bt,
                        index_l1, index_p, index_q,
                        X1, X2_, X3_,
                        &bispectrum[index_p][index_q],
                        NULL),
            pbi->error_message,
            pbi->error_message);

          filled[index_p][index_q] = _TRUE_;

        }


          // /* If q belongs to our l-sampling, there is no need for interpolation */
          //
          // if (q == q_left) {
          //
          //   bispectrum[index_p][q] = b_left;
          //
          // }
          //
          // else if (!is_triangular_int(l1,p,q)) {}
          //
          // /* Extract closest node to the right of q */
          //
          // else {
          //
          //   int q_right = pbi->l[index_q_left+1];
          //   double b_right = 0;
          //
          //   printf ("q_right = %d\n", q_right);
          //
          //   class_call (bispectra_at_l1l2l3 (
          //                 pbi,
          //                 index_bt,
          //                 index_l1, index_p, index_q_left+1,
          //                 X1, X2_, X3_,
          //                 &b_right,
          //                 NULL),
          //     pbi->error_message,
          //     pbi->error_message);
          //
          //
          //   /* Perform interpolation */
          //
          //   double a = (q_right-q)/(double)(q_right-q_left);
          //   bispectrum[index_p][q] = a*b_left + (1-a)*b_right;
          //
          // }

        *result += C_PP *
                   alternating_sign * threej_p[p-p_min_3j] * threej_q[q-q_min_3j] * sixj_q[q-q_min_6j] *
                   delta_p[index_p-index_p_min] * delta_q[index_q-index_q_min] *
                   bispectrum[index_p][index_q];

        count++;

      } // for index_q

      if (q_size > 0)
        free (delta_q);

    } // for index_p

    if (p_size > 0)
      free (delta_p);

  } // for l
  
  /* Debug - Print how many bispectrum configurations we computed in this (l,p,q) loop */
  if ((index_l1 == (pbi->l_size-1)) && (index_l2 == index_l1) && (index_l3 == index_l2))
    printf ("%s_%s: correction[l3=%4d] = %g\n",
      pbi->bt_labels[index_bt], pbi->bfff_labels[X1][X2][X3], l3, *result);
  // if ((index_l1 == (pbi->l_size-1)) && (index_l2 == index_l1))
  //   printf ("count(l1=%d,l2=%d,l3=%4d) = %ld\n", l1, l2, l3, count);
  
  return _SUCCESS_;
  
}


// int bispectra_lensing_convolution (
//      struct precision * ppr,
//      struct spectra * psp,
//      struct lensing * ple,
//      struct bispectra * pbi,
//      int l1, int l2, int l3,
//      int X1, int X2, int X3,
//      int index_bt,
//      double * result
//      )
// {
//
//   class_test (!is_triangular_int(l1,l2,l3),
//     pbi->error_message,
//     "l1=%d, l2=%d and l3=%d do not form a triangle", l1, l2, l3);
//
//   /* Sice the bispectrum does not depend on l, we store it in a (p,q) array */
//   int count = 0;
//   double correction[pbi->full_l_size+2][pbi->full_l_size+2];
//
//   for (int l=2; l <= pbi->l_max; ++l) {
//
//     double C_PP = pbi->cls[psp->index_ct_pp][l-2];
//
//     /* The lensing potential enters as an overall factor; if it vanishes, the whole
//     term vanishes. This happens whenever pbi->lmax_lensing_corrT for temperature or
//     pbi->lmax_lensing_corrE for polarisation are smaller than pbi->l_max. */
//     if (C_PP == 0)
//       continue;
//
//     /* Enforce the triangular inequality on l,l2,p */
//     int p_min = MAX (abs(l-l2), 2);
//     int p_max = MIN (l+l2, pbi->l_max);
//
//     /* Enforce the triangular inequality on l,q,l3 */
//     int q_min = MAX (abs(l-l3), 2);
//     int q_max = MIN (l+l3, pbi->l_max);
//
//     for (int p=p_min; p <= p_max; ++p) {
//
//       /* Enforce the triangular inequality on l1,q,p */
//       int q_min = MAX (q_min, abs(l1-p));
//       int q_max = MIN (q_max, l1+p);
//
//       for (int q=q_min; q <= q_max; ++q) {
//
//         if (count == 0)
//           correction[p][q] = 1234;
//
//
//
//       } // for p
//
//     } // for q
//
//     count++;
//
//   } // for l
//
//
//   // int index_l[4] = {0, index_l1, index_l2, index_l3};
//   // int order[4];
//   //
//   // class_call (ordering_int (index_l, order, pbi->error_message),
//   //   pbi->error_message, pbi->error_message);
//   //
//   // int index_1 = index_l[order[1]]; /* Smallest l */
//   // int index_2 = index_l[order[2]]; /* Mid l */
//   // int index_3 = index_l[order[3]]; /* Largest l */
//   //
//   // /* Index of the current (l1,l2,l3) configuration */
//   // int index_1_max = MIN (index_2, pbi->index_l_triangular_max[index_3][index_2]);
//   // long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_3][index_3-index_2][index_1_max-index_1];
//   //
//   // /* Extract the bispectrum in (l1,l2,l3), using the fact that a permutation of
//   // (l1,l2,l3) is cancelled by the same permutation of (X,Y,Z). For example:
//   // b^TTE(l1,l2,l3) = b^TET(l1,l3,l2). */
//   //
//   // int XYZ[4] = {0, X, Y, Z};
//   //
//   // for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
//   //
//   //   bispectrum[index_bt] = pbi->bispectra[index_bt][XYZ[order[3]]][XYZ[order[2]]][XYZ[order[1]]][index_l1_l2_l3];
//
//
//   return _SUCCESS_;
//
// }
//


/**
 * Prepare the bispectrum to be interpolated in any (l1,l2,l3) configuration.
 */

// int bispectra_interpolation (
//     struct precision * ppr,
//     struct background * pba,
//     struct thermo * pth,
//     struct perturbs * ppt,
//     struct bessels * pbs,
//     struct transfers * ptr,
//     struct primordial * ppm,
//     struct spectra * psp,
//     struct lensing * ple,
//     struct bispectra * pbi
//     )
// {
//
//   // ====================================================================================
//   // =                               Interpolation arrays                               =
//   // ====================================================================================
//
//   /* By default we do never interpolate the analytic bispectra, because they can be
//   computed at any time without spending too much computational time. For testing
//   purposes, however, one might want to interpolate the analytical bispectra as well. */
//   if (pbi->always_interpolate_bispectra == _TRUE_) {
//     for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
//       pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
//     pbi->n[analytical_bispectrum] = 0;
//   }
//
//   /* Determine which is the first bispectrum for which we should compute the interpolation
//   mesh. This is equivalent to the first non-analytical bispectrum, because the analytical
//   bispectra are not interpolated at all */
//   pbi->first_non_analytical_index_ft = 0;
//   while (pbi->bispectrum_type[pbi->index_bt_of_ft[pbi->first_non_analytical_index_ft]] == analytical_bispectrum
//   && (pbi->first_non_analytical_index_ft<pbi->ft_size))
//     pbi->first_non_analytical_index_ft++;
//
//   /* No need to compute meshes if all bispectra are analytical, as interpolation is not
//   needed in this case */
//   pbi->has_only_analytical_bispectra = _FALSE_;
//   if (pbi->first_non_analytical_index_ft == pbi->ft_size)
//     pbi->has_only_analytical_bispectra = _TRUE_;
//
//
//   if (pbi->interpolation_method == mesh_interpolation) && (pbi->has_only_analytical_bispectra == _FALSE_)) {
//
//     // -------------------------------------------------------------------------------
//     // -                        Mesh interpolation parameters                        -
//     // -------------------------------------------------------------------------------
//
//     /* We set the linking length of the fine mesh to roughly the smallest distance
//     in the grid */
//     pbi->link_lengths[0] = 2 * (pbi->l[1] - pbi->l[0]);
//     pbi->group_lengths[0] = 0.1 * (pbi->l[1] - pbi->l[0]);
//     pbi->soft_coeffs[0] = 0.5;
//
//     /* We set the linking length of the coarse mesh to roughly the largest distance
//     in the grid */
//     int linstep = ppr->l_linstep;
//
//     /* If ppr->l_linstep=1, then all l are used and there is effectively no interpolation.
//     However, if also ppr->compute_only_even_ls == _TRUE_, then the actual step between one
//     multipole and the other is doubled, as the odd l are skipped. */
//     if ((l_linstep==1) && ((ppr->compute_only_even_ls==_TRUE_) || (ppr->compute_only_odd_ls==_TRUE_)))
//       linstep = 2:
//
//     pbi->link_lengths[1] = 0.5/sqrt(2) * l_linstep;
//     pbi->group_lengths[1] = 0.1 * l_linstep;
//     pbi->soft_coeffs[1] = 0.5;
//
//     for (int index_mesh=0; index_mesh < pbi->n_meshes; ++index_mesh)
//       class_test (pbi->link_lengths[index_mesh] <= pbi->group_lengths[index_mesh],
//         pbi->error_message,
//         "the linking length must be larger than the grouping length.");
//
//
//
//     // -------------------------------------------------------------------------------
//     // -                            Determine turnover point                         -
//     // -------------------------------------------------------------------------------
//
//     /* Never use the fine grid if the large step is used since the beginning */
//     if ((pbi->l[1]-pbi->l[0]) >= l_linstep) {
//       pbi->l_turnover[0] = pbi->l[0];
//     }
//
//     /* Never use the coarse grid if the large step is never used. The MAX() is needed
//     because the last point is fixed and does not depend on the grid spacing. */
//     else if (MAX(pbi->l[pbi->l_size-1]-pbi->l[pbi->l_size-2], pbi->l[pbi->l_size-2]-pbi->l[pbi->l_size-3]) < l_linstep) {
//       pbi->l_turnover[0] = pbi->l[pbi->l_size-1] + 1;
//     }
//
//     /* Otherwise, find the turnover multipole as the point in the grid where we
//     start using the large step */
//     else {
//       int index_l = 0;
//       while ((index_l < pbi->l_size-1) && ((pbi->l[index_l+1] - pbi->l[index_l]) < l_linstep))
//         index_l++;
//       pbi->l_turnover[0] = pbi->l[index_l-1];
//     }
//
//     printf_log_if (pbi->fisher_verbose, 1,
//       "     * mesh_interpolation: l_turnover=%d, n_boxes=[%d,%d], linking lengths=[%g,%g], grouping lengths=[%g,%g]\n",
//       pbi->l_turnover[0],
//       (int)ceil(pbi->l_turnover[0]/ (pbi->link_lengths[0]*(1.+pbi->soft_coeffs[0]))),
//       (int)ceil(pbi->l[pbi->l_size-1] / (pbi->link_lengths[1]*(1.+pbi->soft_coeffs[1]))),
//       pbi->link_lengths[0], pbi->link_lengths[1],
//       pbi->group_lengths[0], pbi->group_lengths[1]);
//
//   } // if(mesh_interpolation)
//
//
//
//   // ====================================================================================
//   // =                               Assign window functions                            =
//   // ====================================================================================
//
//   /* Assign to each bispectrum a window function suitable for its interpolation. See
//   header file for mode details. Comment out a bispectrum to have it interpolated without
//   a window function. It seems that using window functions that cross the zero might
//   generate some issues. This has yet to be verified rigorously, though. */
//   for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
//
//     pbi->window_function[index_bt] = NULL;
//
//     /* Previously, SONG used to interpolate the bispectra in the two smallest l-multipoles.
//     That was not a good idea, and as a result we needed to use window functions to obtain
//     precise results. Now we interpolate the two largest l's instead, and we don't need
//     window functions anymore. You can ignore what follows in this block. */
//     continue;
//
//     /* In absence of reionisation and for polarisation, better results are obtained without
//     a window  function. The reason is that the C_l's for polarisation are tiny at small l's
//     without reionisation, so multiplying and dividing by the window functions creates numerical
//     instabilities. The patology is stronger for bispectra that peak in the squeezed limit,
//     as in that case the large scales are those where most of the signal comes from. */
//     /* TODO: the above statement has to be corrected in terms of the new interpolation,
//     which relies on interpolating the two largest multipoles */
//     if ((pth->reio_parametrization != reio_none)
//       && (pbi->has_bispectra_e == _TRUE_) && (pbi->bf_size > 1)) {
//       continue;
//     }
//     else {
//
//       if ((pbi->has_local_model == _TRUE_) && (index_bt == pbi->index_bt_local))
//         pbi->window_function[index_bt] = bispectra_local_window_function;
//
//       else if ((pbi->has_intrinsic == _TRUE_) && (index_bt == pbi->index_bt_intrinsic))
//         pbi->window_function[index_bt] = bispectra_local_window_function;
//
//       else if ((pbi->has_intrinsic_squeezed == _TRUE_) && (index_bt == pbi->index_bt_intrinsic_squeezed))
//         pbi->window_function[index_bt] = bispectra_local_window_function;
//     }
//
//     /* For the non-squeezed bispectra, we always use the window function. */
//     if ((pbi->has_equilateral_model == _TRUE_) && (index_bt == pbi->index_bt_equilateral))
//       pbi->window_function[index_bt] = bispectra_local_window_function;
//
//     else if ((pbi->has_orthogonal_model == _TRUE_) && (index_bt == pbi->index_bt_orthogonal))
//       pbi->window_function[index_bt] = bispectra_local_window_function;
//
//     else if ((pbi->has_galileon_model==_TRUE_) && (index_bt == pbi->index_bt_galileon_gradient))
//       pbi->window_function[index_bt] = bispectra_local_window_function;
//
//     else if ((pbi->has_galileon_model==_TRUE_) && (index_bt == pbi->index_bt_galileon_time))
//       pbi->window_function[index_bt] = bispectra_local_window_function;
//
//   }
//
//   return _SUCCESS_;
//
// }




/**
 * Allocate the array of pointers 'meshes' given as an input. This points to two
 * mesh workspaces, one finely sampled and the other coarsely sampled, per each of the
 * bispectra to be interpolated.
 */  

// int fisher_allocate_interpolation_mesh(
//         struct precision * ppr,
//         struct spectra * psp,
//         struct lensing * ple,
//         struct bispectra * pbi,
//         struct fisher * pfi,
//         struct interpolation_mesh ******* meshes
//         )
// {
//
//   /* Allocate one mesh interpolation workspace per type of bispectrum */
//   class_alloc ((*meshes),
//     pfi->ft_size*sizeof(struct interpolation_mesh *****),
//     pfi->error_message);
//
//   /* Allocate worspaces and intialize counters */
//   for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
//
//     /* The analytical bispectra are never interpolated */
//     if (pbi->bispectrum_type[pfi->index_bt_of_ft[index_ft]] == analytical_bispectrum)
//       continue;
//
//     class_alloc ((*meshes)[index_ft],
//       pfi->ff_size*sizeof(struct interpolation_mesh ****),
//       pfi->error_message);
//
//       for (int X = 0; X < pfi->ff_size; ++X) {
//
//         class_alloc ((*meshes)[index_ft][X],
//           pfi->ff_size*sizeof(struct interpolation_mesh ***),
//           pfi->error_message);
//
//         for (int Y = 0; Y < pfi->ff_size; ++Y) {
//
//           class_alloc ((*meshes)[index_ft][X][Y],
//             pfi->ff_size*sizeof(struct interpolation_mesh **),
//             pfi->error_message);
//
//           for (int Z = 0; Z < pfi->ff_size; ++Z) {
//
//             class_alloc ((*meshes)[index_ft][X][Y][Z],
//               pfi->n_meshes*sizeof(struct interpolation_mesh *),
//               pfi->error_message);
//
//             for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
//
//               class_alloc ((*meshes)[index_ft][X][Y][Z][index_mesh],
//                 sizeof(struct interpolation_mesh),
//                 pfi->error_message);
//
//               /* Set the maximum l. The fine grid does not need to go up to l_max */
//               if (index_mesh == 0)
//                 (*meshes)[index_ft][X][Y][Z][index_mesh]->l_max = (double)pfi->l_turnover;
//               else
//                 (*meshes)[index_ft][X][Y][Z][index_mesh]->l_max = (double)pbi->l[pbi->l_size-1];
//
//               /* Initialise the structure with the interpolation parameters */
//               (*meshes)[index_ft][X][Y][Z][index_mesh]->link_length = pfi->link_lengths[index_mesh];
//               (*meshes)[index_ft][X][Y][Z][index_mesh]->group_length = pfi->group_lengths[index_mesh];
//               (*meshes)[index_ft][X][Y][Z][index_mesh]->soft_coeff = pfi->soft_coeffs[index_mesh];
//
//           } // for(index_mesh)
//         } // for(Z)
//       } // for(Y)
//     } // for(X)
//   } // for(index_ft)
//
//   return _SUCCESS_;
//
// }



/**
 * Deallocate the array of pointers 'meshes' given as an input.
 */

// int fisher_free_interpolation_mesh(
//         struct bispectra * pbi,
//         struct fisher * pfi,
//         struct interpolation_mesh ******* meshes
//         )
// {
//
//   for (int index_ft=(pfi->ft_size-1); index_ft >= pfi->first_non_analytical_index_ft; --index_ft) {
//
//     if (pbi->bispectrum_type[pfi->index_bt_of_ft[index_ft]] == analytical_bispectrum)
//       continue;
//
//     for (int X = (pfi->ff_size-1); X >= 0; --X) {
//       for (int Y = (pfi->ff_size-1); Y >= 0; --Y) {
//         for (int Z = (pfi->ff_size-1); Z >= 0; --Z) {
//           for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
//
//             if (((index_mesh == 0) && (pfi->l_turnover <= pbi->l[0]))
//             || ((index_mesh == 1) && (pfi->l_turnover > pbi->l[pbi->l_size-1]))) {
//               free ((*meshes)[index_ft][X][Y][Z][index_mesh]);
//               continue;
//             }
//             else {
//               if (pfi->bispectra_interpolation == mesh_interpolation) {
//                 class_call (mesh_2D_free ((*meshes)[index_ft][X][Y][Z][index_mesh]),
//                   pfi->error_message,
//                   pfi->error_message);
//               }
//               else {
//                 class_call (mesh_3D_free ((*meshes)[index_ft][X][Y][Z][index_mesh]),
//                   pfi->error_message,
//                   pfi->error_message);
//               }
//             }
//           } free ((*meshes)[index_ft][X][Y][Z]);
//         } free ((*meshes)[index_ft][X][Y]);
//       } free ((*meshes)[index_ft][X]);
//     } free ((*meshes)[index_ft]);
//   } free ((*meshes));
//
//   return _SUCCESS_;
//
// }





/**
 * Build the two-dimensional mesh for the interpolation of the bispectra.
 */

// int fisher_create_2D_interpolation_mesh(
//         struct precision * ppr,
//         struct spectra * psp,
//         struct lensing * ple,
//         struct bispectra * pbi,
//         struct fisher * pfi,
//         int index_l1,
//         struct interpolation_mesh ****** meshes
//         )
// {
//
//   printf_log_if (pfi->fisher_verbose, 2,
//     " -> preparing 2D interpolation mesh for index_l1=%d\n", index_l1);
//
//   int l1 = pbi->l[index_l1];
//
//   /* Count the number of (l2,l3) nodes in the considered l1-slice of the bispectrum */
//   long int n_points = 0;
//   for (int index_l2=index_l1; index_l2 < pbi->l_size; ++index_l2) {
//     int index_l3_min = MAX (index_l2, pbi->index_l_triangular_min[index_l1][index_l2]);
//     int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
//     for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
//       n_points ++;
//   }
//
//   printf_log_if (pfi->fisher_verbose, 3,
//     "     * the l1=%d slice has %ld point%s\n", l1, n_points, ((n_points==1)?"":"s"));
//
//   /* Allocate the array that will contain the re-arranged bispectra */
//   double ** values;
//   class_alloc (values, n_points*sizeof(double *), pfi->error_message);
//   for (int index_l2_l3=0; index_l2_l3 < n_points; ++index_l2_l3)
//     class_calloc (values[index_l2_l3], 3, sizeof(double), pfi->error_message);
//
//   for (int index_ft=0; index_ft < pfi->ft_size; ++index_ft) {
//
//     int index_bt = pfi->index_bt_of_ft[index_ft];
//
//     /* The analytical bispectra are never interpolated */
//     if (pbi->bispectrum_type[index_bt] == analytical_bispectrum)
//       continue;
//
//     for (int X = 0; X < pfi->ff_size; ++X) {
//       for (int Y = 0; Y < pfi->ff_size; ++Y) {
//         for (int Z = 0; Z < pfi->ff_size; ++Z) {
//
//           /* Corresponding field indices in the bispectrum module */
//           int X_ = pfi->index_bf_of_ff[X];
//           int Y_ = pfi->index_bf_of_ff[Y];
//           int Z_ = pfi->index_bf_of_ff[Z];
//
//           // ---------------------------------------------------------------------------
//           // -                          Rearrange the bispectra                        -
//           // ---------------------------------------------------------------------------
//
//           /* Index of the considered (l2,l3) node */
//           long int index_l2_l3 = 0;
//
//           for (int index_l2=index_l1; index_l2 < pbi->l_size; ++index_l2) {
//
//             int l2 = pbi->l[index_l2];
//
//             /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
//             int index_l3_min = MAX (index_l2, pbi->index_l_triangular_min[index_l1][index_l2]);
//             int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
//
//             for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {
//
//               int l3 = pbi->l[index_l3];
//
//               /* Determine the index of this particular index_l1_l2_l3, taking into account that now
//               l1 is the smallest multipole */
//               int index_l1_max = MIN (index_l2, pbi->index_l_triangular_max[index_l2][index_l3]);
//               long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l3][index_l3-index_l2][index_l1_max-index_l1];
//
//               /* By default, assume that no window function is needed */
//               double inverse_window = 1;
//
//               /* Compute the value of the window function for this (l1,l2,l3) */
//               if (pbi->window_function[index_bt] != NULL) {
//
//                 double dump;
//
//                 class_call ((*pbi->window_function[index_bt]) (
//                               ppr, psp, ple, pbi,
//                               l3, l2, l1, /* smallest one goes in third position  */
//                               X_, Y_, Z_,
//                               dump,
//                               dump,
//                               dump,
//                               &inverse_window),
//                   pbi->error_message,
//                   pfi->error_message);
//
//               }
//
//               /* The 1st element of 'values' is the function to be interpolated. The natural scaling for the bispectrum
//               (both the templates and the second-order one) is given by Cl1*Cl2 + Cl2*Cl3 + Cl1*Cl3, which we adopt
//               as a window function. When available, we always use the C_l's for the temperature, as they are approximately
//               the same order of magnitude for all l's, contrary to those for polarization which are very small for l<200. */
//               values[index_l2_l3][0] = pbi->bispectra[index_bt][X_][Y_][Z_][index_l1_l2_l3] / inverse_window;
//
//               /* Debug - print the values at the nodes, quick check */
//               // if ((l1==1000)) {
//               // // if ((l1==63)&&(l2==38)&&(l3==27)) {
//               //   printf ("AT THE NODE: %s_%s(%d,%d,%d) = %g\n",
//               //     pbi->bt_labels[index_bt], pbi->bfff_labels[X_][Y_][Z_], l1, l2, l3,
//               //     values[index_l2_l3][0]);
//               // }
//
//               /* Debug - print to file the effect of the window function, as a function of (l2,l3).
//               Plot with splot [:] "interp.dat" u 1:2:(($3)) lw 2 */
//               // if ((l1==6) && (X_==0) && (Y_==0) && (Z_==0)) {
//               //   fprintf (stderr, "%5d %5d %12.4g %12.4g %12.4g %s %s %s %s\n",
//               //     l2, l3, pbi->bispectra[index_bt][X_][Y_][Z_][index_l1_l2_l3], inverse_window, values[index_l2_l3][0],
//               //     pbi->bt_labels[index_bt], pbi->bf_labels[X_], pbi->bf_labels[Y_], pbi->bf_labels[Z_]);
//               // }
//
//               /* Debug - print to file the effect of the window function, as a function of (l1,l2)
//               Plot with splot [:] "interp.dat" u 1:2:(($3)) lw 2 */
//               // if ((l3==1770) && (X_==1) && (Y_==1) && (Z_==0)) {
//               //   fprintf (stderr, "%5d %5d %12.4g %12.4g %12.4g %s %s %s %s\n",
//               //     l1, l2, pbi->bispectra[index_bt][X_][Y_][Z_][index_l1_l2_l3], inverse_window, values[index_l2_l3][0],
//               //     pbi->bt_labels[index_bt], pbi->bf_labels[X_], pbi->bf_labels[Y_], pbi->bf_labels[Z_]);
//               // }
//
//               /* 2nd to 3rd arguments are the coordinates */
//               values[index_l2_l3][1] = (double)(l2);
//               values[index_l2_l3][2] = (double)(l3);
//
//               /* Go to the next (l2,l3) pair */
//               index_l2_l3++;
//
//             } // for(index_l3)
//           } // for(index_l2)
//
//           // ---------------------------------------------------------------------------
//           // -                            Generate the meshes                          -
//           // ---------------------------------------------------------------------------
//
//           for (int index_mesh=0; index_mesh < pfi->n_meshes; ++index_mesh) {
//
//             /* Number of (l2,l3) nodes where the bispectrum is known */
//             meshes[index_ft][X][Y][Z][index_mesh]->n_points = n_points;
//
//             /* Since the grid is shared between different bispectra, it needs to be computed only for the first one.
//             The first one is not necessarily index_ft=0 because analytical bispectra don't need to be interpolated
//             at all */
//             if ((index_ft == pfi->first_non_analytical_index_ft) && (X == 0) && (Y == 0) && (Z == 0)) {
//               meshes[index_ft][X][Y][Z][index_mesh]->compute_grid = _TRUE_;;
//             }
//             else {
//               meshes[index_ft][X][Y][Z][index_mesh]->compute_grid = _FALSE_;
//               meshes[index_ft][X][Y][Z][index_mesh]->grid_2D =
//                 meshes[pfi->first_non_analytical_index_ft][0][0][0][index_mesh]->grid_2D;
//             }
//
//             printf_log_if (pfi->fisher_verbose, 3,
//               "     * computing mesh for bispectrum %s_%s(%d)\n",
//               pbi->bt_labels[index_bt], pfi->ffff_labels[X][Y][Z], index_mesh);
//
//             /* Skip the fine grid if the turnover point is smaller than than the smallest multipole */
//             if ((index_mesh == 0) && (pfi->l_turnover <= pbi->l[0])) {
//               printf_log_if (pfi->fisher_verbose, 3,
//                 "      \\ fine grid not needed because l_turnover <= l_min (%d <= %d)\n",
//                 pfi->l_turnover, pbi->l[0]);
//               continue;
//             }
//
//             /* Skip the coarse grid if the turnover point is larger than than the largest multipole */
//             if ((index_mesh == 1) && (pfi->l_turnover > pbi->l[pbi->l_size-1])) {
//               printf_log_if (pfi->fisher_verbose, 3,
//                 "      \\ coarse grid not needed because l_turnover > l_max (%d > %d)\n",
//                 pfi->l_turnover, pbi->l[pbi->l_size-1]);
//               continue;
//             }
//
//             /* Generate the mesh. This function is already parallelised. */
//             class_call (mesh_2D_sort (
//                           meshes[index_ft][X][Y][Z][index_mesh],
//                           values),
//               pfi->error_message,
//               pfi->error_message);
//
//             printf_log_if (pfi->fisher_verbose, 3,
//               "      \\ allocated (grid,mesh)=(%g,%g) MBs\n",
//               meshes[index_ft][X][Y][Z][index_mesh]->n_allocated_in_grid*8/1e6,
//               meshes[index_ft][X][Y][Z][index_mesh]->n_allocated_in_mesh*8/1e6);
//
//           } // for(index_mesh)
//         } // for(Z)
//       } // for(Y)
//     } // for(X)
//   } // for(index_ft)
//
//   /* We do not need the node values anymore */
//   for (long int index_l2_l3=0; index_l2_l3 < n_points; ++index_l2_l3)
//     free (values[index_l2_l3]);
//
//   free (values);
//
//   return _SUCCESS_;
//
// }



/**
 * Interpolate in a specific (l1,l2,l3) configuration the bispectrum corresponding to the
 * row of the Fisher matrix 'index_ft' and to the fields X,Y,Z, using 2D interpolation.
 *
 * Note that the returned value needs to be divided by the window function. We do not do it
 * inside the interpolate function for optimization purposes.
 *
 */
 
// int fisher_interpolate_bispectrum_mesh_2D (
//     struct bispectra * pbi,
//     struct fisher * pfi,
//     int index_ft,
//     double l3, double l2, double l1,
//     int X, int Y, int Z,
//     struct interpolation_mesh ** mesh,
//     double * interpolated_value
//     )
// {
//
//   /* Use the fine mesh when all of the multipoles are small. */
//   if ((l1<pfi->l_turnover) && (l2<pfi->l_turnover) && (l3<pfi->l_turnover)) {
//     class_call (mesh_2D_int (mesh[0], l2, l3, interpolated_value),
//     pfi->error_message, pfi->error_message);
//   }
//   /* Use the coarse mesh when any of the multipoles is large. */
//   else {
//     class_call (mesh_2D_int (mesh[1], l2, l3, interpolated_value),
//     pfi->error_message, pfi->error_message);
//   }
//
// #ifdef DEBUG
//   /* Check for nan's and crazy values. A value is crazy when it is much larger than
//   the characteristic scale for a bispectrum, A_s*A_s~1e-20 */
//   if (isnan(*interpolated_value) || (fabs(*interpolated_value)>1) )
//     printf ("@@@ WARNING: Interpolated b(%g,%g,%g) = %g for bispectrum %s_%s!!!\n",
//     l1, l2, l3, *interpolated_value,
//     pfi->ft_labels[index_ft],
//     pfi->ffff_labels[X][Y][Z]);
// #endif // DEBUG
//
//   /* Debug - print the interpolated value of the intrinsic bispectrum */
//   // if ((pbi->has_intrinsic) || (pfi->index_ft_intrinsic==_TRUE_)) {
//   //   printf ("%8g %8g %8g %16.7g\n", l1, l2, l3, *interpolated_value);
//   // }
//
//   return _SUCCESS_;
//
// }









/**
 * Load a bispectrum from disk. It is important to keep the index_bt dependence in the argument list,
 * as we might want to selectively load only the secondary bispectra.
 */
int bispectra_load_from_disk(
    struct bispectra * pbi,
    int index_bt
    )
{

  /* Open file for reading */
  class_open (pbi->bispectra_files[index_bt], pbi->bispectra_paths[index_bt], "rb", pbi->error_message);

  /* Print some debug */
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * reading bispectra from disk for index_bt=%d on'%s'\n",
    index_bt, pbi->bispectra_paths[index_bt]);

  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {

        int n_to_read = pbi->n_independent_configurations;
  
        /* Read a chunk with all the independent (l1,l2,l3) triplets to pbi->bispectra[index_bt] */
        int n_read = fread(
                pbi->bispectra[index_bt][X][Y][Z],
                sizeof(double),
                n_to_read,
                pbi->bispectra_files[index_bt]);

        class_test(n_read != n_to_read,
          pbi->error_message,
          "Could not read in '%s' file, read %d entries but expected %d",
            pbi->bispectra_paths[index_bt], n_read, n_to_read);        
      }
    }
  }
  
  /* Close file */
  fclose(pbi->bispectra_files[index_bt]); 

  return _SUCCESS_;
  
}



/**
 * Bispectrum produced by the pi_dot * grad_pi^2 term in the Galileon Lagrangian, 
 * taken from eq. 22 of http://arxiv.org/abs/0905.3746.
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
 * Bispectrum produced by the pi_dot^3 term in the Galileon Lagrangian, 
 * taken from eq. 22 of http://arxiv.org/abs/0905.3746.
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
  
  double pk_one_third_1 = pow(pk_1, one_third);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, one_third);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, one_third);
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
  
  double pk_one_third_1 = pow(pk_1, one_third);
  double pk_two_thirds_1 = pk_one_third_1*pk_one_third_1;

  double pk_one_third_2 = pow(pk_2, one_third);
  double pk_two_thirds_2 = pk_one_third_2*pk_one_third_2;

  double pk_one_third_3 = pow(pk_3, one_third);
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


/*
 * Compute the CMB lensing bispectrum including polarisation.
 *
 * We implement the formula in Eq. 4.5 of Lewis, Challinor & Hanson 2011
 * (http://uk.arxiv.org/abs/1101.2234). This formula is non-perturbative
 * in the sense of Sec. 3.2 (ibidem) and is valid for all (l1,l2,l3)
 * configurations, including non-squeezed ones.
 *
 * With respect to eq. 4.5, in SONG we have i->X1, j->X2, k->X3, and l1<->l3.
 * We have not included the imaginary term in square brackets which is only
 * needed when including B-mode polarisation (not supported yet, TODO).
 */

int bispectra_cmb_lensing_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------
  
  /* When only including temperature, the CMB-lensing bispectrum reduces to
    b^TTT_l1l2l3 = 1/2 * [l2(l2+1) + l3(l3+1) - l1(l1+1)] C_l2^{T\psi} C_l3^{TT} + 5 permutations,
  where \psi is the lensing potential. */
  
  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    /* Temperature-lensing potential correlation */
    double C_l1_Tp = pbi->cls[psp->index_ct_tp][l1-2];
    double C_l2_Tp = pbi->cls[psp->index_ct_tp][l2-2];
    double C_l3_Tp = pbi->cls[psp->index_ct_tp][l3-2];
                        
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
    double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
                  
    /* Use lensed temperature C_l's if available */
    if (pbi->has_lensed_bispectra == _TRUE_) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
      C_l3 = pbi->lensed_cls[ple->index_lt_tt][l3-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l2*(l2+1) + l3*(l3+1) - l1*(l1+1) ) * C_l2_Tp * C_l3
      + ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * C_l3_Tp * C_l2
      + ( l1*(l1+1) + l3*(l3+1) - l2*(l2+1) ) * C_l1_Tp * C_l3
      + ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * C_l3_Tp * C_l1
      + ( l1*(l1+1) + l2*(l2+1) - l3*(l3+1) ) * C_l1_Tp * C_l2
      + ( l2*(l2+1) + l1*(l1+1) - l3*(l3+1) ) * C_l2_Tp * C_l1
    );
  
    /* Debug - print temperature-lensing potential C_l's */
    // if ((l1==pbi->l_max) && (l2==pbi->l_max))
    //   fprintf (stderr, "%4d %10g %10g %10g\n", l3, C_l3_TT, C_l3, C_l3_Tp);
  
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // end of temperature only

  
  // -------------------------------------------------------------------------------
  // -                        Determine field coefficients                         -
  // -------------------------------------------------------------------------------
  
  /* Set the amplitude of the geometrical factor in Eq. 4.4 of Lewis et al. 2011.
  This is either 2 (if the two 3j-symbols add up) or 0 (if they cancel). Whether
  they add up or cancel depends on the parity of the considered field (even for
  T and E and odd for B). */

  double S_X1=2, S_X2=2, S_X3=2;

  /* Uncomment the following to set the quadratic bispectrum to zero based on the
  parity of the involved fields. By default we skip this part, because we want our
  reduced bispectrum to be continuous in order to facilitate its interpolation and
  plotting. */

  // int L = l3-l1-l2;
  //
  // if (pbi->has_bispectra_t == _TRUE_) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_t) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_e == _TRUE_) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b == _TRUE_) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  //   if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
  // }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "CMB-lensing bispectrum for B-modes not implemented yet.");
  

  // -------------------------------------------------------------------------------
  // -                              Obtain the C_l                                 -
  // -------------------------------------------------------------------------------
  
  /* Get the C_l's involving the lensing potential \phi. When implementing the B-modes, remember
  to set them to zero. */
  double C_l1_X1_p = pbi->cls[pbi->index_ct_of_phi_bf[ X1 ]][l1-2];
  double C_l2_X2_p = pbi->cls[pbi->index_ct_of_phi_bf[ X2 ]][l2-2];
  double C_l3_X3_p = pbi->cls[pbi->index_ct_of_phi_bf[ X3 ]][l3-2];
  
  /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
  double C_l3_X1_X3 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X3 ]][l3-2];
  double C_l2_X1_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X2 ]][l2-2];
  double C_l3_X2_X3 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X3 ]][l3-2];
  double C_l1_X2_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X1 ]][l1-2];
  double C_l2_X3_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ X2 ]][l2-2];
  double C_l1_X3_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ X1 ]][l1-2];
    
  /* Use lensed temperature C_l's if available */
  if (pbi->has_lensed_bispectra == _TRUE_) {
    C_l3_X1_X3 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X3 ]][l3-2];
    C_l2_X1_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X2 ]][l2-2];
    C_l3_X2_X3 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X3 ]][l3-2];
    C_l1_X2_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X1 ]][l1-2];
    C_l2_X3_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X3 ][ X2 ]][l2-2];
    C_l1_X3_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X3 ][ X1 ]][l1-2];
  }
  
  
  // -------------------------------------------------------------------------------
  // -                              Get the right 3j                               -
  // -------------------------------------------------------------------------------
    
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];
  int F_X3 = pbi->field_spin[X3];
  
  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l1_l2_l3_FX1_0_mFX1 = 1, threej_l1_l3_l2_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    threej_l1_l2_l3_FX1_0_mFX1 = threej_ratio_20m2;
    threej_l1_l3_l2_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l2_l3_l1_FX2_0_mFX2 = 1, threej_l2_l1_l3_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    threej_l2_l3_l1_FX2_0_mFX2 = threej_ratio_m220;
    threej_l2_l1_l3_FX2_0_mFX2 = threej_ratio_0m22;
  }
    
  double threej_l3_l1_l2_FX3_0_mFX3 = 1, threej_l3_l2_l1_FX3_0_mFX3 = 1;
  if (F_X3==2) {
    threej_l3_l1_l2_FX3_0_mFX3 = threej_ratio_0m22;
    threej_l3_l2_l1_FX3_0_mFX3 = threej_ratio_20m2;
  }

  /* Obtain the geometric factor F^+s_l1l2l3 */
  double F_l1_l2_l3_X1 = 0.25 * ( l2*(l2+1) + l3*(l3+1) - l1*(l1+1) ) * S_X1 * threej_l1_l2_l3_FX1_0_mFX1; /* 1-2-3 */
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l1_l3_X2 = 0.25 * ( l1*(l1+1) + l3*(l3+1) - l2*(l2+1) ) * S_X2 * threej_l2_l1_l3_FX2_0_mFX2; /* 2-1-3 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
  double F_l3_l1_l2_X3 = 0.25 * ( l1*(l1+1) + l2*(l2+1) - l3*(l3+1) ) * S_X3 * threej_l3_l1_l2_FX3_0_mFX3; /* 3-1-2 */
  double F_l3_l2_l1_X3 = 0.25 * ( l2*(l2+1) + l1*(l1+1) - l3*(l3+1) ) * S_X3 * threej_l3_l2_l1_FX3_0_mFX3; /* 3-2-1 */
  
  
  // -------------------------------------------------------------------------------
  // -                              Bispectrum formula                             -
  // -------------------------------------------------------------------------------
  
  /* CMB-lensing bispectrum formula, from Eq. 4.5 of Lewis, Challinor & Hanson 2011. */
  *result = 
      C_l2_X2_p * C_l3_X1_X3 * F_l1_l2_l3_X1   /* 1-2-3 */
    + C_l3_X3_p * C_l2_X1_X2 * F_l1_l3_l2_X1   /* 1-3-2 */
    + C_l1_X1_p * C_l3_X2_X3 * F_l2_l1_l3_X2   /* 2-1-3 */
    + C_l3_X3_p * C_l1_X2_X1 * F_l2_l3_l1_X2   /* 2-3-1 */
    + C_l1_X1_p * C_l2_X3_X2 * F_l3_l1_l2_X3   /* 3-1-2 */
    + C_l2_X2_p * C_l1_X3_X1 * F_l3_l2_l1_X3;  /* 3-2-1 */
                  
  /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "CMB-lensing bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
     l1, l2, l3, *result, exact, diff);
  }

  /* Debug - Print the 3j-symbols */
  // printf ("threej_l3_l2_l1_FX3_0_mFX3 = %g\n", threej_l3_l2_l1_FX3_0_mFX3);
  // printf ("threej_l3_l1_l2_FX3_0_mFX3 = %g\n", threej_l3_l1_l2_FX3_0_mFX3);
  // printf ("threej_l2_l1_l3_FX2_0_mFX2 = %g\n", threej_l2_l1_l3_FX2_0_mFX2);
  // printf ("threej_l2_l3_l1_FX2_0_mFX2 = %g\n", threej_l2_l3_l1_FX2_0_mFX2);
  // printf ("threej_l1_l3_l2_FX1_0_mFX1 = %g\n", threej_l1_l3_l2_FX1_0_mFX1);
  // printf ("threej_l1_l2_l3_FX1_0_mFX1 = %g\n", threej_l1_l2_l3_FX1_0_mFX1);
  // printf("\n");
  
  /* Debug - Print the C_l's */
  // printf ("C_l1_X1_p = %g\n", C_l1_X1_p);
  // printf ("C_l2_X2_p = %g\n", C_l2_X2_p);
  // printf ("C_l3_X3_p = %g\n", C_l3_X3_p);
  // printf ("\n");

  return _SUCCESS_;
  
}


/*
 * Compute the CMB lensing bispectrum kernel including polarisation in the squeezed limit
 * (l3<<l1 and l3<<l2).
 *
 * Here we code just the kernel in Eq. 5.20 of Lewis, Challinor & Hanson 2011
 * (http://uk.arxiv.org/abs/1101.2234); the full squeezed bispectrum is given by
 * the product between the kernel and C_l3^{X3\phi}.
 *
 * With respect to Eq. 5.20 (ibidem) in SONG we use i->X1, j->X2, k->X3 and then we
 * perform the substitution (X1,l1)<->(X3,l3), as in our convention l3 rather than
 * l1 is the smallest multipole. Therefore, this is what we code here:
 *
 *   A^{X1,X2}_{l1l2l3} = \tilde{C}^{X2,X1}_{l1} * F^{X1}_{l1l3l2}
 *                      + \tilde{C}^{X1,X2}_{l2} * F^{X2}_{l2l3l1}
 *
 * We have not included the imaginary term which is only needed when including
 * B-mode polarisation.

 *
 */

int bispectra_cmb_lensing_squeezed_kernel (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------
  
  /* When only including temperature, the CMB-lensing bispectrum reduces to
    b^TTT_l1l2l3 = 1/2 * [l2(l2+1) + l3(l3+1) - l1(l1+1)] C_l2^{T\psi} C_l3^{TT} + 5 permutations,
  where \psi is the lensing potential. */
  
  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
                  
    /* Use lensed temperature C_l's if available */
    if (pbi->has_lensed_bispectra == _TRUE_) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * C_l2
      + ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * C_l1
    );

    /* Uncomment to compute the bispectrum rather than the kernel */
    // ttt *= pbi->cls[psp->index_ct_tp][l3-2];                      
    
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // end of temperature only
  
  
  // -------------------------------------------------------------------------------
  // -                        Determine field coefficients                         -
  // -------------------------------------------------------------------------------
  
  /* Set the amplitude of the geometrical factor in Eq. 4.4 of Lewis et al. 2011.
  This is either 2 (if the two 3j-symbols add up) or 0 (if they cancel). Whether
  they add up or cancel depends on the parity of the considered field (even for
  T and E and odd for B). */

  double S_X1=2, S_X2=2;

  /* Uncomment the following to set the quadratic bispectrum to zero based on the
  parity of the involved fields. By default we skip this part, because we want our
  reduced bispectrum to be continuous in order to facilitate its interpolation and
  plotting. */

  // int L = l3-l1-l2;
  //
  // if (pbi->has_bispectra_t == _TRUE_) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_e == _TRUE_) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b == _TRUE_) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  // }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "CMB-lensing squeezed bispectrum for B-modes not implemented yet.");
  
  // -------------------------------------------------------------------------------
  // -                              Obtain the C_l                                 -
  // -------------------------------------------------------------------------------
  
  /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
  double C_l2_X1_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X2 ]][l2-2];
  double C_l1_X2_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X1 ]][l1-2];
    
  /* Use lensed temperature C_l's if available */
  if (pbi->has_lensed_bispectra == _TRUE_) {
    C_l2_X1_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X1 ][ X2 ]][l2-2];
    C_l1_X2_X1 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X1 ]][l1-2];
  }
  
  // -------------------------------------------------------------------------------
  // -                             Get the right 3j                                -
  // -------------------------------------------------------------------------------
  
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];

  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l1_l2_l3_FX1_0_mFX1 = 1;
  double threej_l1_l3_l2_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    threej_l1_l2_l3_FX1_0_mFX1 = threej_ratio_20m2;
    threej_l1_l3_l2_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l2_l3_l1_FX2_0_mFX2 = 1;
  double threej_l2_l1_l3_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    threej_l2_l3_l1_FX2_0_mFX2 = threej_ratio_m220;
    threej_l2_l1_l3_FX2_0_mFX2 = threej_ratio_0m22;
  }

  /* Obtain the geometric factor F^+s_l1l2l3 */
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1) + l2*(l2+1) - l1*(l1+1) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1) + l1*(l1+1) - l2*(l2+1) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
      

  // -------------------------------------------------------------------------------
  // -                              Bispectrum formula                             -
  // -------------------------------------------------------------------------------
  
  /* Kernel of the CMB-lensing bispectrum in the squeezed limit, from Eq. 5.20 of Lewis,
  Challinor & Hanson 2011. This is simply the general formula with C_l1_X1_p=C_l2_X2_p=0 
  (see bispectra_cmb_lensing_bispectrum()) */
  *result = C_l2_X1_X2 * F_l1_l3_l2_X1
          + C_l1_X2_X1 * F_l2_l3_l1_X2;
                  
  /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "CMB-lensing squeezed bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
     l1, l2, l3, *result, exact, diff);
  }

  return _SUCCESS_;
  
}

/*
 * Compute the CMB lensing bispectrum including polarisation in the squeezed limit, valid
 * when l3<<l1 and l3<<l2. This is given by the kernel computed in 'bispectra_cmb_lensing_squeezed_kernel'
 * times the cross-correlation with the lensing potential, C_l3^{X3\phi} (see Eq. 5.20 of Lewis,
 * Challinor & Hanson 2011 (http://uk.arxiv.org/abs/1101.2234).
 *
 * This formula is non-perturbative (in the sense of Sec. 3.2, ibidem) and is valid only
 * for squeezed configurations, where l3<<l1 and l3<<l2; it is obtained from the general
 * formula (Eq. 4.5, ibidem) by setting C_l1_X1_p=C_l2_X2_p=0. It is an excellent approximation
 * for the CMB-lensing bispectrum, as most of its signal is in these squeezed configurations.
 *
 * Look at the comment above 'bispectra_intrinsic_squeezed_bispectrum' for details
 * about the symmetry properties of this special squeezed bispectrum.
 * 
 */
int bispectra_cmb_lensing_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Compute the kernel */
  class_call (bispectra_cmb_lensing_squeezed_kernel (
                ppr, psp, ple, pbi,
                l1, l2, l3,
                X1, X2, X3,
                threej_ratio_20m2,
                threej_ratio_m220,
                threej_ratio_0m22,
                result),
    pbi->error_message,
    pbi->error_message);

  /* Obtain the bispectrum in the squeezed limit by multiplication with C_l3^{X3\phi} */
  *result *= pbi->cls[pbi->index_ct_of_phi_bf[ X3 ]][l3-2];

  return _SUCCESS_;

}


/** 
 * Squeezed approximation for the local bispectrum with f_nl=1.
 *
 * This is the generalisation to polarisation of the temperature approximation first
 * described in Gangui et al. 1994. See also Komatsu & Spergel 2001 and eq. 2.16 of
 * http://arxiv.org/abs/1201.1010.
 *
 */
int bispectra_local_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only
  sensible choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double cl1_XY = pbi->cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double cl2_XY = pbi->cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Use lensed temperature C_l's if available */
  if (pbi->has_lensed_bispectra == _TRUE_) {
    cl1_XY = pbi->lensed_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];    
    cl2_XY = pbi->lensed_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];    
  }

  /* Obvious generalisation of the TTT result in Eq. 2.16 of http://arxiv.org/abs/1201.1010 */
  *result = 6/5. * cl3_Zz * (cl1_XY + cl2_XY);

  /* TTT result */
  // *result = 6/5. * (
  //     pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
  //   + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
  //   + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2]
  // );      

  /* EEE result */
  // *result = 6/5. * (
  //     pbi->cls[psp->index_ct_ee][l1-2] * pbi->cls[psp->index_ct_ee][l2-2]
  //   + pbi->cls[psp->index_ct_ee][l2-2] * pbi->cls[psp->index_ct_ee][l3-2]
  //   + pbi->cls[psp->index_ct_ee][l3-2] * pbi->cls[psp->index_ct_ee][l1-2]
  // );      

  return _SUCCESS_;

}



/** 
 * Compute the squeezed-limit approximation for the intrinsic bispectrum, as in
 * eq. 4.1 and 4.2 of Lewis 2012 (http://arxiv.org/abs/1204.5018). See also Creminelli
 * et al. 2004 (http://arxiv.org/abs/astro-ph/0405428), Creminelli et al. 2011
 * (http://arxiv.org/abs/1109.1822), Bartolo et al. 2012 (http://arxiv.org/abs/1109.2043).
 * This is the general formula that includes polarisation. With respect to Lewis' formula,
 * (i,l1)->(Z,l3), (j,l2)->(X,l1), (k,l3)->(Y,l2) and \zeta -> z.
 * 
 * ~~~ CONSIDERATIONS THAT APPLY TO ALL "SQUEEZED" BISPECTRA ~~~
 *
 * In this and in the other "squeezed" approximation functions, l3 is taken to be the long-wavelength
 * mode, l1 and l2 the short ones. This choice is preferred because in SONG we loop over (l1,l2,l3)
 * configurations that satisfy the condition l1>=l2>=l3.
 *
 * Therefore, it is crucial that that, in the formulas below, the multipole l3 is associated with the
 * C_l that correlates with the comoving curvature perturbation zeta (for the CMB lensing bispectrum
 * it would be the lensing potential phi) because l3 is the smallest multipole, which describes the
 * long wavelength mode. Also, l3 must be associated with Z (or X3) because our convention
 * for the bispectrum is <X_l1 Y_l2 Z_l3>. If you give l3 to another field, then the Fisher matrix
 * estimator will associate to that field the wrong covariance matrix, and the result will change
 * drastically.
 *
 * Because of the special role played by l3, the bispectrum computed in this and the other 'squeezed'
 * approximation functions is NOT symmetric with respect to an exchange of (l1,X) <-> (l3,Z) or of
 * or (l2,Y) <-> (l3,Z), contrary to the other bispectra, which are computed as <X_l1 Y_l2 Z_l3>.
 * Therefore, for this bispectrum one cannot obtain the configurations outside l1>=l2>=l3 by
 * permuting the XYZ indices, as it is done, for example, in the 'print_bispectra' function.
 * 
 */

int bispectra_intrinsic_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* Uncomment to exclude the polarisation when it is assigned the large-scale mode */
  // if ((pbi->has_bispectra_e==_TRUE_) && (Z==pbi->index_bf_e)) {
  //   *result = 0;
  //   return _SUCCESS_;
  // }

  /* Uncomment to restrict the approximation to squeezed configurations, as in Sec. 6 of Creminelli,
  Pitrou & Vernizzi 2011 (http://arxiv.org/abs/1109.1822). Note that in our case the smallest mode is l3,
  while in that paper it is l1. IMPORTANT: setting a bispectrum to zero might screw up some matrix
  inversions done in the Fisher module, especially when computing the lensing variance of the 
  bispectrum (pfi->include_lensing_effects==_TRUE_). */
  // if (!((l3<=100) && (l2>=10*l3))) {
  //   *result = 0;
  //   return _SUCCESS_;
  // }

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
  choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Use lensed temperature C_l's if available */
  if (pbi->has_lensed_bispectra == _TRUE_) {
    dcl1_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l1-2];    
    dcl2_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l2-2];    
  }
          
  /* Ricci focussing in Lewis 2012 (eq. 4.1) */
  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);

  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
  double bolometric_T_lewis_redshift = 0;
          
  if (pbi->has_bispectra_t == _TRUE_) {
    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf_bf[Z][pbi->index_bf_t]][l3-2];
    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf_bf[X][pbi->index_bf_t]][l1-2];
    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    if (pbi->has_lensed_bispectra == _TRUE_) {
      if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[X][pbi->index_bf_t]][l1-2];
      if (X == pbi->index_bf_t) cl2_Yt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    }
    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
  }
          
  class_test ((pbi->has_bispectra_e == _TRUE_) &&
    ((X == pbi->index_bf_e) && (Y == pbi->index_bf_e))
    && (bolometric_T_lewis_redshift!=0),
    pbi->error_message,
    "anisotropic redshifting should vanish (eq. 4.2 of Lewis 2012)", bolometric_T_lewis_redshift);
          
  /* Sum of Ricci focussing and redshift modulation */
  *result = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

  return _SUCCESS_;

}


/**
 * Same as above, but using the unlensed power spectrum, which gives more signal as the
 * lensed C_l's oscillate less. See 'bispectra_intrinsic_squeezed_bispectrum' for details.
 */

int bispectra_intrinsic_squeezed_unlensed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Test that l3 is the smallest multipole */
  class_test ((l1<l3) || (l2<l3),
    pbi->error_message,
    "in all squeezed approximations, make sure l3 is the smallest multipole");

  /* We take l3 to be the long wavelength and l1 and l2 the short ones. This is the only sensible
  choice as the l-loop we are into is constructed to have l1>=l2>=l3. */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Ricci focussing in Lewis 2012 (eq. 4.1) */
  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);

  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
  double bolometric_T_lewis_redshift = 0;
          
  if (pbi->has_bispectra_t == _TRUE_) {
    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf_bf[Z][pbi->index_bf_t]][l3-2];
    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf_bf[X][pbi->index_bf_t]][l1-2];
    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
  }
          
  /* Sum of Ricci focussing and redshift modulation */
  *result = bolometric_T_lewis_ricci + bolometric_T_lewis_redshift;

  return _SUCCESS_;

}


/** 
 * Compute the quadratic correction to the bispectrum.
 *
 * At second order, the CMB temperature is obtained from the photon distribution
 * function by adding a correction quadratic in the first-order distribution (see
 * Eq. 3.12 of arXiv:1401.3296).
 *
 * At the bispectrum level this extra term translates to a four-point function
 * in the first-order perturbations, which in turn is expressed in terms of
 * products of C_l.
 *
 * The same correction term appears for the variable transformation that is used
 * to treat the redshift therm (see eq. 3.5 ibidem), the so-called
 * delta_tilde transformation.
 *
 * In both cases, the correction has the following form (eq. 3.6, 3.7 and 3.9):
 * 
 *   QC_{l1 l2 l3} = S_X3 * i^(L+V_X3) * ( l'  l''  |    l3 ) * (  l'  l'' | l3 )  
 *                                       ( 0   F_X3 | -F_X3 )   (  m'  m'' | m3 ) 
 *
 *                   * <a^X1_l1m1 * a^X2_l2m2 * a^I_l'm' * a^T_X3_l''m''>
 *
 *                   + 2 permutations (1->2->3)
 *   
 * where L=l3-l'-l'' and the a_lm are first-order. The permutations go over 1->2->3
 * and refer to the position of the second-order perturbation in the bispectrum; in
 * the formula above, the positioning is <a^(1)X1_l1m1 * a^(1)X2_l2m2 * a^(2)X3_l3m3>.
 *
 * The coefficients take different values according to which field is X3:
 * 
 *   X3=I -> F=0, S=1, V_X3=0,  T_X3=I
 *   X3=E -> F=2, S=2, V_X3=0,  T_X3=E
 *   X3=B -> F=2, S=2, V_X3=-1, T_X3=E
 * 
 * For X3=E, the sum over l' and l'' only includes EVEN values of l3-l'-l'', for X3=B
 * it only includes ODD values of l3-l'-l''.
 * 
 * By employing Wick's theorem, the above can be expressed in terms of the angular
 * power spectrum of the CMB, C_l = <a_lm a^*_lm> (note that the a_lm's already
 * include the 1/4 factor):
 * 
 *   QC_{l1 l2 l3} = 4 * B_{l1 l2 l3} * i^V_X3 * G^m1m2m3_l1l2l3 * S_X3
 *
 *                     * [ 
 *                           ( l1    l2     l3 ) * C^{X1,I}_l1 * C^{X2,X3}_l2
 *                           ( 0   F_X3  -F_X3 )
 *
 *                         + (   l1  l2     l3 ) * C^{X1,X3}_l1 * C^{X2,I}_l2
 *                           ( F_X3   0  -F_X3 )
 *
 *                       ] + 2 permutations (1->2->3) ,
 * 
 * where
 * 
 *   B_{l1 l2 l3} = sqrt((2*l1+1)*(2*l2+1)*(2*l3+1)/(4*_PI_)) * ( l1    l2   l3 )
 *                                                              ( m1    m2   m3 ) .
 * 
 * Below we compute the quadratic correction term using this general formula and store
 * it in pbi->bispectra[pbi->index_bt_quadratic].
 * 
 * For intensity (X1=X2=X3=I), the formula reduces simply to
 *
 *    8 * [ C_l1*C_l2 + C_l1*C_l2 + C_l1*C_l2 ]
 * 
 * while any contribution from the second-order B-modes has V_X=-1 and is therefore
 * purely imaginary. Furthermore, the B-mode contribution only includes one term,
 * the one where a^B_lm is second order, as we ignore the first-order B-modes. In
 * any case, we do not compute the quadratic term for a B-mode spectrum.
 */

int bispectra_quadratic_correction (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{
  
  // --------------------------------------------------------------------------------------
  // -                              Temperature-only formula                              -
  // --------------------------------------------------------------------------------------

  /* The quadratic (reduced) bispectrum for TTT is very simple: 
    b^TTT_l1l2l3 = 8 * (C_l1*C_l2 + 2 perms) */

  double ttt;
  
  if (pbi->has_bispectra_t == _TRUE_) {
  
    double C_l1_TT = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2_TT = pbi->cls[psp->index_ct_tt][l2-2];
    double C_l3_TT = pbi->cls[psp->index_ct_tt][l3-2];
  
    ttt = 8 * (C_l1_TT*C_l2_TT + C_l1_TT*C_l3_TT + C_l2_TT*C_l3_TT);
    
    /* If only temperature is requested, skip what follows exit from the 'quadratic' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      goto end_quadratic; 
    }
  }

  // -------------------------------------------------------------------------------
  // -                           Determine field coefficients                      -
  // -------------------------------------------------------------------------------
  
  /* Determine the fields T_X1, T_X2 and T_X3 that go in the C_l (T_I=I, T_E=E, T_B=E).
  and the indices of the cross power spectra between X1, X2, X3 and I. These have
  to be set by hand, rather than using the array pbi->index_ct_of_bf_bf because
  pbi->index_ct_of_bf_bf only contains information on the fields appearing in one of
  the requested bispectrum. For example, if you only request EEE, then
  pbi->index_ct_of_bf_bf does not contain information about <ET>, because no
  bispectrum containing T is requested */
  
  int T_X1=0, T_X2=0, T_X3=0;
  int index_ct_X1_I, index_ct_X2_I, index_ct_X3_I;  
  double S_X1=0, S_X2=0, S_X3=0;

  if (pbi->has_bispectra_t == _TRUE_) {
    if (X1 == pbi->index_bf_t) {S_X1 = 1; T_X1 = pbi->index_bf_t; index_ct_X1_I = psp->index_ct_tt;}
    if (X2 == pbi->index_bf_t) {S_X2 = 1; T_X2 = pbi->index_bf_t; index_ct_X2_I = psp->index_ct_tt;}
    if (X3 == pbi->index_bf_t) {S_X3 = 1; T_X3 = pbi->index_bf_t; index_ct_X3_I = psp->index_ct_tt;}
  }
  if (pbi->has_bispectra_e == _TRUE_) {
    if (X1 == pbi->index_bf_e) {S_X1 = 2; T_X1 = pbi->index_bf_e; index_ct_X1_I = psp->index_ct_te;}
    if (X2 == pbi->index_bf_e) {S_X2 = 2; T_X2 = pbi->index_bf_e; index_ct_X2_I = psp->index_ct_te;}
    if (X3 == pbi->index_bf_e) {S_X3 = 2; T_X3 = pbi->index_bf_e; index_ct_X3_I = psp->index_ct_te;}
  }
  if (pbi->has_bispectra_b == _TRUE_) { /* Note that <TB> vanishes, hence the negative values */
    if (X1 == pbi->index_bf_b) {S_X1 = 2; T_X1 = pbi->index_bf_e; index_ct_X1_I = -1;}
    if (X2 == pbi->index_bf_b) {S_X2 = 2; T_X2 = pbi->index_bf_e; index_ct_X2_I = -1;}
    if (X3 == pbi->index_bf_b) {S_X3 = 2; T_X3 = pbi->index_bf_e; index_ct_X3_I = -1;}
  }

  /* Uncomment the following to set the quadratic bispectrum to zero based on the
  parity of the involved fields. By default we skip this part, because we want our
  reduced bispectrum to be continuous in order to facilitate its interpolation and
  plotting. */

  // int L = l3-l1-l2;
  //
  // if (pbi->has_bispectra_t == _TRUE_) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?1:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?1:0;
  //   if (X3 == pbi->index_bf_t) S_X3 = (L%2==0)?1:0;
  // }
  // if (pbi->has_bispectra_e == _TRUE_) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b == _TRUE_) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  //   if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
  // }
              
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b == _TRUE_,
    pbi->error_message,
    "quadratic correction for B-mode bispectrum not implemented yet.");
    
  
  // -------------------------------------------------------------------------------
  // -                               Obtain the C_l                                -
  // -------------------------------------------------------------------------------
  
  /* Get the C_l's. When implementing the B-modes, remember to 
  set the C_l's to zero, so that the only contribution comes from the
  bispectrum with the second-order a^B_lm.  */
  /* TODO: do we need the lensed C_l's? */
  /* TODO: CHECK FACTORS 4!!!! */
    
  double C_l1_X1_I   = pbi->cls[index_ct_X1_I][l1-2];
  double C_l1_X1_TX2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ T_X2 ]][l1-2];
  double C_l1_X1_TX3 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ T_X3 ]][l1-2];
              
  double C_l2_X2_I   = pbi->cls[index_ct_X2_I][l2-2];
  double C_l2_X2_TX1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ T_X1 ]][l2-2];
  double C_l2_X2_TX3 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ T_X3 ]][l2-2];
  
  double C_l3_X3_I   = pbi->cls[index_ct_X3_I][l3-2];                  
  double C_l3_X3_TX1 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ T_X1 ]][l3-2];
  double C_l3_X3_TX2 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ T_X2 ]][l3-2];
    
  
  // -------------------------------------------------------------------------------
  // -                              Get the right 3j                               -
  // -------------------------------------------------------------------------------
  
  /* Spin of the fields */
  int F_X1 = pbi->field_spin[X1];
  int F_X2 = pbi->field_spin[X2];
  int F_X3 = pbi->field_spin[X3];
  
  /* Obtain the needed 3j ratios (TODO: what about B-modes? */
  double threej_l2_l3_l1_0_FX1_mFX1 = 1, threej_l2_l3_l1_FX1_0_mFX1 = 1;
  if (F_X1==2) {
    threej_l2_l3_l1_0_FX1_mFX1 = threej_ratio_20m2;
    threej_l2_l3_l1_FX1_0_mFX1 = threej_ratio_m220;
  }

  double threej_l3_l1_l2_0_FX2_mFX2 = 1, threej_l3_l1_l2_FX2_0_mFX2 = 1;
  if (F_X2==2) {
    threej_l3_l1_l2_0_FX2_mFX2 = threej_ratio_m220;
    threej_l3_l1_l2_FX2_0_mFX2 = threej_ratio_0m22;
  }
  
  double threej_l1_l2_l3_0_FX3_mFX3 = 1, threej_l1_l2_l3_FX3_0_mFX3 = 1;
  if (F_X3==2) {
    threej_l1_l2_l3_0_FX3_mFX3 = threej_ratio_0m22;
    threej_l1_l2_l3_FX3_0_mFX3 = threej_ratio_20m2;
  }

  // -------------------------------------------------------------------------------
  // -                              Bispectrum formula                             -
  // -------------------------------------------------------------------------------

  /* The sum includes three terms, corresponding to the three possible types
  of second-order perturbations. The first term, involving T_X1, for example,
  refers to the term where a^X1_lm is second order. */
  *result = 4 * S_X1 * (threej_l2_l3_l1_0_FX1_mFX1*C_l2_X2_I*C_l3_X3_TX1 + threej_l2_l3_l1_FX1_0_mFX1*C_l2_X2_TX1*C_l3_X3_I)
          + 4 * S_X2 * (threej_l3_l1_l2_0_FX2_mFX2*C_l3_X3_I*C_l1_X1_TX2 + threej_l3_l1_l2_FX2_0_mFX2*C_l3_X3_TX2*C_l1_X1_I)
          + 4 * S_X3 * (threej_l1_l2_l3_0_FX3_mFX3*C_l1_X1_I*C_l2_X2_TX3 + threej_l1_l2_l3_FX3_0_mFX3*C_l1_X1_TX3*C_l2_X2_I);

  /* Check that for <TTT> the correction is equal to 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) */
  if ((pbi->has_bispectra_t==_TRUE_) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "quadratic bispectrum correction does not simplify to 8 (C_l1*C_l2 + perm); b=%g, exact=%g, diff=%g",
     *result, exact, diff);
  }

  /* Debug. Print the quadratic contribution */
  // printf ("B_quadratic_%3s[%3d,%3d,%3d] = %g\n",
  //   pbi->bfff_labels[X][Y][Z], l1, l2, l3, *result);

  end_quadratic: ;
  
  return _SUCCESS_;

}


/**
 * Return an l-dependent normalisation suitable to visualise the bispectra,
 * especially the intrinsic bispectrum.
 *
 * This function will be called in bispectra_output() to include the 
 * normalisation as a column in the output files.
 *
 * The intrinsic bispectrum divided by this factor is of the same order of 
 * magnitude as fnl_bias, that is, the bias induced by the intrinsic
 * bispectrum on a measurement of the local bispectrum (fnl_bias). This
 * is the same normalisation we used to plot the results in
 * http://arxiv.org/abs/1406.2981.
 *
 * We define the normalisation as
 *
 *      6/5 * cl1_Xz * (cl2_YY + cl3_ZZ)
 *    + 6/5 * cl2_Yz * (cl3_ZZ + cl1_XX)
 *    + 6/5 * cl3_Zz * (cl1_XX + cl2_YY);
 *
 * where z is the curvature perturbation zeta and the C_l are unlensed. By using
 * quadratic combinations of C_l, we factor out from the bispectra the dependence
 * on the primordial power spectrum.
 * 
 * For the TTT and EEE bispectra, the normalisation is the symmetrised version of
 * the squeezed limit of the local bispectrum (which we have coded in
 * bispectra_local_squeezed_bispectrum()). It follows that the normalised local
 * bispectrum should be of order unity, and the normalised intrinsic bispectrum
 * should be of order fnl_bias.
 *
 * This normalisation has the disadvantage of crossing the zero whenever the
 * zeta C_l do, that is, several times after l~60 (Fig. 3 of arxiv:1204.5018).
 * In bispectra_normalisation_positive() we have coded an alternative formula
 * that is strictly positive; both normalisations will be included in the output
 * files.
 */

int bispectra_normalisation (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double not_used_1,
     double not_used_2,
     double not_used_3,
     double * result
     )
{

  /* Extract the C_l */
  double cl1_Xz = pbi->cls[pbi->index_ct_of_zeta_bf[ X ]][l1-2];
  double cl2_Yz = pbi->cls[pbi->index_ct_of_zeta_bf[ Y ]][l2-2];
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];

  double cl1_XX = pbi->cls[pbi->index_ct_of_bf_bf[X][X]][l1-2];
  double cl2_YY = pbi->cls[pbi->index_ct_of_bf_bf[Y][Y]][l2-2];
  double cl3_ZZ = pbi->cls[pbi->index_ct_of_bf_bf[Z][Z]][l3-2];

  /* Build the normalisation */
  *result =   6/5. * cl1_Xz * (cl2_YY + cl3_ZZ)
            + 6/5. * cl2_Yz * (cl3_ZZ + cl1_XX)
            + 6/5. * cl3_Zz * (cl1_XX + cl2_YY);  
  
  return _SUCCESS_;

}


/**
 * Return an l-dependent, strictly positive normalisation factor suitable to
 * visualise the bispectrum.
 *
 * This function will be called in bispectra_output() to include the 
 * normalisation as a column in the output files.
 *
 * We define the normalisation as
 *
 *  -24/5. * (cl1_XX*cl2_YY + cl2_YY*cl3_ZZ + cl3_ZZ*cl1_XX)
 *
 * where the C_l are taken to be unlensed. By using quadratic combinations of
 * C_l, we factor out from the bispectra the dependence on the primordial power
 * spectrum.
 * 
 * This choice of normalisation is particularly advantageous because it never
 * crosses the zero. Furthermore, for TTT and for l1 << 200, the formula
 * is the symmetrised version of the squeezed limit of the local bispectrum, (which
 * we have coded in bispectra_local_squeezed_bispectrum()). It follows that the
 * normalised local bispectrum should be of order unity, and the normalised
 * intrinsic bispectrum should be of order fnl_bias.
 *
 * In bispectra_normalisation() we have coded an alternative normalisation
 * factor that approximates the squeezed local bispectrum also for EEE.
 * Both normalisations will be included in the output
 * files.
 */

int bispectra_normalisation_positive (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double not_used_1,
     double not_used_2,
     double not_used_3,
     double * result
     )
{

  /* Extract the C_l */
  double cl1_XX = pbi->cls[pbi->index_ct_of_bf_bf[X][X]][l1-2];
  double cl2_YY = pbi->cls[pbi->index_ct_of_bf_bf[Y][Y]][l2-2];
  double cl3_ZZ = pbi->cls[pbi->index_ct_of_bf_bf[Z][Z]][l3-2];

  /* Build the normalisation */
  *result = -24/5. * (cl1_XX*cl2_YY + cl2_YY*cl3_ZZ + cl3_ZZ*cl1_XX);
  
  return _SUCCESS_;

}



/** 
 * Simple bispectrum used for testing purposes:
 *  b_l1l2l3 = 6 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) * cos((l1+l2+l3)/50).
 * It is only defined for temperature.
 */
int bispectra_cosine_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{                   
  
  /* Get the squeezed approximations for the local bispectrum */
  class_call (bispectra_local_squeezed_bispectrum (
                ppr, psp, ple, pbi,
                l1, l2, l3,
                X1, X2, X3,
                threej_ratio_20m2,
                threej_ratio_m220,
                threej_ratio_0m22,
                result),
    pbi->error_message,
    pbi->error_message);

  *result *= cos((l1+l2+l3)/50.);

  return _SUCCESS_;
  
}




/** 
 * Window function for the local (l1,l2,l3) bispectrum.
 *
 * The natural scaling for the bispectrum (both the templates and the second-order one) is given by a product 
 * of two power spectra. When available, we always use the C_l's for the temperature as, when multiplied by l^2,
 * they are approximately the same order of magnitude for all l's, contrary to those for polarisation which are
 * very small for l<200.
 *
 */
int bispectra_local_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  /* Interpolate all bispectra with temperature's C_l1*C_l2 */
  if (pbi->has_bispectra_t == _TRUE_) {
    *result = 
        pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
      + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
      + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  }
  else if (pbi->bf_size == 1) {
    *result =
        pbi->cls[0][l1-2] * pbi->cls[0][l2-2]
      + pbi->cls[0][l2-2] * pbi->cls[0][l3-2]
      + pbi->cls[0][l3-2] * pbi->cls[0][l1-2];
  }
  else {
    class_stop (pbi->error_message, "bispectra with two non-temperature fields are not implemented yet.\n");
  }

  /* Interpolate only temperature */
  // if ((pbi->has_bispectra_t == _TRUE_) && (X==pbi->index_bf_t) && (X==Y) && (X==Z)) {
  //   *result = 
  //       pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
  //     + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
  //     + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  // }
  // else {
  //   *result = 1;
  // }

  return _SUCCESS_;

}



/** 
 * Window function for the intrinsic (l1,l2,l3) bispectrum.
 *
 */
int bispectra_intrinsic_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
     double * result
     )
{

  if (pbi->has_bispectra_t == _TRUE_) {
    *result = 
        pbi->cls[psp->index_ct_tt][l1-2] * pbi->cls[psp->index_ct_tt][l2-2]
      + pbi->cls[psp->index_ct_tt][l2-2] * pbi->cls[psp->index_ct_tt][l3-2]
      + pbi->cls[psp->index_ct_tt][l3-2] * pbi->cls[psp->index_ct_tt][l1-2];
  }
  else if (pbi->bf_size == 1) {
    *result =
        pbi->cls[0][l1-2] * pbi->cls[0][l2-2]
      + pbi->cls[0][l2-2] * pbi->cls[0][l3-2]
      + pbi->cls[0][l3-2] * pbi->cls[0][l1-2];
  }
  else {
    class_stop (pbi->error_message, "bispectra with two non-temperature fields are not implemented yet.\n");
  }

  return _SUCCESS_;

}





/**
 * Initialise the interpolation mesh for the input bispectrum.
 *
 * The mesh created by this function can be used to interpolate the bispectrum
 * in the bidimensional slice with fixed l1. The interpolated dimensions, l2
 * and l3, need to be larger than l1.
 *
 * We choose to interpolate the two largest multipoles in the bispectrum
 * because they are the slowest varying directions.
 *
 * The mesh must be already allocated with the compute_grid variables set.
 */
int bispectra_init_interpolation_mesh_2D (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int index_bt,
        int X, int Y, int Z,
        int index_l1, /**< Input: fixed multipole in the interpolation mesh */
        int ** mesh_grid, /**< Input: interpolation grid; set to NULL to recompute it */
        struct interpolation_mesh * mesh /**< Output: the initialised interpolation mesh. It must be
                                         already allocated with the compute_grid variables set. */
        )
{

  int l1 = pbi->l[index_l1];

  /* Count the number of (l2,l3) nodes in the considered l1-slice of the bispectrum */
  mesh->n_points = 0;
  for (int index_l2=index_l1; index_l2 < pbi->l_size; ++index_l2) {
    int index_l3_min = MAX (index_l2, pbi->index_l_triangular_min[index_l1][index_l2]);
    int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
    for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3)
      mesh->n_points++;
  }

  printf_log_if (pbi->bispectra_verbose, 3,
    "     * the l1=%d slice has %ld point%s\n",
    l1, mesh->n_points, ((mesh->n_points==1)?"":"s"));

  /* Allocate the array that will contain the re-arranged bispectra */
  double ** values;
  class_alloc (values, mesh->n_points*sizeof(double *), pbi->error_message);
  for (int index_l2_l3=0; index_l2_l3 < mesh->n_points; ++index_l2_l3)
    class_calloc (values[index_l2_l3], 3, sizeof(double), pbi->error_message);


  // -------------------------------------------------------------------------------
  // -                           Rearrange the bispectrum                          -
  // -------------------------------------------------------------------------------

  /* Identifier for the current (l2,l3) node */
  long int index_l2_l3 = 0;

  for (int index_l2=index_l1; index_l2 < pbi->l_size; ++index_l2) {

    int l2 = pbi->l[index_l2];
    int index_l3_min = MAX (index_l2, pbi->index_l_triangular_min[index_l1][index_l2]);
    int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];

    for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

      int l3 = pbi->l[index_l3];

      /* Determine the index of this particular index_l1_l2_l3, taking into account that now
      l1 is the smallest multipole */
      int index_l1_max = MIN (index_l2, pbi->index_l_triangular_max[index_l2][index_l3]);
      long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l3][index_l3-index_l2][index_l1_max-index_l1];

      /* By default, assume that no window function is needed */
      double inverse_window = 1;
      double dump;

      /* Compute the value of the window function for this (l1,l2,l3) */
      if (pbi->window_function[index_bt] != NULL) {

        class_call ((*pbi->window_function[index_bt]) (
                      ppr, psp, ple, pbi,
                      l3, l2, l1, /* smallest one goes in third position  */
                      X, Y, Z,
                      dump,
                      dump,
                      dump,
                      &inverse_window),
          pbi->error_message,
          pbi->error_message);
      }

      /* The 1st element of the values array must be the function to be interpolated */
      values[index_l2_l3][0] = pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] / inverse_window;

      /* Debug - print the values at the nodes, quick check */
      // if ((l1==63)&&(l2==38)&&(l3==27)) {
      //   printf ("AT THE NODE: %s_%s(%d,%d,%d) = %g\n",
      //     pbi->bt_labels[index_bt], pbi->bfff_labels[X_][Y_][Z_], l1, l2, l3,
      //     values[index_l2_l3][0]);
      // }

      /* Debug - print to file the effect of the window function, as a function of
      (l2,l3). Plot with splot [:] "interp.dat" u 1:2:(($3)) lw 2 */
      // if ((l1==6) && (X_==0) && (Y_==0) && (Z_==0)) {
      //   fprintf (stderr, "%5d %5d %12.4g %12.4g %12.4g %s %s %s %s\n",
      //     l2, l3, pbi->bispectra[index_bt][X_][Y_][Z_][index_l1_l2_l3], inverse_window, values[index_l2_l3][0],
      //     pbi->bt_labels[index_bt], pbi->bf_labels[X_], pbi->bf_labels[Y_], pbi->bf_labels[Z_]);
      // }

      /* Debug - print to file the effect of the window function, as a function of (l1,l2)
      Plot with splot [:] "interp.dat" u 1:2:(($3)) lw 2 */
      // if ((l3==1770) && (X_==1) && (Y_==1) && (Z_==0)) {
      //   fprintf (stderr, "%5d %5d %12.4g %12.4g %12.4g %s %s %s %s\n",
      //     l1, l2, pbi->bispectra[index_bt][X_][Y_][Z_][index_l1_l2_l3], inverse_window, values[index_l2_l3][0],
      //     pbi->bt_labels[index_bt], pbi->bf_labels[X_], pbi->bf_labels[Y_], pbi->bf_labels[Z_]);
      // }

      /* 2nd to 3rd arguments are the coordinates */
      values[index_l2_l3][1] = (double)(l2);
      values[index_l2_l3][2] = (double)(l3);

      /* Go to the next (l2,l3) pair */
      index_l2_l3++;

    } // for(index_l3)
  } // for(index_l2)


  // -------------------------------------------------------------------------------
  // -                               Generate the mesh                             -
  // -------------------------------------------------------------------------------

  /* Compute the grid if the user didn't provide one */
  mesh->compute_grid = (mesh_grid==NULL?_TRUE_:_FALSE_);
  if (mesh->compute_grid == _FALSE_)
    mesh->grid_2D = mesh_grid;

  /* Generate the mesh. This function is already parallelised. */
  class_call (mesh_2D_sort (
                mesh,
                values),
    pbi->error_message,
    pbi->error_message);

  printf_log_if (pbi->bispectra_verbose, 3,
    "      \\ allocated (grid,mesh)=(%g,%g) MBs\n",
    mesh->n_allocated_in_grid*8/1e6, mesh->n_allocated_in_mesh*8/1e6);

  /* We do not need the node values anymore */
  for (long int index_l2_l3=0; index_l2_l3 < mesh->n_points; ++index_l2_l3)
    free (values[index_l2_l3]);

  free (values);

  return _SUCCESS_;

}




