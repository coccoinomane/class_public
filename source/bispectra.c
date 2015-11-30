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
 * If the user specified 'store_bispectra=yes', the module will save all the 
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
 * configurations. This strategy would require interpolating the transfer
 * functions in l3, though.
 *
 * Created by Guido W. Pettinari on 19.07.2012.
 * Last modified by Guido W. Pettinari on 23.10.2015
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
 *    bispectra_store().
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
  if (!ppt->has_cmb_bispectra) {

    pbi->has_bispectra = _FALSE_;

    printf_log_if (pbi->bispectra_verbose, 0,
      "No bispectra requested. Bispectra module skipped.\n");

    pbi->bt_size = 0;

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

  class_call (bispectra_indices (
                ppr,pba,pth,ppt,
                pbs,ptr,ppm,psp,
                ple,pbi),
    pbi->error_message,
    pbi->error_message);



  // =====================================================================================
  // =                                 Compute bispectra                                 =
  // =====================================================================================
  
  /* Compute the three types of CMB bispectra: analytic, separable and non-separable */
  
  class_call (bispectra_harmonic (
                ppr,pba,pth,ppt,
                pbs,ptr,ppm,psp,
                ple,pbi),
    pbi->error_message,
    pbi->error_message);



  // ====================================================================================
  // =                          Prepare interpolation mesh                              =
  // ====================================================================================

  /* Compute the 3D interpolation mesh if needed. Given a bispectrum b_l1_l2_l3,
  the mesh is basically its interpolation table conveniently binned in l space.
  Note that we cannot compute the 3D mesh at the last moment, like we do for the
  2D mesh, because it would create a race condition in fisher_compute_matrix(). 
  See comment in bispectra_mesh_interpolate() for more details. */

  if (pbi->interpolation_method == mesh_interpolation_3D) {

    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

      if (pbi->interpolate_me[index_bt]) {

        /* The array pbi->meshes contains the table required to interpolate the bispectra
        with the 2D and 3D mesh methods. Its first level, pbi->meshes[index_l1], is used
        only for the 2D mesh interpolation. For 3D interpolation, this level is redundant
        because a single 3D mesh is enough to interpolate the whole (l1,l2,l3) space.
        Therefore, for 3D interpolation we have only one entry: pbi->meshes[0]. */
        class_call (bispectra_mesh_create_for_all_probes(
                      ppr, psp, ple, pbi,
                      index_bt,
                      -1,
                      pbi->meshes[0]),
          pbi->error_message,
          pbi->error_message);

      }
    }
  }

  

  // =====================================================================================
  // =                                  Compute lensing                                  =
  // =====================================================================================

  /* Compute the three types of CMB bispectra: analytic, separable and non-separable */
  
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    if (pbi->lens_me[index_bt] && pbi->lens_me_brute_force[index_bt]) {

      class_call (bispectra_lensing (
                    ppr,pba,pth,ppt,
                    pbs,ptr,ppm,psp,
                    ple,pbi,
                    index_bt),
        pbi->error_message,
        pbi->error_message);

    }
  }


  /* Apart from pbi->bispectra, at this point all the arrays in the module have
  been completely filled. If the user requested to load the bispectra from disk,
  we can stop the execution of this module now without regrets. */
    
  if (ppr->load_bispectra) {

    printf_log_if (pbi->bispectra_verbose, 0, 
      " -> the intrinsic and non-separable bispectra will be read from disk\n");

    goto output_and_exit;

  }
  
  

  // =====================================================================================
  // =                              Store bispectra to disk                              =
  // =====================================================================================
  
  /* Save the bispectra to disk if requested. By default, we save only the non-separable
  bispectra because they take much more time to compute. */

  if (ppr->store_bispectra) {

    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
      if (pbi->bispectrum_type[index_bt] == non_separable_bispectrum)
        class_call (bispectra_store (
                      ppr,
                      pbi,
                      index_bt),
          pbi->error_message,
          pbi->error_message);

  }
  
  /* Check that we correctly filled the bispectra array (but only if there are no
  intrinsic bispectra left to be computed)*/
  if (pbi->n[intrinsic_bispectrum] < 1)
    class_warning (pbi->count_allocated_bispectra != pbi->count_memorised_bispectra,
      "there is a mismatch between allocated (%ld) and used (%ld) space!",
      pbi->count_allocated_bispectra, pbi->count_memorised_bispectra);



  output_and_exit:

  // =====================================================================================
  // =                                  Produce output                                   =
  // =====================================================================================

  /* Create output files containing the bispectra. If the intrinsic bispectrum is
  requested, we wait until the bispectra2.c module is executed. */
  
  if (!pbi->has_intrinsic) {

    class_call (bispectra_output (
                  ppr,pba,pth,ppt,
                  pbs,ptr,ppm,psp,
                  ple,pbi),
      pbi->error_message,
      pbi->error_message);
      
  }


  return _SUCCESS_;

}



/**
 * Evaluate the reduced bispectrum at a given (l1,l2,l3) configuration belonging
 * to SONG l-sampling.
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

int bispectra_at_node (
    struct bispectra * pbi,
    int index_bt,
    int index_l1, int index_l2, int index_l3,
    int X, int Y, int Z,
    double * bispectrum,           /**< Output: the bispectrum in (l1,l2,l3) */
    double * bispectrum_unlensed   /**< Output: the unlensed bispectrum in (l1,l2,l3) */
    )
{

  int l1 = pbi->l[index_l1];
  int l2 = pbi->l[index_l2];
  int l3 = pbi->l[index_l3];

  /* Debug: print arguments */
  // printf ("(%d^%s,%d^%s,%d^%s)\n", l1, pbi->bf_labels[X],
  //   l2, pbi->bf_labels[Y], l3, pbi->bf_labels[Z]);
  
#ifdef DEBUG
  class_test (!is_triangular_int(l1,l2,l3),
    pbi->error_message,
    "(l1=%d, l2=%d, l3=%d) is not a triangular configuration", l1, l2, l3);
#endif // DEBUG
  

  /* Find the ordering of (l1,l2,l3) */
  
  int index_l[4] = {0, index_l1, index_l2, index_l3};

  int order[4];

  class_call (ordering_int (index_l, order, pbi->error_message),
    pbi->error_message,
    pbi->error_message);


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

  *bispectrum = pbi->bispectra[index_bt]
                              [XYZ[order[3]]]
                              [XYZ[order[2]]]
                              [XYZ[order[1]]]
                              [index_l1_l2_l3];
  
  if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt])) {

    double correction = pbi->lensing_correction[index_bt]
                                               [XYZ[order[3]]]
                                               [XYZ[order[2]]]
                                               [XYZ[order[1]]]
                                               [index_l1_l2_l3];

    *bispectrum_unlensed = *bispectrum - correction;

  }
    

  /* The squeezed approximation in SONG are computed assuming that l3>=l2>=l1 */ 

  if (!pbi->is_symmetric[index_bt] && ((l2>l3) || (l1>l3) || (l1>l2))) {

      *bispectrum = 0;

      if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt]))
        *bispectrum_unlensed = 0;
  }
  
  return _SUCCESS_;
  
}



/**
 * Interpolate linearly the reduced bispectrum in l3 for a (l1,l2) pair belonging to
 * SONG l-sampling.
 */

int bispectra_at_l3_linear (
    struct transfers * ptr,
    struct bispectra * pbi,
    int index_bt,
    int index_l1, int index_l2, int l3,
    int X, int Y, int Z,
    int extrapolate,               /**< Input: if true, extrapolate the bispectrum in case l3 is outside the region
                                   computed by SONG for this (l1,l2) pair. If false, use the closest neighbour. */
    double * bispectrum,           /**< Output: the bispectrum in (l1,l2,l3) */
    double * bispectrum_unlensed   /**< Output: the unlensed bispectrum in (l1,l2,l3) */
    )
{

  int l1 = pbi->l[index_l1];
  int l2 = pbi->l[index_l2];

  class_test (!pbi->bispectra_available[index_bt],
    pbi->error_message,
    "bispectrum not available! maybe you should load it from disk?");

  class_test ((l3<abs(l1-l2)) || (l3>(l1+l2)),
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d) does not satisfy the triangular condition",
    l1, l2, l3);

  class_test (l3 < pbi->l[0],
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d), l3 is smaller than l_min=%d", l1, l2, l3, pbi->l[0]);

  class_test (l3 > pbi->l_max,
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d), l3 is larger than l_max=%d", l1, l2, l3, pbi->l_max);

  /* First node for this (l1,l2) pair */
  int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
  int l3_min = pbi->l[index_l3_min];
  
  /* Last node for this (l1,l2) pair */
  int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
  int l3_max = pbi->l[index_l3_max];


  // ====================================================================================
  // =                                   Special cases                                  =
  // ====================================================================================

  /* Special case A: if for this (l1,l2) there are no l3 nodes, then return zero
  and tell the user the sampling was insufficient to interpolate. With the current
  sampling scheme, this is never going to happen because there is at least one
  point for each (l1,l2) pair; eg. the slice l1=100, l2=2 has the point (100,2,100). */

  if (pbi->l_triangular_size[index_l1][index_l2] <= 0) {
        
    *bispectrum = 0;
    
    if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt]))
      *bispectrum_unlensed = 0;
    
    class_stop (pbi->error_message, "found l3_size=%d for index_l1=%d, index_l2=%d",
      pbi->l_triangular_size[index_l1][index_l2], index_l1, index_l2);
      
    return _SUCCESS_;    
 
  }

  /* Special case B: if for this (l1,l2) there is only one l3 node then return
  its value */
  
  else if (pbi->l_triangular_size[index_l1][index_l2] == 1) {
    
    int index_l3 = pbi->index_l_triangular_min[index_l1][index_l2];
    
    class_call (bispectra_at_node (
                  pbi,
                  index_bt,
                  index_l1, index_l2, index_l3,
                  X, Y, Z,
                  bispectrum,
                  bispectrum_unlensed),
      pbi->error_message,
      pbi->error_message);
      
    return _SUCCESS_;    
    
  }

  /* Special case C: if the requested l3 is smaller than the first node, extrapolate
  backward linearly using the slope of the first two nodes, unless extrapolate==_FALSE,
  in which case return the value at the first node. */

  else if (l3 < l3_min) {

    if (!extrapolate) {

      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_min,
                    X, Y, Z,
                    bispectrum,
                    bispectrum_unlensed),
        pbi->error_message,
        pbi->error_message);
    }
   
    else {
   
      int l3_first = l3_min;
      double b_first, b_first_unlensed;
   
      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_min,
                    X, Y, Z,
                    &b_first,
                    &b_first_unlensed),
        pbi->error_message,
        pbi->error_message);
   
      int l3_second = pbi->l[index_l3_min+1];
      double b_second, b_second_unlensed;
   
      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_min+1,
                    X, Y, Z,
                    &b_second,
                    &b_second_unlensed),
        pbi->error_message,
        pbi->error_message);
   
      double slope = (b_second-b_first)/(double)(l3_second-l3_first);
      *bispectrum = b_first - (l3_first-l3) * slope;
   
      if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt])) {
        double slope = (b_second_unlensed-b_first_unlensed)/(double)(l3_second-l3_first);
        *bispectrum_unlensed = b_first_unlensed - (l3_first-l3) * slope;
      }
    }

    return _SUCCESS_; 
    
  }

  /* Special case D: if the requested k3 is larger than the last node, extrapolate
  forward linearly using the slope of the first two nodes, unless extrapolate==_FALSE,
  in which case return the value at the last node. */
  
  else if (l3 > l3_max) {
    
    if (!extrapolate) {
    
      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_max,
                    X, Y, Z,
                    bispectrum,
                    bispectrum_unlensed),
        pbi->error_message,
        pbi->error_message);
    }
    
    else {
    
      int l3_last = l3_max;
      double b_last, b_last_unlensed;
    
      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_max,
                    X, Y, Z,
                    &b_last,
                    &b_last_unlensed),
        pbi->error_message,
        pbi->error_message);
    
      int l3_penultimate = pbi->l[index_l3_max-1];
      double b_penultimate, b_penultimate_unlensed;
    
      class_call (bispectra_at_node (
                    pbi,
                    index_bt,
                    index_l1, index_l2, index_l3_max-1,
                    X, Y, Z,
                    &b_penultimate,
                    &b_penultimate_unlensed),
        pbi->error_message,
        pbi->error_message);
    
      double slope = (b_last-b_penultimate)/(double)(l3_last-l3_penultimate);
      *bispectrum = b_last + (l3-l3_last) * slope;
    
      if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt])) {
        double slope = (b_last_unlensed-b_penultimate_unlensed)/(double)(l3_last-l3_penultimate);
        *bispectrum_unlensed = b_last_unlensed + (l3-l3_last) * slope;
      }
    }
  
    return _SUCCESS_; 
        
  }

  /* Special case E: if l3 is a node (ie. it belongs to pbi->l), then just return the
  bispectrum at the tabulated value */

  else if (ptr->index_l[l3] > 0) {
    
    class_call (bispectra_at_node (
                  pbi,
                  index_bt,
                  index_l1, index_l2, ptr->index_l[l3],
                  X, Y, Z,
                  bispectrum,
                  bispectrum_unlensed),
      pbi->error_message,
      pbi->error_message);
      
    return _SUCCESS_;
    
  }


  // ====================================================================================
  // =                                   Interpolation                                  =
  // ====================================================================================

  /* Index in pbi->l preceding l3 */
  int index_l3_left = ptr->index_l_left[l3];
  int l3_left = pbi->l[index_l3_left];

  /* Index in pbi->l following l3; this index must be a node because of special case D */
  int index_l3_right = index_l3_left + 1;
  class_test (index_l3_right > index_l3_max,
    pbi->error_message,
    "could not bracket l3");
  int l3_right = pbi->l[index_l3_right];
  
  /* Extract bispectrum in l3_left */
  double b_left, b_left_unlensed;
  
  class_call (bispectra_at_node (
                pbi,
                index_bt,
                index_l1, index_l2, index_l3_left,
                X, Y, Z,
                &b_left,
                &b_left_unlensed),
    pbi->error_message,
    pbi->error_message);

  /* Extract bispectrum in l3_right */
  double b_right, b_right_unlensed;
  
  class_call (bispectra_at_node (
                pbi,
                index_bt,
                index_l1, index_l2, index_l3_right,
                X, Y, Z,
                &b_right,
                &b_right_unlensed),
    pbi->error_message,
    pbi->error_message);

  /* Perform interpolation */
  double a = (l3_right - l3)/(double)(l3_right - l3_left);
  
  *bispectrum = a*b_left + (1-a)*b_right;

  if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt]))
    *bispectrum_unlensed = a*b_left_unlensed + (1-a)*b_right_unlensed;

  /* Debug: print interpolation quantities for a specific (l1,l2,l3) */
  // if ((l1==1500) && (l2==1200) && (l3 == 302)) {
  //   printf ("l3_left = %d(%d)\n", l3_left, index_l3_left);
  //   printf ("l3_right = %d(%d)\n", l3_right, index_l3_right);
  //   printf ("b_left = %g\n", b_left);
  //   printf ("b_right = %g\n", b_right);
  //   printf ("a = %g\n", a);
  // }

  return _SUCCESS_;

}



/**
 * Interpolate the reduced bispectrum in (l2,l3) with l1 belonging to the
 * SONG's l-sampling, using bilinear interpolation.
 *
 * Depending on the requested point, we adopt three different kinds of
 * interpolations:
 *
 * - For points that are well into the triangular borders (eg. l1=200,
 *   l2=200, l3=100), we interpolate along the l2 and l3 directions
 *   (RECTANGULAR interpolation).
 *
 * - For points that are close to the triangular border or right onto it
 *   (eg. l1=200, l2=300, l3=100), we use a linear interpolation along the
 *   l3 direction and a linear interpolation along the direction parallel
 *   to the border of the triangular region (TRIANGULAR interpolation).
 *
 * - For squeezed configurations where l3 is much smaller than the other
 *   two multipoles (eg. l1=200, l2=195, l3=3), we use a linear interpolation
 *   along the l3 direction and, for the l2 direction, either a closest 
 *   neighbour interpolation (extrapolate == 0) or linear extrapolation
 *   (extrapolate == 1). This mode is never used if l1<=l2<=l3, ie.
 *   in the computation of the Fisher matrix.
 *
 * Ideally, we would interpolate using the two closest nodes in (l2,l3)
 * space. Searching for these two nodes takes time however; our algorithm
 * relies on searching first on the horizontal & vertical lines, and then,
 * if nothing is found, on the diagonal lines. This makes sense because
 * the (l2,l3) domain is determined by horizontal & verical lines (the 
 * l_min and l_max limits on harmonic space) and by diagonal lines (the
 * triangular condition on (l1,l2,l3)), and it is also faster.
 */

int bispectra_at_l2l3_bilinear (
    struct transfers * ptr,
    struct bispectra * pbi,
    int index_bt,
    int index_l1, int l2, int l3,
    int X, int Y, int Z,
    int extrapolate,               /**< Input: if true, extrapolate linearly the bispectrum in the l3 direction
                                   when needed; if false, use the closest neighbour. */
    double * bispectrum,           /**< Output: the bispectrum in (l1,l2,l3) */
    double * bispectrum_unlensed  /**< Output: the unlensed bispectrum in (l1,l2,l3) */
    )
{

  int l1 = pbi->l[index_l1];

  class_test ((l3<abs(l1-l2)) || (l3>(l1+l2)),
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d) does not satisfy the triangular condition",
    l1, l2, l3);

  class_test (l2 < pbi->l[0],
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d), l2 is smaller than l_min=%d", l1, l2, l3, pbi->l[0]);

  class_test (l2 > pbi->l_max,
    pbi->error_message,
    "(l1,l2,l3)=(%d,%d,%d), l2 is larger than l_max=%d", l1, l2, l3, pbi->l_max);


  /* Before even considering interpolation, check whether l2 is a node. If this the
  case, then interpolate only along l3 */

  if (ptr->index_l[l2] > 0) {
    
    class_call (bispectra_at_l3_linear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, ptr->index_l[l2], l3,
                  X, Y, Z,
                  extrapolate,
                  bispectrum,
                  bispectrum_unlensed),
      pbi->error_message,
      pbi->error_message);
      
    return _SUCCESS_;

  }


  // ====================================================================================
  // =                              What strategy to use?                               =
  // ====================================================================================
  
  /* We define four types of interpolation strategies: rectangular, triangular,
  left node and right node. Refer to the documentation below for details. */
  enum {RECTANGULAR, TRIANGULAR, LEFT_NODE, RIGHT_NODE} interpolation;

  /* Index in pbi->l preceding l2 */
  int index_l2_left = ptr->index_l_left[l2];
  int l2_left = pbi->l[index_l2_left];

  /* Index in pbi->l following l2. This index cannot exceed the size of pbi->l because
  we made sure that l2 is strictly smaller than pbi->l_max */
  int index_l2_right = index_l2_left + 1;
  class_test (index_l2_right >= pbi->l_size,
    pbi->error_message,
    "could not bracket l2");
  int l2_right = pbi->l[index_l2_right];

  /* Determine triangular limits for l2_left and l2_right */
  int l2_min = MAX (abs(l1-l3), 2);
  int l2_max = MIN (l1+l3, pbi->l_max);

  /* Variables needed by the triangular method */
  int l3_left, l3_right;
  
  /* Debug: print l2 values and limits */
  // printf ("l1=%d,l2=%d,l3=%d,l2_left=%d,l2_right=%d,l2_min=%d,l2_max=%d\n",
  //   l1,l2,l3,l2_left,l2_right,l2_min,l2_max);
  // fflush (stdout);


  // -------------------------------------------------------------------------------
  // -                                 Rectangular?                                -
  // -------------------------------------------------------------------------------

  /* Are the left and right l2 nodes usable for interpolating in l3? A node is
  usable if its l2 value satisfy the triangular inequality with l1 and l3. If
  both the left and right nodes are usable, it means that (l2,l3) is far from
  the triangular limit and can be interpolated linearly along l2. The diagonal
  limits set by the triangular condition are problematic for this kind of
  interpolation, and are better dealt with the triangular interpolation, below. .*/

  short left_is_usable = (l2_left >= l2_min);
  short right_is_usable = (l2_right <= l2_max);

  if (left_is_usable && right_is_usable) {
    
    interpolation = RECTANGULAR;

  }
  

  // -------------------------------------------------------------------------------
  // -                                  Triangular?                                -
  // -------------------------------------------------------------------------------

  /* If interpolation along l2 is not possible, let's see if we can interpolate
  diagonally. The (l2,l3) dominion is a triangle filled with lines parallel to
  the triangular border. In the triangular method, we interpolate along these
  lines. This kind of interpolation is efficient where the rectangular interpolation
  fails (dealing with the diagonal limits set by the triangular condition), and it
  fails where the rectangular interpolation is ideal (dealing with the hard limit
  on l2 and l3 set by l_min=2 and l_max=pbi->l_max). The best scenario
  for the triangular interpolation is a very small l1, so that the (l2,l3) domain
  is mostly diagonal. The worst scenario is for l1~l_max, where the shape of the
  (l2,l3) domain is basically a right triangle. */

  else {

    /* Determine whether the requested point is in the upper or lower
    part of the (l2,l3) triangle */

    int offset_from_top = (l1+l2) - l3;
    int offset_from_bottom = l3 - abs(l1-l2);

    /* If the point is in the upper part of the triangle, interpolate its
    value along the direction parallel to the upper side of the triangle */
    if (offset_from_top < offset_from_bottom) {

      /* If l1 is equal to the upper limit in multipole space, then the whole upper
      part of the triangle is not accessible, so how can (l1,l2,l3) be there? */
      class_test (l1==pbi->l_max,
        pbi->error_message,
        "(l1=%d, l2=%d, l3=%d) cannot be in the upper triangle!", l1, l2, l3);

      l3_left = (l1+l2_left) - offset_from_top;
      l3_right = (l1+l2_right) - offset_from_top;
    }

    /* If the point is in the lower part of the triangle, interpolate its
    value along the direction parallel to the lower side of the triangle */
    else {
      l3_left = abs(l1-l2_left) + offset_from_bottom;
      l3_right = abs(l1-l2_right) + offset_from_bottom;
    }

    /* If the left node is very distant from l2, it is possible that
    it does not have an l3 value that intersects the diagonal passing
    through (l2,l3). In other words, the triangular condition cannot
    be met because the node is too far. In these (rare) cases, we assume
    for the function the value of the closest node to the right */
    if (!is_triangular_int (l1, l2_left, l3_left)) {

      interpolation = RIGHT_NODE;

    }

    /* It is now safe to use the triangular interpolation */
    else {
      
      interpolation = TRIANGULAR;
      
    }

    /* Check that the right node satisfies the triangular condition.
    Will we ever enter here? The right node should always satisfy
    the triangular condition, because as we increase l2 the triangle 
    never shrinks... */
    class_test (!is_triangular_int (l1, l2_right, l3_right),
      pbi->error_message,
      "l1=%d, l2_right=%d, l3_right=%d not triangular", l1, l2_right, l3_right);


    // -------------------------------------------------------------------------------
    // -                                 Squeezed?                                   -
    // -------------------------------------------------------------------------------

    /* The squeezed corner with l1~l2 is problematic. There, a small step in l3
    (say from 4 to 2) can result in abrupt changes in the bispectrum, especially
    for those bispectra peaked on squeezed configurations, like the local and
    CMB-lensing bispectra. This is not much of a problem for the interpolation along
    l3 (ie. the bispectra_at_l3_linear() function), because our l sampling has a logarithmic
    leg that ensures that the small values (2,3,4,5...) are densely sampled.

    The problem is for the l2 interpolation when using the triangular method. In
    fact, the diagonal method which interpolates along the diagonal border of
    the triangular region; but this is a direction along which the bispectrum varies
    quickly because it goes straight towards the ultra squeezed pole with l1=l2 and
    l3=2. If l1 and l2 are larger than ~ 100, the sampling along this diagonal
    direction is sparse, usually with a linear step of ~25. As a result, a mildy
    squeezed configuration (say l1=300, l2=250, l3=50) could be interpolated using
    another mildy squeezed node (say l1=300, l2=220, l3=70) and the ultra squeezed pole
    (l1=300, l2=300, l3=2). For a bispectrum peaked on squeezed configurations, the
    latter node would dominate and make the interpolated bispectrum in (l2,l3) much
    larger than what it really is.

    Note that using the rectangular interpolation to solve this issue is not possible,
    because, these squeezed configurations are at the limit of the triangular
    condition (eg. l1=l2=1000, l3=10).

    To recap, the triangular method interpolates in a given (l2,l3) point using nodes
    with different values of l3. This is a problem when the interpolated function
    varies a lot with l3, that is, for squeezed configurations where l3 is small and
    l1~l2.

    To solve this problem, we assign to those (l2,l3) points close to a squeezed
    node the value at the squeezed node with the same l3; that is, we assign
    b(l2,l3) = b(l2_left,l3) if the squeezed node is the left one, and
    b(l2,l3) = b(l2_right, l3) if the squeezed node is the right one. In this way
    we avoid using the fastest varying direction, l3. We lose in accuracy because
    we use only one point rather than two, but this loss is more than balanced
    by avoiding varying the l3 direction. */

    /* Determine whether we are dealing with a squeezed node */
    double squeezed_limit = 5;
    short left_is_squeezed = (l2_left/(double)l3_left > squeezed_limit);
    short right_is_squeezed = (l2_right/(double)l3_right > squeezed_limit);

    /* If one of the nodes is a squeezed configuration with small l3, we use it
    for interpolating the bispectrum rather than using the triangular method.
    An exception is when neither the left nor the right l2 values satisfy the
    triangular condition with l3, in which case the diagonal interpolation is
    the only option even if it is not very accurate */

    if (left_is_usable && left_is_squeezed)
      interpolation = LEFT_NODE;

    else if (right_is_usable && right_is_squeezed)
      interpolation = RIGHT_NODE;

  }



  // ====================================================================================
  // =                             Rectangular interpolation                            =
  // ====================================================================================

  /* If (l2,l3) is an internal point, interpolate along l2 */

  if (interpolation == RECTANGULAR) {

    /* If the l3-value of either the left or the right node is outside the triangular
    region, adjust it so that it fits. This should not happen because we made sure
    that we enter here only if both the left and right nodes are usable; we include
    this modifications for debug purposes. TODO: Ideally, if l3_left is changed, then
    also l3_right should be changed, because otherwise the line between l3_left and
    l3_right won't pass through l3. We should therefore define two l3 values rather
    than one. */

    int l3_rectangular = l3;

    int l3_left_min = MAX (abs(l2_left - l1), 2);
    int l3_left_max = MIN (l2_left+l1, pbi->l_max);

    int l3_right_min = MAX (abs(l2_right - l1), 2);
    int l3_right_max = MIN (l2_right+l1, pbi->l_max);

    l3_rectangular = MAX (l3_rectangular, MAX (l3_left_min, l3_right_min));
    l3_rectangular = MIN (l3_rectangular, MIN (l3_left_max, l3_right_max));
    
    /* Extract bispectrum in l2_left and l2_right */
    double b_left, b_right;
    double b_left_unlensed, b_right_unlensed;

    class_call (bispectra_at_l3_linear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, index_l2_left, l3_rectangular,
                  X, Y, Z,
                  extrapolate,
                  &b_left,
                  &b_left_unlensed),
      pbi->error_message,
      pbi->error_message);

    class_call (bispectra_at_l3_linear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, index_l2_right, l3_rectangular,
                  X, Y, Z,
                  extrapolate,
                  &b_right,
                  &b_right_unlensed),
      pbi->error_message,
      pbi->error_message);

    /* Interpolate along the l2 direction */
    double a = (l2_right - l2)/(double)(l2_right - l2_left);
    *bispectrum = a*b_left + (1-a)*b_right;

    if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt]))
      *bispectrum_unlensed = a*b_left_unlensed + (1-a)*b_right_unlensed;

  } // if (RECTANGULAR)


  // ====================================================================================
  // =                              Triangular interpolation                            =
  // ====================================================================================

  /* If (l2,l3) is close to the triangular limit and is not too squeezed,
  interpolate along the diagonal border of the triangular region. The 
  advantage of this interpolation is that any point can be interpolated,
  even if it is right on the triangular border or if has no usable nodes at 
  all. */

  else if (interpolation == TRIANGULAR) {

    /* If the l3-value of either the left or the right node is smaller than l_min or
    larger than l_max, adjust it so that it fits. This might happen because the
    diagonal interpolation is not compatible with lines parallel to the l2 and l3
    axes, such as those corresponding to the l_min and l_max limits. The worst
    situation is for very large l1, where the (l2,l3) domain is basically a right
    triangle. In these cases, we just interpolate using the closest node that is in
    the limits. TODO: Ideally, if l3_left is changed, then also l3_right should be
    changed, because otherwise the line between l3_left and l3_right won't pass
    through l3.  */

    l3_left = MAX (MIN (l3_left, pbi->l_max), 2);
    l3_right = MAX (MIN (l3_right, pbi->l_max), 2);


    // -------------------------------------------------------------------------------
    // -                            Bispectrum at the nodes                          -
    // -------------------------------------------------------------------------------

    double b_left, b_left_unlensed;

    class_call (bispectra_at_l3_linear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, index_l2_left, l3_left,
                  X, Y, Z,
                  extrapolate,
                  &b_left,
                  &b_left_unlensed),
      pbi->error_message,
      pbi->error_message);

    double b_right, b_right_unlensed;

    class_call (bispectra_at_l3_linear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, index_l2_right, l3_right,
                  X, Y, Z,
                  extrapolate,
                  &b_right,
                  &b_right_unlensed),
      pbi->error_message,
      pbi->error_message);


    // -------------------------------------------------------------------------------
    // -                                Interpolate                                  -
    // -------------------------------------------------------------------------------

    /* Distance squared between left and right nodes */
    double h_squared = (l2_right-l2_left)*(l2_right-l2_left) + (l3_right-l3_left)*(l3_right-l3_left);

    /* Distance squared between point and right node */
    double d_squared = (l2_right-l2)*(l2_right-l2) + (l3_right-l3)*(l3_right-l3);

    /* Proximity factor */
    double a = sqrt( d_squared/h_squared );

    /* Linear interpolation */
    *bispectrum = a*b_left + (1-a)*b_right; 

    if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt]))
      *bispectrum_unlensed = a*b_left_unlensed + (1-a)*b_right_unlensed;

    /* Debug: print l2 values and limits */
    // if (l1==10 && l2==15) {
    //
    //   printf ("l1=%d,l2=%d,l3=%d,l2_left=%d,l2_right=%d,l2_min=%d,l2_max=%d\n",
    //     l1,l2,l3,l2_left,l2_right,l2_min,l2_max);
    //
    //   printf ("\tl3_left = %d\n", l3_left);
    //   printf ("\tl3_right = %d\n", l3_right);
    //   printf ("\tb_left = %g\n", b_left);
    //   printf ("\tb_right = %g\n", b_right);
    //   printf ("\ta = %g\n", a);
    // }

  } // if (TRIANGULAR)


  // ====================================================================================
  // =                        Closest neighbour interpolation                           =
  // ====================================================================================

  /* Interpolate using either of the left or right node. The advantage of this
  approach is that it can be used even if a point cannot be bracketed by two valid
  nodes. */
  
  else if (interpolation == LEFT_NODE) {

    /* By default, we just use the previous node, ie. l2_left. If the user asked for
    extrapolation, we try to use backward extrapolation, that is, we use the two previous
    nodes (ie. the left node and the node at its left) and extrapolate forward to l2. We
    do so only if the node with the smallest l2 satisfies the triangular condition. */

    short backward_extrapolation = (
      (extrapolate) && /* user asked for extrapolation */
      ((index_l2_left-1) > 0) && /* left node is not the first one */
      (pbi->l[index_l2_left-1] >= l2_min)); /* node with the smallest l2 satisfies triangular condition */

    if (!backward_extrapolation) {

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_left, l3,
                    X, Y, Z,
                    extrapolate,
                    bispectrum,
                    bispectrum_unlensed),
        pbi->error_message,
        pbi->error_message);

    }

    else {

      /* First l2 node that is smaller than l2 */
      int l2_last = l2_left;
      double b_last, b_last_unlensed;

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_left, l3,
                    X, Y, Z,
                    extrapolate,
                    &b_last,
                    &b_last_unlensed),
        pbi->error_message,
        pbi->error_message);

      /* Second l2 node that is smaller than l2 */
      int l2_penultimate = pbi->l[index_l2_left-1];
      double b_penultimate, b_penultimate_unlensed;

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_left-1, l3,
                    X, Y, Z,
                    extrapolate,
                    &b_penultimate,
                    &b_penultimate_unlensed),
        pbi->error_message,
        pbi->error_message);

      double slope = (b_last-b_penultimate)/(double)(l2_last-l2_penultimate);
      *bispectrum = b_last + (l2-l2_last) * slope;

      if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt])) {
        double slope = (b_last_unlensed-b_penultimate_unlensed)/(double)(l2_last-l2_penultimate);
        *bispectrum_unlensed = b_last_unlensed + (l2-l2_last) * slope;
      }

    } // if backward extrapolation

  } // if (LEFT_NODE)
  
  else if (interpolation == RIGHT_NODE) {

    /* By default, we just use the next node, ie. l2_right. If the user asked for
    extrapolation, we try to use forward extrapolation, that is, we use the next two
    nodes (ie. the right node and the node at its right) and extrapolate back to l2.
    We do so only if the node with the largest l2 satisfies the triangular condition. */

    short forward_extrapolation =
      (extrapolate) && /* user asked for forward extrapolation */
      ((index_l2_right+1) < pbi->l_size) && /* right node is not the last one */
      (pbi->l[index_l2_right+1] <= l2_max); /* node with the largest l2 satisfies triangular condition */

    if (!forward_extrapolation) {

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_right, l3,
                    X, Y, Z,
                    extrapolate,
                    bispectrum,
                    bispectrum_unlensed),
        pbi->error_message,
        pbi->error_message);

    }

    else {

      /* First l2 node that is larger than l2 */
      int l2_first = l2_right;
      double b_first, b_first_unlensed;

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_right, l3,
                    X, Y, Z,
                    extrapolate,
                    &b_first,
                    &b_first_unlensed),
        pbi->error_message,
        pbi->error_message);

      /* Second l2 node that is larger than l2 */
      int l2_second = pbi->l[index_l2_right+1];
      double b_second, b_second_unlensed;

      class_call (bispectra_at_l3_linear (
                    ptr,
                    pbi,
                    index_bt,
                    index_l1, index_l2_right+1, l3,
                    X, Y, Z,
                    extrapolate,
                    &b_second,
                    &b_second_unlensed),
        pbi->error_message,
        pbi->error_message);

      double slope = (b_second-b_first)/(double)(l2_second-l2_first);
      *bispectrum = b_first - (l2_first-l2) * slope;

      if ((bispectrum_unlensed != NULL) && (pbi->lens_me[index_bt])) {
        double slope = (b_second_unlensed-b_first_unlensed)/(double)(l2_second-l2_first);
        *bispectrum_unlensed = b_first_unlensed - (l2_first-l2) * slope;
      }

    } // if forward extrapolation

  } // if (RIGHT_NODE)


  return _SUCCESS_;

}






/**
 * Interpolate the reduced bispectrum in (l2,l3) for a l1 multipole belonging to
 * SONG l-sampling.
 *
 * The algorithm used for the interpolation is determined via the
 * pbi->interpolation_method flag.
 *
 * So far, the following interpolation methods are supported:
 *
 * - bilinear
 * - mesh_interpolation_2D
 * - mesh_interpolation_3D
 */

int bispectra_at_l2l3 (
    struct precision * ppr,
    struct transfers * ptr,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi,
    int index_bt,
    int index_l1,
    int l1, int l2, int l3,
    int X, int Y, int Z,
    double * bispectrum,         /**< Output: the bispectrum in (l1,l2,l3) */
    double * bispectrum_unlensed /**< Output: the unlensed bispectrum in (l1,l2,l3), enabled only for bilinear interpolation */
    )
{

#ifdef DEBUG
  class_test ((index_l1>=0) && l1!=pbi->l[index_l1],
    pbi->error_message,
    "inconsistent input: (l1=%d) != (pbi->l[%d]=%d)", l1, index_l1, pbi->l[index_l1]);
    
  class_test (bispectrum_unlensed != NULL &&
              pbi->lens_me[index_bt] &&
              pbi->interpolation_method != bilinear_interpolation,
    pbi->error_message,
    "to interpolate the unlensed bispectrum use bilinear interpolation (TODO: extend to mesh)");
#endif // DEBUG


  if (pbi->interpolation_method == bilinear_interpolation) {

    class_call (bispectra_at_l2l3_bilinear (
                  ptr,
                  pbi,
                  index_bt,
                  index_l1, l2, l3,
                  X, Y, Z,
                  _TRUE_,
                  bispectrum,
                  bispectrum_unlensed),
      pbi->error_message,
      pbi->error_message);

  }
  
  else if (pbi->interpolation_method == mesh_interpolation_2D) {

    class_test (!pbi->interpolate_me[index_bt],
      pbi->error_message,
      "cannot interpolate %s bispectrum - mesh wasn't computed",
      pbi->bt_labels[index_bt]);

    class_call (bispectra_mesh_interpolate(
                  ppr, psp, ple, pbi,
                  index_bt,
                  index_l1,
                  l1, l2, l3,
                  X, Y, Z,
                  pbi->meshes[index_l1][index_bt][X][Y][Z][0],
                  pbi->meshes[index_l1][index_bt][X][Y][Z][1],
                  bispectrum),
      pbi->error_message,
      pbi->error_message);
             
  }

  else if (pbi->interpolation_method == mesh_interpolation_3D) {

    class_test (!pbi->interpolate_me[index_bt],
      pbi->error_message,
      "cannot interpolate %s bispectrum - mesh wasn't computed",
      pbi->bt_labels[index_bt]);

    /* For 3D interpolation, the index_l1 parameter is ignored because
    l1 is not restricted to the l-sampling in pbi->l. Only l1 is used.
    This is also why we fix index_l1=0 in pbi->meshes[index_l1]: the
    whole l1 level is collapsed onto 0. */

    class_call (bispectra_mesh_interpolate(
                  ppr, psp, ple, pbi,
                  index_bt,
                  -1,
                  l1, l2, l3,
                  X, Y, Z,
                  pbi->meshes[0][index_bt][X][Y][Z][0],
                  pbi->meshes[0][index_bt][X][Y][Z][1],
                  bispectrum),
      pbi->error_message,
      pbi->error_message);
             
  }
          
  return _SUCCESS_;
  
}




/**
 * Save the bispectra in pbi->bispectra to disk for a given type.
 * 
 * See the documentation in perturbations2.h (\ref StorageFiles) for more
 * details.
 */

int bispectra_store (
      struct precision * ppr,
      struct bispectra * pbi,
      int index_bt
      )
{

  class_test (!ppr->store_bispectra,
    pbi->error_message,
    "shouldn't be here");

  /* Print some debug */
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * writing bispectra to disk for index_bt=%d on '%s'\n",
    index_bt, pbi->storage_paths[index_bt]);

  class_test (!pbi->bispectra_allocated[index_bt],
    pbi->error_message,
    "cannot store bispectra if they are not even allocated");

  class_test (!pbi->bispectra_available[index_bt],
    pbi->error_message,
    "cannot store bispectra if they are not available");

  /* Open file for writing */
  class_open (pbi->storage_files[index_bt],
    pbi->storage_paths[index_bt],
    "a+b", pbi->error_message);

  /* Write to file */
  for (int X = 0; X < pbi->bf_size; ++X)
    for (int Y = 0; Y < pbi->bf_size; ++Y)
      for (int Z = 0; Z < pbi->bf_size; ++Z)
        fwrite(
              pbi->bispectra[index_bt][X][Y][Z],
              sizeof(double),
              pbi->n_independent_configurations,
              pbi->storage_files[index_bt]
              );

  /* Close file */
  fclose(pbi->storage_files[index_bt]);
  
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

  if (pbi->has_bispectra) {

    free(pbi->bt_labels);
    free(pbi->bf_labels);
    free(pbi->bfff_labels);

    free(pbi->l);
    free(psp->pk);
    free(psp->pk_pt);

    for(int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {

      for (int index_l2=0; index_l2 <= index_l1; ++index_l2)
        free (pbi->l3[index_l1][index_l2]);

      free(pbi->l_triangular_size[index_l1]);
      free(pbi->index_l_triangular_min[index_l1]);
      free(pbi->index_l_triangular_max[index_l1]);
      free(pbi->l3_size[index_l1]);
      free(pbi->l3[index_l1]);
          
    } // for(index_l1)

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
    free (pbi->bispectra_allocated);
    free (pbi->bispectra_available);
    
    /* Arrays specific to the primordial models */
    if ((pbi->has_local_model) || (pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
        
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
    } // local model

    if ((pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
  
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
    } // equilateral and orthogonal models  
    
    free (pbi->delta_k);
    
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->cls[index_ct]);
    free (pbi->cls);
    for (int index_ct=0; index_ct < psp->ct_size; ++index_ct)
      free (pbi->d_lsq_cls[index_ct]);
    free (pbi->d_lsq_cls);
    if (ppr->extend_lensed_cls) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_d_lsq_cls[index_lt]);
      free (pbi->lensed_d_lsq_cls);
    }
    
    if (ppr->extend_lensed_cls) {
      for (int index_lt=0; index_lt < ple->lt_size; ++index_lt)
        free (pbi->lensed_cls[index_lt]);
      free (pbi->lensed_cls);
    }
      
    /* Free file arrays */
    if (ppr->store_bispectra || ppr->load_bispectra) {
    
      for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt)
        free (pbi->storage_paths[index_bt]);
    
      free (pbi->storage_files);
      free (pbi->storage_paths);
    
    }

    /* Free interpolation meshes */
    
    if (pbi->interpolation_method == mesh_interpolation_2D ||
        pbi->interpolation_method == mesh_interpolation_3D) {

      for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {

        /* There is only one 3D mesh for all l1 values */
        if (pbi->interpolation_method==mesh_interpolation_3D && index_l1!=0)
          continue;

        class_call (bispectra_mesh_free(pbi, &pbi->meshes[index_l1]),
          pbi->error_message,
          pbi->error_message);

      }
    }

    free (pbi->meshes);

  } // if(has_bispectra)
  
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


  if (pbi->lens_me[index_bt]) {
    
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

  pbi->bispectra_allocated[index_bt] = _FALSE_;
  pbi->bispectra_available[index_bt] = _FALSE_;

  return _SUCCESS_;
 
}





/**
 * Allocate the type level of the bispectrum array (pbi->bispectra).
 */

int bispectra_allocate_type_level (
      struct bispectra * pbi,
      int index_bt
      )
{

  /* Allocate memory only if needed */
  if (pbi->bispectra_allocated[index_bt])
    return _SUCCESS_;

  /* Allocate pbi->bispectra[index_bt] */
  class_alloc (pbi->bispectra[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

  for (int X = 0; X < pbi->bf_size; ++X) {
    
    class_alloc (pbi->bispectra[index_bt][X], pbi->bf_size*sizeof(double **), pbi->error_message);

    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      
      class_alloc (pbi->bispectra[index_bt][X][Y], pbi->bf_size*sizeof(double *), pbi->error_message);
      
      for (int Z = 0; Z < pbi->bf_size; ++Z) 
        class_calloc (pbi->bispectra[index_bt][X][Y][Z],
          pbi->n_independent_configurations,
          sizeof(double),
          pbi->error_message);

    }
  }

  #pragma omp atomic
  pbi->count_allocated_bispectra += pbi->n_probes*pbi->n_independent_configurations;


  /* Do the same for the lensing correction */
  
  if (pbi->has_lensed_bispectra) {

    class_alloc (pbi->lensing_correction, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
    for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {

      if (pbi->lens_me[index_bt])
        continue;

      class_alloc (pbi->lensing_correction[index_bt], pbi->bf_size*sizeof(double ***), pbi->error_message);

      for (int X = 0; X < pbi->bf_size; ++X) {
      
        class_alloc (pbi->lensing_correction[index_bt][X], pbi->bf_size*sizeof(double **), pbi->error_message);

        for (int Y = 0; Y < pbi->bf_size; ++Y) {
        
          class_alloc (pbi->lensing_correction[index_bt][X][Y], pbi->bf_size*sizeof(double *), pbi->error_message);
        
          for (int Z = 0; Z < pbi->bf_size; ++Z)
            class_calloc (pbi->lensing_correction[index_bt][X][Y][Z],
              pbi->n_independent_configurations,
              sizeof(double),
              pbi->error_message);

        }
      }
    }
  }

  /* We succesfully allocated the bt level of pbi->bispectra */
  pbi->bispectra_allocated[index_bt] = _TRUE_;

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
 *  -# Store P(k) (power spectrum of phi) in psp->pk for all the k-values in the transfer
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
  

  // ====================================================================================
  // =                              Count bispectra fields                              =
  // ====================================================================================

  /* Find out which kind of bispectra to compute and assign them indices and labels
  Generate indices for the probes (T for temperature, E for E-mode polarisation,
  R for Rayleigh...) that we will use to build the bispectra and the Fisher matrix
  elements. */
  
  int index_bf = 0;
  
  pbi->has_bispectra_t = _FALSE_;
  pbi->has_bispectra_e = _FALSE_;
  pbi->has_bispectra_b = _FALSE_;

  class_calloc (pbi->bf_labels,
    _MAX_NUM_FIELDS_*_MAX_LENGTH_LABEL_,
    sizeof(char),
    pbi->error_message);
    
  class_calloc (pbi->bfff_labels,
    _MAX_NUM_FIELDS_*_MAX_NUM_FIELDS_*_MAX_NUM_FIELDS_*_MAX_LENGTH_LABEL_,
    sizeof(char),
    pbi->error_message);
    
  if (ppt->has_bi_cmb_temperature) {
    pbi->has_bispectra_t = _TRUE_;
    strcpy (pbi->bf_labels[index_bf], "t");
    pbi->field_parity[index_bf] = _EVEN_;
    pbi->field_spin[index_bf] = 0;
    pbi->index_bf_t = index_bf++;
  }
  
  if (ppt->has_bi_cmb_polarization) {
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
    "cannot compute the bispectrum for more than %d fields (e.g. T and E)",
    pbi->bf_size);

  class_test (pbi->bf_size < 1,
    pbi->error_message,
    "no probes requested");

  /* Create labels for the full bispectra */
  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {
        for (int i=0; i < _MAX_LENGTH_LABEL_; ++i)
          pbi->bfff_labels[X][Y][Z][i] = '\0';
        sprintf (pbi->bfff_labels[X][Y][Z], "%s%s%s",
          pbi->bf_labels[X], pbi->bf_labels[Y], pbi->bf_labels[Z]);
      }
    }
  }


  /* Associate to each field X=T,E,... its transfer function, which was computed in
  the transfer.c module, the C_l correlation <X phi> with the lensing potential,
  the C_l correlation <X zeta> with the curvature perturbation and, to each possible
  pair of fields (TT,EE,TE,...), their power spectra, which were computed in the
  spectra.c module. */
  for (int X = 0; X < pbi->bf_size; ++X) {
    
    if ((pbi->has_bispectra_t) && (X == pbi->index_bf_t)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_t;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_tp;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_tt;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_tt;
      if (ppt->has_cl_cmb_zeta)
        pbi->index_ct_of_zeta_bf[X] = psp->index_ct_tz;
      if (ppr->extend_lensed_cls)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_tt;
    }

    if ((pbi->has_bispectra_e) && (X == pbi->index_bf_e)) {
      pbi->index_tt_of_bf[X] = ptr->index_tt_e;
      pbi->index_ct_of_phi_bf[X] = psp->index_ct_ep;
      pbi->index_ct_of_t_bf[X] = psp->index_ct_te;
      pbi->index_ct_of_bf_bf[X][X] = psp->index_ct_ee;
      if (ppt->has_cl_cmb_zeta)
        pbi->index_ct_of_zeta_bf[X] = psp->index_ct_ez;
      if (ppr->extend_lensed_cls)
        pbi->index_lt_of_bf_bf[X][X] = ple->index_lt_ee;
    }

    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      if (((pbi->has_bispectra_t) && (X == pbi->index_bf_t))
       && ((pbi->has_bispectra_e) && (Y == pbi->index_bf_e))) {

        pbi->index_ct_of_bf_bf[X][Y] = pbi->index_ct_of_bf_bf[Y][X] = psp->index_ct_te;
        if (ppr->extend_lensed_cls)
          pbi->index_lt_of_bf_bf[X][Y] = pbi->index_lt_of_bf_bf[Y][X] = ple->index_lt_te;
      }
    }
  }


  // ====================================================================================
  // =                               Count bispectra types                              =
  // ====================================================================================
  
  /* - Initialise variables */
  
  pbi->n[separable_bispectrum] = 0;
  pbi->n[non_separable_bispectrum] = 0;
  pbi->n[analytical_bispectrum] = 0;
  pbi->n[intrinsic_bispectrum] = 0;

  class_calloc (pbi->bt_labels,
    _MAX_NUM_BISPECTRA_*_MAX_LENGTH_LABEL_,
    sizeof(char),
    pbi->error_message);

  for (int index_bt=0; index_bt < _MAX_NUM_BISPECTRA_; ++index_bt) {
    pbi->is_symmetric[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = (pbi->has_lensed_bispectra?_TRUE_:_FALSE_);
    pbi->lens_me_brute_force[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _TRUE_;
  }
  
  int index_bt = 0;

  // -------------------------------------------------------------------------------
  // -                               Separable bispectra                           -
  // -------------------------------------------------------------------------------
  
  if (pbi->has_local_model) {
    pbi->index_bt_local = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_equilateral_model) {
    pbi->index_bt_equilateral = index_bt;
    strcpy (pbi->bt_labels[index_bt], "equilateral");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }

  if (pbi->has_orthogonal_model) {
    pbi->index_bt_orthogonal = index_bt;
    strcpy (pbi->bt_labels[index_bt], "orthogonal");
    pbi->bispectrum_type[index_bt] = separable_bispectrum;
    pbi->n[separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }


  // -------------------------------------------------------------------------------
  // -                             Non-separable bispectra                         -
  // -------------------------------------------------------------------------------

  if (pbi->has_galileon_model) {

    /* Bispectrum induced by pi_dot*pi_grad^2 */
    pbi->index_bt_galileon_gradient = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_grad");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;

    /* Bispectrum induced by pi_dot^3 */
    pbi->index_bt_galileon_time = index_bt;
    strcpy (pbi->bt_labels[index_bt], "galileon_time");
    pbi->bispectrum_type[index_bt] = non_separable_bispectrum;
    pbi->n[non_separable_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }


  // -------------------------------------------------------------------------------
  // -                               Analytic bispectra                            -
  // -------------------------------------------------------------------------------

  if (pbi->has_cmb_lensing) {
    pbi->index_bt_cmb_lensing = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    /* Uncomment the following two lines to make SONG lens and interpolate the
    CMB lensing bispectrum. This is useful to debug the interpolation and lensing
    of the bispectrum, because the CMB lensing bispectrum can be computed and lensed
    analytically. */
    // pbi->interpolate_me[index_bt] = _TRUE_;
    // pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }
  
  if (pbi->has_quadratic_correction) {
    pbi->index_bt_quadratic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "quadratic");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  if (pbi->has_local_squeezed) {
    pbi->index_bt_local_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "local_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->is_symmetric[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_intrinsic_squeezed) {
    pbi->index_bt_intrinsic_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->is_symmetric[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_cmb_lensing_squeezed) {
    pbi->index_bt_cmb_lensing_squeezed = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->is_symmetric[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  /* The kernel for the squeezed CMB-lensing bispectrum is needed to compute
  the lensing contribution to the variance. In the final Fisher matrix, the kernel
  will be multiplied by C_l^{X\phi} to give the actual squeezed bispectrum (see
  eq. 5.20 of http://uk.arxiv.org/abs/1101.2234); it will therefore show up as
  CMB-lensing_sqz rather than kernel_sqz. */
  if (pbi->has_cmb_lensing_kernel) {
    pbi->index_bt_cmb_lensing_kernel = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cmb_lensing_sqz");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->is_symmetric[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }
  
  if (pbi->has_cosine_shape) {
    pbi->index_bt_cosine = index_bt;
    strcpy (pbi->bt_labels[index_bt], "cosine");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me[index_bt] = _FALSE_;
    pbi->interpolate_me[index_bt] = _FALSE_;
    index_bt++;
  }

  if (pbi->has_test_bispectrum) {
    pbi->index_bt_test = index_bt;
    strcpy (pbi->bt_labels[index_bt], "lensing_test");
    pbi->bispectrum_type[index_bt] = analytical_bispectrum;
    pbi->n[analytical_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    pbi->interpolate_me[index_bt] = _TRUE_;
    index_bt++;
  }


  // -------------------------------------------------------------------------------
  // -                             Second-order bispectra                          -
  // -------------------------------------------------------------------------------

  if (pbi->has_intrinsic) {
    pbi->index_bt_intrinsic = index_bt;
    strcpy (pbi->bt_labels[index_bt], "intrinsic");
    pbi->bispectrum_type[index_bt] = intrinsic_bispectrum;
    pbi->n[intrinsic_bispectrum]++;
    pbi->has_reduced_bispectrum[index_bt] = _TRUE_;
    pbi->lens_me_brute_force[index_bt] = _TRUE_;
    index_bt++;
  }

  pbi->bt_size = index_bt;

  class_test (pbi->bt_size > _MAX_NUM_BISPECTRA_,
   "exceeded maximum number of allowed bispectra, increase _MAX_NUM_BISPECTRA_ in common.h",
   pbi->error_message);

#ifndef WITH_SONG2
   class_test (pbi->has_intrinsic,
     pbi->error_message,
     "cannot compute the intrinsic bispectrum without SONG support\n");
#endif // WITH_SONG2

  /* Are the Wigner 3j-symbols needed to compute the requested bispectra? */
  pbi->need_3j_symbols = ((pbi->has_bispectra_e) &&
    ((pbi->has_quadratic_correction) ||
     (pbi->has_cmb_lensing) ||
     (pbi->has_cmb_lensing_squeezed) ||
     (pbi->has_cmb_lensing_kernel)));


  
  
  // ====================================================================================
  // =                               Determine l-sampling                               =
  // ====================================================================================
  
  /* We compute the bispectrum on a mesh where l1>=l2>=l3, with l3 determined by
  the triangular condition. All the multipoles are drawn from pbi->l, which is a
  copy of ptr->l, the sampling determined in transfer_get_l_list(). */

  pbi->l_size = pbs->l_size;
  class_alloc (pbi->l, pbi->l_size*sizeof(int), pbi->error_message);
  for(int index_l=0; index_l<pbi->l_size; ++index_l)
    pbi->l[index_l] = pbs->l[index_l];

  /* Maximum value in pbi->l */
  pbi->l_min = pbi->l[0];
  pbi->l_max = pbi->l[pbi->l_size-1];
  pbi->full_l_size = pbi->l_max - pbi->l_min + 1;
  
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

        /* When the triangular condition is not compatible with index_l3<=index_l2, then
        index_l3_max < index_l3_min+1 will be either zero or negative */ 
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

    } // for(index_l2)
  } // for(index_l1)

  /* Each bispectrum will be computed for the following number of configurations */
  pbi->n_independent_configurations = index_l1_l2_l3;
  

  /* Inform the user on how much her machine will have to suffer */
  if (pbi->bispectra_verbose > 0) {
    printf_log (" -> we shall compute %dx%d=%d bispectr%s for %ld configurations of (l1,l2,l3)\n",
      pbi->bt_size, pbi->n_probes, pbi->bt_size*pbi->n_probes,
      ((pbi->bt_size*pbi->n_probes)!=1?"a":"um"), pbi->n_independent_configurations);
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

  /* Compute the trapezoidal integration measure in k, needed for the integration
  of both the separable and non-separable bispectra. In the latter case, we shall
  define a separate measure for k3, because the grid in k3 is not fixed but it
  depends on k1 and k2 via the triangular condition. */
  class_alloc (pbi->delta_k, k_tr_size * sizeof(double), pbi->error_message);
  
  /* Fill pbi->delta_k */
  class_call (trapezoidal_weights (
                k_tr,
                k_tr_size,
                k_tr[0],
                k_tr[k_tr_size-1],
                _FALSE_,
                pbi->delta_k,
                NULL,
                NULL,
                pbi->error_message),
    pbi->error_message,
    pbi->error_message);
  
  
  
  // ====================================================================================
  // =                            Allocate memory for bispectra                         =
  // ====================================================================================

  /* Allocate and initialize the logical arrays keeping track of the state of
  pbi->bispectra */
  class_calloc (pbi->bispectra_allocated, pbi->bt_size, sizeof(short), pbi->error_message);
  class_calloc (pbi->bispectra_available, pbi->bt_size, sizeof(short), pbi->error_message);

  /* Keep track of memory usage (debug only) */
  pbi->count_allocated_bispectra = 0;
  pbi->count_memorised_bispectra = 0;  

  class_alloc (pbi->bispectra, pbi->bt_size*sizeof(double ****), pbi->error_message);
  
  for (int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
    
    class_call (bispectra_allocate_type_level(pbi, index_bt),
      pbi->error_message,
      pbi->error_message);

  }
  
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * allocated ~ %.3g MB (%ld doubles) for the bispectra array\n",
    pbi->count_allocated_bispectra*sizeof(double)/1e6, pbi->count_allocated_bispectra);


  
  // ====================================================================================
  // =                               Create storage files                               =
  // ====================================================================================
  
  /* Create the files to store the bispectra in */
  if (ppr->store_bispectra || ppr->load_bispectra) {

    /* We are going to store the bispectra in n=bt_size files, one for each requested type of bispectrum */
    class_alloc (pbi->storage_files, pbi->bt_size*sizeof(FILE *), pbi->error_message);
    class_alloc (pbi->storage_paths, pbi->bt_size*sizeof(char *), pbi->error_message);
  
    for(int index_bt=0; index_bt<pbi->bt_size; ++index_bt) {
      
      /* Include the name of the bispectrum in its file */
      class_alloc (pbi->storage_paths[index_bt], _FILENAMESIZE_*sizeof(char), pbi->error_message);
      sprintf (pbi->storage_paths[index_bt], "%s/bispectra_%s.dat", pbi->storage_dir, pbi->bt_labels[index_bt]);
      
    }

    if (ppr->store_bispectra)
      printf_log_if (pbi->bispectra_verbose, 1, 
        "     * will create %d files for the bispectra\n", pbi->bt_size);
      
  }
  
  

  // ====================================================================================
  // =                                 Interpolate P(k)                                 =
  // ====================================================================================
  
  /* Store the primordial power spectrum for the curvature potential Phi inside the
  spectra structure for faster access */
  class_call (spectra_primordial_power_spectrum(
                pba,
                ppt,
                ptr,
                ppm,
                psp),
    psp->error_message,
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



  // ====================================================================================
  // =                                Mesh interpolation                                = 
  // ====================================================================================

  /* Compute the basic parameters for the mesh interpolation. As described in
  bispectra.h, the mesh interpolation is a way to linearly interpolate the
  bispectrum (or any other function) on a mesh rather than on a grid. This 
  method allows to naturally solve the issues related to the triangular (l1,l2,l3) domain
  of the bispectrum. SONG supports both 2D and 3D interpolation using the mesh
  method. */

  if ((pbi->interpolation_method == mesh_interpolation_2D) ||
      (pbi->interpolation_method == mesh_interpolation_3D)) {


    // -------------------------------------------------------------------------------
    // -                                Parameters                                   -
    // -------------------------------------------------------------------------------

    /* We shall use two meshes for each bispectrum: a fine mesh for the small l
    and a coarse mesh for the large l. See bispectra.h for a description of the
    parameters of the two meshes. */
    
    /* Fine mesh parameters. We use a linking length of the same order of magnitude
    of the smallest distance between two points. */
    pbi->mesh_link_lengths[0] = 2 * (pbi->l[1] - pbi->l[0]);
    pbi->mesh_group_lengths[0] = 0.1 * (pbi->l[1] - pbi->l[0]);
    pbi->mesh_soft_coeffs[0] = 0.5;  

    /* Coarse mesh parameters. We use a linking length of the same order of magnitude
    of the largest distance between two points, ppr->l_linstep (the linear step in l) */
    int l_linstep = ppr->l_linstep;

    /* Adjust l_linstep to accomodate for an all evel l-grid. If ppr->l_linstep=1, all
    l are used and there is effectively no interpolation. However, if we use an all even
    l-grid, the actual step between one multipole and the other is doubled, as the odd
    l are skipped. */
    if ((l_linstep==1) && (ppr->compute_only_even_ls))
      l_linstep = 2;

    pbi->mesh_link_lengths[1] = 0.5/sqrt(2) * l_linstep;
    pbi->mesh_group_lengths[1] = 0.1 * l_linstep;
    pbi->mesh_soft_coeffs[1] = 0.5;

    for (int index_mesh=0; index_mesh < 2; ++index_mesh)
      class_test (pbi->mesh_link_lengths[index_mesh] <= pbi->mesh_group_lengths[index_mesh],
        pbi->error_message,
        "the linking length must be larger than the grouping length.");

    /* We restrict the (l1,l2,l3) domain for the interpolation of the bispectra to
    those configurations with l1<=l2<=l3, unless lensing of the bispectrum is
    requested, in which case all configurations are needed */
    pbi->mesh_restrict_l1l2l3 = !pbi->has_lensed_bispectra;


    // -----------------------------------------------------------------------------
    // -                              Turnover point                               -
    // -----------------------------------------------------------------------------

    /* The turnover point is the l-value below which the fine grid is used, and
    above which the coarse grid is used. Setting the turnover to the largerst
    l in the sampling will make SONG always use the fine grid. The result would be
    more or less the same, but the execution time will increase considerably.
    If instead you set the turnver to be smaller than the smallest l, SONG will
    always use the coarse grid. The speed would increase but the result would
    basically ignore the contribution from small l. */

    /* Never use the fine grid if the large step is used since the beginning */
    if ((pbi->l[1]-pbi->l[0]) >= l_linstep) {
      pbi->mesh_l_turnover = pbi->l[0];
    }

    /* Never use the coarse grid if the linstep (ie. the largest possible step)
    is never used. The MAX() is needed because the last point is fixed and does
    not depend on the grid spacing. */
    else if (MAX(pbi->l[pbi->l_size-1]-pbi->l[pbi->l_size-2], pbi->l[pbi->l_size-2]-pbi->l[pbi->l_size-3]) < l_linstep) {
      pbi->mesh_l_turnover = pbi->l[pbi->l_size-1] + 1;
    }

    /* Compute the turnover point in the l-grid where the linear step begins */
    else {
      int index_l = 0;
      while ((index_l < pbi->l_size-1) && ((pbi->l[index_l+1] - pbi->l[index_l]) < l_linstep))
        index_l++;
      pbi->mesh_l_turnover = pbi->l[index_l-1];
    }

    /* Print some info on the interpolation meshes */
    printf_log_if (pbi->bispectra_verbose, 1, 
      "     * mesh_interpolation: l_turnover=%d, n_bins=[%d,%d], linking lengths=[%g,%g], grouping lengths=[%g,%g]\n",
      pbi->mesh_l_turnover,
      (int)ceil(pbi->mesh_l_turnover/ (pbi->mesh_link_lengths[0]*(1+pbi->mesh_soft_coeffs[0]))),
      (int)ceil(pbi->l[pbi->l_size-1] / (pbi->mesh_link_lengths[1]*(1+pbi->mesh_soft_coeffs[1]))),
      pbi->mesh_link_lengths[0], pbi->mesh_link_lengths[1],
      pbi->mesh_group_lengths[0], pbi->mesh_group_lengths[1]);

    
    // -------------------------------------------------------------------------------
    // -                            Allocate meshes                                -
    // -------------------------------------------------------------------------------
    
    /* We allocate a fine grid and a coarse mesh for each bispectrum type & probe.
    If we are dealing with 2D interpolation, we also allocate a mesh for every
    l1 value in our sampling. For 3D interpolation, this is not required because
    a single 3D mesh is enough to interpolate the bispectrum for all values of l1. */
    
    int l1_size;
    
    if (pbi->interpolation_method == mesh_interpolation_2D)
      l1_size = pbi->l_size;
    
    else if (pbi->interpolation_method == mesh_interpolation_3D)
      l1_size = 1;

    class_alloc (pbi->meshes,
      l1_size * sizeof (*pbi->meshes),
      pbi->error_message);

    for (int index_l1=0; index_l1 < l1_size; ++index_l1) {

      class_call (bispectra_mesh_allocate(
                    ppr,psp,ple,pbi,
                    &pbi->meshes[index_l1]),
        pbi->error_message,
        pbi->error_message);
      
    }

  } // if(mesh_interpolation)


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
  if (ppr->extend_lensed_cls) {

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
    if (ppr->extend_lensed_cls) {

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

    if ((ppt->has_cl_cmb_lensing_potential) && (ppt->has_scalars)) {

      pbi->lmax_lensing_corrT = 300;
      pbi->lmax_lensing_corrE = 300;

      if ((l > pbi->lmax_lensing_corrT) && (pbi->has_bispectra_t))
        pbi->cls[psp->index_ct_tp][l-2] = 0;

      if ((l > pbi->lmax_lensing_corrE) && (pbi->has_bispectra_e))
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
    //   /* Debug: print the cls */
    //   // printf ("%d %g %g %g\n", l,
    //   //   pbi->cls[psp->index_ct_tt][l-2],
    //   //   pbi->cls[psp->index_ct_tz][l-2],
    //   //   pbi->d_lsq_cls[psp->index_ct_tt][l-2]);
    
  } // for loop on the l's

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

  } // loop on index_ct
  
  /* Do the same for the lensed C_l's */
  if (ppr->extend_lensed_cls) {
  
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

    } // loop on index_lt
  } // if lensing
  
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
  pbi->count_memorised_bispectra = 0;


  // ====================================================================================
  // =                               Separable bispectra                                =
  // ====================================================================================

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
                  psp,
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



  // ====================================================================================
  // =                               Analytic bispectra                                 =
  // ====================================================================================  
  
  /* Compute the bispectra obtained from simple analytical formulas, such as the
  lensing one */
  
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
  
  
	/* Debug: print some bispectra configurations */
  // for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
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
  //       } // for(index_l3)
  //     } // for(index_l2)
  //   } // for(index_l1)
  // } // for(index_bt)
  

  /* If we are loading the bispectra from disk, nothing else needs to be done */

  if (ppr->load_bispectra) {
        
    return _SUCCESS_;

  }
  
  
  // ====================================================================================
  // =                             Non-separable bispectra                              =
  // ====================================================================================
  
  if (pbi->n[non_separable_bispectrum] > 0) {
  
    struct bispectra_workspace_non_separable * pwb_nonsep;
    class_alloc (pwb_nonsep, sizeof(struct bispectra_workspace_non_separable), pbi->error_message);
  
    /* Compute the non-separable bispectra */
    class_call (bispectra_non_separable_init(
                  ppr, pba, pth, ppt,
                  pbs, ptr, ppm, psp,
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
  if (pbi->bispectra_verbose > 1 && pbi->n[intrinsic_bispectrum] < 1)
    printf(" -> memorised ~ %.3g MB (%ld doubles) in the bispectra array\n",
      pbi->count_memorised_bispectra*sizeof(double)/1e6, pbi->count_memorised_bispectra);
  
  
  
  // ====================================================================================
  // =                          Check bispectra against nan's                           =
  // ====================================================================================

  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

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
   
                class_warning (isnan(bispectrum),
                  "b(%d,%d,%d) = %g for bispectrum %s_%s",
                  pbi->l[index_l1], pbi->l[index_l2], pbi->l[index_l3], bispectrum,
                  pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

              } // for(index_l3)
            } // for(index_l2)
          } // for(index_l1)
        } // for(X)
      } // for(Y)
    } // for(Z)
  } // for(index_bt)
  
  return _SUCCESS_;

}



/**
 * Produce output files for the bispectrum b^XYZ(l1,l2,l3).
 *
 * The output files will be generated only for the multipole values
 * specified by the user with the (l1_out, l2_out, l3_out) parameters.
 * 
 * Three types of files will be created:
 *
 * - The l3 text files containing the bispectra for fixed (l1,l2) as a 
 *   function of l3.
 *
 * - The l2l3 text files containing the bispectra for fixed l1 as a 
 *   function of (l2,l3).
 *
 * - The l1l2l3 binary file containing all bispectra configurations
 *   for all bispectra probes; this file can be read and plotted
 *   using the scripts by Thomas Tram in the python folder.
 *
 * For the l2 and l2l3 outputs, a separate file is generated for each
 * XYZ probe (TTT,TTE,...). On the contrary, the l1l2l3 output will
 * contain all probes.
 * 
 * The binary files include accessory data (such as the l-sampling and 
 * the first-order C_l) and a human-readable ASCII header explaining how
 * to access their content.
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
  
  if (pbi->bispectra_verbose > 0)
    printf (" -> computing output files for the bispectrum\n");
  
  /* Load the bispectra from disk if needed */

  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
    class_call (bispectra_load (
                  ppr,
                  pbi,
                  index_bt),
      pbi->error_message,
      pbi->error_message);


  // ====================================================================================
  // =                            Output in (l2,l3) and l3                              =
  // ====================================================================================

  class_call (bispectra_output_l2l3 (
                ppr,pba,pth,ppt,
                pbs,ptr,ppm,psp,
                ple,pbi),
    pbi->error_message,
    pbi->error_message);

  

  // ====================================================================================
  // =                             Output in (l1,l2,l3)                                 =
  // ====================================================================================

  if (pbi->output_binary_bispectra) {

    class_call (bispectra_output_l1l2l3 (
                  ppr,pba,pth,ppt,
                  pbs,ptr,ppm,psp,
                  ple,pbi),
      pbi->error_message,
      pbi->error_message);

  }
    

  return _SUCCESS_;
  
}



/**
 * Create two kinds of output text files: one for the bispectra as a function
 * of (l2,l3) for the output values in l1_out and tau_out, and one for the
 * bispectra as a function of k3 for the output values in (l1_out,l2_out).
 *
 * A text file will be generated for each probe (TTT, EEE, TTE...) for each
 * kind of output.
 * 
 * The l2l3 output files are suitable to be plotted with gnuplot using the
 * splot function. To produce contour plots and, in general, nicer plots,
 * use the option 'set pm3d'.
 * 
 * See documentation of bispectra_output() for further details.
 */

int bispectra_output_l2l3 (
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

  for (int index_l_out=0; index_l_out < ppr->l_out_size; ++index_l_out) {

    for (int X = 0; X < pbi->bf_size; ++X) {

      for (int Y = 0; Y < pbi->bf_size; ++Y) {

        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          /* Multipole value and index of this output file */
          int l1 = ppr->l1_out[index_l_out];
          int l2 = ppr->l1_out[index_l_out];
          int index_l1 = ppr->index_l1_out[index_l_out];
          int index_l2 = ppr->index_l2_out[index_l_out];

          /* Build filenames */
          int index_probe = X*pbi->bf_size*pbi->bf_size + Y*pbi->bf_size + Z;
          sprintf (ppr->paths_bispectra_l3[index_l_out][index_probe], "%sbispectra_l3_L%03d_%s.txt",
            ppr->paths_bispectra_l3[index_l_out][index_probe], index_l_out, pbi->bfff_labels[X][Y][Z]);
          sprintf (ppr->paths_bispectra_l2l3[index_l_out][index_probe], "%sbispectra_l2l3_L%03d_%s.txt",
            ppr->paths_bispectra_l2l3[index_l_out][index_probe], index_l_out, pbi->bfff_labels[X][Y][Z]);

          /* Open files */
          FILE * file_1D = ppr->files_bispectra_l3[index_l_out][index_probe];
          FILE * file_2D = ppr->files_bispectra_l2l3[index_l_out][index_probe];
          class_open(file_1D, ppr->paths_bispectra_l3[index_l_out][index_probe], "w", pbi->error_message);
          class_open(file_2D, ppr->paths_bispectra_l2l3[index_l_out][index_probe], "w", pbi->error_message);
          

          // -------------------------------------------------------------------------------
          // -                               Print information                             -
          // -------------------------------------------------------------------------------
          
          char line[1024];
          
          /* Write the information header of the 1D and 2D files */
          if (l2 > 0)
            sprintf (line, "CMB reduced bispectra b_l1_l2_l3 tabulated as a function of l3 and bispectrum type for a fixed (l1,l2) pair.");
          else
            sprintf (line, "CMB reduced bispectra b_l1_l_l tabulated as a function of l for a fixed l1 value.");
          fprintf (file_1D, "%s%s\n", _COMMENT_, line);
          sprintf (line, "CMB reduced bispectra b_l1_l2_l3 tabulated as a function of (l2,l3) and bispectrum type for a fixed l1 value.");
          fprintf (file_2D, "%s%s\n", _COMMENT_, line);
          sprintf (line, "This file was generated by SONG %s (%s) on %s.", _SONG_VERSION_, _SONG_URL_, ppr->date);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%s\n", _COMMENT_);
          if (pbi->has_lensed_bispectra) {
            sprintf (line, "The suffix _u denotes unlensed bispectra.");
            fprintf_2way (1, file_1D, 0, file_2D, 0, "%s%s\n", _COMMENT_, line);
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
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%sInformation on the output l:\n", _COMMENT_);
          fprintf_2way (1, file_1D, 0, file_2D, 0, "%sl1 = %d, index_l1 = %d/%d\n", _COMMENT_, l1, index_l1, pbi->l_size-1);
          if (l2 > 0) fprintf (file_1D, "%sl2 = %d, index_l2 = %d/%d\n", _COMMENT_, l2, index_l2, pbi->l_size-1);


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
                
                class_call (bispectra_at_node (
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
                            &normalisation),
                pbi->error_message,
                pbi->error_message);

              class_warning (fabs(normalisation) < _MINUSCULE_,
                "normalisation=%g is small; beware of inf", normalisation);

              double normalisation_positive = 0;

              class_call (bispectra_normalisation_positive (
                            ppr, psp, ple, pbi,
                            l1, l2, l3,
                            X, Y, Z,
                            &normalisation_positive),
                pbi->error_message,
                pbi->error_message);

              class_warning (fabs(normalisation_positive) < _MINUSCULE_,
                "normalisation_positive=%g is small; beware of inf", normalisation_positive);


              // -------------------------------------------------------------------------------
              // -                                   Build row                                 -
              // -------------------------------------------------------------------------------

              /* Arrays containing all the information on the columns to be printed, labels included */
              char (*label)[_MAX_LENGTH_LABEL_] = calloc (_MAX_NUM_COLUMNS_*_MAX_LENGTH_LABEL_, sizeof(char));
              double * value = calloc (_MAX_NUM_COLUMNS_, sizeof(double));
              short * condition = calloc (_MAX_NUM_COLUMNS_, sizeof(short));
  
              /* Initialise column arrays */
              for (int i=0; i < _MAX_NUM_COLUMNS_; ++i)
                condition[i] = _TRUE_;

              /* Initialise column counter  */
              int i = -1;

              /* Multipole l2 (won't be printed on the 1D file) */
              strcpy (label[++i], "l2");
              value[i] = l2;
              
              /* Multipole l3 */
              if (l2 > 0)
                strcpy (label[++i], "l3");
              else
                strcpy (label[++i], "l");
              value[i] = l3;
              
              /* All bispectra types */
              for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

                sprintf (label[++i], "%s", pbi->bt_labels[index_bt]);
                value[i] = bispectrum[index_bt];

                if (pbi->lens_me[index_bt]) {
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

              /* Include a blank line at the beginning of each new l2; this allows to use the
              'set pm3d' option in gnuplot */
              if (index_l3 == index_l3_max)
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
              
              else if ((l2 < 0) && (l2 == l3)) {
                
                /* Write row with labels and append information on l2 to the header */
                int n_columns_1D = 1;
                if (n_rows_1D++ == 0) {
                  fprintf (file_1D, "%swill print b(%d,l,l) as a function of l\n", _COMMENT_, l1);
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

              free (label);
              free (value);
              free (condition);

            } // for l3
          } // for l2

          /* Close the files */
          fclose (file_1D);
          fclose (file_2D);



          // ====================================================================================
          // =                               Debug interpolation                                =
          // ====================================================================================

          int output_interpolated_bispectrum = _FALSE_;
          int extrapolate = _TRUE_;

          /* Output a 2D file where (l2,l3) can take any allowed value, using interpolation */
          
          if (output_interpolated_bispectrum) {
          
            FILE * file_2D_interpolated;
            char * file_path_2D_interpolated;

            class_call (replace_string (
                          ppr->paths_bispectra_l2l3[index_l_out][index_probe],
                          ".txt",
                          "_interpolated.txt",
                          &file_path_2D_interpolated,
                          pbi->error_message),
              pbi->error_message,
              pbi->error_message);

            class_open(file_2D_interpolated, file_path_2D_interpolated, "w", pbi->error_message);
            free (file_path_2D_interpolated);
          
            for (int l2=2; l2 <= pbi->l_max; ++l2) {

              int l3_min = MAX(abs(l1-l2),2);
              int l3_max = MIN((l1+l2),pbi->l_max);

              for (int l3=l3_min; l3 <= l3_max; ++l3) {

                // -------------------------------------------------------------------------------
                // -                              Extract bispectra                              -
                // -------------------------------------------------------------------------------
              
                double bispectrum[pbi->bt_size];
                double bispectrum_unlensed[pbi->bt_size];
  
                for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
                
                  class_call (bispectra_at_l2l3 (
                                ppr,
                                ptr,
                                psp,
                                ple,
                                pbi,
                                index_bt,
                                index_l1,
                                l1, l2, l3,
                                X, Y, Z,
                                &bispectrum[index_bt],
                                &bispectrum_unlensed[index_bt]),
                    pbi->error_message,
                    pbi->error_message);

                }

              
                // -------------------------------------------------------------------------------
                // -                                 Normalisation                               -
                // -------------------------------------------------------------------------------

                double normalisation = 0;

                class_call (bispectra_normalisation (
                              ppr, psp, ple, pbi,
                              l1, l2, l3,
                              X, Y, Z,
                              &normalisation),
                  pbi->error_message,
                  pbi->error_message);

                class_warning (fabs(normalisation) < _MINUSCULE_,
                  "normalisation=%g is small; beware of inf", normalisation);

                double normalisation_positive = 0;

                class_call (bispectra_normalisation_positive (
                              ppr, psp, ple, pbi,
                              l1, l2, l3,
                              X, Y, Z,
                              &normalisation_positive),
                  pbi->error_message,
                  pbi->error_message);

                class_warning (fabs(normalisation_positive) < _MINUSCULE_,
                  "normalisation_positive=%g is small; beware of inf", normalisation_positive);


                // -------------------------------------------------------------------------------
                // -                                   Build row                                 -
                // -------------------------------------------------------------------------------

                /* Arrays containing all the information on the columns to be printed, labels included */
                char (*label)[_MAX_LENGTH_LABEL_] = calloc (_MAX_NUM_COLUMNS_*_MAX_LENGTH_LABEL_, sizeof(char));
                double * value = calloc (_MAX_NUM_COLUMNS_, sizeof(double));
                short * condition = calloc (_MAX_NUM_COLUMNS_, sizeof(short));
  
                /* Initialise column counter  */
                int i = -1;

                /* Multipole l2 */
                strcpy (label[++i], "l2");
                value[i] = l2;
              
                /* Multipole l3 */
                strcpy (label[++i], "l3");
                value[i] = l3;
              
                /* Bispectra */
                for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

                  sprintf (label[++i], "%s", pbi->bt_labels[index_bt]);
                  value[i] = bispectrum[index_bt];

                  if (pbi->lens_me[index_bt]) {
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
                // -                       Print row to interpolation file                       -
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
                    fprintf (file_2D_interpolated, format_label, label[i], n_columns_2D++);
                  fprintf (file_2D_interpolated, "\n");
                }

                /* Write row with data to file */
                for (int i=0; i < n_max_columns; ++i)
                  fprintf (file_2D_interpolated, format_value, value[i]);
                fprintf (file_2D_interpolated, "\n");

                /* Include a blank line at the beginning of each new l2; this allows to use the
                'set pm3d' option in gnuplot */
                if (l3 == l3_max)
                  fprintf (file_2D_interpolated, "\n");
              
              
                free (label);
                free (value);
                free (condition);
              
              } // for l3
            } // for l2

            /* Close the file */
            fclose (file_2D_interpolated);

          } // if output_interpolated_bispectrum

        } // Z
      } // Y
    } // X


  } // for l_out
  
  return _SUCCESS_;

}



/**
 * Output to binary file the bispectra for all (l1,l2,l3) configurations.
 *
 * The binary file will also contain all probes (TTT, EEE, TTE...). The
 * file can can be read and plotted using the scripts by Thomas Tram in
 * the python folder.
 * 
 * See documentation of perturb2_output() for further details.
 */

int bispectra_output_l1l2l3 (
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
                &(ppr->file_bispectra_l1l2l3),
                ppr->path_bispectra_l1l2l3,
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

  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1)
    for (int index_l2=0; index_l2 <= index_l1; ++index_l2)
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

  sprintf (desc, "l3 array: l3[index_l1][index_l2] with index_l1 < pbi->l_size, index_l2 <= index_l1");
  sprintf (name, "pbi->l3");
  int index_l3_block = file->n_blocks;

  for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1)
    for (int index_l2=0; index_l2 <= index_l1; ++index_l2)
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

  sprintf (desc, "number of (l1,l2,l3) configurations in each bispectrum (=%ld)", pbi->n_independent_configurations);
  sprintf (name, "pbi->n_independent_configurations");
  class_call (binary_append_long_int (file, &pbi->n_independent_configurations, 1, desc, name),
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

  sprintf (desc, "unlensed C_l: cls[index_ct][index_l] with index_ct < psp->ct_size, index_l < pbi->full_l_size");
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

  sprintf (desc, "unlensed C_l logarithmic derivative: d_lsq_cls[index_ct][index_l] with index_ct < psp->ct_size, index_l < pbi->full_l_size");
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

  if (ppr->extend_lensed_cls) {

    sprintf (desc, "lensed C_l: cls[index_ct][index_l] with index_ct < psp->ct_size, index_l < pbi->full_l_size");
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

    sprintf (desc, "lensed C_l logarithmic derivative: d_lsq_cls[index_ct][index_l] with index_ct < psp->ct_size, index_l < pbi->full_l_size");
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

    class_test (!pbi->bispectra_available[index_bt],
      pbi->error_message,
      "cannot output bispectrum %d (%s): not available",
      index_bt, pbi->bt_labels[index_bt]);

    sprintf (desc, "%s CMB bispectrum for all values of (l1,l2,l3) and for all fields (X,Y,Z);\
 idexed as [X][Y][Z][index_l1][index_l2][index_l3] with XYZ < pbi->bf_size, index_l1 < pbi->l_size,\
 idex_l2 <= index_l1, index_l3 < pbi->l3_size[index_l1][index_l2]", pbi->bt_labels[index_bt]);
    sprintf (name, "pbi->bispectra[index_bt=%d]", index_bt);
    int index_bispectrum = file->n_blocks;

    for (int X = 0; X < pbi->bf_size; ++X) {
      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          for (int index_l1=0; index_l1 < pbi->l_size; ++index_l1) {        
            for (int index_l2=0; index_l2 <= index_l1; ++index_l2) {

              /* Determine the limits for l3, which come from the triangular inequality |l1-l2| <= l3 <= l1+l2 */
              int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
              int index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]);
              int l3_size = index_l3_max - index_l3_min + 1;

              if (l3_size > 0) {

                /* Index in the bispectrum array with the first l3 value for the current (l1,l2) pair.
                We can do this because all l3 bispectrum values are arranged continuously in memory.  */
                long int index_start = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3_min];

                class_call (binary_add_block (
                              file,
                              &(pbi->bispectra[index_bt][X][Y][Z][index_start]),
                              l3_size,
                              sizeof (double),
                              desc,
                              "double",
                              name,
                              index_bispectrum),
                  file->error_message,
                  pbi->error_message);
               }
                  
            } // for l2
          } // for l1

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

  if (pth->reio_parametrization != reio_none &&
     (pbi->has_bispectra_e || pbi->has_bispectra_b))
    class_warning (*r_min > pba->conformal_age-pth->tau_reio,
      "the r-sampling might be inadequate to compute reionisation");

  /* Check that the r-grid is strictly ascending */
  for (int index_r=0; index_r<(*r_size-1); ++index_r)
    class_test ((*r_grid)[index_r] >= (*r_grid)[index_r+1],
      pbi->error_message,
      "the r grid should be stricty ascending");


  /* Compute the measure for the trapezoidal integration over r */
  class_alloc (*delta_r, *r_size * sizeof(double), pbi->error_message);

  class_call (trapezoidal_weights (
                *r_grid,
                *r_size,
                (*r_grid)[0],
                (*r_grid)[*r_size-1],
                _FALSE_,
                *delta_r,
                NULL,
                NULL,
                pbi->error_message),
    pbi->error_message,
    pbi->error_message);

  
  
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
      struct spectra * psp,
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
  if ((pbi->has_local_model) ||
      (pbi->has_equilateral_model) ||
      (pbi->has_orthogonal_model)) {

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

    if (abort)
      return _FAILURE_;
    
  } // if all models



  // ---------------------------------------------------------
  // -          Specific to equilateral & orthogonal         -
  // ---------------------------------------------------------

  if ((pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {

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

    if (abort)
      return _FAILURE_;
    
  } // if equilateral || orthogonal

  return _SUCCESS_;
  
}




int bispectra_separable_filter_functions (
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct spectra * psp,
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
      if ((pbi->has_local_model) ||
          (pbi->has_orthogonal_model) ||
          (pbi->has_equilateral_model)) {

        for (int index_k=0; index_k < k_size; ++index_k) {
        
          double pk = psp->pk[index_k];
          double tr = transfer[index_k];

          pwb->alpha_integrand[thread][index_k] = tr;
          pwb->beta_integrand[thread][index_k] = tr * pk;
          
          /* Print out the transfer function for a given field */
          // if (ptr->l[index_l]==284)
          //   if ((ppt->has_cl_cmb_polarization) && (pbi->index_tt_of_bf[index_bf] == ptr->index_tt_e))
          //     fprintf (stderr, "%12g %12g\n", ptr->q[index_k], tr);
          
        }

      } // all models
    
      // *** Equilateral and orthogonal models ***
      if ((pbi->has_equilateral_model) ||
          (pbi->has_orthogonal_model)) {
     
        /* Here we basically copy eqs. 15-18 in Creminelli et al. 2006, keeping out the 2/pi factor and the
        Bessel function, and multiplying the rest by k (we use the dimensional power spectrum and we already
        factored out a k^2 factor) */
        for (int index_k=0; index_k < k_size; ++index_k) {      

          double pk_one_third = pow(psp->pk[index_k], 1/3.);
          double pk_two_thirds = pk_one_third*pk_one_third;
          double tr = transfer[index_k];

          pwb->gamma_integrand[thread][index_k] = tr * pk_one_third;
          pwb->delta_integrand[thread][index_k] = tr * pk_two_thirds;

        }
     
      } // equilateral and orthogonal model





      // =========================================================
      // =         Convolve with the Bessel function             =
      // =========================================================
      
      for (int index_r = 0; index_r < pwb->r_size; ++index_r) {
      
        double r = pwb->r[index_r];
      
        printf_log_if (pbi->bispectra_verbose, 4, 
          "       \\ r=%g, index_r=%d\n", r, index_r);
        
        
        // *** All models ***
        if ((pbi->has_local_model) ||
            (pbi->has_equilateral_model) ||
            (pbi->has_orthogonal_model)) {
                      
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
                
        } // local model
        
        // /* Some debug */
        // if (index_r == 70)
        //   fprintf (stderr, "%10d %17.7g %17.7g\n",
        //   pbi->l[index_l], pbi->alpha[index_bf][index_l][index_r], pbi->beta[index_bf][index_l][index_r]);

        
        // *** Equilateral and orthogonal models ***
        if ((pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {

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
      
        } // equilateral and orthogonal models
        
        #pragma omp flush(abort)
      
      } // for(index_r)
  
    } // for(index_l)
    
  } if (abort) return _FAILURE_; /* end of parallel region */

    
  /* Output the filter functions */
  // int index_l = 60;
  // 
  // if ((pbi->has_bispectra_t) && (index_bf == pbi->index_bf_t) ) {
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
      struct spectra * psp,
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
    shared (ppt,pbs,ptr,ppm,psp,pbi,pwb,abort)             \
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

            if ((pbi->has_local_model) && (index_bt == pbi->index_bt_local)) {
  
              /* The primordial bispectrum for the local model has only two extra contributions due to the symmetrization
              P(k1)*P(k2) + P(k1)*P(k3) + P(k2)*P(k3). It is easy to see that each contribution has the same value, but
              we sum them nonetheless to be consistent. */
              integrand = 2 * fnl_R * r*r * (
                  alpha[X][index_l1][index_r] * beta[Y][index_l2][index_r]  * beta[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * beta[Y][index_l2][index_r]  * alpha[Z][index_l3][index_r]
                + beta[X][index_l1][index_r]  * alpha[Y][index_l2][index_r] * beta[Z][index_l3][index_r]
                );
            
            } // local model

            // -------------------------------------------------------------------------
            // -                           Equilateral model                           -
            // -------------------------------------------------------------------------
            
            else if ((pbi->has_equilateral_model) && (index_bt == pbi->index_bt_equilateral)) {
  
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

            } // equilateral model

            // -------------------------------------------------------------------------------
            // -                              Orthogonal model                               -
            // -------------------------------------------------------------------------------
            
            else if ((pbi->has_orthogonal_model) && (index_bt == pbi->index_bt_orthogonal)) {
  
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

            } // orthogonal model

            
            /* Increment the estimate of the integral */
            integral += integrand * pwb->delta_r[index_r];

            /* Debug - Print integrand as a function of r */
            // if (index_bt == pbi->index_bt_local)
            //   if ((X==pbi->index_bf_t) && (Y==pbi->index_bf_t) && (Z==pbi->index_bf_t))
            //     if ((pbi->l[index_l1] == 2) && (pbi->l[index_l2] == 2) && (pbi->l[index_l3] == 2))
            //       fprintf (stderr, "%15.7g %15.7g\n", r, integrand);

  
          } // for(index_r)
  

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3) */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] = integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. In KSW2005, this factor
          was split between the alpha and beta integrals, but from the numerical point of view it is preferable
          to include it at the end of the computation (as in eq. 17 of Fergusson & Shellard 2007). */
          pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3] *= pow(2./_PI_,3);
  
          /* Some debug */
          // if ((index_l1==index_l2) && (index_l2==index_l3))
          //   fprintf (stderr, "%10d %17.7g\n", pbi->l[index_l1], pbi->bispectra[index_bt][index_l1][index_l2][index_l3-index_l3_min]);
    
          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_bispectra++;

        } // for(index_l3)
      } // for(index_l2)
      
      #pragma omp flush(abort)
  
    } // for(index_l1)
  
  } if (abort) return _FAILURE_;  // parallel region
  
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
      struct spectra * psp,
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
                psp,
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
                  psp,
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
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
  
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
                        psp,
                        pbi,
                        index_bt,
                        X,
                        Y,
                        Z,
                        pwb),
            pbi->error_message,
            pbi->error_message);
  
        } // for(X)
      } // for(Y)
    } // for(Z)
    
    /* The separable bispectra are now ready */
    pbi->bispectra_available[index_bt] = _TRUE_;

  } // for(index_bt)
  
  
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

    if ((pbi->has_local_model) || (pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
      free (pwb->alpha_integrand[thread]);
      free (pwb->beta_integrand[thread]);
    }

    if ((pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
      free (pwb->gamma_integrand[thread]);
      free (pwb->delta_integrand[thread]);
    }
  }

  if ((pbi->has_local_model) || (pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
    free (pwb->alpha_integrand);
    free (pwb->beta_integrand);
  }

  if ((pbi->has_equilateral_model) || (pbi->has_orthogonal_model)) {
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

    if ((pbi->has_cmb_lensing) && (index_bt == pbi->index_bt_cmb_lensing))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_bispectrum;

    else if ((pbi->has_cmb_lensing_squeezed) && (index_bt == pbi->index_bt_cmb_lensing_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_bispectrum;

    else if ((pbi->has_cmb_lensing_kernel) && (index_bt == pbi->index_bt_cmb_lensing_kernel))
      pbi->bispectrum_function[index_bt] = bispectra_cmb_lensing_squeezed_kernel;

    else if ((pbi->has_local_squeezed) && (index_bt == pbi->index_bt_local_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_local_squeezed_bispectrum;
    
    else if ((pbi->has_intrinsic_squeezed) && (index_bt == pbi->index_bt_intrinsic_squeezed))
      pbi->bispectrum_function[index_bt] = bispectra_intrinsic_squeezed_bispectrum;
     
    else if ((pbi->has_quadratic_correction) && (index_bt == pbi->index_bt_quadratic))
      pbi->bispectrum_function[index_bt] = bispectra_quadratic_correction;
    
    else if ((pbi->has_test_bispectrum) && (index_bt == pbi->index_bt_test))
      pbi->bispectrum_function[index_bt] = bispectra_test_bispectrum;
    
    else if ((pbi->has_cosine_shape) && (index_bt == pbi->index_bt_cosine))
      pbi->bispectrum_function[index_bt] = bispectra_cosine_bispectrum;

  }
  

  // ===================================================================================
  // =                                   Main loop                                     =
  // ===================================================================================
    
  /* We parallelize the outer loop over 'l1'. */
  #pragma omp parallel for shared (ppt,pbs,ptr,ppm,psp,pbi,abort) private (thread)
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

        if (pbi->need_3j_symbols) {          

          class_call_parallel (threej_ratio_M (l2, l1, l3, 2, &threej_ratio_20m2, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l3, l1, l2, 2, &threej_ratio_m220, pbi->error_message),
            pbi->error_message, pbi->error_message);

          class_call_parallel (threej_ratio_M (l1, l2, l3, 2, &threej_ratio_0m22, pbi->error_message),
            pbi->error_message, pbi->error_message);

        } // 3j computation

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

                /* Compute the bispectrum using the function associated to index_bt. If the
                bispectrum is marked for lensing this function computes the lensed bispectrum,
                unless pbi->lens_me_brute_force is true, in which case the bispectrum will be
                lensed later in the bispectra_lensing() function. */
                
                short lens_me = pbi->lens_me[index_bt] && !pbi->lens_me_brute_force[index_bt];

                class_call_parallel ((*pbi->bispectrum_function[index_bt]) (
                              ppr, psp, ple, pbi,
                              l1, l2, l3,
                              X, Y, Z,
                              lens_me,
                              threej_ratio_20m2,
                              threej_ratio_m220,
                              threej_ratio_0m22,
                              &(pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3])),
                  pbi->error_message,
                  pbi->error_message);
    

                if (lens_me) {
                  
                  /* If we are inside here, then we just stored the lensed bispectrum in
                  pbi->bispectra; now we compute the unlensed bispectrum as well, and
                  store the difference with the lensed one in pbi->lensing_correction */

                  double bispectrum_lensed = pbi->bispectra[index_bt][X][Y][Z][index_l1_l2_l3];
                  double bispectrum_unlensed;
                  
                  class_call_parallel ((*pbi->bispectrum_function[index_bt]) (
                                ppr, psp, ple, pbi,
                                l1, l2, l3,
                                X, Y, Z,
                                _FALSE_, /* give me the unlensed bispectrum */
                                threej_ratio_20m2,
                                threej_ratio_m220,
                                threej_ratio_0m22,
                                &bispectrum_unlensed),
                    pbi->error_message,
                    pbi->error_message);

                  pbi->lensing_correction[index_bt][X][Y][Z][index_l1_l2_l3] = bispectrum_lensed - bispectrum_unlensed;

                }
    
                /* Update the counter */
                #pragma omp atomic
                pbi->count_memorised_bispectra++;

              } // for(Z)
            } // for(Y)
          } // for(X)
        } // for(index_l3)
      } // for(index_l2)
    } // for(index_bt)
    #pragma omp flush(abort)
  } // for(index_l1)

  if (abort)
    return _FAILURE_;

  /* The analytical bispectra are now ready */
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt)
    if (pbi->bispectrum_type[index_bt] == analytical_bispectrum)
      pbi->bispectra_available[index_bt] = _TRUE_;
  

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
      struct spectra * psp,
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
                psp,
                pbi,
                pwb),
    pbi->error_message,
    pbi->error_message);
  

  printf_log_if (pbi->bispectra_verbose, 1, 
    " -> computing non-separable bispectra; r sampled %d times in [%g,%g]\n",
    pwb->r_size, pwb->r_min, pwb->r_max);
  
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

    /* Skip the bispectrum if it not of the non-separable type */
    if (pbi->bispectrum_type[index_bt] != non_separable_bispectrum)
      continue;

    for (int X = 0; X < pbi->bf_size; ++X) {

      pwb->X = X;

      // ==================================================================================================
      // =                              Determine the bispectrum to compute                               =
      // ==================================================================================================
    
      /* Which primordial shape is needed? */
      if ((pbi->has_galileon_model) && (index_bt==pbi->index_bt_galileon_gradient))
        pwb->shape_function = bispectra_galileon_gradient;

      else if ((pbi->has_galileon_model) && (index_bt==pbi->index_bt_galileon_time))
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
                    psp,
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
                      psp,
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
                        psp,
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
                        psp,
                        pbi,
                        pbi->bispectra[index_bt][X][Y][Z],
                        pwb),
            pbi->error_message,
            pbi->error_message);

        } // for(Z)
      } // for(Y)
    } // for(X)

    /* The bispectrum is ready */
    pbi->bispectra_available[index_bt] = _TRUE_;

  } // for(index_bt)
  
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
    struct spectra * psp,
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
    double pk = psp->pk_pt[index_k];
    
    pwb->k_window[index_k] = pow(k,2);

  }

  /* Inverse window function will have ptr->q_size elements */
  class_alloc (pwb->k_window_inverse, ptr->q_size*sizeof(double), pbi->error_message);
  for (int index_k=0; index_k < ptr->q_size; ++index_k) {

    double k = ptr->q[index_k];
    double pk = psp->pk[index_k];

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
  class_alloc (pwb->interpolated_integral, number_of_threads*sizeof(double*), pbi->error_message);
    
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
    class_alloc_parallel (pwb->interpolated_integral[thread], ptr->q_size*sizeof(double), pbi->error_message);
  
  } // parallel region
  
  if (abort) return _FAILURE_;

  return _SUCCESS_;
  
} // bispectra_non_separable_workspace_init







int bispectra_non_separable_integrate_over_k3 (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
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

      } // for(index_k1)
    } // for(index_r)
  } // for(index_l3)
    
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
    shared (ppt,pbs,ptr,ppm,psp,pbi,pwb,abort) private(thread)
  {
  
    #ifdef _OPENMP
    thread = omp_get_thread_num();
    #endif
  
  
    #pragma omp for schedule (dynamic)
    for (int index_k1 = 0; index_k1 < pwb->k_smooth_size; ++index_k1) {

      double k1 = pwb->k_smooth_grid[index_k1];
      double pk_1 = psp->pk_pt[index_k1];
      
      printf_log_if (pbi->bispectra_verbose, 2, 
        "     * computing the k3 integral for k1=%g, index_k1=%d\n",
        pwb->k_smooth_grid[index_k1], index_k1);
  
      /* We only need to consider those k2's that are equal to or larger than k1,
      as the shape function is assumed to be symmetric woth respect to k1<->k2<->k3 */      
      for (int index_k2 = 0; index_k2 <= index_k1; ++index_k2) {
    
        double k2 = pwb->k_smooth_grid[index_k2];
        double pk_2 = psp->pk_pt[index_k2];
    
        /* Set lower and upper limits for the integration over k3 */
        double k_tr_min = k_tr[0];
        double k_tr_max = k_tr[k_tr_size-1];
      
        /* Uncomment to restrict the integration to the k3 values that satisfy the
        triangular condition |k1-k2| <= k3 <= k1+k2. This will give faster but
        imprecise result. The reason is that the Dirac delta function that enforces
        the triangular condition, \delta(k1+k2-k3), is not explicit in the integral.
        In fact, we expanded \delta in Bessel functions using the Rayleigh expansion
        of a plane wave (see eq. 6.30 of my thesis http://arxiv.org/abs/1405.2280 or
        eqs. 9 and 10 of Fergusson and Shellard 2007). Analytically, this does not
        make any difference, but numerically it renders the integral unstable unless
        we take extra oscillations in the k3 direction. This means that to make the
        integral stable, we should extend the range of k3 into the forbidden region,
        similarly to what we do for the second-order bispectrum in bispectra2.c. */
        // k_tr_min = MAX (k_tr_min, fabs(k1-k2));
        // k_tr_max = MIN (k_tr_max, k1+k2);
    
        /* Discard invalid ranges, eg. when CLASS has a different k_max in ppt->k
        and ptr->q (very rare) */
        if (k_tr_max <= k_tr_min) {
          #pragma omp atomic
          pwb->count_memorised_for_integral_over_k3 += pbi->l_size * pwb->r_size;
          continue;
        }

        /* Determine the trapezoidal measure for the integration over k3 */
        int index_k3_lower;
        int index_k3_upper;

        class_call_parallel (trapezoidal_weights (
                               k_tr,
                               k_tr_size,
                               k_tr_min,
                               k_tr_max,
                               _FALSE_,
                               pwb->delta_k3[thread],
                               &index_k3_lower,
                               &index_k3_upper,
                               pbi->error_message),
          pbi->error_message,
          pbi->error_message);

        /* Size of the integration grid. If you chose to enforce the triangular
        condition in k_tr_min and k_tr_max, this might be zero, in which case
        this (k1,k2) does not contribute to the integral */
        int k3_size = index_k3_upper - index_k3_lower + 1;

        /* Compute the shape function for the current (k1,k2) slice */
        for (int index_k3=index_k3_lower; index_k3 <= index_k3_upper; ++index_k3)
          class_call_parallel (pwb->shape_function (
                        ppm,
                        psp,
                        pbi,
                        k1, k2, k_tr[index_k3],
                        pk_1, pk_2, psp->pk[index_k3],
                        &pwb->interpolated_integral[thread][index_k3]),
            pbi->error_message,
            pbi->error_message);
  
        /* We compute the integral over k3 for all possible l-values */
        for (int index_l3 = 0; index_l3 < pbi->l_size; ++index_l3) {

          /* Define the pointer to the first-order transfer functions as a function of k, for this value of l */
          int tt_size = ptr->tt_size[ppt->index_md_scalars];
          int l_size = ptr->l_size[ppt->index_md_scalars];

          double * transfer = &(ptr->transfer
            [ppt->index_md_scalars]
            [((ppt->index_ic_ad * tt_size + index_tt_k3) * l_size + index_l3) * k_tr_size]);

          /* Debug: print transfer function */
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
  
          } // for(index_r)          
        } // for(index_l3)
      } // for(index_k2)
    } // for(index_k1)
  } if (abort) return _FAILURE_; /* end of parallel region */    
  
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
      struct spectra * psp,
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
  
        } // for(index_k1)
      } // for(index_r)
    } // for(index_l2)
    
    printf_log_if (pbi->bispectra_verbose, 2, 
      "     * allocated ~ %.3g MB (%ld doubles) for the k2-integral array\n",
      pwb->count_allocated_for_integral_over_k2*sizeof(double)/1e6,
      pwb->count_allocated_for_integral_over_k2);
  
  } // if(pwb->Y==0)


  // ==============================================================================================================
  // =                                   Compute the INT_l3_l2(r, k1) integral                                    =
  // ==============================================================================================================
  
  /* Initialize counter for the number of integrals computed */
  pwb->count_memorised_for_integral_over_k2 = 0;
    
  abort = _FALSE_;
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,psp,pbi,pwb,abort) private (thread)
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
          Note that we pass the interpolated_integral array rather than accessing it from pwb, because
          it is thread dependent (while pwb isn't). */
          class_call_parallel (bispectra_non_separable_interpolate_over_k2(
                      ppr,
                      ppt,
                      pbs,
                      ptr,
                      ppm,
                      psp,
                      pbi,
                      index_r,
                      index_k1,
                      index_l3,
                      pwb->interpolated_integral[thread],
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

            /* Debug: print transfer function */
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
  
          } // for(index_l2)
        } // for(index_l3)
      } // for(index_k1)
    } // for(index_r)
  } if (abort) return _FAILURE_;  // parallel region
  
  
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
        } // for(index_k1)
        free (pwb->integral_over_k3[index_l2][index_r]);
      } // for(index_r)
      free (pwb->integral_over_k3[index_l2]);
    } // for(index_l2)    
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
      struct spectra * psp,
      struct bispectra * pbi,
      int index_r,
      int index_k1,
      int index_l3,
      double * interpolated_integral,
      struct bispectra_workspace_non_separable * pwb
      )
{

  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;

  /* So far, we always assumed that k1>=k2 because the shape function is symmetric wrt k1<->k2.
  Interpolating an array with such property is complicated, hence we build a temporary array
  where f(index_k1, index_k2) = f(index_k2, index_k1) when index_k1 < index_k2. */

  double * f;
  class_calloc (f, k_pt_size, sizeof(double), pbi->error_message);

  for (int index_k2=0; index_k2 < k_pt_size; ++index_k2) {

    f[index_k2] = (index_k1 > index_k2 ?
      pwb->integral_over_k3[index_l3][index_r][index_k1][index_k2]:
      pwb->integral_over_k3[index_l3][index_r][index_k2][index_k1]);

    /* Multiply by window function */
    f[index_k2] *= pwb->k_window[index_k2];

  }
  
  /* Interpolate f in all k_tr values */
  class_call (interpolate_array (
                k_pt,
                k_pt_size,
                f,
                ppr->transfers_k2_interpolation,
                _SPLINE_EST_DERIV_,
                k_tr,
                k_tr_size,
                interpolated_integral,
                NULL,
                pbi->error_message),
       pbi->error_message,
       pbi->error_message);

  /* Revert the effect of the window function */
  for (int index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  /* Debug: print the original array and the interpolation */
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
  
  free (f);

  return _SUCCESS_;

}






int bispectra_non_separable_integrate_over_k1 (
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct spectra * psp,
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
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,psp,pbi,pwb,abort) private (thread)
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
                        psp,
                        pbi,
                        index_r,
                        index_l3,
                        index_l2,
                        pwb->interpolated_integral[thread],
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
  
          } // for(index_l1)
        } // for(index_l2)
      } // for(index_l3)
    } // for(index_r)
  } if (abort) return _FAILURE_;  // parallel region
  
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
        } // for(index_r)
        free (pwb->integral_over_k2[index_l2][index_l1]);
      } // for(index_l1)
      free (pwb->integral_over_k2[index_l2]);
    } // for(index_l2)
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
      struct spectra * psp,
      struct bispectra * pbi,
      int index_r,
      int index_l3,
      int index_l2,
      double * interpolated_integral,
      struct bispectra_workspace_non_separable * pwb
      )
{
  
  /* Shortcuts */
  int k_pt_size = pwb->k_smooth_size;
  double * k_pt = pwb->k_smooth_grid;
  int k_tr_size = ptr->q_size;
  double * k_tr = ptr->q;
  
  /* Define the function to be interpolated, and multiply it by a window function */

  double * f;
  class_calloc (f, k_pt_size, sizeof(double), pbi->error_message);

  for (int index_k1=0; index_k1 < k_pt_size; ++index_k1) {
    
    f[index_k1] = pwb->integral_over_k2[index_l3][index_l2][index_r][index_k1];
    f[index_k1] *= pwb->k_window[index_k1];
  }
  
  /* Interpolate f in all k_tr values */
  class_call (interpolate_array (
                k_pt,
                k_pt_size,
                f,
                ppr->transfers_k1_interpolation,
                _SPLINE_EST_DERIV_,
                k_tr,
                k_tr_size,
                interpolated_integral,
                NULL,
                pbi->error_message),
       pbi->error_message,
       pbi->error_message);

  /* Revert the effect of the window function */
  for (int index_k_tr = 0; index_k_tr < k_tr_size; ++index_k_tr)
    interpolated_integral[index_k_tr] *= pwb->k_window_inverse[index_k_tr];

  /* Debug: print the original array and the interpolation */
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

  free (f);

  return _SUCCESS_;
  
}









int bispectra_non_separable_integrate_over_r (
    struct precision * ppr,
    struct background * pba,
    struct perturbs * ppt,
    struct bessels * pbs,
    struct transfers * ptr,
    struct primordial * ppm,
    struct spectra * psp,
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
  #pragma omp parallel shared (ppt,pbs,ptr,ppm,psp,pbi,pwb,abort)
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

            /* Debug: output intermediate results on stderr for a custom (l1,l2,l3) configuration */
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
   
          } // for(index_r)

          /* Fill the bispectrum array with the result for this set of (l1,l2,l3) */
          bispectrum[index_l1_l2_l3] = integral;

          /* Account for the overall (2/pi)^3 factor coming from the bispectrum formula. This factor is seen
          (see, for instance, eq. 17 of Fergusson & Shellard 2007). */
          bispectrum[index_l1_l2_l3] *= pow(2./_PI_,3);

          /* Update the counter */
          #pragma omp atomic
          pbi->count_memorised_bispectra++;

          /* Debug: output the integral as a function of r on stderr for a custom (l1,l2,l3) */
          // if ( (l1==l2) && (l2==l3) ) {
          //   fprintf(stderr, "%12d %17.7g\n", l1, pwb->integral_over_r[index_l1][index_l2][index_l3-index_l3_min]);
          // }

        } // for(index_l3)
      } // for(index_l2)
      
      #pragma omp flush(abort)
      
    } // for(index_l1)
  } if (abort) return _FAILURE_;  // parallel region
  
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
    free(pwb->interpolated_integral[thread]);
    
  }  if (abort) return _FAILURE_;
  
  free(pwb->k3_grid);
  free(pwb->delta_k3);
  free(pwb->interpolated_integral);
  
  return _SUCCESS_;
    
}





/**
 * Compute the lensing correction to the CMB bispectra.
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
    R += l*(l+1.0)*(2*l+1.0)/(4*_PI_) * pbi->cls[psp->index_ct_pp][l-2];
  
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

        /* Compute 3j factor needed to convert the bispectrum to the reduced bispectrum */
        double threej[2*pbi->l_max+1];
        int l3_min_3j, l3_max_3j;

        double min_D, max_D;
        class_call_parallel (drc3jj (
                               l1, l2, 0, 0,
                               &min_D, &max_D,
                               threej,
                               (2*pbi->l_max+1),
                               pbi->error_message),
          pbi->error_message,
          pbi->error_message);
        l3_min_3j = (int)(min_D + _EPS_);
        l3_max_3j = (int)(max_D + _EPS_);

        for (int X2=0; X2 < pbi->bf_size; ++X2) {

          int F_X2 = pbi->field_spin[X2];

          for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

            int l3 = pbi->l[index_l3];
            long int index_l1_l2_l3 = pbi->index_l1_l2_l3[index_l1][index_l1-index_l2][index_l3_max-index_l3];

            /* Debug: compute lensing for a reduced set of l */
            // int l_long = 10;
            // if ( ! (l1==l_long || l2==l_long || l3==l_long) )
            // // if ( ! (((l1==l_long) && (l2==l3)) || ((l3==l_long) && (l1==l2))) )
            //   continue;

            for (int X3=0; X3 < pbi->bf_size; ++X3) {

              int F_X3 = pbi->field_spin[X3];


              // ====================================================================================
              // =                                   Overall R factor                               =
              // ====================================================================================

              double R_factor = - 0.5 * R * (
                (l1+F_X1)*(l1-F_X1+1.0) + (l1-F_X1)*(l1+F_X1+1.0) +
                (l2+F_X2)*(l2-F_X2+1.0) + (l2-F_X2)*(l2+F_X2+1.0) +
                (l3+F_X3)*(l3-F_X3+1.0) + (l3-F_X3)*(l3+F_X3+1.0)
              );
                            

              // ====================================================================================
              // =                                 Convolution factor                               =
              // ====================================================================================

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
                
              double convolution_factor = convolution_factor_123 +
                                          convolution_factor_231 +
                                          convolution_factor_312;
                
                
              // -------------------------------------------------------------------------------
              // -                        Convert to reduced bispectrum                        -
              // -------------------------------------------------------------------------------
                
              /* So far we have computed the convolution correction for the angle-averaged
              bispectrum, while SONG expects the one for the reduced bispectrum */

              double conversion_factor = sqrt((2*l1+1.0)*(2*l2+1.0)*(2*l3+1.0)/(4*_PI_)) * threej[l3-l3_min_3j];

              if ((l1+l2+l3)%2==0)
                convolution_factor /= conversion_factor;
              else
                convolution_factor = 0;



              // ====================================================================================
              // =                             Final lensing correction                             =
              // ====================================================================================

              double b = pbi->bispectra[index_bt][X1][X2][X3][index_l1_l2_l3];

              pbi->lensing_correction[index_bt][X1][X2][X3][index_l1_l2_l3] =
                b * R_factor/4 +
                convolution_factor;

            }
          }
        }
      }
    }
  }
  
  if (abort)
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
 * Compute the convolution part of lensing for the input (l1,l2,l3) configuration
 * of the bispectrum.
 *
 * The convolution is a sum over three dummy multipoles (l,p,q) symmetric with respect
 * to p<->q. Each addend involves two cosmological functions, the lensing potential in l
 * and the bispectrum in (l1,p,q), and a bunch of geometrical factors involving 3j
 * and 6j symbols.
 * 
 * Our strategy is to first perform the sum over p and q taking all points, and then
 * integrating over l using only the nodes in pbi->l, with the trapezoidal rule. Ideally
 * we would use the trapezoidal rule also for the p and q directions, thus cutting the
 * computation time drastically. However, the presence of the 3j and 6j symbols makes
 * the integrand function discontinuous and thus this approach would be very imprecise
 * (see bispectra_lensing_convolution_nodes()).
 * 
 * We obtain the bispectrum in (l1,p,q) using bilinear interpolation.
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
  
  // TODO:
  // *) Run lensing on CMB-lensing and see if it matches lensed CMB-lensing
  // *) Look again at interpolation of 3j symbol
  // *) Streamline fisher_compute_matrix() (documentation as well)
  // *) Get rid of fisher_compute_matrix_nodes() 
  // *) Sort out interaction of output functions (bispectra and perturbs2)
  // *) with store_sources and load_sources
  // *) Sort out interaction of bispectra mesh interpolation and lensed
  //    bispectra (ie right now we need to reinterpolate pbi->bispectra
  //    after we lens). Possible solution: generate separate meshes for
  //    lensed and unlensed bispectra.
  // *) Lensing variance is broken right now, because we inverted l1 and l3
  //    in fisher_compute_matrix(). Maybe we need to invert it in the 
  //    fisher_lensing_variance() function?

  int l1 = pbi->l[index_l1], F_X1 = pbi->field_spin[X1];
  int l2 = pbi->l[index_l2], F_X2 = pbi->field_spin[X2];
  int l3 = pbi->l[index_l3], F_X3 = pbi->field_spin[X3];

  class_test (!is_triangular_int(l1,l2,l3),
    pbi->error_message,
    "l1=%d, l2=%d and l3=%d do not form a triangle", l1, l2, l3);

  /* Since the bispectrum does not depend on l, we store it in a (p,q) array
  and reuse it when l changes */
  double (*bispectrum_pq)[pbi->full_l_size] = calloc (pbi->full_l_size*pbi->full_l_size, sizeof(double));
  short (*is_filled_pq)[pbi->full_l_size] = calloc (pbi->full_l_size*pbi->full_l_size, sizeof(short));

  /* Trapezoidal weights for the sum over the l1 direction */
  double * delta_l;
  class_alloc (delta_l, pbi->l_size*sizeof(double), pbi->error_message);

  class_call (trapezoidal_weights_int (
                pbi->l,
                pbi->l_size,
                pbi->l[0],
                pbi->l[pbi->l_size-1],
                delta_l,
                NULL,
                NULL,
                pbi->error_message),
    pbi->error_message,
    pbi->error_message);

  /* Initialise output */
  *result = 0;

  /* Debug: store the integrand function */
  // double (*integrand_lp)[pbi->l_size] = calloc (pbi->l_size*pbi->full_l_size, sizeof (double));
  // double * integrand_l = calloc (pbi->l_size, sizeof (double));
  // int l_long = 10;
  // int l_short = 1500;
  // short print_output = 0;
  // if ( ((l1==l_long) && (l2==l_short) && (l3==l_short)) ||
  //      ((l2==l_long) && (l3==l_short) && (l1==l_short)) ||
  //      ((l3==l_long) && (l1==l_short) && (l2==l_short)))
  //   print_output = 1;


  // ====================================================================================
  // =                                     Loop on l                                    =
  // ====================================================================================

  for (int index_l=0; index_l < pbi->l_size; ++index_l) {

    int l = pbi->l[index_l];

    /* Uncomment to set an upper limit for l */
    // if (l > 300)
    //   continue;

    /* Lensing potential in l */
    double C_PP = pbi->cls[psp->index_ct_pp][l-2];
    
    /* Range of p dictated by the triangular inequality on l,l2,p */
    int p_min = MAX (abs(l-l2), 2);
    int p_max = MIN (l+l2, pbi->l_max);

    /* Compute 3j symbol for all values of p */
    double threej_p[2*pbi->l_max+1];
    int p_min_3j, p_max_3j;

    double min_D, max_D;
    class_call (drc3jj (
                  l, l2, 0, -F_X2,
                  &min_D, &max_D,
                  threej_p,
                  (2*pbi->l_max+1),
                  pbi->error_message),
      pbi->error_message,
      pbi->error_message);
    p_min_3j = (int)(min_D + _EPS_);
    p_max_3j = (int)(max_D + _EPS_);

    for (int p=p_min_3j; p <= p_max_3j; ++p)
      threej_p[p-p_min_3j] *= sqrt((2*l2+1.0)*(2*l+1.0)*(2*p+1.0)/(4*_PI_)) * 0.5 * (-l2*(l2+1.0)+l*(l+1.0)+p*(p+1.0));

    class_test ((p_min != MAX (p_min_3j, 2)) || (p_max != MIN (p_max_3j, pbi->l_max)),
      pbi->error_message,
      "wrong 3j limits");


    // ====================================================================================
    // =                                     Loop on p                                    =
    // ====================================================================================

    for (int p=p_min; p <= p_max; ++p) {

      class_test (!is_triangular_int (l2,p,l),
         pbi->error_message,
         "p not triangular");

      /* Range of q dictated by the triangular inequality on (l,q,l3) and (l1,q,p) */
      int q_min = MAX (MAX (abs(l3-l), abs(l1-p)), 2);
      int q_max = MIN (MIN (l3+l, l1+p), pbi->l_max);

      /* Compute 3j symbol for all values of q */
      double threej_q[2*pbi->l_max+1];
      int q_min_3j, q_max_3j;

      double min_D, max_D;
      class_call (drc3jj (
                    l, l3, 0, -F_X3,
                    &min_D, &max_D,
                    threej_q,
                    (2*pbi->l_max+1),
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      q_min_3j = (int)(min_D + _EPS_);
      q_max_3j = (int)(max_D + _EPS_);

      for (int q=q_min_3j; q <= q_max_3j; ++q)
        threej_q[q-q_min_3j] *= sqrt((2*l3+1.0)*(2*l+1.0)*(2*q+1.0)/(4*_PI_)) * 0.5 * (-l3*(l3+1.0)+l*(l+1.0)+q*(q+1.0));

      class_test ((q_min_3j != abs(l3-l) || (q_max_3j != (l3+l))),
        pbi->error_message,
        "wrong 3j limits");

      /* Compute the 6j symbol for all values of q */
      double sixj_q[2*pbi->l_max+1];
      int q_min_6j, q_max_6j;

      class_call (drc6j (
                    /*q,*/ p, l1, l2, l3, l,
                    &min_D, &max_D,
                    sixj_q,
                    (2*pbi->l_max+1),
                    pbi->error_message
                    ),
        pbi->error_message,
        pbi->error_message);
      q_min_6j = (int)(min_D + _EPS_);
      q_max_6j = (int)(max_D + _EPS_);

      class_test ((q_min != MAX (q_min_6j, 2)) || (q_max != (MIN (q_max_6j, pbi->l_max))),
        pbi->error_message,
        "wrong 3j limits");

      /* Additional 3J factor needed to convert SONG reduced bispectrum to the
      angle averaged one */
      double threej_000[2*pbi->l_max+1];
      int q_min_3j_000, q_max_3j_000;

      /* Before computing the 3j, check whether we have already computed it
      in all the required q values for this p. If so, we shally recycle it. */
      int is_filled = 1;
      for (int q=q_min; q <= q_max; ++q)
        if (!(is_filled = is_filled_pq[p-2][q-2]))
          break;

      if (!is_filled) {
        
        class_call (drc3jj (
                      l1, p, 0, 0,
                      &min_D, &max_D,
                      threej_000,
                      (2*pbi->l_max+1),
                      pbi->error_message),
          pbi->error_message,
          pbi->error_message);
        q_min_3j_000 = (int)(min_D + _EPS_);
        q_max_3j_000 = (int)(max_D + _EPS_);

        for (int q=q_min_3j_000; q <= q_max_3j_000; ++q)
          threej_000[q-q_min_3j_000] *= sqrt((2*l1+1.0)*(2*p+1.0)*(2*q+1.0)/(4*_PI_));

      }
      

      // ====================================================================================
      // =                                     Loop on q                                    =
      // ====================================================================================

      for (int q=q_min; q <= q_max; ++q) {

        class_test (!is_triangular_int (l3,q,l) || !is_triangular_int (l1,p,q),
           pbi->error_message,
           "q not triangular");

        int X2_ = X2;
        int X3_ = X3;


        // --------------------------------------------------------------------------------
        // -                            Extract the bispectrum                            -
        // --------------------------------------------------------------------------------

        if (!is_filled_pq[p-2][q-2]) {

          double reduced_bispectrum;

          class_call (bispectra_at_l2l3 (
                        ppr,
                        ptr,
                        psp,
                        ple,
                        pbi,
                        index_bt,
                        index_l1,
                        l1, p, q,
                        X1, X2_, X3_,
                        &reduced_bispectrum,
                        NULL),
            pbi->error_message,
            pbi->error_message);

          /* Compute the actual bispectrum from the reduced one */
          bispectrum_pq[p-2][q-2] = threej_000[q-q_min_3j_000] * reduced_bispectrum;

          is_filled_pq[p-2][q-2] = 1;

        }

        double integrand = C_PP *
                           ALTERNATING_SIGN (l1+l2+p) * threej_p[p-p_min_3j] * threej_q[q-q_min_3j] * sixj_q[q-q_min_6j] *
                           bispectrum_pq[p-2][q-2] *
                           delta_l[index_l];

        class_test (isnan (integrand),
          pbi->error_message,
          "integrand nan for l1=%d,l2=%d,l3=%d,l=%d,p=%d,q=%d; b=%g, 3j(p)=%g, 3j(q)=%g, 6j(q)=%g",
          l1,l2,l3,l,p,q, bispectrum_pq[p-2][q-2], threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);

        *result += integrand;
        
        /* Debug: save the integrand function */
        // integrand_lp[index_l][p-2] += integrand;

        /* Debug: print the integrand as a function of q */
        // if (print_output && (l==100) && (p==100)) {
        //   fprintf (stderr, "%6d %17g %17g %17g %17g %17g\n", q, integrand, bispectrum_pq[p-2][q-2],
        //     threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);
        //   if (q == q_max) fprintf (stderr, "\n\n");
        // }

      } // for q

      /* Debug: print the integrand as a function of p */
      // if (print_output && (l==28)) {
      //   fprintf (stderr, "%6d %17g %17g\n", p, integrand_lp[index_l][p-2], threej_p[p-p_min_3j]);
      //   if (p == p_max) fprintf (stderr, "\n\n");
      // }

      /* Debug: save the integrand function */
      // integrand_l[index_l] += integrand_lp[index_l][p-2];

    } // for p

    /* Debug: print the integrand as a function of l */
    // if (print_output) {
    //   fprintf (stderr, "%6d %17g\n", l, integrand_l[index_l]);
    //   if (l == pbi->l_max) fprintf (stderr, "\n\n");
    // }

  } // for index_l

  free (bispectrum_pq);
  free (is_filled_pq);
  free (delta_l);

  /* Debug: free the debug arrays */
  // free (integrand_lp);
  // free (integrand_l);
  
  return _SUCCESS_;

}



/**
 * Same as bispectra_lensing_convolution(), but using the trapezoidal rule
 * for the p and q directions as well.
 * 
 * This function implicitly assumes that the integrand function is smooth in
 * p and q. This is not the case because the integrand includes 3j and 6j
 * symbols with alternating signs. Therefore, the only purpose of this
 * function is to give a quick order of magnitude estimate of the effect of
 * lensing.
 */

int bispectra_lensing_convolution_nodes (
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

  /* Since the bispectrum does not depend on l, we store it in a (p,q) array
  and reuse it when l changes */
  double (*bispectrum_pq)[pbi->l_size] = calloc (pbi->l_size*pbi->l_size, sizeof(double));
  short (*is_filled_pq)[pbi->l_size] = calloc (pbi->l_size*pbi->l_size, sizeof(short));

  /* Weights for the trapezoidal sum along the l direction */
  double delta_l[pbi->l_size];
  delta_l[0] = (pbi->l[1] - pbi->l[0] + 1)/2.0;
  for (int index_l=1; index_l < (pbi->l_size-1); ++index_l)
    delta_l[index_l] = (pbi->l[index_l+1] - pbi->l[index_l-1])/2.0;
  delta_l[pbi->l_size-1] = (pbi->l[pbi->l_size-1] - pbi->l[pbi->l_size-2] + 1)/2.0;

  /* Initialise output */
  *result = 0;

  /* Debug: store the integrand function */
  // double (*integrand_lp)[pbi->l_size] = calloc (pbi->l_size*pbi->l_size, sizeof (double));
  // double * integrand_l = calloc (pbi->l_size, sizeof (double));
  // int l_long = 10;
  // int l_short = 1500;
  // short print_output = 0;
  // if ( ((l1==l_long) && (l2==l_short) && (l3==l_short)) ||
  //      ((l2==l_long) && (l3==l_short) && (l1==l_short)) ||
  //      ((l3==l_long) && (l1==l_short) && (l2==l_short)))
  //   print_output = 1;


  // ======================================================================================
  // =                                      Loop on l                                     =
  // ======================================================================================

  for (int index_l=0; index_l < pbi->l_size; ++index_l) {

    int l = pbi->l[index_l];

    /* Uncomment to set an upper limit for l */
    // if (l > 300)
    //   continue;

    /* Lensing potential in l */
    double C_PP = pbi->cls[psp->index_ct_pp][l-2];

    /* Range of p dictated by the triangular inequality on l,l2,p */
    int p_min = MAX (abs(l-l2), 2);
    int p_max = MIN (l+l2, pbi->l_max);

    /* Enforce triagular inequality for p */
    int index_p_min = 0;
    while (pbi->l[index_p_min] < p_min)
      ++index_p_min;

    int index_p_max = pbi->l_size-1;
    while (pbi->l[index_p_max] > p_max)
      --index_p_max;

    int p_size = index_p_max - index_p_min + 1;


    // ---------------------------------------------------------------------------------
    // -                                   3J & weights                                -
    // ---------------------------------------------------------------------------------

    double delta_p[pbi->l_size]; /* will contain the trapezoidal interpolation weights */
    double threej_p[2*pbi->l_max+1]; /* will contain the 3j symbol for all p */
    int p_min_3j, p_max_3j; /* limits in p of the 3j symbol */
    double min_D, max_D; /* temporary variables for the 3j computation */

    if (p_size > 0) {

      /* If only one node falls within p_min and p_max, use closest-neighbour interpolation,
      that is, assume the bispectrum value of that node for the whole p domain */

      if (p_size == 1) {
        delta_p[0] = p_max - p_min + 1;
      }

      /* Otherwise, compute trapezoidal weights */

      if (p_size > 1) {

        /* Build integration grid making sure to include p_min and p_max. This is
        equivalent to assuming closest-neighbour interpolation for those points
        outside the sampled range but inside the triangular region */

        double p[p_size];

        p[0] = p_min;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          p[index_p] = pbi->l[index_p_min+index_p];

        p[p_size-1] = p_max;

        /* Build trapezoidal integration measure */

        delta_p[0] = (p[1] - p[0] + 1)/2.0;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          delta_p[index_p] = (p[index_p+1] - p[index_p-1])/2.0;

        delta_p[p_size-1] = (p[p_size-1] - p[p_size-2] + 1)/2.0;

      }

      /* Compute 3j symbol for all values of p */

      class_call (drc3jj (
                    /*p,*/ l, l2, 0, -F_X2,
                    &min_D, &max_D,
                    threej_p,
                    (2*pbi->l_max+1),
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      p_min_3j = (int)(min_D + _EPS_);
      p_max_3j = (int)(max_D + _EPS_);

      for (int p=p_min_3j; p <= p_max_3j; ++p)
        threej_p[p-p_min_3j] *= sqrt((2*l2+1.0)*(2*l+1.0)*(2*p+1.0)/(4*_PI_)) * 0.5 * (-l2*(l2+1.0)+l*(l+1.0)+p*(p+1.0));

      class_test ((p_min != MAX (p_min_3j, 2)) || (p_max != MIN (p_max_3j, pbi->l_max)),
        pbi->error_message,
        "wrong 3j limits");

    }


    // ======================================================================================
    // =                                      Loop on p                                     =
    // ======================================================================================

    for (int index_p=index_p_min; index_p <= index_p_max; ++index_p) {

      int p = pbi->l[index_p];

      class_test (!is_triangular_int (l2,p,l),
         pbi->error_message,
         "p not triangular");

      /* Range of q dictated by the triangular inequality on (l,q,l3) and (l1,q,p) */
      int q_min = MAX (MAX (abs(l3-l), abs(l1-p)), 2);
      int q_max = MIN (MIN (l3+l, l1+p), pbi->l_max);

      /* Enforce triagular inequality for q */
      int index_q_min = MAX (pbi->index_l_triangular_min[index_l1][index_p], pbi->index_l_triangular_min[index_l3][index_l]);
      int index_q_max = MIN (pbi->index_l_triangular_max[index_l1][index_p], pbi->index_l_triangular_max[index_l3][index_l]);
      int q_size = index_q_max - index_q_min + 1;


      // ---------------------------------------------------------------------------------
      // -                                 3J, 6J & weights                              -
      // ---------------------------------------------------------------------------------

      double delta_q[pbi->l_size]; /* will contain the trapezoidal interpolation weights */
      double threej_q[2*pbi->l_max+1]; /* will contain the 3j symbol for all q */
      int q_min_3j, q_max_3j; /* limits in q of the 3j symbol */
      double sixj_q[2*pbi->l_max+1]; /* will contain the 6j symbol for all q */
      int q_min_6j, q_max_6j; /* limits in q of the 6j symbol */
      double threej_000[2*pbi->l_max+1]; /* will contain 3j symbol to convert the reduced bispectrum */
      int q_min_3j_000, q_max_3j_000; /* limits in q for the conversion 3j */

      if (q_size > 0) {

        /* If only one node falls within q_min and q_max, use closest-neighbour interpolation,
        that is, assume the bispectrum value of that node for the whole q domain */

        if (q_size == 1) {
          delta_q[0] = q_max - q_min + 1;
        }

        /* Otherwise, compute trapezoidal weights */

        if (q_size > 1) {

          /* Build integration grid making sure to include q_min and q_max. This is
          equivalent to assuming closest-neighbour interpolation for those points
          outside the sampled range but inside the triangular region */

          double q[q_size];

          q[0] = q_min;

          for (int index_q=1; index_q < q_size-1; ++index_q)
            q[index_q] = pbi->l[index_q_min+index_q];

          q[q_size-1] = q_max;

          /* Build trapezoidal integration measure */

          delta_q[0] = (q[1] - q[0] + 1)/2.0;

          for (int index_q=1; index_q < q_size-1; ++index_q)
            delta_q[index_q] = (q[index_q+1] - q[index_q-1])/2.0;

          delta_q[q_size-1] = (q[q_size-1] - q[q_size-2] + 1)/2.0;

          /* Debug: print the grid in q */
          // printf ("~~~ l1=%d, l2=%d, l3=%d, l=%d, q=%d, q_size=%d, q_min=%d, q_max=%d, index_q_min=%d, index_q_max=%d\n",
          //   l1, l2, l3, l, q, q_size, q_min, q_max, index_q_min, index_q_max);
          // for (int index_q=0; index_q < q_size; ++index_q)
          //   printf ("q[%d] = %g, delta_q = %g\n", index_q, q[index_q], delta_q[index_q]);

        }

        /* Compute 3j symbol for all values of q */

        class_call (drc3jj (
                      /*q,*/ l, l3, 0, -F_X3,
                      &min_D, &max_D,
                      threej_q,
                      (2*pbi->l_max+1),
                      pbi->error_message),
          pbi->error_message,
          pbi->error_message);
        q_min_3j = (int)(min_D + _EPS_);
        q_max_3j = (int)(max_D + _EPS_);

        for (int q=q_min_3j; q <= q_max_3j; ++q)
          threej_q[q-q_min_3j] *= sqrt((2*l3+1.0)*(2*l+1.0)*(2*q+1.0)/(4*_PI_)) * 0.5 * (-l3*(l3+1.0)+l*(l+1.0)+q*(q+1.0));

        class_test ((q_min_3j != abs(l3-l) || (q_max_3j != (l3+l))),
          pbi->error_message,
          "wrong 3j limits");

        /* Compute the 6j symbol for all values of q */

        class_call (drc6j (
                      /*q,*/ p, l1, l2, l3, l,
                      &min_D, &max_D,
                      sixj_q,
                      (2*pbi->l_max+1),
                      pbi->error_message
                      ),
          pbi->error_message,
          pbi->error_message);
        q_min_6j = (int)(min_D + _EPS_);
        q_max_6j = (int)(max_D + _EPS_);

        class_test ((q_min != MAX (q_min_6j, 2)) || (q_max != (MIN (q_max_6j, pbi->l_max))),
          pbi->error_message,
          "wrong 3j limits");

        /* Before computing the conversion 3j, check whether we have already computed it
        in all the required q values for this p. If so, we shally recycle it. */
        int is_filled = 1;
        for (int index_q=index_q_min; index_q <= index_q_max; ++index_q)
          if (!(is_filled = is_filled_pq[index_p][index_q]))
            break;

        if (!is_filled) {

          class_call (drc3jj (
                        l1, p, 0, 0,
                        &min_D, &max_D,
                        threej_000,
                        (2*pbi->l_max+1),
                        pbi->error_message),
            pbi->error_message,
            pbi->error_message);
          q_min_3j_000 = (int)(min_D + _EPS_);
          q_max_3j_000 = (int)(max_D + _EPS_);

          for (int q=q_min_3j_000; q <= q_max_3j_000; ++q)
            threej_000[q-q_min_3j_000] *= sqrt((2*l1+1.0)*(2*p+1.0)*(2*q+1.0)/(4*_PI_));
        }

      } // if(q_size>0)


      // ======================================================================================
      // =                                      Loop on q                                     =
      // ======================================================================================

      for (int index_q=index_q_min; index_q <= index_q_max; ++index_q) {

        int q = pbi->l[index_q];

        class_test (!is_triangular_int (l3,q,l) || !is_triangular_int (l1,p,q),
           pbi->error_message,
           "q not triangular");

        int X2_ = X2;
        int X3_ = X3;

        if (!is_filled_pq[index_p][index_q]) {

          double reduced_bispectrum;

          class_call (bispectra_at_node (
                        pbi,
                        index_bt,
                        index_l1, index_p, index_q,
                        X1, X2_, X3_,
                        &reduced_bispectrum,
                        NULL),
            pbi->error_message,
            pbi->error_message);

          /* Compute the actual bispectrum from the reduced one */
          bispectrum_pq[index_p][index_q] = threej_000[q-q_min_3j_000] * reduced_bispectrum;

          is_filled_pq[index_p][index_q] = 1;

        }

        double integrand = C_PP *
                           ALTERNATING_SIGN (l1+l2+p) * threej_p[p-p_min_3j] * threej_q[q-q_min_3j] * sixj_q[q-q_min_6j] *
                           delta_l[index_l] * delta_p[index_p-index_p_min] * delta_q[index_q-index_q_min] *
                           bispectrum_pq[index_p][index_q];

        class_test (isnan (integrand),
          pbi->error_message,
          "integrand nan for l1=%d,l2=%d,l3=%d,l=%d,p=%d,q=%d; b=%g, 3j(p)=%g, 3j(q)=%g, 6j(q)=%g",
          l1,l2,l3,l,p,q, bispectrum_pq[index_p][index_q], threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);

        *result += integrand;

        /* Debug: save the integrand function */
        // integrand_lp[index_l][index_p] += integrand;

        /* Debug: print the integrand as a function of q */
        // if (print_output && (l==100) && (p==100)) {
        //   fprintf (stderr, "%6d %17g %17g %17g %17g %17g\n", q, integrand, bispectrum_pq[index_p][index_q],
        //     threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);
        //   if (q == q_max) fprintf (stderr, "\n\n");
        // }

      } // for index_q

      /* Debug: print the integrand as a function of p */
      // if (print_output && (l==28)) {
      //   fprintf (stderr, "%6d %17g %17g\n", p, integrand_lp[index_l][index_p], threej_p[p-p_min_3j]);
      //   if (p == p_max) fprintf (stderr, "\n\n");
      // }

      /* Debug: save the integrand function */
      // integrand_l[index_l] += integrand_lp[index_l][index_p];

    } // for index_p

    /* Debug: print the integrand as a function of l */
    // if (print_output) {
    //   fprintf (stderr, "%6d %17g\n", l, integrand_l[index_l]);
    //   if (l == pbi->l_max) fprintf (stderr, "\n\n");
    // }

  } // for index_l

  free (bispectrum_pq);
  free (is_filled_pq);

  /* Debug: free the debug arrays */
  // free (integrand_lp);
  // free (integrand_l);

  return _SUCCESS_;

}




/**
 * Same as bispectra_lensing_convolution(), but using the trapezoidal rule
 * for the p direction as well and using 1D linear interpolation for the
 * bispectrum.
 *
 * This function implicitly assumes that the integrand function is smooth in p.
 * This is not the case because the integrand includes 3j and 6j symbols with
 * alternating signs. Therefore, the only purpose of this function is to give
 * a quick order of magnitude estimate of the effect of lensing.
 */

int bispectra_lensing_convolution_linear (
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

  /* Since the bispectrum does not depend on l, we store it in a (p,q) array
  and reuse it when l changes */
  double (*bispectrum_pq)[pbi->full_l_size] = calloc (pbi->l_size*pbi->full_l_size, sizeof(double));
  short (*is_filled_pq)[pbi->full_l_size] = calloc (pbi->l_size*pbi->full_l_size, sizeof(short));

  /* Weights for the trapezoidal sum along the l direction */
  double delta_l[pbi->l_size];
  delta_l[0] = (pbi->l[1] - pbi->l[0] + 1)/2.0;
  for (int index_l=1; index_l < (pbi->l_size-1); ++index_l)
    delta_l[index_l] = (pbi->l[index_l+1] - pbi->l[index_l-1])/2.0;
  delta_l[pbi->l_size-1] = (pbi->l[pbi->l_size-1] - pbi->l[pbi->l_size-2] + 1)/2.0;

  /* Initialise output */
  *result = 0;

  /* Debug: store the integrand function */
  // double (*integrand_lp)[pbi->l_size] = calloc (pbi->l_size*pbi->l_size, sizeof (double));
  // double * integrand_l = calloc (pbi->l_size, sizeof (double));
  // int l_long = 10;
  // int l_short = 1500;
  // short print_output = 0;
  // if ( ((l1==l_long) && (l2==l_short) && (l3==l_short)) ||
  //      ((l2==l_long) && (l3==l_short) && (l1==l_short)) ||
  //      ((l3==l_long) && (l1==l_short) && (l2==l_short)))
  //   print_output = 1;


  // ======================================================================================
  // =                                      Loop on l                                     =
  // ======================================================================================

  for (int index_l=0; index_l < pbi->l_size; ++index_l) {

    int l = pbi->l[index_l];

    /* Uncomment to set an upper limit for l */
    // if (l > 300)
    //   continue;

    /* Lensing potential in l */
    double C_PP = pbi->cls[psp->index_ct_pp][l-2];

    /* Range of p dictated by the triangular inequality on l,l2,p */
    int p_min = MAX (abs(l-l2), 2);
    int p_max = MIN (l+l2, pbi->l_max);

    /* Enforce triagular inequality for p */
    int index_p_min = 0;
    while (pbi->l[index_p_min] < p_min)
      ++index_p_min;

    int index_p_max = pbi->l_size-1;
    while (pbi->l[index_p_max] > p_max)
      --index_p_max;

    int p_size = index_p_max - index_p_min + 1;


    // ---------------------------------------------------------------------------------
    // -                                   3J & weights                                -
    // ---------------------------------------------------------------------------------

    double delta_p[pbi->l_size]; /* will contain the trapezoidal interpolation weights */
    double threej_p[2*pbi->l_max+1]; /* will contain the 3j symbol for all p */
    int p_min_3j, p_max_3j; /* limits in p of the 3j symbol */
    double min_D, max_D; /* temporary variables for the 3j computation */

    if (p_size > 0) {

      /* If only one node falls within p_min and p_max, use closest-neighbour interpolation,
      that is, assume the bispectrum value of that node for the whole p domain */

      if (p_size == 1) {
        delta_p[0] = p_max - p_min + 1;
      }

      /* Otherwise, compute trapezoidal weights */

      if (p_size > 1) {

        /* Build integration grid making sure to include p_min and p_max. This is
        equivalent to assuming closest-neighbour interpolation for those points
        outside the sampled range but inside the triangular region */

        double p[p_size];

        p[0] = p_min;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          p[index_p] = pbi->l[index_p_min+index_p];

        p[p_size-1] = p_max;

        /* Build trapezoidal integration measure */

        delta_p[0] = (p[1] - p[0] + 1)/2.0;

        for (int index_p=1; index_p < p_size-1; ++index_p)
          delta_p[index_p] = (p[index_p+1] - p[index_p-1])/2.0;

        delta_p[p_size-1] = (p[p_size-1] - p[p_size-2] + 1)/2.0;

      }

      /* Compute 3j symbol for all values of p */

      class_call (drc3jj (
                    /*p,*/ l, l2, 0, -F_X2,
                    &min_D, &max_D,
                    threej_p,
                    (2*pbi->l_max+1),
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      p_min_3j = (int)(min_D + _EPS_);
      p_max_3j = (int)(max_D + _EPS_);

      for (int p=p_min_3j; p <= p_max_3j; ++p)
        threej_p[p-p_min_3j] *= sqrt((2*l2+1.0)*(2*l+1.0)*(2*p+1.0)/(4*_PI_)) * 0.5 * (-l2*(l2+1.0)+l*(l+1.0)+p*(p+1.0));

      class_test ((p_min != MAX (p_min_3j, 2)) || (p_max != MIN (p_max_3j, pbi->l_max)),
        pbi->error_message,
        "wrong 3j limits");

    }


    // ======================================================================================
    // =                                      Loop on p                                     =
    // ======================================================================================

    for (int index_p=index_p_min; index_p <= index_p_max; ++index_p) {

      int p = pbi->l[index_p];

      class_test (!is_triangular_int (l2,p,l),
         pbi->error_message,
         "p not triangular");

      /* Range of q dictated by the triangular inequality on (l,q,l3) and (l1,q,p) */
      int q_min = MAX (MAX (abs(l3-l), abs(l1-p)), 2);
      int q_max = MIN (MIN (l3+l, l1+p), pbi->l_max);

      /* Enforce triagular inequality for q */
      int index_q_min = MAX (pbi->index_l_triangular_min[index_l1][index_p], pbi->index_l_triangular_min[index_l3][index_l]);
      int index_q_max = MIN (pbi->index_l_triangular_max[index_l1][index_p], pbi->index_l_triangular_max[index_l3][index_l]);
      int q_size = index_q_max - index_q_min + 1;

      /* If there are no nodes in the q range, we skip this p value; how would be able
      to interpolate the bispectrum in q otherwise? */
      if (q_size <= 0)
        continue;

      // ---------------------------------------------------------------------------------
      // -                                   3J and 6J                                   -
      // ---------------------------------------------------------------------------------

      double threej_q[2*pbi->l_max+1]; /* will contain the 3j symbol for all q */
      int q_min_3j, q_max_3j; /* limits in q of the 3j symbol */
      double sixj_q[2*pbi->l_max+1]; /* will contain the 6j symbol for all q */
      int q_min_6j, q_max_6j; /* limits in q of the 6j symbol */
      double threej_000[2*pbi->l_max+1]; /* will contain 3j symbol to convert the reduced bispectrum */
      int q_min_3j_000, q_max_3j_000; /* limits in q for the conversion 3j */


      /* Compute 3j symbol for all values of q */

      class_call (drc3jj (
                    /*q,*/ l, l3, 0, -F_X3,
                    &min_D, &max_D,
                    threej_q,
                    (2*pbi->l_max+1),
                    pbi->error_message),
        pbi->error_message,
        pbi->error_message);
      q_min_3j = (int)(min_D + _EPS_);
      q_max_3j = (int)(max_D + _EPS_);

      for (int q=q_min_3j; q <= q_max_3j; ++q)
        threej_q[q-q_min_3j] *= sqrt((2*l3+1.0)*(2*l+1.0)*(2*q+1.0)/(4*_PI_)) * 0.5 * (-l3*(l3+1.0)+l*(l+1.0)+q*(q+1.0));

      class_test ((q_min_3j != abs(l3-l) || (q_max_3j != (l3+l))),
        pbi->error_message,
        "wrong 3j limits");

      /* Compute the 6j symbol for all values of q */

      class_call (drc6j (
                    /*q,*/ p, l1, l2, l3, l,
                    &min_D, &max_D,
                    sixj_q,
                    (2*pbi->l_max+1),
                    pbi->error_message
                    ),
        pbi->error_message,
        pbi->error_message);
      q_min_6j = (int)(min_D + _EPS_);
      q_max_6j = (int)(max_D + _EPS_);

      class_test ((q_min != MAX (q_min_6j, 2)) || (q_max != (MIN (q_max_6j, pbi->l_max))),
        pbi->error_message,
        "wrong 3j limits");

      /* Before computing the conversion 3j, check whether we have already computed it
      in all the required q values for this p. If so, we shally recycle it. */
      int is_filled = 1;
      for (int q=q_min; q <= q_max; ++q)
        if (!(is_filled = is_filled_pq[index_p][q-2]))
          break;
      
      if (!is_filled) {

        class_call (drc3jj (
                      l1, p, 0, 0,
                      &min_D, &max_D,
                      threej_000,
                      (2*pbi->l_max+1),
                      pbi->error_message),
          pbi->error_message,
          pbi->error_message);
        q_min_3j_000 = (int)(min_D + _EPS_);
        q_max_3j_000 = (int)(max_D + _EPS_);

        for (int q=q_min_3j_000; q <= q_max_3j_000; ++q)
          threej_000[q-q_min_3j_000] *= sqrt((2*l1+1.0)*(2*p+1.0)*(2*q+1.0)/(4*_PI_));
      }


      // ======================================================================================
      // =                                      Loop on q                                     =
      // ======================================================================================

      for (int q=q_min; q <= q_max; ++q) {

        class_test (!is_triangular_int (l3,q,l) || !is_triangular_int (l1,p,q),
           pbi->error_message,
           "q not triangular");

        int X2_ = X2;
        int X3_ = X3;

        if (!is_filled_pq[index_p][q-2]) {

          double reduced_bispectrum;

          class_call (bispectra_at_l3_linear (
                        ptr,
                        pbi,
                        index_bt,
                        index_l1, index_p, q,
                        X1, X2_, X3_,
                        _TRUE_,
                        &reduced_bispectrum,
                        NULL),
            pbi->error_message,
            pbi->error_message);

          /* Compute the actual bispectrum from the reduced one */
          bispectrum_pq[index_p][q-2] = threej_000[q-q_min_3j_000] * reduced_bispectrum;

          is_filled_pq[index_p][q-2] = 1;

        }

        double integrand = C_PP *
                           ALTERNATING_SIGN (l1+l2+p) * threej_p[p-p_min_3j] * threej_q[q-q_min_3j] * sixj_q[q-q_min_6j] *
                           delta_l[index_l] * delta_p[index_p-index_p_min] *
                           bispectrum_pq[index_p][q-2];

        class_test (isnan (integrand),
          pbi->error_message,
          "integrand nan for l1=%d,l2=%d,l3=%d,l=%d,p=%d,q=%d; b=%g, 3j(p)=%g, 3j(q)=%g, 6j(q)=%g",
          l1,l2,l3,l,p,q, bispectrum_pq[index_p][q-2], threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);

        *result += integrand;

        /* Debug: save the integrand function */
        // integrand_lp[index_l][index_p] += integrand;

        /* Debug: print the integrand as a function of q */
        // if (print_output && (l==100) && (p==100)) {
        //   fprintf (stderr, "%6d %17g %17g %17g %17g %17g\n", q, integrand, bispectrum_pq[index_p][q-2],
        //     threej_p[p-p_min_3j], threej_q[q-q_min_3j], sixj_q[q-q_min_6j]);
        //   if (q == q_max) fprintf (stderr, "\n\n");
        // }

      } // for q

      /* Debug: print the integrand as a function of p */
      // if (print_output && (l==28)) {
      //   fprintf (stderr, "%6d %17g %17g\n", p, integrand_lp[index_l][index_p], threej_p[p-p_min_3j]);
      //   if (p == p_max) fprintf (stderr, "\n\n");
      // }

      /* Debug: save the integrand function */
      // integrand_l[index_l] += integrand_lp[index_l][index_p];

    } // for index_p

    /* Debug: print the integrand as a function of l */
    // if (print_output) {
    //   fprintf (stderr, "%6d %17g\n", l, integrand_l[index_l]);
    //   if (l == pbi->l_max) fprintf (stderr, "\n\n");
    // }


  } // for index_l

  free (bispectrum_pq);
  free (is_filled_pq);

  /* Debug: free the debug arrays */
  // free (integrand_lp);
  // free (integrand_l);

  return _SUCCESS_;

}




/**
 * Load a bispectrum from disk.
 *
 * See the documentation in perturbations2.h (\ref StorageFiles) for more
 * details.
 */

int bispectra_load (
      struct precision * ppr,
      struct bispectra * pbi,
      int index_bt
      )
{

  /* Load only if needed */
  if (pbi->bispectra_available[index_bt])
    return _SUCCESS_;

  /* If needed, it must be on disk */
  class_test (!(ppr->store_bispectra || ppr->load_bispectra),
    pbi->error_message,
    "shouldn't be here");

  /* Load only the bispectra that require a lot of time to compute,
  that is, the non-separable and intrinsic bispectra */
  if (!(pbi->bispectrum_type[index_bt] == non_separable_bispectrum ||
        pbi->bispectrum_type[index_bt] == intrinsic_bispectrum))
    return _SUCCESS_;

  /* Print some info */
  printf_log_if (pbi->bispectra_verbose, 2, 
    "     * reading bispectra from disk for index_bt=%d on'%s'\n",
    index_bt, pbi->storage_paths[index_bt]);

  /* Complain if there is no file to load */
  struct stat st;
  class_test (stat (pbi->storage_paths[index_bt], &st) != 0,
		pbi->error_message,
		"cannot load bispectra for index_bt=%d, file '%s' does not exist",
    index_bt, pbi->storage_paths[index_bt]);

  /* Make space */
  class_call (bispectra_allocate_type_level(pbi, index_bt),
    pbi->error_message,
    pbi->error_message);

  /* Open file for reading */
  class_open (pbi->storage_files[index_bt],
    pbi->storage_paths[index_bt],
    "rb", pbi->error_message);

  /* Read from file */
  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {

        int n_to_read = pbi->n_independent_configurations;
  
        /* Read a chunk with all the independent (l1,l2,l3) triplets to pbi->bispectra[index_bt] */
        int n_read = fread(
                pbi->bispectra[index_bt][X][Y][Z],
                sizeof(double),
                n_to_read,
                pbi->storage_files[index_bt]);

        class_test(n_read != n_to_read,
          pbi->error_message,
          "Could not read in '%s' file, read %d entries but expected %d",
            pbi->storage_paths[index_bt], n_read, n_to_read);        
      }
    }
  }
  
  /* Close file */
  fclose(pbi->storage_files[index_bt]); 

  /* Bispectrum is now ready to use */
  pbi->bispectra_available[index_bt] = _TRUE_;

  return _SUCCESS_;
  
}



/**
 * Bispectrum produced by the pi_dot * grad_pi^2 term in the Galileon Lagrangian, 
 * taken from eq. 22 of http://arxiv.org/abs/0905.3746.
 */
int bispectra_galileon_gradient (
      struct primordial * ppm,
      struct spectra * psp,
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
      struct spectra * psp,
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
      struct spectra * psp,
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
      struct spectra * psp,
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
      struct spectra * psp,
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
     int lens_me,
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
  
  if (pbi->has_bispectra_t) {
  
    /* Temperature-lensing potential correlation */
    double C_l1_Tp = pbi->cls[psp->index_ct_tp][l1-2];
    double C_l2_Tp = pbi->cls[psp->index_ct_tp][l2-2];
    double C_l3_Tp = pbi->cls[psp->index_ct_tp][l3-2];
                        
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
    double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
                  
    /* Use lensed temperature C_l's if available */
    if (lens_me) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
      C_l3 = pbi->lensed_cls[ple->index_lt_tt][l3-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l2*(l2+1.0) + l3*(l3+1.0) - l1*(l1+1.0) ) * C_l2_Tp * C_l3
      + ( l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0) ) * C_l3_Tp * C_l2
      + ( l1*(l1+1.0) + l3*(l3+1.0) - l2*(l2+1.0) ) * C_l1_Tp * C_l3
      + ( l3*(l3+1.0) + l1*(l1+1.0) - l2*(l2+1.0) ) * C_l3_Tp * C_l1
      + ( l1*(l1+1.0) + l2*(l2+1.0) - l3*(l3+1.0) ) * C_l1_Tp * C_l2
      + ( l2*(l2+1.0) + l1*(l1+1.0) - l3*(l3+1.0) ) * C_l2_Tp * C_l1
    );
  
    /* Debug - print temperature-lensing potential C_l's */
    // if ((l1==pbi->l_max) && (l2==pbi->l_max))
    //   fprintf (stderr, "%4d %10g %10g %10g\n", l3, C_l3_TT, C_l3, C_l3_Tp);
  
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // temperature only

  
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
  // if (pbi->has_bispectra_t) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_t) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_e) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  //   if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
  // }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b,
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
  if (lens_me) {
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
  double F_l1_l2_l3_X1 = 0.25 * ( l2*(l2+1.0) + l3*(l3+1.0) - l1*(l1+1.0) ) * S_X1 * threej_l1_l2_l3_FX1_0_mFX1; /* 1-2-3 */
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l1_l3_X2 = 0.25 * ( l1*(l1+1.0) + l3*(l3+1.0) - l2*(l2+1.0) ) * S_X2 * threej_l2_l1_l3_FX2_0_mFX2; /* 2-1-3 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1.0) + l1*(l1+1.0) - l2*(l2+1.0) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
  double F_l3_l1_l2_X3 = 0.25 * ( l1*(l1+1.0) + l2*(l2+1.0) - l3*(l3+1.0) ) * S_X3 * threej_l3_l1_l2_FX3_0_mFX3; /* 3-1-2 */
  double F_l3_l2_l1_X3 = 0.25 * ( l2*(l2+1.0) + l1*(l1+1.0) - l3*(l3+1.0) ) * S_X3 * threej_l3_l2_l1_FX3_0_mFX3; /* 3-2-1 */
  
  
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
  if ((pbi->has_bispectra_t) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
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
 * Compute the CMB lensing bispectrum kernel including polarisation in the squeezed
 * limit (l3<<l1 and l3<<l2).
 *
 * Here we code just the kernel in Eq. 5.20 (eq. 45 arXiv) of Lewis, Challinor &
 * Hanson 2011 (http://uk.arxiv.org/abs/1101.2234); the full squeezed bispectrum
 * is given by the product between the kernel and C_l3^{X3\phi}.
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
 */

int bispectra_cmb_lensing_squeezed_kernel (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     int lens_me,
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
  
  if (pbi->has_bispectra_t) {
  
    /* By default take unlensed temperature C_l's */
    double C_l1 = pbi->cls[psp->index_ct_tt][l1-2];
    double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
                  
    /* Use lensed temperature C_l's if available */
    if (lens_me) {
      C_l1 = pbi->lensed_cls[ple->index_lt_tt][l1-2];
      C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
    }
  
    /* CMB lensing bispectrum formula for TTT */
    ttt = 0.5 * (
      + ( l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0) ) * C_l2
      + ( l3*(l3+1.0) + l1*(l1+1.0) - l2*(l2+1.0) ) * C_l1
    );

    /* Uncomment to compute the bispectrum rather than the kernel */
    // ttt *= pbi->cls[psp->index_ct_tp][l3-2];                      
    
    /* If only temperature is requested, skip what follows exit from the 'cmb_lensing' if block */
    if (pbi->bf_size == 1) {
      *result = ttt;
      return _SUCCESS_;
    }
    
  } // temperature only
  
  
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
  // if (pbi->has_bispectra_t) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_e) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  // }
  
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b,
    pbi->error_message,
    "CMB-lensing squeezed bispectrum for B-modes not implemented yet.");
  
  // -------------------------------------------------------------------------------
  // -                              Obtain the C_l                                 -
  // -------------------------------------------------------------------------------
  
  /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
  double C_l2_X1_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X1 ][ X2 ]][l2-2];
  double C_l1_X2_X1 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X1 ]][l1-2];
    
  /* Use lensed temperature C_l's if available */
  if (lens_me) {
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
  double F_l1_l3_l2_X1 = 0.25 * ( l3*(l3+1.0) + l2*(l2+1.0) - l1*(l1+1.0) ) * S_X1 * threej_l1_l3_l2_FX1_0_mFX1; /* 1-3-2 */
  double F_l2_l3_l1_X2 = 0.25 * ( l3*(l3+1.0) + l1*(l1+1.0) - l2*(l2+1.0) ) * S_X2 * threej_l2_l3_l1_FX2_0_mFX2; /* 2-3-1 */
      

  // -------------------------------------------------------------------------------
  // -                              Bispectrum formula                             -
  // -------------------------------------------------------------------------------
  
  /* Kernel of the CMB-lensing bispectrum in the squeezed limit, from Eq. 5.20 of Lewis,
  Challinor & Hanson 2011. This is simply the general formula with C_l1_X1_p=C_l2_X2_p=0 
  (see bispectra_cmb_lensing_bispectrum()) */
  *result = C_l2_X1_X2 * F_l1_l3_l2_X1
          + C_l1_X2_X1 * F_l2_l3_l1_X2;
                  
  /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
  if ((pbi->has_bispectra_t) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
    double exact = ttt;
    double diff = fabs (1-*result/exact);
    class_test (diff > _SMALL_,
     pbi->error_message,
     "CMB-lensing squeezed bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
     l1, l2, l3, *result, exact, diff);
  }

  return _SUCCESS_;
  
}





// int bispectra_cmb_lensing_squeezed_kernel (
//      struct precision * ppr,
//      struct spectra * psp,
//      struct lensing * ple,
//      struct bispectra * pbi,
//      int l1, int l2, int l3,
//      int X1, int X2, int X3,
//      int lens_me,
//      double threej_ratio_20m2,
//      double threej_ratio_m220,
//      double threej_ratio_0m22,
//      double * result
//      )
// {
//
//   /* Test that l1 is the smallest multipole */
//   class_test ((l3<l1) || (l2<l1),
//     pbi->error_message,
//     "(l1=%d,l2=%d,l3=%d) :in all squeezed approximations, make sure l1 is the smallest multipole",
//     l1, l2, l3);
//
//   // --------------------------------------------------------------------------------------
//   // -                              Temperature-only formula                              -
//   // --------------------------------------------------------------------------------------
//
//   /* When only including temperature, the CMB-lensing bispectrum reduces to
//     b^TTT_l1l2l3 = 1/2 * [l1(l1+1) + l2(l2+1) - l3(l3+1)] C_l1^{TT} C_l2^{T\psi} + 5 perm,
//   where \psi is the lensing potential. */
//
//   double ttt;
//
//   if (pbi->has_bispectra_t) {
//
//     /* By default take unlensed temperature C_l's */
//     double C_l2 = pbi->cls[psp->index_ct_tt][l2-2];
//     double C_l3 = pbi->cls[psp->index_ct_tt][l3-2];
//
//     /* Use lensed temperature C_l's if available */
//     if (lens_me) {
//       C_l2 = pbi->lensed_cls[ple->index_lt_tt][l2-2];
//       C_l3 = pbi->lensed_cls[ple->index_lt_tt][l3-2];
//     }
//
//     /* CMB lensing bispectrum formula for TTT */
//     ttt = 0.5 * (
//       + ( l1*(l1+1.0) + l2*(l2+1.0) - l3*(l3+1.0) ) * C_l2
//       + ( l1*(l1+1.0) + l3*(l3+1.0) - l2*(l2+1.0) ) * C_l3
//     );
//
//     /* Uncomment to compute the bispectrum rather than the kernel */
//     // ttt *= pbi->cls[psp->index_ct_tp][l1-2];
//
//     /* If only temperature is requested, skip what follows exit from the cmb_lensing
//     if block */
//     if (pbi->bf_size == 1) {
//       *result = ttt;
//       return _SUCCESS_;
//     }
//
//   } // temperature only
//
//
//   // -------------------------------------------------------------------------------
//   // -                        Determine field coefficients                         -
//   // -------------------------------------------------------------------------------
//
//   /* Set the amplitude of the geometrical factor in Eq. 4.4 of Lewis et al. 2011.
//   This is either 2 (if the two 3j-symbols add up) or 0 (if they cancel). Whether
//   they add up or cancel depends on the parity of the considered field (even for
//   T and E and odd for B). */
//
//   double S_X2=2, S_X3=2;
//
//   /* Uncomment the following to set the bispectrum to zero based on the parity of
//   the involved fields. By default we skip this part, because we want our reduced
//   bispectrum to be continuous in order to facilitate its interpolation and
//   plotting. */
//
//   // int L = l1-l2-l3;
//   //
//   // if (pbi->has_bispectra_t) {
//   //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?2:0;
//   //   if (X3 == pbi->index_bf_t) S_X3 = (L%2==0)?2:0;
//   // }
//   // if (pbi->has_bispectra_e) {
//   //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
//   //   if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
//   // }
//   // if (pbi->has_bispectra_b) {
//   //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
//   //   if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
//   // }
//
//   /* B-mode bispectrum not implemented yet */
//   class_test (pbi->has_bispectra_b,
//     pbi->error_message,
//     "CMB-lensing squeezed bispectrum for B-modes not implemented yet.");
//
//   // -------------------------------------------------------------------------------
//   // -                              Obtain the C_l                                 -
//   // -------------------------------------------------------------------------------
//
//   /* Get the C_l's involving the fields. By default take unlensed temperature C_l's */
//   double C_l2_X3_X2 = pbi->cls[pbi->index_ct_of_bf_bf[ X3 ][ X2 ]][l2-2];
//   double C_l3_X2_X3 = pbi->cls[pbi->index_ct_of_bf_bf[ X2 ][ X3 ]][l3-2];
//
//   /* Use lensed temperature C_l's if available */
//   if (lens_me) {
//     C_l2_X3_X2 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X3 ][ X2 ]][l2-2];
//     C_l3_X2_X3 = pbi->lensed_cls[pbi->index_lt_of_bf_bf[ X2 ][ X3 ]][l3-2];
//   }
//
//   // -------------------------------------------------------------------------------
//   // -                             Get the right 3j                                -
//   // -------------------------------------------------------------------------------
//
//   /* Spin of the fields */
//   int F_X3 = pbi->field_spin[X3];
//   int F_X2 = pbi->field_spin[X2];
//
//   /* Obtain the needed 3j ratios (TODO: what about B-modes?) */
//   double threej_l3_l2_l1_FX3_0_mFX3 = 1;
//   double threej_l3_l1_l2_FX3_0_mFX3 = 1;
//   if (F_X3==2) {
//     threej_l3_l2_l1_FX3_0_mFX3 = threej_ratio_20m2;
//     threej_l3_l1_l2_FX3_0_mFX3 = threej_ratio_m220;
//   }
//
//   double threej_l2_l1_l3_FX2_0_mFX2 = 1;
//   double threej_l2_l3_l1_FX2_0_mFX2 = 1;
//   if (F_X2==2) {
//     threej_l2_l1_l3_FX2_0_mFX2 = threej_ratio_m220;
//     threej_l2_l3_l1_FX2_0_mFX2 = threej_ratio_0m22;
//   }
//
//   /* Obtain the geometric factor F^+s_l3l2l1 */
//   double F_l3_l1_l2_X3 = 0.25 * ( l1*(l1+1.0) + l2*(l2+1.0) - l3*(l3+1.0) ) * S_X3 * threej_l3_l1_l2_FX3_0_mFX3; /* 1-3-2 */
//   double F_l2_l1_l3_X2 = 0.25 * ( l1*(l1+1.0) + l3*(l3+1.0) - l2*(l2+1.0) ) * S_X2 * threej_l2_l1_l3_FX2_0_mFX2; /* 2-3-1 */
//
//
//   // -------------------------------------------------------------------------------
//   // -                              Bispectrum formula                             -
//   // -------------------------------------------------------------------------------
//
//   /* Kernel of the CMB-lensing bispectrum in the squeezed limit, from Eq. 5.20 of Lewis,
//   Challinor & Hanson 2011. This is simply the general formula with C_l3_X3_p=C_l2_X2_p=0
//   (see bispectra_cmb_lensing_bispectrum()) */
//   *result = C_l2_X3_X2 * F_l3_l1_l2_X3
//           + C_l3_X2_X3 * F_l2_l1_l3_X2;
//
//   /* Check that for <TTT> the bispectrum is equal to the one computed with the simpler formula */
//   if ((pbi->has_bispectra_t) && (X1==pbi->index_bf_t) && (X2==X3) && (X1==X3)) {
//     double exact = ttt;
//     double diff = fabs (1-*result/exact);
//     class_test (diff > _SMALL_,
//      pbi->error_message,
//      "CMB-lensing squeezed bispectrum for TTT does not reduce to simple formula; l=(%d,%d,%d), b=%g, exact=%g, diff=%g",
//      l1, l2, l3, *result, exact, diff);
//   }
//
//   return _SUCCESS_;
//
// }




/*
 * Compute the CMB lensing bispectrum including polarisation in the squeezed limit, valid
 * when l3<<l1 and l3<<l2. This is given by the kernel computed in bispectra_cmb_lensing_squeezed_kernel()
 * times the cross-correlation with the lensing potential, C_l3^{X3\phi} (see Eq. 5.20 of Lewis,
 * Challinor & Hanson 2011 (http://uk.arxiv.org/abs/1101.2234).
 *
 * This formula is non-perturbative (in the sense of Sec. 3.2, ibidem) and is valid only
 * for squeezed configurations, where l3<<l1 and l3<<l2; it is obtained from the general
 * formula (Eq. 4.5, ibidem) by setting C_l1_X1_p=C_l2_X2_p=0. It is an excellent approximation
 * for the CMB-lensing bispectrum, as most of its signal is in these squeezed configurations.
 *
 * Look at the documentation of above bispectra_intrinsic_squeezed_bispectrum() for
 * details about the symmetry properties of this special squeezed bispectrum.
 * 
 */
int bispectra_cmb_lensing_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     int lens_me,
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
                lens_me,
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
     int lens_me,
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
  if (lens_me) {
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
 * Squeezed-limit approximation for the intrinsic bispectrum, as in eq. 7 of
 * Pettinari et al. http://arxiv.org/abs/1406.2981.
 *
 * This is the general formula that includes polarisation, which was first
 * derived in eq. 4.1 and 4.2 of Lewis 2012 (http://arxiv.org/abs/1204.5018).
 * With respect to that reference, we switch (i,l1)->(Z,l3), (j,l2)->(X,l1),
 * (k,l3)->(Y,l2).
 *
 * Note that the lensed version of this bispectrum is considerably smaller than
 * the unlensed one, due to the effect explained below eq. 7 of Pettinari et al.
 * (http://arxiv.org/abs/1406.2981).
 *
 * See also:
 * - Creminelli et al. 2004 (http://arxiv.org/abs/astro-ph/0405428)
 * - Creminelli et al. 2011 (http://arxiv.org/abs/1109.1822)
 * - Bartolo et al. 2012 (http://arxiv.org/abs/1109.2043).
 * 
 * CONSIDERATIONS THAT APPLY TO ALL "SQUEEZED" BISPECTRA
 *
 * In this and in the other "squeezed" approximation functions, l3 is taken to be the
 * long-wavelength mode and l1 and l2 the short ones. This ordering is preferred because
 * in SONG we often loop over (l1,l2,l3) configurations that satisfy the condition
 * l1>=l2>=l3.
 *
 * Therefore, it is crucial that that, in the formulas below, the multipole l3 is
 * associated with the C_l that correlates with the comoving curvature perturbation zeta
 * (for the CMB lensing bispectrum it would be the lensing potential phi) because l3 is
 * the smallest multipole, which describes the long wavelength mode. Also, l3 must be
 * associated with Z (or X3) because our convention for the bispectrum is <X_l1 Y_l2 Z_l3>.
 * If you give l3 to another field, then the Fisher matrix estimator will associate to that
 * field the wrong covariance matrix, and the result will change drastically.
 *
 * Because of the special role played by l3, the bispectrum computed in this and the
 * other 'squeezed' approximation functions is NOT symmetric with respect to an exchange
 * of (l1,X) <-> (l3,Z) or of or (l2,Y) <-> (l3,Z), contrary to the other bispectra,
 * which are computed as <X_l1 Y_l2 Z_l3>. Therefore, for these bispectra one cannot
 * obtain the configurations outside l1>=l2>=l3 by permuting the XYZ indices, as it is
 * done in the bispectra_at_node() function.
 */

int bispectra_intrinsic_squeezed_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     int lens_me,
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
  // if ((pbi->has_bispectra_e) && (Z==pbi->index_bf_e)) {
  //   *result = 0;
  //   return _SUCCESS_;
  // }

  /* We take l3 to be the long wavelength and l1 and l2 the short ones */
  double cl3_Zz = pbi->cls[pbi->index_ct_of_zeta_bf[ Z ]][l3-2];
  double dcl1_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l1-2];
  double dcl2_XY = pbi->d_lsq_cls[pbi->index_ct_of_bf_bf[X][Y]][l2-2];
  
  /* Use lensed temperature C_l's if available */
  if (lens_me) {
    dcl1_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l1-2];    
    dcl2_XY = pbi->lensed_d_lsq_cls[pbi->index_lt_of_bf_bf[X][Y]][l2-2];    
  }
          
  /* Ricci focussing in Lewis 2012 (eq. 4.1) */
  double bolometric_T_lewis_ricci = - 0.5 * cl3_Zz * (dcl1_XY/l1 + dcl2_XY/l2);

  /* Redshift modulation in Lewis 2012 (eq. 4.2). This exists only if Y=Z=temperature */
  double bolometric_T_lewis_redshift = 0;
          
  if (pbi->has_bispectra_t) {
    double cl3_Zt_long = 0; double cl1_Xt_short = 0; double cl2_Yt_short = 0;
    cl3_Zt_long = pbi->cls[pbi->index_ct_of_bf_bf[Z][pbi->index_bf_t]][l3-2];
    if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->cls[pbi->index_ct_of_bf_bf[X][pbi->index_bf_t]][l1-2];
    if (X == pbi->index_bf_t) cl2_Yt_short = pbi->cls[pbi->index_ct_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    if (lens_me) {
      if (Y == pbi->index_bf_t) cl1_Xt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[X][pbi->index_bf_t]][l1-2];
      if (X == pbi->index_bf_t) cl2_Yt_short = pbi->lensed_cls[pbi->index_lt_of_bf_bf[Y][pbi->index_bf_t]][l2-2];
    }
    bolometric_T_lewis_redshift = cl3_Zt_long * (cl1_Xt_short + cl2_Yt_short);
  }
          
  class_test ((pbi->has_bispectra_e) &&
    ((X == pbi->index_bf_e) && (Y == pbi->index_bf_e))
    && (bolometric_T_lewis_redshift!=0),
    pbi->error_message,
    "anisotropic redshifting should vanish (eq. 4.2 of Lewis 2012)", bolometric_T_lewis_redshift);
          
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
     int lens_me,
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
  
  if (pbi->has_bispectra_t) {
  
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

  if (pbi->has_bispectra_t) {
    if (X1 == pbi->index_bf_t) {S_X1 = 1; T_X1 = pbi->index_bf_t; index_ct_X1_I = psp->index_ct_tt;}
    if (X2 == pbi->index_bf_t) {S_X2 = 1; T_X2 = pbi->index_bf_t; index_ct_X2_I = psp->index_ct_tt;}
    if (X3 == pbi->index_bf_t) {S_X3 = 1; T_X3 = pbi->index_bf_t; index_ct_X3_I = psp->index_ct_tt;}
  }
  if (pbi->has_bispectra_e) {
    if (X1 == pbi->index_bf_e) {S_X1 = 2; T_X1 = pbi->index_bf_e; index_ct_X1_I = psp->index_ct_te;}
    if (X2 == pbi->index_bf_e) {S_X2 = 2; T_X2 = pbi->index_bf_e; index_ct_X2_I = psp->index_ct_te;}
    if (X3 == pbi->index_bf_e) {S_X3 = 2; T_X3 = pbi->index_bf_e; index_ct_X3_I = psp->index_ct_te;}
  }
  if (pbi->has_bispectra_b) { /* Note that <TB> vanishes, hence the negative values */
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
  // if (pbi->has_bispectra_t) {
  //   if (X1 == pbi->index_bf_t) S_X1 = (L%2==0)?1:0;
  //   if (X2 == pbi->index_bf_t) S_X2 = (L%2==0)?1:0;
  //   if (X3 == pbi->index_bf_t) S_X3 = (L%2==0)?1:0;
  // }
  // if (pbi->has_bispectra_e) {
  //   if (X1 == pbi->index_bf_e) S_X1 = (L%2==0)?2:0;
  //   if (X2 == pbi->index_bf_e) S_X2 = (L%2==0)?2:0;
  //   if (X3 == pbi->index_bf_e) S_X3 = (L%2==0)?2:0;
  // }
  // if (pbi->has_bispectra_b) {
  //   if (X1 == pbi->index_bf_b) S_X1 = (L%2!=0)?2:0;
  //   if (X2 == pbi->index_bf_b) S_X2 = (L%2!=0)?2:0;
  //   if (X3 == pbi->index_bf_b) S_X3 = (L%2!=0)?2:0;
  // }
              
  /* B-mode bispectrum not implemented yet */
  class_test (pbi->has_bispectra_b,
    pbi->error_message,
    "quadratic correction for B-mode bispectrum not implemented yet.");
    
  
  // -------------------------------------------------------------------------------
  // -                               Obtain the C_l                                -
  // -------------------------------------------------------------------------------
  
  /* Get the C_l's. When implementing the B-modes, remember to 
  set the C_l's to zero, so that the only contribution comes from the
  bispectrum with the second-order a^B_lm.  */
    
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
  
  /* Obtain the needed 3j ratios (TODO: what about B-modes?) */
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
  *result = 4.0 * S_X1 * (threej_l2_l3_l1_0_FX1_mFX1*C_l2_X2_I*C_l3_X3_TX1 + threej_l2_l3_l1_FX1_0_mFX1*C_l2_X2_TX1*C_l3_X3_I)
          + 4.0 * S_X2 * (threej_l3_l1_l2_0_FX2_mFX2*C_l3_X3_I*C_l1_X1_TX2 + threej_l3_l1_l2_FX2_0_mFX2*C_l3_X3_TX2*C_l1_X1_I)
          + 4.0 * S_X3 * (threej_l1_l2_l3_0_FX3_mFX3*C_l1_X1_I*C_l2_X2_TX3 + threej_l1_l2_l3_FX3_0_mFX3*C_l1_X1_TX3*C_l2_X2_I);

  /* Check that for <TTT> the correction is equal to 8 * (C_l1*C_l2 + C_l1*C_l3 + C_l2*C_l3) */
  if ((pbi->has_bispectra_t) && (X1==pbi->index_bf_t) && (X1==X2) && (X1==X3)) {
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
 * Simple C_l x C_l bispectrum used for testing purposes.
 *
 * This bispectrum has to be symmetric in (l1,X1) <-> (l2,X2) <-> (l3,X3)
 * and to involve quadratic products of the angular power spectrum C_l.
 */

int bispectra_test_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     int lens_me,
     double threej_ratio_20m2,
     double threej_ratio_m220,
     double threej_ratio_0m22,
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

  /* Build the bispectrum */
  *result =   6/5. * cl1_Xz * (cl2_YY + cl3_ZZ)
            + 6/5. * cl2_Yz * (cl3_ZZ + cl1_XX)
            + 6/5. * cl3_Zz * (cl1_XX + cl2_YY);  
  
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
 * Simple bispectrum used for testing purposes, given by the squeezed 
 * limit of the local bispectrum modulated by a cosine function.
 *
 * The bispectrum reads:
 *
 *   b_l1l2l3 = 6/5 * cl3_Zz * (cl1_XY + cl2_XY) * cos((l1+l2+l3)/50).
 *
 */
int bispectra_cosine_bispectrum (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X1, int X2, int X3,
     int lens_me,
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
                lens_me,
                threej_ratio_20m2,
                threej_ratio_m220,
                threej_ratio_0m22,
                result),
    pbi->error_message,
    pbi->error_message);

  *result *= cos((l1+l2+l3)/50.0);

  return _SUCCESS_;
  
}




/** 
 * Window function for the local (l1,l2,l3) bispectrum (DISABLED).
 *
 * The natural scaling for the bispectrum (both the templates and the second-order
 * one) is given by a product of two power spectra. When available, we always use
 * the C_l's for the temperature as, when multiplied by l^2, they are approximately
 * the same order of magnitude for all l's, contrary to those for polarisation which
 * are very small for l<200.
 */

int bispectra_local_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double * result
     )
{

  if (pbi->has_bispectra_t) {
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
  // if ((pbi->has_bispectra_t) && (X==pbi->index_bf_t) && (X==Y) && (X==Z)) {
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
 * Window function for the intrinsic (l1,l2,l3) bispectrum (DISABLED).
 */

int bispectra_intrinsic_window_function (
     struct precision * ppr,
     struct spectra * psp,
     struct lensing * ple,
     struct bispectra * pbi,
     int l1, int l2, int l3,
     int X, int Y, int Z,
     double * result
     )
{

  if (pbi->has_bispectra_t) {
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
 * Interpolate the bispectrum in the (l1,l2,l3) configuration using 2D mesh
 * interpolation; l1 must belong to pbi->l.
 *
 * This function will interpolate the bispectrum corresponding to index_bt for
 * the fields XYZ (eg. TTT, EEE, EET, ..., EEE).
 *
 * The multipoles must satisfy the condition l1<=l2<=l3<=pbi->l_max. To obtain
 * the bispectrum in a generic configuration, use bispectra_at_l2l3_bilinear_mesh() 
 * instead.
 *
 * If you asked for a window function (ie. pbi->window_function[index_bt] != NULL),
 * remember to multiply the returned value by the window function. We do not do it
 * inside the interpolate function for optimization purposes.
 */
 
int bispectra_mesh_interpolate (
      struct precision * ppr,
      struct spectra * psp,
      struct lensing * ple,
      struct bispectra * pbi,
      int index_bt,
      int index_L1, /**< Input: l1 multipole for which the interpolation will be performed; must be the
                    smallest of the three multipoles; negative ignored for 3D mesh interpolation */
      double l1, double l2, double l3,
      int X, int Y, int Z,
      struct interpolation_mesh * fine_mesh, /**< Input: mesh table to use for the small multipoles */
      struct interpolation_mesh * coarse_mesh, /**< Input: mesh table to use for the large multipoles */
      double * result /** Output: value of the bispectrum in (l1,l2,l3) */
      )
{

#ifdef DEBUG
  class_test ((index_L1>=0) && (int)(l1+_EPS_)!=pbi->l[index_L1],
    pbi->error_message,
    "inconsistent input: (l1=%g) != (pbi->l[%d]=%d)", l1, index_L1, pbi->l[index_L1]);

  class_test (pbi->mesh_restrict_l1l2l3 && (l3 < l2 || l2 < l1),
    pbi->error_message,
    "mesh was computed only for l1<=l2<=l3");
#endif // DEBUG


  // ====================================================================================
  // =                                   2D or 3D?                                      =
  // ====================================================================================

  /* Determine whether the user wants to compute a 2D or a 3D mesh */
  short is_2D = (index_L1 >= 0);
  double x, y, z;
  int use_fine_mesh;
  
  if (is_2D) {
    use_fine_mesh = (l2<pbi->mesh_l_turnover) && (l3<pbi->mesh_l_turnover);
    x = l2;
    y = l3;
    z = 0;
  }
  else {
    use_fine_mesh = (l1<pbi->mesh_l_turnover) && (l2<pbi->mesh_l_turnover) && (l3<pbi->mesh_l_turnover);
    x = l1;
    y = l2;
    z = l3;
  }


  // ====================================================================================
  // =                                Fine or coarse?                                   =
  // ====================================================================================

  /* Choose mesh parameters depending on whether we use the fine or coarse grid */
  int l_max;
  double link_length;
  double group_length;
  double soft_coeff;
  struct interpolation_mesh * mesh;
  int index_mesh;

  if (use_fine_mesh) {
    
    l_max = pbi->mesh_l_turnover;
    link_length = pbi->mesh_link_lengths[0];
    group_length = pbi->mesh_group_lengths[0];
    soft_coeff = pbi->mesh_soft_coeffs[0];
    mesh = fine_mesh;
    index_mesh = 0;

  }
  
  else {
    
    l_max = pbi->l_max;
    link_length = pbi->mesh_link_lengths[1];
    group_length = pbi->mesh_group_lengths[1];
    soft_coeff = pbi->mesh_soft_coeffs[1];
    mesh = coarse_mesh;
    index_mesh = 1;

  }

  
  // ====================================================================================
  // =                                 Create mesh                                      =
  // ====================================================================================

  /* If the mesh is not ready, then create it. If you are calling this function
  inside a parallel loop in (l1,l2,l3), beware of race conditions! Two threads 
  might attempt to initialise the same mesh at the same time, thus leading to
  unexpected behaviour. For 2D mesh interpolation, this is not a problem as long
  as l1 is the only parallelised loop. For 3D mesh interpolation, this is always
  a problem, because the same mesh is used for all (l1,l2,l3) configurations. 
  Wrapping the next block of code in a #pragma omp critical directive fixes the
  isse, but slows down the interpolation considerably. */

  if (!mesh->ready) {

    /* Check for race conditions */

    #ifdef _OPENMP
    int abort = _FALSE_;
    #pragma omp parallel
    class_test_parallel (omp_get_num_threads()>1 && pbi->interpolation_method==mesh_interpolation_3D,
      pbi->error_message,
      "stopping to prevent race condition");
    if (abort)
      return _FAILURE_;
    #endif
    
    class_call (bispectra_mesh_create(
                  ppr, psp, ple, pbi,
                  index_bt,
                  X, Y, Z,
                  pbi->window_function[index_bt],
                  index_L1,
                  l_max,
                  link_length,
                  group_length,
                  soft_coeff,
                  NULL,
                  NULL,
                  mesh),
      pbi->error_message,
      pbi->error_message);
  }



  // =====================================================================================
  // =                                 Interpolate                                       =
  // =====================================================================================

  /* Interpolate the bispectrum in (l1,l2,l3) */

  class_call (mesh_interpolate (
                mesh,
                x,y,z,
                result),
    mesh->error_message,
    pbi->error_message);


#ifdef DEBUG
  /* Check for nan's and crazy values. A value is crazy when it is much larger than
  the characteristic scale for a bispectrum, A_s*A_s~1e-20 */
  class_warning (isnan(*result) || fabs(*result)>1,
    "Interpolated b(%g,%g,%g) = %g for bispectrum %s_%s!!!",
    l1, l2, l3, *result,
    pbi->bt_labels[index_bt],
    pbi->bfff_labels[X][Y][Z]);
#endif // DEBUG

  /* Debug: print the interpolated value of the intrinsic bispectrum */
  // if ((pbi->has_intrinsic) || (pbi->index_bt_intrinsic)) {
  //   printf ("%8g %8g %8g %16.7g\n", l1, l2, l3, *result);
  // }

  return _SUCCESS_;
  
}



/**
 * Allocate memory for an array of interpolation meshes with the following
 * structure: meshes[index_bt][X][Y][Z][2].
 */

int bispectra_mesh_allocate(
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct interpolation_mesh ******* meshes
        )
{

  /* Allocate one mesh interpolation workspace per type of bispectrum */
  class_alloc ((*meshes),
    pbi->bt_size*sizeof(struct interpolation_mesh *****),
    pbi->error_message);

  /* Allocate worspaces and intialize counters */
  for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {
    
    /* The analytical bispectra are never interpolated */
    if (!pbi->interpolate_me[index_bt])
      continue;

    class_alloc ((*meshes)[index_bt],
      pbi->bf_size*sizeof(struct interpolation_mesh ****),
      pbi->error_message);

    for (int X = 0; X < pbi->bf_size; ++X) {

      class_alloc ((*meshes)[index_bt][X],
        pbi->bf_size*sizeof(struct interpolation_mesh ***),
        pbi->error_message);

      for (int Y = 0; Y < pbi->bf_size; ++Y) {
        
        class_alloc ((*meshes)[index_bt][X][Y],
          pbi->bf_size*sizeof(struct interpolation_mesh **),
          pbi->error_message);

        for (int Z = 0; Z < pbi->bf_size; ++Z) {

          class_alloc ((*meshes)[index_bt][X][Y][Z],
            2*sizeof(struct interpolation_mesh *),
            pbi->error_message);
    
          for (int index_mesh=0; index_mesh < 2; ++index_mesh) {

            class_alloc ((*meshes)[index_bt][X][Y][Z][index_mesh],
              sizeof(struct interpolation_mesh),
              pbi->error_message);
            
            (*meshes)[index_bt][X][Y][Z][index_mesh]->ready = _FALSE_;
            
          } // for(index_mesh)
        } // for(Z)
      } // for(Y)
    } // for(X)
  } // for(index_bt)
  
  return _SUCCESS_;

}



/**
 * Deallocate the input array of interpolation meshes
 */

int bispectra_mesh_free(
      struct bispectra * pbi,
      struct interpolation_mesh ******* meshes
      )
{
  
  class_call (bispectra_mesh_empty (pbi, meshes),
    pbi->error_message,
    pbi->error_message);

  for (int index_bt=(pbi->bt_size-1); index_bt >= 0; --index_bt) {
  
    if (!pbi->interpolate_me[index_bt])
      continue;
  
    for (int X = (pbi->bf_size-1); X >= 0; --X) {
      for (int Y = (pbi->bf_size-1); Y >= 0; --Y) {
        for (int Z = (pbi->bf_size-1); Z >= 0; --Z) {          
          for (int index_mesh=0; index_mesh < 2; ++index_mesh) {
            free ((*meshes)[index_bt][X][Y][Z][index_mesh]);
          }
          free ((*meshes)[index_bt][X][Y][Z]);
        }
        free ((*meshes)[index_bt][X][Y]);
      }
      free ((*meshes)[index_bt][X]);
    }
    free ((*meshes)[index_bt]);
  }

  free (*meshes);
    
  return _SUCCESS_;

}




/**
 * Empty all grids and meshes in the input mesh structures
 */

int bispectra_mesh_empty(
      struct bispectra * pbi,
      struct interpolation_mesh ******* meshes
      )
{

  for (int index_bt=(pbi->bt_size-1); index_bt >= 0; --index_bt) {
  
    if (!pbi->interpolate_me[index_bt])
      continue;
  
    for (int X = (pbi->bf_size-1); X >= 0; --X)
      for (int Y = (pbi->bf_size-1); Y >= 0; --Y)
        for (int Z = (pbi->bf_size-1); Z >= 0; --Z)
          for (int index_mesh=0; index_mesh < 2; ++index_mesh)
            mesh_free ((*meshes)[index_bt][X][Y][Z][index_mesh]);

  }
  
  return _SUCCESS_;

}



/**
 * Initialise an interpolation mesh for the input bispectrum at l1.
 *
 * The mesh thus initialised can be used to interpolate the bispectrum in
 * the l1 bidimensional slice, via the function bispectra_at_l2l3_bilinear_mesh().
 *
 * If the l1 index provided is negative, the function will initialise a 3D
 * mesh, ie. it will consider all (l1,l2,l3) points rather than considering
 * only a 2D slice with fixed l1.
 * 
 * If mesh_restrict_l1l2l3 is true, only those configurations with l1<=l2<=l3
 * will be considered. This is an optimisation flag that will reduce the
 * memory usage for 2D interpolation. We choose l1<=l2<=l3 and not viceversa
 * because the two largest multipoles correspond to the slowest varying
 * directions in the bispectrum, and are therefore easier to interpolate.
 */

int bispectra_mesh_create (
      struct precision * ppr,
      struct spectra * psp,
      struct lensing * ple,
      struct bispectra * pbi,

      int index_bt, /**< Input: the index corresponding to the bispectrum to interpolate */
      int X, int Y, int Z,
      window_function_type * window_function, /**< Input: window function to apply to the bispectrum; set to NULL to ignore */
      int index_L1, /**< Input: multipole for which the 2D mesh will be computed; if negative, fill a 3D mesh instead */

      int l_max, /**< Input: the largest multipole to consider for this mesh */
      double link_length, /**< Input: linking length for this mesh (see bispectra.h) */
      double group_length, /**< Input: grouping length for this mesh (see bispectra.h) */
      double soft_coeff, /**< Input: softening of the linking length for this mesh (see bispectra.h) */
      int *** grid, /**< Input: interpolation grid for the mesh; set to NULL to recompute it */
      int **** id, /**< Input: id array for the mesh; set to NULL to recompute it */
      struct interpolation_mesh * mesh /**< Output: the interpolation mesh to initialise */
      )
{

  /* Determine whether the user wants to compute a 2D or a 3D mesh */
  short is_2D = (index_L1 >= 0);
  int n_dim = (is_2D?2:3);

  /* Count the number of nodes in the function to be interpolated. For 3D interpolation,
  this is equal to the total number of (l1,l2,l3) sampling points in the bispectrum, minus
  those with any of l1, l2 or l3 larger than l_max. For 2D interpolation, it is equal to
  the number of (l2,l3) sampling points in the considered l1-plane, minus those with either
  l2 or l3 larger than l_max. */
  long int n_nodes = 0;

  for(int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

    /* For 2D interpolation, we only need a single l1 value */
    if (is_2D && index_l1!=index_L1)
      continue;

    /* Restrict the interpolation to l1<=l2 if requested */
    int index_l2_min = 0;
    if (pbi->mesh_restrict_l1l2l3)
      index_l2_min = index_l1;

    for (int index_l2=index_l2_min; index_l2 < pbi->l_size; ++index_l2) {

      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];
      
      /* Restrict the interpolation to l2<=l3 if requested */
      if (pbi->mesh_restrict_l1l2l3)
        index_l3_min = MAX (index_l2, index_l3_min);
      
      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

        int l1 = pbi->l[index_l1];
        int l2 = pbi->l[index_l2];
        int l3 = pbi->l[index_l3];

        if (l2>l_max || l3>l_max || (l1>l_max && !is_2D))
          continue;
        else
          n_nodes++;

      }
    }
  }

  if (is_2D) {
    printf_log_if (pbi->bispectra_verbose, 2,
      "     * the l1=%d slice has %ld node%s\n",
      pbi->l[index_L1], n_nodes, ((n_nodes==1)?"":"s"));
  }
  else {
    printf_log_if (pbi->bispectra_verbose, 2,
      "     * mesh has %ld node%s\n",
      n_nodes, ((n_nodes==1)?"":"s"));
  }


  // ====================================================================================
  // =                             Rearrange the bispectrum                             =
  // ====================================================================================

  /* Counter for the (l1,l2,l3) nodes */
  int counter_l1l2l3 = 0;

  /* Allocate the array that will contain the rearranged bispectrum */  
  double (*values)[4] = calloc (4*n_nodes, sizeof(double));

  int abort = _FALSE_;

  for(int index_l1 = 0; index_l1 < pbi->l_size; ++index_l1) {

    int l1 = pbi->l[index_l1];

    /* The bidimensional mesh will only contain one value of l1 */
    if (is_2D && index_l1!=index_L1)
      continue;

    /* Counter for the (l2,l3) nodes */
    int counter_l2l3 = 0;

    /* Restrict the interpolation to l1<=l2 if requested */
    int index_l2_min = 0;
    if (pbi->mesh_restrict_l1l2l3)
      index_l2_min = index_l1;

    for (int index_l2=index_l2_min; index_l2 < pbi->l_size; ++index_l2) {

      int l2 = pbi->l[index_l2];
      int index_l3_min = pbi->index_l_triangular_min[index_l1][index_l2];
      int index_l3_max = pbi->index_l_triangular_max[index_l1][index_l2];

      /* Restrict the interpolation to l2<=l3 if requested */
      if (pbi->mesh_restrict_l1l2l3)
        index_l3_min = MAX (index_l2, index_l3_min);

      for (int index_l3=index_l3_min; index_l3<=index_l3_max; ++index_l3) {

        int l3 = pbi->l[index_l3];

        if (l2>l_max || l3>l_max || (l1>l_max && !is_2D))
          continue;

        /* By default, assume that no window function is needed */
        double inverse_window = 1;

        /* Compute the value of the window function for this (l1,l2,l3) */
        if (window_function != NULL) {
          class_call_parallel ((*window_function) (
                                ppr, psp, ple, pbi,
                                l1, l2, l3,
                                X, Y, Z,
                                &inverse_window),
            pbi->error_message,
            pbi->error_message);
        }

        /* Compute the bispectrum in (l1,l2,l3) */
        double bispectrum;

        class_call_parallel (bispectra_at_node (
                               pbi,
                               index_bt,
                               index_l1, index_l2, index_l3,
                               X, Y, Z,
                               &bispectrum,
                               NULL),
          pbi->error_message,
          pbi->error_message);

        /* The first element is the value of the function to be interpolated at the
        nodes, the other arguments are the coordinates of the nodes */

        if (is_2D) {          
          values[counter_l2l3][0] = bispectrum/inverse_window;
          values[counter_l2l3][1] = (double)(l2);
          values[counter_l2l3][2] = (double)(l3);
          values[counter_l2l3][3] = 0;
        }
        else {
          values[counter_l1l2l3][0] = bispectrum/inverse_window;
          values[counter_l1l2l3][1] = (double)(l1);
          values[counter_l1l2l3][2] = (double)(l2);
          values[counter_l1l2l3][3] = (double)(l3);
        }
 
        counter_l2l3++;
        counter_l1l2l3++;
 
      } // for(index_l3)
    } // for(index_l2)
  } // for(index_l1)

  if (abort)
    return _FAILURE_;
  


  // ====================================================================================
  // =                                 Recycle the grid                                 =
  // ====================================================================================

  /* All bispectra share the same (l1,l2,l3) nodes. We use this fact to recycle the
  interpolation grid, ie. the binning of the nodes, across the different bispectra
  types. */

  if ((grid == NULL) || (id == NULL)) {

    int index_l1 = is_2D ? index_L1:0;

    for (int index_bt=0; index_bt < pbi->bt_size; ++index_bt) {

      if (!pbi->interpolate_me[index_bt])
        continue;

      for (int X=0; X < pbi->bf_size; ++X) {
        for (int Y=0; Y < pbi->bf_size; ++Y) {
          for (int Z=0; Z < pbi->bf_size; ++Z) {

            for (int index_mesh=0; index_mesh < 2; ++index_mesh) {

              struct interpolation_mesh * m = pbi->meshes[index_l1][index_bt][X][Y][Z][index_mesh];

              if (m->ready &&
                  m->n_nodes == n_nodes &&
                  (int)(m->max+_EPS_) == l_max &&
                  m->link_length == link_length) {

                grid = m->grid;
                id = m->id;
                goto create_mesh;
                
              }
            }
          }
        }
      }
    }
  }



  // ====================================================================================
  // =                                Generate the mesh                                 =
  // ====================================================================================

  create_mesh:

  /* Generate the mesh calling the mesh_init() function; the function is already
  parallelised. We slightly increase l_max to make sure that the last point in
  the l_sampling (eg. lmax=2000) is included in the mesh. */
  class_call (mesh_init (
                n_dim,
                n_nodes,
                values,
                l_max * (1+ppr->smallest_allowed_variation),
                link_length,
                group_length,
                soft_coeff,
                grid,
                id,
                mesh),
    mesh->error_message,
    pbi->error_message);    

  printf_log_if (pbi->bispectra_verbose, 2,
    "      \\ allocated (grid,mesh)=(%d,%d)=(%g,%g) MBs\n",
    mesh->n_allocated_in_grid, mesh->n_allocated_in_mesh,
    mesh->n_allocated_in_grid*sizeof(double)/1e6,
    mesh->n_allocated_in_mesh*sizeof(double)/1e6);

  /* We do not need the node values anymore */
  free (values);

  return _SUCCESS_;

}




/**
 * Initialise interpolation meshes for each probe (TTT, EEE, EET...) of
 * the bispectrum corresponding to index_bt.
 *
 * This function will initialise two meshes for each probe of the considered
 * bispectrum (probe=TTT, EEE, EET...). The first mesh is the so called fine
 * mesh and will be used to interpolate the configurations where (l1,l2,l3)
 * are all smaller than pbi->mesh_l_turnover. The second mesh is the coarse mesh
 * and will be used to interpolate the other configurations.
 *
 * If the l1 index provided is negative, the function will create meshes for
 * the 3D interpolation of the bispectrum; if it is positive, the function
 * will create 2D meshes for the interpolation of the (l2,l3) plane with
 * fixed l1 = pbi->l[index_L1].
 *
 * This function is basically a wrapper to bispectra_mesh_create(); refer
 * to the documentation there for more details.
 */

int bispectra_mesh_create_for_all_probes (
      struct precision * ppr,
      struct spectra * psp,
      struct lensing * ple,
      struct bispectra * pbi,
      int index_bt, /**< Input: bispectrum type for which the meshes will be computed */
      int index_L1, /**< Input: multipole for which the 2D meshes will be computed; if negative, fill 3D meshes instead */
      struct interpolation_mesh ****** meshes
      )
{
  
  /* Determine whether the user wants to compute 2D or 3D meshes */
  short is_2D = (index_L1 >= 0);
  short is_3D = !is_2D;

  /* If we are computing a 3D mesh, it might take some time */
  if (is_3D)  
    printf_log_if (pbi->bispectra_verbose, 1,
      " -> preparing the interpolation mesh for the %s bispectrum\n",
      pbi->bt_labels[index_bt]);

  for (int X = 0; X < pbi->bf_size; ++X) {
    for (int Y = 0; Y < pbi->bf_size; ++Y) {
      for (int Z = 0; Z < pbi->bf_size; ++Z) {

        if (is_3D)
          printf_log_if (pbi->bispectra_verbose, 0,
            "     * computing fine mesh for bispectrum %s_%s\n",
            pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

        /* Initialise the fine grid; note that we set l_max=l_turnover */
        class_call (bispectra_mesh_create (
                      ppr, psp, ple, pbi,
                      index_bt,
                      X, Y, Z,
                      pbi->window_function[index_bt],
                      index_L1,
                      pbi->mesh_l_turnover,
                      pbi->mesh_link_lengths[0],
                      pbi->mesh_group_lengths[0],
                      pbi->mesh_soft_coeffs[0],
                      NULL,
                      NULL,
                      meshes[index_bt][X][Y][Z][0]),
          pbi->error_message,
          pbi->error_message);

        if (is_3D)
          printf_log_if (pbi->bispectra_verbose, 0,
            "     * computing coarse mesh for bispectrum %s_%s\n",
            pbi->bt_labels[index_bt], pbi->bfff_labels[X][Y][Z]);

        /* Initialise the coarse grid */
        class_call (bispectra_mesh_create (
                      ppr, psp, ple, pbi,
                      index_bt,
                      X, Y, Z,
                      pbi->window_function[index_bt],
                      index_L1,
                      pbi->l_max,
                      pbi->mesh_link_lengths[1],
                      pbi->mesh_group_lengths[1],
                      pbi->mesh_soft_coeffs[1],
                      NULL,
                      NULL,
                      meshes[index_bt][X][Y][Z][1]),
          pbi->error_message,
          pbi->error_message);

      } // for(Z)
    } // for(Y)
  } // for(X)

  return _SUCCESS_;

}


