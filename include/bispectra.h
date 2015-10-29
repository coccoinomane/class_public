/** @file spectra.h Documented includes for bispectra module */

#ifndef __BISPECTRA__
#define __BISPECTRA__

#include "bessel.h"
#include "lensing.h"
#include "spectra.h"
#include "slatec_3j_C.h"
#include "binary.h"
#include "mesh_interpolation.h"


/**
 * Categories of bispectra that SONG can compute.
 *
 * For details, see documentation for the bispectra.c file.
 */
enum bispectra_types {
  separable_bispectrum,     /**< bispectra derived from a separable shape function S(k1,k2,k3) (e.g. local, equilateral, orthogonal)*/
  non_separable_bispectrum, /**< bispectra derived from an arbitrary shape function S(k1,k2,k3) (e.g. galileon) */
  analytical_bispectrum,    /**< bispectra described by a closed form (analytical limits, cmb-lensing bispectrum, etc) */
  intrinsic_bispectrum      /**< bispectra that require second-order transfer functions */
};



/**
 * Possible interpolation methods for the bispectrum.
 */
enum bispectra_interpolation_method {
  smart_interpolation,
  bilinear_interpolation,
  trilinear_interpolation,
  mesh_interpolation_2D,
  mesh_interpolation_3D
};


struct bispectra;


/**
 * Type for the bispectrum window functions
 */
typedef int window_function_type(
  struct precision * ppr,
  struct spectra * psp,
  struct lensing * ple,
  struct bispectra * pbi,
  int l1, int l2, int l3,
  int X, int Y, int Z,
  double *result
  );


/**
 * Type for the analytical bispectra functions
 *
 * The three threej variables are the following 3j-symbols, respectively:
 *
 *   ( l1, l2, l3 )   ( l1, l2, l3 )   ( l1, l2, l3 )
 *   (  2,  0, -2 )   ( -2,  2,  0 )   (  0, -2,  2 )
 *
 * They are needed only for the quadratic and CMB-lensing bispectra.
 */
typedef int analytical_function_type(
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
  double *result
  );



struct bispectra {

  // ====================================================================================
  // =                                     Flags                                        =
  // ====================================================================================

  short has_bispectra; /**< Should we compute any bispectra at all? If this flag is set
                       to _FALSE_, the bispectrum module will be skipped altogether. */

  /**
   * Should we compute the lensing of the bispectrum?
   *
   * Gravitational lensing due to intervening matter bends the trajectory of
   * photons from the last scattering surface to us. This effect alters the
   * shape of the CMB bispectra and is called lensing of the bispectrum.
   *
   * We parameterise the lensing of the bispectrum as an additive factor to
   * the unlensed bispectrum.
   *
   * If this flag is turned on, the lensing correction will be stored in
   * pbi->bispectra_lensing_correction and added to the unlensed bispectrum
   * in pbi->bispectra. Otherwise, pbi->bispectra will contain the unlensed
   * bispectrum.
   *
   * For some analytical bispectra (squeezed approximations, CMB-lensing...)
   * the lensing is computed by just substituting the unlensed C_l in their
   * definition with the lensed C_l.
   *
   * For arbitrary bispectra, however, the lensing correction requires
   * convolving the unlensed bispectrum with the lensing potential C_l^PP.
   * We do so in bispectra_lensing().
   */          
  short has_lensed_bispectra;

  short need_3j_symbols; /**< Are the Wigner 3j-symbols needed to compute the requested bispectra?  
                         This is the case for the CMB-lensing bispectrum and for the quadratic 
                         bispectrum. */



  // ====================================================================================
  // =                             Indices of bispectra                                 =
  // ====================================================================================

  /* What bispectra types should be compute? */
  short has_local_model;              /* Local model */
  short has_equilateral_model;        /* Equilateral model */  
  short has_orthogonal_model;         /* Orthogonal model */
  short has_galileon_model;           /* Galileon inflation (arXiv:1108.0305) */
  short has_intrinsic;                /* Bispectrum induced by non-linear dynamics */
  short has_intrinsic_squeezed;       /* Squeezed limit for the 2nd-order bispectrum (Creminelli et al. 2012, Bartolo et al. 2012, Lewis 2012) */
  short has_local_squeezed;           /* Squeezed limit for local model (Guangui et al. 1994) */
  short has_cosine_shape;             /* The above, multiplied by an oscillating function in l */
  short has_cmb_lensing;              /* CMB-lensing bispectrum (Eq. 4.5 of arxiv:1101.2234) */
  short has_cmb_lensing_squeezed;     /* Squeezed limit of the CMB-lensing bispectrum (Eq. 5.20 of arxiv:1101.2234) */
  short has_cmb_lensing_kernel;       /* Squeezed limit of the CMB-lensing bispectrum (kernel only, Eq. 5.20 of arxiv:1101.2234) */
  short has_quadratic_correction;     /* Four-point contribution to the bispectrum */
  short has_test_bispectrum;          /* Simple C_l x C_l bispectrum suitable to test lensing */

  /* Indices for the above bispectra types */
  int index_bt_local;                 /* Index for the bispectrum for a local model */
  int index_bt_equilateral;           /* Index for the bispectrum for a equilateral model */
  int index_bt_orthogonal;            /* Index for the bispectrum for a orthogonal model */
  int index_bt_galileon_gradient;     /* Index for the bispectrum for the pi_dot*pi_grad^2 term in Galileon inflation */
  int index_bt_galileon_time;         /* Index for the bispectrum for the pi_dot^3 term in Galileon inflation */
  int index_bt_intrinsic;             /* Index for the bispectrum induced by nonlinear dynamics */
  int index_bt_intrinsic_squeezed;    /* Index for the intrinsic bispectrum in the squeezed limit */  
  int index_bt_intrinsic_squeezed_unlensed; /* Index for the unlensed intrinsic bispectrum in the squeezed limit */  
  int index_bt_local_squeezed;        /* Index for the local-model bispectrum in the squeezed limit */  
  int index_bt_cosine;                /* Index for the oscillating bispectrum */  
  int index_bt_cmb_lensing;           /* Index for the bispectrum of CMB-lensing */
  int index_bt_cmb_lensing_squeezed;  /* Index for the bispectrum of CMB-lensing in the squeezed limit */
  int index_bt_cmb_lensing_kernel;    /* Index for the bispectrum of CMB-lensing in the squeezed limit (kernel only) */
  int index_bt_quadratic;             /* Index for the bispectrum induced by a quadratic correction to the distribution function */
  int index_bt_test;                  /* Index for the test bispectrum */
  int bt_size;                        /* Total number of bispectra types requested */

  /* Array of strings that contain the text labels of the various bispectra */
  char bt_labels[_MAX_NUM_BISPECTRA_][_MAX_LENGTH_LABEL_];

  /* Array that identifies whether a bispectrum is separable, non-separable, analytic or intrinsic. */
  enum bispectra_types bispectrum_type[_MAX_NUM_BISPECTRA_];
  
  /* How many bispectra of each type shall we compute? Indexed as pbi->n[bispectrum type], where
  the bispectrum type can be separable, non_separable, analytical, intrinsic */
  int n[_MAX_NUM_BISPECTRA_];
  
  /* Given a generic bispectrum B_l1l2l3 with even parity (such as those involving only
  temperature and E-mode polarisation), its reduced form b_l1l2l3 is defined by
    B_l1l2l3 = sqrt((2*(l1+1)+2*(l2+1)+2*(l3+1))/(4pi)) * (l1,l2,l3)(0,0,0) * b_l1l2l3,
  where the 3j-symbol (l1,l2,l3)(0,0,0) enforces the parity invariance of B_l1l2l3.
  This implies that the bispectrum B_l1l2l3 vanishes for configurations where
  l1+l2+l3 is odd, while the reduced bispectrum is non-zero for all (l1,l2,l3).
  In the Fisher matrix estimator, only B_l1l2l3 enters, therefore the odd configurations
  do not matter. However, computing b_l1l2l3 would make the interpolation more stable as the
  odd configurations fill the gaps in the (l1,l2,l3) space. Therefore, we choose to store
  the reduced bispectrum b_l1l2l3 rather than B_l1l2l3.
  
  However, for certain bispectra it is not possible to analytically extract the
  3j-symbol (l1,l2,l3)(0,0,0). Examples of such bispectra are the |m|>0 intrinsic
  bispectra (see chapter 6 of my thesis) and the CMB-lensing/quadratic bispectra
  for polarisation (see Sec. 4 of http://arxiv.org/abs/1101.2234). For these
  bispectra, the reduced bispectrum b_l1l2l3 has to be computed from B_l1l2l3
  by numerically dividing the latter by (l1,l2,l3)(0,0,0). Therefore, it not
  possible to obtain the value of b_l1l2l3 for odd l1+l2+3. We take note of these
  non-reducible bispectra using the following array. */
  int has_reduced_bispectrum[_MAX_NUM_BISPECTRA_];

  /**
   * Logical array to mark the bispectra that need to be lensed.
   */
  int lens_me[_MAX_NUM_BISPECTRA_];

  /**
   * Logical array to mark the bispectra that need to be interpolated.
   */
  int interpolate_me[_MAX_NUM_BISPECTRA_];
  
  /**
   * Logical array to mark the bispectra that are symmetric with respect to
   * permutations of (l1,X1), (l2,X2) and (l3,X3).
   *
   * By construction, all bispectra must be symmetric with respect to any
   * simultaneous permutation of a multipole index (l=2,3,4...) and the
   * associated field index (X=T,E,B...). This follows from the very
   * definition of a bispectrum: b = <a_l1m1^X1 a_l2m2^X2 a_l3m3^X3>.
   *
   * In SONG we rely on this symmetry to compute the bispectra: we store
   * them only for the configurations where l1>=l2>=l3 and then obtain
   * other combinations by permuting the field indices (X1, X2, X3).
   *
   * Some bispectra however are defined in such a way that the three
   * multipoles are not on the same ground. For example, in the squeezed
   * limit approximations (eg. for the CMB-lensing, local and intrinsic
   * bispectra) one of the multipoles is assigned the large wavelength
   * mode, ie. the mode with small l.
   * 
   * For these bispectra, the 1<->2<->3 symmetry is lost. SONG computes
   * them only in those configurations that satisfy l1>=l2>=l3 and set
   * their value to zero for the other ones.
   */
  int is_symmetric[_MAX_NUM_BISPECTRA_];
  
  /**
   * Logical array to mark the bispectra that will be lensed with the convolution
   * method in the bispectra_lensing() function.
   *
   * The other bispectra will be lensed analytically on a one-by-one case. For
   * example the CMB-lensing bispectrum and the squeezed approximations will be
   * lensed by sustituting the unlensed C_l in their definition with the lensed
   * ones.
   */
  int lens_me_brute_force[_MAX_NUM_BISPECTRA_];
  
  /**
   * Add the bolometric and redshift corrections to the second-order bispectra?
   *
   * Not used unless you requested the intrinsic bispectrum. Refer to
   * bispectra2_add_quadratic_correction() for documentation.
   */
  short add_quadratic_correction;

  /**
   * Functions to use for the computation of the analytical bispectra.
   *
   * This array of functions associates to a bispectrum the function used to
   * compute it.
   *
   * It is defined only for the analytical bispectra; non-analytical bispectra
   * will have NULL entries.
   */
  analytical_function_type * bispectrum_function[_MAX_NUM_BISPECTRA_];

  /**
   * Window functions to use to interpolate each bispectrum.
   *
   * The interpolation in (l1,l2,l3) of a given bispectrum might be easier if it is
   * divided by a window function beforehand. For example, bispectra peaked on
   * squeezed configurations might benefit from a window function that peaks itself
   * on squeezed scales, in order to make the interpolated function smoother.
   *
   * This array contains the window functions associated to each bispectrum. A NULL
   * entry implies that the corresponding bispectrum does not have a window function.
   */
  window_function_type * window_function[_MAX_NUM_BISPECTRA_];



  // ====================================================================================
  // =                                Indices of fields                                 =
  // ====================================================================================

  /* The fields considered in this module and in the Fisher one are denoted by indices 'bf'.
  These can be temperature (T), E-polarisation (E), B-polarisation (B), Rayleigh (R)...  */
  short has_bispectra_t;
  short has_bispectra_e;
  short has_bispectra_b;

  int index_bf_t;
  int index_bf_e;
  int index_bf_b;
  int bf_size;

  /* Actual number of bispectra to be computed (usually 'bf_size' to the power of 3, it is 8 for
  T and E: TTT, TTE, TET, ETT, EET, ETE, TEE, EEE) */
  int n_probes;

  /* Parity of the considered field. Even parity (T,E) is represented by zero, odd parity
  by 1. Indexed as field_parity[index_bf] */
  int field_parity[_MAX_NUM_FIELDS_];
  
  /* Spin of the considered field; this is equal to 2 for polarisation and 0 for
  intensity (see, e.g., Hu & White 1997) */
  int field_spin[_MAX_NUM_FIELDS_];

  /* Array that relates the bispectrum field indices (T,E,B) to the index of the transfer functions */
  int index_tt_of_bf[_MAX_NUM_FIELDS_];
  int index_tt2_of_bf[_MAX_NUM_FIELDS_];

  /* Array that relates each pair of fields (TT,EE,TE...) to their C_l power spectrum, as computed
  in the spectra.c module. This is used to build the cross-power spectrum in the Fisher module.
  For example,
  index_ct_of_bf_bf[pbi->index_bf_t][pbi->index_bf_t] = psp->index_ct_tt
  index_ct_of_bf_bf[pbi->index_bf_t][pbi->index_bf_e] = psp->index_ct_te  */
  int index_ct_of_bf_bf[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_];

  /* Associate to each field X=T,E,... its C_l correlation <X phi> with the lensing potential and
  its C_l correlation <X zeta> with the curvature perturbation. */
  int index_ct_of_phi_bf[_MAX_NUM_FIELDS_];
  int index_ct_of_zeta_bf[_MAX_NUM_FIELDS_];
  int index_ct_of_t_bf[_MAX_NUM_FIELDS_];
  
  /* Same but for the lensed C_l's */
  int index_lt_of_bf_bf[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_];

  /* Array of strings that contain the text labels of the various fields */
  char bf_labels[_MAX_NUM_FIELDS_][_MAX_LENGTH_LABEL_]; /* T,E... */
  char bfff_labels[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_][_MAX_LENGTH_LABEL_]; /* TTT,EEE,TTE... */


  // ====================================================================================
  // =                                    Sampling in l                                 =
  // ====================================================================================

  int * l;                /**< 1D list of multipole values, used as a grid to determine
                          the (l1,l2,l3) configurations where the bispectrum is computed */
  int l_size;             /**< Size of l */
  int l_min;              /**< Minimum value in pbi->l (=2 in most cases) */
  int l_max;              /**< Maximum value in pbi->l */
  int full_l_size;        /**< l_max - l_min + 1 */
  
  /* The bispectrum depends on 3 l-values (l1,l2,l3), which must satisfy the triangular inequality:
  |l1-l2| <= l3 <= l1+l2.  We determine the range of l3 within pbi->l using the following two arrays. */
  int ** l_triangular_size;       /* l_triangular_size[index_l1][index_l2] is the number of allowed l3 values in pbi->l given l1 and l2 */
  int ** index_l_triangular_min;  /* index_l_triangular_min[index_l1][index_l2] is the index closest to abs(l1-l2) in pbi->l */
  int ** index_l_triangular_max;  /* index_l_triangular_max[index_l1][index_l2] is the index closest to l1+l2 in pbi->l */

  int ** l3_size; /* l3_size[index_l1][index_l2-index_l1] is the number of l3 multipoles considered for bispectra computation. 
                  The indices refer to the l values contained in pbi->l. With respect to l_triangular_size, it has index_l1<=index_l2. */
  int *** l3; /* l3[index_l1][index_l2] is the array of l3 multipoles considered for bispectra computation, of
              size l3_size[index_l1][index_l2-index_l1]. The indices refer to the l values contained in pbi->l. 
              With respect to l_triangular_size, it has index_l1<=index_l2. */

  /* Array of indices associated to all (l1,l2,l3) configurations, of size 
  n_independent_configurations. Indices are associated to the (l1,l2,l3) configurations
  belonging to pbi->l where l3 satisfies the triangular condition,and l1>=l2>=l3. The
  array should be indexed as:

      pbi->index_l1_l2_l3[index_l1]
                         [index_l1 - index_l2]
                         [index_l3_max - index_l3]
                         
  where index_l1, index_l2, index_l3 are the indices for l1,l2,l3 in the pbi->l array;
  index_l3 has to be smaller than index_l3_max, which is given by:
   
    index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]) */
  long int *** index_l1_l2_l3;
  
  /* Number of (l1,l2,l3) configurations that will be stored for each bispectrum, and
  size of index_l1_l2_l3. */
  long int n_independent_configurations;

  /**
   * Array where the bispectra will be stored.
   *
   * The bispectra array is indexed in the following way:
   * 
   * bispectra [index_bt][X][Y][Z][index_l1_l2_l3]
   *
   * - index_bt is the bispectrum type index (local, equilateral, intrinsic...) and 
   *   it can be any of the active pbi->index_bt_XXX indices. 
   * - X, Y and Z are field indices (temperature, E-polarisation, B-polarisation...)
   *   and can be any of the active pbi->index_bf_XXX indices.
   * - index_l1_l2_l3 is the index pointing to the (l1,l2,l3) configurations; see
   *   documentation for index_l1_l2_l3 for more information.
   *
   * If the flag has_lensed_bispectra is turned on, the bispectra stored in this array
   * will be lensed.
   */
  double ***** bispectra;

  
  /**
   * Array where the correction to the bispectrum due to lensing will be stored.
   *
   * Gravitational lensing due to intervening matter bends the trajectory of
   * photons from the last scattering surface to us. This effect alters the
   * shape of the CMB bispectra and is called lensing of the bispectrum.
   *
   * We parameterise the lensing of the bispectrum as an additive factor to
   * the unlensed bispectrum.
   */
  double ***** lensing_correction;
  
  /* In computing the CMB-lensing bispectrum, turn the lensing potential C_l to zero on small scales,
  where we cannot trust them. This won't change the result because these C_l's are very small for
  large l's. In CAMB, Antony sets C_l^TP=0 for l>300 and C_l^EP=0 for l>40. */
  int lmax_lensing_corrT;
  int lmax_lensing_corrE;


  // ====================================================================================
  // =                                Filter functions                                  =
  // ====================================================================================

  /* In the local model, the primordial bispectrum is simply given by
        
    B = 2 * f_NL^local * (P(k1)*P(k2) + P(k1)*P(k3) + P(k2)*P(k3))
      
  which means that two of the separated integrals are equal. Here we follow the notation by
  Komatsu and define two separable integrals, one with the primordial spectrum (beta) and one
  without it (alpha). We use the notation of Komatsu, Spergel & Wandelt 2005; note that in
  Komatsu et al 2001 alpha -> b_NL and beta -> b_L.
  
  Each fiter function has an 'index_bf' level that goes over T,E,B...
  */
  double *** alpha;              /* alpha[index_bf][index_l][index_r] for temperature */
  double *** beta;               /* beta[index_bf][index_l][index_r] for temperature */


  /* In the equilateral model (Creminelli et al. 2006), the primordial bispectrum is given by:
      
    B = 6 * ( - P(k1)*P(k2) - P(k1)*P(k3) - P(k2)*P(k3)
              - 2 (P(k1)*P(k2)*P(k3))^(2/3)   + 5 symm.
              + (P(k1)*P(k2)^2*P(k3)^3)^(1/3) + 5 symm. ).
           
  There are two extra separable integrals with respect to the local case, which we name 'gamma'
  and 'delta' following Creminelli et al. 2006 (see eqs. 15-18).
     
  The same integrals are needed for the orthogonal shape (Senatore et al. 2010, eq. 3.2):

     B = 6 * ( - 3*P(k1)*P(k2) - 3*P(k1)*P(k3) - 3*P(k2)*P(k3)
               - 8*(P(k1)*P(k2)*P(k3))^(2/3)   + 5 symm.
               + 3*(P(k1)*P(k2)^2*P(k3)^3)^(1/3) + 5 symm. ).
   */
  double *** gamma;              /* gamma[index_l][index_r] for temperature */
  double *** delta;              /* delta[index_l][index_r] for temperature */



  // ====================================================================================
  // =                                 Power spectra                                    =
  // ====================================================================================

  /* First-order primordial power spectrum as a function of k. It is indexed as pbi->pk[index_k]
  where 'index_k' indexes ptr->k. Its size is ptr->k_size[ppt->index_md_scalars]. */
  double * pk;

  /* Same as above, but now 'index_k' indexes ppt->k. Its size is ppt->k_size[ppt->index_md_scalars]. */
  double * pk_pt;
  
  /* Array that contains the first-order Cl's. It is indexed as cls[index_ct][l-2], where index_ct
  can be psp->index_ct_tt, psp->index_ct_te, etc. 
  IMPORTANT: these arrays are computed for all l's between 2 and pbi->l_max, hence they have to be
  indexed as cls[index_ct][l-2] */
  double ** cls;
  double ** d_lsq_cls;        /* Logarithmic derivative of the C_l's, as in Lewis 2012 */
  double ** lensed_cls;
  double ** lensed_d_lsq_cls;

  /* Array that contains the measure for the trapezoidal integration over k of the first-order transfer
  functions. It is defined as ptr->k[i+1] - ptr->k[i-1]. */
  double * delta_k;



  // ====================================================================================
  // =                            Bispectrum interpolation                              =
  // ====================================================================================
  
  /**
   * What method should we use for interpolating the bispectra (bilinear,
   * mesh_2D, mesh_3D...)
   */
  enum bispectra_interpolation_method interpolation_method;


  /**
   * Parameters for the mesh interpolation of the bispectra.
   *
   * The main difficulty in interpolating the bispectrum is that it is not defined on
   * a cubic grid. In fact, the triangular condition |l2-l3| <= l1 <= l2+l3, results
   * in a domain for (l1,l2,l3) that has the shape of a "tetrapyd", the union of two
   * triangular pyramids through the base (see Fig. 2 of Fergusson et al. 2012).
   *
   * The problem can be circumvented by deforming the allowed region to a cube via a
   * geometrical transformation and then using trilinear interpolation (Fergusson et al.
   * 2009). While viable, this approach would force us to discard the points that do
   * not fall in the transformed grid, thus requiring a finer l-sampling.
   *
   * Rather than relying on a cubic grid, we devise a general interpolation technique
   * that is valid on any mesh. We first define a correlation length L and divide the
   * tetrapyd domain in boxes of side L. To compute the interpolation in an arbitrary
   * triplet (l1,l2,l3), we consider the values of all the nodes in the box where
   * the triplet lives and in the adjacent ones. To each node, we assign a weight that
   * is inversely proportional to its distance from (l1,l2,l3).
   *
   * The problem with this approach is that, the mesh being inhomogeneous, there might
   * be a group of close nodes in one direction that influences the interpolated value
   * in (l1,l2,l3) much more than a closer point in the opposite direction. In order
   * to prevent this, we weight down the nodes that have a high local density within a
   * certain distance from them.
   *
   * This mesh interpolation technique relies on two free parameters:
   *
   * - The linking length L which sets the size of the local region influencing
   *   the interpolation. In a homogenous grid, it should correspond roughly to the
   *   largest distance of two neighbouring points.
   *
   * - The grouping length, that is the distance below which many close nodes are
   *   considered as a single one. It is used to avoid the interpolation being
   *   determined by a bunch of close nodes in one direction. In a homogenous grid,
   *   The grouping length should roughly correspond to the shortest distance between
   *   two points.
   *
   * Since the bispectra sampling is denser for small l and sparser for high l, we use
   * two meshes: a fine mesh with a small linking length, to interpolate the bispectrum
   * in the region where all l are small, and a coarse mesh with a larger linking length,
   * to interpolate it in the remaining region. The turnover point between the two
   * regions is determined based on the l sampling in pbi->l, and is memorised in
   * pbi->mesh_l_turnover.
   */
  //@{

  /** 
   * Data needed to interpolate the bispectra with the mesh technique.
   *
   * The meshes array contains two interpolation meshes per bispectrum:
   *
   * - A fine mesh to interpolate the bispectrum for l < l_turnover.
   * - A coarse mesh to interpolate the bispectrum for l >= l_turnover.
   *
   * The fine and coarse mesh have a linking length roughly equal to the
   * distance between the two closest and furthest neighbouring points,
   * respectively.
   * 
   * Why two meshes? Using only a fine mesh all the way to l_max would
   * imply more memory usage, because a small linking length means more
   * interpolation boxes to keep track of. On the other hand, the coarse
   * mesh is quick to compute but takes much more time to use for
   * interpolation, because larger interpolation boxes mean more nodes
   * to include and the interpolation algorithm goes as n_nodes^2 or
   * n_nodes^3.
   *
   * The meshes array is indexed in the same way as pbi->bispectra, but
   * with two extra levels:
   *
   *  pbi->meshes[index_l1][index_bt][X][Y][Z][index_mesh]
   *
   * For 3D interpolation, the index_l1 level is collapsed to a single
   * element, because a single 3D mesh is enough to interpolate the bispectrum
   * for all values of l1. For 2D interpolation, it has size pbi->l_size.
   *
   * The level indexed by index_mesh has only two elements, corresponding
   * to the fine mesh (index_mesh=0) and to the coarse mesh (index_mesh=1).
   */
  struct interpolation_mesh ******* meshes;
  
  double mesh_link_lengths[2];  /**< Linking lengths of the fine & coarse interpolation meshes */
  double mesh_group_lengths[2]; /**< Grouping lengths of the fine & coarse interpolation meshes */
  double mesh_soft_coeffs[2];   /**< Softening of the linking length of the fine & coarse interpolation meshes, fixed to 0.5 */
  int mesh_l_turnover;  /**< Turnover multipole between the fine and coarse interpolation meshes;
                        determined based on the l-sampling in pbi->l */

  /**
   * Should we restrict the (l1,l2,l3) domain for the mesh interpolation?
   *
   * When building the interpolation mesh for the bispectrum, we can choose to
   * restrict the domain to those configurations with l1<=l2<=l3, which results
   * in less memory usage. This is fine to build the Fisher matrix, which is
   * obtained via a sum over l1<=l2<=l3, but it is not ok to compute the lensing
   * of the bispectrum, where the sum is unrestricted.
   */
  short mesh_restrict_l1l2l3;

  //@}



  // ====================================================================================
  // =                                  Disk storage                                    =
  // ====================================================================================


  char bispectra_dir[_FILENAMESIZE_];  /**< Directory containing the bispectra. If it already exists, and
                                       ppr->load_bispectra_from_disk==_TRUE_, the bispectra will be read from this folder 
                                       into the array pbi->bispectra. If it does not exist, and ppr->store_bispectra_to_disk==_TRUE_,
                                       the bispectra will be first computed and then written to this folder from the array
                                       pbi->bispectra. Either way, the directory contains one binary file for each bispectrum type,
                                       for a total of pbi->bt_size files. The file corresponding to index_bt is located at
                                       pbi->bispectra_paths[index_bt]; its stream is in pbi->bispectra_files[index_bt]. */


  char ** bispectra_paths; /**< bispectra_paths[index_bt] is the path to the file containing the bispectrum type
                           corresponding to index_bt. Used only if ppr->store_bispectra_to_disk==_TRUE_ or
                           ppr->load_bispectra_from_disk==_TRUE_. */

  FILE ** bispectra_files; /**< bispectra_paths[index_bt] is the pointer to the file containing the bispectrum type
                           corresponding to index_bt. Used only if ppr->store_bispectra_to_disk==_TRUE_ or
                           ppr->load_bispectra_from_disk==_TRUE_. */

  FILE * bispectra_status_file;                      /**< NOT IMPLEMENTED YET */
  char bispectra_status_path[_FILENAMESIZE_];        /**< NOT IMPLEMENTED YET */


  // ===========================================================================================
  // =                                   Technical parameters                                  =
  // ===========================================================================================
  
  short bispectra_verbose;                   /**< Flag regulating the amount of information sent to standard output (none if set to zero) */                                                  
  long int n_total_configurations;           /**< Number of (l1,l2,l3) configurations compatible with the triangular condition */
  ErrorMsg error_message;                    /**< Zone for writing error messages */
  long int count_allocated_for_bispectra;    /**< Number of elements allocated in pbi->bispectra */
  long int count_memorised_for_bispectra;    /**< Number of elements memorised in pbi->bispectra */
  short output_binary_bispectra;             /**< Should we write the all bispectra and accessory infromation to fiel? */
  short always_interpolate_bispectra;        /**< Should we interpolate also the analytical bispectra? */

};


/**
 * Workspace that contains the intermediate results for the integration of bispectra which 
 * have separable shape functions.
 *
 * The technique we emply consists in integrating separately the k3,k2,k1 and r levels
 * by using the trapezoidal rule and exploiting the fact that the shape function is separable:
 *
 * B(k1,k2,k3) = B1(k1)*B2(k2)*B3(k3) + symm
 *
 * For a brief review of these shapese, look at Komatsu et al. WMAP7 paper, eqs. 60, 63, 64.
 *
 * In the case of a local model, where we have B = 2 * f_NL * (P(k1)*P(k2) + P(k1)*P(k3) + P(k2)*P(k3)),
 * the separation is of the type:
 *
 * B1(k) = B2(k) = (2*f_NL)^1/3 * P(k)
 * B3(k) = (2*f_NL)^1/3
 *
 *
 */
struct bispectra_workspace_separable {
  
  
  /* Grid in the integration variable 'r'.  This is the parameter that stems from the expansion into
    spherical Bessels of the Dirac Delta \delta(\vec{k1}+\vec{k2}+\vec{k3}) */
  double r_min;
  double r_max;
  int r_size;
  double * r;

  /* Measure for the trapezoidal rule: pwb->r[i+1] - pwb->r[i-1] */
  double * delta_r;

  /* Arrays that will contain the integrand functions for the separable. For the simple local model, one needs
    only two of them (alpha, beta), while for the more complicated models (e.g. equilateral and orthogonal),
    one needs more. */
  double ** alpha_integrand;    /* alpha_integrand[thread][index_r] */
  double ** beta_integrand;     /* beta_integrand[thread][index_r] */
  double ** gamma_integrand;    /* gamma_integrand[thread][index_r] */
  double ** delta_integrand;    /* delta_integrand[thread][index_r] */


};



/**
 * Workspace that contains the intermediate results for the integration of the non-separable
 * bispectra.
 *
 * The technique consists in integrating sequentially the k3,k2,k1 and r levels
 * by using the trapezoidal rule.
 *
 */
struct bispectra_workspace_non_separable {

  /* Pointer to the desired shape function */
  int (*shape_function) (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3, 
    double pk_1, double pk_2, double pk_3, 
    double * out
    );

  /* Considered bispectrum XYZ (e.g. TTT, TEE, ...) */
  int X;
  int Y;
  int Z;
  
  /* The shape function does not need to be sampled as finely as the transfer functions, as it generally does not
  oscillate. We shall use interpolation in order to obtain the value in the integration grid. */
  double * k_smooth_grid;
  int k_smooth_size;

  /* Temporary pointer to store arrays that have pwb->k_smooth_size elements (one per thread) */
  double ** f;

  /* Window function for the interpolation of the k-space bispectrum. It is indexed as pwb->k_window[index_k]
  where 'index_k' indexes pwb->k_smooth_grid. Its size is pwb->k_smooth_size. */
  double * k_window;

  /* Inverse window function for the interpolation of the k-space bispectrum. It is indexed as pwb->k_window[index_k]
  where 'index_k' indexes ptr->k. Its size is ptr->k_size[ppt->index_md_scalars]. */
  double * k_window_inverse;


  /* Grid in the integration variable 'r'.  This is the parameter that stems from the Rayleigh expansion 
  of the Dirac Delta \delta(\vec{k1}+\vec{k2}+\vec{k3}) */
  double r_min;
  double r_max;
  int r_size;
  double * r;


  /* Array to contain the integral over k3:
  
                             /
     INT_l3 (r,k1,k2)   =   |  dk3  k3^2 * j_l3(r*k3) * T_l3(k3,k1,k2)
                            /

    The array is indexed as pbi->integral_over_k3[index_l3][index_r][index_k1][index_k2]  */
  double **** integral_over_k3;
                           
  /* Integration grid in k3 for a given k1 and k3, one for each thread: k3_grid[thread][index_k] */
  double ** k3_grid;

  /* Maximum size of the k3_grid */
  int k3_size_max;
  
  /* Limits in the k3 integration as a function of k1 and k3: k3_lower[index_k1][index_k2-index_k1] */
  int ** index_k3_lower;
  int ** index_k3_upper;
  int ** k3_grid_size;


  /* Array to contain the integral over k2:  
  
                             /
     INT_l2_l3 (r,k1)   =   |  dk2  k2^2 * j_l2(r*k2) * INT_l3 (r,k1,k2) * T_l2(k2)
                            /

    The array is indexed as pbi->integral_over_k2[index_l2][index_l3][index_r][index_k1] */
  double **** integral_over_k2;



  /* Array to contain the integral over k1:
  
                             /
     INT_l1_l2_l3 (r)   =   |  dk1  k1^2 * j_l1(r*k1) * INT_l2_l3 (r,k1) * T_l1(k1)
                            /
  
    The array is indexed as pbi->integral_over_k1[index_l1_l2_l3][index_r]. We use index_l1_l2_l3 because the integral over k1
    yields a function of (l1,l2,l3) that is symmetric under permutations of (l1,l2,l3).  Hence, we compute and store it only for
    l1>=l2>=l3 configurations (that satisfy the triangular condition) */
  double ** integral_over_k1;



  /* Array to contain the integral over r:
  
                          /
     INT_l1_l2_l3    =   |  dr  r^2 * INT_l1_l2_l3(r)
                         /

    The array is indexed as pbi->integral_over_r[index_l1][index_l2][index_l3-index_l_triangular_min], and is basically
    the unsymmetrised bispectrum. */
  double *** integral_over_r;

  

  /* Array that contains the interpolated values of the above integrals in ptr->k. Each thread has one.
  Indexed as integral_splines[thread][index_k] and interpolated_integral[thread][index_k], where
  index_k belongs to ptr->k. */
  double ** integral_splines;
  double ** interpolated_integral;
  
  /* Same as above, but for the k3 integration grid (one per thread) */
  double ** delta_k3;

  /* Array that contains the pwb->r[i+1] - pwb->r[i-1] values needed for the trapezoidal rule by 
    the function 'bispectra_smooth_integration_over_r' */
  double * delta_r;

  /* Memory counters for the various arrays */
  long int count_allocated_for_integral_over_k1;
  long int count_allocated_for_integral_over_k2;
  long int count_allocated_for_integral_over_k3;
  long int count_allocated_for_integral_over_r;

  long int count_memorised_for_integral_over_k1;
  long int count_memorised_for_integral_over_k2;
  long int count_memorised_for_integral_over_k3;
  long int count_memorised_for_integral_over_r;
  
  
};







/*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif


  int bispectra_init(
       struct precision * ppr,
       struct background * pba,
       struct thermo *pth,
       struct perturbs * ppt,
       struct bessels * pbs,
       struct transfers * ptr,
       struct primordial * ppm,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi
       );

  int bispectra_at_node (
      struct bispectra * pbi,
      int index_bt,
      int index_l1, int index_l2, int index_l3,
      int X, int Y, int Z,
      double * bispectrum,
      double * bispectrum_unlensed
      );

  int bispectra_at_l3_linear (
      struct transfers * ptr,
      struct bispectra * pbi,
      int index_bt,
      int index_l1, int index_l2, int l3,
      int X, int Y, int Z,
      int interpolate,
      double * bispectrum,
      double * bispectrum_unlensed
      );
  
  int bispectra_at_l2l3_bilinear (
      struct transfers * ptr,
      struct bispectra * pbi,
      int index_bt,
      int index_l1, int l2, int l3,
      int X, int Y, int Z,
      int extrapolate,
      double * bispectrum,
      double * bispectrum_unlensed
      );
  
  int bispectra_at_l2l3_mesh (
      struct bispectra * pbi,
      int index_bt,
      int index_l1,
      int l1, int l2, int l3,
      int X, int Y, int Z,
      double * bispectrum
      );
  
  int bispectra_free(
       struct precision * ppr,
       struct perturbs * ppt,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi
       );

  int bispectra_free_type_level(
       struct bispectra * pbi,
       int index_bt
       );

  int bispectra_indices(
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
          );

  int bispectra_primordial_power_spectrum (
      struct background * pba,
      struct perturbs * ppt,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi
      );

  int bispectra_cls (
      struct precision * ppr,
      struct perturbs * ppt,
      struct spectra * psp,
      struct lensing * ple,
      struct bispectra * pbi
      );


  int bispectra_harmonic(
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
      );
      
      
  int bispectra_output(
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
      );
      
      
  int bispectra_get_r_grid (
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct bispectra * pbi,
      double * tau_sampling,
      int tau_size,
      double ** r_grid,
      int * r_size,
      double * r_min,
      double * r_max,
      double ** delta_r
      );


  int bispectra_separable_init(
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      struct bispectra_workspace_separable * pwb
      );
  
  
  int bispectra_separable_workspace_init(
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      struct bispectra_workspace_separable * pwb
      );
  
  int bispectra_separable_filter_functions(
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      int index_bf,
      struct bispectra_workspace_separable * pwb
      );


  int bispectra_separable_integrate_over_r(
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      int index_bt,
      int index_bf_1,
      int index_bf_2,
      int index_bf_3,
      struct bispectra_workspace_separable * pwb
      );



  int bispectra_separable_workspace_free(
      struct bispectra * pbi,
      struct bispectra_workspace_separable * pwb
      );


  int bispectra_analytical_init(
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct spectra * psp,
      struct lensing * ple,
      struct bispectra * pbi
      );



  int bispectra_non_separable_init(
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      struct bispectra_workspace_non_separable * pwb
      );


  int bispectra_non_separable_workspace_init(
      struct precision * ppr,
      struct background * pba,
      struct thermo * pth,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      struct bispectra_workspace_non_separable * pwb
      );

  int bispectra_non_separable_workspace_free(
      struct bispectra * pbi,
      struct bispectra_workspace_non_separable * pwb
      );


  int bispectra_non_separable_integrate_over_k3(
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
      );


  int bispectra_non_separable_integrate_over_k2(
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
      );

  int bispectra_non_separable_interpolate_over_k2(
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
      );

  int bispectra_non_separable_integrate_over_k1(
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
      );


  int bispectra_non_separable_interpolate_over_k1(
      struct precision * ppr,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      int index_r,
      int index_l2,
      int index_l3,
      double * integral_splines,
      double * interpolated_integral,
      double * f,
      struct bispectra_workspace_non_separable * pwb
      );


  int bispectra_non_separable_integrate_over_r(
      struct precision * ppr,
      struct background * pba,
      struct perturbs * ppt,
      struct bessels * pbs,
      struct transfers * ptr,
      struct primordial * ppm,
      struct bispectra * pbi,
      double * bispectrum,  /* Output, bispectrum[index_l1l2l3] */
      struct bispectra_workspace_non_separable * pwb
      );


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
       );

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
       );

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
       );

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
       );

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
       );

  int bispectra_store_to_disk(
      struct bispectra * pbi,
      int index_bt
      );

  int bispectra_load_from_disk(
      struct bispectra * pbi,
      int index_bt
      );

  int bispectra_galileon_gradient (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3,
    double pk_1, double pk_2, double pk_3,
    double * out
    );

  int bispectra_galileon_time (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3,
    double pk_1, double pk_2, double pk_3,
    double * out
    );

  int bispectra_local_model (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3, 
    double pk_1, double pk_2, double pk_3, 
    double * out
    );

  int bispectra_equilateral_model (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3, 
    double pk_1, double pk_2, double pk_3, 
    double * out
    );

  int bispectra_orthogonal_model (
    struct primordial * ppm,
    struct bispectra * pbi,
    double k1, double k2, double k3, 
    double pk_1, double pk_2, double pk_3, 
    double * out
    );

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
       );

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
       );

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
       );


  int bispectra_local_squeezed_bispectrum (
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
       );
     
     
  int bispectra_intrinsic_squeezed_bispectrum (
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
       );


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
       ); 

  int bispectra_test_bispectrum (
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
       ); 

  int bispectra_normalisation (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       int l1, int l2, int l3,
       int X, int Y, int Z,
       double * result
       );

  int bispectra_normalisation_positive (
       struct precision * ppr,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       int l1, int l2, int l3,
       int X, int Y, int Z,
       double * result
       );

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
        );


  int bispectra_local_window_function (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int l1, int l2, int l3,
        int X1, int X2, int X3,
        double * result
        );


  int bispectra_intrinsic_window_function (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int l1, int l2, int l3,
        int X, int Y, int Z,
        double * result
        );


  int bispectra_mesh_interpolate (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int index_bt,
        int index_L1,
        double l3, double l2, double l1,
        int X, int Y, int Z,
        struct interpolation_mesh * fine_mesh,
        struct interpolation_mesh * coarse_mesh,
        double * result
        );

  int bispectra_mesh_allocate(
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct interpolation_mesh ******* meshes
        );


  int bispectra_mesh_empty(
        struct bispectra * pbi,
        struct interpolation_mesh ******* meshes
        );


  int bispectra_mesh_free(
        struct bispectra * pbi,
        struct interpolation_mesh ******* meshes
        );


  int bispectra_mesh_create (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int index_bt,
        int X, int Y, int Z,
        window_function_type * window_function,
        int index_l1,
        int l_max,
        double link_length,
        double group_length,
        double soft_coeff,
        int *** grid,
        int **** id,
        struct interpolation_mesh * mesh
        );


  int bispectra_mesh_create_for_all_probes (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        int index_bt,
        int index_l1,
        struct interpolation_mesh ****** meshes
        );


#ifdef __cplusplus
}
#endif

#endif
