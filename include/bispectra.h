/** @file spectra.h Documented includes for bispectra module */

#ifndef __BISPECTRA__
#define __BISPECTRA__

#include "bessel.h"
#include "lensing.h"
#include "spectra.h"
#include "slatec_3j_C.h"

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


struct bispectra {

  // ===========================================================================================
  // =                                         Flags                                           =
  // ===========================================================================================

  /* Should we compute any bispectra at all? */
  short has_bispectra;

  /* Should we include the lensing effects in the computation of the bispectra? This affects
  only analytical bispectra such as the CMB-lensing and the squeezed approximation for the
  intrinsic bispectrum */
  short include_lensing_effects;
  
  /* Flag that overrides 'include_lensing_effects' for the intrinsic bispectrum. If _FALSE_,
  never use the lensed version of the intrinsic bispectrum or of its approximation  */ 
  short lensed_intrinsic;

  /* Are the Wigner 3j-symbols needed to compute the requested bispectra? */
  short need_3j_symbols;

  // ===========================================================================================
  // =                                  Indices of bispectra                                   =
  // ===========================================================================================

  /* What bispectra types should be compute? */
  short has_local_model;              /* Local model */
  short has_equilateral_model;        /* Equilateral model */  
  short has_orthogonal_model;         /* Orthogonal model */
  short has_galileon_model;           /* Galileon inflation (arXiv:1108.0305) */
  short has_intrinsic;                /* Bispectrum induced by non-linear dynamics */
  short has_intrinsic_squeezed;       /* Squeezed limit for the 2nd-order bispectrum (Creminelli et al. 2012, Bartolo et al. 2012, Lewis 2012) */
  short has_intrinsic_squeezed_unlensed; /* Same as above, but unlensed */
  short has_local_squeezed;           /* Squeezed limit for local model (Guangui et al. 1994) */
  short has_cosine_shape;             /* The above, multiplied by an oscillating function in l */
  short has_cmb_lensing;              /* CMB-lensing bispectrum (Eq. 4.5 of arxiv:1101.2234) */
  short has_cmb_lensing_squeezed;     /* Squeezed limit of the CMB-lensing bispectrum (Eq. 5.20 of arxiv:1101.2234) */
  short has_cmb_lensing_kernel;       /* Squeezed limit of the CMB-lensing bispectrum (kernel only, Eq. 5.20 of arxiv:1101.2234) */
  short has_quadratic_correction;     /* Four-point contribution to the bispectrum */

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
  
  /* Add the bolometric and redshift corrections (in the form of C_l*C_l) to the intrinsic
  bispectrum? Ignore if you are only interested in first-order CLASS, and see function
  bispectra2_add_quadratic_correction() for a proper documentation. */
  short add_quadratic_correction;

  /* Array of functions that associates to a bispetrum index its analytical function. Defined only for
  the analytical bispectra. The four 'threej' variables must be respectively the following 3j-symbols:
  ( l1, l2, l3 )      ( l1, l2, l3 )      ( l1, l2, l3 )    ( l1, l2, l3 )
  (  0,  0,  0 )      (  0, -2,  2 )      (  2,  0, -2 )    ( -2,  2,  0 )
  They are needed only for the 'quadratic' and CMB-lensing bispectra. */
  int (*bispectrum_function[_MAX_NUM_BISPECTRA_]) (
    struct precision * ppr,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi,
    int l1, int l2, int l3,
    int X, int Y, int Z,
    double threej_ratio_20m2,
    double threej_ratio_m220,
    double threej_ratio_0m22,
    double *result
  );

  /* To obtain the Fisher matrix involving non-analytical bispectra, these will need to be interpolated
  in (l1,l2,l3). For bispectra with power on squeezed configurations, like the local and intrinisc ones,
  it is preferable to use a window function for the interpolation. This array contains the window
  function associated to each bispectrum. A NULL function implies that no window function will be used. */
  int (*window_function[_MAX_NUM_BISPECTRA_]) (
    struct precision * ppr,
    struct spectra * psp,
    struct lensing * ple,
    struct bispectra * pbi,
    int l1, int l2, int l3,
    int X, int Y, int Z,
    double threej_ratio_20m2,
    double threej_ratio_m220,
    double threej_ratio_0m22,
    double *result
  );


  // ======================================================================================
  // =                                 Indices of fields                                  =
  // ======================================================================================

  /* The fields considered in this module and in the Fisher one are denoted by indices 'bf'.
  These can be temperature (T), E-polarisation (E), B-polarisation (B), Rayleigh (R)...  */
  short has_bispectra_t;
  short has_bispectra_e;
  short has_bispectra_b;
  short has_bispectra_r;

  int index_bf_t;
  int index_bf_e;
  int index_bf_b;
  int index_bf_r;
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


  // ===========================================================================================
  // =                                       Sampling in l                                     =
  // ===========================================================================================

  int * l;                                /* List of multipole values pbi->l[index_l] */
  int l_size;                             /* Number of l's where we compute the bispectrum */
  int l_max;                              /* Maximum value in pbi->l */
  int full_l_size;                        /* Total number of l's, given by l_max - 2 + 1 */
  
  /* The bispectrum depends on 3 l-values (l1,l2,l3), which must satisfy the triangular inequality:
  |l1-l2| <= l3 <= l1+l2.  We determine the range of l3 within pbi->l using the following two arrays. */
  int ** l_triangular_size;       /* l_triangular_size[index_l1][index_l2] is the number of allowed l3 values in pbi->l given l1 and l2 */
  int ** index_l_triangular_min;  /* index_l_triangular_min[index_l1][index_l2] is the index closest to abs(l1-l2) in pbi->l */
  int ** index_l_triangular_max;  /* index_l_triangular_max[index_l1][index_l2] is the index closest to l1+l2 in pbi->l */

  /* Index associated to a given (l1,l2,l3) configuration. This is obtained assuming l3 satisfies the triangular condition,
    and that l1>=l2>=l3. The array should be indexed as:

      pbi->index_l1_l2_l3[index_l1]
                         [index_l1 - index_l2]
                         [index_l3_max - index_l3]
                         
   where index_l1, index_l2, index_l3 are the indices for l1,l2,l3 in the pbi->l array; 'index_l3' has to be smaller
   than index_l3_max, which is given by:
   
   index_l3_max = MIN (index_l2, pbi->index_l_triangular_max[index_l1][index_l2]) */
  long int *** index_l1_l2_l3;
  
  /* Number of (l1,l2,l3) configurations that will be stored for each bispectrum */
  long int n_independent_configurations;

  /* Multi-array where the bispectra will be stored.  It must be indexed as
    bispectra [index_bt][X][Y][Z][index_l1_l2_l3]
  where X,Y,Z=T,E,B are indexed by the 'pbi->index_bf' indices, and
  'index_l1_l2_l3' is described above. */          
  double ***** bispectra;
  
  /* In computing the CMB-lensing bispectrum, turn the lensing potential C_l to zero on small scales,
  where we cannot trust them. This won't change the result because these C_l's are very small for
  large l's. In CAMB, Antony sets C_l^TP=0 for l>300 and C_l^EP=0 for l>40. */
  int lmax_lensing_corrT;
  int lmax_lensing_corrE;

  // ===========================================================================================
  // =                                   Filter functions                                      =
  // ===========================================================================================

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



  // ==========================================================================================
  // =                                    Other arrays                                        =
  // ==========================================================================================

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


  // ==========================================================================================
  // =                                        Disk IO                                         =
  // ==========================================================================================

  /* Files where the bispectra will be stored (one file for each bispectra type) */
  char bispectra_dir[_FILENAMESIZE_];
  FILE ** bispectra_files;
  char ** bispectra_paths;

  /* File that will keep track how how many bispectra have been succesfully written */
  FILE * bispectra_status_file;
  char bispectra_status_path[_FILENAMESIZE_];


  // ===========================================================================================
  // =                                   Technical parameters                                  =
  // ===========================================================================================
  
  short bispectra_verbose;                   /* Flag regulating the amount of information sent to standard output (none if set to zero) */                                                  
  long int n_total_configurations;           /* Number of (l1,l2,l3) configurations compatible with the triangular condition */
  ErrorMsg error_message;                    /* Zone for writing error messages */
  long int count_allocated_for_bispectra;    /* Number of elements allocated in pbi->bispectra */
  long int count_memorised_for_bispectra;    /* Number of elements memorised in pbi->bispectra */


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
  // double **** integral_over_k1;



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

  int bispectra_free(
       struct precision * ppr,
       struct perturbs * ppt,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi
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

  int bispectra_store_to_disk(
      struct bispectra * pbi,
      int index_bt
      );

  int bispectra_load_from_disk(
      struct bispectra * pbi,
      int index_bt
      );


  int bispectra_interpolate (
      struct bispectra * pbi,
      int index_bt,
      double l1,
      double l2,
      double l3,
      double * interpolated_value
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
       double threej_ratio_20m2,
       double threej_ratio_m220,
       double threej_ratio_0m22,
       double * result
       );


  int bispectra_intrinsic_squeezed_unlensed_bispectrum (
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
       );


  int bispectra_quadratic_bispectrum (
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
       ); 


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
       );


  int bispectra_local_window_function (
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
       );


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
       );
       

  int bispectra_local_window_function_test (
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
       );

#ifdef __cplusplus
}
#endif

#endif
