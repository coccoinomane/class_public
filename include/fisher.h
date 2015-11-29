/** @file fisher.h Documented includes for Fisher module */

#ifndef __FISHER__
#define __FISHER__

#include "bispectra.h"


/* Maximum number of frequency bands of the experiment, for the purpose of Fisher matrix computation. */
#define _N_FREQUENCY_CHANNELS_MAX_ 100



/**
 *
 * Note that the f_NL we use is the one for the gravitational potential during matter domination, which is
 * related to the comoving curvature perturbation by a factor -3/5 (fnl_R = -3/5 fnl_psi). 
 *
 */

struct fisher {


  // =============================================================================
  // =                             Flags and indices                             =
  // =============================================================================

  short has_fisher; /**< Should we compute the Fisher matrix at all? If _TRUE_,
                    the Fisher module will be skipped altogether. */

  /**
   * Should we include the lensing effects in the Fisher matrix estimator?
   *
   * Gravitational lensing degrades the bispectrum signal via the lensing noise
   * and the lensing variance. If this flag is turned on, we include these effects
   * in the Fisher matrix. 
   *
   * By lensing noise we refer to using the lensed C_l rather than the unlensed
   * ones in the inverse covariance matrix (ie. the C_l at the denominator of
   * the Fisher matrix estimator).
   *
   * By lensing variance we refer to the effect described in Sec. 5 of
   * http://uk.arxiv.org/abs/1101.2234.
   */
  short include_lensing_effects;

  /**
   * Should we compute the Fisher matrix for all the l up to l_max, rather than
   * only for l_max?
   *
   * When the lensing variance is not requested the flag is not needed, as the
   * computation of the Fisher matrix naturally yields F(l) up to l=l_max. With
   * lensing variance, however, computing F(l) requires running the Fisher module
   * separately for each l, because the lensing variance algorithm used in SONG
   * (Lewis et al 2011, Sec. 5) includes a sum over l1 (smallest multipole in
   * the Fisher summation) rather than over l3 (larger multipole in the Fisher
   * summation).
   */
  short compute_lensing_variance_lmax;

  // ===========================================================================================
  // =                               Indices of Fisher matrix                                  =
  // ===========================================================================================
  
  /* Indices of the bispectra types to be included in the Fisher matrix */
  int index_ft_local;                 /* Index for the bispectrum for a local model */
  int index_ft_equilateral;           /* Index for the bispectrum for a equilateral model */
  int index_ft_orthogonal;            /* Index for the bispectrum for a orthogonal model */
  int index_ft_galileon_gradient;     /* Index for the bispectrum for the pi_dot*pi_grad^2 term in Galileon inflation */
  int index_ft_galileon_time;         /* Index for the bispectrum for the pi_dot^3 term in Galileon inflation */
  int index_ft_intrinsic;             /* Index for the bispectrum induced by nonlinear dynamics */
  int index_ft_intrinsic_squeezed;    /* Index for the intrinsic bispectrum in the squeezed limit */  
  int index_ft_intrinsic_squeezed_unlensed; /* Index for the unlensed intrinsic bispectrum in the squeezed limit */  
  int index_ft_local_squeezed;        /* Index for the local-model bispectrum in the squeezed limit */  
  int index_ft_cosine;                /* Index for the oscillating bispectrum */  
  int index_ft_cmb_lensing;           /* Index for the bispectrum of CMB-lensing */  
  int index_ft_cmb_lensing_squeezed;  /* Index for the bispectrum of CMB-lensing in the squeezed limit */  
  int index_ft_cmb_lensing_kernel;    /* Index for the bispectrum of CMB-lensing in the squeezed limit (kernel only) */
  int ft_size;                    /* Total number of bispectra types requested */
  
  /* Correspondence between rows of the Fisher matrix and the bispectra stored in pbi->bispectra[index_bt] */
  int index_bt_of_ft[_MAX_NUM_BISPECTRA_];

  /* Arrays that relate the bispectra in the Fisher module to quantities computed in the bispectrum module
  See documentation for the equivalent arrays in bispectra.h */
  char ft_labels[_MAX_NUM_BISPECTRA_][_MAX_LENGTH_LABEL_];

  // ======================================================================================
  // =                                 Indices of fields                                  =
  // ======================================================================================

  /* Which fields to include in the bispectrum analysis? By default use all that were computed
  in the bispectrum module, unless they are explicitly ignored by the pfi->ignore_bf flags. */
  short has_fisher_t;
  short has_fisher_e;
  short has_fisher_b;

  /* Should we include in the Fisher matrix analysis all fields (T,E,B...) that were 
  computed in the bispectrum module? This is the case by default, but sometimes
  it is useful to restrict the analysis to only to one of them, for example to see the separate
  effect of T and E in a run that has both. */
  short ignore_t;
  short ignore_e;
  short ignore_b;

  /* The fields considered in this module and in the Fisher one are denoted by indices 'ff'.
  These can be temperature (T), E-polarisation (E), B-polarisation (B), Rayleigh (R)...  */
  int index_ff_t;
  int index_ff_e;
  int index_ff_b;
  int ff_size;

  /* Actual number of bispectra to be included in the Fisher matrix (usually 'ff_size' to
  the power of 3, it is 8 for T and E: TTT, TTE, TET, ETT, EET, ETE, TEE, EEE) */
  int n_probes;

  /* Correspondence between the fields included in the Fisher matrix and the bispectra stored in
  pbi->bispectra[index_bt][index_bf_X][index_bf_Y][index_bf_Z] */
  int index_bf_of_ff[_MAX_NUM_FIELDS_];

  /* Arrays that relate the fields in the Fisher module (T,E,B...) to quantities computed throughout
  the code. See documentation for the equivalent arrays in bispectra.h */
  int index_ct_of_phi_ff[_MAX_NUM_FIELDS_];
  int index_ct_of_ff_ff[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_];
  int index_lt_of_ff_ff[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_];
  char ff_labels[_MAX_NUM_FIELDS_][_MAX_LENGTH_LABEL_]; /* T,E... */
  char ffff_labels[_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_][_MAX_NUM_FIELDS_][_MAX_LENGTH_LABEL_]; /* TTT,EEE,TTE... */
  

  // ===============================================================================
  // =                                 Arrays                                      =
  // ===============================================================================
  
  int l_min_estimator;   /**< Minimum l in the estimator sum, default is pbi->l_min */
  int l_max_estimator;   /**< Maximum l in the estimator sum, default is pbi->l_max */

  /**
   * Limits of the Fisher matrix sum for the l1, l2 and l3 multipoles.
   *
   * These limits are never used unless the user changes them manually
   * in fisher_indices(), ot the l_min_estimator and l_max_estimator
   * variables are set.
   */
  int l1_min_global;
  int l2_min_global;
  int l3_min_global;
  int l1_max_global;
  int l2_max_global;
  int l3_max_global;

  /**
   * Grid in l over which we shall compute the Fisher matrix sum.
   *
   * The Fisher matrix is obtained as a sum over (l1,l2,l3). Having computed the
   * bispectrum on a grid in (l1,l2,l3), we need to resort to some interpolation or
   * integration scheme in order to evaluate the sum. Depending on the adopted
   * method, we consider different grids for the three l-directions, and store
   * them in pfi->l1, pfi->l2 and pfi->l3.
   */
  //@{
  int * l1;           /**< Integration grid for the Fisher matrix in the l1 direction */
  int l1_is_full;     /**< If true, l1 runs from pbi->l_min to pbi->l_max; if false, it only contains the nodes in pbi->l */
  int l1_size;        /**< Number of grid points in the l1 direction */

  int * l2;           /**< Integration grid for the Fisher matrix in the l2 direction */
  int l2_is_full;     /**< If true, l2 runs from pbi->l_min to pbi->l_max; if false, it only contains the nodes in pbi->l */
  int l2_size;        /**< Number of grid points in the l2 direction */

  int * l3;           /**< Integration grid for the Fisher matrix in the l3 direction */
  int l3_is_full;     /**< If true, l3 runs from pbi->l_min to pbi->l_max; if false, it only contains the nodes in pbi->l */
  int l3_size;        /**< Number of grid points in the l3 direction */
  //@}


    
  /* Contribution to the Fisher matrix coming from a given l1 and for a given XYZ bispectrum,
  where XYZ=TTT,TTE,TET, etc. This is the sum over l2, l3, A, B, C of
  b^XYZ(l1,l2,l3) * b^ABC(l1,l2,l3) * cov^XYZABC(l1,l2,l3), with l1>=l2>=l3.
  Indexed as pfi->fisher_matrix_XYZ_largest[X][Y][Z][index_l1][index_bt_1][index_bt_2],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ****** fisher_matrix_XYZ_largest;

  /* Same as above, but for l3, the smallest multipole, and l3 belonging to pfi->l3[index_l3]. */
  double ****** fisher_matrix_XYZ_smallest;

  /* Same as fisher_matrix_XYZ_largest, but summed over XYZ */
  double *** fisher_matrix_largest;

  /* Same as fisher_matrix_XYZ_smallest, but summed over XYZ */
  double *** fisher_matrix_smallest;

  /* Fisher matrix for the considered experiment, as a function of the angular resolution
  and for a given bispectrum XYZ. This is obtained as sum_{lmin<=l1<=lmax} fisher_matrix_XYZ_largest,
  with lmin fixed (=2) and lmax varying.
  Indexed as pfi->fisher_matrix_XYZ_lmax[X][Y][Z][index_l1][index_bt_1][index_bt_2]
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ****** fisher_matrix_XYZ_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double ****** fisher_matrix_XYZ_lmin; 

  /* Fisher matrix for the considered experiment, as a function of the angular resolution.
  It is obtained as sum_{lmin<=l1<=l_max,XYZ} fisher_matrix_XYZ_largest.
  Indexed as pbi->fisher_matrix_lmax[index_l1][index_bt_1][index_bt_2],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double *** fisher_matrix_lmax;
  double *** inverse_fisher_matrix_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double *** fisher_matrix_lmin;        
  double *** inverse_fisher_matrix_lmin;

  /* Array that contains 1/sqrt(F^ii), with i=1,..,pbi->bt_size. For a given bispectrum type,
  it corresponds to the value of f_NL that could be detected by an experiment with a
  resolution of lmin<=l<=l_max, with lmin fixed (=2) and lmax varying.
  Indexed in the same way as fisher_matrix_lmax. 
  Indexed as pfi->sigma_fnl_lmax[index_l1][index_bt],
  where index_l1 refers to the multipole pfi->l1[index_l1]. */
  double ** sigma_fnl_lmax;
  
  /* Same as above, but with lmin varying and lmax fixed, and l belonging to pfi->l3[index_l3] */
  double ** sigma_fnl_lmin;
    
    
  /* Array containing the quantity I_l1_l2_l3 in eq. 13 of Komatsu, Spergel & Wandelt (2005):
  
       I_l1_l2_l3 = sqrt( (2L1+1)(2L2+1)(2L3+1)/4*pi ) *  (L1 L2 L3)
                                                          (0  0  0 )
  
  which is needed to compute the Fisher matrix. (This is the factor that converts a reduced
  bispectrum to an angular averaged one.)  */
  double * I_l1_l2_l3;


  /* Cross-power spectrum of the C_l's, and its inverse; it is needed to compute the covariance
  matrix between the various bispectra. For example, if we consider temperature and polarisation
  bispectra, the full covariance matrix is a 8x8 matrix and the cross-power spectrum is given by
  C = ( C_l^TT C_l^TE )
      ( C_l^TE C_l^EE ).
  The C and inverse_C arrays are indexed as pfi->C[l-2][pfi->index_fp_X][pfi->index_fp_Y],
  where X and Y are the considered probes (T and E for temperature and polarisation) */
  double *** C;
  double *** inverse_C;


  // ==========================================================================================
  // =                                   Lensing variance                                     =
  // ==========================================================================================

  /* Same as the other arrays defined above, but used to contain the full result including lensing
  variance */
  double *** fisher_matrix_lensvar_smallest;
  double *** fisher_matrix_lensvar_lmin;
  double *** inverse_fisher_matrix_lensvar_lmin;
  double ** sigma_fnl_lensvar_lmin;
  double *** fisher_matrix_lensvar_lmax; /* Indexed with index_l1 rather than index_l3! */
  double *** inverse_fisher_matrix_lensvar_lmax; /* Indexed with index_l1 rather than index_l3! */
  double ** sigma_fnl_lensvar_lmax; /* Indexed with index_l1 rather than index_l3! */
  
  /* Same as fisher_matrix_XYZ_smallest, but keeping track of the Z and C field indices instead.
  This is needed to compute the lensing variance, and corresponds to \bar{F}_{l_1 i p}
  in Eq. 5.25 of http://uk.arxiv.org/abs/1101.2234. The indexing of this array is slightly
  different from the others, because we will need to invert it with respect to the last two
  levels (see Eq. 5.35 ibidem):
  fisher_matrix_CZ_smallest[index_l3]
                           [index_ft_1*field_size+index_field_C]
                           [index_ft_2*field_size+index_field_Z] */
  double *** fisher_matrix_CZ_smallest;

  /* Same as above but contains also an l3-level:
  fisher_matrix_CZ_smallest_largest[index_l1]
                                   [index_l3]
                                   [index_ft_1*field_size+index_field_C]
                                   [index_ft_2*field_size+index_field_Z] */
  double **** fisher_matrix_CZ_smallest_largest;


  // ==========================================================================================
  // =                                       Noise model                                      =
  // ==========================================================================================

  /* Beam for each frequency band of the considered experiment With respect to Table I of astro-ph/0506396v2,
    'beam[index_channel]' is theta_fwhm in radians for that frequency channel. */
  int n_channels;
  double beam[_N_FREQUENCY_CHANNELS_MAX_];
  
  /* Amplitude of the noise. With respect to Table I of astro-ph/0506396v2, 'noise' is 'sigma' */
  double noise_t[_N_FREQUENCY_CHANNELS_MAX_];
  double noise_e[_N_FREQUENCY_CHANNELS_MAX_];

  /* Total noise as a function of l, defined by C_l_experiment = C_l_theory + N_l. This includes co-added
  contributions from all frequency channels, as explained in eq. 29 of astro-ph/0506396v2.
  It is indexed as pfi->N_l[pbi->index_bf][l-2], where pbi->index_bt=T,E... */
  double ** N_l;

  /* Sky coverage of the experiment. Equal to 1 for a full-sky experiment. */
  double f_sky;

  /* When larger than one, forces SONG to only consider squeezed triangles when computing the
  Fisher matrix for all bispectra. The squeezed triangles are chosen so that
  L2/L1 > squeezed_ratio. Since L1<=L2<=L3, this forces also L3/L1 to be larger than 'squeezed_ratio'.
  When smaller than minus 1, SONG will only consider equilateral triangles with
  L3/L1 < abs(squeezed_ratio). Therefore, the larger and more positive 'squeezed_ratio' is,
  the more squeezed the considered triangles are. When 'squeezed_ratio' is negative, the triangles
  will be more and more equilateral as 'squeezed_ratio' approaches -1. */
  double squeezed_ratio;


  // ===========================================================================================
  // =                                    Technical parameter                                  =
  // ===========================================================================================
  
  short fisher_verbose;            /* Flag regulating the amount of information sent to standard output (none if set to zero) */                                                    
  char info[8192];                 /* Store Fisher matrix information to be printed out to screen and saved to file */
  char info_lensvar[8192];         /* Same, also accounting for lensing-induced noise */
  ErrorMsg error_message;          /* Zone for writing error messages */

};




/**
 * Variables and arrays needed to compute the f_NL estimator given a pair of bispectra.
 *
 */
struct fisher_workspace {
  
  /* Weights for the linear interpolation. Used for the l1 and l2 sums in the
  Fisher matrix estimator. */
  double * delta_l;
  
  /* Weights for the linear interpolation along the triangular direction (l3).
  Indexed as pfi->delta_l3[thread][index_l3] */
  double ** delta_l3;

  /* Temporary arrays needed to store the 3j-symbols */
  double ** threej_000;
  double ** threej_m220;
  double ** threej_0m22;
  double ** threej_20m2;

};








// /*************************************************************************************************************/

/*
 * Boilerplate for C++ 
 */
#ifdef __cplusplus
extern "C" {
#endif


  int fisher_init(
       struct precision * ppr,
       struct background * pba,
       struct thermo *pth,
       struct perturbs * ppt,
       struct bessels * pbs,
       struct transfers * ptr,
       struct primordial * ppm,
       struct spectra * psp,
       struct lensing * ple,
       struct bispectra * pbi,
       struct fisher * pfi
       );

  int fisher_free(
       struct bispectra * pbi,
       struct fisher * pfi
       );

  int fisher_indices(
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
          );

  int fisher_cross_cls(
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
          );

  int fisher_noise(
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
          );
          
          
  int fisher_compute(
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
        );

  int fisher_compute_matrix (
        struct precision * ppr,
        struct transfers * ptr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        struct fisher_workspace * pw
        );

  int fisher_compute_matrix_nodes (
        struct precision * ppr,
        struct spectra * psp,
        struct lensing * ple,
        struct bispectra * pbi,
        struct fisher * pfi,
        struct fisher_workspace * pw
        );
        
  int fisher_sky_coverage (
         struct fisher * pfi
         );

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
          struct fisher_workspace * pw
          );

#ifdef __cplusplus
}
#endif

#endif
