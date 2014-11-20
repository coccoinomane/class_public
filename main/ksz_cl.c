/** @file ksz_cl.c
 * Guido W. Pettinari, 22.09.2014
 * Compute angular power spectrum for perpendicular kinetic Sunyaev-Zeldovic,
 * using eq. 25 of Munshi, Pettinari, Dixon, Iliev, Coles (to be published in
 * 2014), which is equivalent to Eq. 4 of Ma & Fry (2002). The result is output to
 * stderr. To dump to file, run the program as ./ksz_cl 2> output/ksz_cl.dat.
 */

#include "class.h"

int main(int argc, char **argv) {

  struct precision pr;        /* for precision parameters */
  struct background ba;       /* for cosmological background */
  struct thermo th;           /* for thermodynamics */
  struct perturbs pt;         /* for source functions */
  struct transfers tr;        /* for transfer functions */
  struct primordial pm;       /* for primordial spectra */
  struct spectra sp;          /* for output spectra */
  struct nonlinear nl;        /* for non-linear spectra */
  struct lensing le;          /* for lensed spectra */
  struct output op;           /* for output files */
  ErrorMsg errmsg;            /* for error messages */


  if (input_init_from_arguments(argc, argv,&pr,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le,&op,errmsg) == _FAILURE_) {
    printf("\n\nError running input_init_from_arguments \n=>%s\n",errmsg);
    return _FAILURE_;
  }

  if (background_init(&pr,&ba) == _FAILURE_) {
    printf("\n\nError running background_init \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_init(&pr,&ba,&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_init \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (perturb_init(&pr,&ba,&th,&pt) == _FAILURE_) {
    printf("\n\nError in perturb_init \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (primordial_init(&pr,&pt,&pm) == _FAILURE_) {
    printf("\n\nError in primordial_init \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (nonlinear_init(&pr,&ba,&th,&pt,&pm,&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_init \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (transfer_init(&pr,&ba,&th,&pt,&nl,&tr) == _FAILURE_) {
    printf("\n\nError in transfer_init \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  /* Before computing the kSZ power spectra, check that we have computed
  everything that is needed. Note that we need the CMB because otherwise the
  time sampling of the power spectrum P(k,z) would be too coarse to sample
  reionisation. */
  if (
    (pt.has_cmb == _FALSE_) ||
    (pt.has_pk_ksz == _FALSE_) ||
    (pt.has_scalars == _FALSE_)) {
    printf ("\nERROR: To compute kSZ C_l's, we need both CMB and kSZ! Set 'output=kPk,tCl' in input file.\n");
    return _FAILURE_;
  }

  /* We need to have computed the kSZ from before the start of reionisation and all the way to today */
  if (sp.z_max_pk < 1.1*th.z_reio) {
    printf ("\nERROR: to compute kSZ C_l's, we need P(k,z) from before reionisation; set a z_max_pk>1.1*%g\n", th.z_reio);
    return _FAILURE_;
  }
  
  /* Compute power spectra, including the P(k,z) for kSZ */
  if (spectra_init(&pr,&ba,&pt,&pm,&nl,&tr,&sp) == _FAILURE_) {
    printf("\n\nError in spectra_init \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  /* Compute kSZ angular C_l */
  double * C_l = calloc (tr.l_size[pt.index_md_scalars], sizeof(double));

  if (spectra_compute_cl_ksz(&pr,&ba,&th,&pt,&pm,&nl,&tr,&sp,C_l) == _FAILURE_) {
    printf("\n\nError in spectra_ksz \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }
  
  /* Print kSZ C_l's to stderr */
  fprintf (stderr, "# %15s %17.7s %17.7s\n", "l", "C_l", "D_l");
  for (int index_l=0; index_l < tr.l_size[pt.index_md_scalars]; ++index_l) {
    int l = tr.l[index_l];
    double Cl = C_l[index_l];
    fprintf (stderr, "%17d %17.7g %17.7g\n", l, Cl, l*(l+1)/(2*_PI_)*pow(ba.T_cmb*1e6,2)*Cl);
  }
  
  /* Resume the standard CLASS workflow */
  if (lensing_init(&pr,&pt,&sp,&nl,&le) == _FAILURE_) {
    printf("\n\nError in lensing_init \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (output_init(&ba,&th,&pt,&pm,&tr,&sp,&nl,&le,&op) == _FAILURE_) {
    printf("\n\nError in output_init \n=>%s\n",op.error_message);
    return _FAILURE_;
  }




  /****** all calculations done, now free the structures ******/

  if (lensing_free(&le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le.error_message);
    return _FAILURE_;
  }

  if (spectra_free(&sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp.error_message);
    return _FAILURE_;
  }

  if (transfer_free(&tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr.error_message);
    return _FAILURE_;
  }

  if (nonlinear_free(&nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl.error_message);
    return _FAILURE_;
  }

  if (primordial_free(&pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm.error_message);
    return _FAILURE_;
  }

  if (perturb_free(&pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt.error_message);
    return _FAILURE_;
  }

  if (thermodynamics_free(&th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th.error_message);
    return _FAILURE_;
  }

  if (background_free(&ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba.error_message);
    return _FAILURE_;
  }

  return _SUCCESS_;

}
