/** @file bessel.c
 *
 * This module computes the spherical Bessel functions and stores them for 
 * later interpolation.
 *
 * The spherical Bessel functions are important in cosmology because they project
 * quantities defined in 3D space (such as perturbations at the time of
 * recombination) to a 2D spherical surface (such as our sky today). The
 * basical mathematical reason for this lies in the Rayleigh expansion of the plane
 * wave as a sum of spherical waves (http://en.wikipedia.org/wiki/Plane_wave_expansion).
 *
 * In particular, in CLASS and in SONG the Bessel functions are used to build
 * the projection functions that, once convolved with the line of sight sources,
 * allow to compute today's value of the transfer functions. This process goes 
 * under the name of line of sight formalism; see Seljak and Zaldarriaga 1996 for
 * the original paper and sec. 5.5 of http://arxiv.org/abs/1405.2280 for a general
 * derivation. The Bessel functions are also used to project a primordial
 * bispectrum to today's sky (see eqs. 9 and 10 of Fergusson and Shellard 2007).
 *
 * This module is no longer included in CLASS since v2. In newer versions of CLASS,
 * the (hyper)spherical Bessel functions are computed directly in the transfer.c module
 * using a more general procedure suitable for curved universes. In SONG, we have
 * preferred to keep the old bessel.c module because it is more transparent. However,
 * in order to not disrupt CLASS, the l-list is still computed in the transfer module;
 * that's why this module requires ptr as an input. Note that, due to the aforementioned
 * changes in CLASS, now the Bessel functions in this module are used only in SONG's
 * bispectra.c module to solve the bispectrum integral.
 *
 * The following functions can be called from other modules:
 *
 * -# bessel_init() to run the module, anytime after transfer_init()
 * -# bessel_at_x() for computing a value j_l(x) at any x by spline interpolation
 * -# bessel_at_x_linear() for computing a value j_l(x) at any x by linear interpolation
 * -# bessel_convolution() for convolving a function f(k) with j_l(k*r) for a given r
 * -# bessel_free() at the end to free the structure and its content
 *
 * Created by Julien Lesgourgues on the 26.08.2010.
 * Last modified by Guido W. Pettinari on 01.06.2015
 */

#include "bessel.h"


/**
 * Compute and store the spherical Bessel functions in pbs->j.
 *
 * In detail, this function does:
 *
 * -# Determine the list of multipole values pbs->l where to compute the Bessel functions.
 *    Since CLASS v2, the list is copied directly from the transfer module (see comment on
 *    top of this file).
 *
 * -# Allocate the Bessel array pbs->j and its second derivative pbs->ddj in view of
 *    spline interpolation.
 *
 * -# For each multipole value in pbs->l, call bessel_j_for_l() to compute and store
 *    the spherical Bessel functions.
 *
 * The computation of the Bessel function is determined by the following parameters:
 *
 * -# pbs->l[index_l]: list of size pbs->l_size with the multipole values l where we shall
 *    compute the Bessel functions; it is determined in transfer structure.
 * -# pbs->x_step: step dx for sampling Bessel functions j_l(x)
 * -# pbs->x_max: maximum value of x where to compute the Bessel functions
 * -# pbs->j_cut: value of j_l(x) below which it is approximated by zero in the region x<<l
 *
 * You need to call the background_init() and transfer_init() before calling this function.
 */

int bessel_init(
		struct precision * ppr,
    struct background * pba,
    struct thermo * pth,
		struct transfers * ptr,
		struct bessels * pbs
		)
{

  if ((pbs->l_max==0) || !ptr->has_cls) {

    if (pbs->bessels_verbose > 0)
      printf_log("Bessel module skipped.\n");

    return _SUCCESS_;
  }
  else {
    printf_log_if (pbs->bessels_verbose, 0, 
      "Computing bessels\n");
  }


  /** - define local variables */

  /* index for l (since first value of l is always 2, l=index_l+2) */
  int index_l;
  int num_j;
  int abort = _FALSE_;

  /** - update the value of pbs->x_max of j_l(x) */
  
  /* In the bispectrum integral the upper limit in the time variable is r_max rather
  than tau0. Here we scale the Bessels domain so that it includes the integration domain
  of the bispectrum. */
  if (pbs->has_cmb_bispectra) {

    double r_max;

    if (ppr->bispectra_r_sampling == custom_r_sampling)
      r_max = ppr->r_max;
    else if ((ppr->bispectra_r_sampling == centred_r_sampling) || (ppr->bispectra_r_sampling == sources_r_sampling))
      r_max = (pba->conformal_age-pth->tau_rec) + ppr->r_right*pth->tau_rec;

    pbs->x_max = 1.01 * MAX (pbs->x_max, pbs->x_max*r_max/pba->conformal_age);

  }

  /** - copy l values from the transfer module and set x_max */

  class_call(bessel_get_l_list(ppr,ptr,pbs),
	     pbs->error_message,
	     pbs->error_message);

  /** - check x_step, x_max and j_cut from precision parameters and from l_max */

  class_test(pbs->x_step <= 0.,
	     pbs->error_message,
	     "x_step=%e, stop to avoid segmentation fault",pbs->x_step);

  class_test(pbs->x_max <= 0.,
	     pbs->error_message,
	     "x_max=%e, stop to avoid segmentation fault",pbs->x_max);

  pbs->j_cut = ppr->bessel_j_cut;

  /** - do we need to store also j_l'(x) ? */
  // added for new version
  pbs->has_dj = _TRUE_;

  /** - compute Bessel functions */

  class_alloc(pbs->x_size,pbs->l_size*sizeof(int*),pbs->error_message);
  class_alloc(pbs->buffer,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->x_min,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->j,pbs->l_size*sizeof(double*),pbs->error_message);
  class_alloc(pbs->ddj,pbs->l_size*sizeof(double*),pbs->error_message);
  if (pbs->has_dj) {
    class_alloc(pbs->dj,pbs->l_size*sizeof(double*),pbs->error_message);
    class_alloc(pbs->dddj,pbs->l_size*sizeof(double*),pbs->error_message);
  }

  /* initialize error management flag */
  abort = _FALSE_;
  
  /* beginning of parallel region */

  #pragma omp parallel				\
  shared(ppr,pbs,abort)				\
  private(index_l)
  
  {

    #pragma omp for schedule (dynamic)

    /** (a) loop over l and x values, compute \f$ j_l(x) \f$ for each of them */
    for (index_l = 0; index_l < pbs->l_size; index_l++) {

      class_call_parallel(bessel_j_for_l(ppr,pbs,index_l),
			  pbs->error_message,
			  pbs->error_message);
	
      #pragma omp flush(abort)
		  
    } /* end of loop over l */

  } /* end of parallel region */

  if (abort) return _FAILURE_;

  pbs->x_size_max=0;
  for (index_l=0; index_l < pbs->l_size; index_l++)
    if (pbs->x_size[index_l] > pbs->x_size_max)
      pbs->x_size_max=pbs->x_size[index_l];

  return _SUCCESS_;

}


/** 
 * Bessel function for arbitrary argument x using spline interpolation.
 *
 * Evaluates the spherical Bessel function x at a given value of x by
 * interpolating in the pre-computed table.  This function can be
 * called from whatever module at whatever time, provided that
 * bessel_init() has been called before, and bessel_free() has not
 * been called yet.
 */
int bessel_at_x(
		struct bessels * pbs, /**< input, pointer to bessels structure */
		double x, /**< argument x where to interpolate j_l(x) */
		int index_l, /**< order of the Bessel functions from pbs->l[index_l] */
		double * j /* output, spherical Bessel function j_l(x) */
		) {
  
  /** Summary: */

  /** - define local variables */

  int index_x;          /* index in the interpolation table */
  double a;         /* quantities for the splint interpolation formula */

  /** - if index_l is too large to be in the interpolation table, return  an error */

  class_test(index_l > pbs->l_size,
	     pbs->error_message,
	     "index_l=%d>l_size=%d; increase l_max.",index_l,pbs->l_size);

  /** - if x is too small to be in the interpolation table, return 0 */

  if (x < *(pbs->x_min[index_l])) {
    *j=0;
    return _SUCCESS_;
  }
  else {

    /** - if x is too large to be in the interpolation table, return an error (this should never occur since x_max in the table should really be the highest value needed by the code, given the precision parameters) */

    class_test(x > pbs->x_max,
	       pbs->error_message,
	       "x=%e>x_max=%e in bessel structure",x,pbs->x_max);

    /** - otherwise, interpolation is needed: */

    /** (a) find index_x, i.e. the position of x in the table; no complicated algorithm needed, since values are regularly spaced with a known step and known first value */

    index_x = (int)((x-*(pbs->x_min[index_l]))/pbs->x_step);

    /** (b) find result with the splint algorithm (equivalent to the one in numerical recipies, although terms are rearranged differently to minimize number of operations) */
    a = (*(pbs->x_min[index_l])+pbs->x_step*(index_x+1) - x)/pbs->x_step;
    *j= a * pbs->j[index_l][index_x] 
      + (1.-a) * ( pbs->j[index_l][index_x+1]
          - a * ((a+1.) * pbs->ddj[index_l][index_x]
           +(2.-a) * pbs->ddj[index_l][index_x+1]) 
          * pbs->x_step * pbs->x_step / 6.0);

    /* This would be linear interpolation */
    // *j= a * pbs->j[index_l][index_x] 
    //   + (1.-a) * pbs->j[index_l][index_x+1];

  }

  /* Debug - test interpolation */
  // int l = pbs->l[index_l];
  // printf("from table: j(%d, %20.15g) = %g\n", l, (*pbs->x_min[index_l]) + (index_x)*pbs->x_step, pbs->j[index_l][index_x]);
  // printf("from table: j(%d, %20.15g) = %g\n", l, (*pbs->x_min[index_l]) + (index_x+1)*pbs->x_step, pbs->j[index_l][index_x+1]);
  // printf("interp    : j(%d, %20.15g) = %g\n", l, x, *j);  

  return _SUCCESS_;

}



/** 
 * Bessel function for arbitrary argument x using linear interpolation.
 * 
 * Same as bessel_at_x(), but using linear interpolation instead of cubic
 * spline interpolation.
 */
int bessel_at_x_linear(
		struct bessels * pbs, /**< input, pointer to bessels structure */
		double x, /**< argument x where to interpolate j_l(x) */
		int index_l, /**< order of the Bessel functions from pbs->l[index_l] */
		double * j /* output, spherical Bessel function j_l(x) */
    ) {

  /** Summary: */

  /** - define local variables */

  int index_x;          /* index in the interpolation table */
  double a;         /* quantities for the splint interpolation formula */

  /** - if index_l is too large to be in the interpolation table, return  an error */

  class_test(index_l > pbs->l_size,
       pbs->error_message,
       "index_l=%d>l_size=%d; increase l_max.",index_l,pbs->l_size);

  /* Shortcuts */
  double x_min = *(pbs->x_min[index_l]);
  double x_step = pbs->x_step;

  /** - if x is too small to be in the interpolation table, return 0 */

  if (x < x_min) {
    *j = 0;
    return _SUCCESS_;
  }
  else {

    /** - if x is too large to be in the interpolation table, return an error (this should never occur since x_max in the table should really be the highest value needed by the code, given the precision parameters) */

    class_test(x > pbs->x_max,
         pbs->error_message,
         "x=%e>x_max=%e in bessel structure",x,pbs->x_max);

    /** - otherwise, interpolation is needed: */

    /** (a) find index_x, i.e. the position of x in the table; no complicated algorithm needed, since values are regularly spaced with a known step and known first value */

    index_x = (int)((x-x_min)/x_step);

    /** (b) find result with the splint algorithm (equivalent to the one in numerical recipies, although terms are rearranged differently to minimize number of operations) */
    a = (x_min + x_step*(index_x+1) - x)/x_step;

    double * j_x = &(pbs->j[index_l][index_x]);
    *j= a * (*j_x)    +     (1.-a) *  (*(j_x+1));

  }

  /* Some debug */
  // int l = pbs->l[index_l];
  // printf("from table: j(%d, %20.15g) = %g\n", l, (*pbs->x_min[index_l]) + (index_x)*pbs->x_step, pbs->j[index_l][index_x]);
  // printf("from table: j(%d, %20.15g) = %g\n", l, (*pbs->x_min[index_l]) + (index_x+1)*pbs->x_step, pbs->j[index_l][index_x+1]);
  // printf("interp    : j(%d, %20.15g) = %g\n", l, x, *j);

  return _SUCCESS_;

}



/**
 * Given the integration domain kk and the array f[index_k], compute the following integral
 * using the trapezoidal rule:
 * 
 *     /
 *    |  dk k^2 f[k] * g[k] * j_l[k*r]
 *    /
 *
 * The Bessel functions are interpolated from the table in the Bessel structure pbs. The
 * delta_k array should contain the trapezoidal measure around a given k[i], that is
 * delta_k[i]=k[i+1]-k[i-1], with delta_k[0]=k[1]-k[0] and delta_k[k_size-1]=k[k_size-1]-k[k_size-2].
 *
 * If you give a NULL pointer for the g function, then it is assumed to be unity.
 */

int bessel_convolution(
    struct precision * ppr, /**< pointer to precision structure */
    struct bessels * pbs, /**< pointer to Bessel structure, should be already initiated
                          with bessel_init() */
    double * kk, /**< array with the integration grid in k */
    double * delta_kk, /**< trapezoidal measure, compute as delta_k[i]=k[i+1]-k[i-1],
                       and delta_k[0]=k[1]-k[0], delta_k[k_size-1]=k[k_size-1]-k[k_size-2] */
    int k_size, /**< size of the integration grid in k */
    double * f, /**< array with the integrand function f, of size k_size */
    double * g, /**< array with the integrand function g, of size k_size; pass NULL
                to automatically set it to unity */
    int index_l, /**< order of the Bessel function, taken from the multipole array pbs->l */
    double r, /**< frequency of the Bessel function */
    double * integral, /**< output, estimate of the integral */
    ErrorMsg error_message /**< string to write error message */
    )
{

  /* Test that the Bessel functions have been computed for the requested
  multipole index (index_l) and argument (x=k*r) */
  class_test ((index_l<0) || (index_l>=pbs->l_size),
    error_message,
    "index_l=%d out of bounds (l_size=%d)", index_l, pbs->l_size);

  class_test ((kk[k_size-1]*r)>pbs->x_max,
    error_message,
    "r*k_max=%g is larger than x_max=%g (index_l1=%d, k_max=%g, r=%g)",
    kk[k_size-1]*r, pbs->x_max, index_l, kk[k_size-1], r);

  /* Initialize the integral */
  *integral = 0;

  /* Find the value of l from the Bessel structure */
  int l = pbs->l[index_l];
     
  /* Loop over the integration grid */
  for (int index_k = 0; index_k < k_size; ++index_k) {

    /* Value of the function f in k */
    double f_in_k = f[index_k];

    /* If the function f vanishes, do not bother computing the Bessel function,
    and jump to the next iteration without incrementing the integral.  Note that this
    is important because CLASS transfer functions at first-order are set to zero above
    a certain value of k. */
    if (f_in_k == 0.)
      continue;

    /* Same for the g function, if it is defined */
    double g_in_k = 1.;
    
    if (g != NULL) {

      g_in_k = g[index_k];

      if (g_in_k == 0.)
        continue;
    }

    /* Value of the considered k */
    double k = kk[index_k];

    /* - Bessel interpolation j_l(k*r) */
    double j;
    double x = k*r;

    /* j_l(x) vanishers for x < x_min(l) */
    if (x < *(pbs->x_min[index_l]))
      continue;
    
    int index_x = (int)((x-*(pbs->x_min[index_l]))/pbs->x_step);
    double a = (*(pbs->x_min[index_l])+pbs->x_step*(index_x+1) - x)/pbs->x_step;

    /* Store in 'j' the value of j_l(r*k) */
    if (ppr->bessels_interpolation == linear_interpolation) {
      
      j = a * pbs->j[index_l][index_x] 
          + (1.-a) * pbs->j[index_l][index_x+1];

      /* Uncomment the following if you prefer to use the built-in function (slightly slower,
      but performs cheks on the interpolation bounds) */
      // class_call (bessel_at_x_linear (
      //               pbs,
      //               k*r,
      //               index_l,
      //               &j),
      //   pbs->error_message,
      //   error_message);
    }
    else if (ppr->bessels_interpolation == cubic_interpolation) {

      j = a * pbs->j[index_l][index_x]
          + (1.-a) * ( pbs->j[index_l][index_x+1]
            - a * ((a+1.) * pbs->ddj[index_l][index_x]
            +(2.-a) * pbs->ddj[index_l][index_x+1]) 
            * pbs->x_step * pbs->x_step / 6.0);
 
      /* Uncomment the following if you prefer to use the built-in function (slightly slower,
      but performs cheks on the interpolation bounds) */
      // class_call (bessel_at_x (
      //               pbs,
      //               k*r,
      //               index_l,
      //               &j),
      //   pbs->error_message,
      //   error_message);
    }

    /* Value of the integrand function */
    double integrand = k * k * j * f_in_k * g_in_k;

    /* Increment the estimate of the integral */
    *integral += integrand * delta_kk[index_k];
    
  } // end of for(index_k)
   
   
  /* Divide the integral by a factor 1/2 to account for the trapezoidal rule */
  (*integral) *= 0.5;
   
  return _SUCCESS_;
  
} // end of bessel_convolution



/**
 * Free all memory space allocated by bessel_init().
 *
 * To be called at the end of each run.
 *
 * @param pbs Input : Initialized bessel structure 
 * @return the error status
 */

int bessel_free(
		struct bessels * pbs) {

  int index_l;

  if (pbs->l_max > 0) {

    for (index_l = 0; index_l < pbs->l_size; index_l++) {
      free(pbs->buffer[index_l]);
    }
    free(pbs->buffer);
    free(pbs->x_min);
    free(pbs->j);
    free(pbs->ddj);
    if (pbs->has_dj) {
      free(pbs->dj);
      free(pbs->dddj);
    }
    free(pbs->x_size);
    free(pbs->l);
    
  }
  
  return _SUCCESS_; 
}

/**
 * Define number and values of mutipoles l. This is just
 * copied from ptr->l, which in turn is set in 
 * 'transfer_get_l_list'.
 */

int bessel_get_l_list(
		      struct precision * ppr,
		      struct transfers * ptr,
		      struct bessels * pbs
		      )
{

  // ===============================================================================
  // =                            Write list of multipoles                         =
  // ===============================================================================

  /* We shall copy the multipole list from the transfer module. Before doing so, we
  check that that pbs->l_max, which we computed in the input.c module and is used 
  in input2.c, coincides with the l_max from the transfer module. This should be the
  case because pbs->l_max was computed in input_read_parameters() using the same code
  snippet used in transfer_get_l_list(). */
  class_warning (pbs->l_max != ptr->l[ptr->l_size_max-1],
    "l_max computed in input.c is different from l_max in transfer.c");

  pbs->l_size = ptr->l_size_max;

  class_alloc (pbs->l, pbs->l_size*sizeof(int), pbs->error_message);

  for (int index_l=0; index_l < pbs->l_size; ++index_l)
    pbs->l[index_l] = ptr->l[index_l];

  pbs->l_max = pbs->l[pbs->l_size-1];

  return _SUCCESS_;

}


/**
 * Get spherical Bessel functions for given value of l.
 *
 * Find the first value x_min(l) at which the function is not negligible
 * (for large l values, Bessel functions are very close to zero nearly
 * until x=l). Then, sample it with step x_step till x_max, using the
 * function bessel_j().
 */

int bessel_j_for_l(
		   struct precision * ppr,
		   struct bessels * pbs,
		   int index_l
		   ){
  
  /** Summary: */

  /** - define local variables */

  /* index for x and value x=x_min[index_l]+x_step*index_x */
  int index_x;
  double x;

  /* value of j_l(x), j_{l-1}(x) returned by bessel_j(); plus j_l'(x) */
  double j,jm,jprime;

  /* for computing x_min */
  double x_min_up;
  double x_min_down;
  double x_min;

  int num_j;

  index_x=0;
  j = 0.;

  /** - find x_min[index_l] by bisection */

  x_min_up=(double)pbs->l[index_l]+0.5;
  x_min_down=0.;

  class_call(bessel_j(pbs,
               pbs->l[index_l], /* l */
               x_min_up, /* x */
               &j),  /* j_l(x) */
    pbs->error_message,
    pbs->error_message);
  
  class_test(j < pbs->j_cut,
    pbs->error_message,
    "in bisection, wrong initial guess for x_min_up.");
 
  while ((x_min_up-x_min_down)/x_min_down > ppr->bessel_tol_x_min) {
      
    class_test((x_min_up-x_min_down) < ppr->smallest_allowed_variation,
	       pbs->error_message,
	       "(x_min_up-x_min_down) =%e < machine precision : maybe kmin=%e is too small",
	       (x_min_up-x_min_down),ppr->bessel_tol_x_min);
    
    class_call(bessel_j(pbs,
          			pbs->l[index_l], /* l */
          			0.5 * (x_min_up+x_min_down), /* x */
          			&j),  /* j_l(x) */
      pbs->error_message,
      pbs->error_message);
    
    if (j >= pbs->j_cut) 
      x_min_up=0.5 * (x_min_up+x_min_down);
    else
      x_min_down=0.5 * (x_min_up+x_min_down);

  }
  
  x_min = x_min_down;

  class_call(bessel_j(pbs,
    		      pbs->l[index_l], /* l */
    		      x_min, /* x */
    		      &j),  /* j_l(x) */
    pbs->error_message,
    pbs->error_message);

  /** - define number of x values to be stored (one if all values of j_l(x) were negligible for this l) */
  
  if (x_min >= pbs->x_max)
    pbs->x_size[index_l] = 1;
  else
    pbs->x_size[index_l] = (int)((pbs->x_max-x_min)/pbs->x_step) + 1;

  /** - allocate memory for x_min[index_l], j[index_l], ddj[index_l] in such way that they stand in a contiguous memory location */

  if (pbs->has_dj) {
    num_j = 4;
  }
  else {
    num_j = 2;
  }

  class_alloc(pbs->buffer[index_l],
    (1+num_j*pbs->x_size[index_l])*sizeof(double),
    pbs->error_message);
    
  pbs->x_min[index_l] = pbs->buffer[index_l];
  pbs->j[index_l] = pbs->buffer[index_l]+1;
  pbs->ddj[index_l] = pbs->j[index_l] + pbs->x_size[index_l];
  if (pbs->has_dj) {
    pbs->dj[index_l] = pbs->ddj[index_l] + pbs->x_size[index_l];
    pbs->dddj[index_l] = pbs->dj[index_l] + pbs->x_size[index_l];
  }

  /** - case when all values of j_l(x) were negligible for this l*/

  if (x_min >= pbs->x_max) {
    
    *(pbs->x_min[index_l]) = pbs->x_max;
    pbs->j[index_l][0]=0.;
    pbs->ddj[index_l][0]=0.;
    if (pbs->has_dj) {
      pbs->dj[index_l][0]=0.;
      pbs->dddj[index_l][0]=0.;
    }
  }

  /** -otherwise, write first non-negligible value and then loop over x */
  else {

    *(pbs->x_min[index_l]) = x_min;

    pbs->j[index_l][0] = j;

    class_call(bessel_j(pbs,
			pbs->l[index_l]-1, /* l-1 */
			x_min, /* x */
			&jm),  /* j_{l-1}(x) */
	       pbs->error_message,
	       pbs->error_message);

    jprime = jm - (pbs->l[index_l]+1)*j/x_min; /* j_l'=j_{l-1}-(l+1)j_l/x */
    
    pbs->ddj[index_l][0] = - 2./x_min*jprime 
      + (pbs->l[index_l]*(pbs->l[index_l]+1)/x_min/x_min-1.)*j; /* j_l'' = -2/x j_l' + (l(l+1)/x/x-1)*j */

    if (pbs->has_dj) {

      pbs->dj[index_l][0] = jprime; 

      pbs->dddj[index_l][0] = - 2./x_min*pbs->ddj[index_l][0]
      	+ ((pbs->l[index_l]*(pbs->l[index_l]+1)+2)/x_min/x_min-1.)*jprime
      	- 2.*pbs->l[index_l]*(pbs->l[index_l]+1)/x_min/x_min/x_min*j;
    }

    /* loop over other non-negligible values */
    for (index_x=1; index_x < pbs->x_size[index_l]; index_x++) {

      x = *(pbs->x_min[index_l])+index_x*pbs->x_step;

      class_call(bessel_j(pbs,
                   pbs->l[index_l], /* l */
                   x, /* x */
                   &j),  /* j_l(x) */
        pbs->error_message,
        pbs->error_message);

      class_call(bessel_j(pbs,
          			  pbs->l[index_l]-1, /* l-1 */
          			  x, /* x */
          			  &jm),  /* j_{l-1}(x) */
        pbs->error_message,
        pbs->error_message);

      jprime = jm - (pbs->l[index_l]+1)*j/x; /* j_l'=j_{l-1}-(l+1)j_l/x */
    
      pbs->j[index_l][index_x] = j;

      pbs->ddj[index_l][index_x] = - 2./x*jprime 
        + (pbs->l[index_l]*(pbs->l[index_l]+1)/x/x-1.)*j; /* j_l'' = -2/x j_l' + (l(l+1)/x/x-1)*j */

      if (pbs->has_dj) {

        pbs->dj[index_l][index_x] = jprime; 

        	pbs->dddj[index_l][index_x] = - 2./x*pbs->ddj[index_l][index_x]
        	  + ((pbs->l[index_l]*(pbs->l[index_l]+1)+2)/x/x-1.)*jprime
        	  - 2.*pbs->l[index_l]*(pbs->l[index_l]+1)/x/x/x*j;
      }
    }
  }
  
  return _SUCCESS_;
}
    

/**
 * Compute spherical Bessel function \f$ j_l(x) \f$ for a given l and x.
 *
 * Inspired from Numerical Recipies.
 * 
 * @param pbs Input : pointer to bessel structure (used only for error mess ge)
 * @param l   Input: l value
 * @param x   Input: x value
 * @param jl  Output: \f$ j_l(x) \f$ value
 * @return the error status
 */

int bessel_j(
	     struct bessels * pbs,
	     int l,
	     double x,
	     double * jl
	     ) {
    
  double nu,nu2,beta,beta2;
  double x2,sx,sx2,cx;
  double cotb,cot3b,cot6b,secb,sec2b;
  double trigarg,expterm,fl;
  double l3,cosb;
  
  class_test(l < 0,
	     pbs->error_message,
	     " ");

  class_test(x < 0,
	     pbs->error_message,
	     " ");

  fl = (double)l;

  x2 = x*x;

  /************* Use closed form for l<7 **********/

  if (l < 7) {

    sx=sin(x);
    cx=cos(x);

    if(l == 0) {
      if (x > 0.1) *jl=sx/x;
      else *jl=1.-x2/6.*(1.-x2/20.);
      return _SUCCESS_;
    }

    if(l == 1) {
      if (x > 0.2) *jl=(sx/x -cx)/x;
      else *jl=x/3.*(1.-x2/10.*(1.-x2/28.));
      return _SUCCESS_;
    }

    if (l == 2) {
      if (x > 0.3) *jl=(-3.*cx/x-sx*(1.-3./x2))/x;
      else *jl=x2/15.*(1.-x2/14.*(1.-x2/36.));
      return _SUCCESS_;
    }

    if (l == 3) {
      if (x > 0.4) *jl=(cx*(1.-15./x2)-sx*(6.-15./x2)/x)/x;
      else *jl=x*x2/105.*(1.-x2/18.*(1.-x2/44.));
      return _SUCCESS_;
    }

    if (l == 4) {
      if (x > 0.6) *jl=(sx*(1.-45./x2+105./x2/x2) +cx*(10.-105./x2)/x)/x;
      else *jl=x2*x2/945.*(1.-x2/22.*(1.-x2/52.));
      return _SUCCESS_;
    }
    
    if (l == 5) {
      if (x > 1) *jl=(sx*(15.-420./x2+945./x2/x2)/x -cx*(1.0-105./x2+945./x2/x2))/x;
      else *jl=x2*x2*x/10395.*(1.-x2/26.*(1.-x2/60.));
      return _SUCCESS_;
    }

    if (l == 6) {
      if (x > 1) *jl=(sx*(-1.+(210.-(4725.-10395./x2)/x2)/x2)+
		      cx*(-21.+(1260.-10395./x2)/x2)/x)/x;
      else *jl=x2*x2*x2/135135.*(1.-x2/30.*(1.-x2/68.));
      return _SUCCESS_;
    }

  }

  else {

    if (x <= 1.e-40) {
      *jl=0.0;
      return _SUCCESS_;
    }

    nu= fl + 0.5;
    nu2=nu*nu;

    if ((x2/fl) < 0.5) {
      *jl=exp(fl*log(x/nu/2.)+nu*(1-log(2.))-(1.-(1.-3.5/nu2)/nu2/30.)/12./nu)
	/nu*(1.-x2/(4.*nu+4.)*(1.-x2/(8.*nu+16.)*(1.-x2/(12.*nu+36.))));
      return _SUCCESS_;
    }

    if ((fl*fl/x) < 0.5) {

      beta = x - _PI_/2.*(fl+1.);
      *jl = (cos(beta)*(1.-(nu2-0.25)*(nu2-2.25)/8./x2*(1.-(nu2-6.25)*(nu2-12.25)/48./x2))
	     -sin(beta)*(nu2-0.25)/2./x* (1.-(nu2-2.25)*(nu2-6.25)/24./x2*(1.-(nu2-12.25)*(nu2-20.25)/80./x2)) )/x;
      
      return _SUCCESS_;

    }

    l3 = pow(nu,0.325);

    if (x < nu-1.31*l3) {

      cosb=nu/x;
      sx=sqrt(nu2-x2);
      cotb=nu/sx;
      secb=x/nu;
      beta=log(cosb+sx/x);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      expterm=((2.+3.*sec2b)*cot3b/24.
	       - ((4.+sec2b)*sec2b*cot6b/16.
		  + ((16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b/5760.
		     + (32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b/128./nu)*cot6b/nu)/nu)/nu;
      *jl=sqrt(cotb*cosb)/(2.*nu)*exp(-nu*beta+nu/cotb-expterm);

      return _SUCCESS_;

    }

    if (x > nu+1.48*l3) {

      cosb=nu/x;
      sx=sqrt(x2-nu2);
      cotb=nu/sx;
      secb=x/nu;
      beta=acos(cosb);
      cot3b=cotb*cotb*cotb;
      cot6b=cot3b*cot3b;
      sec2b=secb*secb;
      trigarg=nu/cotb-nu*beta-_PI_/4.
	-((2.+3.*sec2b)*cot3b/24.
	  +(16.-(1512.+(3654.+375.*sec2b)*sec2b)*sec2b)*cot3b*cot6b/5760./nu2)/nu;
      expterm=((4.+sec2b)*sec2b*cot6b/16.
	       -(32.+(288.+(232.+13.*sec2b)*sec2b)*sec2b)*sec2b*cot6b*cot6b/128./nu2)/nu2;
      *jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);

      return _SUCCESS_;
    }
    
    /* last possible case */

    beta=x-nu;
    beta2=beta*beta;
    sx=6./x;
    sx2=sx*sx;
    secb=pow(sx,1./3.);
    sec2b=secb*secb;
    *jl=(_GAMMA1_*secb + beta*_GAMMA2_*sec2b
	 -(beta2/18.-1./45.)*beta*sx*secb*_GAMMA1_
	 -((beta2-1.)*beta2/36.+1./420.)*sx*sec2b*_GAMMA2_
	 +(((beta2/1620.-7./3240.)*beta2+1./648.)*beta2-1./8100.)*sx2*secb*_GAMMA1_
	 +(((beta2/4536.-1./810.)*beta2+19./11340.)*beta2-13./28350.)*beta*sx2*sec2b*_GAMMA2_
	 -((((beta2/349920.-1./29160.)*beta2+71./583200.)*beta2-121./874800.)*
	   beta2+7939./224532000.)*beta*sx2*sx*secb*_GAMMA1_)*sqrt(sx)/12./sqrt(_PI_);

    return _SUCCESS_;

  }

  class_test(0==0,
	     pbs->error_message,
	     "value of l=%d or x=%e out of bounds",l,x);

}

