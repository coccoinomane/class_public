/**
 * Commented standalone version of the 'cisi' function of Numerical
 * Recipes, third edition, adapted for CLASS and C. I have quickly tested the
 * function with Mathematica up to the 15th digit in the range x=[1e-6,20],
 * and it works fine. To convert the C++ limits in the template 'numeric_limits'
 * to the equivalent C limits in 'float.h', I have used the C++/C correspondances
 * in http://www.cplusplus.com/reference/limits/numeric_limits/.
 * 
 * Guido W. Pettinari, 1.12.2014
 */

#include "cisi.h"

/* Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is
 * returned as a large negative number and no error message is generated.
 * For x<0 the routine returns Ci(-x) and you must supply the -􏱉i*pi􏰸 yourself.
 */

int cisi (
      double x,
      double * ci,
      double * si,
      ErrorMsg errmsg)
{

  /* Here EULER is Euler’s constant􏰷; PIBY2 is pi/2; TMIN is the dividing line
  between using the series and continued fraction; EPS is the relative error,
  or absolute error near a zero of Ci(x); FPMIN is a number close to the smallest
  representable floating-point number; and BIG is a number near the machine
  overflow limit. */
	double EULER=0.577215664901533, PIBY2=1.570796326794897,
         TMIN=2.0, EPS=DBL_EPSILON,
         FPMIN=DBL_MIN*4.0,
         BIG=DBL_MAX*EPS;
	int MAXIT=200;
	int i,k;
	short odd;
	double a,err,fact,sign,sum,sumc,sums,t,term;
	double complex h,b,c,d,del,cs;
  t = fabs(x);
	class_test (t == 0.0, errmsg, "Ci(0) diverges!");
	if (t > TMIN) { /* Evaluate continued fraction by modified Lentz’s method */
		b=1.0+I*t;
		c=BIG;
		d=h=1.0/b;
		for (i=1;i<MAXIT;i++) {
			a= -i*i;
			b += 2.0;
			d=1.0/(a*d+b); /* Denominators cannot be zero */
			c=b+a/c;
			del=c*d;
			h *= del;
			if (fabs(creal(del)-1.0)+fabs(cimag(del)) <= EPS) break;
		}
		class_test (i >= MAXIT, errmsg, "cf failed in cisi");
		h=(cos(t)-I*sin(t))*h;
		cs= -conj(h)+I*PIBY2;
	} else { /* Evaluate both series simultaneously */
		if (t < sqrt(FPMIN)) { /* Special case: Avoid failure of convergence test because of underflow. */
			sumc=0.0;
			sums=t;
		} else {
			sum=sums=sumc=0.0;
			sign=fact=1.0;
			odd=_TRUE_;
			for (k=1;k<=MAXIT;k++) {
				fact *= t/k;
				term=fact/k;
				sum += sign*term;
				err=term/fabs(sum);
				if (odd==_TRUE_) {
					sign = -sign;
					sums=sum;
					sum=sumc;
				} else {
					sumc=sum;
					sum=sums;
				}
				if (err < EPS) break;
				odd=!odd;
			}
			class_test (k > MAXIT, errmsg, "maxits exceeded in cisi");
		}
		cs=sumc+log(t)+EULER+I*sums;
	}
	if (x < 0.0) cs = conj(cs);
  *ci = creal(cs);
  *si = cimag(cs);

  return _SUCCESS_;

}
