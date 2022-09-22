#include <math.h> /* Not necessary for ran1 function but it is necessary for the */
                  /* gamma function and the binomial deviates function */

#define PI 3.141592654   /* Only necessary for the Binomial Deviates Function */

/**********************/
/* The ran1 Function  */
/**********************/

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum) {
/* "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added */
/* safeguards.  Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint */
/* values).  Call with idum a negative integer to initialize; thereafter, do not alter idum between */
/* successive deviates in a sequence. RNMX should approximate the largest doubleing value that is less than 1.  */
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   double temp;
   if (*idum <= 0 || !iy) {              /* Initialize */
      if (-(*idum) < 1) *idum=1;         /* be sure to prevent idum = 0                 */
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0;j--) {         /* Load the shuffle table (after 8 warm-ups).  */
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if (*idum < 0) *idum += IM;
         if (j < NTAB) iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;                         /* Start here when not initializing.              */
   *idum=IA*(*idum-k*IQ)-IR*k;           /* Compute idum=(IA*idum) % IM without over-      */
   if (*idum < 0) *idum += IM;           /*         flows by Schrage's method.             */
   j=iy/NDIV;                            /* Will be in the range 0..NTAB-1.                */
   iy=iv[j];                             /* Output previously stored value and refill the  */
   iv[j] = *idum;                        /*         shuffle table.                         */
   if ((temp=AM*iy) > RNMX) return RNMX; /* Because users don't expect endpoint values     */
   else return temp;
}

/*****************************************************************************/
/* The Gamma Function                                                        */
/* From: Nummerical Recipes in C: The Art of Scientific Computing. pg. 214   */
/*****************************************************************************/

double gammln(double xx){      /* returns the value ln[gamma(xx)] for xx>0 */

/* Internal arithmetic will be done in double precision, a nicety that you */
/* can omit if five-figure accuracy is good enough. */

   double x, y, tmp, ser;
   static double cof[6]={76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};
   int j;

   y=x=xx;
   tmp=x+5.5;
   tmp -= (x+0.5)*log(tmp);
   ser=1.000000000190015;
   for (j=0; j<=5; j++) ser += cof[j]/++y;
   return -tmp+log(2.5066282746310005*ser/x);
}

/*****************************************************************************/
/* The Binomial Deviates Function                                            */
/* From: Nummerical Recipes in C: The Art of Scientific Computing. pg. 295   */
/*****************************************************************************/

double bnldev(double pp, int n, long *idum){

/* returns as a doubleing-point number an integer value that is a random deviate drawn */ 
/* form a binomial distribution of n trials each of probability pp, using ran1(idum) */
/* as a source of uniform random deviates. */

        int j;
        static int nold=(-1);
        double am,em,g,angle,p,bnl,sq,t,y;
        static double pold=(-1.0),pc,plog,pclog,en,oldg;
 
        p=(pp <= 0.5 ? pp : 1.0-pp);

        /* The binomial distribution is invariant under changing pp to 1-pp, if we*/
        /* also change the answer to n minus itself; we'll remember to do this below*/
         

        am=n*p;          /* This is the mean of the deviate to be produced. */
        if (n < 25) {    /* Use the direct method while n is not too large. */
                         /* This can require up to 25 calls to ran1. */
                bnl=0.0;
                for (j=1;j<=n;j++)
                        if (ran1(idum) < p) ++bnl;
        } else if (am < 1.0) {     /* If fewer than one event is expected out of 25 */
                                   /* or more trials, then the distribution is quite*/
                                    /* accurately Poisson. Use direct Poisson method.*/
                g=exp(-am);
                t=1.0;
                for (j=0;j<=n;j++) {
                        t *= ran1(idum);
                        if (t < g) break;
                }
                bnl=(j <= n ? j : n);
        } else {                   /* Use the rejection method. */
                if (n != nold) {   /* If n has changed, then compute useful quantities */
                        en=n;
                        oldg=gammln(en+1.0);
                        nold=n;
               } if (p != pold) {  /* If p has changed, then compute useful quantities */
                        pc=1.0-p;
                        plog=log(p);
                        pclog=log(pc);
                        pold=p;
                }
                sq=sqrt(2.0*am*pc);   /* The rejection method with Lorentzian comparison */
                                      /* function. */
                do {
                        do {
                                angle=PI*ran1(idum);
                                y=tan(angle);
                                em=sq*y+am;
                        } while (em < 0.0 || em >= (en+1.0));   /* Reject */
                        em=floor(em);     /* Trick for integer-valued distribution*/            
                        t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
                                -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
                } while (ran1(idum) > t);   /* Reject. This happens about 1.5 times per
                                             * deviate, on average. */
                bnl=em;
        }
        if (p != pp) bnl=n-bnl;        /* Remember to undo the symmetry transformation */
        return bnl;
}


/*****************************************************************************/
/* The Exponential Deviate Function                                          */
/*                                                                           */
/*****************************************************************************/

double expdev(long *idum, double lambda) {
    /* Returns an exponentially distributed, positive, random
     * deviate of mean lambda, using ran1(idum) as the source of
     * uniform deviates */

    double dum;
    do {
        dum=ran1(idum);
    } while (dum == 0.0);
    return -log(1.0 - dum)*lambda;
}

/*****************************************************************************/
/* The Poisson Deviate Function                                              */
/*                                                                           */
/*****************************************************************************/

// #include <math.h>
// #define PI 3.141592654

double poidev(double xm, long *idum) {
	//  double gammln(double xx);
	// double ran1(long *idum);
	static double sq,alxm,g,oldm=(-1.0);
	double em,t,y;

	if (xm < 12.0) {
		if (xm != oldm) {
			oldm=xm;
			g=exp(-xm);
		}
		em = -1;
		t=1.0;
		do {
			++em;
			t *= ran1(idum);
		} while (t > g);
	} else {
		if (xm != oldm) {
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-gammln(xm+1.0);
		}
		do {
			do {
				y=tan(PI*ran1(idum));
				em=sq*y+xm;
			} while (em < 0.0);
			em=floor(em);
			t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
		} while (ran1(idum) > t);
	}
	return em;
}
// #undef PI
