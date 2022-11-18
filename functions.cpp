#include <math.h>
#include <vector>

//This file includes functions from Numerical Recipes, Third
//Edition, we will try it as it is updated from ran1


/* RAN FUNCTION to generate all the random deviates you'll ever need
 * this is our never any doubt generator and meets
 * all the guidelines for avoiding wasteful, overengineered
 * methods - the constructor Ran() takes a single integer argument
 * which becomes the seed for the sequence generated */

struct Ran {

    unsigned long long int u,v,w;
    Ran(unsigned long long int j) : v(4101842887655102017LL), w(1) {
        //Constructor. Call with any integer seed (except v above)
        u = j ^ v; int64();
        v = u; int64();
        w = v; int64();
    }
    inline unsigned long long int int64() {
        //Return 64-bit random integer
        u = u * 2862933555777941757LL + 7046029254386353087LL;
        v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
        w = 4294957665U*(w & 0xffffffff) + (w >> 32);
        unsigned long long int x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
        return (x+v) ^ w;
    }
    inline double doub() { return 5.42101086242752217E-20 * int64() ; }
    //Return random double-precision floating value in the range 0-1
    inline unsigned int int32() {return (unsigned int)int64(); }
    //Return 32-bit random integer.
};

/* DEVIATES FROM OTHER DISTRIBUTIONS NECESSARY FOR PROJECT
 * EXPONENTIAL DEVIATES, NORMAL DEVIATES, POISSON DEVIATES
 * EXPONENTIAL */

// Constructor arguments are beta and a random sequence seed
struct Expondev : Ran {
	double beta;
	Expondev(double bbeta, unsigned long long int i) : Ran(i), beta(bbeta) {}
	double dev() {
		double u;
		do u = doub(); while (u == 0.);
		return -log(u)/beta;
	}
};

// This code is for generating normal deviates by the Box-Muller
// method - consider it for pedagogical use only. Below is a
// significantly faster method for generating normal deviates
struct Normaldev_BM : Ran {
	double mu,sig;
	double storedval;
	Normaldev_BM(double mmu, double ssig, unsigned long long int i)
	: Ran(i), mu(mmu), sig(ssig), storedval(0.) {}
	double dev() {
		double v1,v2,rsq,fac;
		if (storedval == 0.) {
			do {
				v1=2.0*doub()-1.0;
				v2=2.0*doub()-1.0;
				rsq=v1*v1+v2*v2;
			} while (rsq >= 1.0 || rsq == 0.0);
			fac=sqrt(-2.0*log(rsq)/rsq);
			storedval = v1*fac;
			return mu + sig*v2*fac;
		} else {
			fac = storedval;
			storedval = 0.;
			return mu + sig*fac;
		}
	}
};

// Normal Deviates by Ratio-of-Uniforms
// faster method than above

struct Normaldev : Ran {
    double mu, sig;
    Normaldev(double mmu, double ssig, unsigned long long int i)
    : Ran(i), mu(mmu), sig(ssig){}
    double dev() {
        double u,v,x,y,q;
        do {
            u = doub();
            v = 1.7156*(doub()-0.5);
            x = u - 0.449871;
            y = abs(v) + 0.386595;
            q = sqrt(x) + y*(0.19600*y-0.25472*x);
        } while (q > 0.27597 && (q > 0.27846 || sqrt(v) > -4.*log(u)*sqrt(u)));
        return mu + sig*v/u;
    }
};

//CONSTRUCTOR ARGUMENTS ARE LAMBDA AND A RANDOM SEED SEQUENCE
struct Poissondev : Ran {
	double lambda, sqlam, loglam, lamexp, lambold;
    std::vector<double> logfact;
	Poissondev(double llambda, unsigned long long int i) : Ran(i), lambda(llambda),
		logfact(1024,-1.), lambold(-1.) {}
	int dev() {
		double u,u2,v,v2,p,t,lfac;
		int k;
		if (lambda < 5.) {
			if (lambda != lambold) lamexp=exp(-lambda);
			k = -1;
			t=1.;
			do {
				++k;
				t *= doub();
			} while (t > lamexp);
		} else {
			if (lambda != lambold) {
				sqlam = sqrt(lambda);
				loglam = log(lambda);
			}
			for (;;) {
				u = 0.64*doub();
				v = -0.68 + 1.28*doub();
				if (lambda > 13.5) {
					v2 = sqrt(v);
					if (v >= 0.) {if (v2 > 6.5*u*(0.64-u)*(u+0.2)) continue;}
					else {if (v2 > 9.6*u*(0.66-u)*(u+0.07)) continue;}
				}
				k = int(floor(sqlam*(v/u)+lambda+0.5));
				if (k < 0) continue;
				u2 = sqrt(u);
				if (lambda > 13.5) {
					if (v >= 0.) {if (v2 < 15.2*u2*(0.61-u)*(0.8-u)) break;}
					else {if (v2 < 6.76*u2*(0.62-u)*(1.4-u)) break;}
				}
				if (k < 1024) {
					if (logfact[k] < 0.) logfact[k] = gammal(k+1.);
					lfac = logfact[k];
				} else lfac = gammal(k+1.);
				p = sqlam*exp(-lambda + k*loglam - lfac);
				if (u2 < p) break;
			}
		}
		lambold = lambda;
		return k;
	}
	int dev(double llambda) {
		lambda = llambda;
		return dev();
	}
};
