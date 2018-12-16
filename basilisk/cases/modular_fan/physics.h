#if dimension == 3
	#include "SGS.h"
#endif 

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .2 	// Kelvin per meter

//#define STRAT(s) gCONST*(INVERSION + gCONST/CP)*s/TREF 	// Stratification 
#define WIND(s) def.wind/0.41*log((s-0.075)/0.1) 	// log Wind profile TODO

//#define strat(s) gCONST*5.*log(sqrt(s)+1)/TREF

#define QFLX (-0.005)
#define BSURF (1.5*b[] - 0.5*b[0, 1])
#define GFLX (-Lambda*(BSURF - bd))
double bd = 0, Lambda = 0.0125;
#define STRAT(s) log(s + a1 + 1.) - log(a1 + 1.) + (QFLX/Lambda + bd)
double a1 = 1;


double crho = 1.;
scalar b[];
scalar * tracers = {b};
	
face vector av[]; 
struct sCase def;

struct sCase {
	double wind;
	double wphi;
};

void init_physics(){
 	def.wind = -0.1;
        def.wphi = 0.;

	b.nodump = false; // TODO

	if(def.wind == 0.){    	
	    u.n[bottom] = dirichlet(0.);
	    u.t[bottom] = dirichlet(0.);

	    u.n[top] = dirichlet(0.);
 	    u.t[top] = neumann(0.);

            periodic (left);

	} else if(fabs(def.wind) > 0.) {
	    //if(def.wind > 0.) {
	    //} else if(def.wind<0.) {
	    //}

            u.n[right] = dirichlet(WIND(y));

	    u.n[left] = dirichlet(WIND(y));
	    
	    u.n[bottom] = dirichlet(0.);
	    u.t[bottom] = dirichlet(0.); //neumann(0.);

	    u.n[top] = dirichlet(0.);
 	    u.t[top] = neumann(0.);

            b[left] = dirichlet(STRAT(y));
            b[right] = dirichlet(STRAT(y));
        }
	
	b[bottom] = BSURF;
	b[top] = dirichlet(STRAT(y));

	#if dimension == 3
		u.r[bottom] = dirichlet(0.); //	neumann(0.);
		u.r[top] = neumann(0.); 
		
		Evis[bottom] = dirichlet(0.);
        	Evis[top] = dirichlet(0.);

		periodic(front);
	#endif  
	foreach() {
		b[] = STRAT(y);
		if(fabs(def.wind) > 0.) {
       	            u.x[] = WIND(y);
		}
	}

    while(adapt_wavelet((scalar *){u,b},(double []){eps,eps,eps,3.*9.81/273},maxlevel,minlevel).nf) {
	foreach() {
		b[] = STRAT(y);
       	        if(fabs(def.wind) > 0.) {
		    u.x[] = WIND(y);
	        }
	}
    }


}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

mgstats mgb;
//face vector muz;
/* Diffusion */
event tracer_diffusion(i++){
    scalar r[];
    /*
    foreach_face(){
	if(y < Y0 + 1E-10){
            muz.x[] = 0.;
   	} else {
	    muz.x[] = 1.*mu.x[];
        }
    }*/
    
    foreach() {
        r[] = 0;
        if (y < Delta)
            r[] = (QFLX + GFLX)/Delta;
	}
    double flx = 0, bt = 0;
    foreach_boundary(bottom){
        flx += (QFLX + GFLX) * Delta;
         bt += BSURF * Delta;
    }
    bt /= L0;
    flx /= L0;
    printf("%g %g %g %d\n", t, flx, bt, i);  
    
    mgb = diffusion(b, dt, mu, r = r);

}
