#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h" // For 3D
#include "view.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"

/* Global variables */
int minlevel, maxlevel;         // Grid depths
double meps, eps;		// Maximum error and error in u fields

char sim_ID[] = "angles";

#include "physics.h"
#include "fan.h"
#include "diagnostics.h"

/* Initialisation */
int main() {
	for(/*CONDITION*/){
    	init_grid(2<<5);
   	L0 = 100.;
   	X0 = Y0 = Z0 = 0.;
	a = av; 

	u.x.refine = refine_linear; 
	u.y.refine = refine_linear;
	#if dimension == 3
		u.z.refine = refine_linear;
	#endif

	fan.prolongation = fraction_refine;
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; // Flux limiter 

  	minlevel = 3; 
  	maxlevel = 9;
  	meps = 10.;
	DT = 10E-5;
        TOLERANCE=10E-6;
	CFL = 0.5;
	
	sim_dir_create();
	out.sim_i++;
    	run();
	}
}

/* Initialisation */
event init(t = 0){
	init_physics();
	init_rotor();
	fan.prolongation=fraction_refine;
	refine (fan[] > 0. && level < maxlevel);
	eps = min(meps, 0.07*rot.cu);
}

event init_change(i=10){
	TOLERANCE=10E-3;
	DT = 0.05;
}

/* Adaptivity function called */
event adapt(i++) {
	adapt_wavelet((scalar *){fan, u},(double []){0.,eps,eps,eps},maxlevel,minlevel);
}

/* Progress event */
event progress(t+=2.) {
	fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}

event end(t=90){
}
