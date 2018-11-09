/** Include required libraries */ 
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "grid/octree.h" 		// For 3D
#include "view.h"			// For bview
#include "navier-stokes/centered.h"     // Navier stokes 
#include "tracer.h"			// Tracers
#include "diffusion.h"			// Diffusion 

/** Global variables */
int minlevel, maxlevel;         	// Grid depths
double meps, eps;			// Maximum error and error in u fields

char sim_ID[] = "rotation";		// Simulation identifier
char sim_var[] = "None";  		// Notes if a variable is varied over runs

#include "physics.h"			// Physics of the simulation 
#include "fan.h"			// Include a fan
#include "diagnostics.h"		// Perform diagnostics

/** Initialisation */
int main() {	
	minlevel = 5;
  	maxlevel = 9;

   	L0 = 100.;
   	X0 = Y0 = Z0 = 0.;

	// Possibility to run for variable changes
	for(minlevel = 5; minlevel<6; minlevel+=1)
	{
    	init_grid(2<<(minlevel-1));
	a = av; 

	u.x.refine = refine_linear; 			// Momentum conserved 
	u.y.refine = refine_linear;
	#if dimension == 3
		u.z.refine = refine_linear;
	#endif

	fan.prolongation = fraction_refine;		// Fan is a volume fraction
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; 				// Flux limiter 

  	meps = 10.;					// Maximum adaptivity criterion
	DT = 10E-5;					// For poisson solver 
        TOLERANCE=10E-6;				// For poisson solver 
	CFL = 0.5;					// CFL condition

	sim_dir_create();				// Create relevant dir's
	out.sim_i++;					// Simulation iteration 
    	run();						// Start simulation 
	}
}

/** Initialisation */
event init(t = 0){
	rot.rotate = true;
	init_physics();
	init_rotor();
	fan.prolongation=fraction_refine;
	refine (fan[] > 0. && level < maxlevel);
	eps = min(meps, 0.06*rot.cu);
}

/** Return to standard tolerances and DTs for poisson solver */ 
event init_change(i=10){
	TOLERANCE=10E-3;
	DT = 0.05;
}

/** Adaptivity */
event adapt(i++) {
	adapt_wavelet((scalar *){fan, u},(double []){0.,eps,eps,eps},maxlevel,minlevel);
}

/** Progress event */
event progress(t+=2.) {
	fprintf(stderr, "i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
}

/** End the simulation */
event end(t=360){
}
