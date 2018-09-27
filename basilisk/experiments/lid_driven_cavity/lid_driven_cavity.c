#include "navier-stokes/centered.h"

// Physics
double Re; 		// Reynolds number
double vis;		// Dynamic face viscocity

// Grid details
int minlevel, maxlevel; // Min and max grid level 2^n
double err;

// Data analysis 
vertex scalar omega[]; 	// Vorticity 

int main() {
	// Adaptivity 
	minlevel = 3;
	maxlevel = 7;
	err= 0.005;
	
	// Grid initialization 
	init_grid(1<<3);
	L0 = 1.;
	origin (-0.5, -0.5);

	run();
}

// Boundary conditions
u.t[top] = dirichlet(1.);
u.n[top] = dirichlet(0.);
u.t[bottom] = dirichlet(0);
u.n[bottom] = dirichlet(0.);
u.t[right] = dirichlet(0);
u.n[right] = dirichlet(0.);
u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);

// Initialisation
event init (t = 0) {

	// Constants
	Re = 500.;
	vis = 1./Re;
	
	// Grid + adaptation
	astats adapting;
	do {
		foreach () {
			u.x[] = 0;
		}
		boundary ((scalar *) {u});
		adapting = adapt_wavelet((scalar *){u},(double[]){err,err},maxlevel,minlevel);
	} while (adapting.nf);

	// Couple solver to our variables 
	const face vector muc[]={vis,vis};
	mu = muc;
}

event adapt (i++) {
	adapt_wavelet((scalar *){u},(double []){err,err},maxlevel,minlevel);
}

// Output visuals
event movies (t += 0.1) {
	scalar lev[];	 	// Grid depth 
	foreach () {
		omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
		lev[] = level;
	}
	boundary ({lev});
	output_ppm (omega, file = "ppm2mp4 vort.mp4", n = 512, linear = true); 
	output_ppm (lev, file = "pp2mp4 grid_depth.mp4", n = 512, min = minlevel, max = maxlevel);
}

// Final event
event end (t += 2; t <= 15) {
	printf("i = %d t = %g\n", i, t);
}


