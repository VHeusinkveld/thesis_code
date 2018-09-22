#include "navier-stokes/centered.h"

// Physics
double Re; 		// Reynolds number
double vis;		// Dynamic face viscocity

// Grid details
int minlevel, maxlevel; // Min and max grid level 2^n
double err;

// Rotor
double rTdamp;		// Time to start up rotor
double rP;		// Power
double rD, rA;		// Diameter, Area
double rx0, ry0;	// Origin	

// Data analysis 
vertex scalar omega[]; 	// Vorticity
 

// FUNCTIONS
// =================================================================================
 

// Check if gridcell is in fan domain 
bool in_rotor(double x, double y, double Delta, double rx0, double ry0, double rD) {
	bool x_dom, y_dom;

	x_dom = (x + Delta/2. > rx0 - rD/2.) && (x - Delta/2. < rx0 + rD/2.);
	y_dom = (y - Delta/2. < ry0) && (y + Delta/2. > ry0);

	return x_dom && y_dom; 
}


// Code
// =================================================================================
int main() {
	// Adaptivity 
	minlevel = 3;
	maxlevel = 9;
	err= 0.005;
	
	// Grid initialization 
	init_grid(1<<7);
	L0 = 10.;
	origin (0, 0);

	// Rotor details
	rTdamp = 0;
	rP = 100.;
	rD = 1.89+ 0.1*sqrt(2);  
	rx0 = 5.;
	ry0 = 4.89 + 0.1*sqrt(2);
	rA = rD;
	printf("%g\n", ry0);
	run();
}

// Initialisation
event init(t = 0) {

	// Constants
	Re = 500.;
	vis = 1./Re;
	
	// Boundary conditions
	//periodic (right);
	//periodic (bottom);

	// Couple solver to our variables 
	const face vector muc[]={vis,vis};
	mu = muc;
}

event rotor(i=1; i++) {

	refine(in_rotor(x, y, Delta, rx0, ry0, rD) && level < maxlevel);

	foreach () {
		// Checks if gridcell is in the rotor
		if (in_rotor(x, y, Delta, rx0, ry0, rD)) {
			
			// Takes care of rotor edges 
			double d_start = fabs(x - (rx0 - rD/2.));
			double d_end = fabs(x - (rx0 + rD/2.));
			double c=0, c1=0, c2=0;
			
			// TODO solve rotor edges 
			if (d_start < Delta/2.) {
				if ((x - (rx0 - rD/2)) < 0){
					c1 = (Delta/2. - d_start)/Delta;
				} else {
					c1 = (Delta/2. + d_start)/Delta;
				}
			}
			
			if (d_end < Delta/2.){
				if (((rx0 + rD/2) - x) < 0){
					c2 = (Delta/2. - d_end)/Delta;
				} else {
					c2 = (Delta/2. + d_end)/Delta;
				}
			}
				
			c = 1.; //(c1 + c2 > 0) ? c1 + c2 : 1.; // Includes the option of both ends in one cell 

			// Calculate actual addition to the kinetic energy
			double damp = t < rTdamp ? t/rTdamp : 1.; // Linear rotor startup
			
			double temp = pow(u.y[], 3.) - damp*c*2.*rP*dt/rA*Delta/Delta;
			
			if (temp < 0.) {
				u.y[] = -pow(abs(temp), 1./3.);
			} else {
				u.y[] = pow(temp, 1./3.);
			
			}
		}
	} 
}

event adapt(i++) {
	adapt_wavelet((scalar *){u},(double []){err,err},maxlevel,minlevel);
}

// Output visuals
event movies(t += 0.1) {
	scalar lev[];	 	// Grid depth 
	foreach () {
		omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
		lev[] = level;
	}
	boundary ({lev});
	output_ppm (u.x, file = "ppm2mp4 vel_x.mp4", n = 512, linear = true, min = -1, max = 1);
	output_ppm (u.y, file = "ppm2mp4 vel_y.mp4", n = 512, linear = true, min = -1, max = 1);
	output_ppm (omega, file = "ppm2mp4 vort.mp4", n = 512, linear = true); 
	output_ppm (lev, file = "pp2mp4 grid_depth.mp4", n = 512, min = minlevel, max = maxlevel);
}

// Final event
event end(t += 2; t <= 25) {
	printf("i = %d t = %g\n", i, t);
}




