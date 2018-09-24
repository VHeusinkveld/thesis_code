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
double e_old;


// FUNCTIONS
// =================================================================================
 

// Check if gridcell is in fan domain 
bool near_rotor(double x, double y, double Delta, double rx0, double ry0, double rD) {
	bool x_dom, y_dom;

	x_dom = (x + Delta/2. >= rx0 - rD/2.) && (x - Delta/2. <= rx0 + rD/2.);
	y_dom = (y - Delta/2. <= ry0) && (y + Delta/2. >= ry0);

	return x_dom && y_dom; 
}



// Code
// =================================================================================
int main() {
	// Adaptivity 
	minlevel = 3;
	maxlevel = 8;
	err= 0.005;
	
	// Grid initialization 
	init_grid(1<<6);
	L0 = 10.;
	origin (0, 0);

	// Rotor details
	rTdamp = 5.;
	rD = 1. + 0.00001*sqrt(2);  
	rx0 = 5.;
	ry0 = 5. + 0.00001*sqrt(2);
	rA = rD;
	rP = rA*0.5;
	
	DT = 0.1;

	run();
}

// Initialisation
event init(t = 0) {

	// Constants
	//Re = 30000.;
	//vis = 1./Re;
	
	// Boundary conditions
	periodic (right);
	periodic (bottom);

	// Couple solver to our variables 
	const face vector muc[]={vis,vis};
	mu = muc;
}

event diag(t=0; t+=0.1) {
	double e = 0.;
	double rE = rP*(t < rTdamp ? t/rTdamp : (t-(rTdamp*0.5)));

	foreach(reduction(+:e)) {
		e += 0.5*(u.x[]*u.x[] + u.y[]*u.y[])*sq(Delta);
	}
	
	printf("de=%g, drE=%g, de/drE=%g  \n", e-e_old, rP*0.1, (e-e_old)/(rP*0.1));
	e_old = 1.*e;
}

event rotor(i=1; i++) {
		
	foreach () {
		// Checks if gridcell is in the rotor
		if (near_rotor(x, y, Delta, rx0, ry0, rD)) {
			
			// Takes care of rotor edges 
			double d_start = fabs(x - (rx0 - rD/2.));
			double d_end = fabs(x - (rx0 + rD/2.));
			double c=0, c1=0, c2=0;
			
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
		
			// Includes the options of 2 ends in one cell					
			c = (c1 + c2 > 0) ? ((c1 > 0 && c2 > 0) ? c1 + c2 - Delta : c1 + c2) : 1.; 

			// Calculate actual addition to the kinetic energy
			double damp = (t < rTdamp) ? t/rTdamp : 1.; // Linear rotor startup
			
			double temp = pow(u.y[], 3.) - damp*c*2.*rP/rA;
				
			if (temp < 0.) {
				u.y[] = -pow(fabs(temp), 1./3.);
			} else {
				u.y[] = pow(temp, 1./3.);
			
			}
		}
	} 
	// TODO is this needed?	
	boundary((scalar *){u});

}

event adapt(i++) {
	adapt_wavelet((scalar *){u},(double []){err,err},maxlevel,minlevel);
	refine(near_rotor(x, y, 4*Delta, rx0, ry0, rD) && level < maxlevel);

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
event end(t += 2; t <= 30) {
	printf("i = %d t = %g\n", i, t);
}




