#include "navier-stokes/centered.h"

// Physics
double Re; 		// Reynolds number
double vis;		// Dynamic face viscocity

// Grid details
int minlevel, maxlevel; // Min and max grid level 2^n
double err;

// Data analysis 
vertex scalar omega[]; 	// Vorticity
double e_old;


// FUNCTIONS
// =================================================================================
struct sRotor {	
	double rampT;			// Time to start up rotor
	double P;			// Power
	double D, L, A, V;		// Diameter, Thickness, Area ,Volume
	double x0, y0;			// Origin
	double theta, phi;		// Polar and Azimuthal angle 
	double n[3];			// Normal vector 
};

struct sRotor_edge {
	double cx1, cx2, cy1, cy2, cz1, cz2;
};

// Rotor
struct sRotor r;
struct sRotor_edge re;

void rotor_init() {
	// Time to ramp up to full power 
	r.rampT = 2.;

	// Dimensions
	r.D = 1. + 0.00001*sqrt(2);
	r.L = 0.2;
	r.A = r.D;
	r.V = r.A*r.L;

	// Location and orientation
	r.x0 = 5.;
	r.y0 = 5. + 0.00001*sqrt(2);
	r.theta = 0.;  			// Polar angle
	r.phi = M_PI/3.;		// Azimuthal angle 

	r.n[0] = sin(r.phi);
	r.n[1] = cos(r.phi);
	
	// Power density
	r.P = r.V*0.1;
}

bool rotor_domain(double x, double y, double Delta, struct sRotor r) {
	bool x_dom, y_dom;

	double xt = (x - r.x0)*cos(r.phi) - (y - r.y0)*sin(r.phi);
	double yt = (x - r.x0)*sin(r.phi) + (y - r.y0)*cos(r.phi);

	//printf("xt=%g, yt=%g\n", xt, yt);
	x_dom = (xt + Delta/2. >= -r.D/2.) && (xt - Delta/2. <= r.D/2.);
	y_dom = (yt + Delta/2. >= -r.L/2.) && (yt - Delta/2. <= r.L/2.);
		

	return x_dom && y_dom; 
}
/*
void rotor_edges(struct sRotor_edge re) {

}

double rotor_velocity(vector u, struct sRotor r){
	// Calculate actual addition to the kinetic energy
	double damp = (t < r.rampT) ? t/r.rampT : 1.; // Linear rotor startup
	double temp; 

	temp = pow(u.x, 3.) - damp*c*2.*r.P*sq(r.n[0])/r.V;

	if (temp < 0.) {
		u.x[] = -pow(fabs(temp), 1./3.);
	} else {
		u.x[] = pow(temp, 1./3.);
	}

	temp = pow(u.y[], 3.) - damp*c*2.*r.P*sq(r.n[1])/r.V;

	if (temp < 0.) {
		u.y[] = -pow(fabs(temp), 1./3.);
	} else {
		u.y[] = pow(temp, 1./3.);
	}
}*/


	
// Code
// =================================================================================
int main() {
	// Adaptivity 
	minlevel = 3;
	maxlevel = 8;
	err= 0.01;
	
	// Grid initialization 
	init_grid(1<<8);
	L0 = 10.;
	origin (0, 0);

	rotor_init();
	DT = 0.1;

	run();
}

// Initialisation
event init(t = 0) {

	// Constants
	Re = 30000.;
	vis = 1./Re;
	
	// Boundary conditions
	periodic (right);
	periodic (bottom);

	// Couple solver to our variables 
	const face vector muc[]={vis,vis};
	mu = muc;
}

event diag(t=0; t+=0.5) {
	double e = 0.;
	double rE = r.P*(t < r.rampT ? t/r.rampT : (t-(r.rampT*0.5)));

	foreach(reduction(+:e)) {
		e += 0.5*(u.x[]*u.x[] + u.y[]*u.y[])*sq(Delta);
	}
	
	printf("de=%g, drE=%g, de/drE=%g  \n", e-e_old, r.P*0.5, (e-e_old)/(r.P*0.5));
	e_old = 1.*e;
}

event forcing(i=1; i++) {

	foreach () {
		// Checks if gridcell is in the rotor
		if (rotor_domain(x, y, Delta, r)) {
			
			//rotor_edges();
			//rotor_velocity(u, r);
			double damp = (t < r.rampT) ? t/r.rampT : 1.; // Linear rotor startup
			double temp; 

			temp = pow(u.x[], 3.) + damp*2.*r.P*sq(r.n[0])*Delta/r.V;

			if (temp < 0.) {
				u.x[] = -pow(fabs(temp), 1./3.);
			} else {
				u.x[] = pow(temp, 1./3.);
			}

			temp = pow(u.y[], 3.) + damp*2.*r.P*sq(r.n[1])*Delta/r.V;

			if (temp < 0.) {
				u.y[] = -pow(fabs(temp), 1./3.);
			} else {
				u.y[] = pow(temp, 1./3.);
			}
		
		}
	} 

}

event adapt(i++) {
	adapt_wavelet((scalar *){u},(double []){err,err},maxlevel,minlevel);
	refine(rotor_domain(x, y, 4*Delta, r) && level < maxlevel);

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
event end(t += 2; t <= 10) {
	printf("i = %d t = %g\n", i, t);
}




