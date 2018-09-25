#include "navier-stokes/centered.h"

// Physics
double Re; 		// Reynolds number
double vis;		// Dynamic face viscocity

// Grid details
int minlevel, maxlevel; // Min and max grid level 2^n
double err;

// Rotor
struct sRotor r;
struct sRotor_edge re;

// Data analysis 
vertex scalar omega[]; 	// Vorticity
double e_old;


// FUNCTIONS
// =================================================================================
struct sRotor {	
	double rampT;			// Time to start up rotor
	double P;			    // Power
	double D, L, A, V;		// Diameter, Thickness, Area ,Volume
	double x0, y0;			// Origin
	double theta, phi;		// Polar and Azimuthal angle 
	double n[3];			// Normal vector 
};

struct sRotor_edge {
	double cx1, cx2, cy1, cy2, cz1, cz2;
};

void rotor_init() {
	// Time to ramp up to full power 
	r.rampT = 5.;

	// Dimensions
	r.D = 1. + 0.00001*sqrt(2);
	r.L = 0.3.;
	r.A = r.D;
	r.V = r.A*r.L;

	// Location and orientation
	r.x0 = 5.;
	r.y0 = 5. + 0.00001*sqrt(2);
	r.theta = 0.;  					// Polar angle
	r.phi = M_PI/2.;		     	// Azimuthal angle 

	r.n[0] = cos(r.phi);
	r.n[1] = sin(r.phi);
	
	// Power density
	r.P = r.V*0.5;
}

bool rotor_domain(double x, double y, double z, double Delta, struct Rotor Rot) {
	bool x_dom, y_dom;

	double xt = x*cos(r.phi) - y*sin(r.phi);
	double yt = x*sin(r.phi) + y*cos(r.phi);

	x_dom = (xt + Delta/2. >= r.x0 - r.D/2.) && (xt - Delta/2. <= r.x0 + r.D/2.);
	y_dom = (yt + Delta/2. >= r.y0 - r.L/2.) && (yt - Delta/2. <= r.y0 + r.L/2.);
		

	return x_dom && y_dom; 
}
/*
void rotor_edges(struct sRotor_edge re) {

}*/

double rotor_velocity(struct sRotor r){
	// Calculate actual addition to the kinetic energy
	double damp = (t < r.rampT) ? t/r.rampT : 1.; // Linear rotor startup
	double temp_u[2] = {u.x[], u.y[]};

	for (i = 0; i <= 1; i++) {
		double temp = pow(temp_u[i], 3.) - damp*c*2.*r.P*sq(r.n[i])/r.V;

		if (temp < 0.) {
			temp_u[i] = -pow(fabs(temp), 1./3.);
		} else {
			temp_u[i] = pow(temp, 1./3.);
		}
	}

	u.x[] = temp_u[0];
	u.y[] = temp_u[1];
}


	
// Code
// =================================================================================
int main() {
	// Adaptivity 
	minlevel = 3;
	maxlevel = 7;
	err= 0.01;
	
	// Grid initialization 
	init_grid(1<<6);
	L0 = 10.;
	origin (0, 0);

	rotor_init()
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
		if (rotor_domain(x, y, Delta, Rot)) {
			//rotor_edges();
			rotor_velocity();
		}
	} 

}

event adapt(i++) {
	adapt_wavelet((scalar *){u},(double []){err,err},maxlevel,minlevel);
	refine(rotor_domain(x, y, 4*Delta, Rot) && level < maxlevel);

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




