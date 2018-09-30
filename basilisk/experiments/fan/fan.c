//#include "grid/octree.h" // For 3D
#include "navier-stokes/centered.h"
#include "fractions.h"

/*
============================================================================
Global variables 
============================================================================
*/

int minlevel, maxlevel;
double eps;
scalar fan[];

/*
============================================================================
Data structures
============================================================================
*/

struct sRotor {	
	double rampT;			// Time to start up rotor
	double P, Prho;			// Power, powerdensity 
	double R, W, A, V;		// Diameter, Thickness, Area ,Volume
	double x0, y0, z0;		// Origin of rotor
	double theta, phi;		// Polar and Azimuthal angle 
	double nf[3], nr[3];	        // Normal vector fan, rotation 
};

/*
============================================================================
Functions
============================================================================
*/

/* Function returning the sRotor structure */
struct sRotor rotor_init() {
    
	struct sRotor r;
	
    // Set variables 
    r.rampT = 2.;
	r.R = 1.  + 0.00001*sqrt(2);     
	r.W = 0.2 + 0.00001*sqrt(2);                      
    r.Prho = 0.1;
    
    r.x0 = r.y0 = r.z0 = 5.;
	r.theta = M_PI/2.;          // Polar angle
	r.phi = 0.;		            // Azimuthal angle 

    // Set normal vectors 
    r.nr[0] = sin(r.theta)*cos(r.phi);
	r.nr[1] = sin(r.theta)*sin(r.phi);
	r.nr[2] = cos(r.theta);

	r.nf[0] = sin(r.theta)*cos(-r.phi + M_PI/2.);
    r.nf[1] = sin(r.theta)*sin(-r.phi + M_PI/2.);
    r.nf[2] = cos(r.theta);

    // Calculate consequences
    r.A = 1.*r.R;                      
	r.V = r.A*r.W;
	r.P = r.V*r.Prho;

	return r;
}

/* Function returning the volume fractions of a fan object */
scalar rotor_coord(struct sRotor r) {

    scalar fan[], sph[], plnu[], plnd[];
    fan.prolongation = fraction_refine; // Tell basilisk it is a volume field

    fraction(sph, -sq((x - r.x0)) - sq((y - r.y0)) - sq((z - r.z0)) + sq(r.R));
    fraction(plnu, r.nr[0]*(x - r.x0) + r.nr[1]*(y - r.y0) + r.nr[2]*(z - r.z0) + r.W/2.);
    fraction(plnd, -r.nr[0]*(x - r.x0) + -r.nr[1]*(y - r.y0) + -r.nr[2]*(z - r.z0) + r.W/2.);
    
    foreach () {
      fan[] = sph[] * plnu[] * plnd[];
    }
    
	return fan;
}


/*
============================================================================
Main Code, Events
============================================================================
*/

/* Initialisation function */
int main() {

    // Grid variables 
    init_grid(2<<7);
    double L0 = 1.;
    X0 = Y0 = Z0 = 0.;

    // Initialize physics 
    struct sRotor ro = rotor_init(); 
    fan = rotor_coord(ro);
    //mu = {0, 0, 0}

    // Adaptivity
    minlevel = 4; 
    maxlevel = 8;
    eps = 0.05;

    foreach_dimension() {
        periodic (right);
    }
}

/* Forcing by the rotor */
event forcing(i = 1; i++) {
    foreach() {
        if(fan[] > 0.) {
            double ugoal = 1.;
            u.y[] = u.y[] + (ugoal - u.y[])*dt;
        }
    }
}

/* Progress event */
event end(t += 2; t <= 10) {
	printf("i = %d t = %g\n", i, t);
}

/* Adaptivity function called */
event adapt(i++) {
	adapt_wavelet((scalar *){u},(double []){eps,eps},maxlevel,minlevel);
}

/* Visualisation */
event movies(t += 0.1) {
    vertex scalar omega[];  // Vorticity
	scalar lev[];	 	    // Grid depth 
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

