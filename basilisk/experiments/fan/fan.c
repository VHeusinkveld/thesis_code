//#include "grid/octree.h" // For 3D
#include "navier-stokes/centered.h"
#include "utils.h" 
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
    
	struct sRotor rot;
	
    // Set variables 
    rot.rampT = 0.;
	rot.R = 0.25  + 0.00001*sqrt(2);     
	rot.W = 0.1   + 0.00001*sqrt(2);                      
    rot.Prho = 1.;
    
   	rot.x0 = rot.y0 = rot.z0 = L0/2.;
	rot.z0 = 0.;
	rot.theta = 0.; // Polar angle
	rot.phi = 0.;	// Azimuthal angle 

   	// Set normal vectors 
   	rot.nr[0] = sin(rot.theta)*cos(rot.phi);
	rot.nr[1] = sin(rot.theta)*sin(rot.phi);
	rot.nr[2] = cos(rot.theta);

	rot.nf[0] = sin(rot.theta)*cos(-rot.phi + M_PI/2.);
    rot.nf[1] = sin(rot.theta)*sin(-rot.phi + M_PI/2.);
    rot.nf[2] = cos(rot.theta);

    // Calculate consequences
    rot.A = 1.*rot.R;                      
	rot.V = rot.A*rot.W;
	rot.P = rot.V*rot.Prho;

	return rot;
}

/* Function returning the volume fractions of a fan object */
void rotor_coord(struct sRotor rot) {
    scalar sph[], plnu[], plnd[];
    fan.prolongation = fraction_refine; // Tell basilisk it is a volume field
    fraction(sph, -sq((x - rot.x0)) - sq((y - rot.y0)) - sq((z - rot.z0)) + sq(rot.R));
    fraction(plnu, rot.nr[0]*(x - rot.x0) + rot.nr[1]*(y - rot.y0) + rot.nr[2]*(z - rot.z0) + rot.W/2.);
    fraction(plnd, -rot.nr[0]*(x - rot.x0) + -rot.nr[1]*(y - rot.y0) + -rot.nr[2]*(z - rot.z0) + rot.W/2.);
    
    foreach () {
    	fan[] = sph[] * plnu[] * plnd[];
    }
        
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
   	struct sRotor rot = rotor_init(); 
	//TODO incompatible types when init type double using type struct sRotor
	rotor_coord(rot);
  	//mu = {0, 0, 0}

  	// Adaptivity
  	minlevel = 7; 
  	maxlevel = 8;
  	eps = 0.05;

  	foreach_dimension() {
   		periodic (right);
    }

    run();
}

/* Forcing by the rotor */
event forcing(i = 1; i++) {
	double ugoal = 1.;
	foreach() {
        if(fan[] > 0.) {
			printf(" yes");
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
    vertex scalar omega[]; 	// Vorticity
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

