//#include "grid/octree.h" // For 3D
#include "utils.h" 
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "fractions.h"

/*
============================================================================
Declarations
============================================================================
*/

/* Global variables */
int minlevel, maxlevel; // Grid depths
double eps;				// Error in u fields
struct sRotor rot;  	// Rotor details structure 
struct sDiag dia; 		// Diagnostics
scalar fan[];			// Fan volume fraction

scalar b[];		// Buoyancy
scalar * tracers = {b};
face vector av[]; 

/* Functions */
void rotor_init(); 
void rotor_update();
void rotor_coord();
void rotor_forcing();

/* Structures */
struct sRotor {	
	double rampT;			// Time to start up rotor
	double P, Prho;			// Power, powerdensity 
	double R, W, A, V;		// Diameter, Thickness, Area ,Volume
	double x0, y0, z0;		// Origin of rotor
	double theta, phi;		// Polar and Azimuthal angle 
	coord nf, nr;			// Normal vector fan, rotation 
};

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double Wdone;			// Track work done
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Track real rotor volume
};

/*
============================================================================
Main Code, Events
============================================================================
*/

/* Initialisation */
int main() {
    	// Grid variables 
    	init_grid(2<<7);
   	L0 = 5.;
   	X0 = Y0 = Z0 = 0.;
	
	//# Such that momentum is better conserved 
	u.x.refine = refine_linear; 
	u.y.refine = refine_linear;
	#if dimension > 2
		u.z.refine = refine_linear;
	#endif

   	// Initialize physics 
   	rotor_init(); 
	rotor_coord();
	//# Tell basilisk it is a volume field
	fan.prolongation = fraction_refine;
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; // Flux limiter 
	const face vector muc[] = {1/10000, 1/10000};
	mu = muc;
	a = av; // Link acceleration

  	// Adaptivity
  	minlevel = 4; 
  	maxlevel = 8;
  	eps = 0.005;

	// Set boundary conditions
	periodic (left);
    	
	u.t[bottom] = dirichlet(0.);
	b[bottom] = dirichlet(1.);
	b[top] = neumann(1.);

	// Limit maximum time step 
	DT = 0.05;

	// Run the simulation
    	run();
}

/* Initialisation */
event init(t=0){
	foreach() {
		b[] = y + 0.001*noise();
	}
}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
}

/* Diffusion */
event tracer_diffusion(i++){
	diffusion(b, dt, mu);
}

/* Forcing by the rotor */
event forcing(i = 1; i++) {
	rotor_coord();
	rotor_forcing();
}

/* Rotate the rotor */
event rotate(t = rot.rampT; t+=10 ) {
	// Change center  
	rot.x0 += 0;
	rot.y0 += 0;
	rot.z0 += 0;

	// Change angles 
	rot.theta += 0;
	rot.phi += 0;

	rotor_update();
}

/* Progress event */
event end(t += 2; t <= 30) {
	printf("i = %d t = %g\n", i, t);
}

/* Adaptivity function called */
event adapt(i++) {
	adapt_wavelet((scalar *){fan, u},(double []){0.01,eps,eps,eps},maxlevel,minlevel);
}

/* Visualisation */ 
event movies(t += 0.1) {	
    	vertex scalar omega[]; 	// Vorticity
	scalar lev[];	 	// Grid depth
	scalar ekinRho[]; 		// Kinetic energy

	foreach() {
		omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
		ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
		lev[] = level;
	}

	boundary ({lev, omega, ekinRho});

	output_ppm (b, file = "ppm2mp4 buoyancy.mp4", n = 512);
	output_ppm (fan, file = "ppm2mp4 coord_fan.mp4", n = 512, max = 1, min = 0);
	output_ppm (ekinRho, file = "ppm2mp4 ekin.mp4", n = 512);
	output_ppm (u.x, file = "ppm2mp4 vel_x.mp4", n = 512, linear = true, min = -1, max = 1);
	output_ppm (u.y, file = "ppm2mp4 vel_y.mp4", n = 512, linear = true, min = -1, max = 1);
	output_ppm (omega, file = "ppm2mp4 vort.mp4", n = 512, linear = true); 
	output_ppm (lev, file = "pp2mp4 grid_depth.mp4", n = 512, min = minlevel, max = maxlevel);
}

/* Sanity checks */
event sanity (t += 1){
	
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0;
	double tempEkin = 0;	

	foreach(reduction(+:tempVol) reduction(+:tempEkin)) {
		ekin[] = 0.5*rho[]*sq(Delta)*(sq(u.x[]) + sq(u.y[]));
		tempEkin += ekin[];
		#if dimension > 1			
			tempVol += sq(Delta)*fan[];
		#endif
		#if dimension > 2
			tempVol = cube(Delta)*fan[];
		#endif
	}

	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;

	printf("V=%g, Vr=%g, ",rot.V, dia.rotVol);
	printf("Energy: Ek=%g, W=%g, Ek/W=%g, dEk/dW=%g\n", 
		dia.Ekin, dia.Wdone, dia.Ekin/dia.Wdone, 
		(dia.Ekin-dia.EkinOld)/(dia.Wdone-dia.WdoneOld));

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*dia.Wdone;
}


/*
============================================================================
Functions
============================================================================
*/

/* Function returning the sRotor structure */
void rotor_init() {
    
	// Set variables 
    	rot.rampT = 1.;
	rot.R = 0.05;     
	rot.W = 0.01;                      
    	rot.Prho = 50.;
    
   	rot.x0 = L0/2.;
	rot.y0 = 3*L0/4.;
	rot.z0 = 0.;
	
	rot.theta = M_PI/2.;	// Polar angle
	rot.phi = -M_PI/2.;		// Azimuthal angle 

	rotor_update();
}

/* Updating relevant rotor vars */
void rotor_update() {

   	// Set normal vectors 
   	rot.nf.x = sin(rot.theta)*cos(rot.phi);
	rot.nf.y = sin(rot.theta)*sin(rot.phi);
	rot.nf.z = cos(rot.theta);

	rot.nr.x = sin(rot.theta)*cos(rot.phi);
    	rot.nr.y = sin(rot.theta)*sin(rot.phi);
    	rot.nr.z = cos(rot.theta);

    	// Calculate consequences
	#if dimension > 1	
		rot.A = 2*rot.R;
	#endif
	#if dimension > 2    	
		rot.A = sq(rot.R)*M_PI;      
	#endif
               
	rot.V = rot.A*rot.W;
	rot.P = rot.V*rot.Prho;
}


/* Function returning the volume fractions of a fan object */
void rotor_coord() {

      	scalar sph[], plnu[], plnd[];

    	fraction(sph, -sq((x - rot.x0)) - sq((y - rot.y0)) - sq((z - rot.z0)) + sq(rot.R));
    	fraction(plnu,  rot.nr.x*(x - rot.x0) + rot.nr.y*(y - rot.y0) + rot.nr.z*(z - rot.z0) + rot.W/2.);
    	fraction(plnd, -rot.nr.x*(x - rot.x0) - rot.nr.y*(y - rot.y0) - rot.nr.z*(z - rot.z0) + rot.W/2.);	

	foreach () {
    		fan[] = sph[] * plnu[] * plnd[];
   	}
	boundary({fan});
}

 void rotor_forcing(){
	double tempW = 0.;
	double w, wsgn, damp, usgn, utemp;

	foreach(reduction(+:tempW)) {		
		if(fan[] > 0.) {
			foreach_dimension() {

			// Work in respective direction 
			wsgn = sign(rot.nf.x*u.x[]) + (sign(rot.nf.x*u.x[]) == 0)*sign(rot.nf.x);
			damp = rot.rampT > t ? t/rot.rampT : 1.;
			w = wsgn*fan[]*damp*sq(rot.nf.x)*(2./rho[])*(rot.P/rot.V)*dt;
			tempW += 0.5*w*sq(Delta);

			// New kinetic energy
			utemp = sq(u.x[]) + w;

			usgn = 1.*(u.x[] >= 0)*(utemp > 0) +
			    	     -1.*(u.x[] >= 0)*(utemp < 0) +
		 		      1.*(u.x[] <  0)*(utemp < 0) +
				     -1.*(u.x[] <  0)*(utemp > 0); 

			u.x[] = usgn*sqrt(fabs(utemp));
		}
		}
	}
	
	dia.Wdone += tempW;

}


