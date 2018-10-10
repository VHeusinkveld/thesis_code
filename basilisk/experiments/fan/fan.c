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

/* Functions */
void rotor_init(struct sRotor*); 
void rotor_update(struct sRotor*);
void rotor_coord(struct sRotor*, scalar*);
void rotor_forcing(scalar*, vector*, struct sDiag*, struct sRotor*);

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
   	rotor_init(&rot); 
	rotor_coord(&rot, &fan);

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
  	eps = 0.01;

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
		b[] = 9.81*(0.5*y + 273 - 273)/273 + 0.001*noise();
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
	rotor_coord(&rot, &fan);
	rotor_forcing(&fan, &u, &dia, &rot);
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

	rotor_update(&rot);
}

/* Progress event */
event end(t += 2; t <= 30) {
	printf("i = %d t = %g\n", i, t);
}

/* Adaptivity function called */
event adapt(i++) {
	adapt_wavelet((scalar *){fan, u},(double []){0.0001,eps,eps,eps},maxlevel,minlevel);
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
void rotor_init(struct sRotor* sR) {
    
	// Set variables 
    	sR.rampT = 1.;
	sR.R = 0.05;     
	sR.W = 0.01;                      
    	sR.Prho = 5.;
    
   	sR.x0 = L0/2.;
	sR.y0 = 3*L0/4.;
	sR.z0 = 0.;
	
	sR.theta = M_PI/2.;	// Polar angle
	sR.phi = -M_PI/2.;		// Azimuthal angle 

	sR = rotor_update(sR);
}

/* Updating relevant rotor vars */
void rotor_update(struct sRotor* sR) {

   	// Set normal vectors 
   	sR.nf.x = sin(sR.theta)*cos(sR.phi);
	sR.nf.y = sin(sR.theta)*sin(sR.phi);
	sR.nf.z = cos(sR.theta);

	sR.nr.x = sin(sR.theta)*cos(sR.phi);
    	sR.nr.y = sin(sR.theta)*sin(sR.phi);
    	sR.nr.z = cos(sR.theta);

    	// Calculate consequences
	#if dimension > 1	
		sR.A = 2*sR.R;
	#endif
	#if dimension > 2    	
		sR.A = sq(sR.R)*M_PI;      
	#endif
               
	sR.V = sR.A*sR.W;
	sR.P = sR.V*sR.Prho;
}


/* Function returning the volume fractions of a fan object */
void rotor_coord(struct sRotor* sR, scalar* obj) {

      	scalar sph[], plnu[], plnd[];

    	fraction(sph, -sq((x - sR.x0)) - sq((y - sR.y0)) - sq((z - sR.z0)) + sq(sR.R));
    	fraction(plnu,  sR.nr.x*(x - sR.x0) + sR.nr.y*(y - sR.y0) + sR.nr.z*(z - sR.z0) + sR.W/2.);
    	fraction(plnd, -sR.nr.x*(x - sR.x0) - sR.nr.y*(y - sR.y0) - sR.nr.z*(z - sR.z0) + sR.W/2.);	

	foreach () {
    		obj[] = sph[] * plnu[] * plnd[];
   	}
	boundary({obj});

	void
}

 void rotor_forcing(scalar* obj, vector* uu, struct sDiag* dia, struct sRotor* sR){
	double tempW = 0.;
	double w, wsgn, damp, usgn, utemp;

	foreach(reduction(+:tempW)) {		
		if(obj[] > 0.) {
			foreach_dimension() {

			// Work in respective direction 
			wsgn = sign(sR.nf.x*uu.x[]) + (sign(sR.nf.x*uu.x[]) == 0)*sign(sR.nf.x);
			damp = sR.rampT > t ? t/sR.rampT : 1.;
			w = wsgn*obj[]*damp*sq(sR.nf.x)*(2./rho[])*(sR.P/sR.V)*dt;
			tempW += 0.5*w*sq(Delta);

			// New kinetic energy
			utemp = sq(uu.x[]) + w;

			usgn = 1.*(uu.x[] >= 0)*(utemp > 0) +
			    	     -1.*(uu.x[] >= 0)*(utemp < 0) +
		 		      1.*(uu.x[] <  0)*(utemp < 0) +
				     -1.*(uu.x[] <  0)*(utemp > 0); 

			uu.x[] = usgn*sqrt(fabs(utemp));
		}
		}
	}
	
	dia.Wdone += tempW;

}


