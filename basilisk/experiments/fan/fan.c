#include "grid/octree.h" // For 3D
#include "lambda2.h"
#include "view.h"
#include "utils.h" 
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "diffusion.h"
#include "fractions.h"
#include "profile5b.h"

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
struct sKar kar;
scalar fan[];			// Fan volume fraction

double vphi;

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

struct sKar {
	double sL;
	double sV;
	double sT;	
};

/* Functions */
void rotor_init(); 
void rotor_update();
void rotor_coord();
void rotor_forcing();

/*
============================================================================
Main Code, Events
============================================================================
*/

/* Initialisation */
int main() {
	
    	// Grid variables 
    	init_grid(2<<5);
   	L0 = 50.;
   	X0 = Y0 = Z0 = 0.;
	
	//# Such that momentum is better conserved 
	u.x.refine = refine_linear; 
	u.y.refine = refine_linear;
	#if dimension > 2
		u.z.refine = refine_linear;
	#endif

   	// Initialize physics 
   	rotor_init(); 
	const face vector muc[] = {0.*1./3000., 0.*1./3000.};
	mu = muc;
	a = av; // Link acceleration
	kar.sV = pow(4*rot.Prho*rot.W,1./3.);

	// Tell basilisk it is a volume field
	fan.prolongation = fraction_refine;
	p.refine = p.prolongation = refine_linear;
	b.gradient = minmod2; // Flux limiter 

  	// Adaptivity
  	minlevel = 2; 
  	maxlevel = 7;
  	eps = 0.05;

	// Set boundary conditions
	periodic (left);
    	
	u.t[bottom] = dirichlet(0.);
	b[bottom] = dirichlet(0.);
	b[top] = neumann(1.);

	// Limit maximum time step 
	DT = 0.05;
	CFL = 0.5;

	// Run the simulation
    	run();
}

/* Initialisation */
event init(t=0){
	foreach() {
		b[] = 9.81*(y)/273. + 0.001*noise();
	}
	rotor_coord();
	refine (fan[] > 0. && level < maxlevel);
	view(tx = 0., ty = -0.6);
	vphi = -M_PI/6.;
	view(theta=M_PI/4., phi=vphi);
}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
} 

mgstats mgb;

/* Diffusion */
event tracer_diffusion(i++){
	mgb = diffusion(b, dt, mu);
}

/* Forcing by the rotor */
event forcing(i = 1; i++) {
	rotor_coord();
	rotor_forcing();
}

/* Rotate the rotor */
event rotate(t = rot.rampT) {
	// Change center  
	rot.x0 += 0;
	rot.y0 += 0;
	rot.z0 += 0;

	// Change angles 
	rot.theta += 0;
	rot.phi += 0;

	rotor_update();
}

/* Adaptivity function called */
event adapt(i++) {
	adapt_wavelet((scalar *){fan, u},(double []){0.,eps,eps,eps},maxlevel,minlevel);
}

/* Profiler 
event profiler(t += 1.) {
	char name[0x100];
	snprintf(name, sizeof(name), "./output/bout%05d", i);
	profile({b}, name);
} */

/* Visualisation */ 
event movies(t += 0.1) {
	#if dimension < 3
		vertex scalar omega[]; 	// Vorticity
		scalar lev[];	 	// Grid depth
		scalar ekinRho[]; 		// Kinetic energy
		vector db[];
		scalar m[];
		scalar bfy[];
		scalar bfx[];
  
		boundary({b});
		gradients ({b}, {db});
		foreach() {
			m[] = 0;
			foreach_dimension() {
				m[] += sq(db.x[]);
			}
			if (m[] > 0) {
				m[] = log(sqrt(m[])+1.);
			}
			bfy[] = b[]*u.y[];
			bfx[] = b[]*u.x[];
			omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
			ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
			lev[] = level;
		}

		boundary ({m, bfx, bfy, lev, omega, ekinRho});
		output_ppm (m, file = "mfield.mp4", n = 1<<maxlevel, min = 0, max = 1, linear = true);
		output_ppm (bfy, file = "bfluxy.mp4", n = 1<<maxlevel, linear = true);
		output_ppm (bfx, file = "bfluxx.mp4", n = 1<<maxlevel, linear = true);
		output_ppm (b, file = "buoyancy.mp4", n = 1<<maxlevel, linear = true);
		output_ppm (fan, file = "coord_fan.mp4", n = 1<<maxlevel, max = 1, min = 0);
		output_ppm (ekinRho, file = "ekin.mp4", n = 1<<maxlevel, min = 0, max = 0.5*sq(kar.sV));
		output_ppm (omega, file = "vort.mp4", n = 1<<maxlevel, linear = true); 
		output_ppm (lev, file = "grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
	#elif dimension > 2
		scalar bfy[];
		scalar l2[];
		lambda2(u,l2);
	
		foreach(){
			bfy[] = b[]*u.y[];
		}

		boundary({bfy, l2, fan});
		
		if(vphi < M_PI/36){
			vphi += M_PI/(36/7*70);
		}

		clear();
		// translate(-rot.x0/L0, -rot.y0/L0, -rot.z0/L0)
		view(theta= M_PI/4., phi = vphi);
		box(notics=false);
		cells(alpha = rot.z0);
		if(t > 100.){
			squares("bfy", n = {0.,1,0.}, alpha=L0/2.);
		}
		isosurface("l2", color = "bfy");
		draw_vof("fan", fc = {1,0,0});
		save("visual_3d.mp4");
	#endif
}

/* Sanity checks */
event sanity (t += 1){
	
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;	
	double bf;
	
	foreach(reduction(+:tempVol) reduction(+:tempEkin) reduction(+:bf)) {
		ekin[] = 0.5*rho[]*sq(Delta)*(sq(u.x[]) + sq(u.y[]));
		tempEkin += ekin[];
		#if dimension > 1			
			tempVol += sq(Delta)*fan[];
			if (y<rot.y0/2.) {
				bf += u.y[]*b[]*Delta;
			}
		#elif dimension > 2
			tempVol = cube(Delta)*fan[];
			bf += u.y[]*b[]*sq(Delta);
		#endif
	}


	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	printf("bf=%g\n", bf);

	/*
	printf("V=%g, Vr=%g, ",rot.V, dia.rotVol);
	printf("Energy: Ek=%g, W=%g, Ek/W=%g, dEk/dW=%g\n", 
		dia.Ekin, dia.Wdone, dia.Ekin/dia.Wdone, 
		(dia.Ekin-dia.EkinOld)/(dia.Wdone-dia.WdoneOld));
	*/

	
	
	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*dia.Wdone;
}

/* Progress event */
event end(t+=2.; t <= 20.) {
	printf("i=%d t=%g p=%d u=%d b=%d \n", i, t, mgp.i, mgu.i, mgb.i);
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
	rot.R = L0/30.;     
	rot.W = rot.R/4.;                      
    	rot.Prho = 50.;
    
   	rot.x0 = L0/2.;
	rot.y0 = 3.*L0/4.;
	#if dimension < 3
		rot.z0 = 0.;
	#elif dimension > 2
		rot.z0 = L0/2.;
	#endif
	rot.theta = M_PI/2.;	// Polar angle
	rot.phi = -M_PI/2.;	// Azimuthal angle 

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


