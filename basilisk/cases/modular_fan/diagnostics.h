#include "utils.h"
#include "lambda2.h"

#include "output_slices.h"
#include "profile5c.h"

/** Define variables and structures to do: diagnostics, ouput data, output movies. */

scalar b_ave[];				// Averaged buoyancy field
scalar Ri[];
scalar bdiff[];
struct sDiag dia; 			// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Diagnosed rotor volum
	double bE;			// Buoyancy energy
	double bEold;			// Track changes in buoyancy energy
	double diss;			// Diagnose dissipation
};

struct sOutput {
	double dtDiag;
	double dtVisual;
	double dtSlices;
	double dtProfile;
	double startAve;
	double dtAve;
	char main_dir[12];
	char dir[30];
	char dir_profiles[60];
	char dir_slices[60];
	int sim_i;
};

struct sbViewSettings {
	double phi; 			// Phi for 3D bview movie
	double theta;
	double sphi;
	double stheta;			// Theta for sliced image

};

/** Initialize structures */
struct sOutput out = {.dtDiag = 0.5, .dtVisual=1., .dtSlices=30., .dtProfile=1., .startAve=0., .dtAve = 30., .main_dir="results", .sim_i=0};
struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

double dissipation(Point point, vector u);

event init(i = 0){
	bvsets.phi = 0.;
	bvsets.theta = -M_PI/6.;
	bvsets.sphi = 0.;
	bvsets.stheta = 0.;
}
/** Profiles in height */
event profiles(t += out.dtProfile) {
	char nameProf[90];
	snprintf(nameProf, 90, "./%s/t=%05g", out.dir_profiles, t);
	field_profile((scalar *){b, Ri, bdiff, u}, nameProf);
}

/** Diagnosing: kinetic energy, diagnosed rotor volume, buoyancy energy, ammount of cells used.*/
event diagnostics (t+=out.dtDiag){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;
	double tempDiss = 0.;
	double maxVel = 0.;
	double turbVol = 0.;
	double bEnergy = 0.;
		

	foreach(reduction(+:n) reduction(+:tempVol) reduction(+:tempEkin) 
		reduction(max:maxVel) reduction(+:bEnergy) reduction(+:turbVol)
        	reduction(+:tempDiss)) {

		tempVol += dv()*fan[];
		bEnergy += dv()*y*(b[] - strat(y));
		
		foreach_dimension() {
			ekin[] += sq(u.x[]);
		}

		tempDiss += dissipation(point, u);

		maxVel = max(maxVel, sq(ekin[]));
		ekin[] *= 0.5*rho[]*dv();	
		tempEkin += ekin[];

		n++;

		Ri[] = ((b[0,1]-b[0,-1])/(2*Delta))/(sq((u.x[0,-1]-uf.x[0,1])/(2*Delta)) + sq((uf.z[0,1]-uf.z[0,-1])/(2*Delta)) + 0.000000000001); 
		turbVol += dv()*(Ri[] < 0.25 ? 1. : 0.);
	}
	dia.diss = 1.*tempDiss;
	dia.bE = 1.*bEnergy;
	rot.diaVol = dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	/** Check if fan volume is within one percent of definition */
	if(fabs(dia.rotVol/rot.V - 1) > 0.05){
		fprintf(stderr, "ERROR Check fan volume, V=%g, Vr=%g\n",rot.V, dia.rotVol);
	}
	
	if (pid() == 0){
	/** Write away data */ 
	char nameOut[90];
	char nameCase[90];
     	snprintf(nameOut, 90, "./%s/output", out.dir);
     	snprintf(nameCase, 90, "./%s/case", out.dir);
	static FILE * fpout = fopen(nameOut, "w");
	static FILE * fpca = fopen(nameCase, "w");

	if(t==0.){
		fprintf(fpca,"L0\tinversion\txr\tyr\tzr\ttheta\tphi\tr\tW\tP\tmaxlvl\tminlvl\teps\n%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\n", 
				L0, INVERSION, rot.x0, rot.y0, rot.z0, rot.theta, rot.phi, rot.R, rot.W, rot.P, maxlevel, minlevel, eps);
		
	        fprintf(stderr,"n\tred\tEkin\tWork\tbE\tturbVol\tdissipation\n");
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\tturbVol\tdissipation\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE, turbVol, dia.diss);
	
	fprintf(stderr, "%d\t%g\t%g\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work,dia.bE,turbVol,dia.diss);
	
	fflush(fpout);
	fflush(fpca);	

	}

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
	dia.bEold = 1.*dia.bE;
}

/** Ouputting movies */

#if dimension == 2
event movies(t += out.dtVisual) {
    vertex scalar omega[]; 	// Vorticity
    scalar lev[];	 	// Grid depth
    scalar ekinRho[]; 		// Kinetic energy

    foreach() {
        omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
        ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
        lev[] = level;
    }

    boundary ({b, lev, omega, ekinRho});
    output_ppm (b, file = "ppm2mp4 ./results/buoyancy.mp4", n = 1<<maxlevel, linear = true);
    output_ppm (ekinRho, file = "ppm2mp4 ./results/ekin.mp4", n = 1<<maxlevel, min = 0, max = 0.5*sq(rot.cu));
    output_ppm (omega, file = "ppm2mp4 ./results/vort.mp4", n = 1<<maxlevel, linear = true); 
    output_ppm (lev, file = "ppm2mp4 ./results/grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
}
#elif dimension == 3
event movies(t += out.dtVisual) {
    scalar l2[];
    scalar bfy[];
	
    foreach(){
	bdiff[] = b[] - strat(y);
	bfy[] = b[]*u.y[];
    }

    lambda2(u,l2);
    boundary({l2, bdiff});
    
    clear();
    view(fov=25, tx = 0., ty = 0., phi=bvsets.phi, theta=bvsets.theta, width = 800, height = 800);
    
    translate(-rot.x0,-rot.y0,-rot.z0) {
        box(notics=false);
        isosurface("l2", v=-0.05, color="b", min=strat(0.), max=strat(L0));
	draw_vof("fan", fc = {1,0,0});
    }
    translate(-rot.z0,-rot.y0, -L0){
      	squares("b", n = {0.,0,1.}, alpha=rot.z0, min=strat(0.), max=strat(L0));
        cells(n = {0.,0.,1.}, alpha = rot.z0);
    }
    
    translate(0.,-rot.y0,-rot.z0){
        squares("b", n = {1.,0,0.}, alpha=rot.x0, min=strat(0.), max=strat(L0));
    }

    /** Save file with certain fps*/
    char nameVid1[90];
    snprintf(nameVid1, 90, "ppm2mp4 -r %g ./%s/visual_3d.mp4", 10., out.dir);
    save(nameVid1);
}

/** Relevant field slicse */
event slices(t+=out.dtSlices) {
    char nameSlice[90];
    coord slice = {1., 0., 1.};

    for(double yTemp = rot.y0-25.; yTemp<=rot.y0+25.; yTemp+=5.) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = {b}, fp = fpsli, n = 99, plane=slice);
    	fclose(fpsli);
    }

    slice.y = 1.;
    slice.z = rot.z0/L0;
    	
    snprintf(nameSlice, 90, "%st=%05gz=%03g", out.dir_slices, t, rot.z0);
    FILE * fpsli = fopen(nameSlice, "w");
    output_slice(list = {b}, fp = fpsli, n = 99, plane=slice);
    fclose(fpsli);

}
#endif

double dissipation(Point point, vector u) {
    
    double dis = 0.;
    foreach_dimension() {
    dis +=  sq(((u.x[1] - u.x[-1])/(2*Delta))) +
 	    sq(((u.x[0,1] - u.x[0,-1])/(2*Delta))) +
	    sq(((u.x[0,0,1] - u.x[0,0,-1])/(2*Delta)));
    }
    dis *= dv()*Evis[];
 
    return dis;
}

/** Checks if required folders exists, if not they get created. */
void sim_dir_create(){
   
    sprintf(out.dir, "./%s/%s%02d", out.main_dir, sim_ID, out.sim_i);
    sprintf(out.dir_profiles, "%s/profiles/", out.dir);
    sprintf(out.dir_slices, "%s/slices/", out.dir);
 
    if (pid() == 0){
    struct stat st = {0};
    if (stat(out.main_dir, &st) == -1) {
        mkdir(out.main_dir, 0777);
    }
   
    if (stat(out.dir, &st) == -1) {
        mkdir(out.dir, 0777);
    }
    if (stat(out.dir_slices, &st) == -1) {
        mkdir(out.dir_slices, 0777);
    }
    if (stat(out.dir_profiles, &st) == -1) {
        mkdir(out.dir_profiles, 0777);
    }        
    }
}

