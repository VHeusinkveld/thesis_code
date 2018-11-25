#include "utils.h"
#if dimension == 3
	#include "lambda2.h"
#endif
#include "output_slices.h"
#include "profile5c.h"
#include "equi_data.h"

/** Define variables and structures to do: diagnostics, ouput data, output movies. */

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

struct sEquiDiag {
    int level;			// Level for which diagnostics should be done
    int ii; 			// Keep track of how many additions are done
    double dtDiag;
    double startDiag;
    double endDiag;
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
	double theta;			// Theta for 3D bview movie
	double sphi; 			// Polar angle for sliced image RED
	double stheta;			// Azimuthal angle for sliced image RED

};

/** Initialize structures */
struct sOutput out = {.dtDiag = 1., .dtVisual=2., .dtSlices=1000000., .dtProfile=30., .main_dir="results", .sim_i=0};

struct sEquiDiag ediag = {.level = 6, .ii = 0, .startDiag = 60., .dtDiag = 2.};

struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

event init(i = 0){
	bvsets.phi = 0.;
	bvsets.theta = -M_PI/6.;
	bvsets.sphi = 0.;
	bvsets.stheta = 0.;
}

/** Diagnosing: kinetic energy, diagnosed rotor volume, buoyancy energy, ammount of cells used.*/
event diagnostics (t+=out.dtDiag){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy field 
	double tempVol = 0.;    // Temp volume 
	double tempEkin = 0.;   // Temp kinetic energy
	double tempDiss = 0.;   // Temp dissipation 
	double maxVel = 0.;     // Maximum velocity in fan
	double bEnergy = 0.;    // Buoyant energy
		
	/** Loop over cells to get diagnostics */ 
	foreach(reduction(+:n) reduction(+:tempVol) reduction(+:tempEkin) 
		reduction(max:maxVel) reduction(+:bEnergy) reduction(+:tempDiss)) {

		tempVol += dv()*fan[];
		bEnergy += dv()*y*(b[] - strat(y));
		
		foreach_dimension() {
			ekin[] += sq(u.x[]);
		}
		maxVel = max(maxVel, sq(ekin[]));
		ekin[] *= 0.5*rho[]*dv();	
		tempEkin += ekin[];
		n++;
	}
	/** Assign values to respective global sturcture vars */ 
	dia.diss = 1.*tempDiss;
	dia.bE = 1.*bEnergy;
	rot.diaVol = dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	/** Check if fan volume is within one percent of definition */
	if(fabs(dia.rotVol/rot.V - 1) > 0.20){
		fprintf(stderr, "ERROR Check fan volume, V=%g, Vr=%g\n",rot.V, dia.rotVol);
	}
	
	if (pid() == 0){
	/** Write away simulation data and case setup for main thread */ 
	char nameOut[90];
	char nameCase[90];
     	snprintf(nameOut, 90, "./%s/output", out.dir);
     	snprintf(nameCase, 90, "./%s/case", out.dir);
	static FILE * fpout = fopen(nameOut, "w");
	static FILE * fpca = fopen(nameCase, "w");

	if(t==0.){
		fprintf(fpca,"L0\tinversion\tTref\txr\tyr\tzr\ttheta\tphi\tr\tW\tP\tcu\trampT\tmaxlvl\tminlvl\teps\n");
		fprintf(fpca, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\n", 
				L0, INVERSION, TREF, rot.x0, rot.y0, rot.z0, rot.theta, rot.phi, rot.R, rot.W, rot.P, rot.cu, rot.rampT, maxlevel, minlevel, eps);
		
	        fprintf(stderr,"n\tred\tEkin\tWork\tbE\n");
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE);
	
	fprintf(stderr, "%d\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*dimension))/n),dia.Ekin,rot.Work,dia.bE);
	
	fflush(fpout);
	fflush(fpca);	

	}

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
	dia.bEold = 1.*dia.bE;
}

#if dimension == 3
/** Profiles in height */
event profiles(t += out.dtProfile) {
	char nameProf[90];
	snprintf(nameProf, 90, "./%s/t=%05g", out.dir_profiles, t);
	field_profile((scalar *){b}, nameProf, n=128);
}

/** Average in time */
event equidiags(t = ediag.startDiag; t += ediag.dtDiag) {
	ediag.ii = equi_diag(b, ediag.level, ediag.ii);
	fprintf(stderr, "%d\n", ediag.ii);
}
event end(TEND){
    char nameEquif[90];
    snprintf(nameEquif, 90, "%s/%s", out.dir, "equifield");
    static FILE * fped = fopen(nameEquif, "w");
    equi_output(b, fped, ediag.level, ediag.ii);
    fclose(fped);
}

#endif

/** Ouputting movies in 2- or 3-D*/

#if dimension == 2
event movies(t += 0.5) {
    vertex scalar omega[]; 	// Vorticity
    scalar lev[];	 	// Grid depth
    scalar ekinRho[]; 		// Kinetic energy

    foreach() {
        omega[] = ((u.y[1,0] - u.y[-1,0]) - (u.x[0,1] - u.x[0,-1]))/(2*Delta); // Curl(u) 
        ekinRho[] = 0.5*rho[]*(sq(u.x[]) + sq(u.y[]));
        lev[] = level;
    }

    boundary ({b, lev, omega, ekinRho});
    output_ppm (b, file = "ppm2mp4 ./results/buoyancy.mp4", n = 1<<maxlevel, linear = true, max=strat(L0), min=strat(0.));
    output_ppm (ekinRho, file = "ppm2mp4 ./results/ekin.mp4", n = 1<<maxlevel, min = 0, max = 0.5*sq(rot.cu));
    output_ppm (omega, file = "ppm2mp4 ./results/vort.mp4", n = 1<<maxlevel, linear = true); 
    output_ppm (lev, file = "ppm2mp4 ./results/grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
}
#elif dimension == 3
event movies(t += out.dtVisual) {
    scalar l2[];

    lambda2(u,l2);
    boundary({l2});
    
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

/** Take relevant field slices and write away */
event slices(t+=out.dtSlices) {
    char nameSlice[90];
    coord slice = {1., 0., 1.};

    for(double yTemp = rot.y0-30.; yTemp<=rot.y0+20.; yTemp+=5.) {
	slice.y = yTemp/L0;
   
    	snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
    	FILE * fpsli = fopen(nameSlice, "w");
    	output_slice(list = (scalar *){b}, fp = fpsli, n = 64, linear = true, plane=slice);
    	fclose(fpsli);
    }

    slice.y = 1.;
    slice.z = 1.;

    for(double xTemp = rot.x0-10.; xTemp<=rot.x0+25.; xTemp+=5.) {
	slice.x = xTemp/L0;
	
	snprintf(nameSlice, 90, "%st=%05gx=%03g", out.dir_slices, t, xTemp);
	FILE * fpsli = fopen(nameSlice, "w");
	output_slice(list = (scalar *){b}, fp = fpsli, n = 64, linear = true, plane=slice);
	fclose(fpsli);
	}

    slice.x = 1.;
    slice.y = 1.;
    slice.z = rot.z0/L0;
    	
    snprintf(nameSlice, 90, "%st=%05gz=%03g", out.dir_slices, t, rot.z0);
    FILE * fpsli = fopen(nameSlice, "w");
    output_slice(list = (scalar *){b}, fp = fpsli, n = 64, linear = true, plane=slice);
    fclose(fpsli);

}
#endif

/** Usefull functions */ 

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


