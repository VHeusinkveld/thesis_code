#include "utils.h"
#include "lambda2.h"
#include "profile5b.h"

/** Define variables and structures to do: diagnostics, ouput data, output movies. */

scalar b_ave[];				// Averaged buoyancy field
scalar Ri[];
scalar bdiff[];
struct sDiag dia; 			// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Diagnosed rotor volume
	double bE;			// Buoyancy energy
	double bEold;			// Track changes in buoyancy energy
	double diss;			// Diagnose dissipation
};

struct sOutput {
	double dtVisual;
	double dtSlices;
	double dtProfile;
	double startBave;
	char main_dir[12];
	char dir[20];
	char dir_profiles[50];
	char dir_slices[50];
	int sim_i;
};

struct sbViewSettings {
	double phi; 			// Phi for 3D bview movie
	double theta;
	double sphi;
	double stheta;			// Theta for sliced image

};


/** Initialize structures */
struct sOutput out = {.dtVisual=0.2, .dtSlices=1., .dtProfile=1., .startBave=20., .main_dir="results", .sim_i=0};
struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

event init(i = 0){
	bvsets.phi = 0.;
	bvsets.theta = -M_PI/6.;
	bvsets.sphi = 0.;
	bvsets.stheta = 0.;	
}

/** Profiles for the buoyancy */
event field_profiles(t += out.dtProfile) {
	char nameProf[80];
	snprintf(nameProf, 80, "./%s/t=%05g", out.dir_profiles, t);
	field_profile(list = {b, Ri, bdiff}, fname = nameProf);
}

/** Diagnosing: kinetic energy, diagnosed rotor volume, buoyancy energy, ammount of cells used.*/
event diagnostics (t+=0.2){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;
	double max_vel = 0.;
	double turbVol = 0.;
	double bEnergy = 0.;
	
	foreach(reduction(+:n) reduction(+:tempVol) reduction(+:tempEkin) 
		reduction(max:max_vel) reduction(+:bEnergy) reduction(+:turbVol)) {
		tempVol += dv()*fan[];
		bEnergy += dv()*y*(b[] - strat(y));
		
		foreach_dimension() {
			ekin[] += sq(u.x[]);
		}
		max_vel = max(max_vel, sq(ekin[]));
		ekin[] *= 0.5*rho[]*dv();	
		tempEkin += ekin[];
		n++;
		Ri[] = ((b[0,1]-b[0,-1])/(2*Delta))/(sq((u.x[0,-1]-uf.x[0,1])/(2*Delta)) + sq((uf.z[0,1]-uf.z[0,-1])/(2*Delta)) + 0.000000000001) < 0.25 ? 1. : 0.; 
		turbVol += dv()*Ri[];
	}
	
	dia.bE = 1.*bEnergy;
	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	/** Check if fan volume is within one percent of definition */
	if(fabs(dia.rotVol/rot.V - 1) > 0.01){
		fprintf(stderr, "ERROR Check fan volume, V=%g, Vr=%g\n",rot.V, dia.rotVol);
	}
	
	/** Write away data */ 
	char nameOut[80];
	char nameCase[80];
     	snprintf(nameOut, 80, "./%s/output", out.dir);
     	snprintf(nameCase, 80, "./%s/case", out.dir);
	static FILE * fpout = fopen(nameOut, "w");
	static FILE * fpca = fopen(nameCase, "w");

	if(t==0.){
		fprintf(fpca,"L0\n %g\n",L0);	
	
		fprintf(fpca,"inversion\tr\tW\tP\tmaxlvl\tminlvl\teps\n%g\t%g\t%g\t%g\t%d\t%d\t%g\n", 
				INVERSION, rot.R, rot.W, rot.P, maxlevel, minlevel, eps);
		
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\tturbVol\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE, turbVol);
	
	fprintf(stderr, "%d\t%g\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work,dia.bE,turbVol);

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
	dia.bEold = 1.*dia.bE;
}

/** Ouputting movies and slices of relevant fields */

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
        if(t > out.startBave){
    	    b_ave[] = (((t-out.startBave)/out.dtVisual)*b_ave[] + b[])/(1. + (t-out.startBave)/out.dtVisual);
        }
    }

    lambda2(u,l2);
    boundary({l2, bdiff});
    
    clear();
    view(fov=25, tx = 0., ty = 0., phi=bvsets.phi, theta=bvsets.theta, width = 1200, height = 1200);
    
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
    char nameVid1[80];
    snprintf(nameVid1, 80, "ppm2mp4 -r %g ./%s/visual_3d.mp4", 1./out.dtVisual, out.dir);
    save(nameVid1);
}

event slices(t+=out.dtSlices) {

    double sx=rot.x0, sy=rot.y0, sz=rot.z0;
    bool XY=false, XZ=false, YZ=false;

    output_field(list = {b}, n = 128);

}

#endif

/** Checks if folders exists, if not they get created. */
void sim_dir_create(){
    if (pid() == 0){
    struct stat st = {0};
   
    sprintf(out.dir, "./%s/%s%02d", out.main_dir, sim_ID, out.sim_i);
    sprintf(out.dir_profiles, "%s/profiles/", out.dir);
    sprintf(out.dir_slices, "%s/slices/", out.dir);
 
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

