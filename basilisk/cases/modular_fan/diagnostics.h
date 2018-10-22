#include "utils.h"
#include "lambda2.h"
#include "profile5b.h"

double vphi, stheta;
scalar T_ave[];
struct sDiag dia; 			// Diagnostics

struct sDiag {
	double Ekin;			// Total kinetic energy
	double EkinOld;			// Track changes in kin energy 
	double WdoneOld;		// Track changes in work done 
	double rotVol;			// Track real rotor volume
};

struct sOutput {
	double dtVisual;
	double dtProfiles;
};

struct sOutput out = {.dtVisual=0.1, .dtProfiles=0.1};

event init(i = 0){
	vphi = -M_PI/6.;
	stheta = 0.;
}

/* Profiler 
event profiler(t += 1.) {
	char name[0x100];
	snprintf(name, sizeof(name), "./output/bout%05d", i);
	profile({b}, name);
} */

/** Some diagnostic function*/
event diagnostics (t+=0.2){
	int n = 0.;
	scalar ekin[]; 		// Kinetic energy
	double tempVol = 0.;
	double tempEkin = 0.;	
	double max_vel = 0.;
	
	foreach(reduction(+:n) reduction(+:tempVol) 
		reduction(+:tempEkin) reduction(max:max_vel)) {
		
		#if dimension == 2			
			ekin[] = 0.5*rho[]*cube(Delta)*(sq(u.x[]) + sq(u.y[])));
			tempVol += sq(Delta)*fan[];

		#elif dimension == 3
			ekin[] = 0.5*rho[]*cube(Delta)*(sq(u.x[]) + sq(u.y[]) + sq(u.z[]));
			tempVol += cube(Delta)*fan[];
		#endif
		tempEkin += ekin[];
		max_vel = max(max_vel, sqrt(sq(u.x[]) + sq(u.y[]) + sq(u.z[])));
		n++;
	}	

	dia.rotVol = 1.*tempVol;
	dia.Ekin = 1.*tempEkin;
	
	if(fabs(dia.rotVol/rot.V - 1) > 0.05){
		printf("ERROR Check fan volume, V=%g, Vr=%g\n",rot.V, dia.rotVol);
	}
	
	static FILE * fpout = fopen("./results/output","w");
	static FILE * fpca = fopen("./results/case","w");

	if(t==0.){
		fprintf(fpca,"L0\n %g\n",L0);	
	
		fprintf(fpca,"r\tW\tP\tmaxlvl\tminlvl\teps\n%g\t%g\t%g\t%d\t%d\t%g\n", 
				rot.R, rot.W, rot.P, maxlevel, minlevel, eps);
		
		fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\n");
	}
	fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\n",
		i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work);
	printf("%d\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work);

	dia.EkinOld = 1.*dia.Ekin;
	dia.WdoneOld = 1.*rot.Work;
}

/** Ouputting some movies, either in 2d or 3d.*/

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
    output_ppm (b, file = "./results/buoyancy.mp4", n = 1<<maxlevel, linear = true);
    output_ppm (ekinRho, file = "./results/ekin.mp4", n = 1<<maxlevel, min = 0, max = 0.5*sq(rot.cu));
    output_ppm (omega, file = "./results/vort.mp4", n = 1<<maxlevel, linear = true); 
    output_ppm (lev, file = "./results/grid_depth.mp4", n = 1<<maxlevel, min = minlevel, max = maxlevel);
}
#elif dimension == 3
event movies(t += out.dtVisual) {
    scalar bfy[];
    scalar l2[];
    lambda2(u,l2);
    scalar vxz[];
   // if(t>30.){
   // foreach(){
   //     T_ave[] = (((t-30.)/0.1)*T_ave[] + b[])/(1. + (t-30.)/0.1);
   // }
    }
    boundary({bfy, l2, vxz});
    
    if(vphi < -M_PI/12.){
        vphi += M_PI/(12*70);
    }
    
    clear();
    view(fov=15, tx = 0., ty = 0.2);
    // translate(-rot.x0/L0, -rot.y0/L0, -rot.z0/L0)
    view(theta= -stheta+M_PI/2., phi = 0., width = 1000, height = 1000);
    translate(-rot.x0,-rot.y0,-rot.z0) {
    //box(notics=false);
    cells(alpha = rot.z0);
    //squares("b", n = {1.,0,0.}, alpha=rot.x0);
    //squares("b", n = {0.,0,1.}, alpha=rot.z0);
    //isosurface("l2", color="bfy");
	
    double rn[3] = {cos(stheta), 0, sin(stheta)};
    double fac = 0.*L0;
    double alp = rn[0]*(rot.x0 + fac*rn[0]) + rn[1]*(rot.y0 + fac*rn[1]) + rn[2]*(rot.z0 + fac*rn[2]);
    draw_vof("fan", fc = {1,0,0});
    squares("b", n = {rn[0],rn[1],rn[2]}, alpha = alp);
    }
    save("./results/visual_3d.mp4");
    stheta += 0.01*M_PI;

/*
    if(t > 25.){
        clear();
	coord nslice;
	double sliceTheta;
	double slice Phi;
	
	nclice.x =;
	nslice.y =; 
	nslice.z =;
        view(relative=false,theta=0, phi=0., tx=0., ty=-slice/L0-0.1);
        box(notics=false);
        squares("T_ave", n = {0.,1.,0.}, alpha=slice, min=0., max=9.81*50./273.);
        save("results/temp_slab.mp4");
    }*/
	
}
#endif
