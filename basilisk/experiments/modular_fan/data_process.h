#include "utils.h"
#include "lambda2.h"
#include "profile5b.h"

double vphi;
scalar T_ave[];

struct sOutput {
	double dtVisual;
	double dtProfiles;
};

struct sOutput out = {.dtVisual=0.1, .dtProfiles=0.1};

event init(i = 0){
	vphi = -M_PI/6.;
}

/* Profiler 
event profiler(t += 1.) {
	char name[0x100];
	snprintf(name, sizeof(name), "./output/bout%05d", i);
	profile({b}, name);
} */


#if dimension < 3
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

    foreach(){
	printf("b=%g\n",b[]);
        bfy[] = b[]*u.y[];
        T_ave[] = ((t/0.1)*T_ave[] + b[])/(1. + t/0.1);
    }

    boundary({bfy, l2, vxz});
    
    if(vphi < -M_PI/12.){
        vphi += M_PI/(12*70);
    }
    
    clear();
    view(tx = 0., ty = -0.5);
    // translate(-rot.x0/L0, -rot.y0/L0, -rot.z0/L0)
    view(theta= M_PI/4., phi = vphi);
    box(notics=false);
    cells(alpha = rot.z0);
    squares("b", n = {1.,0,0.}, alpha=rot.x0);
    squares("b", n = {0.,0,1.}, alpha=rot.z0);
    isosurface("l2", color="bfy", min=0., max=2.);
    draw_vof("fan", fc = {1,0,0});
    save("./results/visual_3d.mp4");
/*
    if(t > 25.){
        clear();
        double slice = fmod(fabs(rot.y0+2.5-(t-25.)),rot.y0+2.5);
        view(relative=false,theta=0, phi=0., tx=0., ty=-slice/L0-0.1);
        box(notics=false);
        squares("T_ave", n = {0.,1.,0.}, alpha=slice, min=0., max=9.81*50./273.);
        save("results/temp_slab.mp4");
    }*/	
	
}
#endif
