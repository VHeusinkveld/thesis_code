#include "grid/octree.h"
#include "view.h"
#include "fractions.h"	

#include "equi_data.h"

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .2 	// Kelvin per meter

#define strat(s) gCONST*(INVERSION + gCONST/CP)*s/TREF // Stratification 

scalar h[];
scalar f[];
vertex scalar G[];

scalar vol[];

int diaglvl = 7;

int main() {
    init_grid(1<<7);
    L0 = 200;
    X0 = Y0 = Z0 = 0.;
    equi_import_binary(h, "equifieldt=00300", diaglvl);
    int ll = 0;
    foreach() {
        f[] = 273/9.81*(1.*h[] - strat(y));
        //if(h[] > 0.01) {
	  //  ll++;
            //printf("%g, %d\n", h[], ll);
	//}
    }

    foreach_vertex(){
        G[] = -(f[] + f[-1] + f[0,-1] + f[-1,-1] +
	       f[0,0,-1] + f[-1,0,-1] + f[0,-1,-1] + f[-1,-1,-1])/8. + 0.4; 
	//fprintf(stderr, "%g\n", G[]);
	/*
	if(G[] > 0.) {
	    G[] = 100.;
 	} else {
	    G[] = 0.;
	}*/
    }
    
    fractions(G, vol);
	
    view(fov=30, phi=0.3, theta=0.4);
    double temp = 0;
    for(int ii = 0; ii < 1; ii++){
	temp += 0.4;
        cells(n = {1,0,0});
        squares("f", n = {0,0,1}, alpha=100, min=-1., max=1.);
	draw_vof("vol");	
        //isosurface("f", v=temp, color="f");
        //isosurface("h", 0.3);
        //save("ppm2mp4 -r 5 anim.mp4");
        save("image.png");
	clear();
	fprintf(stderr, "i=%d\n", ii);
    }
}
