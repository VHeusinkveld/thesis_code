#include "grid/octree.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include "equi_data.h" 


scalar b[];
scalar * tracers = {b};

int main(){
    init_grid(8);
    L0 = 100.;

    run();	
}

event init(t=0) {
    foreach(){
	u.x[] = 0.;
	u.y[] = 1;
        u.z[] = 0.;
 	b[] = y;
    }
}

event diag(t += 0.5){
    equi_diag(&b, 8, 0.5, 0, t+0.1);
}

event end(t=50){
    equi_output(8, equifield);
}
