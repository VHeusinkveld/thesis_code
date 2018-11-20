#include "grid/octree.h"
#include "run.h"
#include "equi_data.h" 
int minlevel = 5;
int maxlevel = 6;
scalar b[];
scalar lev[];


int main(){    
    DT = 0.1;

    init_grid(1<<minlevel);
    L0 = 100.;

    run();	
}

event init(i=0){
	refine(y <= L0/2. && level < maxlevel);
}

event update(i++) {
    dt = dtnext(DT);
    if(i>2){
    unrefine(y <= L0/2. && level > minlevel-1);
    }

    foreach() {

 	if(i < 5){
    	    b[] = 1.;
        } else {
            b[] = 0.;
        }

        lev[]=level;
    }
    static FILE * fp1 = popen ("ppm2mp4 lev.mp4", "w");
    output_ppm(lev, fp1, n=512, min=minlevel, max=maxlevel);
}
    



event diag(i++){
    fprintf(stderr, "time=%g\n",t);
    equi_diag((scalar *){b});
}

event end(i=9){
    FILE * fp2 = fopen("output", "w");
    equi_output((scalar *){b}, fp2);
    fclose(fp2);
}
