#include "grid/octree.h"
#include "run.h"
#include "equi_data.h" 

int minlevel = 2;
int maxlevel = 3;
scalar b;



int main(){    
    DT = 0.1;

    init_grid(1<<minlevel);
    L0 = 100.;

    run();	
}

event update(i++) {
    dt = dtnext(DT);

    foreach() {
 	if(i <= 50){
    	    b[] = 1.;
        } else {
            b[] = 0.;
        }
    }
}

event diag(i++){
    fprintf(stderr, "time=%g\n",t);
    equi_diag((scalar *){b});
}

event end(i=100){
    //FILE * fp = fopen("output", "w");
    //equi_output((scalar *){b}, fp, Ndia);
    //fclose(fp);
}
