#include "grid/quadtree.h"
#include "run.h"

int maxlevel = 3;
int diaglevel = 5;

int main() {
    N = 1<<1;
    L0 = 10;

    DT = 0.1;

    init_grid(N);

    run();
}

event init(i=0) {
   refine(y < L0/2. && level < maxlevel);
}

event update(i++) {
    int ii = 0;
    int jj = 0;
    int kk = 0;

    foreach_level_or_leaf(diaglevel){
        printf("%g, %g, %g\n", x, y, Delta);
	ii++;
    }

    foreach() {
	if(level==diaglevel){
	    jj++;
	}    
	if(level==maxlevel){
	    kk++;
	}
    }

    printf("%d, %d, %d\n", ii, jj, kk);
}


event end(i = 3){}
