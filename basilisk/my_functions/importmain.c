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
int diaglvl = 7;

int main() {
    init_grid(1<<7);

    equi_import_binary(h, "equifield", diaglvl);
    //int ll = 0;
    foreach() {
        h[] = h[] - strat(y);
       // if(h[] > 0.01) {
	//    ll++;
        //    printf("%g, %d\n", h[], ll);
	//}
    }
    view(phi=0.3, theta=0.3);
    cells(n = {1,0,0});
    isosurface("h", 0.5);
    isosurface("h", 0.3);
    save("image.png");

}
