#if dimension == 3
	#include "SGS.h"
#endif 

#define CP 1005.	// C_p for air 
#define gCONST 9.81	// Gravitational constant
#define TREF 273.	// Kelvin
#define INVERSION .1 	// Kelvin per meter

#define strat(s) gCONST*(INVERSION + gCONST/CP)*s/TREF // Stratification 

double crho = 1.;
scalar b[];		
scalar * tracers = {b};
face vector av[]; 
struct sCase def;

struct sCase {
	double ugeo;
	double vgeo;
	double corf;
};

void init_physics(){
 	def.ugeo = 0.;
	def.vgeo = 0.;
	def.corf = pow(10.,-4.);
	b.nodump = true; 
    	
        u.n[bottom] = dirichlet(0.);
	u.t[bottom] = dirichlet(0.);
	u.n[top] = dirichlet(0.);
	u.t[top] = dirichlet(def.ugeo);


	b[bottom] = neumann(0.);
	b[top] = neumann(0.);

	periodic (left);

	#if dimension == 3
		u.r[bottom] = dirichlet(0.);
		u.r[top] = dirichlet(def.vgeo);
		
		Evis[bottom] = dirichlet(0.);
        	Evis[top] = dirichlet(0.);

		periodic(front);
	#endif  

	foreach() {
		b[] = strat(y);
		u.x[] = 1.*def.ugeo;
		#if dimension == 3
			u.z[] = 1.*def.vgeo;
		#endif
	}
}

/* Gravity forcing */
event acceleration(i++){
	foreach_face(y){
		av.y[] = (b[] + b[0,-1])/2.;
	}
	if(def.ugeo>0.){
		foreach_face(x){
			av.x[] = (def.ugeo - uf.x[])*def.corf;
		}
	}
	#if dimension == 3
	if(def.vgeo>0.){
		foreach_face(z){
			av.z[] = (def.vgeo - uf.z[])*def.corf;
		}
	} 
	#endif
}

mgstats mgb;

/* Diffusion */
event tracer_diffusion(i++){
	mgb = diffusion(b, dt, mu);
}
