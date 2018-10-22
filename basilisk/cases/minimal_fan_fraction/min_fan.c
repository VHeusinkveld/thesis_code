/** We investigate the volume error caused by fraction.h of a sliced sphere (a 'fan').*/
#include "grid/octree.h"
#include "utils.h"
#include "fractions.h"

/** Defining a starting and an ending level, and the global scalar field fan[], since this field is iteratively used to refine where the fan is. */
int startlevel = 1;
int endlevel = 7;
int ilevel;
scalar fan[];

/** For the volume of a sliced sphere one can visit [wikipedia](https://en.wikipedia.org/wiki/Spherical_cap) */

int main() {
	init_grid(2<<startlevel);
	for(ilevel = startlevel; ilevel<=endlevel;ilevel++){
  
	double VolEst = 0.;	// Estimated volume from fractions.h
	double VolCalc = 0.; 	// Calculated volume 

  	double L0 = 1.;
  	X0 = Y0 = Z0 = 0.;
 
  	double R = L0/4.;
  	double w = R/5.;
	
	VolCalc = 4./3.*M_PI*pow(R, 3.) - 2*M_PI*pow(R-w/2., 2.)/3.*(3*R - (R-w/2.));
  
  	double xf = L0/2.;
  	double yf = L0/2.;
  	double zf = L0/2.;

  	double rTheta = 0.;  
  	double rPhi = 0.;
	
    	coord fn = {sin(rTheta)*cos(-rPhi + M_PI/2.), sin(rTheta)*sin(-rPhi + M_PI/2.), cos(rTheta)};

	scalar sph[], planeup[], planedown[];
	fan.prolongation = fraction_refine;
	refine(fan[]>0.000001 && level < ilevel); // Refined where needed, 0.0001 is arbitrarily chosen since otherwise due to truncation errors it will refine the whole domain.

    	fraction(sph, -sq((x - xf)) - sq((y - yf)) - sq((z - zf)) + sq(R));
    	fraction(planeup, fn.x*(x-xf) + fn.y*(y-yf) + fn.z*(z-zf) + w/2.);
    	fraction(planedown, -fn.x*(x-xf) - fn.y*(y-yf) - fn.z*(z-zf) + w/2.);
    	foreach () {
      		fan[] = sph[]*planeup[]*planedown[];
	}
	int nn = 0.;
	int nc = 0.;
	foreach (reduction(+:VolEst) reduction(+:nn) reduction(+:nc)) {
		VolEst += dv()*fan[];
		nc++;
		if(fan[]>0.){
			nn++;
		}
	}
	// printing restults
 	printf("VolCalc=%g\tVolEst=%g\terr=%.2g%%\tnn=%d\tnc=%d\tlvl=%d\tred=%g \n", VolCalc, VolEst, 100.*(VolCalc-VolEst)/VolCalc, nn, nc, ilevel, (double)((2<<(ilevel*3))/nc));
	
/** Also we would like to see how the error varies with fan radius to delta ratio */ 
	if(ilevel == endlevel){

  	double L0 = 1.;
	double startR = L0/4.;
	double endR = L0/100.;
	
	for(R=startR; R>= endR; R-=L0/100.){

	double VolEst = 0.;	// Estimated volume from fractions.h
	double VolCalc = 0.; 	// Calculated volume 

  	X0 = Y0 = Z0 = 0.;
 
  	double w = R/5.;
	
	VolCalc = 4./3.*M_PI*pow(R, 3.) - 2*M_PI*pow(R-w/2., 2.)/3.*(3*R - (R-w/2.));
  
  	double xf = L0/2.;
  	double yf = L0/2.;
  	double zf = L0/2.;

  	double rTheta = 0.;  
  	double rPhi = 0.;
	
    	coord fn = {sin(rTheta)*cos(-rPhi + M_PI/2.), sin(rTheta)*sin(-rPhi + M_PI/2.), cos(rTheta)};

	scalar sph[], planeup[], planedown[];
	fan.prolongation = fraction_refine;

    	fraction(sph, -sq((x - xf)) - sq((y - yf)) - sq((z - zf)) + sq(R));
    	fraction(planeup, fn.x*(x-xf) + fn.y*(y-yf) + fn.z*(z-zf) + w/2.);
    	fraction(planedown, -fn.x*(x-xf) - fn.y*(y-yf) - fn.z*(z-zf) + w/2.);
    	foreach () {
      		fan[] = sph[]*planeup[]*planedown[];
	}
	int nn = 0.;
	int nc = 0.;
	foreach (reduction(+:VolEst) reduction(+:nn) reduction(+:nc)) {
		VolEst += dv()*fan[];
		nc++;
		if(fan[]>0.){
			nn++;
		}
	}
	// printing restults
 	printf("VolCalc=%g\tVolEst=%g\terr=%.2g%%\tnn=%d\tlvl=%d\tratio=%g\n", VolCalc, VolEst, 100.*(VolCalc-VolEst)/VolCalc, nn, ilevel, (double)(R/(L0/(2<<ilevel))));
	}
	}
	}	
}

