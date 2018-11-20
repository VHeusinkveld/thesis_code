int diaglevel = 3;
int diagii = 0;
double *equifield = NULL;

struct sEquiDia {
    scalar * list;
};

struct sEquiOut {
    scalar * list;
    FILE * fp;
};



void equi_diag (struct sEquiDia p){  
    int len = list_len(p.list);
    int n = 1<<diaglevel;

    if(!equifield) {
        equifield = calloc(len*n*n*n, sizeof(double)); // Init with zeros
    } 	

    double dDelta = 0.9999999*L0/n;  

    restriction(p.list); // make sure that coarse level values exist
    boundary(p.list); // update boundaries

    foreach_level_or_leaf(diaglevel) {
   	int l = 0; 
        for(scalar s in p.list){
  	if(level < diaglevel){
	    int pn = 1<<(diaglevel-level);
	    for(int i = 0; i < pn; i++){
		for(int j = 0; j < pn; j++){
		    for(int k = 0; k < pn; k++){
			
			double ctrans = dDelta/2. - Delta/2.; 
			double xx = x + ctrans + i*dDelta;
			double yy = y + ctrans + j*dDelta;
 			double zz = z + ctrans + k*dDelta;

		        int ii = round((xx - dDelta/2.)/dDelta);
	         	int jj = round((yy - dDelta/2.)/dDelta);
            	  	int kk = round((zz - dDelta/2.)/dDelta);

	    		int place = ii*n*n + jj*n + kk*len + l;
    	   		//fprintf(stderr, "%g, %g, %g, %g, %g, %g, %g, %g\n", x, y, z, xx, yy, zz, dDelta, Delta);	
			double temp = interpolate(s, xx, yy, zz);
 			equifield[place] += temp;
			fprintf(stderr, "%g\n", temp);	 
		    }
		}
	    }

        } else {
 
            int ii = round((x - dDelta/2.)/dDelta);
       	    int jj = round((y - dDelta/2.)/dDelta);
     	    int kk = round((z - dDelta/2.)/dDelta);

	    int place = ii*n*n + jj*n + kk*len + l;
	    fprintf(stderr, "%d, %g, %g, %g, %g\n", level, x, y, z, s[]);
	    equifield[place] += s[];
	}
        }
	l++;
    }
    diagii++;
    //fprintf(stderr, "%d, %d, %d, %d\n", aa, bb, cc, (int)sq(1<<diaglevel));

}
/*    
    if (pid() == 0) {
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*p.n*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
@endif

    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*p.n*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
@endif

*/
void equi_output (struct sEquiOut p){

    int n = 1<<diaglevel;
	
    if(pid() == 0){

@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif

        int len = list_len(p.list);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
		int l = 0;
		for (scalar s in p.list){
		    int place = i*n*n + j*n + k*len + l;
		    fprintf(p.fp, "%g\t", equifield[place]/diagii);
		    l++;
   		}
	    }
	}    
    }
    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*p.n*p.n*p.n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
}

