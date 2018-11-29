/** NOTE: Very much in progress.

These functions are geared towards equidistant diagnostics. The level for which diagnostics are done can be specified. Values will be linearly interpolated if necesarry. 

So far this function works with MPI and only works in 3D. 


Short description of what they do:

//void equi_diag(scalar s){
    adds values to the equifield array
}

//void equi_output(scalar s, FILE * fp){
    outputs equifield to file fp in the format: 
    place = x*n*n + y*n + z*n
}


TODO implement 2D
TODO walking average since values in equifield get big, or divide before summing?
*/

double *equifield = NULL; 	// The diagnostics field

int equi_diag (scalar s, int lvl, int diagi){  
    int len = 1;		// TODO implement multiple scalar field averaging
    int n = 1<<lvl;
    if(!equifield) {
        equifield = calloc(len*n*n*n, sizeof(double)); // Init with zeros
    } 	

    double dDelta = 0.9999999*L0/n;  

    restriction({s}); 		// Make sure that coarse level values exist
    boundary({s}); 		// Update boundaries
    
    // Iterate over diag level of leaf cell
    foreach_level_or_leaf(lvl) {
  	if(level < lvl){
	    int pn = 1<<(lvl-level);
	    for(int i = 0; i < pn; i++){
		for(int j = 0; j < pn; j++){
		    for(int k = 0; k < pn; k++){
			
			double ctrans = dDelta/2. - Delta/2.;	// Translation to corner of cell
			double xx = x + ctrans + i*dDelta;   	// Coordinates in equidistant grid
			double yy = y + ctrans + j*dDelta;
 			double zz = z + ctrans + k*dDelta;

		        int ii = round((xx - dDelta/2.)/dDelta); // Iteration in equidistant grid
	         	int jj = round((yy - dDelta/2.)/dDelta);
            	  	int kk = round((zz - dDelta/2.)/dDelta);

	    		int place = ii*n*n + jj*n + kk*len;	// Location in equidistant grid 
			double temp = interpolate(s, xx, yy, zz); // Linear interpolation 

 			equifield[place] += temp;
		    }
		}
	    }

        } else {
 
            int ii = round((x - dDelta/2.)/dDelta); // x adaptive is x equidistant 
       	    int jj = round((y - dDelta/2.)/dDelta); // Derive place in equidistant grid
     	    int kk = round((z - dDelta/2.)/dDelta);

	    int place = ii*n*n + jj*n + kk*len;
	    equifield[place] += s[];
	}
    }
    return diagi + 1;
}

void equi_output (scalar s, FILE * fp, int lvl, int diagi){
    int len = 1;
    int n = 1<<lvl;

    if(pid() == 0){
// Make sure all thread values are summed up using MPI
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
    int len = 1;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
		int place = i*n*n + j*n + k*len;
		double temp = (equifield[place]/((double)diagi)); // Divide by number of additions
	        fprintf(fp, "%g\t", temp);
	    }
	}    
    }
    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
@endif
}


/**
Started on a function that would import data, so far did not make further progress since implementation in python is considered*/

/*
void equi_load (struct sEquiLoad p){
    int len = 1;
    int n = 1<<diaglevel;
    double dDelta = 0.9999999*L0/n;  

@if _MPI
    MPI_Allreduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
@endif

    foreach(){
	s[] = 0.;
    }
    restriction(p.list);
    boundary(p.list);

    foreach_level(diaglevel) {
	int l =0;
	for(scalar s in p.list) {
       	    int i = round((x - dDelta/2.)/dDelta);
       	    int j = round((y - dDelta/2.)/dDelta);
       	    int k = round((z - dDelta/2.)/dDelta);
    
	    int place = i*n*n + j*n + k*len + l;
	    double temp = equifield[place]/diagii;
            fprintf(stderr, "%g", temp);
	    s[] = temp;
	    l++;
        }
    }
    restriction(p.list);
    boundary(p.list);
    foreach_level(diaglevel+1){
	for(scalar s in p.list) {
	    fprintf(stderr, "%g", s[]);
	}
    } 

}
*/



