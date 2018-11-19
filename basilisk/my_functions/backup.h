struct sEquiDia {
    scalar * list;
    int n;
    double timestep;
    double starttime;
    double t;
};

struct sEquiOut {
    scalar * list;
    FILE * fp;
    int n;
};


double *equifield = NULL;

void equi_diag (struct sEquiDia p){  
    int len = list_len(p.list);
    if(!equifield) {
        equifield = malloc(len*sizeof(double)*p.n*p.n*p.n);
    } 	
    double Delta = 0.9999999*L0/(p.n - 1);    
   
    for(int i = 0; i < p.n; i++){
	double x = i*Delta + X0;
        for(int j = 0; j < p.n; j++){
	    double y = j*Delta + Y0;
            for(int k = 0; k < p.n; k++){
                double z = k*Delta + Z0;
		double w1 = (p.t - p.starttime - p.timestep)/(p.t-p.starttime) ;
	        double w2 = p.timestep/(p.t-p.starttime);
  		int l = 0;
                for (scalar s in p.list){
		    int place = i*p.n*p.n + j*p.n + k*len + l;
		    if (p.t > p.starttime){
		        equifield[place] = w1*equifield[place] + w2*interpolate(s, x, y, z);
		    } else {
		    	equifield[place] = 0.;
		    }
		    l++;
 		}
	    }
	}    
    }
    
    if (pid() == 0) {
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*p.n*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
@endif

    }
@if _MPI
    else // slave
    MPI_Reduce (equifield, NULL, len*p.n*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
@endif
    
}

void equi_output (struct sEquiOut p){
   if(pid() == 0){
   int len = list_len(p.list);
   int ja = 0.;
    for(int i = 0; i < p.n; i++){
        for(int j = 0; j < p.n; j++){
            for(int k = 0; k < p.n; k++){
		int l = 0;
		for (scalar s in p.list){
		    int place = i*p.n*p.n + j*p.n + k*len + l;
		    fprintf(p.fp, "%g\t", equifield[place]);
		    l++;
		    ja++;
   		}
	    }
	}    
    }
    fprintf(stderr, "%d\n", ja);
}
}



