struct sEquiDia {
    scalar * s;
    int n;
    double timestep;
    double starttime;
    double t;
};

double *equifield = NULL;

void equi_diag (struct sEquiDia p){  
    fprintf(stderr, "hoi");
    if(!equifield) {
        equifield = malloc(sizeof(double)*p.n*p.n*p.n);
	fprintf(stderr, "jaja");
    } 	
    fprintf(stderr, "hoi");

    double Delta = 0.9999999*L0/(p.n - 1);    

    for(int i = 0; i < p.n; i++){
	double x = i*Delta + X0;
        for(int j = 0; j < p.n; j++){
	    double y = j*Delta + Y0;
            for(int k = 0; k < p.n; k++){
                double z = k*Delta + Z0;

		double w1 = (p.t*p.starttime - p.timestep)/(p.t-p.starttime) ;
	        double w2 = p.timestep/(p.t-p.starttime);
		int place = i*p.n*p.n + j*p.n + k;
		equifield[place] = w1*equifield[place] + w2*interpolate(p.s, x, y, z);
	    }
	}    
    }
/*
    if (pid() == 0) {

@if _MPI
    MPI_Reduce (MPI_INPLACE, equifiedl[0], p.n*p.n*p.n, MPI_FLOAT, MPI_MIN, 0, MPI_COMMWORLD);
@endif

    }
@if _MPI
    else // slave
    MPI_Reduce (equifield[0], NULL, p.n*p.n*p.n, MPI_FLOAT, MPI_MIN, 0, MPI_COMMWORLD);
@endif*/
    
}

void equi_output (int n, double *equifield){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            for(int k = 0; k < n; k++){
		int place = i + j*n + k*n*n;
		fprintf(stderr, "%g", equifield[place]);
	    }
	}    
    }
}



