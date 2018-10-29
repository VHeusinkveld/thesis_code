struct sProf {
  scalar * list; // list of scalar field. The default is `all` 
  char * fname;  // Optional file name 
  double ym;     // lower y coordinate  default is Y0
  double h;      // upper y coordinate. Default is Y0+L0
  double rf;     // reduction factor of query heights. Default is 1
  FILE * fp;     // File pointer, if `fname` is not provided. The default is `stdout`
};

void field_profile(struct sProf p){
  /**
     Default values are set in case they are not provided by the user.
  */
  if(!p.list)
    p.list=all;
  if (!p.ym)
    p.ym=Y0;
  if (!p.h)
    p.h=Y0+L0;
  if (!p.rf)
    p.rf=1;
  if (!p.fname && !p.fp)
    p.fp=stdout;
  
  int len = list_len(p.list);
  boundary(p.list);
  FILE * fp = p.fp;
  char * file = p.fname;
/**
Our favorite worker is tasked with the file writing.
*/
  if (pid()==0){ 
    if (file && (fp = fopen (file, "w")) == NULL) {
      perror (file);
      exit (1);
    }
    assert (fp);
    /**
       For reference, a header is printed.
    */
    fprintf(fp,"y\t");
    for(scalar s in p.list)
      fprintf(fp,"%s\t",s.name);
    fprintf(fp,"\n");
  }
  
  double yp = p.ym;
  double aver[len];

  while (yp<=p.h){
    for (int i=0;i<len;i++)
      aver[i]=0.;
    int m=0;
    double a=0;
    foreach(reduction(+:a) reduction(+:m)){
      if ((fabs(y-yp)<=(Delta/2))){
	m++;
	int k = 0;
#if dimension==2
	a+=Delta;
	for (scalar s in p.list)
	  aver[k++] += point.level >= 0 ? interpolate (s, x, yp,z)*Delta  : 0;
#elif dimension==3
	a+=sq(Delta);
	for (scalar s in p.list)
	  aver[k++] += point.level >= 0 ? interpolate (s, x, yp,z)*sq(Delta)  : 0;
#endif
      }
    }
    // MPI reduction for favorite worker
#if _MPI
    if (pid() == 0){
      MPI_Reduce (MPI_IN_PLACE, &aver[0], len, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
#endif
      int k = 0;
      for (scalar s in p.list){
	if (k ==0)
	  fprintf(fp,"%g\t%g",yp,aver[k]/a);
	if (k > 0)
	  fprintf(fp,"\t%g",aver[k]/a);
	k++;
	if (k == len)
	  fprintf(fp,"\n");
      }
      
      // MPI reduction for slaves
#if _MPI
    }else{ 
      MPI_Reduce (&aver[0], NULL, len, MPI_DOUBLE, MPI_SUM, 0,
		  MPI_COMM_WORLD);
    }
#endif
#if dimension == 2
    yp += p.rf*a/(double)m;
#elif dimension == 3
    yp += p.rf*sqrt(a/(double)m);
#endif
  }
  if (pid()==0)
    fflush(fp);
}

