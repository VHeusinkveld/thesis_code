/**
# An improved version of the profiling function under `profile5b.h`

This profiling function uses the same profiling strategy as the `void
profile()` function presented [here](profil5.h). However, this version
of the same function is subjectively more user friendly. Especially
since it mimics the user interface of other output functions in
Basilisk. Also its behaviour when using MPI and/or 2D grids, is as one
may naively expect. An extra function is deterministic profiles with
at a fixed number of levels (`n`). Furthermore, it addresses the bug
with "`too many open files`", that was found by V. Heusinkveld. 

Main features:  

- Consistent vertical profiles on tree grids, without upsampling.  
- The required computational effort scales approximately with the number of (leaf) grid cells.  
- Parrallel MPI compatible and scaleable.
 
## input

The function takes a structure as argument, this faciliates optional arguments.

You may run *full default* by calling `profile(NULL)` to obtain a
profile of all scalar fields in the standard output (e.g. your
terminal).

At the cost/price of vertical resolution, the *rf* argument can be
used to (linearly) reduce/increase the computional effort this
function requires.
*/

struct prof {
  scalar * list; // list of scalar field. The default is `all` 
  char * fname;  // Optional file name 
  double ym;     // lower y coordinate  default is the first "model level"
  double h;      // upper y coordinate. Default is Y0+L0
  double rf;     // reduction factor of query heights. Default is 1
  FILE * fp;     // File pointer, if `fname` is not provided. The default is `stdout`
  int n;         // Number of levels for equidistant profiles
};

double average_over_yp(scalar * list, double * v, double yp){
  int m = 0;
  double a = 0;
  int len = list_len(list);
  foreach(reduction(+:a) reduction(+:m)){
    if ((fabs(y-yp) <= (Delta/2))){
      m++;
      int k = 0;
#if dimension==2
      a += Delta;
      for (scalar s in list)
	v[k++] += point.level >= 0 ? interpolate (s, x, yp, z)*Delta : 0;
#elif dimension==3
      a += sq(Delta);
      for (scalar s in list)
	v[k++] += point.level >= 0 ? interpolate (s, x, yp, z)*sq(Delta) : 0;
#endif
    }
  }
#if _MPI
  if (pid() != 0){
    MPI_Reduce (&v[0], NULL, len, MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
  }else{
    MPI_Reduce (MPI_IN_PLACE, &v[0], len, MPI_DOUBLE, MPI_SUM, 0,
		MPI_COMM_WORLD);
  }
#endif
  int g = 0;
  for (scalar s in list)
    v[g++] /= a;
#if dimension == 2
  return a/(double)m;
#elif dimension == 3
  return sqrt(a/(double)m);
#endif
}

void field_profile(struct prof p){
  /**
     Default values are set in case they are not provided by the user.
  */
  if(!p.list)
    p.list = all;
  if (!p.ym) 
    p.ym = Y0;
  if (!p.h)
    p.h = Y0+L0;
  if (!p.rf)
    p.rf = 1;
  if (!p.fname && !p.fp)
    p.fp = stdout;
  double dzn;
  if (p.n)
    dzn = (p.h - p.ym) / ((double)p.n - 0.9999999);
  int len = list_len(p.list);
  boundary(p.list);
  FILE * fp = p.fp;
  char * file = p.fname;
  /**
     The favorite worker is tasked with the file writing.
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
  /**
     Here a loop starts that iteratively cycles over different y-coordinates. The vertical-step size is governed by the average grid spacing (at that height) and the reduction factor *rf*, unless `n` is passed by the user. 
  */
  int m = 1;
  while (yp <= p.h){
    double aver[len];
    memset(&aver, 0, sizeof(aver)); 
    double dz = average_over_yp(p.list, aver, yp);
    if (pid() == 0){
      int k = 0;
      for (scalar s in p.list){
	if (k == 0)
	  fprintf(fp,"%g\t%g",yp,aver[k]);
	if (k > 0)
	  fprintf(fp,"\t%g",aver[k]);
	k++;
	if (k == len)
	  fprintf(fp,"\n");
      }
    }
    if (p.n)
      dz = dzn/p.rf;
    yp += p.rf*dz;
  }
  if (pid()==0){
    fflush(fp);
    if (fp != stdout) // Not close stdout
      fclose(fp);
  }
}


