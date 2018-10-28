/**
# Output functions

## *output_field()*: Multiple fields interpolated on a regular grid (text format)

This function interpolates a *list* of fields on a *n+1 x n+1* regular
grid. The resulting data are written in text format in the file
pointed to by *fp*. The correspondance between column numbers and
variables is summarised in the first line of the file. The data are
written row-by-row and each row is separated from the next by a blank
line. This format is compatible with the *splot* command of *gnuplot*
i.e. one could use something like

~~~bash
gnuplot> set pm3d map
gnuplot> splot 'fields' u 1:2:4
~~~

The arguments and their default values are:

*list*
: list of fields to output. Default is *all*.

*fp*
: file pointer. Default is *stdout*.

*n*
: number of points along each dimension. Default is *N*.

*linear*
: use first-order (default) or bilinear interpolation. */

struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};

trace
void output_field (struct OutputField p)
{
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = stdout;
  p.n++;
  
  int len = list_len(p.list);
  double ** field = (double **) matrix_new (p.n, p.n, len*sizeof(double));
  
  double Delta = 0.999999*L0/(p.n - 1);
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + X0;
    for (int j = 0; j < p.n; j++) {
      double y = Delta*j + Y0;
      if (p.linear) {
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = interpolate (s, x, y);
      }
      else {
	Point point = locate (x, y);
	int k = 0;
	for (scalar s in p.list)
	  field[i][len*j + k++] = point.level >= 0 ? s[] : nodata;
      }
    }
  }

  if (pid() == 0) { // master
@if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    for (scalar s in p.list)
      fprintf (p.fp, " %d:%s", i++, s.name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + X0;
      for (int j = 0; j < p.n; j++) {
	double y = Delta*j + Y0;
	//	map (x, y);
	fprintf (p.fp, "%g %g", x, y);
	int k = 0;
	for (scalar s in p.list)
	  fprintf (p.fp, " %g", field[i][len*j + k++]);
	fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
@if _MPI
  else // slave
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
		MPI_COMM_WORLD);
@endif

  matrix_free (field);
}

