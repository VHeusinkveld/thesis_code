#line 1 "main-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "main-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if 1
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif


#line 1 "/home/vinlinux/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if 201307
# include <omp.h>
#elif _MPI
# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe
#endif

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if _MPI
static FILE * qstderr (void);
static FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 72

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if 1
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  assert (d != NULL);
  d->id = pmfunc_index(func, file, line);
  d->size = size;
  pmfunc * f = &pmfuncs[d->id - 1];
  f->total += size;
  if (f->total > f->max)
    f->max = f->total;
  pmtrace.total += size;
  pmtrace.overhead += sizeof(pmdata);
  if (pmtrace.total > pmtrace.max) {
    pmtrace.max = pmtrace.total;
    pmtrace.maxoverhead = pmtrace.overhead;
  }
  pmfunc_trace (f, c);
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", qstderr());
    if (d->size == 0)
      fputs (", possible double free()", qstderr());
    else
      fputs (", not traced?", qstderr());
    fputs (", aborting...\n", qstderr());
    abort();
    return ptr;
  }
  else {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (qstderr(), "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
    return d;
  }
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (qstderr(),
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (qstderr(), "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (qstderr(),
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (qstderr(), "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", qstderr());
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if 0 == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif 0

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if _MPI
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if _MPI
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if _MPI
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if _MPI
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if 201307

#define OMP(x) _Pragma(#x)
#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#elif _MPI

#define OMP(x)

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  assert (!in_prof); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 645

#define prof_stop()\
  assert (in_prof); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 650


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/vinlinux/basilisk/src/common.h", 659);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/vinlinux/basilisk/src/common.h", 660);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/vinlinux/basilisk/src/common.h", 661); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 669

#define mpi_all_reduce_double(v,op) {\
  prof_start ("mpi_all_reduce");\
  double global, tmp = v;\
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);\
  v = global;\
  prof_stop();\
}\

#line 677


#endif

#define QFILE FILE

static FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

static FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL){
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = qstderr();
      fout = qstdout();
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define OMP(x)
#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#endif

void init_solver()
{
#if _CADNA
  cadna_init (-1);
#endif
#if _MPI
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define val(a,k,l,m) data(k,l,m)[a.i]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 803 "/home/vinlinux/basilisk/src/common.h"
#if (1 || __APPLE__) && !201307 && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;
#line 882 "/home/vinlinux/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 885

    norm += sq(n->x);
#line 885

    norm += sq(n->y);
#line 885

    norm += sq(n->z);}
  norm = sqrt(norm);
  {
#line 888

    n->x /= norm;
#line 888

    n->y /= norm;
#line 888

    n->z /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 911 "/home/vinlinux/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define dirichlet(x) (2.*(x) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(x) (Delta*(x) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#line 1 "/home/vinlinux/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/vinlinux/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 927 "/home/vinlinux/basilisk/src/common.h"



#include "_attributes.h"






















int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  if (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  if (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    if (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    if (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  if (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  if (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  if (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  if (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 1033

      if (w.x.i != v.x.i)
 id = false;
#line 1033

      if (w.y.i != v.y.i)
 id = false;
#line 1033

      if (w.z.i != v.z.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    if (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1056
 {
      assert (s->i >= 0);
      v.x = *s++;
    }
#line 1056
 {
      assert (s->i >= 0);
      v.y = *s++;
    }
#line 1056
 {
      assert (s->i >= 0);
      v.z = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  if (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1087
 {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
#line 1087
 {
      assert (v->y.i >= 0);
      t.y = *v++;
    }
#line 1087
 {
      assert (v->z.i >= 0);
      t.z = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;



scalar (* init_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if _MPI
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if _MPI
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}};
vector unityf= {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
scalar unity= {_NVARMAX + 6};
scalar zeroc= {_NVARMAX + 7};



 vector fm = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar cm = {(_NVARMAX + 6)};



static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}
#line 13 "main-cpp.c"
#include "_boundarydecl.h"
#line 1 "main.c"

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#line 1 "grid/octree.h"
#line 1 "/home/vinlinux/basilisk/src/grid/octree.h"


#line 1 "grid/tree.h"
#line 1 "/home/vinlinux/basilisk/src/grid/tree.h"
#line 1 "grid/mempool.h"
#line 1 "/home/vinlinux/basilisk/src/grid/mempool.h"





typedef struct _Pool Pool;

struct _Pool {
  Pool * next;
};

typedef struct {
  char * first, * lastb;
  size_t size;
  size_t poolsize;
  Pool * pool, * last;
} Mempool;

typedef struct {
  char * next;
} FreeBlock;

Mempool * mempool_new (size_t poolsize, size_t size)
{

  assert (poolsize % 8 == 0);
  assert (size >= sizeof(FreeBlock));


  poolsize = min(1 << 20, poolsize + sizeof(Pool));
  Mempool * m = ((Mempool *) pcalloc (1, sizeof(Mempool),__func__,__FILE__,__LINE__));
  m->poolsize = poolsize;
  m->size = size;
  return m;
}

void mempool_destroy (Mempool * m)
{
  Pool * p = m->pool;
  while (p) {
    Pool * next = p->next;
    pfree (p,__func__,__FILE__,__LINE__);
    p = next;
  }
  pfree (m,__func__,__FILE__,__LINE__);
}

void * mempool_alloc (Mempool * m)
{
  if (!m->first) {

    Pool * p = (Pool *) pmalloc (m->poolsize,__func__,__FILE__,__LINE__);
    p->next = NULL;
    if (m->last)
      m->last->next = p;
    else
      m->pool = p;
    m->last = p;
    m->first = m->lastb = ((char *)m->last) + sizeof(Pool);
    FreeBlock * b = (FreeBlock *) m->first;
    b->next = NULL;
  }
  void * ret = m->first;
  FreeBlock * b = (FreeBlock *) ret;
  char * next = b->next;
  if (!next) {
    m->lastb += m->size;
    next = m->lastb;
    if (next + m->size > ((char *) m->last) + m->poolsize)
      next = NULL;
    else {
      FreeBlock * b = (FreeBlock *) next;
      b->next = NULL;
    }
  }
  m->first = next;
#if TRASH
  double * v = (double *) ret;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  return ret;
}

void * mempool_alloc0 (Mempool * m)
{
  void * ret = mempool_alloc (m);
  memset (ret, 0, m->size);
  return ret;
}

void mempool_free (Mempool * m, void * p)
{
#if TRASH
  double * v = (double *) p;
  for (int i = 0; i < m->size/sizeof(double); i++)
    v[i] = undefined;
#endif
  FreeBlock * b = (FreeBlock *) p;
  b->next = m->first;
  m->first = (char *) p;
}
#line 2 "/home/vinlinux/basilisk/src/grid/tree.h"
#line 10 "/home/vinlinux/basilisk/src/grid/tree.h"
# define BGHOSTS 1
#line 22 "/home/vinlinux/basilisk/src/grid/tree.h"
typedef struct {
  unsigned short flags;

  unsigned short neighbors;
  int pid;
} Cell;

enum {
  active = 1 << 0,
  leaf = 1 << 1,
  border = 1 << 2,
  vertex = 1 << 3,
  user = 4,

  face_x = 1 << 0

  , face_y = 1 << 1


  , face_z = 1 << 2

};

#define is_active(cell) ((cell).flags & active)
#define is_leaf(cell) ((cell).flags & leaf)
#define is_coarse() ((cell).neighbors > 0)
#define is_border(cell) ((cell).flags & border)
#define is_local(cell) ((cell).pid == pid())



typedef struct {
  int i;

  int j;


  int k;

} IndexLevel;

typedef struct {
  IndexLevel * p;
  int n, nm;
} CacheLevel;

typedef struct {
  int i;

  int j;


  int k;

  int level, flags;
} Index;

typedef struct {
  Index * p;
  int n, nm;
} Cache;




static char * new_refarray (size_t len, size_t size) {
  return (char *) pcalloc (len + 1, size,__func__,__FILE__,__LINE__);
}

static void refarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)++;
}

static bool unrefarray (void * p, size_t len, size_t size) {
  int * refcount = (int *)(((char *)p) + len*size);
  (*refcount)--;
  if (*refcount == 0) {
    pfree (p,__func__,__FILE__,__LINE__);
    return true;
  }
  return false;
}




typedef struct {





  char **** m;

  Mempool * pool;
  int nc;
  int len;
} Layer;

static size_t _size (size_t depth)
{
  return (1 << depth) + 2*2;
}

static size_t poolsize (size_t depth, size_t size)
{






  return cube(_size(depth))*size;

}


#line 139

static inline
void assign_periodic_x (void ** m, int i, int nl, void * b)
{
  m[i] = b;
  if (Period.x) {
    for (int j = i; j < nl + 2*2; j += nl)
      m[j] = b;
    for (int j = i - nl; j >= 0; j -= nl)
      m[j] = b;
  }
}
#line 139

static inline
void assign_periodic_y (void ** m, int i, int nl, void * b)
{
  m[i] = b;
  if (Period.y) {
    for (int j = i; j < nl + 2*2; j += nl)
      m[j] = b;
    for (int j = i - nl; j >= 0; j -= nl)
      m[j] = b;
  }
}
#line 139

static inline
void assign_periodic_z (void ** m, int i, int nl, void * b)
{
  m[i] = b;
  if (Period.z) {
    for (int j = i; j < nl + 2*2; j += nl)
      m[j] = b;
    for (int j = i - nl; j >= 0; j -= nl)
      m[j] = b;
  }
}

static Layer * new_layer (int depth)
{
  Layer * l = ((Layer *) pmalloc ((1)*sizeof(Layer),__func__,__FILE__,__LINE__));
  l->len = _size (depth);
  if (depth == 0)
    l->pool = NULL;
  else {
    size_t size = sizeof(Cell) + datasize;


    l->pool = mempool_new (poolsize (depth, size), (1 << 3)*size);
  }





  l->m = ((char *** *) pcalloc (l->len, sizeof(char ***),__func__,__FILE__,__LINE__));

  l->nc = 0;
  return l;
}

static void destroy_layer (Layer * l)
{
  if (l->pool)
    mempool_destroy (l->pool);
  pfree (l->m,__func__,__FILE__,__LINE__);
  pfree (l,__func__,__FILE__,__LINE__);
}
#line 199 "/home/vinlinux/basilisk/src/grid/tree.h"
static void layer_add_row (Layer * l, int i, int j)
{
  if (!l->m[i]) {
    assign_periodic_x ((void **) l->m, i, l->len - 2*2,
         (void *) new_refarray (l->len, sizeof (char *)));
    l->nc++;
  }
  refarray (l->m[i], l->len, sizeof(char *));

  if (!l->m[i][j])
    assign_periodic_y ((void **) l->m[i], j, l->len - 2*2,
         (void *) new_refarray (l->len, sizeof (char *)));
  refarray (l->m[i][j], l->len, sizeof(char *));

}

static bool layer_remove_row (Layer * l, int i, int j)
{

  if (unrefarray (l->m[i][j], l->len, sizeof (char *)))
    assign_periodic_y ((void **) l->m[i], j, l->len - 2*2, NULL);

  if (unrefarray (l->m[i], l->len, sizeof (char *))) {
    assign_periodic_x ((void **) l->m, i, l->len - 2*2, NULL);
    if (--l->nc == 0) {
      destroy_layer (l);
      return true;
    }
    assert (l->nc >= 0);
  }
  return false;
}




typedef struct {
  Grid g;
  Layer ** L;

  Cache leaves;
  Cache faces;
  Cache vertices;
  Cache refined;
  CacheLevel * active;
  CacheLevel * prolongation;
  CacheLevel * boundary;

  CacheLevel * restriction;

  bool dirty;
} Tree;



struct _Point {

  int i;

  int j;


  int k;

  int level;
};
static Point last_point;



static void cache_level_append (CacheLevel * c, Point p)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (IndexLevel *) prealloc (c->p, (c->nm)*sizeof(IndexLevel),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->n++;
}

static void cache_level_shrink (CacheLevel * c)
{
  if (c->nm > (c->n/128 + 1)*128) {
    c->nm = (c->n/128 + 1)*128;
    assert (c->nm > c->n);
    c->p = (IndexLevel *) prealloc (c->p, sizeof (Index)*c->nm,__func__,__FILE__,__LINE__);
  }
}

static void cache_append (Cache * c, Point p, unsigned short flags)
{
  if (c->n >= c->nm) {
    c->nm += 128;
    c->p = (Index *) prealloc (c->p, (c->nm)*sizeof(Index),__func__,__FILE__,__LINE__);
  }
  c->p[c->n].i = p.i;

  c->p[c->n].j = p.j;


  c->p[c->n].k = p.k;

  c->p[c->n].level = p.level;
  c->p[c->n].flags = flags;
  c->n++;
}

void cache_shrink (Cache * c)
{
  cache_level_shrink ((CacheLevel *)c);
}
#line 371 "/home/vinlinux/basilisk/src/grid/tree.h"
#define allocated(a,l,n) (point.i+a >= 0 &&\
         point.i+a < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a] &&\
         point.j+l >= 0 &&\
         point.j+l < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l] &&\
         point.k+n >= 0 &&\
         point.k+n < (1 << point.level) + 2*2 &&\
         ((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
         [point.k+n])\

#line 381


#define NEIGHBOR(a,l,n) (((Tree *)grid)->L[point.level]->m[point.i+a][point.j+l]\
                                          [point.k+n])\

#line 385

#define PARENT(a,l,n) (((Tree *)grid)->L[point.level-1]->m[(point.i+2)/2+a]\
    [(point.j+2)/2+l][(point.k+2)/2+n])\

#line 388

#define allocated_child(a,l,n) (level < depth() &&\
         point.i > 0 && point.i <= (1 << level) + 2 &&\
         point.j > 0 && point.j <= (1 << level) + 2 &&\
         point.k > 0 && point.k <= (1 << level) + 2 &&\
         ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
   && ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
   && ((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a][2*point.j-2 +l]\
         [2*point.k-2 +n])\

#line 397

#define CHILD(a,l,n) (((Tree *)grid)->L[point.level+1]->m[2*point.i-2 +a]\
      [2*point.j-2 +l][2*point.k-2 +n])\

#line 400


#define CELL(m) (*((Cell *)(m)))


#define depth() (grid->depth)
#define aparent(k,l,n) CELL(PARENT(k,l,n))
#define child(k,l,n) CELL(CHILD(k,l,n))


#define cell CELL(NEIGHBOR(0,0,0))
#define neighbor(k,l,n) CELL(NEIGHBOR(k,l,n))
#define neighborp(l,m,n) (Point) {\
    point.i + l,\
\
    point.j + m,\
\
\
    point.k + n,\
\
    point.level }\

#line 421



#define data(k,l,n) ((double *) (NEIGHBOR(k,l,n) + sizeof(Cell)))
#define fine(a,k,l,n) ((double *) (CHILD(k,l,n) + sizeof(Cell)))[a.i]
#define coarse(a,k,l,n) ((double *) (PARENT(k,l,n) + sizeof(Cell)))[a.i]

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
\
\
\
\
\
\
\
  struct { int x, y, z; } child = {\
    2*((point.i+2)%2)-1, 2*((point.j+2)%2)-1, 2*((point.k+2)%2)-1\
  };\
\
  NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
\
  parent.j = (point.j + 2)/2;\
\
\
  parent.k = (point.k + 2)/2;\
\
\
\
\
\
\

#line 457


#line 1 "grid/foreach_cell.h"
#line 1 "/home/vinlinux/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/vinlinux/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = { .level = 0 };\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
      for (root.k = 2*Period.z; root.k <= 2*(2 - Period.z); root.k++)\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point;\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
\
\
    Point root = {2,2,2,0};\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = { .level = 0 };\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 460 "/home/vinlinux/basilisk/src/grid/tree.h"
#line 495 "/home/vinlinux/basilisk/src/grid/tree.h"
#define foreach_child() {\
  int _i = 2*point.i - 2, _j = 2*point.j - 2, _k = 2*point.k - 2;\
  point.level++;\
  for (int _l = 0; _l < 2; _l++) {\
    point.i = _i + _l;\
    for (int _m = 0; _m < 2; _m++) {\
      point.j = _j + _m;\
      for (int _n = 0; _n < 2; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 505

#define end_foreach_child()\
      }\
    }\
  }\
  point.i = (_i + 2)/2;point.j = (_j + 2)/2;point.k = (_k + 2)/2;\
  point.level--;\
}\

#line 513

#define foreach_child_break() _l = _m = _n = 2
#line 523 "/home/vinlinux/basilisk/src/grid/tree.h"
#define is_refined_check() ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) &&\
    point.i > 0 && point.i < (1 << level) + 2*2 - 1\
\
    && point.j > 0 && point.j < (1 << level) + 2*2 - 1\
\
\
    && point.k > 0 && point.k < (1 << level) + 2*2 - 1\
\
    )\

#line 532


#define foreach_cache(_cache) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  int _k; unsigned short _flags; NOT_UNUSED(_flags);\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    point.level = _cache.p[_k].level;\
    _flags = _cache.p[_k].flags;\
    POINT_VARIABLES;\

#line 557

#define end_foreach_cache() } } }

#define foreach_cache_level(_cache,_l) {\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
\
\
\
\
\
  Point point = {2,2,2,0};\
\
  point.level = _l;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 0; _k < _cache.n; _k++) {\
    point.i = _cache.p[_k].i;\
\
    point.j = _cache.p[_k].j;\
\
\
    point.k = _cache.p[_k].k;\
\
    POINT_VARIABLES;\

#line 582

#define end_foreach_cache_level() } } }

#define foreach_boundary_level(_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _boundary = ((Tree *)grid)->boundary[_l];\
    foreach_cache_level (_boundary,_l)\

#line 590

#define end_foreach_boundary_level() end_foreach_cache_level(); }}



#define foreach_boundary(_b) {\
  for (int _l = depth(); _l >= 0; _l--)\
    foreach_boundary_level(_l) {\
      if ((- cell.pid - 1) == _b)\
 for (int _d = 0; _d < 3; _d++) {\
   for (int _i = -1; _i <= 1; _i += 2) {\
     if (_d == 0) ig = _i; else if (_d == 1) jg = _i; else kg = _i;\
     if (allocated(-ig,-jg,-kg) &&\
  is_leaf (neighbor(-ig,-jg,-kg)) &&\
  !(neighbor(-ig,-jg,-kg).pid < 0) &&\
  is_local(neighbor(-ig,-jg,-kg))) {\
       point.i -= ig; x -= ig*Delta/2.;\
\
       point.j -= jg; y -= jg*Delta/2.;\
\
\
       point.k -= kg; z -= kg*Delta/2.;\
\

#line 613

#define end_foreach_boundary()\
       point.i += ig; x += ig*Delta/2.;\
\
       point.j += jg; y += jg*Delta/2.;\
\
\
       point.k += kg; z += kg*Delta/2.;\
\
            }\
   }\
   ig = jg = kg = 0;\
 }\
    } end_foreach_boundary_level(); }\

#line 627


#define foreach_halo(_name,_l) {\
  if (_l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _cache = ((Tree *)grid)->_name[_l];\
    foreach_cache_level (_cache, _l)\

#line 634

#define end_foreach_halo() end_foreach_cache_level(); }}

#line 1 "grid/neighbors.h"
#line 1 "/home/vinlinux/basilisk/src/grid/neighbors.h"
#line 35 "/home/vinlinux/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 638 "/home/vinlinux/basilisk/src/grid/tree.h"

static inline bool has_local_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 640 "/home/vinlinux/basilisk/src/grid/tree.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}

static inline void cache_append_face (Point point, unsigned short flags)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 648 "/home/vinlinux/basilisk/src/grid/tree.h"

  Tree * q = ((Tree *)grid);
  cache_append (&q->faces, point, flags);
#line 662 "/home/vinlinux/basilisk/src/grid/tree.h"
  {
#line 662

    if (flags & face_x)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(0,i,j).flags & vertex)) {
     cache_append (&q->vertices, neighborp(0,i,j), 0);
     neighbor(0,i,j).flags |= vertex;
   }
#line 662

    if (flags & face_y)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(j,0,i).flags & vertex)) {
     cache_append (&q->vertices, neighborp(j,0,i), 0);
     neighbor(j,0,i).flags |= vertex;
   }
#line 662

    if (flags & face_z)
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   if (!(neighbor(i,j,0).flags & vertex)) {
     cache_append (&q->vertices, neighborp(i,j,0), 0);
     neighbor(i,j,0).flags |= vertex;
   }}

}



void check_periodic (Tree * q)
{
#line 702 "/home/vinlinux/basilisk/src/grid/tree.h"
}

static void update_cache_f (void)
{
  Tree * q = ((Tree *)grid);

  check_periodic (q);

   { foreach_cache (q->vertices){

#line 710 "/home/vinlinux/basilisk/src/grid/tree.h"

    if (level <= depth() && allocated(0,0,0))
      cell.flags &= ~vertex; } end_foreach_cache(); }


  q->leaves.n = q->faces.n = q->vertices.n = 0;
  for (int l = 0; l <= depth(); l++)
    q->active[l].n = q->prolongation[l].n =
      q->boundary[l].n = q->restriction[l].n = 0;

  const unsigned short fboundary = 1 << user;
   { foreach_cell(){

#line 721 "/home/vinlinux/basilisk/src/grid/tree.h"
 {



    if (is_local(cell) && is_active(cell)) {


      cache_level_append (&q->active[level], point);
    }
#line 745 "/home/vinlinux/basilisk/src/grid/tree.h"
    if (!(cell.pid < 0)) {

       { foreach_neighbor (BGHOSTS)
 if (allocated(0,0,0) && (cell.pid < 0) && !(cell.flags & fboundary)) {
   cache_level_append (&q->boundary[level], point);
   cell.flags |= fboundary;
 } end_foreach_neighbor(); }
    }

    else if (level > 0 && is_local(aparent(0,0,0)))
      cache_level_append (&q->restriction[level], point);

    if (is_leaf (cell)) {
      if (is_local(cell)) {
 cache_append (&q->leaves, point, 0);

 unsigned short flags = 0;
 {
#line 762

   if ((neighbor(-1,0,0).pid < 0) || (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) ||
       is_leaf(neighbor(-1,0,0)))
     flags |= face_x;
#line 762

   if ((neighbor(0,-1,0).pid < 0) || (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) ||
       is_leaf(neighbor(0,-1,0)))
     flags |= face_y;
#line 762

   if ((neighbor(0,0,-1).pid < 0) || (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) ||
       is_leaf(neighbor(0,0,-1)))
     flags |= face_z;}
 if (flags)
   cache_append (&q->faces, point, flags);
 {
#line 768

   if ((neighbor(1,0,0).pid < 0) || (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) ||
       (!is_local(neighbor(1,0,0)) && is_leaf(neighbor(1,0,0))))
     cache_append (&q->faces, neighborp(1,0,0), face_x);
#line 768

   if ((neighbor(0,1,0).pid < 0) || (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) ||
       (!is_local(neighbor(0,1,0)) && is_leaf(neighbor(0,1,0))))
     cache_append (&q->faces, neighborp(0,1,0), face_y);
#line 768

   if ((neighbor(0,0,1).pid < 0) || (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) ||
       (!is_local(neighbor(0,0,1)) && is_leaf(neighbor(0,0,1))))
     cache_append (&q->faces, neighborp(0,0,1), face_z);}

 for (int i = 0; i <= 1; i++)

   for (int j = 0; j <= 1; j++)


     for (int k = 0; k <= 1; k++)

       if (!(neighbor(i,j,k).flags & vertex)) {
  cache_append (&q->vertices, neighborp(i,j,k), 0);
  neighbor(i,j,k).flags |= vertex;
       }

        if (cell.neighbors > 0)
   cache_level_append (&q->prolongation[level], point);
      }
      else if (!(cell.pid < 0) || is_local(aparent(0,0,0))) {

 unsigned short flags = 0;
 {
#line 791

   if (allocated(-1,0,0) &&
       is_local(neighbor(-1,0,0)) && (!is_leaf(neighbor(-1,0,0)) && !neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     flags |= face_x;
#line 791

   if (allocated(0,-1,0) &&
       is_local(neighbor(0,-1,0)) && (!is_leaf(neighbor(0,-1,0)) && !neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     flags |= face_y;
#line 791

   if (allocated(0,0,-1) &&
       is_local(neighbor(0,0,-1)) && (!is_leaf(neighbor(0,0,-1)) && !neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     flags |= face_z;}
 if (flags)
   cache_append_face (point, flags);
 {
#line 797

   if (allocated(1,0,0) && is_local(neighbor(1,0,0)) &&
       (!is_leaf(neighbor(1,0,0)) && !neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     cache_append_face (neighborp(1,0,0), face_x);
#line 797

   if (allocated(0,1,0) && is_local(neighbor(0,1,0)) &&
       (!is_leaf(neighbor(0,1,0)) && !neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     cache_append_face (neighborp(0,1,0), face_y);
#line 797

   if (allocated(0,0,1) && is_local(neighbor(0,0,1)) &&
       (!is_leaf(neighbor(0,0,1)) && !neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     cache_append_face (neighborp(0,0,1), face_z);}
      }

      continue;

    }
  } } end_foreach_cell(); }


  cache_shrink (&q->leaves);
  cache_shrink (&q->faces);
  cache_shrink (&q->vertices);
  for (int l = 0; l <= depth(); l++) {
    cache_level_shrink (&q->active[l]);
    cache_level_shrink (&q->prolongation[l]);
    cache_level_shrink (&q->boundary[l]);
    cache_level_shrink (&q->restriction[l]);
}

  q->dirty = false;


  for (int l = depth(); l >= 0; l--)
     { foreach_boundary_level (l){

#line 823 "/home/vinlinux/basilisk/src/grid/tree.h"

      cell.flags &= ~fboundary; } end_foreach_boundary_level(); }



  grid->n = q->leaves.n;

#if !_MPI
  grid->tn = grid->n;
  grid->maxdepth = grid->depth;
#endif
}

#define foreach() { if (((Tree *)grid)->dirty) update_cache_f(); }; foreach_cache(((Tree *)grid)->leaves)
#define end_foreach() end_foreach_cache()

#define foreach_face_generic()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->faces) 
#line 840

#define end_foreach_face_generic() end_foreach_cache()

#define is_face_x() (_flags & face_x)

#define is_face_y() (_flags & face_y)


#define is_face_z() (_flags & face_z)


#define foreach_vertex()\
  { if (((Tree *)grid)->dirty) update_cache_f(); };\
  foreach_cache(((Tree *)grid)->vertices) {\
    x -= Delta/2.;\
\
    y -= Delta/2.;\
\
\
    z -= Delta/2.;\
\

#line 862

#define end_foreach_vertex() } end_foreach_cache()
#line 874 "/home/vinlinux/basilisk/src/grid/tree.h"
#define foreach_level(l) {\
  if (l <= depth()) {\
    { if (((Tree *)grid)->dirty) update_cache_f(); };\
    CacheLevel _active = ((Tree *)grid)->active[l];\
    foreach_cache_level (_active,l)\

#line 879

#define end_foreach_level() end_foreach_cache_level(); }}

#define foreach_coarse_level(l) foreach_level(l) if (!is_leaf(cell)) {
#define end_foreach_coarse_level() } end_foreach_level()

#define foreach_level_or_leaf(l) {\
  for (int _l1 = l; _l1 >= 0; _l1--)\
    foreach_level(_l1)\
      if (_l1 == l || is_leaf (cell)) {\

#line 889

#define end_foreach_level_or_leaf() } end_foreach_level(); }

#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Tree * q = ((Tree *)grid);

  for (int l = 0; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = 0; i < L->len; i++)
      if (L->m[i])





 for (int j = 0; j < L->len; j++)
   if (L->m[i][j])





          for (int k = 0; k < L->len; k++)
     if (L->m[i][j][k])
       if (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10)
           if (!is_constant(s))
    ((double *)(L->m[i][j][k] + sizeof(Cell)))[s.i] = val;


  }
}

#define cache_level_resize(name, a)\
{\
  for (int i = 0; i <= depth() - a; i++)\
    pfree (q->name[i].p,__func__,__FILE__,__LINE__);\
  pfree (q->name,__func__,__FILE__,__LINE__);\
  q->name = ((CacheLevel *) pcalloc (depth() + 1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));\
}\

#line 935


static void update_depth (int inc)
{
  Tree * q = ((Tree *)grid);
  grid->depth += inc;
  q->L = &(q->L[-1]);
  q->L = (Layer * *) prealloc (q->L, (grid->depth + 2)*sizeof(Layer *),__func__,__FILE__,__LINE__);
  q->L = &(q->L[1]);
  if (inc > 0)
    q->L[grid->depth] = new_layer (grid->depth);
  cache_level_resize (active, inc);
  cache_level_resize (prolongation, inc);
  cache_level_resize (boundary, inc);
  cache_level_resize (restriction, inc);
}

static void alloc_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 953 "/home/vinlinux/basilisk/src/grid/tree.h"

  if (point.level == grid->depth)
    update_depth (+1);
  else if (allocated_child(0,0,0))
    return;


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  size_t len = sizeof(Cell) + datasize;
  char * b = (char *) mempool_alloc0 (L->pool);
  int nl = L->len - 2*2;
  int i = 2*point.i - 2;
  for (int k = 0; k < 2; k++, i++) {
#line 980 "/home/vinlinux/basilisk/src/grid/tree.h"
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      layer_add_row (L, i, j);
      int m = 2*point.k - 2;
      for (int n = 0; n < 2; n++, m++) {
 assert (!L->m[i][j][m]);
 assign_periodic_z ((void **) L->m[i][j], m, nl, (void *) b);
 b += len;
      }
    }

  }

  int pid = cell.pid;
   { foreach_child() {
    cell.pid = pid;
#if TRASH
    if (all) for (scalar s = *all, *_i11 = all; ((scalar *)&s)->i >= 0; s = *++_i11)
      val(s,0,0,0) = undefined;
#endif
  } end_foreach_child(); }
}
#line 1038 "/home/vinlinux/basilisk/src/grid/tree.h"
static void free_children (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1039 "/home/vinlinux/basilisk/src/grid/tree.h"


  Layer * L = ((Tree *)grid)->L[point.level + 1];
  int i = 2*point.i - 2, nl = L->len - 2*2;
  assert (L->m[i][2*point.j - 2][2*point.k - 2]);
  mempool_free (L->pool, L->m[i][2*point.j - 2][2*point.k - 2]);
  for (int k = 0; k < 2; k++, i++) {
    int j = 2*point.j - 2;
    for (int l = 0; l < 2; l++, j++) {
      int m = 2*point.k - 2;
      for (int n = 0; n < 2; n++, m++)
 assign_periodic_z ((void **) L->m[i][j], m, nl, NULL);
      if (layer_remove_row (L, i, j)) {
 assert (point.level + 1 == grid->depth);
 update_depth (-1);
      }
    }
  }
}


void increment_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1061 "/home/vinlinux/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
  if (cell.neighbors++ == 0)
    alloc_children (point);
   { foreach_neighbor (2/2)
    if (cell.neighbors++ == 0)
      alloc_children (point); end_foreach_neighbor(); }
  cell.neighbors--;
}

void decrement_neighbors (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1072 "/home/vinlinux/basilisk/src/grid/tree.h"

  ((Tree *)grid)->dirty = true;
   { foreach_neighbor (2/2)
    if (allocated(0,0,0)) {
      cell.neighbors--;
      if (cell.neighbors == 0)
 free_children (point);
    } end_foreach_neighbor(); }
  if (cell.neighbors) {
    int pid = cell.pid;
     { foreach_child() {
      cell.flags = 0;
      cell.pid = pid;
    } end_foreach_child(); }
  }
}

static void apply_periodic_elem (char ** m, int len)
{
  if (m) {
    int end = len - 2;
    for (int k = 0; k < 2; k++) {
      m[k] = m[k + end - 2];
      m[end + k] = m[k + 2];
    }
  }
}

static void apply_periodic (Tree * q)
{
#line 1121 "/home/vinlinux/basilisk/src/grid/tree.h"
  if (Period.z) {
    for (int i = 0; i < q->L[0]->len; i++)
      for (int j = 0; j < q->L[0]->len; j++)
 for (int k = 0; k < q->L[0]->len; k++)
   q->L[0]->m[i][j][k] = q->L[0]->m[i][j][2];
    for (int l = 1; l <= depth(); l++) {
      Layer * L = q->L[l];
      for (int i = 0; i < L->len; i++)
 if (L->m[i])
   for (int j = 0; j < L->len; j++)
     apply_periodic_elem (L->m[i][j], L->len);
    }
  }

}

void realloc_scalar (void)
{

  Tree * q = ((Tree *)grid);
  size_t newlen = sizeof(Cell) + datasize;
  size_t oldlen = newlen - sizeof(double);

  Layer * L = q->L[0];
  int len = L->len;
  for (int i = Period.x*2; i < len - Period.x*2; i++) {



    for (int j = Period.y*2; j < len - Period.y*2; j++) {



      for (int k = Period.z*2; k < len - Period.z*2; k++)
 L->m[i][j][k] = (char *) prealloc (L->m[i][j][k], (newlen)*sizeof(char),__func__,__FILE__,__LINE__);

    }

  }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    int len = L->len;
    Mempool * oldpool = L->pool;
    L->pool = mempool_new (poolsize (l, newlen), (1 << 3)*newlen);
    for (int i = Period.x*2; i < len - Period.x*2; i += 2)
      if (L->m[i]) {
#line 1176 "/home/vinlinux/basilisk/src/grid/tree.h"
 for (int j = Period.y*2; j < len - Period.y*2; j += 2)
   if (L->m[i][j]) {
#line 1187 "/home/vinlinux/basilisk/src/grid/tree.h"
     for (int k = Period.z*2; k < len - Period.z*2; k += 2)
       if (L->m[i][j][k]) {
  char * new = (char *) mempool_alloc (L->pool);
  for (int l = 0; l < 2; l++)
    for (int m = 0; m < 2; m++)
      for (int n = 0; n < 2; n++) {
        memcpy (new, L->m[i+l][j+m][k+n], oldlen);
        L->m[i+l][j+m][k+n] = new;
        new += newlen;
      }
       }

   }

      }
    mempool_destroy (oldpool);
  }
  apply_periodic (q);
  check_periodic (q);
}



#define VN v.x
#define VT v.y
#define VR v.z




#if _MPI
# define disable_fpe_for_mpi() disable_fpe (FE_DIVBYZERO|FE_INVALID)
# define enable_fpe_for_mpi() enable_fpe (FE_DIVBYZERO|FE_INVALID)
#else
# define disable_fpe_for_mpi()
# define enable_fpe_for_mpi()
#endif

static inline void no_restriction (Point point, scalar s);

static bool normal_neighbor (Point point, scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1228 "/home/vinlinux/basilisk/src/grid/tree.h"

  for (int k = 1; k <= BGHOSTS; k++)
    {
#line 1230

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
   Point neighbor = neighborp(i,0,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.x,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);

     scalar vt = VT;
     val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);


     scalar vr = VR;
     val(v.z,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.z);

   }
   return true;
 }
#line 1230

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
   Point neighbor = neighborp(0,i,0);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.y,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);

     scalar vt = VT;
     val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);


     scalar vr = VR;
     val(v.x,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.x);

   }
   return true;
 }
#line 1230

      for (int i = -k; i <= k; i += 2*k)
 if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
   Point neighbor = neighborp(0,0,i);
   int id = (- cell.pid - 1);
   if (scalars) for (scalar s = *scalars, *_i12 = scalars; ((scalar *)&s)->i >= 0; s = *++_i12)
     val(s,0,0,0) = _attribute[s.i].boundary[id](neighbor, point, s);
   if (vectors) for (vector v = *vectors, *_i13 = vectors; ((scalar *)&v)->i >= 0; v = *++_i13) {
     scalar vn = VN;
     val(v.z,0,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.z);

     scalar vt = VT;
     val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);


     scalar vr = VR;
     val(v.y,0,0,0) = _attribute[vr.i].boundary[id](neighbor, point, v.y);

   }
   return true;
 }}
  return false;
}

static bool diagonal_neighbor_2D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1256 "/home/vinlinux/basilisk/src/grid/tree.h"


  for (int k = 1; k <= BGHOSTS; k++)

    {
#line 1260


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(i,j,0) && (allocated(i,j,0) && !(neighbor(i,j,0).pid < 0)) &&
       allocated(i,0,0) && (neighbor(i,0,0).pid < 0) &&
       allocated(0,j,0) && (neighbor(0,j,0).pid < 0)) {
     Point n = neighborp(i,j,0),
       n1 = neighborp(i,0,0), n2 = neighborp(0,j,0);
     int id1 = (- neighbor(i,0,0).pid - 1), id2 = (- neighbor(0,j,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,i,j,0));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.x) +
         _attribute[vn.i].boundary[id2](n,n2,v.x) -
         val(v.x,i,j,0));
       val(v.y,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.y) +
         _attribute[vt.i].boundary[id2](n,n2,v.y) -
         val(v.y,i,j,0));

       scalar vr = VR;
       val(v.z,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.z) +
         _attribute[vr.i].boundary[id2](n,n2,v.z) -
         val(v.z,i,j,0));

     }
     return true;
   }
#line 1260


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(0,i,j) && (allocated(0,i,j) && !(neighbor(0,i,j).pid < 0)) &&
       allocated(0,i,0) && (neighbor(0,i,0).pid < 0) &&
       allocated(0,0,j) && (neighbor(0,0,j).pid < 0)) {
     Point n = neighborp(0,i,j),
       n1 = neighborp(0,i,0), n2 = neighborp(0,0,j);
     int id1 = (- neighbor(0,i,0).pid - 1), id2 = (- neighbor(0,0,j).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,0,i,j));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.y) +
         _attribute[vn.i].boundary[id2](n,n2,v.y) -
         val(v.y,0,i,j));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.z) +
         _attribute[vt.i].boundary[id2](n,n2,v.z) -
         val(v.z,0,i,j));

       scalar vr = VR;
       val(v.x,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.x) +
         _attribute[vr.i].boundary[id2](n,n2,v.x) -
         val(v.x,0,i,j));

     }
     return true;
   }
#line 1260


      for (int i = -k; i <= k; i += 2*k)
 for (int j = -k; j <= k; j += 2*k)
   if (allocated(j,0,i) && (allocated(j,0,i) && !(neighbor(j,0,i).pid < 0)) &&
       allocated(0,0,i) && (neighbor(0,0,i).pid < 0) &&
       allocated(j,0,0) && (neighbor(j,0,0).pid < 0)) {
     Point n = neighborp(j,0,i),
       n1 = neighborp(0,0,i), n2 = neighborp(j,0,0);
     int id1 = (- neighbor(0,0,i).pid - 1), id2 = (- neighbor(j,0,0).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i14 = scalars; ((scalar *)&s)->i >= 0; s = *++_i14)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n,n1,s) + _attribute[s.i].boundary[id2](n,n2,s) -
       val(s,j,0,i));
     if (vectors) for (vector v = *vectors, *_i15 = vectors; ((scalar *)&v)->i >= 0; v = *++_i15) {
       scalar vt = VT, vn = VN;
       val(v.z,0,0,0) = (_attribute[vt.i].boundary[id1](n,n1,v.z) +
         _attribute[vn.i].boundary[id2](n,n2,v.z) -
         val(v.z,j,0,i));
       val(v.x,0,0,0) = (_attribute[vn.i].boundary[id1](n,n1,v.x) +
         _attribute[vt.i].boundary[id2](n,n2,v.x) -
         val(v.x,j,0,i));

       scalar vr = VR;
       val(v.y,0,0,0) = (_attribute[vr.i].boundary[id1](n,n1,v.y) +
         _attribute[vr.i].boundary[id2](n,n2,v.y) -
         val(v.y,j,0,i));

     }
     return true;
   }}

  return false;
}

static bool diagonal_neighbor_3D (Point point,
      scalar * scalars, vector * vectors)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1296 "/home/vinlinux/basilisk/src/grid/tree.h"


  for (int n = 1; n <= BGHOSTS; n++)
    for (int i = -n; i <= n; i += 2*n)
      for (int j = -n; j <= n; j += 2*n)
 for (int k = -n; k <= n; k += 2*n)
   if ((allocated(i,j,k) && !(neighbor(i,j,k).pid < 0)) &&
       (neighbor(i,j,0).pid < 0) &&
       (neighbor(i,0,k).pid < 0) &&
       (neighbor(0,j,k).pid < 0)) {
     Point
       n0 = neighborp(i,j,k),
       n1 = neighborp(i,j,0),
       n2 = neighborp(i,0,k),
       n3 = neighborp(0,j,k);
     int
       id1 = (- neighbor(i,j,0).pid - 1),
       id2 = (- neighbor(i,0,k).pid - 1),
       id3 = (- neighbor(0,j,k).pid - 1);
     if (scalars) for (scalar s = *scalars, *_i16 = scalars; ((scalar *)&s)->i >= 0; s = *++_i16)
       val(s,0,0,0) = (_attribute[s.i].boundary[id1](n0,n1,s) +
       _attribute[s.i].boundary[id2](n0,n2,s) +
       _attribute[s.i].boundary[id3](n0,n3,s) -
       2.*val(s,i,j,k));
     if (vectors) for (vector v = *vectors, *_i17 = vectors; ((scalar *)&v)->i >= 0; v = *++_i17) {
       scalar vt = VT, vn = VN, vr = VR;
       val(v.x,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.x) +
         _attribute[vt.i].boundary[id2](n0,n2,v.x) +
         _attribute[vn.i].boundary[id3](n0,n3,v.x) -
         2.*val(v.x,i,j,k));
       val(v.y,0,0,0) = (_attribute[vt.i].boundary[id1](n0,n1,v.y) +
         _attribute[vn.i].boundary[id2](n0,n2,v.y) +
         _attribute[vt.i].boundary[id3](n0,n3,v.y) -
         2.*val(v.y,i,j,k));
       val(v.z,0,0,0) = (_attribute[vn.i].boundary[id1](n0,n1,v.z) +
         _attribute[vr.i].boundary[id2](n0,n2,v.z) +
         _attribute[vr.i].boundary[id3](n0,n3,v.z) -
         2.*val(v.z,i,j,k));
     }
     return true;
   }

  return false;
}



#line 1342

static Point tangential_neighbor_x (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1344 "/home/vinlinux/basilisk/src/grid/tree.h"

  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(-1,j,0) && !(neighbor(-1,j,0).pid < 0))) {
 *zn = false;
 return neighborp(0,j,0);
      }


      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(-1,0,j) && !(neighbor(-1,0,j).pid < 0))) {
 *zn = true;
 return neighborp(0,0,j);
      }

    }
  return (Point){.level = -1};
}
#line 1342

static Point tangential_neighbor_y (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1344 "/home/vinlinux/basilisk/src/grid/tree.h"

  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(0,0,j) && !(neighbor(0,0,j).pid < 0)) || (allocated(0,-1,j) && !(neighbor(0,-1,j).pid < 0))) {
 *zn = false;
 return neighborp(0,0,j);
      }


      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,-1,0) && !(neighbor(j,-1,0).pid < 0))) {
 *zn = true;
 return neighborp(j,0,0);
      }

    }
  return (Point){.level = -1};
}
#line 1342

static Point tangential_neighbor_z (Point point, bool * zn)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1344 "/home/vinlinux/basilisk/src/grid/tree.h"

  for (int k = 1; k <= BGHOSTS; k++)
    for (int j = -k; j <= k; j += 2*k) {
      if ((allocated(j,0,0) && !(neighbor(j,0,0).pid < 0)) || (allocated(j,0,-1) && !(neighbor(j,0,-1).pid < 0))) {
 *zn = false;
 return neighborp(j,0,0);
      }


      if ((allocated(0,j,0) && !(neighbor(0,j,0).pid < 0)) || (allocated(0,j,-1) && !(neighbor(0,j,-1).pid < 0))) {
 *zn = true;
 return neighborp(0,j,0);
      }

    }
  return (Point){.level = -1};
}


static inline bool is_boundary_point (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1363 "/home/vinlinux/basilisk/src/grid/tree.h"

  return (cell.pid < 0);
}

static void box_boundary_level (const Boundary * b, scalar * list, int l)
{
  disable_fpe_for_mpi();
  scalar * scalars = NULL;
  vector * vectors = NULL, * faces = NULL;
  if (list) for (scalar s = *list, *_i18 = list; ((scalar *)&s)->i >= 0; s = *++_i18)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i) {
 if (_attribute[s.i].face)
   faces = vectors_add (faces, _attribute[s.i].v);
 else
   vectors = vectors_add (vectors, _attribute[s.i].v);
      }
      else if (_attribute[s.i].v.x.i < 0)
 scalars = list_add (scalars, s);
    }

   { foreach_boundary_level (l){

#line 1384 "/home/vinlinux/basilisk/src/grid/tree.h"
 {
    if (!normal_neighbor (point, scalars, vectors) &&
 !diagonal_neighbor_2D (point, scalars, vectors) &&
 !diagonal_neighbor_3D (point, scalars, vectors)) {

      if (scalars) for (scalar s = *scalars, *_i19 = scalars; ((scalar *)&s)->i >= 0; s = *++_i19)
 val(s,0,0,0) = undefined;
      if (vectors) for (vector v = *vectors, *_i20 = vectors; ((scalar *)&v)->i >= 0; v = *++_i20)
 {
#line 1392

   val(v.x,0,0,0) = undefined;
#line 1392

   val(v.y,0,0,0) = undefined;
#line 1392

   val(v.z,0,0,0) = undefined;}
    }
    if (faces) {
      int id = (- cell.pid - 1);
      {
#line 1397

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(i,0,0) && !(neighbor(i,0,0).pid < 0))) {
     Point neighbor = neighborp(i,0,0);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.x,(i + 1)/2,0,0) = _attribute[vn.i].boundary[id](neighbor, point, v.x);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_x (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(-1,0,0).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.x,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.x);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.x,0,0,0) = 0.;
   }

 }
#line 1397

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,i,0) && !(neighbor(0,i,0).pid < 0))) {
     Point neighbor = neighborp(0,i,0);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.y,0,(i + 1)/2,0) = _attribute[vn.i].boundary[id](neighbor, point, v.y);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_y (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,-1,0).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.y,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.y);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.y,0,0,0) = 0.;
   }

 }
#line 1397

 for (int i = -1; i <= 1; i += 2) {

   if ((allocated(0,0,i) && !(neighbor(0,0,i).pid < 0))) {
     Point neighbor = neighborp(0,0,i);
     if (faces) for (vector v = *faces, *_i21 = faces; ((scalar *)&v)->i >= 0; v = *++_i21) {
       scalar vn = VN;
       if (_attribute[vn.i].boundary[id])
  val(v.z,0,0,(i + 1)/2) = _attribute[vn.i].boundary[id](neighbor, point, v.z);
     }
   }

   else if (i == -1) {

     bool zn;
     Point neighbor = tangential_neighbor_z (point, &zn);
     if (neighbor.level >= 0) {
       int id = is_boundary_point (neighbor) ?
  (- neighbor(0,0,-1).pid - 1) : (- cell.pid - 1);
       if (faces) for (vector v = *faces, *_i22 = faces; ((scalar *)&v)->i >= 0; v = *++_i22) {



  scalar vt = zn ? VT : VR;

  val(v.z,0,0,0) = _attribute[vt.i].boundary[id](neighbor, point, v.z);
       }
     }
     else

       if (faces) for (vector v = *faces, *_i23 = faces; ((scalar *)&v)->i >= 0; v = *++_i23)
  val(v.z,0,0,0) = 0.;
   }

 }}
    }
  } } end_foreach_boundary_level(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (vectors,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
  enable_fpe_for_mpi();
}



#undef VN
#undef VT
#define VN _attribute[s.i].v.x
#define VT _attribute[s.i].v.y

static double masked_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1449 "/home/vinlinux/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (!(cell.pid < 0) && val(s,0,0,0) != nodata)
      sum += val(s,0,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}


#line 1457

static double masked_average_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/vinlinux/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.x < 0 && (!(cell.pid < 0) || !(neighbor(1,0,0).pid < 0)) &&
 val(s,1,0,0) != nodata)
      sum += val(s,1,0,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}
#line 1457

static double masked_average_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/vinlinux/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.y < 0 && (!(cell.pid < 0) || !(neighbor(0,1,0).pid < 0)) &&
 val(s,0,1,0) != nodata)
      sum += val(s,0,1,0), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}
#line 1457

static double masked_average_z (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 1459 "/home/vinlinux/basilisk/src/grid/tree.h"

  double sum = 0., n = 0.;
   { foreach_child()
    if (child.z < 0 && (!(cell.pid < 0) || !(neighbor(0,0,1).pid < 0)) &&
 val(s,0,0,1) != nodata)
      sum += val(s,0,0,1), n++; end_foreach_child(); }
  return n ? sum/n : nodata;
}

static void masked_boundary_restriction (const Boundary * b,
      scalar * list, int l)
{
  scalar * scalars = NULL;
  vector * faces = NULL;
  if (list) for (scalar s = *list, *_i24 = list; ((scalar *)&s)->i >= 0; s = *++_i24)
    if (!is_constant(s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].v.x.i == s.i && _attribute[s.i].face)
 faces = vectors_add (faces, _attribute[s.i].v);
      else
 scalars = list_add (scalars, s);
    }

   { foreach_halo (restriction, l){

#line 1481 "/home/vinlinux/basilisk/src/grid/tree.h"
 {
    if (scalars) for (scalar s = *scalars, *_i25 = scalars; ((scalar *)&s)->i >= 0; s = *++_i25)
      val(s,0,0,0) = masked_average (parent, s);
    if (faces) for (vector v = *faces, *_i26 = faces; ((scalar *)&v)->i >= 0; v = *++_i26)
      {
#line 1485
 {
 double average = masked_average_x (parent, v.x);
 if ((neighbor(-1,0,0).pid < 0))
   val(v.x,0,0,0) = average;
 if ((neighbor(1,0,0).pid < 0))
   val(v.x,1,0,0) = average;
      }
#line 1485
 {
 double average = masked_average_y (parent, v.y);
 if ((neighbor(0,-1,0).pid < 0))
   val(v.y,0,0,0) = average;
 if ((neighbor(0,1,0).pid < 0))
   val(v.y,0,1,0) = average;
      }
#line 1485
 {
 double average = masked_average_z (parent, v.z);
 if ((neighbor(0,0,-1).pid < 0))
   val(v.z,0,0,0) = average;
 if ((neighbor(0,0,1).pid < 0))
   val(v.z,0,0,1) = average;
      }}
  } } end_foreach_halo(); }

  pfree (scalars,__func__,__FILE__,__LINE__);
  pfree (faces,__func__,__FILE__,__LINE__);
}
#line 1521 "/home/vinlinux/basilisk/src/grid/tree.h"
static void free_cache (CacheLevel * c)
{
  for (int l = 0; l <= depth(); l++)
    pfree (c[l].p,__func__,__FILE__,__LINE__);
  pfree (c,__func__,__FILE__,__LINE__);
}

void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Tree * q = ((Tree *)grid);
  pfree (q->leaves.p,__func__,__FILE__,__LINE__);
  pfree (q->faces.p,__func__,__FILE__,__LINE__);
  pfree (q->vertices.p,__func__,__FILE__,__LINE__);
  pfree (q->refined.p,__func__,__FILE__,__LINE__);


  Layer * L = q->L[0];
#line 1557 "/home/vinlinux/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++) {
    for (int j = Period.y*2; j < L->len - Period.y*2; j++) {
      for (int k = Period.z*2; k < L->len - Period.z*2; k++)
 pfree (L->m[i][j][k],__func__,__FILE__,__LINE__);
      pfree (L->m[i][j],__func__,__FILE__,__LINE__);
    }
    pfree (L->m[i],__func__,__FILE__,__LINE__);
  }

  for (int l = 1; l <= depth(); l++) {
    Layer * L = q->L[l];
    for (int i = Period.x*2; i < L->len - Period.x*2; i++)
      if (L->m[i]) {
 for (int j = Period.y*2; j < L->len - Period.y*2; j++)
   if (L->m[i][j])
     pfree (L->m[i][j],__func__,__FILE__,__LINE__);
 pfree (L->m[i],__func__,__FILE__,__LINE__);
      }
  }

  for (int l = 0; l <= depth(); l++)
    destroy_layer (q->L[l]);
  q->L = &(q->L[-1]);
  pfree (q->L,__func__,__FILE__,__LINE__);
  free_cache (q->active);
  free_cache (q->prolongation);
  free_cache (q->boundary);
  free_cache (q->restriction);
  pfree (q,__func__,__FILE__,__LINE__);
  grid = NULL;
}

static void refine_level (int depth);


void init_grid (int n)
{ trace ("init_grid", "/home/vinlinux/basilisk/src/grid/tree.h", 1593);

  assert (sizeof(Cell) % 8 == 0);

  free_grid();
  int depth = 0;
  while (n > 1) {
    if (n % 2) {
      fprintf (qstderr(), "tree: N must be a power-of-two\n");
      exit (1);
    }
    n /= 2;
    depth++;
  }
  Tree * q = ((Tree *) pcalloc (1, sizeof(Tree),__func__,__FILE__,__LINE__));
  grid = (Grid *) q;
  grid->depth = 0;


  q->L = ((Layer * *) pmalloc ((2)*sizeof(Layer *),__func__,__FILE__,__LINE__));

  q->L[0] = NULL; q->L = &(q->L[1]);

  Layer * L = new_layer (0);
  q->L[0] = L;
#line 1651 "/home/vinlinux/basilisk/src/grid/tree.h"
  for (int i = Period.x*2; i < L->len - Period.x*2; i++)
    for (int j = Period.y*2; j < L->len - Period.y*2; j++) {
      layer_add_row (L, i, j);
      for (int k = Period.z*2; k < L->len - Period.z*2; k++)
 L->m[i][j][k] = (char *) pcalloc (1, sizeof(Cell) + datasize,__func__,__FILE__,__LINE__);
    }
  apply_periodic (q);
  CELL(L->m[2][2][2]).flags |= leaf;
  if (pid() == 0)
    CELL(L->m[2][2][2]).flags |= active;
  for (int k = - 2*(1 - Period.x); k <= 2*(1 - Period.x); k++)
    for (int l = -2*(1 - Period.y); l <= 2*(1 - Period.y); l++)
      for (int n = -2*(1 - Period.z); n <= 2*(1 - Period.z); n++)
 CELL(L->m[2 +k][2 +l][2 +n]).pid = (k > 0 ? -1 - right :
       k < 0 ? -1 - left :
       l > 0 ? -1 - top :
       l < 0 ? -1 - bottom :
       n > 0 ? -1 - front :
       n < 0 ? -1 - back :
       0);
  CELL(L->m[2][2][2]).pid = 0;

  q->active = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->prolongation = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->boundary = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->restriction = ((CacheLevel *) pcalloc (1, sizeof(CacheLevel),__func__,__FILE__,__LINE__));
  q->dirty = true;
  N = 1 << depth;
#if _MPI
  void mpi_boundary_new();
  mpi_boundary_new();
#endif

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  b->restriction = masked_boundary_restriction;
  add_boundary (b);
  refine_level (depth);
  reset (all, 0.);
  { if (((Tree *)grid)->dirty) update_cache_f(); };
 end_trace("init_grid", "/home/vinlinux/basilisk/src/grid/tree.h", 1691); }
#line 1727 "/home/vinlinux/basilisk/src/grid/tree.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  for (int l = depth(); l >= 0; l--) {
    Point point = { .level = l };
    int n = 1 << point.level;
    point.i = (p.x - X0)/L0*n + 2;

    point.j = (p.y - Y0)/L0*n + 2;


    point.k = (p.z - Z0)/L0*n + 2;

    if (point.i >= 0 && point.i < n + 2*2

 && point.j >= 0 && point.j < n + 2*2


 && point.k >= 0 && point.k < n + 2*2

 ) {
      if (allocated(0,0,0) && is_local(cell) && is_leaf(cell))
 return point;
    }
    else
      break;
  }
  Point point = { .level = -1 };
  return point;
}



bool tree_is_full()
{
  { if (((Tree *)grid)->dirty) update_cache_f(); };
  return (grid->tn == 1L << grid->maxdepth*3);
}

#line 1 "grid/tree-common.h"
#line 1 "/home/vinlinux/basilisk/src/grid/tree-common.h"



#line 1 "grid/multigrid-common.h"
#line 1 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/vinlinux/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (qstderr(), "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      assert (parent < 0);
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/vinlinux/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX);
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level));\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
    double Delta_z = Delta;\
\
  double x = (ig/2. + (point.i - 2) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
  double z = (kg/2. + (point.k - 2) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
    NOT_UNUSED(Delta_z);\
\
  ;\

#line 32


#line 1 "grid/fpe.h"
#line 1 "/home/vinlinux/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static void gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', qstderr());
    fflush (qstderr());
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  system (command);
}

static void caught_abort (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (qstderr(), "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 35 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

scalar new_scalar (const char * name)
{
  int nvar = datasize/sizeof(double);
  scalar s;
  for (s.i = 0; s.i < nvar; s.i++)
    if (!list_lookup (all, s)) {
      init_scalar (s, name);
      trash (((scalar []){s, {-1}}));
      all = list_append (all, s);
      return s;
    }


  assert (nvar < _NVARMAX);
  datasize += sizeof(double); nvar++;
  _attribute = (_Attributes *) prealloc (_attribute, (nvar)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar-1], 0, sizeof (_Attributes));
  s = (scalar){nvar - 1};
  realloc_scalar();
  init_scalar (s, name);
  trash (((scalar []){s, {-1}}));
  all = list_append (all, s);
  return s;
}

scalar new_vertex_scalar (const char * name)
{
  scalar s = new_scalar (name);
  {
#line 66

    _attribute[s.i].d.x = -1;
#line 66

    _attribute[s.i].d.y = -1;
#line 66

    _attribute[s.i].d.z = -1;}
  return s;
}

static vector alloc_vector (const char * name)
{
  vector v;
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  {
#line 76
 {
    sprintf (cname, ext.x, name);
    v.x = new_scalar (cname);
  }
#line 76
 {
    sprintf (cname, ext.y, name);
    v.y = new_scalar (cname);
  }
#line 76
 {
    sprintf (cname, ext.z, name);
    v.z = new_scalar (cname);
  }}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_vector (name);
  init_face_vector (v, NULL);
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 102
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 102
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
#line 102
 {
    sprintf (cname, ext.z, name);
    t.z = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 115
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 115
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }
#line 115
 {
    sprintf (cname, ext.z, name);
    t.z.z = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;


    sprintf (cname, "%s.x.z", name);
    t.x.z = new_scalar(cname);
    t.z.x = t.x.z;
    sprintf (cname, "%s.y.z", name);
    t.y.z = new_scalar(cname);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 159

    init_const_scalar (v.x, name, *val++);
#line 159

    init_const_scalar (v.y, name, *val++);
#line 159

    init_const_scalar (v.z, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 166

    v.x.i = _NVARMAX + i++;
#line 166

    v.y.i = _NVARMAX + i++;
#line 166

    v.z.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar) =
    _attribute[a.i].boundary_homogeneous;
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  if (l) for (scalar s = *l, *_i27 = l; ((scalar *)&s)->i >= 0; s = *++_i27) {
    scalar c = new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  if (list) for (scalar s = *list, *_i28 = list; ((scalar *)&s)->i >= 0; s = *++_i28)
    {
#line 201

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 201

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
#line 201

      if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  if (list) for (scalar f = *list, *_i29 = list; ((scalar *)&f)->i >= 0; f = *++_i29) {
    if (_attribute[f.i].delete)
      _attribute[f.i].delete (f);
    pfree (_attribute[f.i].name,__func__,__FILE__,__LINE__); _attribute[f.i].name = NULL;
    pfree (_attribute[f.i].boundary,__func__,__FILE__,__LINE__); _attribute[f.i].boundary = NULL;
    pfree (_attribute[f.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[f.i].boundary_homogeneous = NULL;
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  if (list) for (scalar f = *list, *_i30 = list; ((scalar *)&f)->i >= 0; f = *++_i30) {
    scalar * s = all;
    for (; s->i >= 0 && s->i != f.i; s++);
    if (s->i == f.i)
      for (; s->i >= 0; s++)
 s[0] = s[1];
  }
}

typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}

void free_solver()
{
  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if 0
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_flux) (vector *);


void boundary (scalar * list)
{ trace ("boundary", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 289);
  if (list == NULL)
    { ; end_trace("boundary", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 291);  return; }
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i31 = list; ((scalar *)&s)->i >= 0; s = *++_i31)
    if (!is_constant(s) && _attribute[s.i].face)
      listf = vectors_add (listf, _attribute[s.i].v);
  if (listf) {
    boundary_flux (listf);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
  boundary_level (list, -1);
 end_trace("boundary", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 301); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_flux (vector * list)
{

}

static double symmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 314 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 319 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  pfree (_attribute[s.i].boundary,__func__,__FILE__,__LINE__);
  pfree (_attribute[s.i].boundary_homogeneous,__func__,__FILE__,__LINE__);

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].name = pname;

  _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
  {
#line 351
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 351
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
#line 351
 {
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }
#line 368
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_scalar (v.z, cname);
    }
    else
      init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 388
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 388
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
#line 388
 {
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }
#line 400
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_vector (t.z, cname);
    }
    else
      init_vector (t.z, NULL);
  }}
#line 424 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
    assert (false);

  return t;
}

void output_cells (FILE * fp)
{
   { foreach(){

#line 431 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
 {
    Delta /= 2.;
#line 443 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

  } } end_foreach(); }
  fflush (fp);
}

static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 494 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells (fp);
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  if (all) for (scalar v = *all, *_i32 = all; ((scalar *)&v)->i >= 0; v = *++_i32)





    fprintf (fp, "x y z %s ", _attribute[v.i].name);

  fputc ('\n', fp);
#line 541 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   if (all) for (scalar v = *all, *_i33 = all; ((scalar *)&v)->i >= 0; v = *++_i33) {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  if (all) for (scalar s = *all, *_i34 = all; ((scalar *)&s)->i >= 0; s = *++_i34) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (qstderr(), "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (qstderr(), _attribute[0].name, name, stencil);
  fflush (qstderr());

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_flux = cartesian_boundary_flux;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 597 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;
#line 614 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"
  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  z = (p.z - z)/Delta - _attribute[v.i].d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 629);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 630 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 632);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 633);  return _ret; }
 end_trace("interpolate", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 634); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 638);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 641 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      if (list) for (scalar s = *list, *_i35 = list; ((scalar *)&s)->i >= 0; s = *++_i35)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      if (list) for (scalar s = *list, *_i36 = list; ((scalar *)&s)->i >= 0; s = *++_i36)
 v[j++] = nodata;
  }
#if _MPI
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/vinlinux/basilisk/src/grid/cartesian-common.h", 660); }



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  if (all) for (scalar s = *all, *_i37 = all; ((scalar *)&s)->i >= 0; s = *++_i37) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  if (all) for (scalar s = *all, *_i38 = all; ((scalar *)&s)->i >= 0; s = *++_i38) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 680

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 680

 _attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
#line 680

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 692 "/home/vinlinux/basilisk/src/grid/cartesian-common.h"

  return nodata;
}

static void periodic_boundary (int d)
{

  if (all) for (scalar s = *all, *_i39 = all; ((scalar *)&s)->i >= 0; s = *++_i39)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  if (all) for (scalar s = *all, *_i40 = all; ((scalar *)&s)->i >= 0; s = *++_i40)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{





    assert (dir <= back);


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}
#line 4 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 27 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3);
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 35 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/val_cm(cm,0,0,0);
 }}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 43 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  }
#line 44
 {







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  }
#line 44
 {







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }}
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 61 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 64 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 67 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--)
     { foreach_coarse_level (l){

#line 76 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
 {
      double sc[1 << 3];
      int c = 0;
       { foreach_child()
 sc[c++] = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
      c = 0;
       { foreach_child() {

 val(w,0,0,0) = sc[c] - val(s,0,0,0);
 val(s,0,0,0) = sc[c++];
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }

   { foreach_level(0){

#line 90 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
 val(w,0,0,0) = 0.; } end_foreach_level(); }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 94 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"








    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 112 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 123 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

#line 138 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
  assert (false);
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 144 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"







  assert (false);
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 157 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 163 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 166

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 169

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 175

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 166

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 166

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 166

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 169

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 169

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 169

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 175

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 175

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 175

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 183 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 189 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 206

    _attribute[v.y.i].restriction = no_restriction;
#line 206

    _attribute[v.z.i].restriction = no_restriction;
#line 206

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 213 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 246 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     if (all) for (scalar v = *all, *_i41 = all; ((scalar *)&v)->i >= 0; v = *++_i41)
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 303 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     if (all) for (scalar v = *all, *_i42 = all; ((scalar *)&v)->i >= 0; v = *++_i42) {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }
     fputc ('\n', fp);
   }
      fprintf (qstderr(), ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (qstderr());
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i43 = list; ((scalar *)&s)->i >= 0; s = *++_i43)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 342

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 342

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 342

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 351 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i44 = listdef; ((scalar *)&s)->i >= 0; s = *++_i44)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i45 = listc; ((scalar *)&s)->i >= 0; s = *++_i45)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { foreach(){

#line 386 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 395 "/home/vinlinux/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 5 "/home/vinlinux/basilisk/src/grid/tree-common.h"






#line 21 "/home/vinlinux/basilisk/src/grid/tree-common.h"
int refine_cell (Point point, scalar * list, int flag, Cache * refined)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 22 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  int nr = 0;


  if (level > 0)
    for (int k = 0; k != 2*child.x; k += child.x)

      for (int l = 0; l != 2*child.y; l += child.y)


 for (int m = 0; m != 2*child.z; m += child.z)

   if (aparent(k,l,m).pid >= 0 && is_leaf(aparent(k,l,m))) {
     Point p = point;


     p.level = point.level - 1;
     p.i = (point.i + 2)/2 + k;
     do { if (p.i < 2) p.i += 1 << p.level; else if (p.i >= 2 + (1 << p.level)) p.i -= 1 << p.level; } while(0);

       p.j = (point.j + 2)/2 + l;
       do { if (p.j < 2) p.j += 1 << p.level; else if (p.j >= 2 + (1 << p.level)) p.j -= 1 << p.level; } while(0);


       p.k = (point.k + 2)/2 + m;
       do { if (p.k < 2) p.k += 1 << p.level; else if (p.k >= 2 + (1 << p.level)) p.k -= 1 << p.level; } while(0);

     nr += refine_cell (p, list, flag, refined);
     aparent(k,l,m).flags |= flag;
   }



  cell.flags &= ~leaf;


  increment_neighbors (point);

  int cflag = is_active(cell) ? (active|leaf) : leaf;
   { foreach_child()
    cell.flags |= cflag; end_foreach_child(); }


  if (list) for (scalar s = *list, *_i46 = list; ((scalar *)&s)->i >= 0; s = *++_i46)
    if (is_local(cell) || _attribute[s.i].face)
      _attribute[s.i].refine (point, s);

#if _MPI
  if (is_border(cell)) {
     { foreach_child() {
      bool bord = false;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   bord = true, foreach_neighbor_break();
 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       bord = true, foreach_child_break(); end_foreach_child(); }
 if (bord)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      if (bord)
 cell.flags |= border;
    } end_foreach_child(); }
    if (refined)
      cache_append (refined, point, cell.flags);
    nr++;
  }
#endif
  return nr;
}





bool coarsen_cell (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 99 "/home/vinlinux/basilisk/src/grid/tree-common.h"




  int pid = cell.pid;
   { foreach_child()
    if (cell.neighbors || (cell.pid < 0 && cell.pid != pid))
      return false; end_foreach_child(); }



  if (list) for (scalar s = *list, *_i47 = list; ((scalar *)&s)->i >= 0; s = *++_i47) {
    _attribute[s.i].restriction (point, s);
    if (_attribute[s.i].coarsen)
      _attribute[s.i].coarsen (point, s);
  }


  cell.flags |= leaf;


  decrement_neighbors (point);

#if _MPI
  if (!is_local(cell)) {
    cell.flags &= ~(active|border);
     { foreach_neighbor(1)
      if (cell.neighbors)
  { foreach_child()
   if (allocated(0,0,0) && is_local(cell) && !is_border(cell))
     cell.flags |= border; end_foreach_child(); } end_foreach_neighbor(); }
  }
#endif

  return true;
}

void coarsen_cell_recursive (Point point, scalar * list)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 137 "/home/vinlinux/basilisk/src/grid/tree-common.h"



   { foreach_child()
    if (cell.neighbors)
       { foreach_neighbor(1)
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
   coarsen_cell_recursive (point, list); end_foreach_neighbor(); } end_foreach_child(); }

  assert (coarsen_cell (point, list));
}

void mpi_boundary_refine (scalar *);
void mpi_boundary_coarsen (int, int);
void mpi_boundary_update (scalar *);

typedef struct {
  int nc, nf;
} astats;

struct Adapt {
  scalar * slist;
  double * max;
  int maxlevel;
  int minlevel;
  scalar * list;
};


astats adapt_wavelet (struct Adapt p)
{ trace ("adapt_wavelet", "/home/vinlinux/basilisk/src/grid/tree-common.h", 167);
  if (p.list == NULL)
    p.list = all;
  if (is_constant(cm))
    restriction (p.slist);
  else {
    scalar * listr = list_concat (((scalar []){cm,{-1}}), p.slist);
    restriction (listr);
    pfree (listr,__func__,__FILE__,__LINE__);
  }

  astats st = {0, 0};
  scalar * listc = NULL;
  if (p.list) for (scalar s = *p.list, *_i48 = p.list; ((scalar *)&s)->i >= 0; s = *++_i48)
    if (!is_constant(s) && _attribute[s.i].restriction != no_restriction)
      listc = list_add (listc, s);


  if (p.minlevel < 1)
    p.minlevel = 1;
  ((Tree *)grid)->refined.n = 0;
  static const int refined = 1 << user, too_fine = 1 << (user + 1);
   { foreach_cell(){

#line 189 "/home/vinlinux/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell)) {
      static const int too_coarse = 1 << (user + 2);
      if (is_leaf (cell)) {
 if (cell.flags & too_coarse) {
   cell.flags &= ~too_coarse;
   refine_cell (point, listc, refined, &((Tree *)grid)->refined);
   st.nf++;
 }
 continue;
      }
      else {
 if (cell.flags & refined) {

   cell.flags &= ~too_coarse;
   continue;
 }

 bool local = is_local(cell);
 if (!local)
    { foreach_child()
     if (is_local(cell))
       local = true, foreach_child_break(); end_foreach_child(); }
 if (local) {
   int i = 0;
   static const int just_fine = 1 << (user + 3);
   if (p.slist) for (scalar s = *p.slist, *_i49 = p.slist; ((scalar *)&s)->i >= 0; s = *++_i49) {
     double max = p.max[i++], sc[1 << 3];
     int c = 0;
      { foreach_child()
       sc[c++] = val(s,0,0,0); end_foreach_child(); }
     _attribute[s.i].prolongation (point, s);
     c = 0;
      { foreach_child() {
       double e = fabs(sc[c] - val(s,0,0,0));
       if (e > max && level < p.maxlevel) {
  cell.flags &= ~too_fine;
  cell.flags |= too_coarse;
       }
       else if ((e <= max/1.5 || level > p.maxlevel) &&
         !(cell.flags & (too_coarse|just_fine))) {
  if (level >= p.minlevel)
    cell.flags |= too_fine;
       }
       else if (!(cell.flags & too_coarse)) {
  cell.flags &= ~too_fine;
  cell.flags |= just_fine;
       }
       val(s,0,0,0) = sc[c++];
     } end_foreach_child(); }
   }
    { foreach_child() {
     cell.flags &= ~just_fine;
     if (!is_leaf(cell)) {
       cell.flags &= ~too_coarse;
       if (level >= p.maxlevel)
  cell.flags |= too_fine;
     }
     else if (!is_active(cell))
       cell.flags &= ~too_coarse;
   } end_foreach_child(); }
 }
      }
    }
    else
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_refine (listc);



  for (int l = depth(); l >= 0; l--) {
     { foreach_cell(){

#line 261 "/home/vinlinux/basilisk/src/grid/tree-common.h"

      if (!(cell.pid < 0)) {
 if (level == l) {
   if (!is_leaf(cell)) {
     if (cell.flags & refined)

       cell.flags &= ~(refined|too_fine);
     else if (cell.flags & too_fine) {
       if (is_local(cell) && coarsen_cell (point, listc))
  st.nc++;
       cell.flags &= ~too_fine;
     }
   }
   if (cell.flags & too_fine)
     cell.flags &= ~too_fine;
   else if (level > 0 && (aparent(0,0,0).flags & too_fine))
     aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
 else if (is_leaf(cell))
   continue;
      } } end_foreach_cell(); }
    mpi_boundary_coarsen (l, too_fine);
  }
  pfree (listc,__func__,__FILE__,__LINE__);

  mpi_all_reduce (st.nf, MPI_INT, MPI_SUM);
  mpi_all_reduce (st.nc, MPI_INT, MPI_SUM);
  if (st.nc || st.nf)
    mpi_boundary_update (p.list);

  { astats _ret =  st; end_trace("adapt_wavelet", "/home/vinlinux/basilisk/src/grid/tree-common.h", 292);  return _ret; }
 end_trace("adapt_wavelet", "/home/vinlinux/basilisk/src/grid/tree-common.h", 293); }
#line 314 "/home/vinlinux/basilisk/src/grid/tree-common.h"
static void refine_level (int depth)
{
  int refined;
  do {
    refined = 0;
    ((Tree *)grid)->refined.n = 0;
     { foreach_leaf(){

#line 320 "/home/vinlinux/basilisk/src/grid/tree-common.h"

      if (level < depth) {
 refine_cell (point, NULL, 0, &((Tree *)grid)->refined);
 refined++;
 continue;
      } } end_foreach_leaf(); }
    mpi_all_reduce (refined, MPI_INT, MPI_SUM);
    if (refined) {
      mpi_boundary_refine (NULL);
      mpi_boundary_update (NULL);
    }
  } while (refined);
}
#line 359 "/home/vinlinux/basilisk/src/grid/tree-common.h"
static void halo_flux (vector * list)
{
  vector * listv = NULL;
  if (list) for (vector v = *list, *_i50 = list; ((scalar *)&v)->i >= 0; v = *++_i50)
    if (!is_constant(v.x))
      listv = vectors_add (listv, v);

  if (listv) {
    for (int l = depth() - 1; l >= 0; l--)
       { foreach_halo (prolongation, l){

#line 368 "/home/vinlinux/basilisk/src/grid/tree-common.h"

 {
#line 369
 {
#line 385 "/home/vinlinux/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.x,0,0,0) = (fine(f.x,0,0,0) + fine(f.x,0,1,0) +
         fine(f.x,0,0,1) + fine(f.x,0,1,1))/4.;
   if ((!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i52 = listv; ((scalar *)&f)->i >= 0; f = *++_i52)
       val(f.x,1,0,0) = (fine(f.x,2,0,0) + fine(f.x,2,1,0) +
   fine(f.x,2,0,1) + fine(f.x,2,1,1))/4.;

      }
#line 369
 {
#line 385 "/home/vinlinux/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.y,0,0,0) = (fine(f.y,0,0,0) + fine(f.y,0,0,1) +
         fine(f.y,1,0,0) + fine(f.y,1,0,1))/4.;
   if ((!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0))
     if (listv) for (vector f = *listv, *_i52 = listv; ((scalar *)&f)->i >= 0; f = *++_i52)
       val(f.y,0,1,0) = (fine(f.y,0,2,0) + fine(f.y,0,2,1) +
   fine(f.y,1,2,0) + fine(f.y,1,2,1))/4.;

      }
#line 369
 {
#line 385 "/home/vinlinux/basilisk/src/grid/tree-common.h"
   if ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0))
     if (listv) for (vector f = *listv, *_i51 = listv; ((scalar *)&f)->i >= 0; f = *++_i51)
       val(f.z,0,0,0) = (fine(f.z,0,0,0) + fine(f.z,1,0,0) +
         fine(f.z,0,1,0) + fine(f.z,1,1,0))/4.;
   if ((!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0))
     if (listv) for (vector f = *listv, *_i52 = listv; ((scalar *)&f)->i >= 0; f = *++_i52)
       val(f.z,0,0,1) = (fine(f.z,0,0,2) + fine(f.z,1,0,2) +
   fine(f.z,0,1,2) + fine(f.z,1,1,2))/4.;

      }} } end_foreach_halo(); }
    pfree (listv,__func__,__FILE__,__LINE__);
  }
}



static scalar tree_init_scalar (scalar s, const char * name)
{
  s = multigrid_init_scalar (s, name);
  _attribute[s.i].refine = _attribute[s.i].prolongation;
  return s;
}


#line 408

static void refine_face_x (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 410 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(-1,0,0)))) {
    double g1 = (val(v.x,0,+1,0) - val(v.x,0,-1,0))/8.;
    double g2 = (val(v.x,0,0,+1) - val(v.x,0,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,0,j,k) = val(v.x,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0) && neighbor(1,0,0).neighbors &&
      (is_local(cell) || is_local(neighbor(1,0,0)))) {
    double g1 = (val(v.x,1,+1,0) - val(v.x,1,-1,0))/8.;
    double g2 = (val(v.x,1,0,+1) - val(v.x,1,0,-1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,2,j,k) = val(v.x,1,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.x,0,+1,0) + val(v.x,1,+1,0) - val(v.x,0,-1,0) - val(v.x,1,-1,0))/16.;
    double g2 = (val(v.x,0,0,+1) + val(v.x,1,0,+1) - val(v.x,0,0,-1) - val(v.x,1,0,-1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.x,1,j,k) = (val(v.x,0,0,0) + val(v.x,1,0,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 408

static void refine_face_y (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 410 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,-1,0)))) {
    double g1 = (val(v.y,0,0,+1) - val(v.y,0,0,-1))/8.;
    double g2 = (val(v.y,+1,0,0) - val(v.y,-1,0,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,0,j) = val(v.y,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0) && neighbor(0,1,0).neighbors &&
      (is_local(cell) || is_local(neighbor(0,1,0)))) {
    double g1 = (val(v.y,0,1,+1) - val(v.y,0,1,-1))/8.;
    double g2 = (val(v.y,+1,1,0) - val(v.y,-1,1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,2,j) = val(v.y,0,1,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.y,0,0,+1) + val(v.y,0,1,+1) - val(v.y,0,0,-1) - val(v.y,0,1,-1))/16.;
    double g2 = (val(v.y,+1,0,0) + val(v.y,+1,1,0) - val(v.y,-1,0,0) - val(v.y,-1,1,0))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.y,k,1,j) = (val(v.y,0,0,0) + val(v.y,0,1,0))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}
#line 408

static void refine_face_z (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 410 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  if (!(!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) &&
      (is_local(cell) || is_local(neighbor(0,0,-1)))) {
    double g1 = (val(v.z,+1,0,0) - val(v.z,-1,0,0))/8.;
    double g2 = (val(v.z,0,+1,0) - val(v.z,0,-1,0))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,0) = val(v.z,0,0,0) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (!(!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0) && neighbor(0,0,1).neighbors &&
      (is_local(cell) || is_local(neighbor(0,0,1)))) {
    double g1 = (val(v.z,+1,0,1) - val(v.z,-1,0,1))/8.;
    double g2 = (val(v.z,0,+1,1) - val(v.z,0,-1,1))/8.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,2) = val(v.z,0,0,1) + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
  if (is_local(cell)) {
    double g1 = (val(v.z,+1,0,0) + val(v.z,+1,0,1) - val(v.z,-1,0,0) - val(v.z,-1,0,1))/16.;
    double g2 = (val(v.z,0,+1,0) + val(v.z,0,+1,1) - val(v.z,0,-1,0) - val(v.z,0,-1,1))/16.;
    for (int j = 0; j <= 1; j++)
      for (int k = 0; k <= 1; k++)
 fine(v.z,j,k,1) = (val(v.z,0,0,0) + val(v.z,0,0,1))/2. + (2*j - 1)*g1 + (2*k - 1)*g2;
  }
}

void refine_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 438 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  vector v = _attribute[s.i].v;
  {
#line 440

    _attribute[v.x.i].prolongation (point, v.x);
#line 440

    _attribute[v.y.i].prolongation (point, v.y);
#line 440

    _attribute[v.z.i].prolongation (point, v.z);}
}

void refine_face_solenoidal (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 445 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  refine_face (point, s);

  if (is_local(cell)) {

    vector v = _attribute[s.i].v;
    double d[1 << 3], p[1 << 3];
    int i = 0;
     { foreach_child() {
      d[i] = 0.;
      {
#line 455

 d[i] += val(v.x,1,0,0) - val(v.x,0,0,0);
#line 455

 d[i] += val(v.y,0,1,0) - val(v.y,0,0,0);
#line 455

 d[i] += val(v.z,0,0,1) - val(v.z,0,0,0);}
      i++;
    } end_foreach_child(); }
#line 469 "/home/vinlinux/basilisk/src/grid/tree-common.h"
    static double m[7][7] = {
      {7./12,5./24,3./8,5./24,3./8,1./4,1./3},
      {5./24,7./12,3./8,5./24,1./4,3./8,1./3},
      {3./8,3./8,3./4,1./4,3./8,3./8,1./2},
      {5./24,5./24,1./4,7./12,3./8,3./8,1./3},
      {3./8,1./4,3./8,3./8,3./4,3./8,1./2},
      {1./4,3./8,3./8,3./8,3./8,3./4,1./2},
      {1./3,1./3,1./2,1./3,1./2,1./2,5./6}
    };
    p[0] = 0.;
    for (int i = 0; i < 7; i++) {
      p[i + 1] = 0.;
      for (int j = 0; j < 7; j++)
 p[i + 1] += m[i][j]*d[j+1];
    }
    for (int k = 0; k <= 1; k++) {
      fine(v.x,1,0,k) += p[4+k] - p[0+k];
      fine(v.x,1,1,k) += p[6+k] - p[2+k];
      fine(v.y,0,1,k) += p[2+k] - p[0+k];
      fine(v.y,1,1,k) += p[6+k] - p[4+k];
    }
    fine(v.z,0,0,1) += p[1] - p[0];
    fine(v.z,0,1,1) += p[3] - p[2];
    fine(v.z,1,0,1) += p[5] - p[4];
    fine(v.z,1,1,1) += p[7] - p[6];

  }

}

vector tree_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 502

    _attribute[v.x.i].restriction = _attribute[v.x.i].refine = no_restriction;
#line 502

    _attribute[v.y.i].restriction = _attribute[v.y.i].refine = no_restriction;
#line 502

    _attribute[v.z.i].restriction = _attribute[v.z.i].refine = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  _attribute[v.x.i].refine = refine_face;
  {
#line 506

    _attribute[v.x.i].prolongation = refine_face_x;
#line 506

    _attribute[v.y.i].prolongation = refine_face_y;
#line 506

    _attribute[v.z.i].prolongation = refine_face_z;}
  return v;
}

static void tree_boundary_level (scalar * list, int l)
{
  int depth = l < 0 ? depth() : l;

  if (tree_is_full()) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, depth); };
    return;
  }

  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i53 = list; ((scalar *)&s)->i >= 0; s = *++_i53)
    if (!is_constant (s)) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 530

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 530

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 530

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, depth); };
    for (int l = depth - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 540 "/home/vinlinux/basilisk/src/grid/tree-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i54 = listdef; ((scalar *)&s)->i >= 0; s = *++_i54)
   restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i55 = listc; ((scalar *)&s)->i >= 0; s = *++_i55)
   _attribute[s.i].restriction (point, s);
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }

  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i56 = list; ((scalar *)&s)->i >= 0; s = *++_i56)
    if (!is_constant (s) && _attribute[s.i].refine != no_restriction) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }

  if (listr || listf) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, 0); };
    for (int i = 0; i < depth; i++) {
       { foreach_halo (prolongation, i){

#line 566 "/home/vinlinux/basilisk/src/grid/tree-common.h"
 {
 if (listr) for (scalar s = *listr, *_i57 = listr; ((scalar *)&s)->i >= 0; s = *++_i57)
          _attribute[s.i].prolongation (point, s);
 if (listf) for (vector v = *listf, *_i58 = listf; ((scalar *)&v)->i >= 0; v = *++_i58)
   {
#line 570

     _attribute[v.x.i].prolongation (point, v.x);
#line 570

     _attribute[v.y.i].prolongation (point, v.y);
#line 570

     _attribute[v.z.i].prolongation (point, v.z);}
      } } end_foreach_halo(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, i + 1); };
    }
    pfree (listr,__func__,__FILE__,__LINE__);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
}

double treex (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 580 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;




  assert (false);
  double i = 0;

  return treex(parent) + i/(1 << 2*(level - 1));
}

double treey (Point point) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 593 "/home/vinlinux/basilisk/src/grid/tree-common.h"

  if (level == 0)
    return 0;
  return treey(parent) + 4./(1 << 2*(level - 1));
}

void output_tree (FILE * fp)
{
   { foreach_cell(){

#line 601 "/home/vinlinux/basilisk/src/grid/tree-common.h"

    if (cell.neighbors)
       { foreach_child()
 if (is_local(cell))
   fprintf (fp, "%g %g\n%g %g\n\n",
     treex(parent), treey(parent), treex(point), treey(point)); end_foreach_child(); }; } end_foreach_cell(); }
}

void tree_check()
{


  long nleaves = 0, nactive = 0;
   { foreach_cell_all(){

#line 614 "/home/vinlinux/basilisk/src/grid/tree-common.h"
 {
    if (is_leaf(cell)) {
      assert (cell.pid >= 0);
      nleaves++;
    }
    if (is_local(cell))
      assert (is_active(cell) || (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0));
    if (is_active(cell))
      nactive++;

    int neighbors = 0;
     { foreach_neighbor(1)
      if (allocated(0,0,0) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 neighbors++; end_foreach_neighbor(); }
    assert (cell.neighbors == neighbors);


    if (!cell.neighbors)
      assert (!allocated_child(0,0,0));
  } } end_foreach_cell_all(); }


  long reachable = 0;
   { foreach_cell(){

#line 637 "/home/vinlinux/basilisk/src/grid/tree-common.h"
 {
    if (is_active(cell))
      reachable++;
    else
      continue;
  } } end_foreach_cell(); }
  assert (nactive == reachable);


  reachable = 0;
   { foreach_cell(){

#line 647 "/home/vinlinux/basilisk/src/grid/tree-common.h"

    if (is_leaf(cell)) {
      reachable++;
      continue;
    } } end_foreach_cell(); }
  assert (nleaves == reachable);
}

static void tree_restriction (scalar * list) {
  if (tree_is_full())
    multigrid_restriction (list);

}

void tree_methods()
{
  multigrid_methods();
  init_scalar = tree_init_scalar;
  init_face_vector = tree_init_face_vector;
  boundary_level = tree_boundary_level;
  boundary_flux = halo_flux;
  restriction = tree_restriction;
}
#line 1768 "/home/vinlinux/basilisk/src/grid/tree.h"


void tree_periodic (int dir)
{
  int depth = grid ? depth() : -1;
  if (grid)
    free_grid();
  periodic (dir);
  if (depth >= 0)
    init_grid (1 << depth);
}


#if _MPI
#line 1 "grid/tree-mpi.h"
#line 1 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

int debug_iteration = -1;

void debug_mpi (FILE * fp1);

typedef struct {
  CacheLevel * halo;
  void * buf;
  MPI_Request r;
  int depth;
  int pid;
  int maxdepth;
} Rcv;

typedef struct {
  Rcv * rcv;
  char * name;
  int npid;
} RcvPid;

typedef struct {
  RcvPid * rcv, * snd;
} SndRcv;

typedef struct {
  Boundary parent;

  SndRcv mpi_level, mpi_level_root, restriction;
  Array * send, * receive;
} MpiBoundary;

static void cache_level_init (CacheLevel * c)
{
  c->p = NULL;
  c->n = c->nm = 0;
}

static void rcv_append (Point point, Rcv * rcv)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 39 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

  if (level > rcv->depth) {
    rcv->halo = (CacheLevel *) prealloc (rcv->halo, (level + 1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__);
    for (int j = rcv->depth + 1; j <= level; j++)
      cache_level_init (&rcv->halo[j]);
    rcv->depth = level;
  }
  cache_level_append (&rcv->halo[level], point);
  if (level > rcv->maxdepth)
    rcv->maxdepth = level;
}

void rcv_print (Rcv * rcv, FILE * fp, const char * prefix)
{
  for (int l = 0; l <= rcv->depth; l++)
    if (rcv->halo[l].n > 0)
       { foreach_cache_level(rcv->halo[l], l){

#line 55 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

 fprintf (fp, "%s%g %g %g %d %d\n", prefix, x, y, z, rcv->pid, level); } end_foreach_cache_level(); }
}

static void rcv_free_buf (Rcv * rcv)
{
  if (rcv->buf) {
    prof_start ("rcv_pid_receive");
    MPI_Wait (&rcv->r, MPI_STATUS_IGNORE);
    pfree (rcv->buf,__func__,__FILE__,__LINE__);
    rcv->buf = NULL;
    prof_stop();
  }
}

static void rcv_destroy (Rcv * rcv)
{
  rcv_free_buf (rcv);
  for (int i = 0; i <= rcv->depth; i++)
    if (rcv->halo[i].n > 0)
      pfree (rcv->halo[i].p,__func__,__FILE__,__LINE__);
  pfree (rcv->halo,__func__,__FILE__,__LINE__);
}

static RcvPid * rcv_pid_new (const char * name)
{
  RcvPid * r = ((RcvPid *) pcalloc (1, sizeof(RcvPid),__func__,__FILE__,__LINE__));
  r->name = pstrdup (name,__func__,__FILE__,__LINE__);
  return r;
}

static Rcv * rcv_pid_pointer (RcvPid * p, int pid)
{
  assert (pid >= 0 && pid < npe());

  int i;
  for (i = 0; i < p->npid; i++)
    if (pid == p->rcv[i].pid)
      break;

  if (i == p->npid) {
    p->rcv = (Rcv *) prealloc (p->rcv, (++p->npid)*sizeof(Rcv),__func__,__FILE__,__LINE__);
    Rcv * rcv = &p->rcv[p->npid-1];
    rcv->pid = pid;
    rcv->depth = rcv->maxdepth = 0;
    rcv->halo = ((CacheLevel *) pmalloc ((1)*sizeof(CacheLevel),__func__,__FILE__,__LINE__));
    rcv->buf = NULL;
    cache_level_init (&rcv->halo[0]);
  }
  return &p->rcv[i];
}

static void rcv_pid_append (RcvPid * p, int pid, Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 108 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

  rcv_append (point, rcv_pid_pointer (p, pid));
}

static void rcv_pid_append_pids (RcvPid * p, Array * pids)
{

  for (int i = 0; i < p->npid; i++) {
    int pid = p->rcv[i].pid, j, * a;
    for (j = 0, a = pids->p; j < pids->len/sizeof(int); j++,a++)
      if (*a == pid)
 break;
    if (j == pids->len/sizeof(int))
      array_append (pids, &pid, sizeof(int));
  }
}

void rcv_pid_write (RcvPid * p, const char * name)
{
  for (int i = 0; i < p->npid; i++) {
    Rcv * rcv = &p->rcv[i];
    char fname[80];
    sprintf (fname, "%s-%d-%d", name, pid(), rcv->pid);
    FILE * fp = fopen (fname, "w");
    rcv_print (rcv, fp, "");
    fclose (fp);
  }
}

static void rcv_pid_print (RcvPid * p, FILE * fp, const char * prefix)
{
  for (int i = 0; i < p->npid; i++)
    rcv_print (&p->rcv[i], fp, prefix);
}

static void rcv_pid_destroy (RcvPid * p)
{
  for (int i = 0; i < p->npid; i++)
    rcv_destroy (&p->rcv[i]);
  pfree (p->rcv,__func__,__FILE__,__LINE__);
  pfree (p->name,__func__,__FILE__,__LINE__);
  pfree (p,__func__,__FILE__,__LINE__);
}

static Boundary * mpi_boundary = NULL;






void debug_mpi (FILE * fp1);

static void apply_bc (Rcv * rcv, scalar * list, vector * listf, int l,
        MPI_Status s)
{
  double * b = rcv->buf;
   { foreach_cache_level(rcv->halo[l], l){

#line 165 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (list) for (scalar s = *list, *_i59 = list; ((scalar *)&s)->i >= 0; s = *++_i59)
      val(s,0,0,0) = *b++;
    if (listf) for (vector v = *listf, *_i60 = listf; ((scalar *)&v)->i >= 0; v = *++_i60) {
      {
#line 169
 {
 val(v.x,0,0,0) = *b++;
 if (allocated(1,0,0))
   val(v.x,1,0,0) = *b;
 b++;
      }
#line 169
 {
 val(v.y,0,0,0) = *b++;
 if (allocated(0,1,0))
   val(v.y,0,1,0) = *b;
 b++;
      }
#line 169
 {
 val(v.z,0,0,0) = *b++;
 if (allocated(0,0,1))
   val(v.z,0,0,1) = *b;
 b++;
      }}
    }
  } } end_foreach_cache_level(); }
  size_t size = b - (double *) rcv->buf;
  pfree (rcv->buf,__func__,__FILE__,__LINE__);
  rcv->buf = NULL;

  int rlen;
  MPI_Get_count (&s, MPI_DOUBLE, &rlen);
  if (rlen != size) {
    fprintf (qstderr(),
      "rlen (%d) != size (%ld), %d receiving from %d at level %d\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      rlen, size, pid(), rcv->pid, l);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -2);
  }
}
#line 215 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
static void mpi_recv_check (void * buf, int count, MPI_Datatype datatype,
       int source, int tag,
       MPI_Comm comm, MPI_Status * status,
       const char * name)
{
#line 250 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
  int errorcode = MPI_Recv (buf, count, datatype, source, tag, comm, status);
  if (errorcode != MPI_SUCCESS) {
    char string[MPI_MAX_ERROR_STRING];
    int resultlen;
    MPI_Error_string (errorcode, string, &resultlen);
    fprintf (qstderr(),
      "ERROR MPI_Recv \"%s\" (count = %d, source = %d, tag = %d):\n%s\n"
      "Calling debug_mpi(NULL)...\n"
      "Aborting...\n",
      name, count, source, tag, string);
    fflush (qstderr());
    debug_mpi (NULL);
    MPI_Abort (MPI_COMM_WORLD, -1);
  }





}


static int mpi_waitany (int count, MPI_Request array_of_requests[], int *indx,
   MPI_Status *status)
{ trace ("mpi_waitany", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 274);
  { int _ret =  MPI_Waitany (count, array_of_requests, indx, status); end_trace("mpi_waitany", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 275);  return _ret; }
 end_trace("mpi_waitany", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 276); }

static void rcv_pid_receive (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_receive");

  int len = list_len (list) + 2*3*vectors_len (listf);

  MPI_Request r[m->npid];
  Rcv * rrcv[m->npid];
  int nr = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);






      MPI_Irecv (rcv->buf, rcv->halo[l].n*len, MPI_DOUBLE, rcv->pid,
   (l), MPI_COMM_WORLD, &r[nr]);
      rrcv[nr++] = rcv;






    }
  }


  if (nr > 0) {
    int i;
    MPI_Status s;
    mpi_waitany (nr, r, &i, &s);
    while (i != MPI_UNDEFINED) {
      Rcv * rcv = rrcv[i];
      assert (l <= rcv->depth && rcv->halo[l].n > 0);
      assert (rcv->buf);
      apply_bc (rcv, list, listf, l, s);
      mpi_waitany (nr, r, &i, &s);
    }
  }

  prof_stop();
}


static void rcv_pid_wait (RcvPid * m)
{ trace ("rcv_pid_wait", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 332);

  for (int i = 0; i < m->npid; i++)
    rcv_free_buf (&m->rcv[i]);
 end_trace("rcv_pid_wait", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 336); }

static void rcv_pid_send (RcvPid * m, scalar * list, vector * listf, int l)
{
  if (m->npid == 0)
    return;

  prof_start ("rcv_pid_send");

  int len = list_len (list) + 2*3*vectors_len (listf);


  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0) {
      assert (!rcv->buf);
      rcv->buf = pmalloc (sizeof (double)*rcv->halo[l].n*len,__func__,__FILE__,__LINE__);
      double * b = rcv->buf;
       { foreach_cache_level(rcv->halo[l], l){

#line 354 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
 if (list) for (scalar s = *list, *_i61 = list; ((scalar *)&s)->i >= 0; s = *++_i61)
   *b++ = val(s,0,0,0);
 if (listf) for (vector v = *listf, *_i62 = listf; ((scalar *)&v)->i >= 0; v = *++_i62)
   {
#line 358
 {
     *b++ = val(v.x,0,0,0);
     *b++ = allocated(1,0,0) ? val(v.x,1,0,0) : undefined;
   }
#line 358
 {
     *b++ = val(v.y,0,0,0);
     *b++ = allocated(0,1,0) ? val(v.y,0,1,0) : undefined;
   }
#line 358
 {
     *b++ = val(v.z,0,0,0);
     *b++ = allocated(0,0,1) ? val(v.z,0,0,1) : undefined;
   }}
      } } end_foreach_cache_level(); }





      MPI_Isend (rcv->buf, (b - (double *) rcv->buf),
   MPI_DOUBLE, rcv->pid, (l), MPI_COMM_WORLD,
   &rcv->r);
    }
  }

  prof_stop();
}

static void rcv_pid_sync (SndRcv * m, scalar * list, int l)
{
  scalar * listr = NULL;
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i63 = list; ((scalar *)&s)->i >= 0; s = *++_i63)
    if (!is_constant(s)) {
      if (_attribute[s.i].face)
 listf = vectors_add (listf, _attribute[s.i].v);
      else
 listr = list_add (listr, s);
    }
  rcv_pid_send (m->snd, listr, listf, l);
  rcv_pid_receive (m->rcv, listr, listf, l);
  rcv_pid_wait (m->snd);
  pfree (listr,__func__,__FILE__,__LINE__);
  pfree (listf,__func__,__FILE__,__LINE__);
}

static void snd_rcv_destroy (SndRcv * m)
{
  rcv_pid_destroy (m->rcv);
  rcv_pid_destroy (m->snd);
}

static void snd_rcv_init (SndRcv * m, const char * name)
{
  char s[strlen(name) + 5];
  strcpy (s, name);
  strcat (s, ".rcv");
  m->rcv = rcv_pid_new (s);
  strcpy (s, name);
  strcat (s, ".snd");
  m->snd = rcv_pid_new (s);
}

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  snd_rcv_destroy (&m->mpi_level);
  snd_rcv_destroy (&m->mpi_level_root);
  snd_rcv_destroy (&m->restriction);
  array_free (m->send);
  array_free (m->receive);
  pfree (m,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_level", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 425);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->mpi_level, list, l);
  rcv_pid_sync (&m->mpi_level_root, list, l);
 end_trace("mpi_boundary_level", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 429); }


static void mpi_boundary_restriction (const Boundary * b, scalar * list, int l)
{ trace ("mpi_boundary_restriction", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 433);
  MpiBoundary * m = (MpiBoundary *) b;
  rcv_pid_sync (&m->restriction, list, l);
 end_trace("mpi_boundary_restriction", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 436); }

void mpi_boundary_new()
{
  mpi_boundary = (Boundary *) ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  mpi_boundary->destroy = mpi_boundary_destroy;
  mpi_boundary->level = mpi_boundary_level;
  mpi_boundary->restriction = mpi_boundary_restriction;
  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;
  snd_rcv_init (&mpi->mpi_level, "mpi_level");
  snd_rcv_init (&mpi->mpi_level_root, "mpi_level_root");
  snd_rcv_init (&mpi->restriction, "restriction");
  mpi->send = array_new();
  mpi->receive = array_new();
  add_boundary (mpi_boundary);
}

static FILE * fopen_prefix (FILE * fp, const char * name, char * prefix)
{
  if (fp) {
    sprintf (prefix, "%s-%d ", name, pid());
    return fp;
  }
  else {
    strcpy (prefix, "");
    char fname[80];
    if (debug_iteration >= 0)
      sprintf (fname, "%s-%d-%d", name, debug_iteration, pid());
    else
      sprintf (fname, "%s-%d", name, pid());
    return fopen (fname, "w");
  }
}

void debug_mpi (FILE * fp1)
{
  void output_cells (FILE * fp);

  char prefix[80];
  FILE * fp;


  if (fp1 == NULL) {
    char name[80];
    sprintf (name, "halo-%d", pid()); remove (name);
    sprintf (name, "cells-%d", pid()); remove (name);
    sprintf (name, "faces-%d", pid()); remove (name);
    sprintf (name, "neighbors-%d", pid()); remove (name);
    sprintf (name, "mpi-level-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-level-root-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-rcv-%d", pid()); remove (name);
    sprintf (name, "mpi-restriction-snd-%d", pid()); remove (name);
    sprintf (name, "mpi-border-%d", pid()); remove (name);
    sprintf (name, "exterior-%d", pid()); remove (name);
    sprintf (name, "depth-%d", pid()); remove (name);
    sprintf (name, "refined-%d", pid()); remove (name);
  }


  fp = fopen_prefix (fp1, "halo", prefix);
  for (int l = 0; l < depth(); l++)
     { foreach_halo (prolongation, l){

#line 499 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

       { foreach_child()
      fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); end_foreach_child(); }; } end_foreach_halo(); }
  if (!fp1)
    fclose (fp);

  if (!fp1) {
    fp = fopen_prefix (fp1, "cells", prefix);
    output_cells (fp);
    fclose (fp);
  }

  fp = fopen_prefix (fp1, "faces", prefix);
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 512
{

#line 512 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 512
{

#line 512 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 512
{

#line 512 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); }  }}  end_foreach_face_generic()
#line 513
 end_foreach_face(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "neighbors", prefix);
   { foreach(){

#line 518 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    int n = 0;
     { foreach_neighbor(1)
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
 n++; end_foreach_neighbor(); }
    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, cell.neighbors);
    assert (cell.neighbors == n);
  } } end_foreach(); }
  if (!fp1)
    fclose (fp);

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;

  fp = fopen_prefix (fp1, "mpi-level-rcv", prefix);
  rcv_pid_print (mpi->mpi_level.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-rcv", prefix);
  rcv_pid_print (mpi->mpi_level_root.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-rcv", prefix);
  rcv_pid_print (mpi->restriction.rcv, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-snd", prefix);
  rcv_pid_print (mpi->mpi_level.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-level-root-snd", prefix);
  rcv_pid_print (mpi->mpi_level_root.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-restriction-snd", prefix);
  rcv_pid_print (mpi->restriction.snd, fp, prefix);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "mpi-border", prefix);
   { foreach_cell(){

#line 562 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (is_border(cell))
      fprintf (fp, "%s%g %g %g %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors, cell.pid);
    else
      continue;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "exterior", prefix);
   { foreach_cell(){

#line 575 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (!is_local(cell))
      fprintf (fp, "%s%g %g %g %d %d %d %d\n",
        prefix, x, y, z, level, cell.neighbors,
        cell.pid, cell.flags & leaf);






  } } end_foreach_cell(); }
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "depth", prefix);
  fprintf (fp, "depth: %d %d\n", pid(), depth());
  fprintf (fp, "======= mpi_level.snd ======\n");
  RcvPid * snd = mpi->mpi_level.snd;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  fprintf (fp, "======= mpi_level.rcv ======\n");
  snd = mpi->mpi_level.rcv;
  for (int i = 0; i < snd->npid; i++)
    fprintf (fp, "%d %d %d\n", pid(), snd->rcv[i].pid, snd->rcv[i].maxdepth);
  if (!fp1)
    fclose (fp);

  fp = fopen_prefix (fp1, "refined", prefix);
   { foreach_cache (((Tree *)grid)->refined){

#line 604 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    fprintf (fp, "%s%g %g %g %d\n", prefix, x, y, z, level); } end_foreach_cache(); }
  if (!fp1)
    fclose (fp);
}

static void snd_rcv_free (SndRcv * p)
{
  char name[strlen(p->rcv->name) + 1];
  strcpy (name, p->rcv->name);
  rcv_pid_destroy (p->rcv);
  p->rcv = rcv_pid_new (name);
  strcpy (name, p->snd->name);
  rcv_pid_destroy (p->snd);
  p->snd = rcv_pid_new (name);
}

static bool is_root (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 622 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

  if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
     { foreach_child()
      if (is_local(cell))
 return true; end_foreach_child(); }
  return false;
}


static bool is_local_prolongation (Point point, Point p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 632 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"




  struct { int x, y, z; } dp = {p.i - point.i, p.j - point.j, p.k - point.k};

  {
#line 638
 {
    if (dp.x == 0 && ((!is_leaf (neighbor(-1,0,0)) && neighbor(-1,0,0).neighbors && neighbor(-1,0,0).pid >= 0) || (!is_leaf (neighbor(1,0,0)) && neighbor(1,0,0).neighbors && neighbor(1,0,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(dp.x,0,0)) && neighbor(dp.x,0,0).neighbors && neighbor(dp.x,0,0).pid >= 0))
      return true;
  }
#line 638
 {
    if (dp.y == 0 && ((!is_leaf (neighbor(0,-1,0)) && neighbor(0,-1,0).neighbors && neighbor(0,-1,0).pid >= 0) || (!is_leaf (neighbor(0,1,0)) && neighbor(0,1,0).neighbors && neighbor(0,1,0).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,dp.y,0)) && neighbor(0,dp.y,0).neighbors && neighbor(0,dp.y,0).pid >= 0))
      return true;
  }
#line 638
 {
    if (dp.z == 0 && ((!is_leaf (neighbor(0,0,-1)) && neighbor(0,0,-1).neighbors && neighbor(0,0,-1).pid >= 0) || (!is_leaf (neighbor(0,0,1)) && neighbor(0,0,1).neighbors && neighbor(0,0,1).pid >= 0)))
      return true;
    if ((!is_leaf (neighbor(0,0,dp.z)) && neighbor(0,0,dp.z).neighbors && neighbor(0,0,dp.z).pid >= 0))
      return true;
  }}
  return false;
}



static void append_pid (Array * pids, int pid)
{
  for (int i = 0, * p = (int *) pids->p; i < pids->len/sizeof(int); i++, p++)
    if (*p == pid)
      return;
  array_append (pids, &pid, sizeof(int));
}

static int locals_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 658 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

  if (is_leaf(cell)) {
    if (is_local(cell)) {
      Point p = point;
       { foreach_neighbor(1) {
 if ((cell.pid >= 0 && cell.pid != pid()) &&
     ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p)))
   append_pid (pids, cell.pid);
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
    { foreach_child()
     if ((cell.pid >= 0 && cell.pid != pid()))
       append_pid (pids, cell.pid); end_foreach_child(); }
      } end_foreach_neighbor(); }
    }
  }
  else
     { foreach_neighbor(1) {
      if ((cell.pid >= 0 && cell.pid != pid()))
 append_pid (pids, cell.pid);
      if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0))
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
     append_pid (pids, cell.pid); end_foreach_child(); }
    } end_foreach_neighbor(); }
  return pids->len/sizeof(int);
}

static int root_pids (Point point, Array * pids)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 686 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if ((cell.pid >= 0 && cell.pid != pid()))
      append_pid (pids, cell.pid); end_foreach_child(); }
  return pids->len/sizeof(int);
}







static void rcv_pid_row (RcvPid * m, int l, int * row)
{
  for (int i = 0; i < npe(); i++)
    row[i] = 0;
  for (int i = 0; i < m->npid; i++) {
    Rcv * rcv = &m->rcv[i];
    if (l <= rcv->depth && rcv->halo[l].n > 0)
      row[rcv->pid] = rcv->halo[l].n;
  }
}

void check_snd_rcv_matrix (SndRcv * sndrcv, const char * name)
{
  int maxlevel = depth();
  mpi_all_reduce (maxlevel, MPI_INT, MPI_MAX);
  int * row = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
  for (int l = 0; l <= maxlevel; l++) {
    int status = 0;
    if (pid() == 0) {


      int ** send = matrix_new (npe(), npe(), sizeof(int));
      int ** receive = matrix_new (npe(), npe(), sizeof(int));
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, &send[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, &receive[0][0], npe(), MPI_INT, 0,
    MPI_COMM_WORLD);

      int * astatus = ((int *) pmalloc ((npe())*sizeof(int),__func__,__FILE__,__LINE__));
      for (int i = 0; i < npe(); i++)
 astatus[i] = 0;
      for (int i = 0; i < npe(); i++)
 for (int j = 0; j < npe(); j++)
   if (send[i][j] != receive[j][i]) {
     fprintf (qstderr(), "%s: %d sends    %d to   %d at level %d\n",
       name, i, send[i][j], j, l);
     fprintf (qstderr(), "%s: %d receives %d from %d at level %d\n",
       name, j, receive[j][i], i, l);
     fflush (qstderr());
     for (int k = i - 2; k <= i + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
     for (int k = j - 2; k <= j + 2; k++)
       if (k >= 0 && k < npe())
  astatus[k] = 1;
   }
      MPI_Scatter (astatus, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      pfree (astatus,__func__,__FILE__,__LINE__);

      matrix_free (send);
      matrix_free (receive);
    }
    else {
      rcv_pid_row (sndrcv->snd, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      rcv_pid_row (sndrcv->rcv, l, row);
      MPI_Gather (row, npe(), MPI_INT, NULL, npe(), MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Scatter (NULL, 1, MPI_INT, &status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    if (status) {
      fprintf (qstderr(),
        "check_snd_rcv_matrix \"%s\" failed\n"
        "Calling debug_mpi(NULL)...\n"
        "Aborting...\n",
        name);
      fflush (qstderr());
      debug_mpi (NULL);
      MPI_Abort (MPI_COMM_WORLD, -3);
    }
  }
  pfree (row,__func__,__FILE__,__LINE__);
}

static bool has_local_child (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 775 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

   { foreach_child()
    if (is_local(cell))
      return true; end_foreach_child(); }
  return false;
}


void mpi_boundary_update_buffers()
{ trace ("mpi_boundary_update_buffers", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 784);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_update_buffers", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 786);  return; }

  prof_start ("mpi_boundary_update_buffers");

  MpiBoundary * m = (MpiBoundary *) mpi_boundary;
  SndRcv * mpi_level = &m->mpi_level;
  SndRcv * mpi_level_root = &m->mpi_level_root;
  SndRcv * restriction = &m->restriction;

  snd_rcv_free (mpi_level);
  snd_rcv_free (mpi_level_root);
  snd_rcv_free (restriction);

  static const unsigned short used = 1 << user;
   { foreach_cell(){

#line 800 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell) && !is_border(cell))



      continue;

    if (cell.neighbors) {

      Array pids = {NULL, 0, 0};
      int n = locals_pids (point, &pids);
      if (n) {
  { foreach_child()
   if (is_local(cell))
     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (mpi_level->snd, *p, point); end_foreach_child(); }
 pfree (pids.p,__func__,__FILE__,__LINE__);
      }

      bool locals = false;
      if (is_leaf(cell)) {
 if ((cell.pid >= 0 && cell.pid != pid())) {
   Point p = point;
    { foreach_neighbor(1)
     if ((is_local(cell) &&
   ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) || is_local_prolongation (point, p))) ||
  is_root(point))
       locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
 }
      }
      else
  { foreach_neighbor(1)
   if (is_local(cell) || is_root(point))
     locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
  { foreach_child()
   if ((cell.pid >= 0 && cell.pid != pid()))
            rcv_pid_append (mpi_level->rcv, cell.pid, point),
       cell.flags |= used; end_foreach_child(); }


      if (!is_leaf(cell)) {

 if (is_local(cell)) {
   Array pids = {NULL, 0, 0};

   int n = root_pids (point, &pids);
   if (n) {
      { foreach_neighbor()
       for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
  if (cell.pid >= 0 && cell.pid != *p)
    rcv_pid_append (mpi_level_root->snd, *p, point); end_foreach_neighbor(); }

     for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
       rcv_pid_append (restriction->snd, *p, point);
     pfree (pids.p,__func__,__FILE__,__LINE__);
   }
 }

 else if ((cell.pid >= 0 && cell.pid != pid())) {
   bool root = false;
    { foreach_child()
     if (is_local(cell))
       root = true, foreach_child_break(); end_foreach_child(); }
   if (root) {
     int pid = cell.pid;
      { foreach_neighbor()
       if ((cell.pid >= 0 && cell.pid != pid()))
  rcv_pid_append (mpi_level_root->rcv, pid, point),
    cell.flags |= used; end_foreach_neighbor(); }

     rcv_pid_append (restriction->rcv, pid, point);
   }
 }
      }
    }


    if (level > 0) {
      if (is_local(cell)) {

 Array pids = {NULL, 0, 0};
 if ((aparent(0,0,0).pid >= 0 && aparent(0,0,0).pid != pid()))
   append_pid (&pids, aparent(0,0,0).pid);
 int n = root_pids (parent, &pids);
 if (n) {
   for (int i = 0, * p = (int *) pids.p; i < n; i++, p++)
     rcv_pid_append (restriction->snd, *p, point);
   pfree (pids.p,__func__,__FILE__,__LINE__);
 }
      }
      else if ((cell.pid >= 0 && cell.pid != pid())) {

 if (is_local(aparent(0,0,0)) || has_local_child (parent))
   rcv_pid_append (restriction->rcv, cell.pid, point);
      }
    }
  } } end_foreach_cell(); }





  static const unsigned short keep = 1 << (user + 1);
  for (int l = depth(); l >= 0; l--)
     { foreach_cell(){

#line 905 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

      if (level == l) {
 if (level > 0 && (cell.pid < 0 || is_local(cell) || (cell.flags & used)))
   aparent(0,0,0).flags |= keep;
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !(cell.flags & keep))
   coarsen_cell (point, NULL);
 cell.flags &= ~(used|keep);
 continue;
      } } end_foreach_cell(); }


  m->send->len = m->receive->len = 0;
  rcv_pid_append_pids (mpi_level->snd, m->send);
  rcv_pid_append_pids (mpi_level_root->snd, m->send);
  rcv_pid_append_pids (mpi_level->rcv, m->receive);
  rcv_pid_append_pids (mpi_level_root->rcv, m->receive);

  prof_stop();
#line 937 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 end_trace("mpi_boundary_update_buffers", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 937); }


void mpi_boundary_refine (scalar * list)
{ trace ("mpi_boundary_refine", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 941);
  prof_start ("mpi_boundary_refine");

  MpiBoundary * mpi = (MpiBoundary *) mpi_boundary;


  Array * snd = mpi->send;
  MPI_Request r[2*snd->len/sizeof(int)];
  int nr = 0;
  for (int i = 0, * dest = snd->p; i < snd->len/sizeof(int); i++,dest++) {
    int len = ((Tree *)grid)->refined.n;
    MPI_Isend (&((Tree *)grid)->refined.n, 1, MPI_INT, *dest,
        (128), MPI_COMM_WORLD, &r[nr++]);
    if (len > 0)
      MPI_Isend (((Tree *)grid)->refined.p, sizeof(Index)/sizeof(int)*len,
   MPI_INT, *dest, (128), MPI_COMM_WORLD, &r[nr++]);
  }



  Array * rcv = mpi->receive;
  Cache rerefined = {NULL, 0, 0};
  for (int i = 0, * source = rcv->p; i < rcv->len/sizeof(int); i++,source++) {
    int len;
    mpi_recv_check (&len, 1, MPI_INT, *source, (128),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE,
      "mpi_boundary_refine (len)");
    if (len > 0) {
      Index p[len];
      mpi_recv_check (p, sizeof(Index)/sizeof(int)*len,
        MPI_INT, *source, (128),
        MPI_COMM_WORLD, MPI_STATUS_IGNORE,
        "mpi_boundary_refine (p)");
      Cache refined = {p, len, len};
       { foreach_cache (refined){

#line 975 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

 if (level <= depth() && allocated(0,0,0)) {
   if (is_leaf(cell)) {
     bool neighbors = false;
      { foreach_neighbor()
       if (allocated(0,0,0) && (is_active(cell) || is_local(aparent(0,0,0))))
  neighbors = true, foreach_neighbor_break(); end_foreach_neighbor(); }

     if (neighbors)
       refine_cell (point, list, 0, &rerefined);
   }
 } } end_foreach_cache(); }
    }
  }


  if (nr)
    MPI_Waitall (nr, r, MPI_STATUSES_IGNORE);


  pfree (((Tree *)grid)->refined.p,__func__,__FILE__,__LINE__);
  ((Tree *)grid)->refined = rerefined;

  prof_stop();



  mpi_all_reduce (rerefined.n, MPI_INT, MPI_SUM);
  if (rerefined.n)
    mpi_boundary_refine (list);
 end_trace("mpi_boundary_refine", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1005); }

static void check_depth()
{
#line 1040 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
}

typedef struct {
  int refined, leaf;
} Remote;




void mpi_boundary_coarsen (int l, int too_fine)
{ trace ("mpi_boundary_coarsen", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1050);
  if (npe() == 1)
    { ; end_trace("mpi_boundary_coarsen", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1052);  return; }

  check_depth();

  assert (sizeof(Remote) == sizeof(double));

  scalar remote= new_scalar("remote");
   { foreach_cell(){

#line 1059 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (is_local(cell)) {
 ((Remote *)&val(remote,0,0,0))->refined = (!is_leaf (cell) && cell.neighbors && cell.pid >= 0);
 ((Remote *)&val(remote,0,0,0))->leaf = is_leaf(cell);
      }
      else {
 ((Remote *)&val(remote,0,0,0))->refined = true;
 ((Remote *)&val(remote,0,0,0))->leaf = false;
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);

   { foreach_cell(){

#line 1076 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (level == l) {
      if (!is_local(cell)) {
 if ((!is_leaf (cell) && cell.neighbors && cell.pid >= 0) && !((Remote *)&val(remote,0,0,0))->refined)
   coarsen_cell_recursive (point, NULL);
 else if (is_leaf(cell) && cell.neighbors && ((Remote *)&val(remote,0,0,0))->leaf) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
      }
      continue;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  check_depth();

  if (l > 0) {
     { foreach_cell(){

#line 1096 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 val(remote,0,0,0) = is_local(cell) ? cell.neighbors : 0;
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
    mpi_boundary_level (mpi_boundary, ((scalar []){remote,{-1}}), l);
     { foreach_cell(){

#line 1105 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l)
 if (!is_local(cell) && is_local(aparent(0,0,0)) && val(remote,0,0,0)) {
   aparent(0,0,0).flags &= ~too_fine;
   continue;
 }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
 delete (((scalar []){remote,{-1}}));  end_trace("mpi_boundary_coarsen", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1115); }

static void flag_border_cells()
{
   { foreach_cell(){

#line 1119 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    if (is_active(cell)) {
      short flags = cell.flags & ~border;
       { foreach_neighbor() {
 if (!is_local(cell) || (level > 0 && !is_local(aparent(0,0,0))))
   flags |= border, foreach_neighbor_break();

 if (is_refined_check())
    { foreach_child()
     if (!is_local(cell))
       flags |= border, foreach_child_break(); end_foreach_child(); }
 if (flags & border)
   foreach_neighbor_break();
      } end_foreach_neighbor(); }
      cell.flags = flags;
    }
    else {
      cell.flags &= ~border;

    }
    if (is_leaf(cell)) {
      if (cell.neighbors) {
  { foreach_child()
   cell.flags &= ~border; end_foreach_child(); }
 if (is_border(cell)) {
   bool remote = false;
    { foreach_neighbor (2/2)
     if (!is_local(cell))
       remote = true, foreach_neighbor_break(); end_foreach_neighbor(); }
   if (remote)
      { foreach_child()
       cell.flags |= border; end_foreach_child(); }
 }
      }
      continue;
    }
  } } end_foreach_cell(); }
}

static int balanced_pid (long index, long nt, int nproc)
{
  long ne = max(1, nt/nproc), nr = nt % nproc;
  int pid = index < nr*(ne + 1) ?
    index/(ne + 1) :
    nr + (index - nr*(ne + 1))/ne;
  return min(nproc - 1, pid);
}



void mpi_partitioning()
{ trace ("mpi_partitioning", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1170);
  prof_start ("mpi_partitioning");

  long nt = 0;
   { foreach(){

#line 1174 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    nt++; } end_foreach(); }


  long i = 0;
  ((Tree *)grid)->dirty = true;
   { foreach_cell_post (is_active (cell)){

#line 1180 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    if (is_active (cell)) {
      if (is_leaf (cell)) {
 cell.pid = balanced_pid (i++, nt, npe());
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else {
 cell.pid = child(0,0,0).pid;
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    } } end_foreach_cell_post(); }

  flag_border_cells();

  prof_stop();

  mpi_boundary_update_buffers();
 end_trace("mpi_partitioning", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1208); }

void restore_mpi (FILE * fp, scalar * list1)
{
  long index = 0, nt = 0, start = ftell (fp);
  scalar size= new_scalar("size"), * list = list_concat (((scalar []){size,{-1}}), list1);;
  long offset = sizeof(double)*list_len(list);


  static const unsigned short set = 1 << user;
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y,fm.z},{{-1},{-1},{-1}}});
   { foreach_cell(){

#line 1219 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    if (balanced_pid (index, nt, npe()) <= pid()) {
      unsigned flags;
      if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
 fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
 exit (1);
      }
      if (list) for (scalar s = *list, *_i64 = list; ((scalar *)&s)->i >= 0; s = *++_i64) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      if (level == 0)
 nt = val(size,0,0,0);
      cell.pid = balanced_pid (index, nt, npe());
      cell.flags |= set;
      if (!(flags & leaf) && is_leaf(cell)) {
 if (balanced_pid (index + val(size,0,0,0) - 1, nt, npe()) < pid()) {
   fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
   index += val(size,0,0,0);
   continue;
 }
 refine_cell (point, listm, 0, NULL);
      }
      index++;
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }


  fseek (fp, start, SEEK_SET);
  index = 0;
   { foreach_cell(){

#line 1255 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (qstderr(), "restore(): error: expecting 'flags'\n");
      exit (1);
    }
    if (cell.flags & set)
      fseek (fp, offset, SEEK_CUR);
    else {
      if (list) for (scalar s = *list, *_i65 = list; ((scalar *)&s)->i >= 0; s = *++_i65) {
 double val;
 if (fread (&val, sizeof(double), 1, fp) != 1) {
   fprintf (qstderr(), "restore(): error: expecting a scalar\n");
   exit (1);
 }
 if (s.i != INT_MAX)
   val(s,0,0,0) = val;
      }
      cell.pid = balanced_pid (index, nt, npe());
      if (is_leaf(cell) && cell.neighbors) {
 int pid = cell.pid;
  { foreach_child()
   cell.pid = pid; end_foreach_child(); }
      }
    }
    if (!(flags & leaf) && is_leaf(cell)) {
      bool locals = false;
       { foreach_neighbor(1)
 if ((cell.flags & set) && (is_local(cell) || is_root(point)))
   locals = true, foreach_neighbor_break(); end_foreach_neighbor(); }
      if (locals)
 refine_cell (point, listm, 0, NULL);
      else {
 fseek (fp, (sizeof(unsigned) + offset)*(val(size,0,0,0) - 1), SEEK_CUR);
 index += val(size,0,0,0);
 continue;
      }
    }
    index++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }


   { foreach_cell_post (is_active (cell)){

#line 1299 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
    cell.flags &= ~set;
    if (is_active (cell)) {
      if (is_leaf (cell)) {
 if (cell.neighbors > 0) {
   int pid = cell.pid;
    { foreach_child()
     cell.pid = pid; end_foreach_child(); }
 }
 if (!is_local(cell))
   cell.flags &= ~active;
      }
      else if (!is_local(cell)) {
 bool inactive = true;
  { foreach_child()
   if (is_active(cell))
     inactive = false, foreach_child_break(); end_foreach_child(); }
 if (inactive)
   cell.flags &= ~active;
      }
    }
  } } end_foreach_cell_post(); }

  flag_border_cells();

  mpi_boundary_update (list);
  pfree (list,__func__,__FILE__,__LINE__);
 delete (((scalar []){size,{-1}})); }
#line 1348 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1350);



  scalar size= new_scalar("size");
  subtree_size (size, leaves);






  double maxi = -1.;
  if (pid() == 0)
     { foreach_level(0){

#line 1364 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

      maxi = val(size,0,0,0) - 1.; } end_foreach_level(); }




   { foreach_level(0){

#line 1370 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"

    val(index,0,0,0) = 0; } end_foreach_level(); }
  for (int l = 0; l < depth(); l++) {
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), l); };
     { foreach_cell(){

#line 1374 "/home/vinlinux/basilisk/src/grid/tree-mpi.h"
 {
      if (level == l) {
 if (is_leaf(cell)) {
   if (is_local(cell) && cell.neighbors) {
     int i = val(index,0,0,0);
      { foreach_child()
       val(index,0,0,0) = i; end_foreach_child(); }
   }
 }
 else {
   bool loc = is_local(cell);
   if (!loc)
      { foreach_child()
       if (is_local(cell))
  loc = true, foreach_child_break(); end_foreach_child(); }
   if (loc) {
     int i = val(index,0,0,0) + !leaves;
      { foreach_child() {
       val(index,0,0,0) = i;
       i += val(size,0,0,0);
     } end_foreach_child(); }
   }
 }
 continue;
      }
      if (is_leaf(cell))
 continue;
    } } end_foreach_cell(); }
  }
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){index,{-1}}), depth()); };

  { double _ret =  maxi; delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1405);  return _ret; }
 delete (((scalar []){size,{-1}}));  end_trace("z_indexing", "/home/vinlinux/basilisk/src/grid/tree-mpi.h", 1406); }
#line 1783 "/home/vinlinux/basilisk/src/grid/tree.h"
#line 1 "grid/balance.h"
#line 1 "/home/vinlinux/basilisk/src/grid/balance.h"


typedef struct {
  short leaf, prolongation;
  int pid;
} NewPid;



#if TRASH
# define is_newpid() (!isnan(val(newpid,0,0,0)) && ((NewPid *)&val(newpid,0,0,0))->pid > 0)
#else
# define is_newpid() (((NewPid *)&val(newpid,0,0,0))->pid > 0)
#endif

Array * linear_tree (size_t size, scalar newpid)
{
  const unsigned short sent = 1 << user, next = 1 << (user + 1);
  Array * a = array_new();

   { foreach_cell_post_all (true){

#line 21 "/home/vinlinux/basilisk/src/grid/balance.h"

    if (level > 0 && (cell.flags & (sent|next)))
      aparent(0,0,0).flags |= next; } end_foreach_cell_post_all(); }

  bool empty = true;
   { foreach_cell_all(){

#line 26 "/home/vinlinux/basilisk/src/grid/balance.h"
 {
    if (cell.flags & sent) {
      array_append (a, &cell, size);
      cell.flags &= ~sent;
      empty = false;
    }
    else {
      if (cell.pid >= 0 && ((NewPid *)&val(newpid,0,0,0))->leaf)
 assert (is_leaf(cell));
      if (is_refined_check()) {


 bool prolo = false;
  { foreach_child()
   if (((NewPid *)&val(newpid,0,0,0))->prolongation)
     prolo = true; end_foreach_child(); }
 if (prolo) {

   cell.flags |= leaf;
   array_append (a, &cell, sizeof(Cell));
   cell.flags &= ~leaf;
 }
 else
   array_append (a, &cell, sizeof(Cell));
      }
      else
 array_append (a, &cell, sizeof(Cell));
    }
    if (cell.flags & next)
      cell.flags &= ~next;
    else
      continue;
  } } end_foreach_cell_all(); }

  if (empty)
    a->len = 0;
  return a;
}

#define foreach_tree(t, size, list)\
{\
  const unsigned short _sent = 1 << user, _next = 1 << (user + 1);\
  scalar * _list = list;\
  char * _i = (char *) (t)->p;\
  foreach_cell_all() {\
    Cell * c = (Cell *) _i;\
    if (c->flags & _sent) {\
      _i += size;\

#line 74


#define end_foreach_tree()\
    }\
    else\
      _i += sizeof(Cell);\
    if (c->flags & _next) {\
      assert (c->neighbors);\
      if (!(c->flags & leaf) && is_leaf(cell) &&\
   (!is_newpid() || !((NewPid *)&val(newpid,0,0,0))->leaf))\
\
 refine_cell (point, _list, 0, NULL);\
      else if (!cell.neighbors)\
\
 alloc_children (point);\
    }\
    else\
      continue;\
  } end_foreach_cell_all();\
}\

#line 94


Array * neighborhood (scalar newpid, int nextpid, FILE * fp)
{
  const unsigned short sent = 1 << user;
   { foreach_cell(){

#line 99 "/home/vinlinux/basilisk/src/grid/balance.h"
 {

    bool root = false;
    if ((!is_local(cell) || ((NewPid *)&val(newpid,0,0,0))->pid - 1 != nextpid) && (!is_leaf (cell) && cell.neighbors && cell.pid >= 0)) {
       { foreach_child()
 if (is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid)
   root = true, foreach_child_break(); end_foreach_child(); }
      if (root && cell.pid != nextpid) {
  { foreach_neighbor()
   if (cell.pid != nextpid && is_newpid()) {
     if (fp)
       fprintf (fp, "%g %g %g %d %d root\n",
         x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
     cell.flags |= sent;
   } end_foreach_neighbor(); }
      }
    }

    if ((is_local(cell) && ((NewPid *)&val(newpid,0,0,0))->pid - 1 == nextpid) || root) {
       { foreach_neighbor(1)
 if (cell.neighbors && cell.pid != nextpid)
    { foreach_child()
     if (cell.pid != nextpid && is_newpid()) {
       if (fp)
  fprintf (fp, "%g %g %g %d %d nextpid\n",
    x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
       cell.flags |= sent;
     } end_foreach_child(); } end_foreach_neighbor(); }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  return linear_tree (sizeof(Cell) + datasize, newpid);
}

static void send_tree (Array * a, int to, MPI_Request * r)
{
  MPI_Isend (&a->len, 1, MPI_LONG, to, (256), MPI_COMM_WORLD, &r[0]);
  if (a->len > 0) {
    MPI_Isend (a->p, a->len, MPI_BYTE, to, (256), MPI_COMM_WORLD, &r[1]);
    ((Tree *)grid)->dirty = true;
  }
}

static void receive_tree (int from, scalar newpid, FILE * fp)
{
  Array a;
  mpi_recv_check (&a.len, 1, MPI_LONG, from, (256),
    MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (len)");
  if (a.len > 0) {
    a.p = pmalloc (a.len,__func__,__FILE__,__LINE__);
    if (fp)
      fprintf (fp, "receiving %ld from %d\n", a.len, from);
    mpi_recv_check (a.p, a.len, MPI_BYTE, from, (256),
      MPI_COMM_WORLD, MPI_STATUS_IGNORE, "receive_tree (p)");

     { foreach_tree (&a, sizeof(Cell) + datasize, NULL){

#line 156 "/home/vinlinux/basilisk/src/grid/balance.h"
 {
      memcpy (((char *)&cell) + sizeof(Cell), ((char *)c) + sizeof(Cell),
       datasize);
      assert (((NewPid *)&val(newpid,0,0,0))->pid > 0);
      if (fp)
 fprintf (fp, "%g %g %g %d %d %d %d %d %d recv\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
   c->flags & leaf,
   cell.flags & leaf, from, ((NewPid *)&val(newpid,0,0,0))->leaf);
    } } end_foreach_tree(); }
    pfree (a.p,__func__,__FILE__,__LINE__);
    ((Tree *)grid)->dirty = true;
  }
}

static void wait_tree (Array * a, MPI_Request * r)
{
  MPI_Wait (&r[0], MPI_STATUS_IGNORE);
  if (a->len > 0)
    MPI_Wait (&r[1], MPI_STATUS_IGNORE);
}

static void check_flags()
{







}

struct {
  int min;
  bool leaves;

  int npe;
} mpi = {
  1,
  true
};


bool balance()
{ trace ("balance", "/home/vinlinux/basilisk/src/grid/balance.h", 201);
  if (npe() == 1)
    { bool _ret =  false; end_trace("balance", "/home/vinlinux/basilisk/src/grid/balance.h", 203);  return _ret; }

  assert (sizeof(NewPid) == sizeof(double));

  check_flags();

  long nl = 0, nt = 0;
   { foreach_cell(){

#line 210 "/home/vinlinux/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      nt++;
      if (is_leaf(cell))
 nl++;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  grid->n = grid->tn = nl;
  grid->maxdepth = depth();
  long nmin = nl, nmax = nl;

  mpi_all_reduce (nmax, MPI_LONG, MPI_MAX);
  mpi_all_reduce (nmin, MPI_LONG, MPI_MIN);
  mpi_all_reduce (grid->tn, MPI_LONG, MPI_SUM);
  mpi_all_reduce (grid->maxdepth, MPI_INT, MPI_MAX);
  if (mpi.leaves)
    nt = grid->tn;
  else
    mpi_all_reduce (nt, MPI_LONG, MPI_SUM);

  long ne = max(1, nt/npe());

  if (ne < mpi.min) {
    mpi.npe = max(1, nt/mpi.min);
    ne = max(1, nt/mpi.npe);
  }
  else
    mpi.npe = npe();

  if (nmax - nmin <= 1)
    { bool _ret =  false; end_trace("balance", "/home/vinlinux/basilisk/src/grid/balance.h", 243);  return _ret; }

  scalar newpid= new_scalar("newpid");
  double zn = z_indexing (newpid, mpi.leaves);
  if (pid() == 0)
    assert (zn + 1 == nt);

  FILE * fp = NULL;
#line 260 "/home/vinlinux/basilisk/src/grid/balance.h"
  bool next = false, prev = false;
   { foreach_cell_all(){

#line 261 "/home/vinlinux/basilisk/src/grid/balance.h"
 {
    if (is_local(cell)) {
      int pid = balanced_pid (val(newpid,0,0,0), nt, mpi.npe);
      pid = clamp (pid, cell.pid - 1, cell.pid + 1);
      if (pid == pid() + 1)
 next = true;
      else if (pid == pid() - 1)
 prev = true;
      ((NewPid *)&val(newpid,0,0,0))->pid = pid + 1;
      ((NewPid *)&val(newpid,0,0,0))->leaf = is_leaf(cell);
      ((NewPid *)&val(newpid,0,0,0))->prolongation = (!is_leaf(cell) && !cell.neighbors && cell.pid >= 0);
      if (fp)
 fprintf (fp, "%g %g %d %d newpid\n", x, y, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
    }
    else
      val(newpid,0,0,0) = 0;
  } } end_foreach_cell_all(); }
  for (int l = 0; l <= depth(); l++)
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, ((scalar []){newpid,{-1}}), l); };
#line 304 "/home/vinlinux/basilisk/src/grid/balance.h"
  Array * anext = next ? neighborhood (newpid, pid() + 1, fp) : array_new();
  Array * aprev = prev ? neighborhood (newpid, pid() - 1, fp) : array_new();

  if (fp)
    fflush (fp);

  check_flags();


  MPI_Request rprev[2], rnext[2];
  if (pid() > 0)
    send_tree (aprev, pid() - 1, rprev);
  if (pid() < npe() - 1)
    send_tree (anext, pid() + 1, rnext);


  if (pid() < npe() - 1)
    receive_tree (pid() + 1, newpid, fp);
  if (pid() > 0)
    receive_tree (pid() - 1, newpid, fp);


  if (pid() > 0)
    wait_tree (aprev, rprev);
  array_free (aprev);
  if (pid() < npe() - 1)
    wait_tree (anext, rnext);
  array_free (anext);

  if (fp)
    fflush (fp);


  int pid_changed = false;
   { foreach_cell_all(){

#line 338 "/home/vinlinux/basilisk/src/grid/balance.h"
 {
    if (cell.pid >= 0) {
      if (is_newpid()) {
 if (fp)
   fprintf (fp, "%g %g %g %d %d %d %d %d new\n",
     x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid,
     is_leaf(cell), cell.neighbors, ((NewPid *)&val(newpid,0,0,0))->leaf);
 if (cell.pid != ((NewPid *)&val(newpid,0,0,0))->pid - 1) {
   cell.pid = ((NewPid *)&val(newpid,0,0,0))->pid - 1;
   cell.flags &= ~(active|border);
   if (is_local(cell))
     cell.flags |= active;
   pid_changed = true;
 }
 if (((NewPid *)&val(newpid,0,0,0))->leaf && !is_leaf(cell) && cell.neighbors)
   coarsen_cell_recursive (point, NULL);
      }
      else if (level > 0 && ((NewPid *)&coarse(newpid,0,0,0))->leaf)
 cell.pid = aparent(0,0,0).pid;
    }

    if (!cell.neighbors && allocated_child(0,0,0)) {
      if (fp)
 fprintf (fp, "%g %g %g %d %d freechildren\n",
   x, y, z, ((NewPid *)&val(newpid,0,0,0))->pid - 1, cell.pid);
      free_children (point);
    }
  } } end_foreach_cell_all(); }

  if (((Tree *)grid)->dirty || pid_changed) {


     { foreach_cell_post (!is_leaf (cell)){

#line 370 "/home/vinlinux/basilisk/src/grid/balance.h"

      if (!is_leaf(cell) && !is_local(cell)) {
 unsigned short flags = cell.flags & ~active;
  { foreach_child()
   if (is_active(cell))
     flags |= active, foreach_child_break(); end_foreach_child(); }
 cell.flags = flags;
      } } end_foreach_cell_post(); }

    flag_border_cells();
    pid_changed = true;
  }

  if (fp)
    fclose (fp);

  mpi_all_reduce (pid_changed, MPI_INT, MPI_MAX);
  if (pid_changed)
    mpi_boundary_update_buffers();

  { bool _ret =  pid_changed; delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/vinlinux/basilisk/src/grid/balance.h", 390);  return _ret; }
 delete (((scalar []){newpid,{-1}}));  end_trace("balance", "/home/vinlinux/basilisk/src/grid/balance.h", 391); }

void mpi_boundary_update (scalar * list)
{
  mpi_boundary_update_buffers();
  grid->tn = 0;
  boundary (list);
  while (balance());
}
#line 1784 "/home/vinlinux/basilisk/src/grid/tree.h"
#else
void mpi_boundary_refine (scalar * list){}
void mpi_boundary_coarsen (int a, int b){}
void mpi_boundary_update (scalar * list) {
  boundary (list);
}
#endif
#line 4 "/home/vinlinux/basilisk/src/grid/octree.h"

void octree_methods() {
  tree_methods();
}
#line 8 "main.c"
#line 1 "view.h"
#line 1 "/home/vinlinux/basilisk/src/view.h"
#line 64 "/home/vinlinux/basilisk/src/view.h"
# include <GL/gl.h>
# include <GL/glu.h>


#include <gl/framebuffer.h>
#include <gl/gl2ps/gl2ps.h>
#include <gl/trackball.h>
#include <gl/utils.h>


#line 1 "utils.h"
#line 1 "/home/vinlinux/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if _MPI
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _n = n; 
#line 69
foreach(){

#line 69 "/home/vinlinux/basilisk/src/utils.h"
 _n++; } end_foreach();OMP(omp critical) n += _n;
mpi_all_reduce_double (n, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if 1
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if _MPI
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Octree"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if _MPI
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; double _avg = avg; double _rms = rms; double _volume = volume; 
#line 135

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/vinlinux/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/vinlinux/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) avg += _avg;
mpi_all_reduce_double (avg, MPI_SUM);
OMP(omp critical) rms += _rms;
mpi_all_reduce_double (rms, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = avg/volume;
  n.rms = sqrt(rms/volume);
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; double _sum2 = sum2; double _volume = volume; double _max = max; double _min = min; 
#line 163

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/vinlinux/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/vinlinux/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);
OMP(omp critical) sum2 += _sum2;
mpi_all_reduce_double (sum2, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);
OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) if (_min < min) min = _min;
mpi_all_reduce_double (min, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 186 "/home/vinlinux/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 212 "/home/vinlinux/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 236 "/home/vinlinux/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
   { foreach(){

#line 239 "/home/vinlinux/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i0 = f; vector * _i1 = g; if (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i0, v = *++_i1) {
      if (_attribute[s.i].gradient)
 {
#line 243

   val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
#line 243

   val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
#line 243

   val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;}
      else
 {
#line 246

   val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 246

   val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 246

   val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    }
  } } end_foreach(); }
  boundary ((scalar *) g);
}
#line 262 "/home/vinlinux/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{

     { foreach(){

#line 265 "/home/vinlinux/basilisk/src/utils.h"

      val(omega,0,0,0) = (val(u.y,1,0,0) - val(u.y,-1,0,0) + val(u.x,0,-1,0) - val(u.x,0,1,0))/(2.*Delta); } end_foreach(); }
    boundary (((scalar []){omega,{-1}}));

}





double change (scalar v, scalar vn)
{
  double max = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; 
#line 278
foreach(){

#line 278 "/home/vinlinux/basilisk/src/utils.h"
 {
    double dv = fabs (val(v,0,0,0) - val(vn,0,0,0));
    if (dv > _max)
      _max = dv;
    val(vn,0,0,0) = val(v,0,0,0);
  } } end_foreach();OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 283
 }
  return max;
}

#line 1 "output.h"
#line 1 "/home/vinlinux/basilisk/src/output.h"
#line 33 "/home/vinlinux/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/vinlinux/basilisk/src/output.h", 42);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
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
 if (p.list) for (scalar s = *p.list, *_i66 = p.list; ((scalar *)&s)->i >= 0; s = *++_i66)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});
      }
      else {
 Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 62 "/home/vinlinux/basilisk/src/output.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i67 = p.list; ((scalar *)&s)->i >= 0; s = *++_i67)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    if (p.list) for (scalar s = *p.list, *_i68 = p.list; ((scalar *)&s)->i >= 0; s = *++_i68)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + X0;
      for (int j = 0; j < p.n; j++) {
 double y = Delta*j + Y0;

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i69 = p.list; ((scalar *)&s)->i >= 0; s = *++_i69)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/vinlinux/basilisk/src/output.h", 102); }
#line 130 "/home/vinlinux/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/vinlinux/basilisk/src/output.h", 139);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 157 "/home/vinlinux/basilisk/src/output.h"

 assert (point.level >= 0);
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/vinlinux/basilisk/src/output.h", 165); }
#line 174 "/home/vinlinux/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  assert (i < 127 - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 311 "/home/vinlinux/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);



    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  assert (pid() == 0);
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  assert (pid() == 0);
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 541 "/home/vinlinux/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/vinlinux/basilisk/src/output.h", 556);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
    p.min = avg - spread; p.max = avg + spread;
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 594 "/home/vinlinux/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 604 "/home/vinlinux/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = qstdout();
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/vinlinux/basilisk/src/output.h", 636); }
#line 668 "/home/vinlinux/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/vinlinux/basilisk/src/output.h", 679);

  if (!p.fp) p.fp = qstdout();
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 715 "/home/vinlinux/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 725 "/home/vinlinux/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/vinlinux/basilisk/src/output.h", 737); }
#line 764 "/home/vinlinux/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/vinlinux/basilisk/src/output.h", 794);
  char * fname = p.file;

#if _MPI



  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = qstdout();
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);


  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if _MPI
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  if (list) for (scalar s = *list, *_i70 = list; ((scalar *)&s)->i >= 0; s = *++_i70)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 872 "/home/vinlinux/basilisk/src/output.h"
 {
#if _MPI
    if (is_local(cell))
#endif
    {
#if _MPI
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :
#line 893 "/home/vinlinux/basilisk/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      if (list) for (scalar s = *list, *_i71 = list; ((scalar *)&s)->i >= 0; s = *++_i71)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != nodata ? - val(s,0,0,0) : (double) DBL_MAX;
     }


     else
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;

   }
   else
     a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if _MPI
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if _MPI
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = qstdout();
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/vinlinux/basilisk/src/output.h", 970); }
#line 990 "/home/vinlinux/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  if (lista) for (scalar s = *lista, *_i72 = lista; ((scalar *)&s)->i >= 0; s = *++_i72)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  if (list) for (scalar s = *list, *_i73 = list; ((scalar *)&s)->i >= 0; s = *++_i73) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !_MPI

void dump (struct Dump p)
{ trace ("dump", "/home/vinlinux/basilisk/src/output.h", 1043);
  FILE * fp = p.fp;
  char * file = p.file;

  if (file && (fp = fopen (file, "w")) == NULL) {
    perror (file);
    exit (1);
  }
  assert (fp);

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 1062 "/home/vinlinux/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    if (list) for (scalar s = *list, *_i74 = list; ((scalar *)&s)->i >= 0; s = *++_i74)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/vinlinux/basilisk/src/output.h", 1080); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/vinlinux/basilisk/src/output.h", 1084);
  FILE * fp = p.fp;
  char * file = p.file;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  FILE * fh = fopen (file, "w");

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };







  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  if (list) for (scalar s = *list, *_i75 = list; ((scalar *)&s)->i >= 0; s = *++_i75)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);

  subtree_size (size, false);

   { foreach_cell(){

#line 1121 "/home/vinlinux/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      fseek (fh, offset, SEEK_SET);
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      if (list) for (scalar s = *list, *_i76 = list; ((scalar *)&s)->i >= 0; s = *++_i76)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/vinlinux/basilisk/src/output.h", 1139); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/vinlinux/basilisk/src/output.h", 1144);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/vinlinux/basilisk/src/output.h", 1148);  return _ret; }
  assert (fp);

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }


  init_grid (1);
   { foreach_cell(){

#line 1159 "/home/vinlinux/basilisk/src/output.h"
 {
    cell.pid = pid();
    cell.flags |= active;
  } } end_foreach_cell(); }
  ((Tree *)grid)->dirty = true;
#line 1184 "/home/vinlinux/basilisk/src/output.h"
  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 if (list) for (scalar s = *list, *_i77 = list; ((scalar *)&s)->i >= 0; s = *++_i77)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }
#line 1259 "/home/vinlinux/basilisk/src/output.h"
  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y,fm.z},{{-1},{-1},{-1}}});



   { foreach_cell(){

#line 1263 "/home/vinlinux/basilisk/src/output.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    if (list) for (scalar s = *list, *_i78 = list; ((scalar *)&s)->i >= 0; s = *++_i78) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (list);

  boundary (listm);

  scalar * other = NULL;
  if (all) for (scalar s = *all, *_i79 = all; ((scalar *)&s)->i >= 0; s = *++_i79)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/vinlinux/basilisk/src/output.h", 1309);  return _ret; }
 end_trace("restore", "/home/vinlinux/basilisk/src/output.h", 1310); }
#line 288 "/home/vinlinux/basilisk/src/utils.h"
#line 75 "/home/vinlinux/basilisk/src/view.h"
#line 1 "input.h"
#line 1 "/home/vinlinux/basilisk/src/input.h"
#line 16 "/home/vinlinux/basilisk/src/input.h"
struct InputPGM {

  scalar s;
  FILE * fp;

  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (qstderr(), "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (qstderr(), "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (qstderr(), "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (qstderr(), "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((width*height)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (qstderr(), "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 73 "/home/vinlinux/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((width*height)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (qstderr(), "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 96 "/home/vinlinux/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (qstderr(), "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 166 "/home/vinlinux/basilisk/src/input.h"

void input_gfs (struct OutputGfs p)
{ trace ("input_gfs", "/home/vinlinux/basilisk/src/input.h", 168);
  not_mpi_compatible();

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  bool input_all = (p.list == all);
  if (p.list == NULL) p.list = all;


  init_grid (1);


  next_char (p.fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (qstderr(), "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (qstderr(), "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (qstderr(), "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    if (p.list) for (scalar s = *p.list, *_i80 = p.list; ((scalar *)&s)->i >= 0; s = *++_i80)
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (qstderr(), "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (qstderr(), "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = ((scalar []){cm,fm.x,fm.y,fm.z,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

   { foreach_cell(){

#line 273 "/home/vinlinux/basilisk/src/input.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (qstderr(), "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (qstderr(), "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    if (input) for (scalar s = *input, *_i81 = input; ((scalar *)&s)->i >= 0; s = *++_i81) {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
 fprintf (qstderr(), "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }


   else
     val(s,0,0,0) = a;

 }
 else
   val(s,0,0,0) = a;
      }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (listm);
  boundary (input);

  pfree (input,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
 end_trace("input_gfs", "/home/vinlinux/basilisk/src/input.h", 328); }
#line 357 "/home/vinlinux/basilisk/src/input.h"
struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  fscanf (p.fp, "%s %d", waste, &nx);
  fscanf (p.fp, "%s %d", waste, &ny);
  fscanf (p.fp, "%s %lf", waste, &XG0);
  fscanf (p.fp, "%s %lf", waste, &YG0);
  fscanf (p.fp, "%s %lf", waste, &DeltaGRD);
  fscanf (p.fp, "%s %lf", waste, &ndv);


  if (!p.nodatavalue)
    p.nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++)
      fscanf (p.fp, "%lf ", &value[j + i*nx]);

  double LGx0 = nx*DeltaGRD;
  double LGy0 = ny*DeltaGRD;
  bool warning = false;

   { foreach(){

#line 409 "/home/vinlinux/basilisk/src/input.h"
 {

    bool incl = (x >= XG0 - DeltaGRD/2. &&
   x <= XG0 + LGx0 + DeltaGRD/2. &&
   y >= YG0 - DeltaGRD/2. &&
   y <= YG0 + LGy0 + DeltaGRD/2.);
    if (incl) {
      double val;

      bool ring = (x >= XG0 &&
     x <= XG0 + LGx0 &&
     y >= YG0 &&
     y <= YG0 + LGy0);
      if (p.linear && ring ) {
 int j = (x - XG0)/DeltaGRD; int i = (y - YG0)/DeltaGRD;
 double dx = x -(j*DeltaGRD + XG0); double dy = y - (i*DeltaGRD + YG0);
 val = value[j + i*nx]
   + dx*(value[j + 1 + i*nx] - value[j + i*nx])/DeltaGRD
   + dy*(value[j + (i + 1)*nx] - value[j + i*nx])/DeltaGRD
   + dx*dy*(value[j + i*nx] + value[j +1 + (i+1)*nx]
     - value[j + (i + 1)*nx]-value[j + 1 + i*nx])/sq(DeltaGRD);
      }
      else {
 int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
 int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
 val = value[j + i*nx];
      }
      val(input,0,0,0) = val;
      if (val == ndv)
 val(input,0,0,0) = p.nodatavalue;
    }
    else {
      val(input,0,0,0) = p.nodatavalue;
      warning = true;
    }
  } } end_foreach(); }
  pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (qstderr(),
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (opened)
    fclose (p.fp);
}
#line 76 "/home/vinlinux/basilisk/src/view.h"






struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;

  bool gfsview;
  bool vector;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  void (* map) (coord *);

  int ni;

  bool active;
};

typedef struct _bview bview;




bview * bview_new() {
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);




  p->bg[0] = 0.3; p->bg[1] = 0.4; p->bg[2] = 0.6;

  p->res = 1.;
  p->lc = 0.001;

  p->vector = false;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view()
{
  assert (_view);
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  double max = 2.;





  gluPerspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);
  glMatrixMode (GL_MODELVIEW);

  glLoadIdentity ();
  glTranslatef (view->tx, view->ty, - (1. + max));

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  return view;
}







typedef void * pointer;



static pointer compose_image (bview * view) { trace ("compose_image", "/home/vinlinux/basilisk/src/view.h", 239);
  { pointer _ret =  framebuffer_image((view)->fb); end_trace("compose_image", "/home/vinlinux/basilisk/src/view.h", 240);  return _ret; }
 end_trace("compose_image", "/home/vinlinux/basilisk/src/view.h", 241); }
#line 339 "/home/vinlinux/basilisk/src/view.h"
#line 1 "draw.h"
#line 1 "/home/vinlinux/basilisk/src/draw.h"




#line 1 "fractions.h"
#line 1 "/home/vinlinux/basilisk/src/fractions.h"
#line 11 "/home/vinlinux/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/vinlinux/basilisk/src/geometry.h"
#line 28 "/home/vinlinux/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 133 "/home/vinlinux/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 237 "/home/vinlinux/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
#line 240
 {
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }}
  return plane_volume (n1, alpha);
}
#line 277 "/home/vinlinux/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
    {
#line 312
 {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    }
#line 312
 {
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    }
#line 312
 {
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }}
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
 {
#line 323
 {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 }
#line 323
 {
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 }
#line 323
 {
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }}
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  if (n.x < 0.) {
    alpha -= n.x;
    n.x = - n.x;
  }
  if (n.y < 0.) {
    alpha -= n.y;
    n.y = - n.y;
  }

  p->x = p->y = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 371

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 371

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 397
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 397
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  if (fabs (m.x) < 1e-4) {
    coord n, q;
    n.x = m.y;
    n.y = m.z;
    double length = line_length_center (n, alpha, &q);
    p->x = 0.;
    p->y = q.x;
    p->z = q.y;
    return sq(length);
  }
  if (fabs (m.y) < 1e-4) {
    coord n, q;
    n.x = m.z;
    n.y = m.x;
    double length = line_length_center (n, alpha, &q);
    p->x = q.y;
    p->y = 0.;
    p->z = q.x;
    return sq(length);
  }
  if (fabs (m.z) < 1e-4) {
    double length = line_length_center (m, alpha, p);
    p->z = 0.;
    return sq(length);
  }

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  {
#line 441

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 441

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
#line 441

    if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }}

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

  {
#line 456
 {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  }
#line 456
 {
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  }
#line 456
 {
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }}

  amax = alpha - amax;
  {
#line 467
 {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  }
#line 467
 {
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  }
#line 467
 {
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }}

  area *= 3.;
  {
#line 478
 {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 478
 {
    if (area) {
      p->y /= area*n.y;
      p->y = clamp (p->y, 0., 1.);
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  }
#line 478
 {
    if (area) {
      p->z /= area*n.z;
      p->z = clamp (p->z, 0., 1.);
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }}

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}
#line 12 "/home/vinlinux/basilisk/src/fractions.h"
#line 20 "/home/vinlinux/basilisk/src/fractions.h"
#line 1 "myc.h"
#line 1 "/home/vinlinux/basilisk/src/myc.h"
#line 16 "/home/vinlinux/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 17 "/home/vinlinux/basilisk/src/myc.h"

  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 21 "/home/vinlinux/basilisk/src/fractions.h"
#line 32 "/home/vinlinux/basilisk/src/fractions.h"
void fraction_refine (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 33 "/home/vinlinux/basilisk/src/fractions.h"






  double cc = val(c,0,0,0);
  if (cc <= 0. || cc >= 1.)
     { foreach_child()
      val(c,0,0,0) = cc; end_foreach_child(); }
  else {




    coord n = mycs (point, c);
    double alpha = plane_alpha (cc, n);






    coord a, b;
    {
#line 57
 {
      a.x = 0.; b.x = 0.5;
    }
#line 57
 {
      a.y = 0.; b.y = 0.5;
    }
#line 57
 {
      a.z = 0.; b.z = 0.5;
    }}

     { foreach_child() {
      coord nc;
      {
#line 63

 nc.x = child.x*n.x;
#line 63

 nc.y = child.y*n.y;
#line 63

 nc.z = child.z*n.z;}
      val(c,0,0,0) = rectangle_fraction (nc, alpha, a, b);
    } end_foreach_child(); }
  }
}











static void alpha_refine (Point point, scalar alpha)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 81 "/home/vinlinux/basilisk/src/fractions.h"

  vector n = _attribute[alpha.i].n;
  double alphac = 2.*val(alpha,0,0,0);
  coord m;
  {
#line 85

    m.x = val(n.x,0,0,0);
#line 85

    m.y = val(n.y,0,0,0);
#line 85

    m.z = val(n.z,0,0,0);}
   { foreach_child() {
    val(alpha,0,0,0) = alphac;
    {
#line 89

      val(alpha,0,0,0) -= child.x*m.x/2.;
#line 89

      val(alpha,0,0,0) -= child.y*m.y/2.;
#line 89

      val(alpha,0,0,0) -= child.z*m.z/2.;}
  } end_foreach_child(); }
}
#line 116 "/home/vinlinux/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
};


void fractions (struct Fractions a)
{ trace ("fractions", "/home/vinlinux/basilisk/src/fractions.h", 124);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
#line 136 "/home/vinlinux/basilisk/src/fractions.h"
  vector p= new_vector("p");
#line 148 "/home/vinlinux/basilisk/src/fractions.h"
   { foreach_vertex(){

#line 148 "/home/vinlinux/basilisk/src/fractions.h"
 {
#line 148
 if (neighbor(1,0,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,1,0,0) < 0.) {






      val(p.x,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 173 "/home/vinlinux/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,1,0,0) > 0.);
  }
#line 148
 if (neighbor(0,1,0).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,1,0) < 0.) {






      val(p.y,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < 0.)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 173 "/home/vinlinux/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,1,0) > 0.);
  }
#line 148
 if (neighbor(0,0,1).flags & vertex) {





    if (val(Phi,0,0,0)*val(Phi,0,0,1) < 0.) {






      val(p.z,0,0,0) = val(Phi,0,0,0)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < 0.)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 173 "/home/vinlinux/basilisk/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > 0. || val(Phi,0,0,1) > 0.);
  }} } end_foreach_vertex(); }
#line 187 "/home/vinlinux/basilisk/src/fractions.h"
  scalar s_x = s.x, s_y = s.y, s_z = s.z;
   { foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 188
{

#line 188 "/home/vinlinux/basilisk/src/fractions.h"






  {
#line 226 "/home/vinlinux/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 228
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 245

 n.x /= nn;
#line 245

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 269 "/home/vinlinux/basilisk/src/fractions.h"
      val(s_z,0,0,0) = ni ? line_area (n.x, n.y, alpha/ni) : max (val(p.x,0,0,0), val(p.y,0,0,0));
    }
  } }  }}  { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 188
{

#line 188 "/home/vinlinux/basilisk/src/fractions.h"






  {
#line 226 "/home/vinlinux/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    }
#line 228
 {
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }}





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      {
#line 245

 n.y /= nn;
#line 245

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0))*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}
#line 269 "/home/vinlinux/basilisk/src/fractions.h"
      val(s_x,0,0,0) = ni ? line_area (n.y, n.z, alpha/ni) : max (val(p.y,0,0,0), val(p.z,0,0,0));
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 188
{

#line 188 "/home/vinlinux/basilisk/src/fractions.h"






  {
#line 226 "/home/vinlinux/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 228
 {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    }
#line 228
 {
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }}





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      {
#line 245

 n.z /= nn;
#line 245

 n.x /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 255

   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0))*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
#line 255

   if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i))*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}
#line 269 "/home/vinlinux/basilisk/src/fractions.h"
      val(s_y,0,0,0) = ni ? line_area (n.z, n.x, alpha/ni) : max (val(p.z,0,0,0), val(p.x,0,0,0));
    }
  } }  }}  end_foreach_face_generic()
#line 271
 end_foreach_face(); }







  boundary_flux (((vector []){{s.x,s.y,s.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 280 "/home/vinlinux/basilisk/src/fractions.h"
 {




    coord n;
    double nn = 0.;
    {
#line 287
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 287
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 287
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    if (nn == 0.)
      val(c,0,0,0) = val(s.x,0,0,0);
    else {
      {
#line 294

 n.x /= nn;
#line 294

 n.y /= nn;
#line 294

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
#line 305

     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j))*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
#line 305

     if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i))*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
#line 305

     if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0))*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      val(c,0,0,0) = ni ? plane_volume (n, alpha/ni) : val(s.x,0,0,0);
    }
  } } end_foreach(); }





  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){p.x,p.y,p.z,{-1}}));  { if (!(a.s).x.i) delete (((scalar []){s.x,s.y,s.z,{-1}})); }  end_trace("fractions", "/home/vinlinux/basilisk/src/fractions.h", 324); }
#line 349 "/home/vinlinux/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 350 "/home/vinlinux/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  assert (3 == 2);
  {
#line 354
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 354
 {
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  }
#line 354
 {
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }}

  if (nn > 0.)
    {
#line 361

      n.x /= nn;
#line 361

      n.y /= nn;
#line 361

      n.z /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 373 "/home/vinlinux/basilisk/src/fractions.h"

  coord n;
  if (s.x.i < 0)
    n = mycs (point, c);
  else {
    double nn = 0.;
    {
#line 379
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 379
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 379
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    assert (nn > 0.);
    {
#line 384

      n.x /= nn;
#line 384

      n.y /= nn;
#line 384

      n.z /= nn;}
  }
  return n;
}
#line 397 "/home/vinlinux/basilisk/src/fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "/home/vinlinux/basilisk/src/fractions.h", 399);
   { foreach(){

#line 400 "/home/vinlinux/basilisk/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 408

 val(n.x,0,0,0) = 0.;
#line 408

 val(n.y,0,0,0) = 0.;
#line 408

 val(n.z,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);

      {
#line 420

 val(n.x,0,0,0) = m.x;
#line 420

 val(n.y,0,0,0) = m.y;
#line 420

 val(n.z,0,0,0) = m.z;}
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 433 "/home/vinlinux/basilisk/src/fractions.h"
  {
#line 433

    _attribute[n.x.i].refine = _attribute[n.x.i].prolongation = refine_injection;
#line 433

    _attribute[n.y.i].refine = _attribute[n.y.i].prolongation = refine_injection;
#line 433

    _attribute[n.z.i].refine = _attribute[n.z.i].prolongation = refine_injection;}




  _attribute[alpha.i].n = n;
  _attribute[alpha.i].refine = _attribute[alpha.i].prolongation = alpha_refine;







  boundary (((scalar []){n.x,n.y,n.z,alpha,{-1}}));
 end_trace("reconstruction", "/home/vinlinux/basilisk/src/fractions.h", 449); }
#line 469 "/home/vinlinux/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "/home/vinlinux/basilisk/src/fractions.h", 477);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = qstdout();
  if (!s.x.i) s.x.i = -1;

   { foreach(){

#line 483 "/home/vinlinux/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (val(c,0,0,0), n);







      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (p.fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', p.fp);

    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "/home/vinlinux/basilisk/src/fractions.h", 505); }








double interface_area (scalar c)
{ trace ("interface_area", "/home/vinlinux/basilisk/src/fractions.h", 515);
  double area = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _area = area; 
#line 517
foreach (){

#line 517 "/home/vinlinux/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      _area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    } } end_foreach();OMP(omp critical) area += _area;
mpi_all_reduce_double (area, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 522
 }
  { double _ret =  area; end_trace("interface_area", "/home/vinlinux/basilisk/src/fractions.h", 523);  return _ret; }
 end_trace("interface_area", "/home/vinlinux/basilisk/src/fractions.h", 524); }
#line 6 "/home/vinlinux/basilisk/src/draw.h"
#line 1 "gl/font.h"
#line 1 "/home/vinlinux/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character );
void gl_StrokeString( const char *string );
float gl_StrokeWidth( int character );
float gl_StrokeLength( const char *string );
GLfloat gl_StrokeHeight( );
#line 7 "/home/vinlinux/basilisk/src/draw.h"




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 46 "/home/vinlinux/basilisk/src/draw.h"
struct _view_set {
  float tx, ty;
  float fov;
  float quat[4];
  float sx, sy, sz;
  unsigned width, height, samples;
  float bg[3];
  float theta, phi, psi;
  bool relative;
  float res;
  char * camera;
  void (* map) (coord *);
  float p1x, p1y, p2x, p2y;
  bview * view;
};

void view (struct _view_set p)
{
  bview * v = p.view ? p.view : get_view();
  if (p.fov) {
    if (p.relative)
      v->fov += (0.1 + 3.*v->fov)*p.fov;
    else
      v->fov = p.fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (p.quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = p.quat[j];
      break;
    }
  if (p.tx) v->tx = p.relative ? v->tx + p.tx*0.02*(0.01 + 3.*v->fov) : p.tx;
  if (p.ty) v->ty = p.relative ? v->ty + p.ty*0.02*(0.01 + 3.*v->fov) : p.ty;
  if (p.sx) v->sx = p.sx;
  if (p.sy) v->sy = p.sy;
  if (p.sz) v->sz = p.sz;
  if (p.bg[0] || p.bg[1] || p.bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = p.bg[i];

  if (p.camera) {
    v->gfsview = false;
    if (strlen(p.camera) >= 4 &&
 !strcmp (&p.camera[strlen(p.camera) - 4], ".gfv")) {
      FILE * fp = fopen (p.camera, "r");
      if (!fp) {
 perror (p.camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (qstderr(), "%s: not a valid gfv file\n", p.camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (p.camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (p.camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, pi, v->quat);
    else if (!strcmp (p.camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, pi/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - pi/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (qstderr(), "view(): unknown camera '%s'\n", p.camera);
      exit (1);
    }
  }
  else if (p.theta || p.phi || p.psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - p.phi, q);
    if (p.relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (p.map)
    v->map = p.map;

  if (p.p1x || p.p1y || p.p2x || p.p2y) {
    float q[4];
    gl_trackball(q, p.p1x, p.p1y, p.p2x, p.p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (p.res)
    v->res = p.res;

  if ((p.width && p.width != v->width) ||
      (p.height && p.height != v->height) ||
      (p.samples && p.samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (p.width) v->width = p.width;
    if (p.height) v->height = p.height;
    if (p.samples) v->samples = p.samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);
    v->fb = framebuffer_new (v->width, v->height);
    init_gl();
  }

  clear();
}







struct _translate {
  float x, y, z;
};

void begin_translate (struct _translate p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glTranslatef (p.x, p.y, p.z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 224 "/home/vinlinux/basilisk/src/draw.h"
struct _mirror {
  coord n;
  double alpha;
};

void begin_mirror (struct _mirror p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  normalize (&p.n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*p.n.x*p.n.x;
  s[1] = - 2.*p.n.x*p.n.y; s[2] = - 2.*p.n.x*p.n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*p.n.y*p.n.y; s[6] = - 2.*p.n.y*p.n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*p.n.z*p.n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = -1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = -1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = -1.; t[11] = 0.;
  t[12] = - 2.*p.n.x*p.alpha;
  t[13] = - 2.*p.n.y*p.alpha;
  t[14] = - 2.*p.n.z*p.alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
}

void end_mirror() {
  end_translate();
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
\
\
  double _r = Delta*0.87;\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (!sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
      < (view)->res) {\
    if (is_active(cell) && is_local(cell)) {\

#line 301

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 308

#line 319 "/home/vinlinux/basilisk/src/draw.h"
#define foreach_visible_plane(view, n1, alpha1)\
coord n = {(n1).x, (n1).y, (n1).z};\
double _alpha = 0.9999999*(alpha1);\
{\
  double norm = sqrt(sq(n.x) + sq(n.y) + sq(n.z));\
  if (!norm)\
    n.z = 1.;\
  else\
    n.x /= norm, n.y /= norm, n.z /= norm, _alpha /= norm;\
}\
glNormal3d (n.x, n.y, n.z);\
foreach_cell() {\
\
  double _r = Delta*0.87, alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;\
  if (fabs(alpha) > 0.87 || !sphere_in_frustum (x, y, z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      sphere_diameter (x, y, z, _r/L0, &(view)->frustum) < (view)->res) {\
    if (is_active(cell) && is_local(cell)) {\

#line 338

#define end_foreach_visible_plane()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 345



static scalar lookup_field (const char * name)
{
  if (name)
    if (all) for (scalar s = *all, *_i82 = all; ((scalar *)&s)->i >= 0; s = *++_i82)
      if (!strcmp (_attribute[s.i].name, name))
 return s;
  return (scalar){-1};
}

static vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    if (all) for (scalar s = *all, *_i83 = all; ((scalar *)&s)->i >= 0; s = *++_i83)
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;
  }
  return (vector){{-1}};
}

static void draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
}

static inline double interp (Point point, coord p, scalar col) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 379 "/home/vinlinux/basilisk/src/draw.h"

  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}
#line 439 "/home/vinlinux/basilisk/src/draw.h"
static void begin_colorized (float fc[3],
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      color c = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = c.r/255.;
      texture[3*i + 1] = c.g/255.;
      texture[3*i + 2] = c.b/255.;
    }
    glTexImage1D (GL_TEXTURE_1D, 0, GL_RGB, 256,0, GL_RGB, GL_FLOAT, texture);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glEnable (GL_TEXTURE_1D);
  }
  glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (GL_TEXTURE_1D);
}
#line 493 "/home/vinlinux/basilisk/src/draw.h"
struct _draw_vof {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
};







static bool cfilter (Point point, scalar c, double cmin)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 514 "/home/vinlinux/basilisk/src/draw.h"

  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    {
#line 517

      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
#line 517

      if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;
#line 517

      if (val(c,0,0,1) >= 1. - cmin1 || val(c,0,0,-1) >= 1. - cmin1)
 return true;}
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    {
#line 523

      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
#line 523

      if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;
#line 523

      if (val(c,0,0,1) <= cmin1 || val(c,0,0,-1) <= cmin1)
 return true;}
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
   { foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 3))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  } end_foreach_neighbor(); }
  return max - min > 0.5;
}
#line 550 "/home/vinlinux/basilisk/src/draw.h"
static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    glVertex3d (p.x, p.y, p.z);
  }
  else
    glVertex3d (x, y, z);
}



bool draw_vof (struct _draw_vof p)
{ trace ("draw_vof", "/home/vinlinux/basilisk/src/draw.h", 563);
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (qstderr(), "draw_vof(): no field named '%s'\n", p.c);
    { bool _ret =  false; end_trace("draw_vof", "/home/vinlinux/basilisk/src/draw.h", 567);  return _ret; }
  }
  vector s = lookup_vector (p.s);

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = lookup_field (p.color); if (col.i < 0) { fprintf (qstderr(), "colorize_args(): no field named '%s'\n", p.color); { bool _ret =  false; end_trace("draw_vof", "/home/vinlinux/basilisk/src/draw.h", 571);  return _ret; } } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { double spread = (p.spread ? p.spread : 5.)*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  double cmin = 1e-3;



  if (_attribute[c.i].prolongation != fraction_refine) {
    _attribute[c.i].prolongation = _attribute[c.i].refine = fraction_refine;
    boundary (((scalar []){c,{-1}}));
  }


  bview * view = draw();
#line 661 "/home/vinlinux/basilisk/src/draw.h"
  double larger =
    p.larger ? p.larger : p.edges || (p.color && !p.linear) ? 1. : 1.1;
  { begin_colorized (p.fc, cmap, !view->vector && p.color && p.linear && col.i >= 0); {
     { foreach_visible (view){

#line 664 "/home/vinlinux/basilisk/src/draw.h"

      if (cfilter (point, c, cmin)) {
 coord n = facet_normal (point, c, s);
 double alpha = plane_alpha (val(c,0,0,0), n);
 coord v[12];
 int m = facets (n, alpha, v, larger);
 if (m > 2) {
   if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); };
   if (view->gfsview)
     glNormal3d (- n.x, - n.y, - n.z);
   else
     glNormal3d (n.x, n.y, n.z);
   glBegin (GL_POLYGON);
   for (int i = 0; i < m; i++) {
     if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, interp (point, v[i], col), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], col); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex3d (view,
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
   }
   glEnd ();
   view->ni++;
 }
      } } end_foreach_visible(); }
  } end_colorized(); }
  if (p.edges) {
    draw_lines (view, p.lc, p.lw);
     { foreach_visible (view){

#line 689 "/home/vinlinux/basilisk/src/draw.h"

      if (cfilter (point, c, cmin)) {
 coord n = facet_normal (point, c, s);
 double alpha = plane_alpha (val(c,0,0,0), n);
 coord v[12];
 int m = facets (n, alpha, v, larger);
 if (m > 2) {
   glBegin (GL_LINE_LOOP);
   for (int i = 0; i < m; i++)
     glvertex3d (view,
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
   glEnd ();
   view->ni++;
 }
      } } end_foreach_visible(); }
    glPopMatrix ();
  }


  { bool _ret =  true; end_trace("draw_vof", "/home/vinlinux/basilisk/src/draw.h", 708);  return _ret; }
 end_trace("draw_vof", "/home/vinlinux/basilisk/src/draw.h", 709); }
#line 722 "/home/vinlinux/basilisk/src/draw.h"
struct _cells {
  coord n;
  double alpha;
  float lc[3], lw;
};


void cells (struct _cells p)
{ trace ("cells", "/home/vinlinux/basilisk/src/draw.h", 730);
  bview * view = draw();
  draw_lines (view, p.lc, p.lw);
#line 744 "/home/vinlinux/basilisk/src/draw.h"
   { foreach_visible_plane (view, p.n, p.alpha){

#line 744 "/home/vinlinux/basilisk/src/draw.h"
 {
    coord v[12];
    int m = facets (n, alpha, v, 1.);
    if (m > 2) {
      glBegin (GL_LINE_LOOP);
      for (int i = 0; i < m; i++)
 glvertex3d (view, x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      glEnd ();
      view->ni++;
    }
  } } end_foreach_visible_plane(); }

  glPopMatrix ();
 end_trace("cells", "/home/vinlinux/basilisk/src/draw.h", 757); }
#line 774 "/home/vinlinux/basilisk/src/draw.h"
struct _squares {
  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];

  coord n;
  double alpha;
};


bool squares (struct _squares p)
{ trace ("squares", "/home/vinlinux/basilisk/src/draw.h", 787);
  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = lookup_field (p.color); if (col.i < 0) { fprintf (qstderr(), "colorize_args(): no field named '%s'\n", p.color); { bool _ret =  false; end_trace("squares", "/home/vinlinux/basilisk/src/draw.h", 788);  return _ret; } } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { double spread = (p.spread ? p.spread : 5.)*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;
  scalar f = col;

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (p.linear) {
    { begin_colorized (p.fc, cmap, !view->vector && p.color && p.linear && col.i >= 0); {
#line 817 "/home/vinlinux/basilisk/src/draw.h"
       { foreach_visible_plane (view, p.n, p.alpha){

#line 817 "/home/vinlinux/basilisk/src/draw.h"

 if (val(f,0,0,0) != nodata) {
   coord v[12];
   int m = facets (n, alpha, v, 1.);
   if (m > 2) {
     coord c = {0,0,0};
     for (int i = 0; i < m; i++)
       {
#line 824

  c.x += v[i].x/m;
#line 824

  c.y += v[i].y/m;
#line 824

  c.z += v[i].z/m;}
     glBegin (GL_TRIANGLE_FAN);
     if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, interp (point, c, f), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, c, f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex3d (view, x + c.x*Delta, y + c.y*Delta, z + c.z*Delta);
     for (int i = 0; i < m; i++) {
       if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, interp (point, v[i], f), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     }
     if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, interp (point, v[0], f), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[0], f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex3d (view,
   x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta);
     glEnd ();
     view->ni++;
   }
 } } end_foreach_visible_plane(); }

    } end_colorized(); }
  }
  else {
#line 858 "/home/vinlinux/basilisk/src/draw.h"
     { foreach_visible_plane (view, p.n, p.alpha){

#line 858 "/home/vinlinux/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 coord v[12];
 int m = facets (n, alpha, v, 1.);
 if (m > 2) {
   if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); };
   glBegin (GL_POLYGON);
   for (int i = 0; i < m; i++)
     glvertex3d (view,
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
   glEnd ();
   view->ni++;
 }
      } } end_foreach_visible_plane(); }

  }
  { bool _ret =  true; end_trace("squares", "/home/vinlinux/basilisk/src/draw.h", 874);  return _ret; }
 end_trace("squares", "/home/vinlinux/basilisk/src/draw.h", 875); }
#line 886 "/home/vinlinux/basilisk/src/draw.h"
struct _box {
  bool notics;
  float lc[3], lw;
};


bool box (struct _box p)
{ trace ("box", "/home/vinlinux/basilisk/src/draw.h", 893);
  bview * view = draw();
  draw_lines (view, p.lc, p.lw);

  float height = 0.5*gl_StrokeHeight();
  float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
  float Z1 = 3 == 2 ? 0. : Z0;
  char label[80];

  glMatrixMode (GL_MODELVIEW);

  if (!p.notics) {
    int nt = 8;
    for (int i = 0; i <= nt; i++) {
      glPushMatrix();
      glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
      glRotatef (-90, 0, 0, 1);
      glScalef (scale, scale, scale);
      sprintf (label, "%g", X0 + i*L0/nt);
      gl_StrokeString (label);
      glPopMatrix();

      glPushMatrix();
      sprintf (label, "%g", Y0 + i*L0/nt);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + width/3.)*scale,
      Y0 + i*L0/nt - height/2.*scale, Z1);
      glScalef (scale, scale, scale);
      gl_StrokeString (label);
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Z0 + i*L0/nt);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + width/3.)*scale,
      Y0, Z0 + i*L0/nt + height/2.*scale);
      glRotatef (-90, 1, 0, 0);
      glScalef (scale, scale, scale);
      gl_StrokeString (label);
      glPopMatrix();

    }

    glPushMatrix();
    sprintf (label, "%g", X0 + L0/2.);
    length = gl_StrokeLength (label);
    glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
    glScalef (2.*scale, 2.*scale, 2.*scale);
    gl_StrokeString ("X");
    glPopMatrix();


    glPushMatrix();
    sprintf (label, "%g", Y0 + L0/2.);
    length = gl_StrokeLength (label);
    glTranslatef (X0 - (length + 4.*width)*scale,
    Y0 + L0/2. - height*scale, Z1);
    glScalef (2.*scale, 2.*scale, 2.*scale);
    gl_StrokeString ("Y");
    glPopMatrix();


    glPushMatrix();
    sprintf (label, "%g", Z0 + L0/2.);
    length = gl_StrokeLength (label);
    glTranslatef (X0 - (length + 4.*width)*scale,
    Y0, Z0 + L0/2. + height*scale);
    glRotatef (-90, 1, 0, 0);
    glScalef (2.*scale, 2.*scale, 2.*scale);
    gl_StrokeString ("Z");
    glPopMatrix();

  }
#line 979 "/home/vinlinux/basilisk/src/draw.h"
   { foreach_level (0){

#line 979 "/home/vinlinux/basilisk/src/draw.h"
 {
    for (int i = -1; i <= 1; i += 2) {
      glBegin (GL_LINE_LOOP);
      glvertex3d (view, x - Delta/2., y - Delta/2., z + i*Delta/2.);
      glvertex3d (view, x + Delta/2., y - Delta/2., z + i*Delta/2.);
      glvertex3d (view, x + Delta/2., y + Delta/2., z + i*Delta/2.);
      glvertex3d (view, x - Delta/2., y + Delta/2., z + i*Delta/2.);
      glEnd ();
      view->ni++;
      glBegin (GL_LINES);
      for (int j = -1; j <= 1; j += 2) {
 glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z - Delta/2.);
 glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z + Delta/2.);
      }
      glEnd ();
      view->ni++;
    }
  } } end_foreach_level(); }


  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  { bool _ret =  true; end_trace("box", "/home/vinlinux/basilisk/src/draw.h", 1001);  return _ret; }
 end_trace("box", "/home/vinlinux/basilisk/src/draw.h", 1002); }
#line 1014 "/home/vinlinux/basilisk/src/draw.h"
struct _isosurface {
  char * f;
  double v;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];
};


bool isosurface (struct _isosurface p)
{ trace ("isosurface", "/home/vinlinux/basilisk/src/draw.h", 1027);

  scalar f = lookup_field (p.f);
  if (f.i < 0) {
    fprintf (qstderr(), "isosurface(): no field named '%s'\n", p.f);
    { bool _ret =  false; end_trace("isosurface", "/home/vinlinux/basilisk/src/draw.h", 1032);  return _ret; }
  }

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = lookup_field (p.color); if (col.i < 0) { fprintf (qstderr(), "colorize_args(): no field named '%s'\n", p.color); { bool _ret =  false; end_trace("isosurface", "/home/vinlinux/basilisk/src/draw.h", 1035);  return _ret; } } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { double spread = (p.spread ? p.spread : 5.)*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  scalar v= new_vertex_scalar("v");
   { foreach_vertex(){

#line 1038 "/home/vinlinux/basilisk/src/draw.h"

    val(v,0,0,0) = (val(f,0,0,0) + val(f,-1,0,0) + val(f,0,-1,0) + val(f,-1,-1,0) +
    val(f,0,0,-1) + val(f,-1,0,-1) + val(f,0,-1,-1) + val(f,-1,-1,-1))/8.; } end_foreach_vertex(); }

  vector n= new_vector("n");
   { foreach(){

#line 1043 "/home/vinlinux/basilisk/src/draw.h"

    {
#line 1044

      val(n.x,0,0,0) = (val(f,1,0,0) - val(f,-1,0,0))/(2.*Delta);
#line 1044

      val(n.y,0,0,0) = (val(f,0,1,0) - val(f,0,-1,0))/(2.*Delta);
#line 1044

      val(n.z,0,0,0) = (val(f,0,0,1) - val(f,0,0,-1))/(2.*Delta);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{n.x,n.y,n.z},{{-1},{-1},{-1}}}));

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  { begin_colorized (p.fc, cmap, !view->vector && p.color && p.linear && col.i >= 0); {
     { foreach_visible (view){

#line 1051 "/home/vinlinux/basilisk/src/draw.h"
 {
      double val[8] = {
 val(v,0,0,0), val(v,1,0,0), val(v,1,0,1), val(v,0,0,1),
 val(v,0,1,0), val(v,1,1,0), val(v,1,1,1), val(v,0,1,1)
      };
      double t[5][3][3];
      int nt = polygonize (val, p.v, t);
      for (int i = 0; i < nt; i++) {
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); };
 glBegin (GL_POLYGON);
 for (int j = 0; j < 3; j++) {
   coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;
   {
#line 1063

     np.x = interp (point, v, n.x);
#line 1063

     np.y = interp (point, v, n.y);
#line 1063

     np.z = interp (point, v, n.z);}
   glNormal3d (np.x, np.y, np.z);
   if (p.color && p.linear && col.i >= 0) { if (view->vector) { color b = colormap_color (cmap, interp (point, v, col), p.min, p.max); glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v, col); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   glvertex3d (view, x + v.x*Delta, y + v.y*Delta, z + v.z*Delta);
 }
 glEnd ();
 view->ni++;
      }
    } } end_foreach_visible(); }
  } end_colorized(); }

  { bool _ret =  true; delete (((scalar []){n.x,n.y,n.z,v,{-1}}));  end_trace("isosurface", "/home/vinlinux/basilisk/src/draw.h", 1075);  return _ret; }
 delete (((scalar []){n.x,n.y,n.z,v,{-1}}));  end_trace("isosurface", "/home/vinlinux/basilisk/src/draw.h", 1076); }
#line 1086 "/home/vinlinux/basilisk/src/draw.h"
struct _travelling {
  double start, end;
  float tx, ty, quat[4], fov;
};




void travelling (struct _travelling p)
{
  static float tx, ty, quat[4], fov;
  static double told = -1.;
  if (told < p.start && t >= p.start) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  if (t >= p.start && t <= p.end)
    view ((struct _view_set){.tx = (!p.tx ? tx : ((t - p.start)*(p.tx) + (p.end - t)*(tx))/(p.end - p.start)), .ty = (!p.ty ? ty : ((t - p.start)*(p.ty) + (p.end - t)*(ty))/(p.end - p.start)),
   .fov = (!p.fov ? fov : ((t - p.start)*(p.fov) + (p.end - t)*(fov))/(p.end - p.start)),
   .quat = {(!p.quat[0] ? quat[0] : ((t - p.start)*(p.quat[0]) + (p.end - t)*(quat[0]))/(p.end - p.start)), (!p.quat[1] ? quat[1] : ((t - p.start)*(p.quat[1]) + (p.end - t)*(quat[1]))/(p.end - p.start)),
           (!p.quat[2] ? quat[2] : ((t - p.start)*(p.quat[2]) + (p.end - t)*(quat[2]))/(p.end - p.start)), (!p.quat[3] ? quat[3] : ((t - p.start)*(p.quat[3]) + (p.end - t)*(quat[3]))/(p.end - p.start))}});
  if (told < p.end && t >= p.end) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  told = t;
}
#line 1133 "/home/vinlinux/basilisk/src/draw.h"
struct _draw_string {
  char * str;
  int pos;
  float size;
  float lc[3], lw;
};


bool draw_string (struct _draw_string p)
{ trace ("draw_string", "/home/vinlinux/basilisk/src/draw.h", 1142);
  bview * view = draw();

  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glColor3f (p.lc[0], p.lc[1], p.lc[2]);
  glLineWidth (view->samples*(p.lw > 0. ? p.lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!p.size)
    p.size = 40;
  float hscale = 2./(p.size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (p.pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (p.pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (p.pos == 2)
    glTranslatef (1. - strlen(p.str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(p.str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (p.str);

  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();

  { bool _ret =  true; end_trace("draw_string", "/home/vinlinux/basilisk/src/draw.h", 1177);  return _ret; }
 end_trace("draw_string", "/home/vinlinux/basilisk/src/draw.h", 1178); }
#line 340 "/home/vinlinux/basilisk/src/view.h"
#line 360 "/home/vinlinux/basilisk/src/view.h"
struct _load {
  FILE * fp;
  char * file;
  Array * buf;
  Array * history;
};

bool load (struct _load p);
#line 411 "/home/vinlinux/basilisk/src/view.h"
struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw;
  int sort, options;

  Array * history;
  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}

static void redraw_feedback (struct _save * p)
{
  bview * view = p->view ? p->view : get_view();
  assert (p->history);
  if (p->history->len) {
    float res = view->res;
    view->res = 0.;
    view->vector = true;
    redraw();


    int list = glGenLists (1);
    glNewList (list, GL_COMPILE);
    load ((struct _load){.buf = p->history});
    glEndList();
    glCallList (list);
    glFinish ();
    glDeleteLists (list, 1);
    enable_fpe (FE_DIVBYZERO|FE_INVALID);
    view->active = false;
    view->vector = false;
    view->res = res;
  }
}




bool save (struct _save p)
{ trace ("save", "/home/vinlinux/basilisk/src/view.h", 459);
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
 dot = s;
 s = strchr (s + 1, '.');
      }
      if (dot)
 p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();

  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    { bool _ret =  true; end_trace("save", "/home/vinlinux/basilisk/src/view.h", 486);  return _ret; }
  }

  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    { bool _ret =  false; end_trace("save", "/home/vinlinux/basilisk/src/view.h", 491);  return _ret; }
  }
  if (!p.fp)
    p.fp = qstdout();

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (p.format, "bv")) {
    assert (p.history);
    fprintf (p.fp,
      "view (fov = %g, quat = {%g,%g,%g,%g}, "
      "tx = %g, ty = %g, "
      "bg = {%g,%g,%g}, "
      "width = %d, height = %d, samples = %d"
      ");\n",
      view->fov,
      view->quat[0], view->quat[1], view->quat[2], view->quat[3],
      view->tx, view->ty,
      view->bg[0], view->bg[1], view->bg[2],
      view->width/view->samples, view->height/view->samples,
      view->samples);
    fwrite (p.history->p, 1, p.history->len, p.fp);
  }

  else if (!strcmp (p.format, "gnu") ||
    !strcmp (p.format, "obj") ||
    !strcmp (p.format, "kml")) {
    int format = (!strcmp (p.format, "gnu") ? FEEDBACK_GNU :
    !strcmp (p.format, "obj") ? FEEDBACK_OBJ :
    !strcmp (p.format, "kml") ? FEEDBACK_KML :
    -1);
    unsigned buffsize = 1 << 24;
    bool done = false;
    while (!done && buffsize <= (1 << 28)) {
      float * f = gl_feedback_begin (buffsize);
      redraw_feedback (&p);
      done = gl_feedback_end (f, p.fp, format);
      buffsize *= 2;
    }
    if (!done)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }

  else if (!strcmp (p.format, "ps") ||
    !strcmp (p.format, "eps") ||
    !strcmp (p.format, "tex") ||
    !strcmp (p.format, "pdf") ||
    !strcmp (p.format, "svg") ||
    !strcmp (p.format, "pgf")) {
    GLint format = (!strcmp (p.format, "ps") ? GL2PS_PS :
      !strcmp (p.format, "eps") ? GL2PS_EPS :
      !strcmp (p.format, "tex") ? GL2PS_TEX :
      !strcmp (p.format, "pdf") ? GL2PS_PDF :
      !strcmp (p.format, "svg") ? GL2PS_SVG :
      !strcmp (p.format, "pgf") ? GL2PS_PGF :
      -1);
    GLint state = GL2PS_OVERFLOW;
    GLint sort = p.sort ? p.sort : GL2PS_SIMPLE_SORT;
    GLint options = p.options ? p.options : (GL2PS_SIMPLE_LINE_OFFSET |
          GL2PS_SILENT |
          GL2PS_BEST_ROOT |
          GL2PS_OCCLUSION_CULL |
          GL2PS_USE_CURRENT_VIEWPORT |
          GL2PS_TIGHT_BOUNDING_BOX);
    unsigned buffsize = 1 << 24;
    while (state == GL2PS_OVERFLOW && buffsize <= (1 << 28)) {
      gl2psBeginPage ("", "bview",
        NULL,
        format, sort, options,
        GL_RGBA, 0, NULL,
        0, 0, 0,
        buffsize, p.fp, "");
      redraw_feedback (&p);
      disable_fpe (FE_DIVBYZERO|FE_INVALID);
      state = gl2psEndPage();
      enable_fpe (FE_DIVBYZERO|FE_INVALID);
      buffsize *= 2;
    }
    if (state == GL2PS_OVERFLOW)
      fprintf (ferr, "save(): error: exceeded maximum feedback buffer size\n");
  }

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    { bool _ret =  false; end_trace("save", "/home/vinlinux/basilisk/src/view.h", 584);  return _ret; }
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

  { bool _ret =  true; end_trace("save", "/home/vinlinux/basilisk/src/view.h", 591);  return _ret; }
 end_trace("save", "/home/vinlinux/basilisk/src/view.h", 592); }







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}

static void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  if (all) for (scalar s = *all, *_i84 = all; ((scalar *)&s)->i >= 0; s = *++_i84)
    fprintf (ferr, " %s", _attribute[s.i].name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  if (all) for (scalar s = *all, *_i85 = all; ((scalar *)&s)->i >= 0; s = *++_i85) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

static void draw_append (char * buf, Array * history, FILE * interactive)
{
  if (interactive) {
    if (history->len)
      load ((struct _load){.buf = history});
    save ((struct _save){.fp = interactive});
  }
  array_append (history, buf, strlen(buf)*sizeof(char));
}






#line 1 "draw_get.h"
#line 1 "/home/vinlinux/basilisk/src/draw_get.h"


bool _view_set_get (struct _view_set * p) {
  Params params[] = {
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"fov", pfloat, &p->fov},
    {"quat", pfloat, p->quat, 4},
    {"sx", pfloat, &p->sx},
    {"sy", pfloat, &p->sy},
    {"sz", pfloat, &p->sz},
    {"width", punsigned, &p->width},
    {"height", punsigned, &p->height},
    {"samples", punsigned, &p->samples},
    {"bg", pfloat, p->bg, 3},
    {"theta", pfloat, &p->theta},
    {"phi", pfloat, &p->phi},
    {"psi", pfloat, &p->psi},
    {"relative", pbool, &p->relative},
    {"res", pfloat, &p->res},
    {"camera", pstring, &p->camera},
    {"p1x", pfloat, &p->p1x},
    {"p1y", pfloat, &p->p1y},
    {"p2x", pfloat, &p->p2x},
    {"p2y", pfloat, &p->p2y},
    {NULL}
  };
  return parse_params (params);
}

bool _translate_get (struct _translate * p) {
  Params params[] = {
    {"x", pfloat, &p->x},
    {"y", pfloat, &p->y},
    {"z", pfloat, &p->z},
    {NULL}
  };
  return parse_params (params);
}

bool _mirror_get (struct _mirror * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_vof_get (struct _draw_vof * p) {
  Params params[] = {
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _cells_get (struct _cells * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _squares_get (struct _squares * p) {
  Params params[] = {
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _box_get (struct _box * p) {
  Params params[] = {
    {"notics", pbool, &p->notics},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _isosurface_get (struct _isosurface * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"v", pdouble, &p->v},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {NULL}
  };
  return parse_params (params);
}

bool _travelling_get (struct _travelling * p) {
  Params params[] = {
    {"start", pdouble, &p->start},
    {"end", pdouble, &p->end},
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"quat", pfloat, p->quat, 4},
    {"fov", pfloat, &p->fov},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_string_get (struct _draw_string * p) {
  Params params[] = {
    {"str", pstring, &p->str},
    {"pos", pint, &p->pos},
    {"size", pfloat, &p->size},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}
#line 647 "/home/vinlinux/basilisk/src/view.h"

static bool process_line (char * line, Array * history, FILE * interactive)
{
  if (line[0] == '\0')
    return true;
  char * buf = pstrdup (line,__func__,__FILE__,__LINE__);
  char * s = strtok (remove_blanks (line), "(");
  if (!s) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      if (!restore ((struct Dump){.file = file, .list = all}))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();

 if (history->len && load ((struct _load){.buf = history}) && interactive)
   save ((struct _save){.fp = interactive});
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ((struct Dump){.file = file});
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ((struct OutputGfs){.file = file, .list = all});
      restriction (all);
      fields_stats();
      clear();

      if (history->len && load ((struct _load){.buf = history}) && interactive)
 save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ((struct _save){.file = file, .history = history});
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file && load ((struct _load){.file = file, .history = history}) && interactive) {
      load ((struct _load){.buf = history});
      save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    _cells_get (&p);
    cells (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    _draw_vof_get (&p);
    if (draw_vof (p))
      draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_translate")) {
    end_translate();
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "end_mirror")) {
    end_mirror();
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    _squares_get (&p);
    squares (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    _isosurface_get (&p);
    isosurface (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    _draw_string_get (&p);
    draw_string (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "display")) {
    if (interactive && history->len && load ((struct _load){.buf = history}))
      save ((struct _save){.fp = interactive});
  }

  else if (!strcmp (s, "clear")) {
    clear();
    if (interactive)
      save ((struct _save){.fp = interactive});
    history->len = 0;
  }

  else if (!strcmp (s, "show")) {
    if (interactive && history->len)
      save ((struct _save){.fp = ferr, .format = "bv", .history = history});
  }

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    _box_get (&p);
    box (p);
    draw_append (buf, history, interactive);
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
    if (p.width || p.height || p.samples) {

      if (history->len && load ((struct _load){.buf = history}) && p.samples && interactive)
 save ((struct _save){.fp = interactive});
    }
  }

  else if (!strcmp (s, "quit")) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return false;
  }

  else if (s[0] != '\n')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  pfree (buf,__func__,__FILE__,__LINE__);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  Array * history = array_new();
  if (p.fp) {
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line, history, NULL));
  }
  else if (p.buf) {
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line, history, NULL);
      }
    }
  }
  if (p.history)
    array_append (p.history, history->p, history->len);
  array_free (history);

  return true;
}
#line 9 "main.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 26 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/vinlinux/basilisk/src/run.h"
#line 9 "/home/vinlinux/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 12 "/home/vinlinux/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/vinlinux/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/vinlinux/basilisk/src/run.h", 37); }
#line 27 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/vinlinux/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _dtmax = dtmax; 
#line 6

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.x,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.y,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/vinlinux/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta*val_cm(cm,0,0,0)/fabs(val(u.z,0,0,0));
      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 10
 end_foreach_face(); }OMP(omp critical) if (_dtmax < dtmax) dtmax = _dtmax;
mpi_all_reduce_double (dtmax, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 10
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 28 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/vinlinux/basilisk/src/bcg.h"
#line 11 "/home/vinlinux/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




   { 
if (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





      double vn = val(uf.y,i,0,0)/val_fm_y(fm.y,i,0,0) + val(uf.y,i,1,0)/val_fm_y(fm.y,i,1,0);
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.z,i,0,0)/val_fm_z(fm.z,i,0,0) + val(uf.z,i,0,1)/val_fm_z(fm.z,i,0,1);
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





      double vn = val(uf.z,0,i,0)/val_fm_z(fm.z,0,i,0) + val(uf.z,0,i,1)/val_fm_z(fm.z,0,i,1);
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.x,0,i,0)/val_fm_x(fm.x,0,i,0) + val(uf.x,1,i,0)/val_fm_x(fm.x,1,i,0);
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/vinlinux/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





      double vn = val(uf.x,0,0,i)/val_fm_x(fm.x,0,0,i) + val(uf.x,1,0,i)/val_fm_x(fm.x,1,0,i);
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(4.*Delta);


      double wn = val(uf.y,0,0,i)/val_fm_y(fm.y,0,0,i) + val(uf.y,0,1,i)/val_fm_y(fm.y,0,1,i);
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(4.*Delta);


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 55
 end_foreach_face(); } }





  boundary_flux (((vector []){{flux.x,flux.y,flux.z},{{-1},{-1},{-1}}}));
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc) {
    scalar zero= new_const_scalar("zero", 8,  0.);
    if (p.tracers) for (scalar s = *p.tracers, *_i86 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i86)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  scalar * _i2 = p.tracers; scalar * _i3 = lsrc; if (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i2, src = *++_i3) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);
     { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 94
foreach(){

#line 94 "/home/vinlinux/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 94
foreach(){

#line 94 "/home/vinlinux/basilisk/src/bcg.h"

      {
#line 95

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 95

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }
   delete (((scalar []){flux.x,flux.y,flux.z,{-1}})); }
  boundary (p.tracers);

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 29 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 1 "./viscosity.h"
#line 1 "/home/vinlinux/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/vinlinux/basilisk/src/poisson.h"
#line 32 "/home/vinlinux/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 54 "/home/vinlinux/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i87 = da; ((scalar *)&s)->i >= 0; s = *++_i87)
   val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 63 "/home/vinlinux/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i88 = da; ((scalar *)&s)->i >= 0; s = *++_i88)
   val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { foreach(){

#line 81 "/home/vinlinux/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i4 = a; scalar * _i5 = da; if (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i4, ds = *++_i5)
      val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
  boundary (a);
}
#line 99 "/home/vinlinux/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
} mgstats;
#line 120 "/home/vinlinux/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    if (p.a) for (scalar s = *p.a, *_i89 = p.a; ((scalar *)&s)->i >= 0; s = *++_i89) {
      scalar r = new_scalar("r");
      res = list_append (res, r);
    }






  for (int b = 0; b < nboundary; b++)
    if (da) for (scalar s = *da, *_i90 = da; ((scalar *)&s)->i >= 0; s = *++_i90)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; 
#line 160
foreach (){

#line 160 "/home/vinlinux/basilisk/src/poisson.h"

    if (p.b) for (scalar s = *p.b, *_i91 = p.b; ((scalar *)&s)->i >= 0; s = *++_i91)
      _sum += val(s,0,0,0); } end_foreach();OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 162
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > TOLERANCE);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data, s.nrelax, 0, grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);







    if (s.resa > TOLERANCE) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }
    resb = s.resa;
  }




  if (s.resa > TOLERANCE)
    fprintf (ferr,
      "WARNING: convergence not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n",
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 237 "/home/vinlinux/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax;
  scalar * res;
};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 272 "/home/vinlinux/basilisk/src/poisson.h"
  scalar c = a;






   { 
if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/vinlinux/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/vinlinux/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/vinlinux/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 279
foreach_level_or_leaf (l){

#line 279 "/home/vinlinux/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 281
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 281
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 281
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
    val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 304 "/home/vinlinux/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
  double maxres = 0.;


  vector g= new_face_vector("g");
   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 321
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.x,0,0,0) = val_alpha_x(alpha.x,0,0,0)*(val(a,0,0,0) - val(a,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.y,0,0,0) = val_alpha_y(alpha.y,0,0,0)*(val(a,0,0,0) - val(a,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 321
{

#line 321 "/home/vinlinux/basilisk/src/poisson.h"

    val(g.z,0,0,0) = val_alpha_z(alpha.z,0,0,0)*(val(a,0,0,0) - val(a,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 322
 end_foreach_face(); } }
  boundary_flux (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 324

if (!is_constant(lambda)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#line 324
foreach (){

#line 324 "/home/vinlinux/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#line 324
foreach (){

#line 324 "/home/vinlinux/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 326

      val(res,0,0,0) += (val(g.x,0,0,0) - val(g.x,1,0,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.y,0,0,0) - val(g.y,0,1,0))/Delta;
#line 326

      val(res,0,0,0) += (val(g.z,0,0,0) - val(g.z,0,0,1))/Delta;}
    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 330
 }
#line 342 "/home/vinlinux/basilisk/src/poisson.h"
  boundary (resl);
  { double _ret =  maxres; delete (((scalar []){g.x,g.y,g.z,{-1}}));  return _ret; }
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }
#line 355 "/home/vinlinux/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i) {
    vector alpha= new_const_vector("alpha", 9, (double []) {1.,1.,1.});
    p.alpha = alpha;
  }
  if (!p.lambda.i) {
    scalar lambda= new_const_scalar("lambda", 12,  0.);
    p.lambda = lambda;
  }




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;
  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax, &p, p.nrelax, p.res});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 416 "/home/vinlinux/basilisk/src/poisson.h"
struct Project {
  vector u;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/vinlinux/basilisk/src/poisson.h", 426);
  vector u = q.u;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { foreach(){

#line 439 "/home/vinlinux/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 441

      val(div,0,0,0) += val(u.x,1,0,0) - val(u.x,0,0,0);
#line 441

      val(div,0,0,0) += val(u.y,0,1,0) - val(u.y,0,0,0);
#line 441

      val(div,0,0,0) += val(u.z,0,0,1) - val(u.z,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 455 "/home/vinlinux/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 461
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 461
{

#line 461 "/home/vinlinux/basilisk/src/poisson.h"

    val(u.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 462
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/vinlinux/basilisk/src/poisson.h", 465);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/vinlinux/basilisk/src/poisson.h", 466); }
#line 2 "/home/vinlinux/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/vinlinux/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/vinlinux/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/vinlinux/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/vinlinux/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/vinlinux/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); } }
#line 85 "/home/vinlinux/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;


  {
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.x,0,0,0) = 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
      (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
      (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.x)) {
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.x)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
      (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
      (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);}
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_y()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.y,0,0,0) = 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_z()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.z,0,0,0) = val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
      (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
      (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.y)) {
#undef val_mu_y
#define val_mu_y(a,k,i,j) val(a,k,i,j)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_z
#define val_mu_z(a,k,i,j) val(a,k,i,j)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,k,i,j)
#undef val_mu_x
#define val_mu_x(a,k,i,j) val(a,k,i,j)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,k,i,j)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,k,i,j)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.y)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.y.i -_NVARMAX], _constant[mu.z.i - _NVARMAX], _constant[mu.x.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_y
#define val_mu_y(a,k,i,j) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,k,i,j) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,k,i,j) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_x()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
      (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
      (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,k,i,j) val(a,k,i,j)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,k,i,j)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,k,i,j)
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,k,i,j) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);}
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }
#line 98
 {
    vector taux= new_face_vector("taux");
     { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 100
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_z()) {
#line 100
{

#line 100 "/home/vinlinux/basilisk/src/viscosity.h"

      val(taux.z,0,0,0) = 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 101
 end_foreach_face(); } }

       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 103
foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_x()) {
#line 103
{

#line 103 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.x,0,0,0) = val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
      (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
      (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 106
 end_foreach_face(); } }


       { 
if (!is_constant(mu.z)) {
#undef val_mu_z
#define val_mu_z(a,j,k,i) val(a,j,k,i)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_x
#define val_mu_x(a,j,k,i) val(a,j,k,i)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,j,k,i)
#undef val_mu_y
#define val_mu_y(a,j,k,i) val(a,j,k,i)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,j,k,i)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,j,k,i)
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); }
if (is_constant(mu.z)) {
const struct { double x, y, z; } _const_mu = {_constant[mu.z.i -_NVARMAX], _constant[mu.x.i - _NVARMAX], _constant[mu.y.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_z
#define val_mu_z(a,j,k,i) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#undef val_mu_x
#define val_mu_x(a,j,k,i) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,j,k,i) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#line 109
foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_y()) {
#line 109
{

#line 109 "/home/vinlinux/basilisk/src/viscosity.h"

 val(taux.y,0,0,0) = val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
      (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
      (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)/Delta; }  }}  end_foreach_face_generic()
#line 112
 end_foreach_face(); } }

    boundary_flux (((vector []){{taux.x,taux.y,taux.z},{{-1},{-1},{-1}}}));
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 115

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,j,k,i) val(a,j,k,i)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,j,k,i)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,j,k,i)
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,j,k,i) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 115
foreach (){

#line 115 "/home/vinlinux/basilisk/src/viscosity.h"
 {
      double d = 0.;
      {
#line 117

 d += val(taux.z,0,0,1) - val(taux.z,0,0,0);
#line 117

 d += val(taux.x,1,0,0) - val(taux.x,0,0,0);
#line 117

 d += val(taux.y,0,1,0) - val(taux.y,0,0,0);}
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) + dt/val_rho(rho,0,0,0)*d/Delta;
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 122
 }
   delete (((scalar []){taux.x,taux.y,taux.z,{-1}})); }}
  boundary (resl);
#line 153 "/home/vinlinux/basilisk/src/viscosity.h"
  return maxres;
}



mgstats viscosity (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
   { foreach(){

#line 161 "/home/vinlinux/basilisk/src/viscosity.h"

    {
#line 162

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 162

      val(r.y,0,0,0) = val(u.y,0,0,0);
#line 162

      val(r.z,0,0,0) = val(u.z,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,mu.z,rho,{-1}}));

  { mgstats _ret =  mg_solve ((struct MGSolve){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res}); delete (((scalar []){r.x,r.y,r.z,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}})); }

mgstats viscosity_explicit (struct Viscosity p)
{
  vector u = p.u, r= new_vector("r");
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}), &p);
   { foreach(){

#line 178 "/home/vinlinux/basilisk/src/viscosity.h"

    {
#line 179

      val(u.x,0,0,0) += val(r.x,0,0,0);
#line 179

      val(u.y,0,0,0) += val(r.y,0,0,0);
#line 179

      val(u.z,0,0,0) += val(r.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  { mgstats _ret =  mg; delete (((scalar []){r.x,r.y,r.z,{-1}}));  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}})); }
#line 30 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
#line 39 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2},{3}}, g= {{4},{5},{6}};
scalar pf= {7};
vector uf= {{8},{9},{10}};
#line 65 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, a = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, alpha = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar rho = {(_NVARMAX + 6)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 79 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; } 







static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; } 


static void _set_boundary4 (void) { _attribute[p.i].boundary[front] = _boundary4; _attribute[p.i].boundary_homogeneous[front] = _boundary4_homogeneous; } 
static void _set_boundary5 (void) { _attribute[p.i].boundary[back] = _boundary5; _attribute[p.i].boundary_homogeneous[back] = _boundary5_homogeneous; } 






static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 100); 
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 119
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 120
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 119
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 119
{

#line 119 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 120
 end_foreach_face(); } }
    boundary ((scalar *)((vector []){{alpha.x,alpha.y,alpha.z},{{-1},{-1},{-1}}}));
  }






  _attribute[uf.x.i].refine = refine_face_solenoidal;

 end_trace("defaults", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 131); } return 0; } 





double dtmax;

static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 139); 
{
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 144
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 143
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 143
{

#line 143 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.; }  }}  end_foreach_face_generic()
#line 144
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));




  event ("properties");





  dtmax = DT;
  event ("stability");
 end_trace("init", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 158); } return 0; } 
#line 167 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 167);  dtmax = DT; end_trace("set_dtmax", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 167);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 169);  {
  dt = dtnext (timestep (uf, dtmax));
 end_trace("stability", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 171); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 179); ; end_trace("vof", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 179);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 180); ; end_trace("tracer_advection", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 180);  return 0; } 
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion (const int i, const double t, Event * _ev) { trace ("tracer_diffusion", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 181); ; end_trace("tracer_diffusion", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 181);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 188);  {
  boundary (((scalar []){alpha.x,alpha.y,alpha.z,mu.x,mu.y,mu.z,rho,{-1}}));
 end_trace("properties", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 190); } return 0; } 
#line 202 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
  {
#line 205
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 205
 {
    scalar s = new_scalar("s");
    du.y = s;
  }
#line 205
 {
    scalar s = new_scalar("s");
    du.z = s;
  }}

  if (_attribute[u.x.i].gradient)
     { foreach(){

#line 211 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      {
#line 212

        val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
#line 212

        val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
#line 212

        val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;}; } end_foreach(); }
  else
     { foreach(){

#line 215 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

      {
#line 216

        val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
#line 216

        val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
#line 216

        val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));

  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 221
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 234
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 221
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);


      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);


      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 221
{

#line 221 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);


      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 234
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));

  delete ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));
}
#line 249 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 249); 
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 256); } return 0; } 







static void correction (double dt)
{
   { foreach(){

#line 266 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    {
#line 267

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 267

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);
#line 267

      val(u.z,0,0,0) += dt*val(g.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
}
#line 279 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 279); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }






  vector af = a;
  trash (((vector []){{uf.x,uf.y,uf.z},{af.x,af.y,af.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 294
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 298
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 294
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(val(u.x,0,0,0) + val(u.x,-1,0,0))/2.;
    if (!is_constant(af.x))
      val(af.x,0,0,0) = 0.;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(val(u.y,0,0,0) + val(u.y,0,-1,0))/2.;
    if (!is_constant(af.y))
      val(af.y,0,0,0) = 0.;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 294
{

#line 294 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
 {
    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(val(u.z,0,0,0) + val(u.z,0,0,-1))/2.;
    if (!is_constant(af.z))
      val(af.z,0,0,0) = 0.;
  } }  }}  end_foreach_face_generic()
#line 298
 end_foreach_face(); } }
 end_trace("viscous_term", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 299); } return 0; } 
#line 314 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 314); 
{
  boundary ((scalar *)((vector []){{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 317
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) += dt*val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) += dt*val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 317
{

#line 317 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) += dt*val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0); }  }}  end_foreach_face_generic()
#line 318
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
 end_trace("acceleration", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 320); } return 0; } 
#line 329 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"
static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 329); 
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});





  vector gf= new_face_vector("gf");
   { 
if (!is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && !is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && !is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (!is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); }
if (is_constant(a.x) && is_constant(alpha.x) && is_constant(fm.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 338
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)/val_fm_x(fm.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)/val_fm_y(fm.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 338
{

#line 338 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)/val_fm_z(fm.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 339
 end_foreach_face(); } }
  boundary_flux (((vector []){{gf.x,gf.y,gf.z},{{-1},{-1},{-1}}}));





  trash (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 347 "/home/vinlinux/basilisk/src/navier-stokes/centered.h"

    {
#line 348

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/2.;
#line 348

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/2.;
#line 348

      val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/2.;}; } end_foreach(); }
  boundary ((scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




  correction (dt);
 delete (((scalar []){gf.x,gf.y,gf.z,{-1}}));  end_trace("projection", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 356); } return 0; } 







static int adapt_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int adapt (const int i, const double t, Event * _ev) { trace ("adapt", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 364);  {
  event ("properties");
 end_trace("adapt", "/home/vinlinux/basilisk/src/navier-stokes/centered.h", 366); } return 0; } 
#line 10 "main.c"
#line 1 "tracer.h"
#line 1 "/home/vinlinux/basilisk/src/tracer.h"
#line 15 "/home/vinlinux/basilisk/src/tracer.h"
extern scalar * tracers;
extern vector uf;
extern double dt;







static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/vinlinux/basilisk/src/tracer.h", 25);  {
  if (tracers) for (scalar s = *tracers, *_i92 = tracers; ((scalar *)&s)->i >= 0; s = *++_i92) {
    _attribute[s.i].refine = refine_linear;
    _attribute[s.i].restriction = restriction_volume_average;
  }
 end_trace("defaults_0", "/home/vinlinux/basilisk/src/tracer.h", 30); } return 0; } 





#line 1 "bcg.h"
#line 37 "/home/vinlinux/basilisk/src/tracer.h"

static int tracer_advection_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection_0 (const int i, const double t, Event * _ev) { trace ("tracer_advection_0", "/home/vinlinux/basilisk/src/tracer.h", 38);  {
  advection ((struct Advection){tracers, uf, dt});
 end_trace("tracer_advection_0", "/home/vinlinux/basilisk/src/tracer.h", 40); } return 0; } 




static int tracer_diffusion_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion_0 (const int i, const double t, Event * _ev) { trace ("tracer_diffusion_0", "/home/vinlinux/basilisk/src/tracer.h", 45); ; end_trace("tracer_diffusion_0", "/home/vinlinux/basilisk/src/tracer.h", 45);  return 0; } 
#line 11 "main.c"
#line 1 "diffusion.h"
#line 1 "/home/vinlinux/basilisk/src/diffusion.h"
#line 25 "/home/vinlinux/basilisk/src/diffusion.h"
#line 1 "poisson.h"
#line 26 "/home/vinlinux/basilisk/src/diffusion.h"
#line 39 "/home/vinlinux/basilisk/src/diffusion.h"
struct Diffusion {

  scalar f;
  double dt;

  vector D;
  scalar r, beta;
  scalar theta;
};


mgstats diffusion (struct Diffusion p)
{ trace ("diffusion", "/home/vinlinux/basilisk/src/diffusion.h", 51);




  if (p.dt == 0.) {
    mgstats s = {0};
    { mgstats _ret =  s; end_trace("diffusion", "/home/vinlinux/basilisk/src/diffusion.h", 58);  return _ret; }
  }




  scalar f = p.f, r = (p.r).i ? (p.r) : new_scalar("r");




  scalar idt= new_const_scalar("idt", 13,  - 1./p.dt);
   scalar theta_idt = p.theta.i ? p.theta : idt;

  if (p.theta.i) {
    scalar theta_idt = p.theta;
     { foreach(){

#line 74 "/home/vinlinux/basilisk/src/diffusion.h"

      val(theta_idt,0,0,0) *= _val_constant(idt,0,0,0); } end_foreach(); }
  }




  if (p.r.i)
     { 
if (!is_constant(theta_idt)) {
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) val(a,i,j,k)
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) fine(a,i,j,k)
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) coarse(a,i,j,k)
#line 82
foreach(){

#line 82 "/home/vinlinux/basilisk/src/diffusion.h"

      val(r,0,0,0) = val_theta_idt(theta_idt,0,0,0)*val(f,0,0,0) - val(r,0,0,0); } end_foreach(); }
if (is_constant(theta_idt)) {
const double _const_theta_idt = _constant[theta_idt.i -_NVARMAX];
NOT_UNUSED(_const_theta_idt);
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) _const_theta_idt
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) _const_theta_idt
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) _const_theta_idt
#line 82
foreach(){

#line 82 "/home/vinlinux/basilisk/src/diffusion.h"

      val(r,0,0,0) = val_theta_idt(theta_idt,0,0,0)*val(f,0,0,0) - val(r,0,0,0); } end_foreach(); } }
  else
     { 
if (!is_constant(theta_idt)) {
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) val(a,i,j,k)
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) fine(a,i,j,k)
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) coarse(a,i,j,k)
#line 85
foreach(){

#line 85 "/home/vinlinux/basilisk/src/diffusion.h"

      val(r,0,0,0) = val_theta_idt(theta_idt,0,0,0)*val(f,0,0,0); } end_foreach(); }
if (is_constant(theta_idt)) {
const double _const_theta_idt = _constant[theta_idt.i -_NVARMAX];
NOT_UNUSED(_const_theta_idt);
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) _const_theta_idt
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) _const_theta_idt
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) _const_theta_idt
#line 85
foreach(){

#line 85 "/home/vinlinux/basilisk/src/diffusion.h"

      val(r,0,0,0) = val_theta_idt(theta_idt,0,0,0)*val(f,0,0,0); } end_foreach(); } }




  scalar lambda = theta_idt;
  if (p.beta.i) {
    scalar beta = p.beta;
     { 
if (!is_constant(theta_idt)) {
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) val(a,i,j,k)
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) fine(a,i,j,k)
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) coarse(a,i,j,k)
#line 94
foreach(){

#line 94 "/home/vinlinux/basilisk/src/diffusion.h"

      val(beta,0,0,0) += val_theta_idt(theta_idt,0,0,0); } end_foreach(); }
if (is_constant(theta_idt)) {
const double _const_theta_idt = _constant[theta_idt.i -_NVARMAX];
NOT_UNUSED(_const_theta_idt);
#undef val_theta_idt
#define val_theta_idt(a,i,j,k) _const_theta_idt
#undef fine_theta_idt
#define fine_theta_idt(a,i,j,k) _const_theta_idt
#undef coarse_theta_idt
#define coarse_theta_idt(a,i,j,k) _const_theta_idt
#line 94
foreach(){

#line 94 "/home/vinlinux/basilisk/src/diffusion.h"

      val(beta,0,0,0) += val_theta_idt(theta_idt,0,0,0); } end_foreach(); } }
    lambda = beta;
  }
  boundary (((scalar []){lambda,{-1}}));




  { mgstats _ret =  poisson ((struct Poisson){f, r, p.D, lambda}); { if (!(p.r).i) delete (((scalar []){r,{-1}})); }  end_trace("diffusion", "/home/vinlinux/basilisk/src/diffusion.h", 103);  return _ret; }
 { if (!(p.r).i) delete (((scalar []){r,{-1}})); }  end_trace("diffusion", "/home/vinlinux/basilisk/src/diffusion.h", 104); }
#line 12 "main.c"


int minlevel, maxlevel;
double meps, eps;
double TEND = 1000;

char sim_ID[] = "krab";
char sim_var[] = "Trot";

#line 1 "physics.h"
#line 1 "./physics.h"

#line 1 "SGS.h"
#line 1 "/home/vinlinux/basilisk/src/SGS.h"






#line 1 "vreman.h"
#line 1 "/home/vinlinux/basilisk/src/vreman.h"
#line 13 "/home/vinlinux/basilisk/src/vreman.h"
void eddyviscosity(double Cs, vector u, double molv, scalar Evis){
  double d1v1, d2v1, d3v1, d1v2, d2v2, d3v2, d1v3, d2v3, d3v3;
  double b11, b12, b13, b22, b23, b33;
  double abeta, bbeta;
   { foreach(){

#line 17 "/home/vinlinux/basilisk/src/vreman.h"
{
    d1v1 = (val(u.x,1,0,0) - val(u.x,-1,0,0))/2/Delta;
    d2v1 = (val(u.x,0,1,0) - val(u.x,0,-1,0))/2/Delta;
    d3v1 = (val(u.x,0,0,1) - val(u.x,0,0,-1))/2/Delta;
    d1v2 = (val(u.y,1,0,0) - val(u.y,-1,0,0))/2/Delta;
    d2v2= (val(u.y,0,1,0) - val(u.y,0,-1,0))/2/Delta;
    d3v2= (val(u.y,0,0,1) - val(u.y,0,0,-1))/2/Delta;
    d1v3= (val(u.z,1,0,0) - val(u.z,-1,0,0))/2/Delta;
    d2v3= (val(u.z,0,1,0) - val(u.z,0,-1,0))/2/Delta;
    d3v3= (val(u.z,0,0,1) - val(u.z,0,0,-1))/2/Delta;
    b11 = sq(Delta)*(d1v1*d1v1 + d2v1*d2v1 + d3v1*d3v1);
    b12 = sq(Delta)*(d1v1*d1v2 + d2v1*d2v2 + d3v1*d3v2);
    b13 = sq(Delta)*(d1v1*d1v3 + d2v1*d2v3 + d3v1*d3v3);
    b22 = sq(Delta)*(d1v2*d1v2 + d2v2*d2v2 + d3v2*d3v2);
    b23 = sq(Delta)*(d1v2*d1v3 + d2v2*d2v3 + d3v2*d3v3);
    b33 = sq(Delta)*(d1v3*d1v3 + d2v3*d2v3 + d3v3*d3v3);
    abeta = sq(d1v1) + sq(d2v1) + sq(d3v1)+
            sq(d1v2) + sq(d2v2) + sq(d3v2)+
            sq(d1v3) + sq(d2v3) + sq(d3v3);
    bbeta = b11*b22 - sq(b12) + b11*b33 - sq(b13) + b22*b33 - sq(b23);

    val(Evis,0,0,0) = (abeta > 10E-5 && bbeta > (abeta/10E6))? 2.5*sq(Cs)*sqrt(bbeta/(abeta)) + molv: molv;
  } } end_foreach(); }
}
#line 8 "/home/vinlinux/basilisk/src/SGS.h"



vector Km= {{11},{12},{13}},Kh= {{14},{15},{16}};
 vector Pr;
scalar Evis= {17};
double molvis;
double Csmag;
scalar * tracers;




static inline void Evisprol(Point point,scalar s){ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 21 "/home/vinlinux/basilisk/src/SGS.h"

   { foreach_child()
    val(Evis,0,0,0)=bilinear(point,Evis)/4.; end_foreach_child(); }
}
static inline void Evisres(Point point,scalar s){ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 25 "/home/vinlinux/basilisk/src/SGS.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/2.;
}




static int defaults_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i=0);   *ip = i; *tp = t;   return ret; } static int defaults_1 (const int i, const double t, Event * _ev) { trace ("defaults_1", "/home/vinlinux/basilisk/src/SGS.h", 35); {
  if (3!=3)
    fprintf(qstdout(),"Warning %dD grid. The used formulations only make sense for 3D turbulence simulations\n",3);
  mu=Kh;
  Pr=unityf;
  molvis=0.;
  Csmag=0.12;
  _attribute[Evis.i].prolongation=Evisprol;
  _attribute[Evis.i].restriction=Evisres;
#line 53 "/home/vinlinux/basilisk/src/SGS.h"
 end_trace("defaults_1", "/home/vinlinux/basilisk/src/SGS.h", 53); } return 0; } 



static int Eddyvis_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int Eddyvis (const int i, const double t, Event * _ev) { trace ("Eddyvis", "/home/vinlinux/basilisk/src/SGS.h", 57); {
  eddyviscosity(Csmag,u,molvis,Evis);
  boundary(((scalar []){Evis,{-1}}));
   { 
if (!is_constant(Pr.x)) {
#undef val_Pr_x
#define val_Pr_x(a,i,j,k) val(a,i,j,k)
#undef fine_Pr_x
#define fine_Pr_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_Pr_x
#define coarse_Pr_x(a,i,j,k) coarse(a,i,j,k)
#undef val_Pr_y
#define val_Pr_y(a,i,j,k) val(a,i,j,k)
#undef fine_Pr_y
#define fine_Pr_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_Pr_y
#define coarse_Pr_y(a,i,j,k) coarse(a,i,j,k)
#undef val_Pr_z
#define val_Pr_z(a,i,j,k) val(a,i,j,k)
#undef fine_Pr_z
#define fine_Pr_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_Pr_z
#define coarse_Pr_z(a,i,j,k) coarse(a,i,j,k)
#line 60
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.x,0,0,0)=(val(Evis,0,0,0)+val(Evis,-1,0,0))/2;
    val(Kh.x,0,0,0)=(val_Pr_x(Pr.x,0,0,0)*(val(Km.x,0,0,0)-molvis))+molvis;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.y,0,0,0)=(val(Evis,0,0,0)+val(Evis,0,-1,0))/2;
    val(Kh.y,0,0,0)=(val_Pr_y(Pr.y,0,0,0)*(val(Km.y,0,0,0)-molvis))+molvis;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.z,0,0,0)=(val(Evis,0,0,0)+val(Evis,0,0,-1))/2;
    val(Kh.z,0,0,0)=(val_Pr_z(Pr.z,0,0,0)*(val(Km.z,0,0,0)-molvis))+molvis;
  } }  }}  end_foreach_face_generic()
#line 63
 end_foreach_face(); }
if (is_constant(Pr.x)) {
const struct { double x, y, z; } _const_Pr = {_constant[Pr.x.i -_NVARMAX], _constant[Pr.y.i - _NVARMAX], _constant[Pr.z.i - _NVARMAX]};
NOT_UNUSED(_const_Pr);
#undef val_Pr_x
#define val_Pr_x(a,i,j,k) _const_Pr.x
#undef fine_Pr_x
#define fine_Pr_x(a,i,j,k) _const_Pr.x
#undef coarse_Pr_x
#define coarse_Pr_x(a,i,j,k) _const_Pr.x
#undef val_Pr_y
#define val_Pr_y(a,i,j,k) _const_Pr.y
#undef fine_Pr_y
#define fine_Pr_y(a,i,j,k) _const_Pr.y
#undef coarse_Pr_y
#define coarse_Pr_y(a,i,j,k) _const_Pr.y
#undef val_Pr_z
#define val_Pr_z(a,i,j,k) _const_Pr.z
#undef fine_Pr_z
#define fine_Pr_z(a,i,j,k) _const_Pr.z
#undef coarse_Pr_z
#define coarse_Pr_z(a,i,j,k) _const_Pr.z
#line 60
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.x,0,0,0)=(val(Evis,0,0,0)+val(Evis,-1,0,0))/2;
    val(Kh.x,0,0,0)=(val_Pr_x(Pr.x,0,0,0)*(val(Km.x,0,0,0)-molvis))+molvis;
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.y,0,0,0)=(val(Evis,0,0,0)+val(Evis,0,-1,0))/2;
    val(Kh.y,0,0,0)=(val_Pr_y(Pr.y,0,0,0)*(val(Km.y,0,0,0)-molvis))+molvis;
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 60
{

#line 60 "/home/vinlinux/basilisk/src/SGS.h"
{
    val(Km.z,0,0,0)=(val(Evis,0,0,0)+val(Evis,0,0,-1))/2;
    val(Kh.z,0,0,0)=(val_Pr_z(Pr.z,0,0,0)*(val(Km.z,0,0,0)-molvis))+molvis;
  } }  }}  end_foreach_face_generic()
#line 63
 end_foreach_face(); } }

  boundary_flux(((vector []){{Km.x,Km.y,Km.z},{Kh.x,Kh.y,Kh.z},{{-1},{-1},{-1}}}));
 end_trace("Eddyvis", "/home/vinlinux/basilisk/src/SGS.h", 66); } return 0; } 
#line 3 "./physics.h"
#line 18 "./physics.h"
scalar b= {18};
scalar * tracers = ((scalar []){{18},{-1}});

double crho = 1.;

vector av= {{19},{20},{21}};
struct sCase def;

struct sCase {
 double wind;
 double wphi;
};

void init_physics(){
  def.wind = (-1.);
        def.wphi = 0.;

 _attribute[b.i].nodump = false;

        _attribute[u.x.i].boundary[bottom] = _boundary6; _attribute[u.x.i].boundary_homogeneous[bottom] = _boundary6_homogeneous;
        _attribute[u.y.i].boundary[bottom] = _boundary7; _attribute[u.y.i].boundary_homogeneous[bottom] = _boundary7_homogeneous;
        _attribute[u.x.i].boundary[top] = _boundary8; _attribute[u.x.i].boundary_homogeneous[top] = _boundary8_homogeneous;
        _attribute[u.y.i].boundary[top] = _boundary9; _attribute[u.y.i].boundary_homogeneous[top] = _boundary9_homogeneous;

        tree_periodic(left);

 _attribute[b.i].boundary[bottom] = _boundary10; _attribute[b.i].boundary_homogeneous[bottom] = _boundary10_homogeneous;
 _attribute[b.i].boundary[top] = _boundary11; _attribute[b.i].boundary_homogeneous[top] = _boundary11_homogeneous;


  _attribute[u.z.i].boundary[bottom] = _boundary12; _attribute[u.z.i].boundary_homogeneous[bottom] = _boundary12_homogeneous;
                _attribute[u.z.i].boundary[top] = _boundary13; _attribute[u.z.i].boundary_homogeneous[top] = _boundary13_homogeneous;
         _attribute[u.y.i].boundary[bottom] = _boundary14; _attribute[u.y.i].boundary_homogeneous[bottom] = _boundary14_homogeneous;
  _attribute[u.y.i].boundary[top] = _boundary15; _attribute[u.y.i].boundary_homogeneous[top] = _boundary15_homogeneous;

  tree_periodic(front);

  { foreach(){

#line 55 "./physics.h"
 {
  val(b,0,0,0) = 9.81*(.5)*y/273.;
          val(u.x,0,0,0) = (-1.);
 } } end_foreach(); }
}


static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "./physics.h", 62); {
  { foreach_face_generic() { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 63
{

#line 63 "./physics.h"
{
  val(av.y,0,0,0) = (val(b,0,0,0) + val(b,0,-1,0))/2.;
 } }  }}  end_foreach_face_generic()
#line 65
 end_foreach_face(); }
 end_trace("acceleration_0", "./physics.h", 66); } return 0; } 

static int inflow_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int inflow (const int i, const double t, Event * _ev) { trace ("inflow", "./physics.h", 68); {
    double sides = 25;
    double relaxtime = dt/25.;
     { foreach(){

#line 71 "./physics.h"
{
 if((x < sides || x > L0-sides ||
     z < sides || z > L0-sides ||
      y < sides || y > L0-sides )) {
     val(u.x,0,0,0) = val(u.x,0,0,0) + ((-1.)-val(u.x,0,0,0))*relaxtime;
      val(b,0,0,0) = val(b,0,0,0) + (9.81*(.5)*y/273. - val(b,0,0,0))*relaxtime;
     val(u.y,0,0,0) = val(u.y,0,0,0) - val(u.y,0,0,0)*relaxtime;
            val(u.z,0,0,0) = val(u.z,0,0,0) - val(u.z,0,0,0)*relaxtime;
 }
    } } end_foreach(); }
 end_trace("inflow", "./physics.h", 81); } return 0; } 

mgstats mgb;

static int tracer_diffusion_1_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion_1 (const int i, const double t, Event * _ev) { trace ("tracer_diffusion_1", "./physics.h", 85); {
    mgb = diffusion((struct Diffusion){b, dt, mu});
 end_trace("tracer_diffusion_1", "./physics.h", 87); } return 0; } 
#line 22 "main.c"
#line 1 "fan.h"
#line 1 "./fan.h"


struct sRotor rot;
scalar fan= {22};

struct sRotor {
    bool fan;
    double start;
    double stop;
    double rampT;
    double P, Prho;
    double R, W, A, V;
    double diaVol;
    double x0, y0, z0;
    double xt, yt, zt;
    double theta, phi;
    double thetat, phit;
    double Work;
    double cu;
    bool rotate;
    coord nf, nr;
};


void init_rotor();
void rotor_update();
void rotor_coord();
void rotor_forcing();


void init_rotor() {
    rot.Work = 0.;
    if(!rot.start)
     rot.start = 0.;
    if(!rot.stop)
     rot.stop = 1E10;
    if(!rot.rampT)
     rot.rampT = 10;
    if(!rot.R)
 rot.R = 2.;
    if(!rot.W)
 rot.W = 0.3;
    if(!rot.Prho)
      rot.Prho = 3000.;
    if(!rot.x0)
     rot.x0 = L0/2.;
    if(!rot.y0)
 rot.y0 = L0/2.;
    if(!rot.z0){



            rot.z0 = L0/2.;

    }
    if(!rot.theta)
     rot.theta = 97*M_PI/180.;
    if(!rot.phi)
     rot.phi = 0.*M_PI/180.;

    if(rot.rotate) {
 if(!rot.phit){
     rot.phit = -2*M_PI/240;
  }
 if(!rot.xt){
     rot.xt = 0;
  }
        if(!rot.yt){
     rot.yt = 0;
  }
 if(!rot.zt){
     rot.zt = 0;
  }
 if(!rot.thetat){
     rot.thetat = 0.;
  }

 } else {
        rot.xt = 0;
        rot.yt = 0;
        rot.zt = 0;
        rot.thetat = 0.;
        rot.phit = 0.;
    }

    rotor_update();
}


static int forcing_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int forcing (const int i, const double t, Event * _ev) { trace ("forcing", "./fan.h", 90);  {
    if(rot.fan && t>rot.start && t<rot.stop) {
 rotor_coord();
 rotor_forcing();
    }
 end_trace("forcing", "./fan.h", 95); } return 0; } 


static int rotate_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int rotate (const int i, const double t, Event * _ev) { trace ("rotate", "./fan.h", 98);  {
    if(t>rot.start && t<rot.stop){

        rot.x0 += rot.xt;
        rot.y0 += rot.yt;
        rot.z0 += rot.zt;


        rot.theta += dt*rot.thetat;
        rot.phi += dt*rot.phit;

        rotor_update();
    }
 end_trace("rotate", "./fan.h", 111); } return 0; } 


void rotor_update() {
     rot.nf.x = sin(rot.theta)*cos(rot.phi);
 rot.nf.z = sin(rot.theta)*sin(rot.phi);
 rot.nf.y = cos(rot.theta);

 rot.nr.x = sin(rot.theta)*cos(rot.phi);
    rot.nr.z = sin(rot.theta)*sin(rot.phi);
     rot.nr.y = cos(rot.theta);




  rot.A = sq(rot.R)*M_PI;






  rot.V = 4.*M_PI*pow(rot.R,3.)/3. -
   2*M_PI*pow(rot.R-rot.W/2., 2.)/3.*(2*rot.R + rot.W/2.);


 rot.P = rot.V*rot.Prho;
     rot.cu = pow(3*rot.Prho*rot.W*crho, 1./3.);
}


void rotor_coord() {
    scalar sph= new_scalar("sph"), plnu= new_scalar("plnu"), plnd= new_scalar("plnd");
    do { scalar phi= new_vertex_scalar("phi");  { foreach_vertex(){

#line 144 "./fan.h"
 val(phi,0,0,0) = -sq((x - rot.x0)) - sq((y - rot.y0)) - sq((z - rot.z0)) + sq(rot.R); } end_foreach_vertex(); } fractions ((struct Fractions){phi, sph});  delete (((scalar []){phi,{-1}})); } while(0);
    do { scalar phi= new_vertex_scalar("phi");  { foreach_vertex(){

#line 145 "./fan.h"
 val(phi,0,0,0) = rot.nr.x*(x - rot.x0) + rot.nr.y*(y - rot.y0) + rot.nr.z*(z - rot.z0) + rot.W/2.; } end_foreach_vertex(); } fractions ((struct Fractions){phi, plnu});  delete (((scalar []){phi,{-1}})); } while(0);
    do { scalar phi= new_vertex_scalar("phi");  { foreach_vertex(){

#line 146 "./fan.h"
 val(phi,0,0,0) = -rot.nr.x*(x - rot.x0) - rot.nr.y*(y - rot.y0) - rot.nr.z*(z - rot.z0) + rot.W/2.; } end_foreach_vertex(); } fractions ((struct Fractions){phi, plnd});  delete (((scalar []){phi,{-1}})); } while(0);

     { foreach (){

#line 148 "./fan.h"
 {
        val(fan,0,0,0) = val(sph,0,0,0)*val(plnu,0,0,0)*val(plnd,0,0,0);
    } } end_foreach(); }
    boundary(((scalar []){fan,{-1}}));
 delete (((scalar []){plnd,plnu,sph,{-1}})); }



void rotor_forcing(){
    double tempW = 0.;
    double tempVol = 0.;
    double w, wsgn, damp, usgn, utemp, corrP;

     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _tempVol = tempVol; 
#line 161

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 161
foreach(){

#line 161 "./fan.h"
{
         _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 161
foreach(){

#line 161 "./fan.h"
{
         _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
    } } end_foreach(); }OMP(omp critical) tempVol += _tempVol;
mpi_all_reduce_double (tempVol, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 163
 }
    rot.diaVol = 1*tempVol;

    Point point = locate((struct _locate){rot.x0, rot.y0, rot.z0});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 166 "./fan.h"


    if(val(fan,0,0,0) == 0.){
 {
#line 169
{

      wsgn = sign(rot.nf.x*val(u.x,0,0,0)) + (sign(rot.nf.x*val(u.x,0,0,0)) == 0)*sign(rot.nf.x);
      damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
      w = wsgn*damp*sq(rot.nf.x)*(2.)*rot.P*dt;
      tempW += fabs(w);


       utemp = sq(val(u.x,0,0,0)) + w;
       usgn = 1.*(val(u.x,0,0,0) >= 0)*(utemp > 0) +
               -1.*(val(u.x,0,0,0) >= 0)*(utemp < 0) +
              1.*(val(u.x,0,0,0) < 0)*(utemp < 0) +
             -1.*(val(u.x,0,0,0) < 0)*(utemp > 0);


              val(u.x,0,0,0) = usgn*min(sqrt(fabs(utemp)), sq(damp)*1.5*rot.cu);

         }
#line 169
{

      wsgn = sign(rot.nf.y*val(u.y,0,0,0)) + (sign(rot.nf.y*val(u.y,0,0,0)) == 0)*sign(rot.nf.y);
      damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
      w = wsgn*damp*sq(rot.nf.y)*(2.)*rot.P*dt;
      tempW += fabs(w);


       utemp = sq(val(u.y,0,0,0)) + w;
       usgn = 1.*(val(u.y,0,0,0) >= 0)*(utemp > 0) +
               -1.*(val(u.y,0,0,0) >= 0)*(utemp < 0) +
              1.*(val(u.y,0,0,0) < 0)*(utemp < 0) +
             -1.*(val(u.y,0,0,0) < 0)*(utemp > 0);


              val(u.y,0,0,0) = usgn*min(sqrt(fabs(utemp)), sq(damp)*1.5*rot.cu);

         }
#line 169
{

      wsgn = sign(rot.nf.z*val(u.z,0,0,0)) + (sign(rot.nf.z*val(u.z,0,0,0)) == 0)*sign(rot.nf.z);
      damp = rot.rampT + rot.start > t ? (t-rot.start)/rot.rampT : 1.;
      w = wsgn*damp*sq(rot.nf.z)*(2.)*rot.P*dt;
      tempW += fabs(w);


       utemp = sq(val(u.z,0,0,0)) + w;
       usgn = 1.*(val(u.z,0,0,0) >= 0)*(utemp > 0) +
               -1.*(val(u.z,0,0,0) >= 0)*(utemp < 0) +
              1.*(val(u.z,0,0,0) < 0)*(utemp < 0) +
             -1.*(val(u.z,0,0,0) < 0)*(utemp > 0);


              val(u.z,0,0,0) = usgn*min(sqrt(fabs(utemp)), sq(damp)*1.5*rot.cu);

         }}
     } else {

      { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _tempW = tempW; 
#line 189

if (!is_constant(rho)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 189
foreach(){

#line 189 "./fan.h"
 {
         if(val(fan,0,0,0) > 0.) {
                {
#line 191
 {

  wsgn = sign(rot.nf.x*val(u.x,0,0,0)) + (sign(rot.nf.x*val(u.x,0,0,0)) == 0)*sign(rot.nf.x);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.x)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.x,0,0,0)) + w;
  usgn = 1.*(val(u.x,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.x,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.x,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.x,0,0,0) < 0)*(utemp > 0);


  val(u.x,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }
#line 191
 {

  wsgn = sign(rot.nf.y*val(u.y,0,0,0)) + (sign(rot.nf.y*val(u.y,0,0,0)) == 0)*sign(rot.nf.y);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.y)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.y,0,0,0)) + w;
  usgn = 1.*(val(u.y,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.y,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.y,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.y,0,0,0) < 0)*(utemp > 0);


  val(u.y,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }
#line 191
 {

  wsgn = sign(rot.nf.z*val(u.z,0,0,0)) + (sign(rot.nf.z*val(u.z,0,0,0)) == 0)*sign(rot.nf.z);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.z)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.z,0,0,0)) + w;
  usgn = 1.*(val(u.z,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.z,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.z,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.z,0,0,0) < 0)*(utemp > 0);


  val(u.z,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }}
 }
    } } end_foreach(); }
if (is_constant(rho)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 189
foreach(){

#line 189 "./fan.h"
 {
         if(val(fan,0,0,0) > 0.) {
                {
#line 191
 {

  wsgn = sign(rot.nf.x*val(u.x,0,0,0)) + (sign(rot.nf.x*val(u.x,0,0,0)) == 0)*sign(rot.nf.x);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.x)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.x,0,0,0)) + w;
  usgn = 1.*(val(u.x,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.x,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.x,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.x,0,0,0) < 0)*(utemp > 0);


  val(u.x,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }
#line 191
 {

  wsgn = sign(rot.nf.y*val(u.y,0,0,0)) + (sign(rot.nf.y*val(u.y,0,0,0)) == 0)*sign(rot.nf.y);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.y)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.y,0,0,0)) + w;
  usgn = 1.*(val(u.y,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.y,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.y,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.y,0,0,0) < 0)*(utemp > 0);


  val(u.y,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }
#line 191
 {

  wsgn = sign(rot.nf.z*val(u.z,0,0,0)) + (sign(rot.nf.z*val(u.z,0,0,0)) == 0)*sign(rot.nf.z);
  damp = rot.rampT > t ? t/rot.rampT : 1.;
  corrP = rot.diaVol > 0. ? rot.V/rot.diaVol : 1.;
  w = wsgn*val(fan,0,0,0)*damp*sq(rot.nf.z)*(2./val_rho(rho,0,0,0))*(corrP*rot.P/rot.V)*dt;
  _tempW += fabs(w);


  utemp = sq(val(u.z,0,0,0)) + w;
  usgn = 1.*(val(u.z,0,0,0) >= 0)*(utemp > 0) +
        -1.*(val(u.z,0,0,0) >= 0)*(utemp < 0) +
     1.*(val(u.z,0,0,0) < 0)*(utemp < 0) +
    -1.*(val(u.z,0,0,0) < 0)*(utemp > 0);


  val(u.z,0,0,0) = usgn*min(sqrt(fabs(utemp)), damp*1.5*rot.cu);
      }}
 }
    } } end_foreach(); }OMP(omp critical) tempW += _tempW;
mpi_all_reduce_double (tempW, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 210
 }
    }
    rot.Work += tempW;
}
#line 23 "main.c"
#line 1 "diagnostics.h"
#line 1 "./diagnostics.h"


#line 1 "lambda2.h"
#line 1 "/home/vinlinux/basilisk/src/lambda2.h"
static void eigsrt (double d[3],
      double v[3][3])
{
  int k, j, i;
  double p;

  for (i = 0; i < 3 - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < 3; j++)
      if (d[j] >= p)
 p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < 3; j++) {
 p = v[j][i];
 v[j][i] = v[j][k];
 v[j][k] = p;
      }
    }
  }
}
#line 37 "/home/vinlinux/basilisk/src/lambda2.h"
void eigenvalues (double a[3][3],
    double d[3],
    double v[3][3])
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

  for (ip = 0; ip < 3; ip++) {
    for (iq = 0; iq < 3; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < 3; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++)
 sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(3*3);
    else
      tresh = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++) {
 g = 100.0*fabs (a[ip][iq]);
 if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
     fabs(d[iq]) + g == fabs(d[iq]))
   a[ip][iq] = 0.0;
 else if (fabs (a[ip][iq]) > tresh) {
   h = d[iq] - d[ip];
   if (fabs(h) + g == fabs(h))
     t = a[ip][iq]/h;
   else {
     theta = 0.5*h/a[ip][iq];
     t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
     if (theta < 0.0) t = -t;
   }
   c = 1.0/sqrt (1 + t*t);
   s = t*c;
   tau = s/(1.0 + c);
   h = t*a[ip][iq];
   z[ip] -= h;
   z[iq] += h;
   d[ip] -= h;
   d[iq] += h;
   a[ip][iq] = 0.0;
   for (j = 0; j <= ip - 1; j++)
     { g=a[j][ip];h=a[j][iq];a[j][ip]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = ip + 1; j <= iq - 1; j++)
     { g=a[ip][j];h=a[j][iq];a[ip][j]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = iq + 1; j < 3; j++)
     { g=a[ip][j];h=a[iq][j];a[ip][j]=g-s*(h+g*tau);a[iq][j]=h+s*(g-h*tau);};
   for (j = 0; j < 3; j++)
     { g=v[j][ip];h=v[j][iq];v[j][ip]=g-s*(h+g*tau);v[j][iq]=h+s*(g-h*tau);};
 }
      }
    }
    for (ip = 0; ip < 3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      fprintf (qstderr(), "%10.3g ", a[i][j]);
    fprintf (qstderr(), "\n");
  }
  assert (false);
}

void lambda2 (const vector u, scalar l2)
{
   { foreach(){

#line 121 "/home/vinlinux/basilisk/src/lambda2.h"
 {
    double JJ[3][3];
    scalar s = u.x;
    int i = 0;
    {
#line 125

      JJ[0][i++] = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 125

      JJ[0][i++] = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 125

      JJ[0][i++] = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    s = u.y; i = 0;
    {
#line 128

      JJ[1][i++] = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 128

      JJ[1][i++] = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 128

      JJ[1][i++] = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    s = u.z; i = 0;
    {
#line 131

      JJ[2][i++] = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
#line 131

      JJ[2][i++] = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
#line 131

      JJ[2][i++] = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);}
    double S2O2[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
 S2O2[i][j] = 0.;
 for (int k = 0; k < 3; k++)
   S2O2[i][j] += JJ[i][k]*JJ[k][j] + JJ[k][i]*JJ[j][k];
      }
    double lambda[3], ev[3][3];
    eigenvalues (S2O2, lambda, ev);
    val(l2,0,0,0) = lambda[1]/2.;
  } } end_foreach(); }
  boundary (((scalar []){l2,{-1}}));
}
#line 4 "./diagnostics.h"

#line 1 "output_slices.h"
#line 1 "/home/vinlinux/basilisk/src/output_slices.h"
#line 39 "/home/vinlinux/basilisk/src/output_slices.h"
struct sOutputSlice {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  coord plane;
};


void output_slice (struct sOutputSlice p)
{ trace ("output_slice", "/home/vinlinux/basilisk/src/output_slices.h", 49);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = qstdout();
  if (!p.plane.x) p.plane.x = 1.;
  if (!p.plane.y) p.plane.y = 1.;
  if (!p.plane.z) p.plane.z = 0.;
p.n++;

  int len = list_len(p.list);
  double ** field = (double **) matrix_new (p.n, p.n, len*sizeof(double));
  double Delta = 0.999999*L0/(p.n - 1);

  for (int i = 0; i < p.n; i++) {
    double varCoord1 = Delta*i;
    bool varX = !(p.plane.x < 1.);
    double x = (!varX ? p.plane.x*L0 : varCoord1) + X0;

    for (int j = 0; j < p.n; j++) {
      double varCoord2 = Delta*j;
      double y = (varX ? (p.plane.y < 1. ? p.plane.y*L0 : varCoord2) : varCoord1) + Y0;
      double z = (p.plane.z < 1. ? p.plane.z*L0 : varCoord2) + Z0;
      if (p.linear) {
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i93 = p.list; ((scalar *)&s)->i >= 0; s = *++_i93)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y, z});
      }
      else {
 Point point = locate ((struct _locate){x, y, z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 77 "/home/vinlinux/basilisk/src/output_slices.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i94 = p.list; ((scalar *)&s)->i >= 0; s = *++_i94)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if _MPI
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

    fprintf (p.fp, "x=%g\ty=%g\tz=%g\tn=%d\tlen=%d\n", p.plane.x*L0, p.plane.y*L0, p.plane.z*L0, p.n, len);
    int k = 0;
    if (p.list) for (scalar s = *p.list, *_i95 = p.list; ((scalar *)&s)->i >= 0; s = *++_i95) {
    fprintf (p.fp, "%s\n", _attribute[s.i].name);
    for (int i = 0; i < p.n; i++) {
      for (int j = 0; j < p.n; j++) {
        fprintf (p.fp, "%g\t", (float) field[i][len*j + k]);
      }
      fputc ('\n', p.fp);
    }
    k++;
    }
    fflush (p.fp);
  }
#if _MPI
  else
    MPI_Reduce (field[0], NULL, len*p.n*p.n, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_slice", "/home/vinlinux/basilisk/src/output_slices.h", 112); }
#line 6 "./diagnostics.h"
#line 1 "profile5c.h"
#line 1 "/home/vinlinux/basilisk/src/profile5c.h"
#line 33 "/home/vinlinux/basilisk/src/profile5c.h"
struct prof {
  scalar * list;
  char * fname;
  double ym;
  double h;
  double rf;
  FILE * fp;
  int n;
};

double average_over_yp(scalar * list, double * v, double yp){
  int m = 0;
  double a = 0;

   { foreach_leaf(){

#line 47 "/home/vinlinux/basilisk/src/profile5c.h"




    {
      if ((fabs(y-yp) < (Delta/2))){
 m++;
 double b = 1./Delta;
 {
#line 55

   b *= Delta;
#line 55

   b *= Delta;
#line 55

   b *= Delta;}
 a += b;
 int k = 0;
 if (list) for (scalar s = *list, *_i96 = list; ((scalar *)&s)->i >= 0; s = *++_i96)
   v[k++] += point.level >= 0 ? interpolate ((struct _interpolate){s, x, yp, z})*b : 0;
      }
    } } end_foreach_leaf(); }
#line 73 "/home/vinlinux/basilisk/src/profile5c.h"
  int g = 0;
  if (list) for (scalar s = *list, *_i97 = list; ((scalar *)&s)->i >= 0; s = *++_i97)
    v[g++] /= a;



  return sqrt(a/(double)m);

}

void field_profile(struct prof p){



  if(!p.list)
    p.list = all;
  if (!p.ym)
    p.ym = Y0 + 0.9999999*(L0 / (double)(1 << ((grid->maxdepth) + 1)));
  if (!p.h)
    p.h = Y0 + L0 - 0.9999999*(L0 / (double)(1 << ((grid->maxdepth) + 1)));
  if (!p.rf)
    p.rf = 1;
  if (!p.fname && !p.fp)
    p.fp = qstdout();
  double dzn;
  if (p.n)
    dzn = 0.9999999*(p.h - p.ym)/((double)p.n - 1.);
  int len = list_len(p.list);
  boundary(p.list);
  FILE * fp = p.fp;
  char * file = p.fname;



  if (pid()==0){
    if (file && (fp = fopen (file, "w")) == NULL) {
      perror (file);
      exit (1);
    }
    assert (fp);



    fprintf(fp,"y\t");
    if (p.list) for (scalar s = *p.list, *_i98 = p.list; ((scalar *)&s)->i >= 0; s = *++_i98)
      fprintf(fp,"%s\t",_attribute[s.i].name);
    fprintf(fp,"\n");
  }
  double yp = p.ym;



  while (yp <= p.h){
    double aver[len];
    memset(&aver, 0, sizeof(aver[0])*len);
    double dz = average_over_yp(p.list, aver, yp);
    if (pid() == 0){
      int k = 0;
      if (p.list) for (scalar s = *p.list, *_i99 = p.list; ((scalar *)&s)->i >= 0; s = *++_i99){
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
    if (fp != qstdout())
      fclose(fp);
  }
}
#line 7 "./diagnostics.h"
#line 1 "equi_data.h"
#line 1 "/home/vinlinux/basilisk/src/equi_data.h"
#line 36 "/home/vinlinux/basilisk/src/equi_data.h"
double *equifield = NULL;

int equi_diag (scalar s, int lvl, int diagi){
    int len = 1;
    int n = 1<<lvl;
    if(!equifield) {
        equifield = pcalloc(len*n*n*n, sizeof(double),__func__,__FILE__,__LINE__);
    }

    double dDelta = 0.9999999*L0/n;

    restriction(((scalar []){s,{-1}}));
    boundary(((scalar []){s,{-1}}));


     { foreach_level_or_leaf(lvl){

#line 51 "/home/vinlinux/basilisk/src/equi_data.h"
 {
   if(level < lvl){
     int pn = 1<<(lvl-level);
     for(int i = 0; i < pn; i++) {
  for(int j = 0; j < pn; j++) {
      for(int k = 0; k < pn; k++) {

   double ctrans = dDelta/2. - Delta/2.;
   double xx = x + ctrans + i*dDelta;
   double yy = y + ctrans + j*dDelta;
    double zz = z + ctrans + k*dDelta;

          int ii = round((xx - dDelta/2.)/dDelta);
           int jj = round((yy - dDelta/2.)/dDelta);
                int kk = round((zz - dDelta/2.)/dDelta);

       int place = ii*n*n + jj*n + kk*len;
   double temp = interpolate((struct _interpolate){s, xx, yy, zz});

    equifield[place] += temp;
      }
  }
     }

        } else {

            int ii = round((x - dDelta/2.)/dDelta);
            int jj = round((y - dDelta/2.)/dDelta);
          int kk = round((z - dDelta/2.)/dDelta);

     int place = ii*n*n + jj*n + kk*len;
     equifield[place] += val(s,0,0,0);
 }
    } } end_foreach_level_or_leaf(); }
    return diagi + 1;
}





void equi_output_binary (scalar s, char * fname, int lvl, int diagi) {
    int len = 1;
    int n = 1<<lvl;

    if(pid() == 0){
    FILE * fpm = fopen(fname, "w");

#if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    int sizef = sizeof(double);
    int len = 1;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
  int place = i*n*n + j*n + k*len;
  double temp = (double)(equifield[place]/((double)diagi));
         fseek(fpm, place*sizef, SEEK_SET);
      fwrite(&temp, sizef, 1, fpm); ;
     }
 }
    }
    fclose(fpm);
    }
#if _MPI
    else
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}




void equi_output_ascii (scalar s, char * fname, int lvl, int diagi) {
    int len = 1;
    int n = 1<<lvl;

    if(pid() == 0){
    FILE * fpm = fopen(fname, "w");

#if _MPI
    MPI_Reduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
    int len = 1;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
  int place = i*n*n + j*n + k*len;
  double temp = (equifield[place]/((double)diagi));
         fprintf(fpm, "%g\t", temp);
     }
 }
    }
    fclose(fpm);
    }
#if _MPI
    else
    MPI_Reduce (equifield, NULL, len*n*n*n, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
#endif
}





void equi_import_binary (scalar s, char * fname, int lvl){
    int len = 1;
    int n = 1<<lvl;

    pfree(equifield,__func__,__FILE__,__LINE__);
    equifield = NULL;

    if(!equifield) {
        equifield = pcalloc(len*n*n*n, sizeof(double),__func__,__FILE__,__LINE__);
    }

    double dDelta = 0.9999999*L0/n;


    if(pid() == 0){
    FILE * fpm = fopen(fname, "r");
    int sizef = sizeof(double);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {
  int place = i*n*n + j*n + k*len;
  fseek(fpm, place*sizef, SEEK_SET);
      fread(&equifield[place], sizef, 1, fpm); ;
     }
 }
    }
    fclose(fpm);
    }

#if _MPI
    MPI_Allreduce (MPI_IN_PLACE, equifield, len*n*n*n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif



     { foreach(){

#line 193 "/home/vinlinux/basilisk/src/equi_data.h"
 {
   if(level < lvl){
      double aveval = 0.;
     int pn = 1<<(lvl-level);
     for(int i = 0; i < pn; i++) {
  for(int j = 0; j < pn; j++) {
      for(int k = 0; k < pn; k++) {

   double ctrans = dDelta/2. - Delta/2.;
   double xx = x + ctrans + i*dDelta;
   double yy = y + ctrans + j*dDelta;
    double zz = z + ctrans + k*dDelta;

          int ii = round((xx - dDelta/2.)/dDelta);
           int jj = round((yy - dDelta/2.)/dDelta);
                int kk = round((zz - dDelta/2.)/dDelta);

       int place = ii*n*n + jj*n + kk*len;

    aveval += equifield[place];
      }
  }
     }
     val(s,0,0,0) = 1.*aveval/cube(pn);

        } else {

            int ii = round((x-dDelta/2)/dDelta);
            int jj = round((y-dDelta/2)/dDelta);
          int kk = round((z-dDelta/2)/dDelta);

     int place = ii*n*n + jj*n + kk*len;
     val(s,0,0,0) = 1.* equifield[place];
 }
    } } end_foreach(); }
    boundary(((scalar []){s,{-1}}));
}
#line 8 "./diagnostics.h"



struct sDiag dia;

struct sDiag {
 double Ekin;
 double EkinOld;
 double WdoneOld;
 double rotVol;
 double bE;
 double bEold;
 double diss;
};

struct sEquiDiag {
    int level;
    int ii;
    double dtDiag;
    double startDiag;
    double endDiag;
    double dtOutput;
};

struct sOutput {
    double dtDiag;
    double dtVisual;
    double dtSlices;
    double dtProfile;
    double startAve;
    double dtAve;
    char main_dir[12];
    char dir[30];
    char dir_profiles[60];
    char dir_slices[60];
    char dir_equifields[60];
    char dir_strvel[60];
    char dir_refvel[60];
    char dir_diffbins[60];
    char dir_dts[60];
    int sim_i;
};

struct sbViewSettings {
    double phi;
    double theta;
    double sphi;
    double stheta;
};


struct sOutput out = {.dtDiag = 1., .dtVisual=1., .dtSlices=5., .dtProfile=60., .main_dir="results", .sim_i=0};

struct sEquiDiag ediag = {.level = 5, .ii = 0, .startDiag = 0., .dtDiag = 0., .dtOutput = 0.};

struct sbViewSettings bvsets = {.phi=0., .theta=0., .sphi=0., .stheta=0.};

static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "./diagnostics.h", 65); {
 bvsets.phi = 0.;
 bvsets.theta = -M_PI/6.;
 bvsets.sphi = 0.;
 bvsets.stheta = 0.;
 end_trace("init_0", "./diagnostics.h", 70); } return 0; } 


static int diagnostics_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t+=out.dtDiag);   *ip = i; *tp = t;   return ret; } static int diagnostics (const int i, const double t, Event * _ev) { trace ("diagnostics", "./diagnostics.h", 73); {
 int n = 0.;
 scalar ekin= new_scalar("ekin");
 double tempVol = 0.;
 double tempEkin = 0.;
 double tempDiss = 0.;
 double maxVel = 0.;
 double bEnergy = 0.;


  { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _n = n; double _tempVol = tempVol; double _tempEkin = tempEkin; double _maxVel = maxVel; double _bEnergy = bEnergy; double _tempDiss = tempDiss; 
#line 83

if (!is_constant(cm) && !is_constant(rho)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 83
foreach(){

#line 84 "./diagnostics.h"
 {

  _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
  if(y + Delta/2. <= rot.y0){
       _bEnergy += (cube(Delta)*val_cm(cm,0,0,0))*y*(val(b,0,0,0) - 9.81*(.5)*y/273.);
  }
  {
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.x,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.y,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.z,0,0,0));
  }}
  _maxVel = max(_maxVel, sq(val(ekin,0,0,0)));
  val(ekin,0,0,0) *= 0.5*val_rho(rho,0,0,0)*(cube(Delta)*val_cm(cm,0,0,0));
  _tempEkin += val(ekin,0,0,0);
  _n++;
 } } end_foreach(); }
if (is_constant(cm) && !is_constant(rho)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#line 83
foreach(){

#line 84 "./diagnostics.h"
 {

  _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
  if(y + Delta/2. <= rot.y0){
       _bEnergy += (cube(Delta)*val_cm(cm,0,0,0))*y*(val(b,0,0,0) - 9.81*(.5)*y/273.);
  }
  {
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.x,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.y,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.z,0,0,0));
  }}
  _maxVel = max(_maxVel, sq(val(ekin,0,0,0)));
  val(ekin,0,0,0) *= 0.5*val_rho(rho,0,0,0)*(cube(Delta)*val_cm(cm,0,0,0));
  _tempEkin += val(ekin,0,0,0);
  _n++;
 } } end_foreach(); }
if (!is_constant(cm) && is_constant(rho)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 83
foreach(){

#line 84 "./diagnostics.h"
 {

  _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
  if(y + Delta/2. <= rot.y0){
       _bEnergy += (cube(Delta)*val_cm(cm,0,0,0))*y*(val(b,0,0,0) - 9.81*(.5)*y/273.);
  }
  {
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.x,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.y,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.z,0,0,0));
  }}
  _maxVel = max(_maxVel, sq(val(ekin,0,0,0)));
  val(ekin,0,0,0) *= 0.5*val_rho(rho,0,0,0)*(cube(Delta)*val_cm(cm,0,0,0));
  _tempEkin += val(ekin,0,0,0);
  _n++;
 } } end_foreach(); }
if (is_constant(cm) && is_constant(rho)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#line 83
foreach(){

#line 84 "./diagnostics.h"
 {

  _tempVol += (cube(Delta)*val_cm(cm,0,0,0))*val(fan,0,0,0);
  if(y + Delta/2. <= rot.y0){
       _bEnergy += (cube(Delta)*val_cm(cm,0,0,0))*y*(val(b,0,0,0) - 9.81*(.5)*y/273.);
  }
  {
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.x,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.y,0,0,0));
  }
#line 90
 {
   val(ekin,0,0,0) += sq(val(u.z,0,0,0));
  }}
  _maxVel = max(_maxVel, sq(val(ekin,0,0,0)));
  val(ekin,0,0,0) *= 0.5*val_rho(rho,0,0,0)*(cube(Delta)*val_cm(cm,0,0,0));
  _tempEkin += val(ekin,0,0,0);
  _n++;
 } } end_foreach(); }OMP(omp critical) n += _n;
mpi_all_reduce_double (n, MPI_SUM);
OMP(omp critical) tempVol += _tempVol;
mpi_all_reduce_double (tempVol, MPI_SUM);
OMP(omp critical) tempEkin += _tempEkin;
mpi_all_reduce_double (tempEkin, MPI_SUM);
OMP(omp critical) if (_maxVel > maxVel) maxVel = _maxVel;
mpi_all_reduce_double (maxVel, MPI_MAX);
OMP(omp critical) bEnergy += _bEnergy;
mpi_all_reduce_double (bEnergy, MPI_SUM);
OMP(omp critical) tempDiss += _tempDiss;
mpi_all_reduce_double (tempDiss, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 97
 }

 dia.diss = 1.*tempDiss;
 dia.bE = 1.*bEnergy;
 rot.diaVol = dia.rotVol = 1.*tempVol;
 dia.Ekin = 1.*tempEkin;







 if (pid() == 0){

 char nameOut[90];
 char nameCase[90];
      snprintf(nameOut, 90, "./%s/output", out.dir);
      snprintf(nameCase, 90, "./%s/case", out.dir);
 static FILE * fpout =NULL; if (!fpout || i == 0) fpout = pid() > 0 ? fopen("/dev/null", "w") :  fopen(nameOut, "w");
 static FILE * fpca =NULL; if (!fpca || i == 0) fpca = pid() > 0 ? fopen("/dev/null", "w") :  fopen(nameCase, "w");

 if(t==0.){
  fprintf(fpca,"L0\tinversion\thubU\tTref\tLambda\txr\tyr\tzr\ttheta\tphit\tr\tW\tP\tcu\trampT\tmaxlvl\tminlvl\teps\n");
  fprintf(fpca, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%g\n",
    L0,273./9.81*9.81*(.5)*rot.y0/273. -9.81*(.5)*1.5/273., (-1.),273., 1., rot.x0, rot.y0, rot.z0, rot.theta, rot.phit, rot.R, rot.W, rot.P, rot.cu, rot.rampT, maxlevel, minlevel, eps);

         fprintf(qstderr(),"n\tred\tEkin\tWork\tbE\n");
  fprintf(fpout,"i\tt\tn\tred\tEkin\tWork\tbE\n");
 }
 fprintf(fpout, "%d\t%g\t%d\t%g\t%g\t%g\t%g\n",
  i,t,n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work, dia.bE);

 fprintf(qstderr(), "%d\t%g\t%g\t%g\t%g\n",n,(double)((1<<(maxlevel*3))/n),dia.Ekin,rot.Work,dia.bE);

 fflush(fpout);
 fflush(fpca);

 }

 dia.EkinOld = 1.*dia.Ekin;
 dia.WdoneOld = 1.*rot.Work;
 dia.bEold = 1.*dia.bE;
 delete (((scalar []){ekin,{-1}}));  end_trace("diagnostics", "./diagnostics.h", 140); } return 0; } 



static int profiles_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += out.dtProfile);   *ip = i; *tp = t;   return ret; } static int profiles (const int i, const double t, Event * _ev) { trace ("profiles", "./diagnostics.h", 144);  {
 char nameProf[90];
 snprintf(nameProf, 90, "./%s/t=%05g", out.dir_profiles, t);
 field_profile((struct prof){(scalar *)((scalar []){b,u.x,u.y,u.z,{-1}}), nameProf, .n=100});
 end_trace("profiles", "./diagnostics.h", 148); } return 0; } 


static int equidiags_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = ediag.startDiag);   *ip = i; *tp = t;   return ret; } static int equidiags_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t += ediag.dtDiag);   *ip = i; *tp = t;   return ret; } static int equidiags (const int i, const double t, Event * _ev) { trace ("equidiags", "./diagnostics.h", 151);  {
 ediag.ii = equi_diag(b, ediag.level, ediag.ii);
 end_trace("equidiags", "./diagnostics.h", 153); } return 0; } 


static int equioutputs_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = (ediag.dtOutput + ediag.startDiag));   *ip = i; *tp = t;   return ret; } static int equioutputs_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t += ediag.dtOutput);   *ip = i; *tp = t;   return ret; } static int equioutputs (const int i, const double t, Event * _ev) { trace ("equioutputs", "./diagnostics.h", 156);  {
    char nameEquif[90];
    snprintf(nameEquif, 90, "%s/%st=%05g", out.dir_equifields, "equifield", t);
    equi_output_binary(b, nameEquif, ediag.level, ediag.ii);

    ediag.ii = 0;

    pfree(equifield,__func__,__FILE__,__LINE__);
    equifield = NULL;
 end_trace("equioutputs", "./diagnostics.h", 165); } return 0; } 
#line 217 "./diagnostics.h"
static int fanvelocity_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t+=1);   *ip = i; *tp = t;   return ret; } static int fanvelocity (const int i, const double t, Event * _ev) { trace ("fanvelocity", "./diagnostics.h", 217);  {
    if(pid() == 0) {
    char nameStrvel[90];
    snprintf(nameStrvel, 90, "%st=%05g", out.dir_strvel, t);
    FILE * fpstr = fopen(nameStrvel, "w");
    fprintf(fpstr, "R,v,vx,vy,vz\n");

    double length = 200;
    int ntot = 200;
    double vels1[ntot];

    double xf0 = rot.x0;
    double yf0 = rot.y0;
    double zf0 = rot.z0;

    for(int n; n < ntot; n++) {
 double dist = length*n/ntot - 100;
 double phi = dist < 0 ? M_PI/360. : 0.;
 double theta = 97*M_PI/180;

 double nfx1 = sin(theta)*cos(phi);
 double nfz1 = sin(theta)*sin(phi);
 double nfy1 = cos(theta);

 double xx1 = xf0 + dist*nfx1;
 double yy1 = yf0 + dist*nfz1;
 double zz1 = zf0 + dist*nfy1;


 double valx1 = interpolate((struct _interpolate){u.x, xx1, yy1, zz1});
 double valy1 = interpolate((struct _interpolate){u.y, xx1, yy1, zz1});
 double valz1 = interpolate((struct _interpolate){u.z, xx1, yy1, zz1});

        vels1[n] = sqrt(sq(valx1) + sq(valy1) + sq(valz1));
 fprintf(fpstr, "%g,%g,%g,%g,%g\n", dist, vels1[n], valx1, valy1, valz1);
    }
    fclose(fpstr);
    }
 end_trace("fanvelocity", "./diagnostics.h", 255); } return 0; } 

static int refvelocity_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t+=1);   *ip = i; *tp = t;   return ret; } static int refvelocity (const int i, const double t, Event * _ev) { trace ("refvelocity", "./diagnostics.h", 257);  {
    if(pid() == 0) {
    char nameRefvel[90];
    snprintf(nameRefvel, 90, "%st=%05g", out.dir_refvel, t);
    FILE * fpstr = fopen(nameRefvel, "w");
    fprintf(fpstr, "x,v,vx,vy,vz\n");

    double length = L0;
    int ntot = L0/2;
    double vels1[ntot];

    double xf0 = 0;
    double yf0 = rot.y0-7;
    double zf0 = L0/2;

    for(int n; n < ntot; n++) {
 double dist = length*n/ntot;
 double xx = xf0 + dist;
 double yy = yf0;
 double zz = zf0;

 double valx1 = interpolate((struct _interpolate){u.x, xx, yy, zz});
 double valy1 = interpolate((struct _interpolate){u.y, xx, yy, zz});
 double valz1 = interpolate((struct _interpolate){u.z, xx, yy, zz});

        vels1[n] = sqrt(sq(valx1) + sq(valy1) + sq(valz1));
 fprintf(fpstr, "%g,%g,%g,%g,%g\n", xx, vels1[n], valx1, valy1, valz1);
    }
    fclose(fpstr);
    }
 end_trace("refvelocity", "./diagnostics.h", 287); } return 0; } 

static int dts_meas_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += 1);   *ip = i; *tp = t;   return ret; } static int dts_meas (const int i, const double t, Event * _ev) { trace ("dts_meas", "./diagnostics.h", 289);  {
    if(pid() == 0) {
 char nameDtshor[90];
     snprintf(nameDtshor, 90, "%shor_t=%05g", out.dir_dts, t);
        FILE * fpstrhor = fopen(nameDtshor, "w");
        fprintf(fpstrhor, "x,y,z,b\n");

 double lengthhor = L0;
        int ntothor = L0;

 double xf0 = rot.x0;
 double yf0 = rot.y0;
 double zf0 = rot.z0;

 for(int n = 0; n <= ntothor; n++) {
     double dist = lengthhor*n/ntothor;

     double xx = xf0;
     double yy = rot.y0-9.;
     double zz = dist;

     double valb = interpolate((struct _interpolate){b, xx, yy, zz});

     fprintf(fpstrhor, "%g,%g,%g,%g\n", xx, yy, zz, valb);
 }
  fclose(fpstrhor);

 char nameDtsver[90];
     snprintf(nameDtsver, 90, "%sver_t=%05g", out.dir_dts, t);
        FILE * fpstrver = fopen(nameDtsver, "w");
        fprintf(fpstrver, "x,y,z,b\n");

 double lengthver = 20;
        int ntotver = 200;

 for(int n = 0; n <= ntotver; n++) {
     double dist = lengthver*n/ntotver;

     double xx = xf0 + 30;
     double yy = dist;
     double zz = zf0;

     double valb = interpolate((struct _interpolate){b, xx, yy, zz});

     fprintf(fpstrver, "%g,%g,%g,%g\n", xx, yy, zz, valb);
 }
  fclose(fpstrver);
    }

 end_trace("dts_meas", "./diagnostics.h", 338); } return 0; } 
#line 364 "./diagnostics.h"
static int movies_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += out.dtVisual);   *ip = i; *tp = t;   return ret; } static int movies (const int i, const double t, Event * _ev) { trace ("movies", "./diagnostics.h", 364);  {
    scalar l2= new_scalar("l2");

    lambda2(u,l2);
    boundary(((scalar []){l2,{-1}}));

    view((struct _view_set){.fov=25, .tx = 0., .ty = 0., .phi=bvsets.phi, .theta=bvsets.theta, .width = 800, .height = 800});

    { begin_translate((struct _translate){-rot.x0,-rot.y0,-rot.z0}); {
        box((struct _box){.notics=false});
        isosurface((struct _isosurface){"l2", .v=-0.02, .color="b", .min=9.81*(.5)*0/273., .max=9.81*(.5)*L0/273.});
 draw_vof((struct _draw_vof){"fan", .fc = {1,0,0}});
    } end_translate(); }
    { begin_translate((struct _translate){-rot.z0,-rot.y0, -L0});{
       squares((struct _squares){"u.x", .n = {0,0,1}, .alpha=rot.z0, .min=-fabs((-1.)), .max=fabs((-1.))});
        cells((struct _cells){.n = {0,0,1}, .alpha = rot.z0});
    } end_translate(); }

    { begin_translate((struct _translate){0.,-rot.y0,-rot.z0});{
        squares((struct _squares){"u.x", .n = {1,0,0}, .alpha=rot.x0, .min=-fabs((-1.)), .max=fabs((-1.))});
    } end_translate(); }


    char nameVid1[90];
    snprintf(nameVid1, 90, "ppm2mp4 -r %g ./%s/visual_3d_vel.mp4", 10., out.dir);
    save((struct _save){nameVid1});
    clear();

    view((struct _view_set){.fov=25, .tx = 0., .ty = 0., .phi=bvsets.phi, .theta=bvsets.theta, .width = 800, .height = 800});

    { begin_translate((struct _translate){-rot.x0,-rot.y0,-rot.z0}); {
        box((struct _box){.notics=false});
        isosurface((struct _isosurface){"l2", .v=-0.02, .color="b", .min=9.81*(.5)*0./273., .max=9.81*(.5)*L0/273.});
 draw_vof((struct _draw_vof){"fan", .fc = {1,0,0}});
    } end_translate(); }
    { begin_translate((struct _translate){-rot.z0,-rot.y0, -L0});{
       squares((struct _squares){"b", .n = {0,0,1}, .alpha=rot.z0, .min=9.81*(.5)*0./273., .max=9.81*(.5)*L0/273.});
        cells((struct _cells){.n = {0,0,1}, .alpha = rot.z0});
    } end_translate(); }

    { begin_translate((struct _translate){0.,-rot.y0,-rot.z0});{
        squares((struct _squares){"b", .n = {1,0,0}, .alpha=rot.x0, .min=9.81*(.5)*0./273., .max=9.81*(.5)*L0/273.});
    } end_translate(); }


    char nameVid2[90];
    snprintf(nameVid2, 90, "ppm2mp4 -r %g ./%s/visual_3d_b.mp4", 10., out.dir);
    save((struct _save){nameVid2});
    clear();
 delete (((scalar []){l2,{-1}}));  end_trace("movies", "./diagnostics.h", 413); } return 0; } 


static int slices_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t=out.dtSlices);   *ip = i; *tp = t;   return ret; } static int slices_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t+=out.dtSlices);   *ip = i; *tp = t;   return ret; } static int slices (const int i, const double t, Event * _ev) { trace ("slices", "./diagnostics.h", 416);  {
    char nameSlice[90];
    coord slice = {1., 0., 1.};
    int res = L0/2;

    for(double yTemp = 6; yTemp<=13; yTemp+=2.) {
 slice.y = (rot.y0 - yTemp)/L0;

     snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
     FILE * fpsli = fopen(nameSlice, "w");
     output_slice((struct sOutputSlice){.list = (scalar *)((scalar []){b,{-1}}), .fp = fpsli, .n = res, .linear = true, .plane=slice});
     fclose(fpsli);
    }

    for(double yTemp = 0; yTemp<=1.; yTemp+=2.) {
 slice.y = (rot.y0 - yTemp)/L0;

     snprintf(nameSlice, 90, "%st=%05gy=%03g", out.dir_slices, t, yTemp);
     FILE * fpsli = fopen(nameSlice, "w");
     output_slice((struct sOutputSlice){.list = (scalar *)((scalar []){b,{-1}}), .fp = fpsli, .n = res, .linear = true, .plane=slice});
     fclose(fpsli);
    }
#line 461 "./diagnostics.h"
 end_trace("slices", "./diagnostics.h", 461); } return 0; } 





void sim_dir_create(){

    sprintf(out.dir, "./%s/%s%02d", out.main_dir, sim_ID, out.sim_i);
    sprintf(out.dir_profiles, "%s/profiles/", out.dir);
    sprintf(out.dir_slices, "%s/slices/", out.dir);
    sprintf(out.dir_equifields, "%s/equifields/", out.dir);
    sprintf(out.dir_strvel, "%s/strvel/", out.dir);
    sprintf(out.dir_refvel, "%s/refvel/", out.dir);
    sprintf(out.dir_diffbins, "%s/diffbins/", out.dir);
    sprintf(out.dir_dts, "%s/dts/", out.dir);


    if (pid() == 0){
    struct stat st = {0};
    if (stat(out.main_dir, &st) == -1) {
        mkdir(out.main_dir, 0777);
    }

    if (stat(out.dir, &st) == -1) {
        mkdir(out.dir, 0777);
    }
    if (stat(out.dir_slices, &st) == -1) {
        mkdir(out.dir_slices, 0777);
    }
    if (stat(out.dir_profiles, &st) == -1) {
        mkdir(out.dir_profiles, 0777);
    }
    if (stat(out.dir_equifields, &st) == -1) {
        mkdir(out.dir_equifields, 0777);
    }