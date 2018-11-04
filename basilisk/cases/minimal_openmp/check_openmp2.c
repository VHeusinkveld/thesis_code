int main(){
  init_grid(1<<5);
  double a = 0;
#if _OPENMP
  printf("Apparently _OPENMP != 0\n");
  foreach_leaf()
#else
  foreach(reduction(+:a))
#endif
    a += dv();
  printf("total: %g\n", a);
}
