int main(){
  int i = 0, j = 0;
#if _OPENMP
  i++;
#endif
@if _OPENMP
  j++;
@endif
  printf("i = %d, j = %d\n", i, j);
}
