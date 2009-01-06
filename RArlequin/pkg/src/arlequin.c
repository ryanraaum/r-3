#include <stdio.h>
#include <R.h>
#include <Rdefines.h>

void 
convolve(double *a, int *na, double *b, int *nb, double *ab)
{
  int i, j, nab = *na + *nb - 1;
     
  for(i = 0; i < nab; i++)
    ab[i] = 0.0;
  for(i = 0; i < *na; i++)
    for(j = 0; j < *nb; j++)
      ab[i + j] += a[i] * b[j];
}

void read_arlequin_file(const char *filename, 
                        int *ploidy,
                        char **genotypes, 
                        int *sample_sizes, 
                        char *groups,
                        char **markers ) 
{
  printf("will try to read '%s'\n", filename);
}

void test(const char *filename, int *a)
{
  a[0] = 2;
}

SEXP test2(SEXP a)
{
  int i, na;
  double *xa, *xres;
  SEXP res;

  PROTECT(a = AS_NUMERIC(a));
  na = LENGTH(a);
  PROTECT(res = NEW_NUMERIC(na));
  xa = NUMERIC_POINTER(a);
  xres = NUMERIC_POINTER(res);
  for (i=0; i<na; i++)
    xres[i] = xa[i];
  UNPROTECT(2);
  return(res);
}
