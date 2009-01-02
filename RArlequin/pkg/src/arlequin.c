#include <stdio.h>
#include "arlequin.h"


void convolve(double *a, int *na, double *b, int *nb, double *ab)
     {
       int i, j, nab = *na + *nb - 1;
     
       for(i = 0; i < nab; i++)
         ab[i] = 0.0;
       for(i = 0; i < *na; i++)
         for(j = 0; j < *nb; j++)
           ab[i + j] += a[i] * b[j];
     }

int read_arlequin_file( const char *filename, 
                        int *ploidy,
                        char **genotypes, 
                        int *sample_sizes, 
                        char *groups,
                        char **markers ) 
{
  printf("will try to read '%s'\n", filename);
  return 0;
}
