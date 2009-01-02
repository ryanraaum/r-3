#ifndef  arlequin_INC
#define  arlequin_INC

#include <R.h>

int read_arlequin_file( const char *filename, 
                        int *ploidy,
                        char **genotypes, 
                        int *sample_sizes, 
                        char *groups,
                        char **markers );

#endif   /* ----- #ifndef arlequin_INC  ----- */
