#include <stdio.h>

#include "arlequin.h"

void usage()
{
  printf("Usage: test <filename>\n");
}

int main( int argc, const char* argv[] )
{
  int ploidy;
  char *genotypes; 
  int sample_sizes; 
  char groups;
  char *markers;  
  
  if (argc != 2) {
    usage();
  }
  else {
    read_arlequin_file(argv[1], &ploidy, &genotypes, &sample_sizes, &groups, &markers); 
	}
}

