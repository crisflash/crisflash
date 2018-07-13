/*******                                 

Copyright 2018 Adrien Jacquin and Margus Lukk, CRUK-CI, University of Cambridge

This file is part of the crisflash software tool.

Crisflash is free software: you can redistribute it
and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

Crisflash is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied   
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Crisflash.  If not, see <http://www.gnu.org/licenses/>.
 
********/
#ifdef __gnu_linux__
        #define _XOPEN_SOURCE >= 500
        #define _GNU_SOURCE
#endif

/* Written by Adrien Jacquin, July 2017 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "vcf.h"

static int vcf_usage()
{  
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: vcf\n");
  fprintf(stderr, "Description: Returns one or two genome sequences constructed by overlaying\n");
  fprintf(stderr, "             the reference sequence with variant data provided in VCF file.\n");
  fprintf(stderr, "             The number of output files depends on whether the VCF file is phased:\n");
  fprintf(stderr, "             two files (one per haplotype) for phased VCF and one file otherwise.\n");
  fprintf(stderr, "Version: %s\n", CRISFLASH_VERSION);
  fprintf(stderr, "Usage:   vcf -g <genome.fa> -v <variants.vcf> [options] \n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "          -g FILE\tFASTA format reference genome.\n");
  fprintf(stderr, "          -v FILE\tVCF file.\n"); /* vcf file */
  fprintf(stderr, "          -o STRING\tPrefix for the output file name(s). Default prefix: 'result'.\n");
  fprintf(stderr, "          -h\t\tPrint help.\n");
  fprintf(stderr, "\n\n");

  exit(1);
  return 1;
}

int main(int argc, char **argv)
{
	char* referenceGenomePath = NULL;
	char* vcfPath = NULL;
	char* prefixName = NULL;
	int index;
	int c;
	while ((c = getopt(argc, argv, "g:v:o:h")) != -1)
	{
		switch (c)
		{
		case 'g':
			referenceGenomePath = optarg;
			break;
		case 'v':
			vcfPath = optarg;
			break;
		case 'o':
			prefixName = optarg;
			break;
		case 'h':
		        vcf_usage();
			break;
		case '?':
			if (isprint(optopt))
			{
				fprintf(stderr, "WARNING: Unknown option `-%c'.\n", optopt);
			}
			else
			{
				fprintf(stderr, "WARNING: Unknown option character `\\x%x'.\n", optopt);
			}
			return 1;
		default:
		  abort();
		}
	}

	for (index = optind; index < argc; index++)
	{
		printf ("\n\nERROR: Non-option argument %s\n", argv[index]);
	}
	if (!vcfPath || !referenceGenomePath)
	{
	  fprintf(stderr,"\n\nERROR: either reference genome ( -g genome.fa ) or VCF file (-v var.vcf) not provided!\n\n\n");
	  vcf_usage();
	  exit(1);
	}
	
	if (!prefixName) {
	  char *gbasename = basename(referenceGenomePath);
	  char *vbasename = basename(vcfPath);
	  prefixName = malloc(sizeof(char)*(strlen(gbasename) + strlen(vbasename) + 6));
	  strcpy(prefixName,gbasename);
	  strcat(prefixName,"_");
	  strcat(prefixName,vbasename);
	}
	fprintf(stderr,"Using output prefix '%s'\n", prefixName);

	int phased;
	char *genome1 = malloc(sizeof(char)*strlen(prefixName) + 5);
	char *genome2 = malloc(sizeof(char)*strlen(prefixName) + 5);
	strcpy(genome1,prefixName);
	strcpy(genome2,prefixName);
	strcat(genome1,"1.fa");
	strcat(genome2,"2.fa");
	phased = VCF_to_genome(referenceGenomePath, vcfPath, genome1, genome2);
	free(genome1);
	free(genome2);
	return 0;
}
