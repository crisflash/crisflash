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

/* 
   Written by Adrien Jacquin, June 2017.
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "read.h"
#include "vcf.h"

int string_ends_with(const char * str, const char * suffix)
{
  int str_len = strlen(str);
  int suffix_len = strlen(suffix);

  return 
    (str_len >= suffix_len) &&
    (0 == strcmp(str + (str_len-suffix_len), suffix));
}

char *add_suffix(char *str, char *str2) {
  char *new_str = malloc(sizeof(char)*(strlen(str)+strlen(str)+1));
  strcpy(new_str, str);
  strcat(new_str, str2);
  return new_str;
}

void file_accessible(char *fname, int mode) {
  // Following modes are supported:
  // F_OK tests existence also (R_OK,W_OK,X_OK). for readable, writeable, executable
  if ( access (fname, mode) != 0 )
    {
      printf("[crisflash] ERROR: %s not accessible!\n", fname);
      exit(1);
    }
}

static int usage()
{
  /*	while ((c = getopt(argc, argv, "g:b:o:s:v:m:@:it:S:hu")) != -1)  */
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: crisflash (A tool for CRISPR/Cas9 sgRNA design and off-target validation)\n");
  fprintf(stderr, "Version: %s\n\n", CRISFLASH_VERSION);
  fprintf(stderr, "Usage:   crisflash -g <genome.fa> -s <input.fa> -o <output> [options]\n\n");
  fprintf(stderr, "Options:\n");
  /* input file option for reference genome */
  fprintf(stderr, "          -g FILE\tFASTA format reference genome.\n");
  /* input file options for candidate area and gRNAs */
  fprintf(stderr, "          -s FILE\tFASTA file containing candidate sequence.\n"); /* input sequence */
  fprintf(stderr, "          -p PAM\tPAM sequence. Default: NGG\n");
  fprintf(stderr, "          -V FILE\tphased VCF file.\n"); /* vcf file */  
  /* boolean options for output file format */
  /*
    -B bam output
    -C cas-offinder output
  */
  fprintf(stderr, "          -B\t\t save output in BED format, with sequence provided on comment field and off-target score on score field. (Default)\n");
  fprintf(stderr, "          -C\t\t save output in cas-offinder format.\n");
  /* output file name */
  fprintf(stderr, "          -o FILE\toutput file name saved in BED format.\n");
  /* options for parameters */
  fprintf(stderr, "          -m INT\tNumber of mismatches allowed. Default: 2.\n");
  fprintf(stderr, "          -t INT\tNumber of threads for off-target scoring. Default: 1.\n");
  fprintf(stderr, "          -u\tExclude low complexity genomic sequences marked lowercase in soft masked genomes.\n");
  /* help and version options */
  fprintf(stderr, "          -h\tPrint help.\n");
  fprintf(stderr, "          -v\tPrint version.\n\n");

  fprintf(stderr, "Examples: crisflash -g genome.fa -s candidate_area.fa -o validated_gRNAs.bed -m 5\n");
  fprintf(stderr, "          crisflash -g genome.fa -s candidate_area.fa -o validated_gRNAs.bed -m 5 -C\n");
  fprintf(stderr, "          crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o validated_gRNAs.bed\n");
  fprintf(stderr, "          crisflash -g genome.fa -o all_unscored_gRNAs.bed\n");
  fprintf(stderr, "          crisflash -s candidate_area.fa -o all_unscored_gRNAs.bed\n");
  fprintf(stderr, "\n\n");
  
  exit(1);
  return 1;
}

int main(int argc, char **argv)
{
	// Main function which will print results according to user inputs
	char* referenceGenomePath = NULL; // path for the reference genome
	char* sequencePath = NULL; // path for the sequence
	char* vcfPath = NULL; // path for the VCF file
	char* prefixName = NULL; // prefix for the output BED files, default = results
	char* outFile = NULL; // output file in BED format for results.
	int maxMismatch = 2; // number of mismatches allowed for matching, default = 2
	int threadFlag = 0; // booelan: use a multi-threaded fucntion for matching itself, default = False
	int nr_of_threads = 1; // number of threads in case of using multi-threaded function, default = 1
	int upper_case_only = 0; // 1 if soft-masked lowercase bases are excluded both from genome indexing to trie as well as from candidate identification
	trie *T;
	int index;
	int c;
	char *bedCandidatesPath = NULL;
	int matchItselfPathBool = 0;
	int inputMatchItselfBool = 0;
	int outFileType = 1; // default is 1 = BED format; 2 = cas-offinder format;
	int printGRNAsOnly = 0;

	// default PAM sequence is 'NGG'
	char * pam =  malloc(sizeof(char)*4);			     
	strcpy(pam,"NGG");
	       
	// while ((c = getopt(argc, argv, "g:b:o:s:v:m:@:it:lS:h")) != -1)
	while ((c = getopt(argc, argv, "BCg:o:s:vV:m:t:lhp:")) != -1)
	  {
	    switch (c)
	      {		
	      case 'g':
		file_accessible(optarg, R_OK);
		if (string_ends_with(optarg,".fa")) { referenceGenomePath = optarg; }
		else {
		  fprintf(stderr,"[crisflash] ERROR: file does not appear to be fasta file (.fa).\n");
		  exit(1);
		}
		break;
	      case 'o':
		outFile = optarg;
		break;
	      case 's':
		file_accessible(optarg, R_OK);
		if (string_ends_with(optarg,".fa")) { sequencePath = optarg; }
		else {
		  fprintf(stderr,"[crisflash] ERROR: file does not appear to be fasta file (.fa).\n");
		  exit(1);
		}
		break;
	      case 'V':
		if (string_ends_with(optarg,".vcf") || string_ends_with(optarg,".VCF")) { vcfPath = optarg; }
		else {
		  fprintf(stderr,"[crisflash] ERROR: file does not appear to be vcf file (.vcf).\n");
		  exit(1);
		}
		vcfPath = optarg;
		break;
	      case 'm':
		maxMismatch = atoi(optarg);
		break;
	      case 't':
		threadFlag = 1;
		nr_of_threads = atoi(optarg);
		break;
	      case 'u':
		upper_case_only = 1; // When set to 1 lower case / low complexity sequences in soft masked genome are ignored.
		break;
	      case 'h':
		usage();
		break;
	      case 'v':
		fprintf(stdout,"crisflash %s\n",CRISFLASH_VERSION);
		exit(0);
		break;
	      case 'B':
		/* Set output file type BED */
		outFileType = 1;
		break;
	      case 'C':
		/* Set output file type cas-offinder */
		outFileType = 2;
		break;
	      case 'p':
		free(pam);
		pam = malloc((sizeof(char)*strlen(optarg))+1);
		strcpy(pam, optarg);
		break;      	
	      case '?':
		if (isprint(optopt))
		  {
		    fprintf(stderr, "[crisflash] ERROR: Unknown option `-%c'.\n", optopt);
		  }
		else
		  {
		    fprintf(stderr, "[crisflash] ERROR: Unknown option character `\\x%x'.\n", optopt);
		  }
		return 1;
	      default:
		fprintf(stderr, "[crisflash] ERROR: Unknown option '-%c'\n", c);
		exit(1);
	      }
	  }
	
	// We need values for at least two arguments.
	int acount = 0; 
	while(argv[++acount] != NULL);
	if (acount == 2) {
	  usage();
	}
	for (index = optind; index < argc; index++)
	  {
	    fprintf(stderr,"[crisflash] ERROR: Non-option argument '%s'. Exiting!\n", argv[index]);
	    exit(1);
	  }

	// if no reference genome or input sequence, display help and exit
	if ((!referenceGenomePath)&&(!sequencePath)) {
	  usage();
	  exit(1);
	}

	// if output file had not been set, display help and exit.
	if (!outFile) {
	  usage();
	  exit(1);
	}

	/*
	  If values for only referenceGenome or sequence path, print all gRNA candidates and exit.
	  No need to build Trie and match candidates.
	*/
	if((referenceGenomePath)&&(!sequencePath)) {
	  readFaToTrie(T, referenceGenomePath, pam, outFile, 0, upper_case_only,1);
	  exit(0);
	}
	if((!referenceGenomePath)&&(sequencePath)) {
	  readFaToTrie(T, sequencePath, pam, outFile, 0, upper_case_only,1);
	  exit(0);
	}

	// We have values for referene genome and sequence for target area. Therefore, we will index reference genome to Trie!
	T = TrieCreate(PROTOSPACER_LENGTH + strlen(pam));
	
	if (!vcfPath) {
	  readFaToTrie(T, referenceGenomePath, pam, outFile, 0, upper_case_only,0);
	}
	else {
	  // variant data has been provided. We create haplotype based genomic sequences first and load them to Trie one by one.
	  char *output_VCF1 = malloc(sizeof(char)*(strlen(referenceGenomePath) + strlen(vcfPath) + 6));
	  char *output_VCF2 = malloc(sizeof(char)*(strlen(referenceGenomePath) + strlen(vcfPath) + 6));	    
	  int phased;
	  int append = 0;
	  strcpy(output_VCF1,referenceGenomePath);
	  strcpy(output_VCF2,referenceGenomePath);
	  strcat(output_VCF1,"_");
	  strcat(output_VCF2,"_");
	  strcat(output_VCF1,vcfPath);
	  strcat(output_VCF2,vcfPath);
	  strcat(output_VCF1,"1.fa");
	  strcat(output_VCF2,"2.fa");
	  fprintf(stdout,"[crisflash] Creating phased genomes for %s.\n", vcfPath);
	  phased = VCF_to_genome(referenceGenomePath, vcfPath, output_VCF1, output_VCF2);
	  // printf("============== 1st HAPLOTYPE ==============\n");
	  readFaToTrie(T, output_VCF1, pam, outFile, append, upper_case_only,0);
	  if (phased == 1)
	    {
	      // printf("============== 2nd HAPLOTYPE ==============\n");
	      append += T->nr_sequences;
	      readFaToTrie(T, output_VCF2, pam, outFile, append, upper_case_only,0);
	    }
	  free(output_VCF1);
	  free(output_VCF2);
	}	  
	
	TrieAMatchSequenceThreads(T, sequencePath, maxMismatch, outFile, outFileType, pam, upper_case_only, nr_of_threads,0);

	// We will not free Trie as the process is slow and the block of memory will be handed back to OS on exit anyway.
	/*
	  T = TrieDestroy(T);
	  free(pam);
	  if(bedCandidatesPath) { free(bedCandidatesPath); }
	  if(outFile) { free(outFile); }
	*/
	return 0;
}
