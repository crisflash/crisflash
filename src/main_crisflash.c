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
  /*	while ((c = getopt(argc, argv, "g:b:o:s:v:m:@:it:lS:h")) != -1)  */
  fprintf(stderr, "\n");
  fprintf(stderr, "Program: crisflash (A tool for CRISPR/Cas9 sgRNA design and off-target validation)\n");
  fprintf(stderr, "Version: %s\n\n", CRISFLASH_VERSION);
  fprintf(stderr, "Usage:   crisflash (-g <genome.fa> | -b <genome_sgRNAs.bed>) [options] \n\n");
  fprintf(stderr, "Options:\n");
  /* input file options for reference genome */
  fprintf(stderr, "          -g FILE\tFASTA format reference genome.\n");
  fprintf(stderr, "          -b FILE\tBED file containing all gRNAs of a genome (the file can be created with -a flag).\n");
  /* input file options for candidate area and gRNAs */
  fprintf(stderr, "          -s FILE\tFASTA file containing candidate sequence.\n"); /* input sequence */
  fprintf(stderr, "          -@ FILE\tBED file containing candidate gRNA coordinates and sequences on comment field.\n");
  fprintf(stderr, "          -V FILE\tphased VCF file.\n"); /* vcf file */  
  /* output file options */
  fprintf(stderr, "          -o FILE\toutput saved in BED format. 'sequence/exact match count/off-target count' in comment field, off-target score in score field.\n");
  fprintf(stderr, "          -a FILE\toutput of all gRNAs in a genome in BED format. No validation of off-targets.\n"); /* vcf file */
  fprintf(stderr, "          -A FILE\toutput of all gRNAs in a genome in BED format with off-target validation. (NB! takes long time even with -t option)\n");
  /* parameter options */
  fprintf(stderr, "          -m INT\tNumber of mismatches allowed. Default: 2.\n");
  fprintf(stderr, "          -t INT\tNumber of threads to use for off-target validation. Default: 1.\n");
  fprintf(stderr, "          -l\tInclude low complexity areas to off-target validation. Default: use only upper case sequences in soft masked genome.\n");
  /* help and version options */
  fprintf(stderr, "          -h\tPrint help.\n");
  fprintf(stderr, "          -v\tPrint version.\n\n");
  fprintf(stderr, "Examples: crisflash -g genome.fa -s candidate_area.fa -o validated_gRNAs.bed -m 5\n");
  fprintf(stderr, "          crisflash -g genome.fa -@ candidate_gRNAs.bed -o validated_gRNAs.bed\n");
  fprintf(stderr, "          crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o validated_gRNAs.bed\n");
  fprintf(stderr, "          crisflash -g genome.fa -a genome_gRNAs_dump.bed\n");
  fprintf(stderr, "          crisflash -b genome_gRNAs_dump.bed -t 40 -A all_gRNAs_validated.bed\n");
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
	char* bedFilePath = NULL; // path for the input BED file
	char* prefixName = NULL; // prefix for the output BED files, default = results
	char* outFile = NULL; // output file in BED format for results.
	int maxMismatch = 2; // number of mismatches allowed for matching, default = 2
	int threadFlag = 0; // booelan: use a multi-threaded fucntion for matching itself, default = False
	int nr_of_threads = 1; // number of threads in case of using multi-threaded function, default = 1
	int upper_case_only = 1; // boolean: find gRNAs in upper case only sequences, default = True (faster)
	trie *T = TrieCreate(LENGTH_SEQ);
	int index;
	int c;
	char *bedCandidatesPath = NULL;
	int matchItselfPathBool = 0;
	int inputMatchItselfBool = 0;
	
	char *allsgRNAsFile = NULL;

	// while ((c = getopt(argc, argv, "g:b:o:s:v:m:@:it:lS:h")) != -1)
	while ((c = getopt(argc, argv, "a:A:g:b:o:s:vV:m:@:t:l:h")) != -1)
	  {
	    switch (c)
	      {
	      case 'a':
		if (string_ends_with(optarg,".bed")) { allsgRNAsFile = optarg; }
		else { allsgRNAsFile = add_suffix(optarg,".bed"); }
		break;
	      case 'A':		
		if (string_ends_with(optarg,".bed")) { outFile = optarg; }
		else { outFile = add_suffix(optarg,".bed"); }
		matchItselfPathBool = 1;
		break;
	      case 'g':
		file_accessible(optarg, R_OK);
		if (string_ends_with(optarg,".fa")) { referenceGenomePath = optarg; }
		else {
		  fprintf(stderr,"[crisflash] ERROR: file does not appear to be fasta file (.fa).\n");
		  exit(1);
		}
		break;
	      case 'b':
		if (string_ends_with(optarg,".bed")) { bedFilePath = optarg; }
		else {
		  fprintf(stderr,"[crisflash] ERROR: file does not appear to be bed file (.bed).\n");
		  exit(1);
		}
		break;
	      case 'o':
		if (string_ends_with(optarg,".bed")) { outFile = optarg; }
		else { outFile = add_suffix(optarg,".bed"); }
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
	      case '@':
		file_accessible(optarg, R_OK);
		if (string_ends_with(optarg,".bed")) { bedCandidatesPath = optarg; }
		else { bedCandidatesPath = add_suffix(optarg,".bed"); }
		inputMatchItselfBool = 1;
		break;      	
	      case 'm':
		maxMismatch = atoi(optarg);
		break;
	      case 't':
		threadFlag = 1;
		nr_of_threads = atoi(optarg);
		break;
	      case 'l':
		upper_case_only = 0; // This means we will search in lower case sequences as well
		break;
	      case 'h':
		usage();
		break;
	      case 'v':
		fprintf(stdout,"crisflash %s\n",CRISFLASH_VERSION);
		exit(0);
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
		exit(1);
	      }
	  }
	
	int acount = 0; 
	while(argv[++acount] != NULL);
	if (acount == 1) {
	  usage();	  
	}
	for (index = optind; index < argc; index++)
	  {
	    fprintf(stderr,"[crisflash] ERROR: Non-option argument '%s'. Exiting!\n", argv[index]);
	    exit(1);
	  }

	// if both reference genome and BED containing all gRNAs missing
	if (!referenceGenomePath && !bedFilePath) {	  
	  fprintf(stderr,"[crisflash] ERROR: No reference genome (option -g) or a BED with all genomic gRNAs (option -b). Exiting!\n\n");
	  exit(1);
	}

	// if vcf, reference genome has to be provided as well.
	if (vcfPath && !referenceGenomePath)
	{
	  fprintf(stderr,"[crisflash] ERROR: reference genome is missing! Use option -g. Exiting!\n");
	  exit(1);
	}

	if ((sequencePath || bedCandidatesPath) && (!outFile)) {
	  fprintf(stderr,"[crisflash] ERROR: output file is missing! Use option -o. Exiting!\n");
	  exit(1);
	}

	if (matchItselfPathBool && inputMatchItselfBool) {
	  fprintf(stderr,"[crisflash] ERROR: use either -A or -@ option, not both at the same time. Exiting!\n");
	  exit(1);
	}

	// If -A
	if (matchItselfPathBool) {      
	  if(!bedFilePath) { 
            char *bname = basename(referenceGenomePath);
            if (!allsgRNAsFile) {
              allsgRNAsFile = malloc(sizeof(char)*(strlen(bname) + 2));
              int i = 0;
              for(;i < strlen(bname); i++) { allsgRNAsFile[i] = bname[i]; }
              allsgRNAsFile[i-2] = 'b';
              allsgRNAsFile[i-1] = 'e';
              allsgRNAsFile[i] = 'd';
              allsgRNAsFile[i+1] = '\0';
            }
	    bedCandidatesPath = allsgRNAsFile;
	  }
	  else { 
	    bedCandidatesPath = bedFilePath;
	  }	  
	}
	
	// If reference genome is the only input output all gRNAs in the genome.
	if (!sequencePath && !bedCandidatesPath && !allsgRNAsFile) {
	  if (referenceGenomePath) {
	    char *bname = basename(referenceGenomePath);
	    if (!allsgRNAsFile) {
	      allsgRNAsFile = malloc(sizeof(char)*(strlen(bname) + 2));
	      int i = 0;
	      for(;i < strlen(bname); i++) { allsgRNAsFile[i] = bname[i]; }
	      allsgRNAsFile[i-2] = 'b';
	      allsgRNAsFile[i-1] = 'e';
	      allsgRNAsFile[i] = 'd';
	      allsgRNAsFile[i+1] = '\0';
	    }
	    fprintf(stderr,"[crisflash] No candidates provided (see flags -s and -@). Saving all sgRNAs in genome without off-target validation to %s\n",allsgRNAsFile);
	  }
	  else {
	    fprintf(stderr,"[crisflash] No candidates provided (see flags -s and -@). Exiting!\n");
	    exit(1);
	  }
	}

	// if reference genome, identify and read gRNAs from reference to memory
	if (referenceGenomePath) {
	  if (!vcfPath) {
	    if (upper_case_only) {
	      readFaToTrie_UpperCaseOnly(T, referenceGenomePath, allsgRNAsFile, 0);
	    }
	    else {
	      readFaToTrie(T, referenceGenomePath, allsgRNAsFile, 0);
	    }
	  }
	  else {
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
	    printf("[crisflash] Creating phased genomes for %s.\n", vcfPath);
	    phased = VCF_to_genome(referenceGenomePath, vcfPath, output_VCF1, output_VCF2);
	    // printf("============== 1st HAPLOTYPE ==============\n");
	    readFaToTrie_UpperCaseOnly(T, output_VCF1, allsgRNAsFile, append);
	    if (phased == 1)
	      {
		// printf("============== 2nd HAPLOTYPE ==============\n");
		append += T->nr_sequences;
		readFaToTrie_UpperCaseOnly(T, output_VCF2, allsgRNAsFile, append); // we add haplotype 2
	      }
	    free(output_VCF1);
	    free(output_VCF2);
	  }	  
	}
	else if (bedFilePath) {
	  fprintf(stdout,"[crisflash] Reading genomic sgRNAs from %s ...\n", bedFilePath);
	  T = Bed2Trie(bedFilePath);
	}

	// If no candidates to validate, exit here
	if (allsgRNAsFile && !sequencePath && !bedCandidatesPath) {
	  fprintf(stdout,"[crisflash] Freeing memory for shutdown ... ");
	  fflush(stdout);
	  T = TrieDestroy(T);
	  fprintf(stdout," DONE!\n\n");
	  return 0;
	}

	fprintf(stdout,"[crisflash] Ready to process candidates.\n");
	fflush(stdout);
	
	/* Identify candidate sgRNAs and score off-targets */
	if (bedCandidatesPath) {
	  if (threadFlag) {
	    TrieAMatchItself_thread(T, maxMismatch, bedCandidatesPath, outFile, nr_of_threads);
	  }
	  else {
	    TrieAMatchItself(T, maxMismatch, bedCandidatesPath, outFile);
	  }
	}
	if (sequencePath) {
	  if (threadFlag) { TrieAMatchSequence_thread(T, sequencePath, maxMismatch, outFile, nr_of_threads); }
	  else { TrieAMatchSequence(T, sequencePath, maxMismatch, outFile); }
	}
	
	/* Remove sgRNA trie from memory for shut down */
	// if (referenceGenomePath || bedFilePath) { T = TrieDestroy(T); }
	fprintf(stdout,"[crisflash] Freeing memory for shutdown ... ");
	fflush(stdout);
	T = TrieDestroy(T);
	/*
	if(allsgRNAsFile) { free(allsgRNAsFile); }
	if(bedCandidatesPath) { free(bedCandidatesPath); }
	if(outFile) { free(outFile); }
	*/
	fprintf(stdout," DONE!\n\n");
	return 0;
}
