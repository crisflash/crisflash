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


#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <ctype.h>
#include <sys/socket.h>
#include <arpa/inet.h>

#include "read.h"
#include "nary_tree.h"
#include "vcf.h"

/* Written by Adrien Jacquin, April 2017
	based on a Margus Lukk script     */

double now()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1000000.;
}

char *invSeq(char* seq, int pos, int slen, int bytes_read, int upper_case_only)
{
	// Transposes DNA sequence, copy seq[pos : pos+slen-1] in seq2 (with the reverse nt order)
	// if upper_case_only == 0, we don't care about distinguishing upper and lower case
	// if upper_case_only == 1, we consider only upper case
	int i=0;
	int j=slen-1;
	char *seq2 = malloc(sizeof(char)*(slen+1));
	seq2[slen] = '\0';
	while(j >= 0)
	{
		if ((pos+slen) > bytes_read) { free(seq2); return NULL; } // if we exit the buffer (right side)
		if (seq[pos + i] == 'C') { seq2[j] = 'G'; i++; j--; continue; }
		if ((seq[pos + i] == 'c') && (upper_case_only == 0)) { seq2[j] = 'G'; i++; j--; continue; }
		if (seq[pos + i] == 'G') { seq2[j] = 'C'; i++; j--; continue; }
		if ((seq[pos + i] == 'g') && (upper_case_only == 0)) { seq2[j] = 'C'; i++; j--; continue; }
		if (seq[pos + i] == 'A') { seq2[j] = 'T'; i++; j--; continue; }
		if ((seq[pos + i] == 'a') && (upper_case_only == 0)) { seq2[j] = 'T'; i++; j--; continue; }
		if (seq[pos + i] == 'T') { seq2[j] = 'A'; i++; j--; continue; }
		if ((seq[pos + i] == 't') && (upper_case_only == 0)) { seq2[j] = 'A'; i++; j--; continue; }
		if (seq[pos + i] == '\n' || seq[pos + i] == '-' || seq[pos + i] == '+') { i++; slen +=1; continue; }
		// if running into any other characters, including Ns etc then break!
		else { free(seq2); seq2 = NULL; break; }
	}
	return seq2;
}

char *fSeq(char* seq, int pos, int slen, int upper_case_only)
{
	// Format DNA sequence, copy seq[pos-slen+1 : pos] in seq2 (with the same nt order)
	// if upper_case_only == 0, we don't care about distinguishing upper and lower case
	// if upper_case_only == 1, we consider only upper case
	int i=0;
	int j=slen-1;
	char *seq2 = (char *) malloc(sizeof(char)*(slen+1));
	seq2[slen] = '\0';
	while(j >=0)
	{
		if (pos - slen < -1) { free(seq2); return NULL; } // if we exit the buffer (left side)
		if (seq[pos + i] == 'C') { seq2[j] = 'C'; i--; j--; continue; }
		if ((seq[pos + i] == 'c') && (upper_case_only == 0)) { seq2[j] = 'C'; i--; j--; continue; }
		if (seq[pos + i] == 'G') { seq2[j] = 'G'; i--; j--; continue; }
		if ((seq[pos + i] == 'g') && (upper_case_only == 0)) { seq2[j] = 'G'; i--; j--; continue; }
		if (seq[pos + i] == 'A') { seq2[j] = 'A'; i--; j--; continue; }
		if ((seq[pos + i] == 'a') && (upper_case_only == 0)) { seq2[j] = 'A'; i--; j--; continue; }
		if (seq[pos + i] == 'T') { seq2[j] = 'T'; i--; j--; continue; }
		if ((seq[pos + i] == 't') && (upper_case_only == 0)) { seq2[j] = 'T'; i--; j--; continue; }
		if (seq[pos + i] == '\n' || seq[pos + i] == '-' || seq[pos + i] == '+') { i--; slen +=1; continue; }
		// if running into any other characters, including Ns etc then break!
		else { free(seq2); seq2 = NULL; break; }
	}
	return seq2;
}

void print_seq(char *seq, int len)
{
	for (int i=0; i < len; i++)
	{
		printf("%c", seq[i]);
		i++;
	}
	printf("\n");
}

void writeInBed(FILE* f, mcontainer *m, char* gRNA, trie* T, int scoreFlag)
{
	// chrom    chromStart    chromEnd    name(sequence gRNA /  genome gRNA)    score(or nrMismatch)    strand
	if (scoreFlag == 0)
	{
		for (int i = 0; i < m->pos; i++)
		{
			// Note: the start position is the first base, and the end position is the base after the last base
			// Note: the position begins at 0 for bed files (this explains the -1)
			for (int j = 0; j < m->nodes[i]->hits; j++)
			{
				fprintf(f, "%s\t%d\t%d\t%s/%s\t%d\t%c\n", T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j]-1, m->nodes[i]->starts[j]+20, gRNA, m->sarray[i], m->marray[i], m->nodes[i]->strands[j]);
			}
		}
	}
	else
	{
		for (int i = 0; i < m->pos; i++)
		{
			// Note: the start position is the first base, and the end position is the base after the last base
			// Note: the position begins at 0 for bed files (this explains the -1)
			for (int j = 0; j < m->nodes[i]->hits; j++)
			{
				fprintf(f, "%s\t%d\t%d\t%s/%s\t%f\t%c\n", T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j]-1, m->nodes[i]->starts[j]+20, gRNA, m->sarray[i], m->off_scores[i], m->nodes[i]->strands[j]);
			}
		}
	}
}

void writeInBed_lighter(FILE* f, mcontainer *m, char* chrom, char* start, char* end, char* gRNA, char strand, trie* T)
{
	// chrom    chromStart    chromEnd    comment(sequence gRNA /  nr of exact matches)    score    strand
	int nr_of_exact_matches = 0;
	int nr_of_off_targets = 0;
	for (int i = 0; i < m->pos; i++)
	{
		if (m->marray[i] == 0) // no mismatch
		{
		  nr_of_exact_matches += m->nodes[i]->hits;
		}
		else {
		  nr_of_off_targets += m->nodes[i]->hits;
		}
	}
	fprintf(f, "%s\t%s\t%s\t%s/%d/%d\t%f\t%c\n", chrom, start, end, gRNA, nr_of_exact_matches, nr_of_off_targets, m->score, strand);
}

int newPAM(char* buffer, int x, int bytes_read, int upper_case_only)
{
	// Return the position of the last G or last C if there is a PAM, 0 otherwise
	int n = 0;
	if (upper_case_only == 1)
	{
		if (x+2 >= bytes_read) { return 0; } // we exit the buffer
		if (buffer[x] == 'G' || buffer[x] == 'C')
		{
			if (buffer[x+1] == buffer[x]) { return x+1; }
			if ((buffer[x+1] == '\n' || buffer[x+1] == '+') && buffer[x+2] == buffer[x]) { return x+2; }
		}
		return 0;
	}
	if (toupper(buffer[x]) == 'G' || toupper(buffer[x]) == 'C')
	{
		if (buffer[x+1] == buffer[x]) { return x+1; }
		if ((buffer[x+1] == '\n' || buffer[x+1] == '+') && buffer[x+2] == buffer[x]) { return x+2; }
	}
	return 0;
}

void readFaToTrie(trie* T, char *fname, char* outputName, int append)
{
	// Read a fasta file and return a trie built with all the gRNAs found
	// Write the gRNAs in a BED file
	// If append = 1, append the gRNAs in the BED file, otherwise erase the existing file before writing
	char *buffer;
	char *chr; // chromosome name
	char* seq;
	int chr_number = 0; // number of current chromosome. Note that chromosome numbers are extacted from fasta header lines starting with '>'. Chromosome numbers run from 0 to n.
	int chr_pos = 0; // position in chromosome.
	int nr_sequences = 0;
	double end_time;
	double total_time;
	int x;
	int n; // position of the last letter in the PAM
	int add; // 1 if sequence added into the tree, 0 otherwise
	ssize_t z = 0; // total bytes read
	double start_time = now();
	ssize_t bytes_read;
	bool bool_chr = false;
	bool bool_last_read = false;
	int identical_seq = 0; // counts the sequences not added in the tree because they were already there

	buffer = calloc(BUFFER_SIZE + 1, sizeof(char));
	buffer[BUFFER_SIZE] = '\0';
	chr = (char *) malloc(128*sizeof(char));
	int chrp = 0; // position in the chromosome name

	FILE* f;

	if (outputName)
	{
	  if (append == 1) { f = fopen(outputName,"a"); }
	  else { f = fopen(outputName,"w"); }
	  if (!f)
	    {
	      fprintf(stderr,"[crisflash] Failed to open %s for writing. Exiting!\n", outputName);
	      exit(1);
	    }
	  fprintf(stdout,"[crisflash] sgRNAs in genome will be written to \'%s\'\n",outputName);
	}

	int fd;
	fd = open(fname, O_RDONLY); // fd is the file descriptor of the file fname

	// Note about the read function: ssize_t read(int fildes, void *buffer, size_t nbyte). Returns the size of the reading.
	// It reads nbyte bytes from the file associated with the open file descriptor, fildes, into the buffer.
	bytes_read = read(fd, buffer, BUFFER_SIZE);
	if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
	{
		bool_last_read = true;
	}

	while(bytes_read > 0)
	{
		for(x = 0; x < bytes_read; x += 1)
		{
			if (buffer[x] == '\n')
                        {
                                if (bool_chr) // end of a chromosome line
                                {
                                        bool_chr = false;
                                        chr[chrp] = '\0';
                                        chrp = 0;
                                        // Add chromosome name to the list of names
					chr_number = addChr(T, chr, strlen(chr)); // we get the chr number in order of appearence in fasta file
                                        chr_pos = 1; // set reading position in chromosome to 1
				}
                                continue;
			}
			else if(!bool_chr)
			{
				n = newPAM(buffer, x, bytes_read, 0); // n is the position of the last letter in the PAM, or 0
				if (n > 0 && toupper(buffer[x]) == 'G')
				{
					seq = fSeq(buffer, n, 23, 0);
					if ((seq != NULL) && (chr_pos-22 > -1))
					{
					  add = TrieAdd(T, seq, LENGTH_SEQ, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
					  // Add the gRNA in the bed file: (1st base = position 0, end position excluded)
					  // chrom    chromStart    chromEnd    name(gRNA)    score    strand
					  if (outputName && add == 1) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-22, chr_pos+1, seq, 0., '+'); }
					  free(seq);
					}
				}
				if (n > 0 && toupper(buffer[x]) == 'C')
				{
					char *seq = invSeq(buffer, x, 23, bytes_read, 0);
					if ((seq != NULL) && (chr_pos-1 > -1))
					{
					  add = TrieAdd(T, seq, LENGTH_SEQ, chr_number, '-', chr_pos+23, 0, &identical_seq);
					  // Add the gRNA in the bed file: (1st base = position 0, end position excluded)
					  // chrom    chromStart    chromEnd    name(gRNA)    score    strand
					  if (outputName && add == 1) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-1, chr_pos+22, seq, 0., '-'); }
					  free(seq);
					}
				}
				else if (buffer[x] == '>')
				{
					bool_chr = true;
				}
			} // end of else if (!bool_chr)
			else // we write the chromosome name in 'chr' :
			{
				chr[chrp] = buffer[x];
				chrp++;
			}
			chr_pos++;
		}
		z += bytes_read - 23;
		chr_pos -= 23; // we will read again the last 23 nt

		bytes_read = pread(fd, buffer, BUFFER_SIZE, z);
		if (bool_last_read) // that was the last read, let's get out of here!
		{
			z += 23;
			break;
		}
		if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
		{
			bool_last_read = true;
		}
	}
	close(fd);
	if (outputName) {
	  fprintf(stdout, "[crisflash] sgRNAs in genome were written to %s.\n", outputName);
	  fclose(f);
	}

	free(buffer);
	free(chr);

	end_time = now();
	total_time = end_time - start_time;

	fprintf(stdout,"[crisflash] %d sgRNA sites in genome. (%f seconds)\n", T->nr_sequences, total_time);
	fflush(stdout);
}

void readFaToTrie_UpperCaseOnly(trie* T, char *fname, char* outputName, int append)
{

        // Read a fasta file and complete trie with all the gRNAs, considering only upper case sequences
	// Write the gRNAs in a BED file
	// If append > 0, append the gRNAs in the BED file, otherwise erase the existing file before writing
	char* seq;
	int chr_number = 0; // number of current chromosome. Note that chromosome numbers are extacted from fasta header lines starting with '>'. Chromosome numbers run from 0 to n.
	int chr_pos = 0; // position in chromosome.
	int nr_sequences = 0;
	double end_time;
	double total_time;
	int n; // position of the last letter in the new found PAM
	int add; // 1 if sequence added into the tree, 0 otherwise
	double start_time = now();
	int identical_seq = 0; // counts the sequences not added in the tree because they were already there

	int chr_positions[CHR_NUMBER_MAX]; // the list of positions of the chromosomes in the genome file
	int chr_sizes[CHR_NUMBER_MAX]; // the list of size of the chromosomes in the genome file
	char* chr_names[CHR_NUMBER_MAX]; // the list of chromosome names in the genome file	
	int m = find_chromosome_positions_ref(fname,chr_positions,chr_sizes,chr_names); // n is the tables size
	
	FILE* f;
	char* output = NULL;

	fflush(stdout);
	
	if (outputName)
	  {
	    if (append > 0) { f = fopen(outputName,"a"); }
	    else { f = fopen(outputName,"w"); }
	    if (!f)
	      {
		fprintf(stderr,"[crisflash] ERROR: failed to open %s for writing. Exiting!\n", outputName);
		exit(1);
	      }
	    fprintf(stdout,"[crisflash] sgRNAs in genome will be written to \'%s\'\n",outputName);
	  }

	FILE* g = fopen(fname,"r"); // genome file

	for (int i = 0; i < m; i++)
	{
	        /* We read entirely the chromosome i */
	        fprintf(stdout,"[crisflash] Processing chromosome %s ...\n", chr_names[i]);
	  
		fseek(g, chr_positions[i], SEEK_SET);
		char* buffer = malloc(sizeof(char)*(chr_sizes[i]+1));
		buffer[chr_sizes[i]] = '\0';
		int read = fread(buffer,chr_sizes[i],1,g);
		if (read == -1) { fprintf(stderr,"[crisflash] ERROR in reading reference genome from %s\n", fname); exit(1); }

		int x = 1;
		char* chr = malloc(sizeof(char)*256);
		while (buffer[x] != '\n') // we know we are at the beginning of a chromosome
		{
			chr[x-1] = buffer[x]; // we store chromosome name without > and \n
			x++;
		}
		chr[x-1] = '\0';
		chr_number = addChr(T, chr, strlen(chr)); // we get the chr number in order of appearence in fasta file
		int chr_pos = 1; // initializing position in chromosome
		x++; // first letter of the genome
		
		/* We read all the nucelotides of the chromosome */
		while (x < chr_sizes[i])
		{
			if (buffer[x] == '\n') { x++; continue; }
			int n = 0; // n is the position of the last letter in the PAM, or 0
			if (buffer[x] == '+')
			{
				if (buffer[x+1] == 'G' || buffer[x+1] == 'C')
				{
					if (x+2 < chr_sizes[i] && (buffer[x+2] == '+' || buffer[x+2] == '\n') && buffer[x+3] == buffer[x+1]) { n = x+3; }
					else if (x+2 < chr_sizes[i] && buffer[x+2] == buffer[x+1]) { n = x+2; }
				}
				x++; // we go on the letter
			}
			else if (x+1 < chr_sizes[i] && (buffer[x] == 'G' || buffer[x] == 'C') && (buffer[x+1] == '-'))
			{
				int i = 1;
				while (buffer[x+i] == '-') { i++; }
				if (buffer[x+i] == buffer[x]) { n = x+i; }
			}
			else { n = newPAM(buffer, x, chr_sizes[i], 1); }
			if (n > 0 && buffer[x] == 'G')
			{
				seq = fSeq(buffer, n, 23, 1);
				if (seq != NULL)
				{
				  add = TrieAdd(T, seq, LENGTH_SEQ, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
				  // Add the gRNA in the bed file: (1st base = position 0, end position excluded)
				  // chrom    chromStart    chromEnd    name(gRNA)    score    strand
				  if (outputName && add == 1) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-22, chr_pos+1, seq, 0., '+'); }
				  free(seq);
				}
			}
			if (n > 0 && buffer[x] == 'C')
			{
				char *seq = invSeq(buffer, x, 23, chr_sizes[i], 1);
				if (seq != NULL)
				{
				  add = TrieAdd(T, seq, LENGTH_SEQ, chr_number, '-', chr_pos+23, 0, &identical_seq);				    
				  // Add the gRNA in the bed file: (1st base = position 0, end position excluded)
				  // chrom    chromStart    chromEnd    name(gRNA)    score    strand
				  if (outputName && add == 1) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-1, chr_pos+22, seq, 0., '-'); }
				  free(seq);
				}
			}
			x++;
			if (buffer[x-2] == '+') { continue; }
			chr_pos++;
		}
		free(buffer);
		free(chr);
		buffer = NULL;
		chr=NULL;
	}
	fclose(g);
	if (outputName) {
	  fprintf(stdout, "[crisflash] sgRNAs in genome were written to %s.\n", outputName);
	  fclose(f);
	}
	for (int i=0; i<m; i++) { free(chr_names[i]); }

	end_time = now();
	total_time = end_time - start_time;

	fprintf(stdout,"[crisflash] %d sgRNA sites in genome. (%f seconds)\n", T->nr_sequences, total_time);
	fflush(stdout);
}

// FIXME
// This function is meant to parse the reference genome and a VCF file in the same time, to build the tree
// For now it's only working with SNPs in positive strands
// It still has to be implemented in negative strands and for indels
trie *readFaToTrie_VCF(char *fname, char *vcfName, char* outputName)
{
	// Returns a tree with every potential protospacers of the reference genome (fname),
	// considering variations of the VCF file:
	// - if the var deletes the PAM ==> 0 branch added
	// - if the variation occurs in the protospacer and is present in both haplotypes, or non of them ==> 1 branch added
	// - if the variation occurs in the protospacer and in one of the two haplotypes ==> 2 branches added
	char *buffer;
	char *chr; // chromosome name
	char* seq;
	int chr_number = 0; // number of current chromosome. Note that chromosome numbers are extacted from fasta header lines starting with '>'. Chromosome numbers run from 0 to n.
	int chr_pos = 0; // position in chromosome.
	int nr_sequences = 0;
	double end_time;
	double total_time;
	int x;
	ssize_t z = 0; // total bytes read
	double start_time = now();
	ssize_t bytes_read;
	bool bool_chr = false;
	bool bool_last_read = false;
	int end_vcf = 0;
	int identical_seq = 0; // counts the sequences not added in the tree because they were already there

	buffer = calloc(BUFFER_SIZE + 1, sizeof(char));
	buffer[BUFFER_SIZE] = '\0';
	chr = malloc(128*sizeof(char));
	int chrp = 0; // position in the chromosome name

	trie* T = TrieCreate(LENGTH_SEQ);
	int fd_ref = open(fname,O_RDONLY);
	FILE* f_vcf = fopen(vcfName,"r");
	FILE* f;

	/* Create output file */
	if (outputName)
	{
	  f = fopen(outputName,"w");
	  if (!f)
	    {
	      fprintf(stderr,"[crisflash] ERROR: failed to open %s for writing. Exiting!\n", outputName);
	      exit(1);	     
	    }
	  fprintf(stdout,"[crisflash] sgRNAs in genome will be written to \'%s\'\n",outputName);
	}
	
	/* Install the VCF structure */
	VCF* vcf = VCF_install(f_vcf); // create vcf and store first line
	int previous_pos = 0; // we memorize the previous pos in VCF file for future PAM

	/* Install the buffer */
	bytes_read = read(fd_ref, buffer, BUFFER_SIZE);
	if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
	{
		bool_last_read = true;
	}

	while(bytes_read > 0)
	{
		for(x = 0; x < bytes_read; x += 1)
		{
			if (buffer[x] == '\n')
			{
				if (bool_chr) // end of a chromosome line
				{
					bool_chr = false;
					chr[chrp] = '\0';
					chrp = 0;
					// Add chromosome name to the list of names
					chr_number = addChr(T, chr, strlen(chr)); // we get the chr number in order of appearence in fasta file
					chr_pos = 1; // set reading position in chromosome to 1
					while (strcmp(vcf->chr,chr) != 0 && end_vcf != -1) { end_vcf = VCF_update(f_vcf,vcf); }
					continue;
				}
				else { continue; }
			}
			else if(!bool_chr)
			{
				if ((buffer[x] == 'G' && buffer[x+1] == 'G') || (buffer[x] == 'G' && buffer[x+1] == '\n' && buffer[x+2] == 'G'))
				{
					// we found a PAM, we store every related variations in a table
					VCF* t_vcf[10];
					int pt_vcf = 0; // position in t_vcf
					while (end_vcf == 0 && vcf->ub < chr_pos-21 && strcmp(vcf->chr,chr) == 0)
					{
						end_vcf = VCF_update(f_vcf,vcf);
						previous_pos = ftell(f_vcf);
					}
					while (end_vcf == 0 && check_var(vcf,chr_pos,buffer,x,chr) == 1)
					{
						t_vcf[pt_vcf] = VCF_copy(vcf);
						VCF_print(t_vcf[pt_vcf]);
						pt_vcf++;
						end_vcf = VCF_update(f_vcf,vcf);
					}
					fseek(f_vcf, previous_pos, SEEK_SET); // we come back to the previous pos for the future PAM
					VCF_destroy(vcf); // we reset vcf
					vcf = VCF_install(f_vcf);
					if (buffer[x+1] == '\n')
					{
						if (pt_vcf > 0) // there is at least one variation in the sequence
						{
							seq = compute_seq(t_vcf,pt_vcf,chr_pos,buffer,x);
							printf("Sequence(s) found and treated with VCF: %s\n", seq);
						}
						else
						{
							seq = fSeq(buffer, x+2, 23, 0);
						}
					}
					else
					{
						if (pt_vcf > 0) // there is at least one variation in the sequence
						{
							seq = compute_seq(t_vcf,pt_vcf,chr_pos,buffer,x);
							printf("Sequence(s) found and treated with VCF: %s\n", seq);
						}
						else
						{
							seq = fSeq(buffer, x+1, 23, 0);
						}
					}
					if ((seq != NULL) && (chr_pos-23 > -1))
					{
					  TrieAdd(T, seq, LENGTH_SEQ, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
					  // Add the gRNA in the bed file: (1st base = position 0, end position excluded)
					  // chrom    chromStart    chromEnd    name(gRNA)    score    strand
					  if (outputName) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-22, chr_pos+1, seq, 0., '+'); }
					  // fprintf(f, "%s\n", seq);
					  free(seq);
					  T->nr_sequences++;
					}
				}
				// if ((buffer[x] == 'C' && buffer[x+1] == 'C') || (buffer[x] == 'C' && buffer[x+1] == '\n' && buffer[x+2] == 'C'))
				// {
				// 	char *seq = invSeq(buffer, x, 23, bytes_read, 1);
				// 	if ((seq != NULL) && (chr_pos-2 > -1))
				// 	{
				// 		// print_seq(seq, 20);
				// 		// print_seq(seq, 23);
				// 		TrieAdd(T, seq, LENGTH_SEQ, chr_number, '-', chr_pos+23, 0, &identical_seq);
				// 		// Add the gRNA in the bed file: (1st base = position 0, end position excluded)
				// 		// chrom    chromStart    chromEnd    name(gRNA)    strand
				// 		if (outputName) { fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", chr, chr_pos-1, chr_pos+22, seq, 0., '-'); }
				// 		// fprintf(f, "%s\n", seq);
				// 		free(seq);
				// 		T->nr_sequences++;
				// 	}
				// }
				else if (buffer[x] == '>')
				{
					bool_chr = true;
				}
			} // end of else if (!bool_chr)
			else // we write the chromosome name in 'chr' :
			{
				chr[chrp] = buffer[x];
				chrp++;
			}
			chr_pos++;
		}
		z += bytes_read - 23;
		chr_pos -= 23; // we will read again the last 23 nt

		bytes_read = pread(fd_ref, buffer, BUFFER_SIZE, z);
		if (bool_last_read) // that was the last read, let's get out of here!
		{
			z += 23;
			break;
		}
		if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
		{
			bool_last_read = true;
		}
	}
	fclose(f);
	fclose(f_vcf);
	close(fd_ref);
	free(buffer);
	free(chr);
	VCF_destroy(vcf);

	end_time = now();
	total_time = end_time - start_time;

	fprintf(stdout,"[crisflash] %d sgRNA sites in genome. (%f seconds)\n", T->nr_sequences, total_time);
	if (outputName)
	{
	  fprintf(stdout, "[crisflash] sgRNAs in genome were written to %s.\n", outputName);
	}
	fflush(stdout);	

	return T;
}

// FIXME
// This function is not finished, it goes with the previous one
int check_var(VCF* vcf, int chr_pos, char* buffer, int x, char* chr)
{
	// Returns 1 if there is a variation in the new found sequence, 0 otherwise
	if (strcmp(vcf->hap,"0|0") == 0 || strcmp(vcf->chr,chr) != 0) { return 0; }
	if (vcf->dif == 0) // variation = SNP
	{
		if (vcf->pos > chr_pos-22 && vcf->pos < chr_pos+2)
		{
			// we still have to check if the reference allele is the same in the reference genome
			int i = chr_pos - vcf->pos;

			int j = 0;
			if (i > j) // we have to take into account the carriage returns when we move in the buffer
			{
				while (j < i)
				{
					if (buffer[x-j] == '\n') { i++; j++; }
					j++;
				}
			}
			else
			{
				while (j > i)
				{
					if (buffer[x-j] == '\n') { i--; j--; }
					j--;
				}
			}
			// printf("%c %c\n", buffer[x-i], vcf->ref[0]);
			if (buffer[x-i] == (vcf->ref)[0]) { return 1; }
			return 0;
		}
		return 0;
	}
	// // variation = indel
	// if (vcf.ub > chr_pos-22 && vcf.pos+1 < chr_pos+2) // we add 1 to vcf.pos because in indels the first position doesn't change
	// {
	// 	// we still have to check if the reference allele is the same in the reference genome
	// 	int i = chr_pos - vcf.pos;
	// 	int j = 0;
	// 	if (i > j) // we have to take into account the carriage returns when we move in the buffer
	// 	{
	// 		while (j < i)
	// 		{
	// 			if (buffer[x-j] == '\n') { i++; j++; }
	// 			j++;
	// 		}
	// 	}
	// 	else
	// 	{
	// 		while (j > i)
	// 		{
	// 			if (buffer[x-j] == '\n') { i--; j--; }
	// 			j--;
	// 		}
	// 	}
	// 	if (buffer[x-i] == (vcf.ref)[0]) { return 1; }
	// 	return 0;
	// }
	return 0;
}

char* compute_seq(VCF** t_vcf, int pt_vcf, int chr_pos, char* buffer, int x)
{
	// Compute the sequence(s) in the positive strand considering the variation(s) in vcf
	int j = x+1; // where we will begin to copy the nucleotides
	int k = 22; // index for new string
	int n = pt_vcf-1;
	char* seq = malloc(sizeof(char)*(24)); // new string
	seq[23] = '\0';
	if (buffer[x+1] == '\n') { j = x+2; }
	while (n >= 0)
	{
		if (t_vcf[n]->dif == 0) // SNP
		{
			int e = t_vcf[n]->pos - (chr_pos - 21); // if e < 0, it means beginning of var is before the sequence
			while (k >= 0)
			{
				if (buffer[j] == '\n') { j--; }
				if (k == e) // there is a variation here
				{
					seq[k] = t_vcf[n]->alt[0];
					j--;
					k--;
					VCF_destroy(t_vcf[n]);
					n--;
					// we made the variation, we have to check if there is another one:
					if (n >= 0) { break; }
				}
				seq[k] = buffer[j];
				j--;
				k--;
			}
			if (k >= 0) // we exited because of a variation and there is another one
			{
				continue; // we go back to the beginning of the while
			}
			else { return seq; }
		}
	}
	return NULL;
}

trie* Bed2Trie(char* fname)
{
	// Return a trie based on a BED files that contains all sequences to add in the trie
	// The goal is to save time
	trie* T = TrieCreate(LENGTH_SEQ);
	int chr_number;
	int start_int;
	int identical_seq = 0; // counts the sequences not added in the tree because they were already there

	FILE *f = fopen(fname,"r");
	if (!f)
	{
	  fprintf(stderr,"[crisflash] ERROR: unabole to open %s!\n", fname);
	  exit(1);
	}
	char* line = NULL;
	size_t len = 0;
	ssize_t read;

	char* chr = malloc(sizeof(char)*128);
	char* start = malloc(sizeof(char)*32);
	char* seq = malloc(sizeof(char)*32);
	int i = 0;
	int j = 0;
	
	double end_time, total_time;
	double start_time = now();	

	fprintf(stdout,"[crisflash] Reading %s ...\n", fname);

	while ((read = getline(&line, &len, f)) != -1)
	  {
	        int i = 0;
		int j = 0;
		
		while (line[i] != '\t')
		{
			chr[i] = line[i];
			i++;
		}
		chr[i] = '\0';
		i++;
		chr_number = addChr(T, chr, i);
		while (line[i] != '\t')
		{
			start[j] = line[i];
			i++;
			j++;
		}
		start[j] = '\0';
		i++;
		j = 0;
		start_int = atoi(start) + 1;
		while (line[i] != '\t') { i++; }
		i++;
		while (line[i] != '\t')
		{
			seq[j] = line[i];
			i++;
			j++;
		}
		seq[j] = '\0';
		i++; // now we are on the score
		while (line[i] != '\t') { i++; }
		i++;
		TrieAdd(T, seq, LENGTH_SEQ, chr_number, line[i], start_int, 0, &identical_seq);
	}

	fclose(f);
	if (line)
	{
		free(line);
	}
	free(seq);
	free(chr);
	free(start);

	end_time = now();
	total_time = end_time - start_time;

	fprintf(stdout,"[crisflash] %d sequences read from file. (%f seconds)\n", T->nr_sequences, total_time);
	fflush(stdout);
	
	return T;
}

int countLines(char* fname)
{
	FILE* f = fopen(fname, "r");
	char ch;
	char prev_ch = ch;
	int notEmpty = 0;
	int cnt = 0;
	if (f)
	{
		while ((ch=getc(f)) != EOF)
		{
			if (notEmpty == 0) { notEmpty = 1; }
			if (ch == '\n') { cnt++; }
			prev_ch = ch;
		}
		if (notEmpty) { cnt++; }
		if (prev_ch == '\n') { cnt--; } // it means there is a useless carriage return just before the end of the file
		fclose(f);
	}
	else
	{
	  fprintf(stderr,"[crisflash] ERROR: unabole to open %s!\n", fname);
	  exit(1);
	}
	return cnt;
}

void TrieAMatchItself(trie* T, int maxMismatch, char* outputName, char* outputMatchResult)
{
	// This function directly reads in the output file of readFaToTrie
	// It reads every gRNA one by one and match them against the whole trie
	// It writes the results in an output BED file
	double end_time;
	double start_time = now();
	int cnt = 0;
	mcontainer* m;
	int i = 0;
	
	fprintf(stdout,"Reading %s for candidate validation.\n",outputName);
	fflush(stdout);
	int nr_of_lines = countLines(outputName);
	if (nr_of_lines == 0)
	  {
	    fprintf(stderr, "[crisflash] ERROR: %s is an empty file!\n", outputName);
	    exit(1);
	  }
		
	FILE* f = fopen(outputName,"r");
	if (!f)
	{
	  fprintf(stderr,"[crisflash] ERROR: failed to open %s!\n",outputName);
	  exit(1);
	}
	
	FILE *f2 = fopen(outputMatchResult,"w");
	if (!f2)
	{
	  fprintf(stderr,"[crisflash] ERROR: failed to open '%s' for writing!\n", outputMatchResult);
	  exit(1);
	}
	fprintf(stdout,"[crisflash] Saving results to '%s'\n", outputMatchResult);
	fflush(stdout);

	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	while ((read = getline(&line, &len, f)) != -1)
	{
		int k = 0;
		char* chrom = malloc(sizeof(char)*128);
		char* start = malloc(sizeof(char)*32);
		char* end = malloc(sizeof(char)*32);
		char* gRNA = malloc(sizeof(char) * (T->readlen + 1));
		char strand;
		for (int j = 0; j < 3; j++)
		{
			int l = 0; // position in the new char*
			while (line[k] != '\t')
			{
				if (j==0) { chrom[l] = line[k]; }
				if (j==1) { start[l] = line[k]; }
				if (j==2) { end[l] = line[k]; }
				k++;
				l++;
			}
			k++;
			if (j==0) { chrom[l] = '\0'; }
			if (j==1) { start[l] = '\0'; }
			if (j==2) { end[l] = '\0'; }
			l = 0;
		}
		for (int j = 0; j < LENGTH_SEQ; j++)
		{
			gRNA[j] = line[k];
			k++;
		}
		gRNA[LENGTH_SEQ] = '\0';
		while (line[k] != '\t') { k++; }
		k++;
		strand = line[k];
		m = TrieAMatch(T, gRNA, T->readlen, maxMismatch);
		if (f2 != NULL)
		{
			writeInBed_lighter(f2, m, chrom, start, end, gRNA, strand, T); // write the results in a .bed file
		}
		// print match results
		// printf("Sequence to match = %s\nSequence(s) which matched:\n", gRNA);
		// mcontainer_print(m);
		// mcontainer_print_pretty(m, maxMismatch); // the i-th element shows the number of sequences with i mismatches
		// mcontainer_free(m);
		cnt++;
		// printf("%d\n", cnt);
		if (cnt % 1000 == 0)
		{
		  fprintf(stdout,"[crisflash] %d gRNAs tested ... [%d%%]\n", cnt, cnt*100/nr_of_lines);
		}
		free(chrom);
		free(start);
		free(end);
		free(gRNA);
		chrom = NULL;
		start = NULL;
		end = NULL;
		gRNA = NULL;
	}
	fclose(f);
	fclose(f2);
	if (line)
	{
		free(line);
	}

	end_time = now();
	fprintf(stdout,"[crisflash] %d sgRNA candidates with up to %d mismatches tested in total. (%f seconds).\n", cnt, maxMismatch, (end_time - start_time));
	fprintf(stdout,"[crisflash] Results saved in %s\n", outputMatchResult);
}

void *thread_worker(void *arg)
{
	// cast the argument into its real nature
	struct arg_struct *args = (struct arg_struct *)arg;

	if (args->gRNA && args->read != -1)
	{
		mcontainer* m;

		// match
		m = TrieAMatch(args->T, args->gRNA, LENGTH_SEQ, args->nr_of_mismatches);
		if (args->output != NULL)
		{
			writeInBed_lighter(args->output, m, args->chrom, args->start, args->end, args->gRNA, args->strand, args->T); // write the results in a file
			// printf("Writing results in a file, thread #%d, counter = %d\n", args->i, args->cnt);
		}

		// free memory
		mcontainer_free(m);
	}

	(void) arg;
	pthread_exit(NULL);
}

void TrieAMatchItself_thread(trie *T, int maxMismatch, char* outputName, char* outputMatchResult, int nr_of_threads)
{
	// Open outputName.bed and match its sequences agaisnt the tree T using multiple threads
	// Write the results in outputMatchResult.bed
	double end_time;
	double start_time = now();
	int cnt = 0; // counter

	// install files
	int i = 0;

	int nr_of_lines = countLines(outputName);
	if (nr_of_lines == 0)
	{
	  fprintf(stderr, "[crisflash] ERROR: %s is an empty file!\n", outputName);
	  exit(1);
	}
	
	FILE* f = fopen(outputName,"r");
	if (!f)
	{
	  fprintf(stderr,"[crisflash] ERROR: failed to open %s!\n", outputName);
	  exit(1);
	}
	fprintf(stdout,"[crisflash] Reading %s for candidate validation.\n",outputName);
	fflush(stdout);

	FILE *f2 = fopen(outputMatchResult,"w");
	if (!f2)
	{
	  fprintf(stderr,"[crisflash] ERROR: failed to open '%s' for writing!\n", outputMatchResult);
	  exit(1);
	}
	fprintf(stdout,"[crisflash] Saving results to '%s'\n", outputMatchResult);
	fflush(stdout);

	// install threads
	pthread_t threads[nr_of_threads];
	// install pthread attributes necessary for keeping track and joining the thread from main program.
	pthread_attr_t attr;
	// initialize and set thread detached attribute
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	// table of parameters
	struct arg_struct args_table[nr_of_threads];

	char * line = NULL;
	size_t len = 0;
	ssize_t read = 0;
	void* status;
	int k = 0;
	int end = nr_of_threads;

	while (cnt < T->nr_sequences && read != -1)
	  {
		for (int i=0; i < nr_of_threads; i++)
		{
			read = getline(&line, &len, f);
			// install the parameters
			args_table[i].read = read; // indicates if the line has been read or not
			args_table[i].T = T;
			args_table[i].i = i;
			args_table[i].chrom = malloc(sizeof(char)*128);
			args_table[i].start = malloc(sizeof(char)*32);
			args_table[i].end = malloc(sizeof(char)*32);
			args_table[i].gRNA = malloc(sizeof(char) * (LENGTH_SEQ + 1));
			k = 0;
			for (int j = 0; j < 3; j++)
			{
				int l = 0; // position in the new char*
				while (line[k] != '\t')
				{
					if (j==0) { args_table[i].chrom[l] = line[k]; }
					if (j==1) { args_table[i].start[l] = line[k]; }
					if (j==2) { args_table[i].end[l] = line[k]; }
					k++;
					l++;
				}
				k++;
				if (j==0) { args_table[i].chrom[l] = '\0'; }
				if (j==1) { args_table[i].start[l] = '\0'; }
				if (j==2) { args_table[i].end[l] = '\0'; }
				l = 0;
			}
			for (int j = 0; j < LENGTH_SEQ; j++)
			{
				args_table[i].gRNA[j] = line[k];
				k++;
			}
			args_table[i].gRNA[LENGTH_SEQ] = '\0';
			while (line[k] != '\t') { k++; } // we reach the end of the gRNA sequence
			k++;
			while (line[k] != '\t') { k++; } // we read the score (without storing it)
			k++;
			args_table[i].strand = line[k]; // strand
			if (read != -1) { cnt++; }
			args_table[i].cnt = cnt;
			args_table[i].output = f2;
			args_table[i].nr_of_mismatches = maxMismatch;

			// create thread
			if(pthread_create(&threads[i], &attr, thread_worker, &args_table[i]))
			{
				perror("pthread_create");
				exit(1);
			}
			if (cnt % 100000 == 0)
			{
			  fprintf(stdout,"[crisflash] %d gRNAs tested ... [%d%%]\n", cnt, cnt*100/nr_of_lines);
			  fflush(stdout);
			}
			if (cnt == T->nr_sequences)
			{
				end = i + 1;
				break;
			}
		}
		// wait for the other threads
		for (int i = 0; i < end; i++)
		{
			if (pthread_join(threads[i], &status))
			{
				perror("pthread_join");
				exit(1);
			}
		}
		for (int i = 0; i < nr_of_threads; i++)
		{
			free(args_table[i].chrom);
			free(args_table[i].start);
			free(args_table[i].end);
			free(args_table[i].gRNA);
			args_table[i].chrom = NULL;
			args_table[i].start = NULL;
			args_table[i].end = NULL;
			args_table[i].gRNA = NULL;
		}
	}
	// free attribute
	pthread_attr_destroy(&attr);

	if (line)
	{
		free(line);
	}
	fclose(f);
	fclose(f2);

	end_time = now();

	fprintf(stdout,"[crisflash] Candidate evaluation complete. %d sgRNAs with up to %d mismatches tested in %f seconds.\n", cnt, maxMismatch, (end_time - start_time));
	fprintf(stdout,"[crisflash] Results saved in %s\n", outputMatchResult);
	fflush(stdout);
}

void TrieAMatchSequence_thread(trie* T, char* fname, int maxMismatch, char* outputName, int nr_of_threads)
{
	// Creates a second tree and a BED file, based on the sequence file.
	// Calls the multiple threads function to match the sequence against the tree.
	// Writes the results in the BED file outputName.
	double end_time;
	double start_time = now();
	int n;
	int append = 0;

	trie* T2 = TrieCreate(LENGTH_SEQ);

	readFaToTrie_UpperCaseOnly(T2, fname, "sequenceTree", append); // Create a tree and a BED file for the sequence
	TrieAMatchItself_thread(T, maxMismatch, "sequenceTree", outputName, nr_of_threads);
	T2 = TrieDestroy(T2);
	n = remove("sequenceTree.bed"); // delete the BED file we needed temporary
	if (n < 0) { perror("[crisflash] ERROR deleting sequence BED file"); }
}

void TrieAMatchSequence(trie* T, char* fname, int maxMismatch, char* outputName)
{
	// Finds all the gRNAs in a fasta file.
	// It matches each one of them in the trie taken as an argument, with up to maxMismatch mismatches
	// Write the results in a .bed file
	char *buffer;
	char* gRNA; // string to store the gRNA
	char* seq;
	int chr_pos = 0; // position in chromosome
	char* chr; // chromosome name
	int chrp = 0;
	int nrRNA = 0;
	int nrRNA_gen = 0;
	double end_time;
	double total_time;
	int x;
	ssize_t z = 0; // total bytes read
	double start_time = now();
	ssize_t bytes_read;
	bool bool_chr = false;
	bool bool_last_read = false;

	chr = malloc(128*sizeof(char));

	FILE *f = fopen(outputName,"w");
	if (!f)
	{
	  fprintf(stderr,"[crisflash] ERROR: failed to open %s for writing!\n", outputName);
	  exit(1);
	}
	
	mcontainer* m;

	buffer = calloc(BUFFER_SIZE + 1, sizeof(char));
	buffer[BUFFER_SIZE] = '\0';

	int fd;
	fd = open(fname, O_RDONLY); // fd is the file descriptor of the file fname

	// FIRST: we have to find all the gRNAs
	bytes_read = read(fd, buffer, BUFFER_SIZE);
	// If the sequence is smaller than the BUFFER_SIZE :
	if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
	{
		bool_last_read = true;
	}

	while(bytes_read > 0)
	{
		for(x = 0; x < bytes_read; x += 1)
		{
			if (buffer[x] == '\n' && bool_chr) // end of a chromosome line
			{
				bool_chr = false;
				chr[chrp] = '\0';
				chrp = 0;
				chr_pos = 1;
			}
			else if(!bool_chr)
			{
				if ((buffer[x] == 'G' && buffer[x+1] == 'G') || (buffer[x] == 'G' && buffer[x+1] == '\n' && buffer[x+2] == 'G'))
				{
					if (buffer[x+1] == '\n') { seq = fSeq(buffer, x+2, 23, 1); }
					else { seq = fSeq(buffer, x+1, 23, 1); }
					if ((seq != NULL) && (chr_pos-23 > -1))
					{
						gRNA = malloc(sizeof(char)*(LENGTH_SEQ + 1));
						gRNA[LENGTH_SEQ] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<LENGTH_SEQ; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, LENGTH_SEQ, maxMismatch);
							for (int i = 0; i < m->pos; i++)
							{
								nrRNA_gen += m->nodes[i]->hits;
							}
							// mcontainer_print(m);
							// mcontainer_print_pretty(m, maxMismatch);
							if (f != NULL)
							{
								char start[32], end[32];
								sprintf(start, "%d", chr_pos - 23);
								sprintf(end, "%d", chr_pos);
								writeInBed_lighter(f, m, chr, start, end, gRNA, '+', T);
							}
							free(gRNA);
							mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
							nrRNA++;
						}
						free(seq);
					}
				}
				if ((buffer[x] == 'C' && buffer[x+1] == 'C') || (buffer[x] == 'C' && buffer[x+1] == '\n' && buffer[x+2] == 'C'))
				{
					char *seq = invSeq(buffer, x, 23, bytes_read, 1);
					if ((seq != NULL) && (chr_pos-2 > -1))
					{
						gRNA = malloc(sizeof(char)*(LENGTH_SEQ + 1));
						gRNA[LENGTH_SEQ] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<LENGTH_SEQ; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, LENGTH_SEQ, maxMismatch);
							for (int i = 0; i < m->pos; i++)
							{
								nrRNA_gen += m->nodes[i]->hits;
							}
							// mcontainer_print(m);
							// mcontainer_print_pretty(m, maxMismatch);
							if (f != NULL)
							{
								char start[32], end[32];
								sprintf(start, "%d", chr_pos - 2);
								sprintf(end, "%d", chr_pos + 21);
								writeInBed_lighter(f, m, chr, start, end, gRNA, '-', T);
							}
							free(gRNA);
							mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
							nrRNA++;
						}
						free(seq);
					}
				}
				else if (buffer[x] == '>')
				{
					bool_chr = true;
				}
			} // end of else if (!bool_chr)
			else // we write the chromosome name in 'chr' :
			{
				chr[chrp] = buffer[x];
				chrp++;
			}
			chr_pos++;
		}
		z += bytes_read - 23;
		chr_pos -= 23; // we will read again the last 23 nt

		bytes_read = pread(fd, buffer, BUFFER_SIZE, z);
		if (bool_last_read) // that was the last read, let's get out of here!
		{
			z += 23;
			break;
		}
		if (bytes_read < BUFFER_SIZE) // the last read is smaller than the others
		{
			bool_last_read = true;
		}
	}
	close(fd);
	fclose(f);

	free(buffer);

	end_time = now();
	
	fprintf(stdout,"[crisflash] %d sgRNA candidates with up to %d mismatches tested. (%f seconds).\n", nrRNA, maxMismatch, (end_time - start_time));
	fprintf(stdout,"[crisflash] Results saved in %s\n", outputName);
}
