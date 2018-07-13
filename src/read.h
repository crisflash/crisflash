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

#ifndef DEF_READ
#define DEF_READ

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "nary_tree.h"

// #define BUFFER_SIZE (1 * 1024 * 1024) // at least > 23
#define BUFFER_SIZE 50
#define ITERATIONS (10 * 1024)
#define LENGTH_SEQ (20)
#define CHR_NUMBER_MAX 1024

/*************************/
/*      STRUCTURES       */
/*************************/

typedef int bool;
#define true 1
#define false 0

// Define structure useful for thread calculations
struct arg_struct
{
	int i;
	trie* T;
	int cnt; // counter
	// For the output line =
	size_t read; // size of the line
	char* chrom;
	char* start;
	char* end;
	char* gRNA;
	char strand;
	//======================
	FILE* output;
	int nr_of_mismatches;
};

/*************************/
/*      FONCTIONS        */
/*************************/

double now();
char *invSeq(char* seq, int pos, int slen, int bytes_read, int upper_case_only);
char *fSeq(char* seq, int pos, int slen, int upper_case_only);
void print_seq(char *seq, int len);
void writeInBed(FILE* f, mcontainer *m, char* gRNA, trie* T, int scoreFlag);
void writeInBed_lighter(FILE* f, mcontainer *m, char* chrom, char* start, char* end, char* gRNA, char strand, trie* T);
int newPAM(char* buffer, int x, int bytes_read, int upper_case_only);
void readFaToTrie(trie* T, char *fname, char* outputName, int append);
void readFaToTrie_UpperCaseOnly(trie* T, char *fname, char* outputName, int append);
trie *readFaToTrie_VCF(char *fname, char *vcfName, char* outputName);
int check_var(VCF* vcf, int chr_pos, char* buffer, int x, char* chr);
char* compute_seq(VCF** vcf, int pt_vcf, int chr_pos, char* buffer, int x);
trie* Bed2Trie(char* fname);
int countLines(char* fname);
void TrieAMatchItself(trie* T, int maxMismatch, char* outputName, char* outputMatchResult);
void *thread_worker(void *arg);
void TrieAMatchItself_thread(trie *T, int maxMismatch, char* outputName, char* outputMatchResult, int nr_of_threads);
void TrieAMatchSequence(trie* T, char* fname, int maxMismatch, char* outputName);
void TrieAMatchSequence_thread(trie* T, char* fname, int maxMismatch, char* outputName, int nr_of_threads);
int ClientServer(char* sequencePath, int maxMismatch, int nr_of_threads);

#endif
