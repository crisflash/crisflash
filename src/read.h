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

#define BUFFER_SIZE 4096
#define ITERATIONS (10 * 1024)
#define CHR_NUMBER_MAX 1024
#define PROTOSPACER_LENGTH 20

/*************************/
/*      STRUCTURES       */
/*************************/

typedef int bool;
#define true 1
#define false 0

typedef struct faread_struct
{
  /*
    The struct is used for reading fasta file one sequence at the time.
    Used by installFastaReader, fastaReader and freed by freeFastaReader.
  */
  long long limit; // maximul length allocated for s and sr at any time
  char *s; // fasta sequence '+' strand
  char *sr; // fasta sequence '-' strand
  long long slen; // string length for s and sr
  int fd; // fasta file handler
  char * buffer; // read buffer, limited to BUFFER_SIZE
  ssize_t bytes_read; // nr of bytes in read buffer
  int i; // character position in read buffer
  char *header; // fasta header, limited to BUFFER_SIZE
  int fasta_header; // 1 if read character is from sequence, 0 otherwise.
  int headerlen; // length of the fasta header
  int nrheaders; // nr of sequences read so far
  int alive; // 1 - if file has not yet reached its end, 0 otherwise.
  int sequence; // 1 - if fasta sequence was read, 0 otherwise.
  int no_low_complexity; // 1 if no low complexity sequence, 0 otherwise.
} faread_struct;

typedef struct grna_list {
  /* The struct is used for storing list of gRNA start positions */
  char *pstrand; // pointer to forward strand of fasta sequence
  char *nstrand; // pointer to reverse strand of fasta sequence (translation from forward)
  char *chr; // chromosome name
  long long length; // number of space/positions allocatd for gRNAs
  long long pos; // nr of gRNas
  long long *starts; // array of gRNA start positions
  char *strands; // array of gRNA strands
} grna_list;

struct arg_struct
{
  /* The struct is used for passing arguments to thread worker responsible for matching gRNA candidate in Trie. */
  int i;
  trie* T;
  int maxMismatch;
  int free;
  FILE *outfh;
  int outFileType;
  char *grna;
  char *chr;
  long long start;
  long long end;
  char strand;
};

/*************************/
/*      FUNCTIONS        */
/*************************/

double now();
void readFaToTrie(trie *T, char *genomefname, char *pam, char* outputfname, int append, int upper_case_only, int printGRNAs);
void TrieAMatchSequenceThreads(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, char *pam, int upper_case_only, int threads, int printOnly);

#endif
