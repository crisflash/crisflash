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
#include <ctype.h>

#include "read.h"
#include "nary_tree.h"
#include "vcf.h"

FILE* open_file(char *fname,char *mode) {
  FILE *fh = fopen(fname, mode);
  if (!fh) {
    fprintf(stderr,"[crisflash] ERROR: failed to open '%s'!\n", fname);
    exit(1);
  }
  return fh;
}

int baseInOther(char seqnt) {
  /* Checks whether char, provided from 'variant edited' haplotype sequence is an editing mark or base pair. Values returned are used for calculating */
  switch(seqnt) {
  case '+':
    return -1; // + in front of a base stands for one extra base in DNA sequence altered by variant data compared to standard reference, hence -1;
    break;
  case '-':
    return +1; // - stands for a base which was deleted in this position in standard reference, hence +1;
    break;
  case '|':
    return 0; // | stands for modified base, compared to standard reference, the number of bases compared to reference is the same.
    break;
  default:
    return 5; // note that this is arbritary value and could, in theory be anything except -1,0,1. However, if changed, corresponding match values from code calling this function should be changed as well. 
    break;
  }
}

int baseMatch(char seqnt, char pamnt) {
  /* Checks whether base in DNA sequence (seqnt) matches with base in PAM pattern (pamnt).
     Assumes uppercase values for pamnt.
     Matching is provided based on  IUPAC representation for a position on a DNA sequence
     and were taken from https://en.wikipedia.org/wiki/Nucleic_acid_notation.
  */

  // first check whether sequence base happens to be one of following: '+', '-', '|'.
  switch(seqnt) {
  case '+':
    return 2;
    break;
  case '-':
    return 2;
    break;
  case '|':
    return 2;
    break;
  default:
    break;
  }

  switch(pamnt)
    {
    case 'N':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	case 'G':
	  return 1;
	  break;
	case 'C':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'G':
      switch(seqnt)
	{
	case 'G':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'A':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	default:
	  return 0;
	  break;	  
	}
      break;
    case 'T':
      switch(seqnt)
	{
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'C':
      switch(seqnt)
	{
	case 'C':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;     
    case 'R':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'G':
	    return 1;
	    break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'Y':
      switch(seqnt)
	{
	case 'C':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'V':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'C':
	  return 1;
	  break;
	case 'G':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;    
    case 'W':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'S':
	switch(seqnt)
	  {
	  case 'G':
	    return 1;
	    break;
	  case 'T':
	    return 1;
	    break;
	  default:
	    return 0;
	    break;
	  }
	break;
    case 'M':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'C':
	  return 1;
	    break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'K':
      switch(seqnt)
	{
	case 'G':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'B':
      switch(seqnt)
	{
	case 'C':
	  return 1;
	  break;
	case 'G':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    case 'D':
      switch(seqnt)
	{
	  case 'A':
	    return 1;
	    break;
	  case 'G':
	    return 1;
	    break;
	  case 'T':
	    return 1;
	    break;
	  default:
	    return 0;
	    break;
	}
      break;
    case 'H':
      switch(seqnt)
	{
	case 'A':
	  return 1;
	  break;
	case 'C':
	  return 1;
	  break;
	case 'T':
	  return 1;
	  break;
	default:
	  return 0;
	  break;
	}
      break;
    default:
      fprintf(stderr,"[crisflash] ERROR: '%c' part of PAM sequence not a valid nucleotide\n", pamnt);
      exit(1);
    }
}

int baseMatchComplement(char seqnt, char pamnt) {
  /*
    Checks whether complement of a base in DNA sequence (seqnt) matches with base in PAM pattern (pamnt).
    Operates in two steps: 1) translates base to its complement, 2) matches complement with base in PAM using baseMatch() function.
  */
  switch(seqnt)
    {
    case 'N':
      return baseMatch('N',pamnt);
      break;
    case 'A':
      return baseMatch('T',pamnt);
      break;
    case 'a':
      return baseMatch('t',pamnt);
      break;
    case 'T':
      return baseMatch('A',pamnt);
      break;
    case 't':
      return baseMatch('a',pamnt);
      break;
    case 'G':
      return baseMatch('C',pamnt);
      break;
    case 'g':
      return baseMatch('c',pamnt);
      break;
    case 'C':
      return baseMatch('G',pamnt);
      break;
    case 'c':
      return baseMatch('g',pamnt);
      break;
    case 'n':
      return baseMatch('n',pamnt);
      break;      
    default:      
      fprintf(stderr,"[crisflash] ERROR: '%c' is not valid DNA base symbol!\n", seqnt);
      exit(1);
      return 0;
      break;
    }
}

char baseComplement(char seqnt) {  
  /* Returns complement of a base. */
  switch(seqnt)
    {
    case 'N':
      return 'N';
      break;
    case 'n':
      return 'n';
      break;   
    case 'A':
      return 'T';
      break;
    case 'a':
      return 't';
      break;     
    case 'T':
      return 'A';
      break;
    case 't':
      return 'a';
      break;      
    case 'G':
      return 'C';
      break;
    case 'g':
      return 'c';
      break;
    case 'C':
      return 'G';
      break;
    case 'c':
      return 'g';
      break;      
    default:
      fprintf(stderr,"[crisflash] ERROR: '%c' not valid DNA base symbol!\n", seqnt);
      exit(0);
      break;
    }
}

char baseComplementUppercaseOnly(char seqnt) {  
  /* Returns complement of a base. */
  switch(seqnt)
    {
    case 'A':
      return 'T';
      break;
    case 'T':
      return 'A';
      break;
    case 'G':
      return 'C';
      break;
    case 'C':
      return 'G';
      break;
    case '+':
      return '+';
      break;
    case '-':
      return '-';
      break;
    case '|':
      return '|';
      break;
    default:
      fprintf(stderr,"[crisflash] ERROR: '%c' not valid DNA base symbol!\n", seqnt);
      exit(0);
      break;
    }
}


char baseUpper(char seqnt) {  
  /* Returns uppercase of a DNA base. */
  switch(seqnt)
    {
    case 'N':
      return 'N';
      break;
    case 'n':
      return 'N';
      break;   
    case 'A':
      return 'A';
      break;
    case 'a':
      return 'A';
      break;     
    case 'T':
      return 'T';
      break;
    case 't':
      return 'T';
      break;      
    case 'G':
      return 'G';
      break;
    case 'g':
      return 'G';
      break;
    case 'C':
      return 'C';
      break;
    case 'c':
      return 'C';
      break;      
    default:
      fprintf(stderr,"[crisflash] ERROR: '%c' not valid DNA base symbol!\n", seqnt);
      exit(0);
      break;
    }
}

char isUppercase(char seqnt) {  
  switch(seqnt)
    {
    case 'A':
      return 1;
      break;
    case 'T':
      return 1;
      break;
    case 'G':
      return 1;
      break;
    case 'C':
      return 1;
      break;
    default:
      return 0;
      break;
    }
}

void updateSizeGRNAs(grna_list *g) {
  g->pos++;
  if(g->pos == g->length) {
    g->length = g->length+300000;
    g->starts = realloc(g->starts,g->length*sizeof(long long));
    g->strands = realloc(g->strands,g->length*sizeof(char));
  }
}

void freeGRNAs(grna_list *g) {
  free(g->starts);
  free(g->strands);
  free(g);
}

grna_list* installGRNAList(int guidelen) {
  grna_list *glist = malloc(sizeof(grna_list));
  glist->starts = malloc(sizeof(long long)*glist->length);
  glist->strands = malloc(sizeof(char)*glist->length);
  glist->mtypes = malloc(sizeof(char)*glist->length);
  glist->sequences = malloc(sizeof(char)*glist->length*guidelen);
  return glist;
}

void updateSizeGRNAList(grna_list *g, int guidelen) {
  g->pos++;
  if(g->pos == g->length) {
    g->length = g->length+300000;
    g->starts = realloc(g->starts,g->length*sizeof(long long));
    g->strands = realloc(g->strands,g->length*sizeof(char));    
    g->mtypes = realloc(g->mtypes,g->length*sizeof(char));
    g->sequences = realloc(g->sequences,g->length*sizeof(char)*guidelen);
  }
}

void freeGRNAList(grna_list *g, int guidelen) {
  free(g->sequences);
  free(g->mtypes);
  free(g->strands);
  free(g->starts);
  free(g);
}

void printGRNAs(grna_list *g, int guidelen, FILE *f) {
  /** Prints gRNAs in grna_list **/
  
  long long i=0;
  int j;
  long long k;
  char *grna = malloc(sizeof(char)*(guidelen+1));
  
  for(i=0;i<g->pos;i++) {
    if(g->strands[i] == '+') {       
      strncpy(grna,g->pstrand+(g->starts[i]),guidelen);
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<guidelen;j++) {
	grna[j]=g->nstrand[k-j];
      }
    }
    grna[guidelen] = '\0';
    // we add 1 to start position as for bed files we want chromosomes to start from position one. Internally all start
    // coordinates are computed starting from 0 for a reason that cas-offinder also does so
    fprintf(f,"%s\t%lld\t%lld\t%s\t0\t%c\n", g->chr, (g->starts[i]+1), (g->starts[i]+guidelen), grna,g->strands[i]);
  }
}

void printGRNAsHaplotypes(grna_list *g1, grna_list *g2, int guidelen, FILE *f) {
  /* Prints grnas for both haplotypes by sorting (to best of its abilities) chromosomal coordinates on fly */
  long long g1pos=0;
  long long g2pos=0;
  long long res;
  char *grna = malloc(sizeof(char)*(guidelen+1));
  char *grna2 = malloc(sizeof(char)*(guidelen+1));

  while( (g1pos < g1->pos) && (g2pos < g2->pos) ) {
    res = g1->starts[g1pos] - g2->starts[g2pos];
    if(res>0) {
      //process g2
      strncpy(grna,g2->sequences+(g2pos*guidelen),guidelen);
      grna[guidelen] = '\0';
      fprintf(f,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna, g2->mtypes[g2pos],g2->strands[g2pos]);
      g2pos++;
    }
    else if(res<0) {
      // process g1
      strncpy(grna,g1->sequences+(g1pos*guidelen),guidelen);
      grna[guidelen] = '\0';
      fprintf(f,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
      g1pos++;
    }
    else {
      // g1 and g2 positions must be equal
      strncpy(grna, g1->sequences+(g1pos*guidelen),guidelen);
      grna[guidelen] = '\0';

      strncpy(grna2, g2->sequences+(g2pos*guidelen),guidelen);
      grna2[guidelen] = '\0';
      // if sequences are equal as well, assume strand and mutation profiles are too. Report gRNA for both haplotypes.
      if(strcmp(grna, grna2) == 0) {
	fprintf(f,"%s\t%lld\t%lld\t%s:3%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
      }
      else {
	fprintf(f,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna, g1->mtypes[g1pos], g1->strands[g1pos]);
	fprintf(f,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna2, g2->mtypes[g2pos],g2->strands[g2pos]);
      }
      g1pos++;
      g2pos++;
    }

    if(g1pos == g1->pos) {
      // deal with rest of g2
      for(;g2pos < g2->pos;g2pos++) {
	strncpy(grna,g2->sequences+(g2pos*guidelen),guidelen);
	grna[guidelen] = '\0';
	fprintf(f,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna, g2->mtypes[g2pos],g2->strands[g2pos]);
      }
      break;
    }
    if(g2pos == g2->pos) {
      // deal with rest of g1
      for(;g1pos<g1->pos;g1pos++) {
	strncpy(grna,g1->sequences+(g1pos*guidelen),guidelen);
	grna[guidelen] = '\0';
	fprintf(f,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
      }
      break;
    }
  }
  free(grna);
  free(grna2);
}

void GRNAsToTrie(trie *T, grna_list *g, int guidelen, int chr_number) {
  /** Adds gRNAs in g to Trie tPrints gRNAs only for one chromosome at the time. **/

  // the loop structure for this function is identical to printGRNAs.
  long long i=0;
  long long k;
  int j;
  char *grna = malloc(sizeof(char)*(guidelen+1));
  int identical_seq = 0; // counts the sequences not added in the tree because they were already there
  
  for(i=0;i<g->pos;i++) {
    if(g->strands[i] == '+') {
      // '0' and '3' stand for 'no variant modification in guide' and for 'guide sequences in both haplotypes being the same', respectively.
      if(!TrieAdd(T, g->pstrand+(g->starts[i]), guidelen, chr_number, g->strands[i], '0', '3',g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
	exit(1);
      }
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<guidelen;j++) {
	grna[j]=g->nstrand[k-j];
      }
      if(!TrieAdd(T, grna, guidelen, chr_number, g->strands[i], '0', '3',g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence in chr=%s, pos=%lld, strand=%c). Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
	exit(1);
      }
    }
  }
  free(grna);
}

void GRNAsToTrieHaplotypes(trie *T, grna_list *g1, grna_list *g2, int guidelen, int chr_number) {
  /** Adds gRNAs in g1 (gRNAs from haplotype 1) and g2 (gRNAs from haplotype 2) from specified sequence/chromosome to Trie. **/
  // the function follows same structure/layout as printGRNAsHaplotypes

  long long g1pos=0;
  long long g2pos=0;
  long long res;
  char *grna = malloc(sizeof(char)*(guidelen+1));
  char *grna2 = malloc(sizeof(char)*(guidelen+1));
  int identical_seq = 0; // counts the sequences not added in the tree because they were already there

  while( (g1pos < g1->pos) && (g2pos < g2->pos) ) {
    res = g1->starts[g1pos] - g2->starts[g2pos];
    if(res>0) {
      //process g2
      // strncpy(grna,g2->sequences+(g2pos*guidelen),guidelen);
      // grna[guidelen] = '\0';
      // fprintf(stdout,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna, g2->mtypes[g2pos],g2->strands[g2pos]);
      if(!TrieAdd(T, g2->sequences+(g2pos*guidelen), guidelen, chr_number, g2->strands[g2pos], g2->mtypes[g2pos], '2', g2->starts[g2pos], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g2->starts[g2pos], g2->strands[g2pos]);
	exit(1);
      }
      g2pos++;
    }
    else if(res<0) {
      // process g1
      // strncpy(grna,g1->sequences+(g1pos*guidelen),guidelen);
      // grna[guidelen] = '\0';
      // fprintf(stdout,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
      if(!TrieAdd(T, g1->sequences+(g1pos*guidelen), guidelen, chr_number, g1->strands[g1pos], g1->mtypes[g1pos], '1', g1->starts[g1pos], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g1->starts[g1pos], g1->strands[g1pos]);
	exit(1);
      }
      g1pos++;
    }
    else {
      // check g1 and g2 positions must be equal
      strncpy(grna, g1->sequences+(g1pos*guidelen),guidelen);
      grna[guidelen] = '\0';
      strncpy(grna2, g2->sequences+(g2pos*guidelen),guidelen);
      grna2[guidelen] = '\0';
      // if sequences are equal, enter only one guide, assume strand and mutation profiles are too. Report gRNA for both haplotypes.
      if(strcmp(grna, grna2) == 0) {
	// fprintf(stdout,"%s\t%lld\t%lld\t%s:3%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
	if(!TrieAdd(T, g1->sequences+(g1pos*guidelen), guidelen, chr_number, g1->strands[g1pos], g1->mtypes[g1pos], '3',g1->starts[g1pos], 0, &identical_seq)) {
	  fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g1->starts[g1pos], g1->strands[g1pos]);
	  exit(1);	
	}
      }
      // g2 and g2 start positions were not equal!
      else {
	// fprintf(stdout,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna, g1->mtypes[g1pos], g1->strands[g1pos]);
	if(!TrieAdd(T, g1->sequences+(g1pos*guidelen), guidelen, chr_number, g1->strands[g1pos], g1->mtypes[g1pos], '1',g1->starts[g1pos], 0, &identical_seq)) {
	  fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g1->starts[g1pos], g1->strands[g1pos]);
	  exit(1);
	}
	// fprintf(stdout,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna2, g2->mtypes[g2pos],g2->strands[g2pos]);
	if(!TrieAdd(T, g2->sequences+(g2pos*guidelen), guidelen, chr_number, g2->strands[g2pos], g2->mtypes[g2pos], '2', g2->starts[g2pos], 0, &identical_seq)) {
	  fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g2->starts[g2pos], g2->strands[g2pos]);
	  exit(1);
	}
      }
      
      g1pos++;
      g2pos++;
    }
    
    // check for end in g1 and g2 gRNA lists
    if(g1pos == g1->pos) {
      // deal with rest of g2
      for(;g2pos < g2->pos;g2pos++) {
	// strncpy(grna,g2->sequences+(g2pos*guidelen),guidelen);
	// grna[guidelen] = '\0';
	// fprintf(stdout,"%s\t%lld\t%lld\t%s:2%c\t0\t%c\n", g2->chr, (g2->starts[g2pos]+1), (g2->starts[g2pos]+guidelen), grna, g2->mtypes[g2pos],g2->strands[g2pos]);
	if(!TrieAdd(T, g2->sequences+(g2pos*guidelen), guidelen, chr_number, g2->strands[g2pos], g2->mtypes[g2pos], '2', g2->starts[g2pos], 0, &identical_seq)) {
          fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g2->starts[g2pos], g2->strands[g2pos]);
          exit(1);
        }
      }
      break;
    }
    if(g2pos == g2->pos) {
      // deal with rest of g1
      for(;g1pos<g1->pos;g1pos++) {
	// strncpy(grna,g1->sequences+(g1pos*guidelen),guidelen);
	// grna[guidelen] = '\0';
	// fprintf(stdout,"%s\t%lld\t%lld\t%s:1%c\t0\t%c\n", g1->chr, (g1->starts[g1pos]+1), (g1->starts[g1pos]+guidelen), grna,g1->mtypes[g1pos],g1->strands[g1pos]);
	if(!TrieAdd(T, g1->sequences+(g1pos*guidelen), guidelen, chr_number, g1->strands[g1pos], g1->mtypes[g1pos], '1', g1->starts[g1pos], 0, &identical_seq)) {
          fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g1->starts[g1pos], g1->strands[g1pos]);
          exit(1);
        }
      }
      break;
    }
  }
  free(grna);
  free(grna2);

}


grna_list *fastaSequenceToGRNAs(char *header, char *s, char *sr, long long spos, char *pam) {

  grna_list *glist = malloc(sizeof(grna_list));
  glist->pstrand = s;
  glist->nstrand = sr;
  glist->chr = header;
  glist->length = 300000;
  glist->pos = 0;
  glist->starts = malloc(sizeof(long long)*glist->length);
  glist->strands = malloc(sizeof(char)*glist->length);
  int nfree;  
  
  // variables for calculating PAMs
  int pamlen = strlen(pam); // length of the PAM motif
  int guidelen = pamlen+PROTOSPACER_LENGTH; // length of the gRNA
  int y, nt, endstart;    
  
  // For first bases up to the length of protospacer, search for PAMs on negative strand only
  for(nt=0; nt<PROTOSPACER_LENGTH; nt++)
    {
      // Following line is a safetynet for avoiding reading over 3' end of the sequence
      if ((nt+guidelen) > spos) { break; }
      y = 0;
      while(baseMatch(sr[nt+y], pam[pamlen-1-y]))
	{	  
	  y++;
	  if (y==pamlen) {
	    // add gRNA to glist only if the sequence is free from 'N'.
	    nfree = 1;
	    for(int j=0;j<guidelen;j++) {
	      if (sr[nt+guidelen-j-1] == 'N') { nfree = 0; break; }
	    }
	    if (nfree) {
	      glist->starts[glist->pos] = nt;
	      glist->strands[glist->pos] = '-';
	      updateSizeGRNAs(glist);
	    }
	    break;
	  }
	}  
    }
  
  // Search for PAM sites on both strands
  for(nt=PROTOSPACER_LENGTH; nt<(spos-guidelen); nt++)
    {
      y = 0;
      // Search PAM sequences on positive strand.		
      while(baseMatch(s[nt+y], pam[y]))
	{
	  y++;	  
	  if (y==pamlen) {
	    // add gRNA to glist only if the sequence is free from 'N'.
	    nfree = 1;
	    for(int j=0;j<(guidelen);j++) {
	      if(s[nt-PROTOSPACER_LENGTH+j] == 'N') { nfree = 0; break; }
	    }
	    if (nfree) {
	      glist->strands[glist->pos] = '+';
	      glist->starts[glist->pos] = nt-PROTOSPACER_LENGTH;
	      updateSizeGRNAs(glist);
	    }
	    break;
	  }
	} // end of while
      
      // Search PAM on negative strand
      y = 0;
      while(baseMatch(sr[nt+y], pam[pamlen-1-y]))
	{
	  y++;
	  if (y==pamlen) {
	    // add gRNA to glist only if the sequence is free from 'N'.
	    nfree = 1;
	    for(int j=0;j<(guidelen);j++) {
	      if (sr[nt+guidelen-j-1] == 'N') { nfree =0; break; }
	    }
	    if (nfree) {	      
	      glist->strands[glist->pos] = '-';
	      glist->starts[glist->pos] = nt;
	      updateSizeGRNAs(glist);
	    }
	    break;
	  }
	} // end of while	
    }
  
  // Search for PAMs only on positive strand for the remaining bit of the sequence
  if((spos-guidelen-PROTOSPACER_LENGTH)<0) {
    endstart = PROTOSPACER_LENGTH;
  }
  else {
    endstart = (spos-guidelen);
  }
  for(nt=endstart; nt<(spos-pamlen+1); nt++)
    {
      // safetynet of exceptionally shot DNA sequence
      if ((nt-PROTOSPACER_LENGTH) < 0) { continue; }			       
      y = 0;
      while(baseMatch(s[nt+y], pam[y]))
	{
	  y++;
	  if (y==pamlen) {
	    // add gRNA to glist only if the sequence is free from 'N'.
	    nfree = 1;
	    for(int j=0;j<(guidelen);j++) {
	      if(s[nt-PROTOSPACER_LENGTH+j] == 'N') { nfree = 0; break; }
	    }
	    if(nfree) {
	      glist->starts[glist->pos] = nt-PROTOSPACER_LENGTH;
	      glist->strands[glist->pos] = '+';
	      updateSizeGRNAs(glist);
	    }
	    break;
	  }
	}
    }
  fprintf(stderr,"[crisflash] %lld gRNA sequences discovered in '%s'.\n",glist->pos, header);
  return glist;
}

long long isPamStartOnPositiveStrand(grna_list *glist, long long nt, char *pam, int pamlen) {
  /* Matches pam sequence of lenght pamlen in sequence s of length spos, starting from position nt. Returns the end position of the pam sequence in s or -1. */
  int iS = 0;
  int iP = 0;
  int res;
 
  while( (res=baseMatch(glist->pstrand[nt+iS], pam[iP]) )) {
    // fprintf(stdout,"Seq[%lld]=%c\tPam[%d]=%c\tiS=%d\tiP=%d\n",nt+iS,s[nt+iS], iP, pam[iP], iS, iP);   
    if(res == 0) {
      return -1;
    }
    else {
      if (res == 1) {
	iP++;
	if (iP == pamlen) {
	  return nt+iS;
	}
      }
      iS++;
      if( (nt+iS) == glist->slen) {
	return -1;
      }
    }
  }
  return -1;
}

long long isPamStartOnNegativeStrand(grna_list *glist, long long nt, char *pam, int pamlen) {
  /* Matches pam sequence of lenght pamlen in reverse in sequence sr of length spos, starting from position nt. Returns the end position of the pam sequence in s or -1. */
  int iS = 0;
  int iP = 0;
  int res;

  while( (res=baseMatch(glist->nstrand[nt+iS], pam[pamlen-1-iP]) )) {
    if(res == 0) {
      return -1;
    }
    else {
      if (res == 1) {
	iP++;
	if (iP == pamlen) {
	  return nt;
	}
      }
      iS++;
      if( (nt+iS) == glist->slen) {
	return -1;
      }
    }
  }
  return -1;
}

int addGrnaOnPositiveStrand(grna_list *glist, long long pamendpos, int guidelen, int pamlen, long long real_pam_start) {
  /** Adds gRNA sequence to the list gRNA sequences in glist->sequences unless the protospacer hangs over an edge of the target (chromosome) sequence. Returns 1 for success and 0 for failure. **/
  /** real_pam_start is a computed hypothetical location of where the PAM would have started if the reference genome would not have been edited with variant and indel info. **/
  long long i = pamendpos;
  int gpos = 0; 
  char mtype = '0';
  int mod;

  while(i > -1) {
    if(glist->pstrand[i] == 'N') {
      // fprintf(stdout,"Sequence not added! %d bases available before running into N\n", gpos);  
      return 0;
    }
    if((mod = baseInOther(glist->pstrand[i])) != 5) {
      i--;      
      if(gpos < pamlen) {
	mtype = 'b';
      }
      else {
	if (mtype == 'b') { mtype = 'f'; }
	else { mtype = 'c'; }
      }
      continue;
    }
    glist->sequences[(glist->pos*guidelen)+(guidelen-gpos-1)] = glist->pstrand[i];
    i--;
    gpos++;
    if(gpos == guidelen) {
      // gRNA start here is reported as location where the gRNA would have started in standard reference genome. All effort has been made to make sure the real_pam_start value would be correct ...
      glist->starts[glist->pos] = real_pam_start - PROTOSPACER_LENGTH;
      glist->strands[glist->pos] = '+';
      glist->mtypes[glist->pos] = mtype;
      updateSizeGRNAList(glist,guidelen);
      return 1;
    }
  }
  // fprintf(stderr,"Sequence not added! %d bases available before running over the edge of the target/chromosome sequence\n", gpos);  
  return 0;
}

int addGrnaOnNegativeStrand(grna_list *glist, long long pamendpos, int guidelen, int pamlen, long long real_pam_start) {
  /** Adds gRNA sequence to the list gRNA sequences in glist->sequences unless the protospacer hangs over an edge of the target (chromosome) sequence. Returns 1 for success and 0 for failure. **/
  /** real_pam_start is a computed hypothetical location of where the PAM would have started if the reference genome would not have been edited with variant and indel info. **/
  long long i = pamendpos;
  int gpos = 0; 
  char mtype = '0';
  int mod;

  while(i < glist->slen) {
    if(glist->nstrand[i] == 'N') {
      // fprintf(stdout,"Sequence not added! %d bases available before running into N\n", gpos);  
      return 0;
    }
    if((mod = baseInOther(glist->nstrand[i])) != 5) {
      i++;
      if(gpos < pamlen) {
	mtype = 'b';
      }
      else {
	if (mtype == 'b') { mtype = 'f'; }
	else { mtype = 'c'; }
      }
      continue;
    }
    glist->sequences[(glist->pos*guidelen)+(guidelen-gpos-1)] = glist->nstrand[i];
    i++;
    gpos++;
    if(gpos == guidelen) {
      // gRNA start here is reported as location where the gRNA would have started in standard reference genome. All effort has been made to make sure the real_pam_start value would be correct ..
      glist->starts[glist->pos] = real_pam_start;
      glist->strands[glist->pos] = '-';
      glist->mtypes[glist->pos] = mtype;
      updateSizeGRNAList(glist,guidelen);
      return 1;
    }
  }
  // fprintf(stderr,"Sequence not added! %d bases available before running over the edge of the target/chromosome sequence\n", gpos);  
  return 0;
}

grna_list *fastaVariantSequenceToGRNAsequences(char *header, char *s, char *sr, long long spos, char *pam) {
  /* Identifies and returns a list of gRNAs from sequence edited for variants and indels. */

  // variables for calculating PAMs
  int pamlen = strlen(pam); // length of the PAM motif
  int guidelen = pamlen+PROTOSPACER_LENGTH; // length of the gRNA
  // grna_list *glist = installGRNAlist(guidelen);
  grna_list *glist = malloc(sizeof(grna_list));
  glist->length = 300000;
  glist->starts = malloc(sizeof(long long)*glist->length);
  glist->strands = malloc(sizeof(char)*glist->length);
  glist->mtypes = malloc(sizeof(char)*glist->length);
  glist->sequences = malloc(sizeof(char)*glist->length*guidelen);
  glist->pstrand = s;
  glist->nstrand = sr;
  glist->slen = spos;
  glist->chr = header;
  glist->pos = 0;
 
  long long nt;
  long long pamendpos;
  long long realpos = 0; // note that chromosomes should normally start on position 1. Here we start from 0 as does cas-offinder.
  int mod;

  for(nt=0; nt < spos; nt++) {
    // fprintf(stdout,"s[%lld]=%c\n",nt,s[nt]);
    if( (mod = baseInOther(s[nt])) != 5) {
      // fprintf(stdout,"realpos modified from %lld to %lld\n",realpos, realpos+mod);
      realpos = realpos + mod;
      continue;
    }    
    if ((pamendpos = isPamStartOnPositiveStrand(glist, nt, pam, pamlen)) != -1) {
      addGrnaOnPositiveStrand(glist,pamendpos, guidelen, pamlen, realpos);
      // fprintf(stdout,"+PAM ends at %lld, glist->pos=%lld\n",pamendpos, glist->pos);
    }
    if ((pamendpos = isPamStartOnNegativeStrand(glist, nt, pam, pamlen)) != -1) {
      addGrnaOnNegativeStrand(glist,pamendpos, guidelen, pamlen, realpos);
      // fprintf(stdout,"-PAM ends at %lld, glist->pos=%lld\n",pamendpos, glist->pos);
    }
    realpos++;
  }
  fprintf(stdout,"[crisflash] Estimated reference length for '%s' excluding indels: %lld\n",header,realpos);
  return glist;
}

faread_struct* installFastaReader(char *genomefname, int no_low_complexity) {
  /* Installs faread_struct for provided file name. */
  faread_struct* fas = malloc(sizeof(faread_struct));
  
  fas->no_low_complexity = no_low_complexity;
  fas->limit = 300000000; // set maximul length for s. Note that this will be adjusted as needed later in the reading process.
  fas->s = malloc(sizeof(char)*fas->limit); // chromosome sequence + strand
  fas->sr = NULL;
  fas->slen = 0;
  fas->buffer = malloc(sizeof(char)*BUFFER_SIZE+1);
  fas->i = 0;
  fas->header = malloc(sizeof(char)*BUFFER_SIZE);
  fas->header[0] = '\0';
  fas->fasta_header = 0; // boolean value to show whether char read is from fasta header or from sequence: 0 - sequence, 1 - header
  fas->headerlen = 0;
  fas->nrheaders = 0;
  
  fas->alive = 1;
  fas->sequence = 0;
  
  // open file and read first chunk (first fasta entry, e.g. chromosome)
  fas->fd = open(genomefname, O_RDONLY);
  fas->bytes_read = read(fas->fd, fas->buffer, BUFFER_SIZE);

  return fas;
}

void fastaReaderImproveSequence(faread_struct *fas) {
  /*
    Pre-processing of fasta sequence before calling gRNAs.
    In case no_low_complexity == 1, converts all low complexity sequence to 'N'.
    Converts all bases to uppercase.
    Creates sequence for negative stand in fas->sr;
  */    
  int nt;
  if(fas->sr != NULL) {
    free(fas->sr);
  }
  fas->sr=malloc(sizeof(char)*fas->slen);  
  if(fas->no_low_complexity) {
    // we will convert all soft masked bases (lower case) to N and others to uppercase.
    for(nt=0; nt<fas->slen; nt++) {
      switch(fas->s[nt]) {
      case 'A':
	fas->sr[nt] = 'T';
	break;
      case 'T':
	fas->sr[nt] = 'A';
	break;
      case 'G':
	fas->sr[nt] = 'C';
	break;
      case 'C':
	fas->sr[nt] = 'G';
	break;
      case '+':
	fas->sr[nt] = '+';
      case '-':
	fas->sr[nt] = '-';
      case '|':
	fas->sr[nt] = '|';
      default:
	fas->s[nt] = 'N';
	fas->sr[nt] = 'N';
      }
    }
  }
  else {
    // convert all bases to uppercase
    for(nt=0; nt<fas->slen; nt++) {
      switch(fas->s[nt]) {
      case 'a':
	fas->s[nt] = 'A';
	fas->sr[nt] = 'T';
	break;
      case 't':
	fas->s[nt] = 'T';
	fas->sr[nt] = 'A';
	break;
      case 'g':
	fas->s[nt] = 'G';
	fas->sr[nt] = 'C';
	break;
      case 'c':
	fas->s[nt] = 'C';
	fas->sr[nt] = 'G';
	break;
      case 'n':
	fas->s[nt] = 'N';
	fas->sr[nt] = 'N';
	break;
      case 'N':
	fas->sr[nt] = 'N';
	break;
      default:
	fas->sr[nt] = baseComplementUppercaseOnly(fas->s[nt]);
	break;
      }
    }
  }
}

void freeFastaReader(faread_struct *fas) {
  /* Frees allocations within faread_struct. */
  free(fas->s);
  if(fas->sr) { free(fas->sr); }
  free(fas->buffer);
  free(fas->header);  
  free(fas);
}

int fastaReader(faread_struct *fas) {
  /*
    Reads fasta file one sequence at the time. Returns 1 if sequence was read successfully and 0 otherwise.
    After each reading round:
    fas->header: contains header of the fasta sequence
    fas->s: contains sequence from + strand
    fas->slen: contains the sequence length
  */
  if(!fas->alive) {
    return 0;
  }
  
  // continue reading from the position
  while(!fas->sequence) {
    while(fas->i < fas->bytes_read) {      
      
      // check for and deal with fasta header lines
      if (fas->buffer[fas->i] == '>') {
	fas->fasta_header = 1;
	fas->i++;
	if(fas->nrheaders > 0) {
	  // identify all gRNAs in the chromosome by PAM sequences. Make sure one does not run over the edge of the chromosome in either end.
	  fas->sequence = 1;
	  break;
	}
	continue;
      }
      if (fas->fasta_header) {
	if ((fas->headerlen+1)==BUFFER_SIZE) {
	  fprintf(stderr,"[crisflash] ERROR: Fasta header too long! Advice: edit BUFFER_SIZE and re-compile! Exiting!\n");
	  exit(1);
	}
	if (fas->buffer[fas->i] == '\n') {
	  fas->header[fas->headerlen] = '\0';
	  fas->fasta_header = 0;
	  fas->headerlen = 0;
	  fas->slen = 0;	    
	  fas->nrheaders++;
	}
	else {
	  fas->header[fas->headerlen] = fas->buffer[fas->i];
	  fas->headerlen++;	  
	}
	fas->i++;
	continue;
      }
      
      // ignoe rewlines while reading chromosome sequence.
      if(fas->buffer[fas->i] == '\n') {
	fas->i++;
	continue;
      }
      
      // add base to chromosome sequence
      fas->s[fas->slen] = fas->buffer[fas->i];
      fas->slen++;
      
     // update chromosome length
      if (fas->slen == fas->limit) {
	fas->limit = fas->limit+300000000;
	fas->s = (char*) realloc(fas->s, fas->limit);
      }
      fas->i++;
    }
    // if end of buffer, read more from file handler
    if (fas->i == fas->bytes_read) {
      fas->bytes_read = read(fas->fd, fas->buffer, BUFFER_SIZE);
      fas->i = 0;
      // if no more bytes were read, shut the reading down
      if(fas->bytes_read == 0) {
	fas->alive = 0;
	fas->sequence = 1;
      }
    }
  }
  fas->sequence = 0;
  if(fas->header[0] == '\0') {
    fprintf(stderr,"[crisflash] ERROR: No fasta header. Exiting!\n");    
    exit(1);
  }
  
  return 1;
}

void readFaToTrie(trie *T, char *genomefname, char *pam, char* outputfname, int uppercaseOnly, int printGRNAsOnly) {
  /* Identifies all gRNAs in reference sequence and loads them to Trie. */
  faread_struct *fas = installFastaReader(genomefname, uppercaseOnly);
  grna_list *g = NULL;
  int chr_number;
  int guidelen = PROTOSPACER_LENGTH+strlen(pam);
  FILE* outfh;
  
  if(printGRNAsOnly) {
    outfh = open_file(outputfname, "w");
  }

  while(fastaReader(fas)) {
    fprintf(stdout,"[crisflash] Processing %s (length %lld) in %s ...\n", fas->header, fas->slen, genomefname);
    fastaReaderImproveSequence(fas);   
    // faSequenceToTrie(T, fas->header, fas->s, fas->sr, fas->slen, pam, outputfname, uppercaseOnly);
    g = fastaSequenceToGRNAs(fas->header, fas->s, fas->sr, fas->slen, pam);
    if(printGRNAsOnly) {
      printGRNAs(g, guidelen, outfh);
    }
    else {
      chr_number = addChr(T, g->chr, strlen(g->chr)); // we get the chr number in order of appearence in fasta file
      GRNAsToTrie(T, g, guidelen,chr_number);
    }
    freeGRNAs(g);
  }
  freeFastaReader(fas);
}

void readFaToTrieVCF(trie *T, char *genomefname1, char *genomefname2, char *pam, char* outputfname, int uppercaseOnly, int printGRNAsOnly) {
  /** Reads haplotype level genome references, identifies gRNAs and loads them to trie. If gRNA is identical in both haplotypes, only one of them is loaded to Trie and labeled accordinly. **/
  faread_struct *fas1 = installFastaReader(genomefname1, uppercaseOnly);
  faread_struct *fas2 = installFastaReader(genomefname2, uppercaseOnly);
  grna_list *g1 = NULL;
  grna_list *g2 = NULL;
  int chr_number;
  int guidelen = PROTOSPACER_LENGTH+strlen(pam);
  FILE* outfh;
  int read_success1, read_success2;
  
  if(printGRNAsOnly) {
    outfh = open_file(outputfname, "w");
  }

  fprintf(stdout,"[crisflash] Starting gRNA discovery.\n");
  // read first fasta sequence (chromosome) in file
  fprintf(stdout,"[crisflash] Reading sequence from %s.\n",genomefname1);
  fflush(stdout);
  read_success1 = fastaReader(fas1);
  fprintf(stdout,"[crisflash] Reading sequence from %s.\n",genomefname1);
  fflush(stdout);
  read_success2 = fastaReader(fas2);
  
  while(read_success1 && read_success2) {

    if (strcmp(fas1->header, fas2->header) != 0) {
      fprintf(stdout,"[crisflash] ERROR! Chr order in %s and %s not the same: %s != %s. Exiting!\n", genomefname1, genomefname2, fas1->header, fas2->header);
      exit(1);
    }
    
    fprintf(stdout,"[crisflash] Processing 1st haploid sequence for '%s;.\n", fas1->header);
    fflush(stdout);
    fastaReaderImproveSequence(fas1);
    fprintf(stdout,"[crisflash] Processing 2nd haploid sequence for '%s'.\n", fas1->header);
    fflush(stdout);
    fastaReaderImproveSequence(fas2);

    fprintf(stdout,"[crisflash] Running gRNA discovery for 1st haploid sequence of '%s'\n",fas1->header);
    fflush(stdout);
    g1 = fastaVariantSequenceToGRNAsequences(fas1->header, fas1->s, fas1->sr, fas1->slen, pam);
    fprintf(stdout,"[crisflash] Running gRNA discovery for 2nd haploid sequence of '%s'\n",fas1->header);
    fflush(stdout);
    g2 = fastaVariantSequenceToGRNAsequences(fas2->header, fas2->s, fas2->sr, fas2->slen, pam);

    if(printGRNAsOnly) {
      fprintf(stdout,"[crisflash] Saving gRNAs in '%s'.\n", fas1->header);
      fflush(stdout);
      printGRNAsHaplotypes(g1, g2, guidelen, outfh);
    }
    else {
      fprintf(stdout,"[crisflash] Indexing gRNAs in %s.\n", fas1->header);    
      fflush(stdout);
      chr_number = addChr(T, g1->chr, strlen(g1->chr)); // we get the chr number in order of appearence in fasta file
      chr_number = addChr(T, g2->chr, strlen(g1->chr)); // we get the chr number in order of appearence in fasta file
      GRNAsToTrieHaplotypes(T, g1, g2, guidelen,chr_number);
    }
    freeGRNAList(g1,guidelen);
    freeGRNAList(g2,guidelen);

    // read next fasta sequence (chromosome)
    fprintf(stdout,"[crisflash] Reading sequence from %s.\n",genomefname1);
    fflush(stdout);
    read_success1 = fastaReader(fas1);
    fprintf(stdout,"[crisflash] Reading sequence from %s.\n",genomefname1);
    fflush(stdout);
    read_success2 = fastaReader(fas2);
  }
  freeFastaReader(fas1);
  freeFastaReader(fas2);
}

double now()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + tv.tv_usec / 1000000.;
}

void writeInBed_lighter(FILE* f, mcontainer *m, char* chrom, long long start, long long end, char* gRNA, char strand, trie* T)
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
  // safety net in case candidate gRNA had 0 exact matches. No exact matches is currently screwing up the scoring formula and the score
  // returned may be out of 0 to 1 bounds. Henece replacing such score with 0.
  if(nr_of_exact_matches == 0) {
    m->score = 0;
  }
  fprintf(f, "%s\t%lld\t%lld\t%s/%d/%d\t%f\t%c\n", chrom, start, end, gRNA, nr_of_exact_matches, nr_of_off_targets, m->score, strand);
}

void writeInCasOffinder(FILE* fh, char *gRNA, mcontainer *m, trie* T, int addVariantHaplotypeInfo)
{
  /* Writes match information in Cas-OFFinder format. If addVariantInfo !=0 adds extra column containing variant info on haplotype level.  */
  
  // if mcontainer is empty, i.e. no matches, return right away.
  if(m->pos == 0) {
    return;
  }
  
  // Columns in the output: sequence+NNN   chr   pos   sequence_with_mismatched_bases_lowercase   strand   nr_of_mismatches
  int i = 0;
  // T->readlen is the same as guidelen or PROTOSPACER_LENGTH+pamlen
  char *gRNA4print = (char *) malloc(sizeof(char)*(T->readlen+1));
  while(i < T->readlen) {
    if(i>=PROTOSPACER_LENGTH) {
      gRNA4print[i] = 'N';
    }
    else {
      gRNA4print[i] = gRNA[i];
    }
    i++;
  }
  gRNA4print[T->readlen] = '\0';

  i=0;  
  char *seq = (char *) malloc(sizeof(char)*(strlen(m->sarray[i])+1));  

  while(i < m->pos)
    {
      strcpy(seq,m->sarray[i]);
      for(char* c=seq; (*c=toupper(*c)); ++c) ;      
      int j = 0;
      while (j < m->marray[i]-1)
	{
	  seq[m->parray[i][j]-1] = tolower(seq[m->parray[i][j]-1]);
	  j++;
	}
      if (m->marray[i] > 0)
	{
	  seq[m->parray[i][j]-1] = tolower(seq[m->parray[i][j]-1]);
	}
      if (addVariantHaplotypeInfo) {
	for(j=0; j < m->nodes[i]->hits;j++)
	  {
	    fprintf(fh,"%s\t%s\t%d\t%s\t%c\t%d\t%c%c\n", gRNA4print, T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j], seq, m->nodes[i]->strands[j], m->marray[i], m->nodes[i]->htypes[i],m->nodes[i]->mtypes[i]);
	  }
      }
      else {
	for(j=0; j < m->nodes[i]->hits;j++)
	  {
	    fprintf(fh,"%s\t%s\t%d\t%s\t%c\t%d\n", gRNA4print, T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j], seq, m->nodes[i]->strands[j], m->marray[i]);
	  }
      }
      i++;
    }
  free(gRNA4print);
  free(seq);
}

void *thread_worker(void *arg)
{
  // cast the argument into its real nature
  struct arg_struct *args_table = (struct arg_struct *)arg;

  // match
  mcontainer* m = TrieAMatch(args_table->T, args_table->grna, args_table->T->readlen, args_table->maxMismatch);
  if(args_table->outFileType == 2) {
    writeInCasOffinder(args_table->outfh, args_table->grna, m, args_table->T,0);
  }
  else if(args_table->outFileType == 3) {
    writeInCasOffinder(args_table->outfh, args_table->grna, m, args_table->T, 1);
  }
  else {
    writeInBed_lighter(args_table->outfh, m, args_table->chr, args_table->start, args_table->end, args_table->grna, args_table->strand, args_table->T);
  }
  // mcontainer_print(m);
  mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
  // label thread being free/available
  args_table->free = 1;

  pthread_exit(NULL);
}

void GRNAsMatch(trie *T, grna_list *g, int guidelen, char *pam, int maxMismatch, FILE *outfh, int outFileType) {
  /** Matches each gRNA in g to trie T. Most of the code strucure is more or less copy/paste from GRNAsToTrie and printGRNAs **/
  
  // the loop structure for this function is identical to printGRNAs.
  long long i=0;
  long long k;
  int j;
  char *grna = malloc(sizeof(char)*(guidelen+1)); // here the gRNA sequence will be gRNA protospacer + PAM, not the sequence as in candidate area!
  mcontainer* m;
  int pamlen = strlen(pam);

  for(i=0;i<g->pos;i++) {

    // create gRNA sequence consisting protospacer sequence + PAM sequence.
    if(g->strands[i] == '+') {
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	grna[j]=g->pstrand[g->starts[i]+j];
      }
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	grna[j]=g->nstrand[k-j];
      }
    }
    // add PAM sequence
    for(j=0;j<pamlen;j++) {
      grna[PROTOSPACER_LENGTH+j]=pam[j];
    }
    grna[guidelen] = '\0';

    fprintf(stdout,"%s\t%lld\t%lld\t%s\t0\t%c\n", g->chr, g->starts[i], (g->starts[i]+guidelen), grna, g->strands[i]);
    
    // match grna
    m = TrieAMatch(T, grna, T->readlen, maxMismatch);
    if(outFileType==2) {
      writeInCasOffinder(outfh, grna, m, T, 0);
    }
    else if(outFileType == 3) {
      writeInCasOffinder(outfh, grna, m, T, 1);
    }
    else {
      // i.e. outFileType==1; and default behaviour
      writeInBed_lighter(outfh, m, g->chr, g->starts[i], g->starts[i]+guidelen, grna, g->strands[i], T);
    }

    mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
  }
  free(grna);
}

void GRNAsMatchThreaded(trie *T, grna_list *g, int guidelen, char *pam, int maxMismatch, FILE *outfh, int outFileType, int nr_of_threads) {
  /** Matches each gRNA in g to trie T. The same code as in GRNAsMatch function except for additional code blocks to support multithreading **/
    
  // the loop structure for this function is identical to printGRNAs.
  long long i=0;
  long long k;
  int j;
  char *grna = NULL;
  grna = malloc(sizeof(char)*(guidelen+1)); // here the gRNA sequence will be gRNA protospacer + PAM, not the sequence as in candidate area!
  grna = calloc((guidelen+1),sizeof(char)); // here the gRNA sequence will be gRNA protospacer + PAM, not the sequence as in candidate area!
  int pamlen = strlen(pam);

  // Install threading related variables and attributes
  pthread_t threads[nr_of_threads]; // list of threads
  pthread_attr_t attr; // thread attribute necessary for thread tracking on joining.
  // initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  void *status;
  int res; // to catch thread exit result 
  struct arg_struct args_table[nr_of_threads];
  int tpos; // thread position, expected to be 0 <= tpos < nr_of_threads
  // install varables in args table
  
  for(tpos=0;tpos<nr_of_threads;tpos++) {  
    args_table[tpos].T = T;
    args_table[tpos].i = tpos; // thread number
    args_table[tpos].maxMismatch = maxMismatch;
    args_table[tpos].free = 2; // 0 - thread in use, 1 - thread was in use and has become free, 2 - thread has is free and has not yet been used.
    args_table[tpos].outfh = outfh;
    args_table[tpos].outFileType = outFileType;
    args_table[tpos].grna = malloc(sizeof(char)*(guidelen+1));
  }
  tpos = 0;
  
  // we will process one gRNA at the time
  for(i=0;i<g->pos;i++) {    
    // create gRNA sequence consisting protospacer sequence + PAM sequence.
    if(g->strands[i] == '+') {
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	grna[j]=g->pstrand[g->starts[i]+j];
      }
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	grna[j]=g->nstrand[k-j];
      }
    }
    // add PAM sequence
    for(j=0;j<pamlen;j++) {
      grna[PROTOSPACER_LENGTH+j]=pam[j];
    }
    
    // find next free thread for matching the found gRNA
    res = 0;
    while(!res) {
      if(tpos == nr_of_threads) { tpos = 0; }
      if(args_table[tpos].free > 0) {
	// thread has been in use in past has completed, try to join
	if(args_table[tpos].free == 1) {
	  if (pthread_join(threads[tpos], &status)) {
	    perror("pthread_join");
	    exit(1);
	  }
	}
	// load_args_table with relevant values
	args_table[tpos].free = 0;
	strcpy(args_table[tpos].grna, grna);
	args_table[tpos].chr = g->chr;
	args_table[tpos].start = g->starts[i];
	args_table[tpos].end = g->starts[i] + T->readlen;
	args_table[tpos].strand = g->strands[i];	
	if(pthread_create(&threads[tpos], &attr, thread_worker, &args_table[tpos])) {
	  perror("pthread_create");
	  exit(1);
	}
	res = 1;
      }
      tpos++;
    } // end of while trying to find free thread to match gRNA candidate
  } // end of for loop for processing all gRNAs

  // join the remaining threads
  for(tpos=0;tpos<nr_of_threads; tpos++) {
    if(args_table[tpos].free != 2) {
      if (pthread_join(threads[tpos], &status)) {
	perror("pthread_join");
	exit(1);
      }
    }
    free(args_table[tpos].grna);
  }
  // free grna
  free(grna);
}

void TrieAMatchSequenceThreads(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, char *pam, int uppercaseOnly, int threads, int printOnly)
{
  faread_struct *fas = installFastaReader(fname, uppercaseOnly);
  grna_list *g;
  int pamlen = strlen(pam);
  int guidelen = PROTOSPACER_LENGTH+pamlen;

  FILE* outfh = open_file(outputName, "w");
  
  while(fastaReader(fas)) {
    fprintf(stdout,"[crisflash] Processing %s (length %lld) in %s ...\n", fas->header, fas->slen, fname);
    fflush(stderr);
    fastaReaderImproveSequence(fas);
    g = fastaSequenceToGRNAs(fas->header, fas->s, fas->sr, fas->slen, pam);   
    if(threads > 1) {
      GRNAsMatchThreaded(T, g, guidelen, pam, maxMismatch, outfh, outFileType, threads);
    }
    else {
      GRNAsMatch(T, g, guidelen, pam, maxMismatch, outfh, outFileType);      
    }
    freeGRNAs(g);
  }
  freeFastaReader(fas);
  fclose(outfh);
}
