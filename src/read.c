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

/* Written by Adrien Jacquin, April 2017
	based on a Margus Lukk script     */

FILE* open_file(char *fname,char *mode) {
  FILE *fh = fopen(fname, mode);
  if (!fh)
        {
          fprintf(stderr,"[crisflash] ERROR: failed to open '%s'!\n", fname);
          exit(1);
        }
  fprintf(stderr,"[crisflash] Saving results to '%s'\n", fname);
  fflush(stderr);
  return fh;
}

int baseMatch(char seqnt, char pamnt) {
  /* Checks whether base in DNA sequence (seqnt) matches with base in PAM pattern (pamnt).
     Assumes uppercase values for pamnt.
     Matching is provided based on  IUPAC representation for a position on a DNA sequence
     and were taken from https://en.wikipedia.org/wiki/Nucleic_acid_notation.
  */
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
      fprintf(stderr,"Error! '%c' part of PAM sequence not a valid nucleotide\n", pamnt);
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
      fprintf(stderr,"Error! '%c' is not valid DNA base symbol!\n", seqnt);
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
      fprintf(stderr,"'%c' not valid DNA base symbol!\n", seqnt);
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
    default:
      fprintf(stderr,"'%c' not valid DNA base symbol!\n", seqnt);
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
      fprintf(stderr,"'%c' not valid DNA base symbol!\n", seqnt);
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

typedef struct grna_list {
  char *pstrand; // pointer to forward strand of fasta sequence
  char *nstrand; // pointer to reverse strand of fasta sequence (translation from forward)
  char *chr; // chromosome name
  long long length; // number of space/positions allocatd for gRNAs
  long long pos; // nr of gRNas
  long long *starts; // array of gRNA start positions
  char *strands; // array of gRNA strands
} grna_list;

void updateSizeGRNAs(grna_list *g) {
  g->pos++;
  if(g->pos == g->length) {
    g->length = g->length+300000;
    g->starts = realloc(g->starts,g->length*sizeof(long long));
    g->strands = realloc(g->strands,g->length*sizeof(char));
  }
}

void printGRNAs(grna_list *g, int guidelen, FILE *f) {
  /** Note that this function prints gRNAs only for one chromosome at the time. **/
  
  long long i=0;
  int j;
  long long k;
  char *grna = malloc(sizeof(char)*(guidelen+1));
  
  for(i=0;i<g->pos;i++) {
    // fprintf(stdout,"%s\t%ld\t%ld\t\t0\t%c\ti=%ld,g->pos=%ld\n", g->chr, g->starts[i], (g->starts[i]+guidelen), g->strands[i],i,g->pos);
    // fflush(stdout);
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
    fprintf(f,"%s\t%lld\t%lld\t%s\t0\t%c\n", g->chr, g->starts[i], (g->starts[i]+guidelen), grna,g->strands[i]);
  }
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
      if(!TrieAddNEW(T, g->pstrand+(g->starts[i]), guidelen, chr_number, g->strands[i], g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
	exit(1);
      }
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<guidelen;j++) {
	grna[j]=g->nstrand[k-j];
      }
      if(!TrieAddNEW(T, grna, guidelen, chr_number, g->strands[i], g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"Failed to add gRNA sequence in chr=%s, pos=%lld, strand=%c). Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
	exit(1);
      }
    }
  }
  free(grna);
}

void freeGRNAs(grna_list *g) {
  free(g->starts);
  free(g->strands);
  free(g);
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
      // fprintf(stdout,"%ld (in first %d)\n",nt, PROTOSPACER_LENGTH);
      // safetynet of avoiding to read over the 3' edge of the sequence
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
      // fprintf(stdout,"MIDDLE %ld\n",nt);
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

/*
void faSequenceToTrie(trie* T, char *header, char *s, char *sr, long long spos, char *pam, char *outputfname, int upper_case_only) {
  // This function has been deprecated we are using faSequenceToGRNA instead

  // variables for calculating PAMs
  int pamlen = strlen(pam); // length of the PAM motif
  int guidelen = pamlen+PROTOSPACER_LENGTH; // length of the gRNA  
  char *gRNA = malloc(sizeof(char)*(20+pamlen+1));
  int identical_seq = 0;
  int add = 0; // to catch return values from TrieAdd. 0 - gRNA not added to trie. 1 otherwise.
  int y, nt;
  int chr_number;

  // add header (chr name) to to list of chromosome names in T.
  chr_number = addChr(T, header, strlen(header));
    
  // For first 20bp, search for PAMs on negative strand only
  for(nt=0; nt<PROTOSPACER_LENGTH; nt++)
    {
      // safetynet of avoiding to read over the 3' edge of the sequence
      if ((nt+guidelen) > spos) { break; }
      y = 0;
      while(baseMatch(sr[nt+y], pam[pamlen-1-y]))
	{
	  y++;
	  if (y==pamlen) {
	    for(int j=0;j<guidelen;j++) {
	      if (sr[nt+guidelen-j-1] == 'N') { break; }
	      gRNA[j]=sr[nt+guidelen-j-1];
	    }
	    //	    add = TrieAdd(T, gRNA, guidelen, chr_number, '-', spos, 0, &identical_seq);
	    // report found gRNA on positive strand
	    fprintf(stdout,"%s\t%d\t%ld\t-\n",gRNA,header,nt);
	    break;
	  }
	}	 		  
    }
  
  // Search for PAM sites on both strands
  for(nt=20; nt<(spos-guidelen); nt++)
    {
      y = 0;
      // Search PAM sequences on positive strand.		
      while(baseMatch(s[nt+y], pam[y]))
	{
	  y++;
	  if (y==pamlen) {
	    // add = TrieAdd(T, s+(nt-PROTOSPACER_LENGTH), guidelen, chr_number, '-', spos, 0, &identical_seq);
	    for(int j=0;j<(guidelen);j++) {
	      if(s[nt-PROTOSPACER_LENGTH+j] == 'N') { break; }
	      gRNA[j]=s[nt-PROTOSPACER_LENGTH+j];
	    }
	    // report found gRNA on positive strand
	    fprintf(stdout,"%s\t%s\t%ld\t+\n",gRNA,header,(nt-PROTOSPACER_LENGTH));
	    break;
	  }
	} // end of while
      
      // Search PAM on negative strand
      y = 0;
      while(baseMatch(sr[nt+y], pam[pamlen-1-y]))
	{
	  y++;
	  if (y==pamlen) {		    
	    for(int j=0;j<(guidelen);j++) {
	      if (sr[nt+guidelen-j-1] == 'N') { break; }
	      gRNA[j]=sr[nt+guidelen-j-1];
	    }
	    // add = TrieAdd(T, gRNA, guidelen, chr_number, '-', spos, 0, &identical_seq);
	    // report found gRNA on negative strand
	    fprintf(stdout,"%s\t%s\t%ld\t-\n",gRNA,header,nt);
	    break;
	  }
	} // end of while	     		
    }
  
  // Search for PAMs only on positive strand for the remaining bit of the sequence
  for(nt=(spos-guidelen); nt<(spos-pamlen); nt++)
    {
      // safetynet of exceptionally shot DNA sequence
      if ((nt-PROTOSPACER_LENGTH) < 0) { continue; }			       
      y = 0;
      while(baseMatch(sr[nt+y], pam[pamlen-1-y]))
	{
	  y++;
	  if (y==pamlen) {
	    // add = TrieAdd(T, s+(nt-PROTOSPACER_LENGTH), guidelen, chr_number, '-', spos, 0, &identical_seq);
	    for(int j=0;j<(guidelen);j++) {
	      if(s[nt-PROTOSPACER_LENGTH+j] == 'N') { break; }
	      gRNA[j]=s[nt-PROTOSPACER_LENGTH+j];
	    }
	    fprintf(stdout,"%s\t%s\t%ld\t+\n",gRNA,header,nt); // <- this is old sanity check
	    // fprintf(f, "%s\t%d\t%d\t%s\t%f\t%c\n", T->chr_names[chr_number], chr_pos-22, chr_pos+1, seq, 0., '+'); }
	    break;
	  }
	}	  
    }
  
  free(gRNA);
}
*/


/* Struct returned by installFastaReader, used by fastaReader and freed by freeFastaReader.*/
typedef struct faread_struct
{
  long long limit; // maximul length allocated for s and sr at any time
  char *s; // fasta sequence + strand
  char *sr; // fasta sequence - strand
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
  
  // open file and read first chunk of data
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
	  fprintf(stderr,"Fasta header too long! Edit BUFFER_SIZE and re-compile! EXITING!\n");
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
	// fprintf(stderr,"Reallocated mem to %ld at character position %ld\n", limit, spos);
      }
      fas->i++;
    }
    // if end of buffer, read more from file handler
    if (fas->i == fas->bytes_read) {
      fas->bytes_read = read(fas->fd, fas->buffer, BUFFER_SIZE);
      // fprintf(stderr,"Bytes read=%d\n",fas->bytes_read);
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
    fprintf(stderr,"[crisflash] No fasta header. Exiting!\n");    
    exit(1);
  }
  
  return 1;
}

void readFaToTrieNEW(trie *T, char *genomefname, char *pam, char* outputfname, int append, int uppercaseOnly, int printGRNAsOnly) {
  faread_struct *fas = installFastaReader(genomefname, uppercaseOnly);
  grna_list *g;
  int chr_number;
  int guidelen = PROTOSPACER_LENGTH+strlen(pam);
  FILE* outfh;
  
  if(printGRNAsOnly) {
    outfh = open_file(outputfname, "w");
  }

  while(fastaReader(fas)) {
    fprintf(stderr,"[crisflash] Processing %s (length %lld) in %s ...\n", fas->header, fas->slen, genomefname);
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

void writeInBed_lighterNEW(FILE* f, mcontainer *m, char* chrom, long long start, long long end, char* gRNA, char strand, trie* T)
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
  fprintf(f, "%s\t%lld\t%lld\t%s/%d/%d\t%f\t%c\n", chrom, start, end, gRNA, nr_of_exact_matches, nr_of_off_targets, m->score, strand);
}

void writeInCasOffinder(FILE* fh, char *gRNA, mcontainer *m, trie* T)
{
  /*************

    TODO:
      - The number of hits found is substantially less than by cas-offinder. Check why! Could it be that cas-offinder is considering masked regions as well?
      - cas-offinder reports PAM sequence as well for candidate area. We don't as it is not saved/stored in trie. This needs to be changed!

   *************/ 

  // in case mcontainer is empty
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
  // gRNA4print[T->readlen] = '\0';

  i=0;  
  char *seq = (char *) malloc(sizeof(char)*(strlen(m->sarray[i])+1));  

  while(i < m->pos)
    {
      strcpy(seq,m->sarray[i]);
      // fprintf(stdout,"Befre: %s\n",seq);
      for(char* c=seq; (*c=toupper(*c)); ++c) ;      
      // fprintf(stdout,"After: %s\n",seq);
      
      // fprintf(stdout,"%s\t[",seq);
      int j = 0;
      while (j < m->marray[i]-1)
	{
	  // fprintf(stdout,"%d,", m->parray[i][j]);
	  seq[m->parray[i][j]-1] = tolower(seq[m->parray[i][j]-1]);
	  j++;
	}
      if (m->marray[i] > 0)
	{
	  seq[m->parray[i][j]-1] = tolower(seq[m->parray[i][j]-1]);
	  // fprintf(stdout,"%d]\t", m->parray[i][j]);
	}
      //      else { fprintf(stdout,"]\t"); }
      // print all hits for a particular gRNA
      //      fprintf(stdout,"Nr of hit locations: %d\n", m->nodes[i]->hits);
      //      fflush(stdout);
      //      fprintf(stdout,"gRNA=%s\n", gRNA);
      //      fflush(stdout);
      for(j=0; j < m->nodes[i]->hits;j++)
	{
	  /*
	  fprintf(stdout,"j=%d\n",j);
	  fprintf(stdout,"Chr=%d\n", T->chr_names[m->nodes[i]->chrs[j]]);
	  fflush(stdout);
	  fprintf(stdout,"Strand=%d\n", m->nodes[i]->starts[j]);
	  fflush(stdout);
	  fprintf(stdout,"Mismatches in sequence=%s\n",seq);
	  fflush(stdout);
	  fprintf(stdout,"Strand=%c\n", m->nodes[i]->strands[j]);
	  fflush(stdout);
	  fprintf(stdout,"mismatches=%d\n",m->marray[i]);
	  fflush(stdout);
	  */
	  fprintf(stdout,"%s\t%s\t%d\t%s\t%c\t%d\n", gRNA4print, T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j], seq, m->nodes[i]->strands[j], m->marray[i]);
	  fprintf(fh,"%s\t%s\t%d\t%s\t%c\t%d\n", gRNA4print, T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j], seq, m->nodes[i]->strands[j], m->marray[i]);
	  // fflush(stdout);
	}
      //      fprintf(stdout,"%f\n", m->off_scores[i]);
      i++;
    }
  free(gRNA4print);
  free(seq);
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
					  add = TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
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
					  add = TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '-', chr_pos+23, 0, &identical_seq);
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
				  add = TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
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
				  add = TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '-', chr_pos+23, 0, &identical_seq);				    
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

	trie* T = TrieCreate(PROTOSPACER_LENGTH);
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
					  TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '+', chr_pos-22, 0, &identical_seq); // 22 because we stopped at the first G
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
				// 		TrieAdd(T, seq, PROTOSPACER_LENGTH, chr_number, '-', chr_pos+23, 0, &identical_seq);
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

void TrieAMatchItself(trie* T, int maxMismatch, char* outputName, char* outputMatchResult, int outFileType, char *pam)
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
		for (int j = 0; j < PROTOSPACER_LENGTH; j++)
		{
			gRNA[j] = line[k];
			k++;
		}
		gRNA[PROTOSPACER_LENGTH] = '\0';
		while (line[k] != '\t') { k++; }
		k++;
		strand = line[k];
		m = TrieAMatch(T, gRNA, T->readlen, maxMismatch);
		if (f2 != NULL)
		{
		  if(outFileType==1) {
		    // write the results in a .bed file
		    writeInBed_lighter(f2, m, chrom, start, end, gRNA, strand, T);
		  }
		  if (outFileType==2) {
		    writeInCasOffinder(f2, gRNA, m, T);
		  }
			// mcontainer_print(m);
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
		m = TrieAMatch(args->T, args->gRNA, PROTOSPACER_LENGTH, args->nr_of_mismatches);
		if (args->output != NULL)
		{
			writeInBed_lighter(args->output, m, args->chrom, args->start, args->end, args->gRNA, args->strand, args->T); // write the results in a file
			// printf("Writing results in a file, thread #%d, counter = %d\n", args->i, args->cnt);
		}

		// free memory
		mcontainer_free(m);
	}

	// (void) arg;
	args->free = 1;
	pthread_exit(NULL);
}

void *thread_workerN(void *arg)
{
  // cast the argument into its real nature
  struct arg_structN *args_tableN = (struct arg_structN *)arg;

  /*
    args_tableN[tpos].T = T;
    args_tableN[tpos].i = tpos; // thread number
    args_tableN[tpos].maxMismatch = maxMismatch;
    args_tableN[tpos].free = 1;
    args_tableN[tpos].outfh = outfh;
    args_tableN[tpos].outFileType = outFileType;
    args_tableM[tpos].free = 0;
    strcpy(args_tableN[tpos].grna, gRNA);
    args_tableN[tpos].chr = T->chr_names[g->chr];
    args_tableN[tpos].start = g->starts[i];
    args_tableN[tpos].end = g->starts[i] + T->readlen;
    args_tableN[tpos].strand = g->strands[i];
  */

  // match
  mcontainer* m = TrieAMatch(args_tableN->T, args_tableN->grna, args_tableN->T->readlen, args_tableN->maxMismatch);
  if(args_tableN->outFileType == 2) {
    writeInCasOffinder(args_tableN->outfh, args_tableN->grna, m, args_tableN->T);
  }
  else {
    writeInBed_lighterNEW(args_tableN->outfh, m, args_tableN->chr, args_tableN->start, args_tableN->end, args_tableN->grna, args_tableN->strand, args_tableN->T);
  }
  // mcontainer_print(m);
  mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
  // label thread being free/available
  args_tableN->free = 1;

  pthread_exit(NULL);
}

void TrieAMatchItself_thread(trie *T, int maxMismatch, char* outputName, char* outputMatchResult, int outFileType, int nr_of_threads, char *pam)
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
			args_table[i].gRNA = malloc(sizeof(char) * (PROTOSPACER_LENGTH + 1));
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
			for (int j = 0; j < PROTOSPACER_LENGTH; j++)
			{
				args_table[i].gRNA[j] = line[k];
				k++;
			}
			args_table[i].gRNA[PROTOSPACER_LENGTH] = '\0';
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

void TrieAMatchSequence_thread(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, int nr_of_threads, char *pam)
{
	// Creates a second tree and a BED file, based on the sequence file.
	// Calls the multiple threads function to match the sequence against the tree.
	// Writes the results in the BED file outputName.
	double end_time;
	double start_time = now();
	int n;
	int append = 0;

	trie* T2 = TrieCreate(PROTOSPACER_LENGTH);

	readFaToTrie_UpperCaseOnly(T2, fname, "sequenceTree", append); // Create a tree and a BED file for the sequence
	TrieAMatchItself_thread(T, maxMismatch, "sequenceTree", outputName, outFileType, nr_of_threads, pam);
	T2 = TrieDestroy(T2);
	n = remove("sequenceTree.bed"); // delete the BED file we needed temporary
	if (n < 0) { perror("[crisflash] ERROR deleting sequence BED file"); }
}

void TrieAMatchSequence(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, char *pam)
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
						gRNA = malloc(sizeof(char)*(PROTOSPACER_LENGTH + 1));
						gRNA[PROTOSPACER_LENGTH] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<PROTOSPACER_LENGTH; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, PROTOSPACER_LENGTH, maxMismatch);
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
								mcontainer_print(m);
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
						gRNA = malloc(sizeof(char)*(PROTOSPACER_LENGTH + 1));
						gRNA[PROTOSPACER_LENGTH] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<PROTOSPACER_LENGTH; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, PROTOSPACER_LENGTH, maxMismatch);
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
								mcontainer_print(m);
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

void GRNAsMatch(trie *T, grna_list *g, int guidelen, char *pam, int maxMismatch, FILE *outfh, int outFileType) {
  /** Matches each gRNA in g to trie T. Most of the code strucure is more or less copy/paste from GRNAsToTrie and printGRNAs **/
  
  // the loop structure for this function is identical to printGRNAs.
  long long i=0;
  long long k;
  int j;
  char *grna = malloc(sizeof(char)*(guidelen+1)); // here the gRNA sequence will be gRNA protospacer + PAM, not the sequence as in candidate area!
  int identical_seq = 0; // counts the sequences not added in the tree because they were already there
  mcontainer* m;
  int pamlen = strlen(pam);

  for(i=0;i<g->pos;i++) {

    // create gRNA sequence consisting protospacer sequence + PAM sequence.
    if(g->strands[i] == '+') {
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	fprintf(stderr,"grna[%d] + '%c'\n",j,g->pstrand[g->starts[i]+j]);
	fflush(stderr);
	grna[j]=g->pstrand[g->starts[i]+j];
      }
    }
    else {
      k=g->starts[i]+guidelen-1-pamlen;
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	fprintf(stderr,"grna[%d] - '%c'\n",j,g->nstrand[k-j]);
	fflush(stderr);
	grna[j]=g->nstrand[k-j];
      }
    }
    // add PAM sequence
    for(j=0;j<pamlen;j++) {
      fprintf(stderr,"grna[%d] pam pos '%c'\n",j+PROTOSPACER_LENGTH,pam[j]);
      fflush(stderr);
      grna[PROTOSPACER_LENGTH+j]=pam[j];
    }
    
    // match grna
    // fprintf(stderr,"Tyring to match %s ...\n",grna);
    fprintf(stderr,"Next will try to match %s, readlen=%d\n",grna,T->readlen);
    fflush(stderr);
    m = TrieAMatch(T, grna, T->readlen, maxMismatch);
    fprintf(stderr,"Match completed!\n");
    fflush(stderr);
    if(outFileType==2) {
      fprintf(stderr,"Writing output in CasOffinder format\n");
      writeInCasOffinder(outfh, grna, m, T);
    }
    else {
      fprintf(stderr,"Writing output in bed format\n");
      // write the results in bed
      writeInBed_lighterNEW(outfh, m, g->chr, g->starts[i], g->starts[i]+guidelen, grna, g->strands[i], T);
    }
    // mcontainer_print(m);
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
  char *grna;
  grna = malloc(sizeof(char)*(guidelen+1)); // here the gRNA sequence will be gRNA protospacer + PAM, not the sequence as in candidate area!
  int identical_seq = 0; // counts the sequences not added in the tree because they were already there
  mcontainer* m;
  int pamlen = strlen(pam);

  // Install threading related variables and attributes
  pthread_t threads[nr_of_threads]; // list of threads
  pthread_attr_t attr; // thread attribute necessary for thread tracking on joining.
  // initialize and set thread detached attribute
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  void *status;
  int res; // to catch thread exit result 
  struct arg_structN args_tableN[nr_of_threads];
  int tpos; // thread position, expected to be 0 <= tpos < nr_of_threads
  // int fre; // 1 if threads is fre, 0 otherwise.
  // install varables in args table
  for(tpos=0;tpos<nr_of_threads;tpos++) {  
    args_tableN[tpos].T = T;
    args_tableN[tpos].i = tpos; // thread number
    args_tableN[tpos].maxMismatch = maxMismatch;
    args_tableN[tpos].free = 1;
    args_tableN[tpos].outfh = outfh;
    args_tableN[tpos].outFileType = outFileType;
  }
  tpos = 0;
  
  // we will process one gRNA at the time
  for(i=0;i<g->pos;i++) {

    // create gRNA sequence consisting protospacer sequence + PAM sequence.
    if(g->strands[i] == '+') {
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	fprintf(stderr,"grna[%d] + '%c'\n",j,g->pstrand[g->starts[i]+j]);
	fflush(stderr);
	grna[j]=g->pstrand[g->starts[i]+j];
      }
    }
    else {
      k=g->starts[i]+guidelen-1-pamlen;
      for(j=0;j<PROTOSPACER_LENGTH;j++) {
	fprintf(stderr,"grna[%d] - '%c'\n",j,g->nstrand[k-j]);
	fflush(stderr);
	grna[j]=g->nstrand[k-j];
      }
    }
    // add PAM sequence
    for(j=0;j<pamlen;j++) {
      fprintf(stderr,"grna[%d] pam pos '%c'\n",j+PROTOSPACER_LENGTH,pam[j]);
      fflush(stderr);
      grna[PROTOSPACER_LENGTH+j]=pam[j];
    }
    
    // find next free thread for matching the found gRNA
    res = 0;
    while(!res) {
      if(tpos == nr_of_threads) { tpos = 0; }
      if(args_tableN[tpos].free == 1) {
	if (pthread_join(threads[tpos], &status)) {
	  perror("pthread_join");
	  exit(1);
	}
	// load_args_tableN with relevant values
	args_tableN[tpos].free = 0;
	strcpy(args_tableN[tpos].grna, grna);
	args_tableN[tpos].chr = g->chr;
	args_tableN[tpos].start = g->starts[i];
	args_tableN[tpos].end = g->starts[i] + T->readlen;
	args_tableN[tpos].strand = g->strands[i];

	// SIIA JAI POOLELI. See SIIN PEAKS T**TAMA AGA POLE PROOVINUD

	if(pthread_create(&threads[i], &attr, thread_workerN, &args_tableN[tpos])) {
	  perror("pthread_create");
	  exit(1);
	}
	res = 1;
      }
      tpos++;
    }

  } // end of processing gRNAs

  // join the remaining threads
  for(tpos=0;tpos<nr_of_threads; tpos++) {
    if(args_tableN[tpos].free != 1) {
      if (pthread_join(threads[tpos], &status)) {
	perror("pthread_join");
	exit(1);
      }
    }
  }
  // free grna
  free(grna);
}

void TrieAMatchSequenceNEWThreads(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, char *pam, int uppercaseOnly, int threads, int printOnly)
{
  faread_struct *fas = installFastaReader(fname, uppercaseOnly);
  grna_list *g;
  int chr_number;
  int pamlen = strlen(pam);
  int guidelen = PROTOSPACER_LENGTH+pamlen;
  mcontainer* m;
  int i;

  FILE* outfh = open_file(outputName, "w");
  
  while(fastaReader(fas)) {
    fprintf(stderr,"[crisflash] Processing %s (length %lld) in %s ...\n", fas->header, fas->slen, outputName);
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

void TrieAMatchSequenceNEWold(trie* T, char* fname, int maxMismatch, char* outputName, int outFileType, char *pam)
{
  fprintf(stdout,"\n\nStarting with NEW matcher\n\n");
  //  exit(0);
  
	// Finds all the gRNAs in candidate fasta file.
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
			  fprintf(stdout,"%c",buffer[x]);
				bool_chr = false;
				chr[chrp] = '\0';
				chrp = 0;
				chr_pos = 1;
			}
			else if(!bool_chr)
			{
			  // THIS IS WHERE THE MAGIC HAPPENS. WE CHECK WHAT THE LETTER FOR THE BASE IS AND ACT ACCORDINGLY.			  
			  /*
				if ((buffer[x] == 'G' && buffer[x+1] == 'G') || (buffer[x] == 'G' && buffer[x+1] == '\n' && buffer[x+2] == 'G'))
				{
					if (buffer[x+1] == '\n') { seq = fSeq(buffer, x+2, 23, 1); }
					else { seq = fSeq(buffer, x+1, 23, 1); }
					if ((seq != NULL) && (chr_pos-23 > -1))
					{
						gRNA = malloc(sizeof(char)*(PROTOSPACER_LENGTH + 1));
						gRNA[PROTOSPACER_LENGTH] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<PROTOSPACER_LENGTH; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, PROTOSPACER_LENGTH, maxMismatch);
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
								mcontainer_print(m);
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
						gRNA = malloc(sizeof(char)*(PROTOSPACER_LENGTH + 1));
						gRNA[PROTOSPACER_LENGTH] = '\0';
						if (gRNA != NULL) // if malloc worked
						{
							for (int i=0; i<PROTOSPACER_LENGTH; i++)
							{
								gRNA[i] = seq[i];
							}
							// printf("gRNA %s found!\nSequence(s) in the trie that matched:\n", gRNA);
							m = TrieAMatch(T, gRNA, PROTOSPACER_LENGTH, maxMismatch);
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
								mcontainer_print(m);
							}
							free(gRNA);
							mcontainer_free(m); // TrieAMatch creates a new container each time, so we have to free it here
							nrRNA++;
						}
						free(seq);
					}
				}
			  */
			        if (buffer[x] == '>')
				{
				  fprintf(stdout,"%c",buffer[x]);
				  bool_chr = true;
				}
			} // end of else if (!bool_chr)
			else // we write the chromosome name in 'chr' :
			{
			  fprintf(stdout,"%c",buffer[x]);
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
	exit(0);
}
