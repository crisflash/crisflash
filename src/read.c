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
  if (!fh) {
    fprintf(stderr,"[crisflash] ERROR: failed to open '%s'!\n", fname);
    exit(1);
  }
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

void printGRNAs(grna_list *g, int guidelen, FILE *f) {
  /** Note that this function prints gRNAs only for one chromosome at the time. **/
  
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
      if(!TrieAdd(T, g->pstrand+(g->starts[i]), guidelen, chr_number, g->strands[i], g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence detected in chr=%s, pos=%lld, strand=%c. Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
	exit(1);
      }
    }
    else {
      k=g->starts[i]+guidelen-1;
      for(j=0;j<guidelen;j++) {
	grna[j]=g->nstrand[k-j];
      }
      if(!TrieAdd(T, grna, guidelen, chr_number, g->strands[i], g->starts[i], 0, &identical_seq)) {
	fprintf(stderr,"[crisflash] ERROR: Failed to add gRNA sequence in chr=%s, pos=%lld, strand=%c). Exiting!",T->chr_names[chr_number],g->starts[i], g->strands[i]);
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

void readFaToTrie(trie *T, char *genomefname, char *pam, char* outputfname, int append, int uppercaseOnly, int printGRNAsOnly) {
  faread_struct *fas = installFastaReader(genomefname, uppercaseOnly);
  grna_list *g;
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
  fprintf(f, "%s\t%lld\t%lld\t%s/%d/%d\t%f\t%c\n", chrom, start, end, gRNA, nr_of_exact_matches, nr_of_off_targets, m->score, strand);
}

void writeInCasOffinder(FILE* fh, char *gRNA, mcontainer *m, trie* T)
{
  /* Writes match information in Cas-OFFinder format */
  
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
      for(j=0; j < m->nodes[i]->hits;j++)
	{
	  fprintf(fh,"%s\t%s\t%d\t%s\t%c\t%d\n", gRNA4print, T->chr_names[m->nodes[i]->chrs[j]], m->nodes[i]->starts[j], seq, m->nodes[i]->strands[j], m->marray[i]);
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
    writeInCasOffinder(args_table->outfh, args_table->grna, m, args_table->T);
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
  int identical_seq = 0; // counts the sequences not added in the tree because they were already there
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
      writeInCasOffinder(outfh, grna, m, T);
    }
    else {
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
  struct arg_struct args_table[nr_of_threads];
  int tpos; // thread position, expected to be 0 <= tpos < nr_of_threads
  // int fre; // 1 if threads is fre, 0 otherwise.
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
  int chr_number;
  int pamlen = strlen(pam);
  int guidelen = PROTOSPACER_LENGTH+pamlen;
  mcontainer* m;
  int i;

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
