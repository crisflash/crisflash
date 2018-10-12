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

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "nary_tree.h"

/* 
   Written by Adrien Jacquin, April 2017
   based on a prototype developed by Margus Lukk.
*/

/***************************************
          BASE MATCHING FUNCTIONS
***************************************/

int baseMismatch(char seqnt, char pamnt) {
  /* Checks whether base in DNA sequence (seqnt) mismatches with base in PAM pattern (pamnt).
     Assumes uppercase values. Returns 1 if mismatch, 0 otherwise.
     Matching is provided based on  IUPAC representation for a position on a DNA sequence
     and were taken from https://en.wikipedia.org/wiki/Nucleic_acid_notation.
  */
  switch(pamnt)
    {
    case 'N':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	case 'G':
	  return 0;
	  break;
	case 'C':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'G':
      switch(seqnt)
	{
	case 'G':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'A':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'T':
      switch(seqnt)
	{
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'C':
      switch(seqnt)
	{
	case 'C':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'R':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'G':
	    return 0;
	    break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'Y':
      switch(seqnt)
	{
	case 'C':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'V':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'C':
	  return 0;
	  break;
	case 'G':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'W':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'S':
	switch(seqnt)
	  {
	  case 'G':
	    return 0;
	    break;
	  case 'T':
	    return 0;
	    break;
	  default:
	    return 1;
	    break;
	  }
	break;
    case 'M':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'C':
	  return 0;
	    break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'K':
      switch(seqnt)
	{
	case 'G':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'B':
      switch(seqnt)
	{
	case 'C':
	  return 0;
	  break;
	case 'G':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    case 'D':
      switch(seqnt)
	{
	  case 'A':
	    return 0;
	    break;
	  case 'G':
	    return 0;
	    break;
	  case 'T':
	    return 0;
	    break;
	  default:
	    return 1;
	    break;
	}
      break;
    case 'H':
      switch(seqnt)
	{
	case 'A':
	  return 0;
	  break;
	case 'C':
	  return 0;
	  break;
	case 'T':
	  return 0;
	  break;
	default:
	  return 1;
	  break;
	}
      break;
    default:
      fprintf(stderr,"[crisflash] ERROR: '%c' part of PAM sequence not a valid nucleotide\n", pamnt);
      exit(1);
    }
}

/***************************************
      MATCHING CONTAINER FONCTIONS
***************************************/

mcontainer* mcontainer_install(int nelem, int calen)
{
  /* Installs mcontainer to keep track of gRNA matching information */
	mcontainer *m = malloc(sizeof(struct mcontainer));
	m->nodes = malloc( sizeof(trieNodeTT*) * nelem);
	m->sarray = malloc( sizeof(char**) * nelem);
	m->marray = malloc(sizeof(int) * nelem);
	m->parray = malloc(sizeof(int*) * nelem);
	m->off_scores = NULL;
	m->score = 0;
	m->nelem = nelem;
	m->pos = 0;
	m->calen = calen;
	m->capos = 0;
	m->matches = 0;
	return m;
}

void mcontainer_free(mcontainer *m)
{
  /* Frees mcontainer. */
	int i = 0;

	while(i < m->pos)
	{
		free(m->sarray[i]);
		m->nodes[i] = NULL;
		if (m->marray[i] > 0)
		{
			free(m->parray[i]);
			m->parray[i] = NULL;
		}
		i++;
	}
	if (m->off_scores)
	{
		free(m->off_scores);
	}
	free(m->sarray);
	free(m->nodes);
	free(m->marray);
	free(m->parray);
	free(m);
}

int mcontainer_add_str(mcontainer *m, int slen)
{
  // Returns the position in m->sarray where the char array with length slen was created.
  char *seq = malloc(sizeof(char)*(slen+1));
  seq[slen] = '\0';
  m->sarray[m->pos] = seq;
  m->marray[m->pos] = 1;
  m->pos++;
  return m->pos-1 ;
}

void mcontainer_add_nt(mcontainer *m, mcontainer *m_new, int mpos, int seqpos, char nt, int mismatch, int max_mismatches, trieNodeTT* nodeT)
{
  /*
    Creates new character array for a particular character array in m->sarray[mpos] and adds new nt to it.
    mpos is the position of sequence in m sarray
    seqpos is position of char in seq array
    nt is the char to be added to a character array
    mismatch: 0 if no mismatch, 1 otherwise.
  */
        char *seq;

	if(m_new->pos >= m_new->nelem)
	{
		fprintf(stderr, "[crisflash] ERROR: Not enough space allocated for new approximate sequences. Exiting!\n");
		exit(1);
	}

	if( (m->marray[mpos] + mismatch) <= max_mismatches) // we check that there are not too much mismatches
	{
		// seq = m->sarray[mpos];
		seq = malloc(sizeof(char) * (m_new->calen + 1)); // seq is initialized with the size of readlen
		seq[m_new->calen] = '\0'; // string terminator
		memcpy(seq, m->sarray[mpos], seqpos); // seqpos here is to indicate the size of the original seq to copy
		m_new->sarray[m_new->pos] = seq; // add seq in the new container
		seq[seqpos] = nt; // add the nt to the sequence
		// char array position for a whole container
		if(m_new->capos != (seqpos+1))
		{
			m_new->capos = seqpos + 1;
		}
		m_new->marray[m_new->pos] = m->marray[mpos] + mismatch;
		m_new->nodes[m_new->pos] = nodeT;
		if (m_new->marray[m_new->pos] > 0) // we have to copy the positions list
		{
			m_new->parray[m_new->pos] = malloc(sizeof(int) * m_new->marray[m_new->pos]);
			int i = 0;
			while (i < m->marray[mpos]) // copy of the former positions
			{
				m_new->parray[m_new->pos][i] = m->parray[mpos][i];
				i++;
			}
			if (mismatch == 1)
			{
				// add the new mismatch position
				m_new->parray[m_new->pos][i] = seqpos + 1; // +1 because we want the first position to be 1, and not 0
			}
		}
		m_new->pos++;
	}
}

void mcontainer_print(mcontainer *m)
{
  /* Prints the content of match container. */
	int i = 0;
	char *seq;
	while(i < m->pos)
	{
		seq = m->sarray[i];
		fprintf(stdout,"%s\t[",seq);		
		int j = 0;
		while (j < m->marray[i]-1)
		{
		  fprintf(stdout,"%d,", m->parray[i][j]);
		  j++;
		}
		if (m->marray[i] > 0) { fprintf(stdout,"%d]\t", m->parray[i][j]); }
		else { fprintf(stdout,"]\t"); }
		fprintf(stdout,"%f\n", m->off_scores[i]);
		i++;
	}
	fprintf(stdout,"gRNA score = %f\n", m->score);
}

void mcontainer_print_pretty(mcontainer *m, int mismatches)
{
  /* Prints the number of gRNAs for each number of mismatches */
	int i = 0;
	int *matches = malloc(sizeof(int) * (mismatches+1));

	// installize matches array
	while(i < (mismatches+1))
	{
		matches[i] = 0;
		i++;
	}

	// for each number of mismatches, find the number of matches
	i = 0;
	while(i < m->pos)
	{
		// count number of different sequences instead of matches
		for (int j = 0; j < m->nodes[i]->hits; j++) // do not forget when there are several hits
		{
			matches[ m->marray[i] ] += 1;
		}
		i++;
	}

	// print match info for each nr of mismatches
	i = 0;
	while(i < (mismatches+1))
	{
		fprintf(stdout,"\t%d", matches[i]);
		i++;
	}
	fprintf(stdout,"\n\n");

	free(matches);
}

void mcontainer_score(mcontainer *m, trie* T)
{
	// Calculate the off-target scores (m->off_scores table) and the gRNA score (m->score)
	// Each off-target score corresponds to a gRNA found on the input genome
	// The gRNA score corresponds to the gRNA found in the sequence and that is being matched in the tree

	// We define the W vector (experimentally-determined):
	float W[20] = {0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583};
	m->off_scores = malloc(sizeof(float) * m->pos);
	float sum = 0;
	for (int i = 0; i < m->pos; i++)
	{
		int pos_max = 0; // declare min and max, in case there is no mismatch
		int pos_min = 0;
		float d = 0;
		m->off_scores[i] = 1;
		// parse the mismatch positions to obtain the max & min positions, and calculate the first term of the product:
		for (int j = 0; j < m->marray[i]; j++)
		{
			if (j == 0) // initialize min and max
			{
				pos_max = m->parray[i][j];
				pos_min = m->parray[i][j];
			}
			if (m->parray[i][j] > pos_max) { pos_max = m->parray[i][j]; }
			if (m->parray[i][j] < pos_min) { pos_min = m->parray[i][j]; }
			for (int k = 0; k < m->nodes[i]->hits; k++) // we must take into account when there are several hits
			{
				m->off_scores[i] *= (1 - W[ m->parray[i][j] - 1 ]);
			}
		}
		if (m->marray[i] > 0)
		{
			d = (pos_max - pos_min)/(m->marray[i]); // this is the mean pairwise distance between mismatches
			m->off_scores[i] /= (((19 - d)/19)*4 + 1); // second term of the product
			m->off_scores[i] /= (m->marray[i] * m->marray[i]); // last term of the product
		}
		// if (m->marray[i] == 0) // there is one exact match in tree
		// {
		// 	trieNodeTT* nt = m->nodes[i];
		// 	for (int k = 0; k < nt->hits; k++) // we must take into account when there are several hits
		// 	{
		// 		printf("The gRNA matches EXACTLY in chromosome %s:%d-%d %c\n", T->chr_names[nt->chrs[k]], nt->starts[k], nt->starts[k]+23, nt->strands[k]);
		// 	}
		// }
		m->off_scores[i] *= 100; // percentage
		sum += m->off_scores[i]; // we sum the off-target scores and the on-target scores
	}

	// After calculating every off-target score, we calculate the gRNA score:
	if (sum != 0) { m->score = 100. / sum; }
	else { m->score = 0; }
}

/**************************
        VCF FONCTIONS
**************************/

VCF* VCF_install(FILE* f_vcf)
{
  /* Installs VCF struct. */
	VCF* vcf = malloc(sizeof(struct VCF));
	char* line = NULL;
	char *str, *ref, *alt, *hap;
	size_t len = 0;
	int read = getline(&line, &len, f_vcf);
	while (read != -1 && line[0] == '#')
	{
		read = getline(&line, &len, f_vcf);
	}
	if (read == -1) { return NULL; }
	// Split the line and store the data:
	str = strtok(line, "\t");
	vcf->chr = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->chr,str);
	for (int i = 0; i < 9; i++)
	{
		str = strtok(NULL, "\t");
		if (i==0) { vcf->pos = atoi(str);}
		if (i==2) { ref = str; }
		if (i==3) { alt = str; }
		if (i==8) { hap = str; }
	}
	str = strtok(ref,",");
	vcf->ref = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->ref,str);
	str = strtok(alt,",");
	vcf->alt = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->alt,str);
	str = strtok(hap,"\n");
	vcf->hap = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->hap,str);
	vcf->dif = strlen(vcf->alt)-strlen(vcf->ref);
	if (vcf->dif < 0) { vcf->ub = vcf->pos - vcf->dif; }
	else { vcf->ub = vcf->pos + vcf->dif; }
	if (line) { free(line); }
	return vcf;
}

void VCF_destroy(VCF* vcf)
{
  /* Frees VCF struct. */
	if (vcf)
	{
		if (vcf->chr) { free(vcf->chr); }
		if (vcf->ref) { free(vcf->ref); }
		if (vcf->alt) { free(vcf->alt); }
		if (vcf->hap) { free(vcf->hap); }
		free(vcf);
	}
}

int VCF_update(FILE* f_vcf, VCF* vcf)
{
  /* Updates VCF struct. */
	/* First we delete the old data */
	free(vcf->chr);
	free(vcf->ref);
	free(vcf->alt);
	free(vcf->hap);
	vcf->chr = NULL;
	vcf->ref = NULL;
	vcf->alt = NULL;
	vcf->hap = NULL;
	char* line = NULL;
	char *str, *ref, *alt, *hap;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f_vcf);
	while (read != -1 && line[0] == '#')
        {       
                read = getline(&line, &len, f_vcf);
        }
	if (read == -1) { return -1; }
	// Split the line and store the data:
	str = strtok(line, "\t");
	vcf->chr = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->chr,str);
	for (int i = 0; i < 9; i++)
	{
		str = strtok(NULL, "\t");
		if (i==0) { vcf->pos = atoi(str);}
		if (i==2) { ref = str; }
		if (i==3) { alt = str; }
		if (i==8) { hap = str; }
	}
	str = strtok(ref,",");
	vcf->ref = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->ref,str);
	str = strtok(alt,",");
	vcf->alt = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->alt,str);
	str = strtok(hap,"\n");
	vcf->hap = malloc(sizeof(char)*(strlen(str)+1));
	strcpy(vcf->hap,str);
	vcf->dif = strlen(vcf->alt)-strlen(vcf->ref);
	if (vcf->dif < 0) { vcf->ub = vcf->pos - vcf->dif; }
	else { vcf->ub = vcf->pos + vcf->dif; }
	if (line) { free(line); }
	return 0;
}

void VCF_print(VCF* vcf)
{
  /* Prints VCF struct. */
	if (vcf)
	{
		fprintf(stdout, "VCF line: chr=%s\tpos=%d\tub=%d\tref=%s\talt=%s\thap=%s\n", vcf->chr, vcf->pos, vcf->ub, vcf->ref, vcf->alt, vcf->hap);
	}
	else { fprintf(stderr, "ERROR reading vcf structure\n"); }
}

VCF* VCF_copy(VCF* vcf)
{
  /* Returns copy of the vcf struct */
	VCF* vcf2 = malloc(sizeof(struct VCF));
	vcf2->chr = malloc(sizeof(char)*(strlen(vcf->chr)+1));
	vcf2->ref = malloc(sizeof(char)*(strlen(vcf->ref)+1));
	vcf2->alt = malloc(sizeof(char)*(strlen(vcf->alt)+1));
	vcf2->hap = malloc(sizeof(char)*(strlen(vcf->hap)+1));
	strcpy(vcf2->chr,vcf->chr);
	strcpy(vcf2->ref,vcf->ref);
	strcpy(vcf2->alt,vcf->alt);
	strcpy(vcf2->hap,vcf->hap);
	vcf2->pos = vcf->pos;
	vcf2->ub = vcf->ub;
	vcf2->dif = vcf->dif;
	return vcf2;
}

/**************************
       TRIE FUNCTIONS
**************************/

trieNodeT* TrieCreateNode(char nt)
{
	// Adds node to trieNodeT
	trieNodeT* a = (trieNodeT*) malloc(sizeof(trieNodeT));
	a->brother = NULL;
	a->child = NULL;
	a->nt = nt;
	return a;
}

trieNodeTT* TrieCreateNodeTerminal(char nt)
{
  /* Adds terminal node to Trie. */
	trieNodeTT* a = (trieNodeTT*) malloc(sizeof(trieNodeTT));
	a->allocated_array_len=10;
	a->chrs = malloc(sizeof(int)*a->allocated_array_len);
	a->starts = malloc(sizeof(int)*a->allocated_array_len);
	a->strands = malloc(sizeof(char)*a->allocated_array_len);
	a->brother = NULL;
	a->hits = 0;
	a->nt = nt;
	a->score = 0;
	return a;
}

trieNodeT* TrieAddNode(trieNodeT* level, char nt)
{
	// Add a new node as the last brother of the level
	// Return this new node address
	if (level == NULL)
	{
		trieNodeT* n = TrieCreateNode(nt);
		return n;
	}
	if (level->brother)
	{
		return TrieAddNode(level->brother, nt);
	}
	trieNodeT* n = TrieCreateNode(nt);
	level->brother = n;
	return n;
}

trieNodeTT* TrieAddNodeTerminal(trieNodeTT* levelT, char nt)
{
	// Add a new terminal node as the last brother of the level's children
	// Return this new node address
	if (levelT == NULL)
	{
		trieNodeTT* n = TrieCreateNodeTerminal(nt);
		return n;
	}
	if (levelT->brother)
	{
		return TrieAddNodeTerminal(levelT->brother, nt);
	}
	trieNodeTT* n = TrieCreateNodeTerminal(nt);
	levelT->brother = n;
	return n;
}

int numberBrotherNode(trieNodeT* n)
{
	// Return the number of little brother(s) of n in the level, including n
	if (n == NULL) { return 0; }
	int i = 1; // we include n in the count
	while (n->brother)
	{
		i++;
		n = n->brother;
	}
	return i;
}

int numberBrotherNodeTerminal(trieNodeTT* n)
{
	// Return the number of little brother(s) of n in the terminal level, including n
	if (n == NULL) { return 0; }
	int i = 1; // we include n in the count
	while (n->brother)
	{
		i++;
		n = n->brother;
	}
	return i;
}


int widthTrie(trie* T)
{
	// Return the width of the trie
	// Remark: the width is the number of nodes in the terminal level
	if (T != NULL)
	{
		if (T->root != NULL)
		{
			return widthSubtrie(T->root, 0, T->readlen);
		}
	}
	return 0;
}

int widthSubtrie(trieNodeT* n, int depth, int readlen)
{
	// Return the width of the subtrie, from the node n
	if (n != NULL)
	{
		int i = 0;
		do
		{
			if (depth + 1 == readlen && n->child)
			{
				i += numberBrotherNodeTerminal((trieNodeTT*)n->child);
			}
			else if (n->child)
			{
				i += widthSubtrie(n->child, depth + 1, readlen);
			}
			n = n->brother;
		} while (n);
		return i;
	}
	return 0;
}

void displayNodeTerminal(trieNodeTT* child)
{
	// Display a terminal node
	if (child != NULL)
	{ 
		fprintf(stdout,"%c\n",child->nt);
		if (child->brother)
		{
			displayNodeTerminal(child->brother);
		}
	}
}

void displayNode(trieNodeT* child, int depth, int readlen)
{
	// Display a node
	if (child)
	{
		fprintf(stdout,"%c\n",child->nt);
		if (depth + 1 == readlen)
		{
			// we want to display a terminal node
			displayNodeTerminal((trieNodeTT*)child->child);
		}
		else if (child->child)
		{
			displayNode(child->child, depth+1, readlen);
		}
		if (child->brother)
		{
			displayNode(child->brother, depth, readlen);
		}
	}
}

trieNodeT* TrieCopyNode(trieNodeT* n)
{
	// Return a copy of the node n
	if (n == NULL) { return NULL; }
	trieNodeT* cp = TrieCreateNode('8');
	cp->child = n->child;
	cp->brother = n->brother;
	cp->nt = n->nt;
	return cp;
}

void TrieRemoveNodeTerminal(trieNodeTT* child)
{
	// Remove a terminal node from the trie
	if (child == NULL) { return ; }
	if (child->brother)
	{
		TrieRemoveNodeTerminal(child->brother);
	}
	if (child->chrs != NULL) { free(child->chrs); }
	if (child->strands != NULL) { free(child->strands); }
	if (child->starts != NULL) { free(child->starts); }
	free(child);
	child = NULL;
}

void TrieRemoveBrotherNodes(trieNodeT* n, int depth, int readlen)
{
	// Remove the brother of the node n, but keep its brothers safe
	if (n == NULL) { return ; }
	if (n->brother)
	{
		//n->brother = TrieRemoveNode(n->brother,depth,readlen);
		if (n->brother->brother)
		{
			TrieRemoveBrotherNodes(n->brother, depth, readlen);
		}
		free(n->brother);
		n->brother = NULL;
	}
}

trieNodeT* TrieRemoveNode(trieNodeT* n, int depth, int readlen)
{
	// Remove the node n, its children and its children brothers
	// Return the n's brother to save or NULL if no brother
	if (n == NULL) { return NULL; }
	trieNodeT* tmp = NULL;
	if (n->child)
	{
		if (depth + 1 == readlen) // we want to remove a terminal node
		{
			//printf("%d is being freed!\n",n->nt);
			TrieRemoveNodeTerminal((trieNodeTT*)n->child); // remove the child
		}
		else
		{
			trieNodeT* output = TrieRemoveNode(n->child, depth+1, readlen); // deletion of the child
			TrieRemoveBrotherNodes(output, depth+1, readlen); // deletion of previous saved brothers
			free(output); // deletion of the output
			output = NULL;
		}
		n->child = NULL;
	}
	if (n->brother)
	{
		tmp = TrieCopyNode(n->brother); // copy of the brother
		tmp->brother = TrieRemoveNode(n->brother, depth, readlen); // deletion of the brother
		// NB: don't froget his saved brother(s)!
		n->brother = NULL;
	}
	//printf("%c is being freed!\n",n->nt);
	free(n);
	n = NULL;
	return tmp;
}

trie* TrieCreate(int readlen)
{
	// create a trie
	trie* trie = malloc(sizeof(struct trie)); // NB: struct trie is important!
	trie->readlen = readlen;
	trie->branches = 0;
	trie->root = TrieCreateNode('0');
	trie->chr_names = NULL; // default is NULL, to be completed with addChr function
	trie->nr_sequences = 0;
	trie->nrchrs = 0;
	return trie;
}

void displayTrie(trie* T)
{
	// display a trie by printing its nucleotides
	if (T == NULL) { return ; }
	fprintf(stdout,"===== TRIE =====\n%c\n",T->root->nt);
	if (T->root->child)
	{
		displayNode(T->root->child,1,T->readlen);
	}
	fprintf(stdout,"===== END =====\n");
}

trie* TrieDestroy(trie* T)
{
	// Free the tree and return a void tree
	int i=0;

	TrieRemoveNode(T->root, 0, T->readlen); // it will detroy the root and
	// recursively all the children and children's brothers.
	while (i<T->nrchrs) // free chromosme records
	{
		free(T->chr_names[i]);
		T->chr_names[i] = NULL;
		i++;
	}
	free(T->chr_names);
	free(T);
	return NULL;
}

int addChr(trie* T, char* chrname, int namelen)
{
	// Add a chromosome in the table chr_names of the trie
	// Return the chromosome number (in order of discovery in fasta file)
	// We check if the chr is stored already
	for (int i=0; i<T->nrchrs; i++)
	{
		if (strcmp(T->chr_names[i],chrname) == 0)
		{
			// fprintf(stdout,"%d\t%s\n", i, T->chr_names[i]);
			return i;
		}
	}
	// it doesn't exist
	if (T->nrchrs == 0) // create the table
	{
		T->chr_names = malloc(sizeof(char**));
	}
	else
	{
		// We increase the size of the table
		T->chr_names = realloc(T->chr_names,sizeof(char**)*(T->nrchrs + 1));
	}
	// copy of the chrname in the table of the trie
	char* seq = malloc(sizeof(char)*(namelen+1));
	seq[namelen] = '\0';
	memcpy(seq, chrname, namelen);
	T->chr_names[T->nrchrs] = seq;
	// fprintf(stdout,"%d\t%s\n", T->nrchrs, T->chr_names[T->nrchrs]);
	T->nrchrs++;
	return T->nrchrs-1;
}

trieNodeT* characterFound(trieNodeT* level, char nt)
{
	// Search for the nucleotide in the brothers of the level
	// Return its location or NULL if it doesn't exist
	if (level == NULL) { return NULL ; }
	if (level->nt == nt) { return level; }
	if (level->brother)
	{
		return characterFound(level->brother,nt);
	}
	return NULL;
}

trieNodeTT* characterFoundTerminal(trieNodeTT* levelT, char nt)
{
	// Search for the nucleotide in the brothers of the terminal level
	// Return its location or NULL if it doesn't exist
	if (levelT == NULL) { return NULL ; }
	if (levelT->nt == nt) { return levelT; }
	if (levelT->brother)
	{
		return characterFoundTerminal(levelT->brother,nt);
	}
	return NULL;
}

int TrieAdd(trie* T, char* array, int alen, int chr, char strand, int start, float score, int* identical_seq)
{
  /* This function is 30% faster compared to TrieAdd. With hg38.fa -l being 30% compared to TrieAdd (i.e. 20m vs 30m). */
  
        // Take as arguments a trie, a sequence, the length of this sequence, the chr number,
	// the strand (+,- or *), the start position, the score (default = 0) and the vcf structure
	// Return 1 if sequence added in the tree, 0 otherwise
	trieNodeT* level = T->root;
	trieNodeTT* levelT = NULL;
	trieNodeT* locationCharacter;
	int i = 0;
	int new_nodes = 0;

	while(i < (alen-1)) {
	  locationCharacter = characterFound(level->child,array[i]);
	  if (locationCharacter == NULL) // the character is not yet in this level
	    {
	      // Add a new standard node and returns its address:
	      //printf("Add the node %c into the trie\n",array[i]);

	      if (level->child == NULL) // if we need to create a child
		{
		  level->child = TrieAddNode(level->child, array[i]);
		  level = level->child; // we move to the next generation
		}
	      else
		{
		  level = TrieAddNode(level->child, array[i]);
		}
	      new_nodes++;
	    }
	  else
	    {
	      level = locationCharacter;
	    }
	  i++;
	}

	// note that in following line, i == (alen-1)
	locationCharacter = (trieNodeT*)characterFoundTerminal((trieNodeTT*)level->child,array[i]);
	if (locationCharacter == NULL) // the character is not yet in this level
	  {
	    
	    // Add a new terminal node and return its address:
	    //printf("Add the node %c into the trie\n",array[i]);
	    if (level->child == NULL) // if we need to create a child
	      {
		level->child = (trieNodeT*) TrieAddNodeTerminal(NULL,array[i]);
		levelT = (trieNodeTT*) level->child;
	      }
	    else
	      {
		levelT = TrieAddNodeTerminal((trieNodeTT*)level->child, array[i]);
	      }
	  }
	else {
	  levelT = (trieNodeTT*)locationCharacter;
	}

	// Lets now deal with additing chr, start and strand info to terminal node
	
	// We check if it already exists (with a difference in the position, due to genomic variations)
	// note that on installation levelT->hits is set to 0, so the check should be fine!
	for (int i=0; i<levelT->hits; i++)
	  {
	    if (levelT->chrs[i] == chr && levelT->strands[i] == strand && levelT->starts[i] == start)
	      {
		*identical_seq = *identical_seq + 1;
		return 0;
	      }
	  }

	// if all slots of levelT->allocated_array_len have been allocated, we will allocate twice the length.
	if (levelT->hits == levelT->allocated_array_len) {
	  levelT->allocated_array_len = levelT->allocated_array_len*2;

	  levelT->chrs = realloc(levelT->chrs, sizeof(int) * (levelT->allocated_array_len));
	  if (levelT->chrs == NULL)
	    {
	      fprintf(stderr,"Failed to allocate memory for %d-th genomic location of gRNA (chrs). Exiting!\n", levelT->hits);
	      exit(1);
	    }
	  levelT->starts = realloc(levelT->starts, sizeof(int) * (levelT->allocated_array_len));
	  if (levelT->starts == NULL)
	    {
	      fprintf(stderr,"Failed to allocate memory for %d-th genomic location of gRNA (starts). Exiting!\n", levelT->hits);
	      exit(1);
	    }
	  levelT->strands = realloc(levelT->strands, sizeof(char) * (levelT->allocated_array_len));
	  if (levelT->strands == NULL)
	    {
	      fprintf(stderr,"Failed to allocate memory for %d-th genomic location of gRNA (strands). Exiting!\n", levelT->hits);
	      exit(1);
	    }
	}
	levelT->chrs[levelT->hits] = chr;
	levelT->starts[levelT->hits] = start;
	levelT->strands[levelT->hits] = strand;
	levelT->score = score;
	// printf("Line %d for seq %s:\n",levelT->hits,array);
	// printf("%d\n",chr);
	// printf("%d\n",start);
	// printf("%f\n", score);
	// printf("%c\n\n",strand);
	levelT->hits++;
	// fprintf(stdout,"%d\t%d\t%c\t%s\t%d\n", chr, start, strand, grna, levelT->hits);
	T->nr_sequences++;
	
	if (new_nodes > 0)
	  {
	    T->branches++;
	  }
	return 1;
}


int TrieMatch(trie *trie, char array[], int alen)
{
	// Return the number of exact matches of the sequence taken as argument in the trie
	// The sequence is in array, its length is given by alen
	trieNodeT *level = trie->root;
	trieNodeTT *levelT = NULL;
	int i = 0;
	int match = 1;

	if(alen != trie->readlen)
	{
    	fprintf(stderr,"Depth of the trie and length of the sequence are not the same!\n");
    	exit(1);
	}

	while(i < alen)
	{
		if (i == (alen-1))
		{
			levelT = characterFoundTerminal((trieNodeTT*)level->child, array[i]);
		}
		else
		{
			level = characterFound(level->child, array[i]);
		}
		if ((i == (alen-1) && levelT == NULL) || (level == NULL)) // we havn't found an exact match
		{
			fprintf(stdout,"%c\n",array[i]);
			match = 0;
			break;
		}
		i++;
	}
	if (match != 0)
	{
		// print found hits:
		return levelT->hits;
	}
	else
	{
		return 0;
	}
}

mcontainer *TrieAMatch(trie *T, char array[], int alen, int nmis)
{
	// Return a container holding information about approximative match of a sequence against the trie
	// The sequence is in array and its length is given by alen
	// nmis is the number of authorized mismatches

	// Define two match containers for holding the approximate match info
	mcontainer *m;
	mcontainer *mn;
	char *str;
	int nelem = 0;
	int pos = 1;
	int i;
	trieNodeT *node;
	trieNodeTT *nodeT;
	int mismatch = 0;

	if(alen != T->readlen)
	{
		fprintf(stderr,"[crisflash] ERROR: Depth of the trie and length of the sequence are not the same!\n");
		exit(1);
	}

	// fprintf(stderr,"Array=%s\n",array);
	
	// install 4 starting strings for the match container
	m = mcontainer_install(4, alen);
	node = T->root->child;
	while (node != NULL)
	{
		str = malloc(sizeof(char)*(alen+1));
		str[alen] = '\0';
		str[0] = node->nt;
		m->sarray[m->pos] = str;
		if (array[0] != node->nt) // there is already a mismatch at the first position
		{
			m->marray[m->pos] = 1; // it's the first mismatch
			m->parray[m->pos] = malloc(sizeof(int));
			m->parray[m->pos][0] = 1; // it's the first position
		}
		else { m->marray[m->pos] = 0; }
		m->nodes[m->pos] = (trieNodeTT*)node;
		m->pos++;
		node = node->brother;
	}
	m->capos = 1;

	// printf("Match container installed. Starting to match moving to the depth of trie one position at the time\n\n");
	while (pos < alen)
	{
	        //fprintf(stdout,"Pos=%d, candidates=%d\n", pos, m->pos);
		mn = mcontainer_install(m->pos * 4, alen); // for each node already registered, we can have 4 new nodes
		nelem = 0;
    
		while (nelem < m->pos)
		{
			if (m->nodes[nelem]) // we start from the stored nodes
			{
			        node = ((trieNodeT*)(m->nodes[nelem]))->child; // the nodes are stored as terminal nodes, here we need the child
				while (node != NULL)
				  {
				    // N.B: The casting into a trieNodeT of a terminal node raises issues on the nucleotide character
				    // But the following code is correct
				    if (pos == alen-1)
				      {
					nodeT = (trieNodeTT*)node;
					while (nodeT != NULL)
					  {					    
					    /*
					    if (array[pos] != nodeT->nt) { mismatch = 1; }
					    else { mismatch = 0; }
					    */
					    mismatch = baseMismatch(nodeT->nt,array[pos]); // returns 1; if mismatch, 0 otherwise.
					    mcontainer_add_nt(m, mn, nelem, pos, nodeT->nt, mismatch, nmis, nodeT);
					    nodeT = nodeT->brother;
					  }
					node = NULL;
				      }
				    else
				      {
					/*
					if (array[pos] != node->nt) { mismatch = 1; }
					else { mismatch = 0; }
					*/
					mismatch = baseMismatch(node->nt,array[pos]); // returns 1; if mismatch, 0 otherwise.
					mcontainer_add_nt(m, mn, nelem, pos, node->nt, mismatch, nmis, (trieNodeTT*)node);
					node = node->brother;
				      }
				  }
			}
			nelem++;
		}
		// fprintf(stdout,"Freeing m and making mn n\n"); fflush(stdout);
		mcontainer_free(m);
		m = mn;

		// mcontainer_print(m);

		pos++;
	}

	m->matches = m->pos;
	mcontainer_score(m, T); // calculate off-target scores and the gRNA (found in the sequence) score

	return m;
}

int TrieAMatchSummary(trie *T, char array[], int alen, int nmis)
{
	mcontainer *m;
	int nrofmatches;

	m = TrieAMatch(T, array, alen, nmis);
	nrofmatches = m->matches;

	// print results
	fprintf(stdout,"[crisflash] Container holds these approximate sequences:\n");
	mcontainer_print(m);

	// free container
	mcontainer_free(m);

	return nrofmatches;
}
