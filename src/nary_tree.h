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

#ifndef DEF_NARY_TREE
#define DEF_NARY_TREE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>


/*************************/
/*      STRUCTURES       */
/*************************/

typedef struct trieNodeTag
{
	// defines a N-ary tree: each node is linked to a brother and a child
	struct trieNodeTag *brother;
	struct trieNodeTag *child;
	char nt;
} trieNodeT;

typedef struct trieNodeTerminal
{
	struct trieNodeTerminal *brother; // a trieNodeTT can still have a brother
        int allocated_array_len; //
        int hits; // number of hits (nr of elements in chrs, strands and starts)
	int *chrs; // array of chromosome names
	char *strands; // array of strands
	int *starts; // array of start positions
        char *mtypes; // array of variant modifications in gRNA
        char *htypes; // array of source haplotypes. 1 for haplotype 1 (i.e. 1|0 in VCF), 2 for (0|1) and 3 for (1|1) or haplotype level info not available
	char nt;
	float score;
} trieNodeTT;

typedef struct trie
{
	trieNodeT *root;
	int readlen;
	unsigned int branches;
	int nr_sequences;
	char **chr_names; // array of pointers to char arrays holding chr names
	int nrchrs; // number of chromosomes
} trie;

typedef struct mcontainer
{
	// This mismatch container is designed to keep matching info while parsing
	// the trie with a gRNA
	trieNodeTT **nodes; // array of pointers to nodes in trie
	char **sarray; // array of pointers to strings
	int *marray; // array of integers containing mismatch info for each string
	int **parray; // array of integer arrays, for listing the mismatch positions in the sequence
	float *off_scores; // off-target scores table, score is calculated thanks to the positions and the number of mismatches
	float score; // score of the gRNA found in the sequence
	int nelem; // length of arrays (sarray and marray)
	int pos; // positions in use in arrays
	int calen; // length of character arrays (len(sarray[i])) (readlen)
	int capos; // last position used in character array
	int matches; // total nr of matches
} mcontainer;

// Define structure useful for reading VCF file in the same time than reference genome
typedef struct VCF
{
	char* chr;
	int pos;
	int dif; // length difference introduced by the variation (0, positive or negative)
	int ub; // last position in the chromosome to vary
	char* ref;
	char* alt;
	char* hap;
} VCF;

/*************************/
/*      FONCTIONS        */
/*************************/

mcontainer* mcontainer_install(int nelem, int calen);
void mcontainer_free(mcontainer *m);
int mcontainer_add_str(mcontainer *m, int slen);
void mcontainer_add_nt(mcontainer *m, mcontainer *m_new, int mpos, int seqpos, char nt, int mismatch, int max_mismatches, trieNodeTT* nodeT);
void mcontainer_print(mcontainer *m);
void mcontainer_print_pretty(mcontainer *m, int mismatches);
void mcontainer_score(mcontainer *m, trie* T);

VCF* VCF_install(FILE* f_vcf);
void VCF_destroy(VCF* vcf);
int VCF_update(FILE* f_vcf, VCF* vcf);
void VCF_print(VCF* vcf);
VCF* VCF_copy(VCF* vcf);

trieNodeT* TrieCreateNode(char nt);
trieNodeTT* TrieCreateNodeTerminal(char nt);
trieNodeT* TrieAddNode(trieNodeT* level, char nt);
trieNodeTT* TrieAddNodeTerminal(trieNodeTT* levelT, char nt);
int numberBrotherNode(trieNodeT* n);
int numberBrotherNodeTerminal(trieNodeTT* n);
int widthTrie(trie* T);
int widthSubtrie(trieNodeT* n, int depth, int readlen);
void displayNodeTerminal(trieNodeTT* child);
void displayNode(trieNodeT* child, int depth, int readlen);
trieNodeT* TrieCopyNode(trieNodeT* n);
void TrieRemoveNodeTerminal(trieNodeTT* child);
void TrieRemoveBrotherNodes(trieNodeT* n, int depth, int readlen);
trieNodeT* TrieRemoveNode(trieNodeT* n, int depth, int readlen);
trie* TrieCreate(int readlen);
void displayTrie(trie* T);
trie* TrieDestroy(trie* T);
int addChr(trie* T, char* chrname, int namelen);
trieNodeT* characterFound(trieNodeT* level, char nt);
trieNodeTT* characterFoundTerminal(trieNodeTT* levelT, char nt);
int TrieAdd(trie* T, char* array, int alen, int chr, char strand, char mtype, char htype, int start, float score, int* identical_seq);
int TrieMatch(trie *trie, char array[], int alen);
mcontainer *TrieAMatch(trie *T, char array[], int alen, int nmis);
int TrieAMatchSummary(trie *T, char array[], int alen, int nmis);

#endif
