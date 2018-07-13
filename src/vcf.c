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
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <ctype.h>

#include "vcf.h"
#include "nary_tree.h"
#include "read.h"

int find_chromosome_positions_ref(char* gName, int* chr_pos, int* chr_size, char** chr_names)
{
	// Parse a first time the reference genome and store the chromosome positions, sizes and names
	// in the tables chr_pos, chr_size and chr_names
	// Return pT, the current position in the tables, or the size of the tables
        FILE* f = fopen(gName,"r");
	if (!f)
	  {
	    fprintf(stderr,"[crisflash] Failed to read %s. Exiting!\n", gName);
	    exit(1);
	  }
	
	int pT = 0;
	int previous_pos = 0;
	int size = 0;
	char* line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);

	long elapsed;
	struct timeval t0, t1;
	gettimeofday(&t0, 0);

	fprintf(stdout,"[crisflash] Analysing %s ...\n", gName);
	
	while (read != -1)
	{
		if (line[0] == '>') // new chromosome
		{
			if (pT > 0) { chr_size[pT-1] = size; }
			size = 0;
			chr_pos[pT] = previous_pos;
			chr_names[pT] = malloc(sizeof(char)*256);
			int i = 1;
			while (line[i] != ' ' && line[i] != '\n')
			{
				chr_names[pT][i-1] = line[i];
				i++;
			}
			chr_names[pT][i-1] = '\0';
			pT++;
			if (pT > CHR_NUMBER_MAX) {
			  fprintf(stderr,"[crisflash] Max number chromosomes reached (%d)! (See CHR_NUMBER_MAX constant in read.h). Exiting!\n\n",CHR_NUMBER_MAX);
			  exit(1);
			}
		}
		previous_pos = ftell(f);
		size += sizeof(char)*strlen(line);
		read = getline(&line, &len, f);
	}
	chr_size[pT-1] = size;
	free(line);

	gettimeofday(&t1, 0);
	elapsed = (t1.tv_sec-t0.tv_sec);
	if(elapsed < 1) { elapsed = 1; }
	fprintf(stderr,"[crisflash] %d chromosomes and chromosomal fragments identified (%ld seconds).\n", pT, elapsed);
	
	return pT;
}

int find_chromosome_positions_vcf(char* vcfName, int* chr_pos, char** chr_names)
{
	// Parse a first time the VCF file and store the chromosome positions and names
	// in the tables chr_pos and chr_names
	// Return pT, the current position in the tables, or the size of the tables
	FILE* f = fopen(vcfName,"r");
	int pT = 0;
	int previous_pos = 0;
	char* line = NULL;
	size_t len = 0;
	ssize_t read = getline(&line, &len, f);
	while (read != -1)
	{
		if (line[0] != '#')
		{
			char* name = malloc(sizeof(char)*256);
			int i = 0;
			while (line[i] != '\t')
			{
				name[i] = line[i];
				i++;
			}
			name[i] = '\0';
			if (pT == 0)
			{
				chr_names[pT] = malloc(sizeof(char)*256);
				strcpy(chr_names[pT],name);
				chr_pos[pT] = previous_pos;
				pT++;
			}
			else
			{
				if (strcmp(name,chr_names[pT-1]) != 0)
				{
					chr_names[pT] = malloc(sizeof(char)*256);
					strcpy(chr_names[pT],name);
					chr_pos[pT] = previous_pos;
					pT++;
				}
			}
			free(name);
			name = NULL;
		}
		previous_pos = ftell(f);
		read = getline(&line, &len, f);
	}
	free(line);
	return pT;
}

int write_variation(VCF* vcf, char* buffer, int length_buffer, int pos, FILE* f_new, FILE* f_new2, int phased, int* var1, int* var2)
{
	// Write the variation in the new fasta files, taking the haplotypes into account
	// Return 1 if the variation was written, 0 otherwise
	if (strcmp(vcf->hap,"0|0") == 0 || strcmp(vcf->hap,"0/0") == 0) { return 0; }
	if (vcf->alt[0] == '<') { return 0; } // there is no recorded alternative nucleotide
	int i = 0;
	int j = 0;
	while (j < strlen(vcf->ref)) // we check if every letters in vcf->ref match the reference genome
	{
		if (pos + i == length_buffer) { return 0; }
                if (buffer[pos+i] == '\n') { i++; continue; }
                if (vcf->ref[j] != toupper(buffer[pos+i])) { return 0; }
		i++;
		j++;
	}
	if (phased == 1)
	{
		if (strcmp(vcf->hap,"1|0") == 0) // variation in haplotype 1
		{
			if (vcf->dif >= 0) // insertion or SNP
			{
				for (int i=0; i<strlen(vcf->alt); i++)
				{
					if (i >= strlen(vcf->ref)) { fputc('+',f_new); }
					fputc(vcf->alt[i],f_new);
				}
				*var1 += 1;
			}
			else // deletion
			{
				for (int i=0; i<strlen(vcf->ref); i++)
				{
					if (i < strlen(vcf->alt)) { fputc(vcf->alt[i],f_new); }
					else { fputc('-',f_new); }
				}
				*var1 += 1;
			}
			for (int i=0; i<strlen(vcf->ref); i++) { fputc(vcf->ref[i],f_new2); } // no variation in haplotype 2
			return 1;
		}
		if (strcmp(vcf->hap,"0|1") == 0) // variation in haplotype 2
		{
			if (vcf->dif >= 0) // insertion or SNP
                        {
                                for (int i=0; i<strlen(vcf->alt); i++)
                                {
                                        if (i >= strlen(vcf->ref)) { fputc('+',f_new2); }
                                        fputc(vcf->alt[i],f_new2);
                                }
                                *var2 += 1;
                        }
                        else // deletion
                        {
                                for (int i=0; i<strlen(vcf->ref); i++)
                                {
                                        if (i < strlen(vcf->alt)) { fputc(vcf->alt[i],f_new2); }
                                        else { fputc('-',f_new2); }
                                }
                                *var2 += 1;
                        }
                        for (int i=0; i<strlen(vcf->ref); i++) { fputc(vcf->ref[i],f_new); } // no variation in haplotype 1
                        return 1;
		}
	}
	// variation in both haplotypes
	if (vcf->dif >= 0) // insertion or SNP
	{
		for (int i=0; i<strlen(vcf->alt); i++)
		{
			if (i >= strlen(vcf->ref)) { fputc('+',f_new); fputc('+',f_new2); }
			fputc(vcf->alt[i],f_new);
			fputc(vcf->alt[i],f_new2);
		}
		*var1 += 1;
		*var2 += 1;
	}
	else // deletion
	{
		for (int i=0; i<strlen(vcf->ref); i++)
		{
			if (i < strlen(vcf->alt)) { fputc(vcf->alt[i],f_new); fputc(vcf->alt[i],f_new2); }
			else { fputc('-',f_new); fputc('-',f_new2); }
		}
		*var1 += 1;
		*var2 += 1;
	}
	return 1;
}

int VCF_to_genome(char* gName, char* vcfName, char* newName, char* newName2)
{
	// Parse a reference genome chromosome by chromosome in the same time than a VCF file
	// Write the individual genome in two new files (phased genome ==> one haplotype per file)
	// Take into account SNPs and indels
	// Return 1 if phased, 0 otherwise
	double start_time = now();
	double end_time;
	double total_time;

	FILE* f_ref = fopen(gName,"r");
	FILE* f_vcf = fopen(vcfName,"r");
	FILE* f_new = fopen(newName,"w");
	FILE* f_new2 = fopen(newName2,"w");

	int phased;
	/* Install the VCF structure */
	VCF* vcf = VCF_install(f_vcf); // create vcf and store first line
	/* Find out if genome is phased or not */
	if (vcf->hap[1] == '|') { phased = 1; } // genome is phased
	else { phased = 0; } // genome is not phased

	int end_vcf = 0;
	int variation;
	int nb_of_var_1 = 0;
	int nb_of_var_2 = 0;

	int chr_pos[CHR_NUMBER_MAX]; // the list of positions of the chromosomes in the reference genome file
	int chr_size[CHR_NUMBER_MAX]; // the list of size of the chromosomes in the reference genome file
	char* chr_names[CHR_NUMBER_MAX]; // the list of chromosome names in the reference genome file
	int n = find_chromosome_positions_ref(gName,chr_pos,chr_size,chr_names); // n is the tables size
	/* Same for the VCF file */
	int chr_pos_vcf[CHR_NUMBER_MAX]; // the list of positions of the chromosomes in the VCF file
	char* chr_names_vcf[CHR_NUMBER_MAX]; // the list of chromosome names in the VCF file
	int m = find_chromosome_positions_vcf(vcfName,chr_pos_vcf,chr_names_vcf); // m is the tables size
	
	for (int i = 0; i < n; i++)
	{
		fseek(f_ref, chr_pos[i], SEEK_SET);
		char* buffer = malloc(sizeof(char)*(chr_size[i]+1));
		buffer[chr_size[i]] = '\0';
		int read = fread(buffer,chr_size[i],1,f_ref);
		if (read == -1) { fprintf(stderr,"ERROR reading reference genome\n"); exit(1); }

		int x = 1;
		char* chr = malloc(sizeof(char)*256);
		while (buffer[x] != ' ' && buffer[x] != '\n') // we know we are at the beginning of a chromosome
		{
			chr[x-1] = buffer[x]; // we store chromosome name without > and \n
			x++;
		}
		chr[x-1] = '\0';
		fprintf(stdout,"[crisflash] Processing chromosome %s ...\n", chr);
		// Initializing position in buffer and chromosome
		int pos = 0;
		int chr_pos = 1;
		int j = 0;
		// we are looking for the chromosome in the VCF file
		while (j<m && strcmp(chr_names_vcf[j],chr_names[i]) != 0) { j++; }
		if (j==m)
		{
		        fprintf(stdout,"[crisflash] WARNING: No variants for chromosome %s!\n", chr);
			int write = fwrite(buffer,chr_size[i],1,f_new);
			if (write == -1) { fprintf(stderr,"[crisflash] ERROR: Writing to %s failed. Exiting!\n", newName2); exit(1); }
			if (phased == 1) { write = fwrite(buffer,chr_size[i],1,f_new2); }
			if (write == -1) { fprintf(stderr,"[crisflash] ERROR: Writing to %s failed. Exiting!\n", newName2); exit(1); }
		}
		else
		{
			// Write the first line
			while(buffer[pos] != '\n')
			{
				fputc(buffer[pos],f_new);
				if (phased == 1) { fputc(buffer[pos],f_new2); }
				pos++;
			}
			fputc(buffer[pos],f_new);
			if (phased == 1) { fputc(buffer[pos],f_new2); }
			pos++;
			// Find the first variation line in VCF file, for this chromosome
			fseek(f_vcf, chr_pos_vcf[j], SEEK_SET);
			end_vcf = VCF_update(f_vcf,vcf);
			while (end_vcf == 0 && strcmp(vcf->chr,chr) == 0 && pos < chr_size[i])
			{
				while (pos < chr_size[i] && (chr_pos < vcf->pos))
				{
					// if (chr_pos == vcf->pos-1) { printf("%c %c %c %c\n",buffer[pos], buffer[pos+1], buffer[pos+2], buffer[pos+3]); }
					fputc(buffer[pos],f_new);
					if (phased == 1) { fputc(buffer[pos],f_new2); }
					pos++;
					if (buffer[pos] == '\n') { continue; } // we don't increment chr_pos if carriage return
					chr_pos++;
				}
				if (pos == chr_size[i]) { continue; } // we will exit the while loop
				// Here, there is a variation
				// printf("%d  %d\n",chr_pos, vcf->pos);
				// if (chr_pos != vcf->pos) { printf("ERROR: There is a lag between reference genome chromosome position and individual genome\n"); exit(1); }
				variation = write_variation(vcf,buffer,chr_size[i],pos,f_new,f_new2,phased,&nb_of_var_1,&nb_of_var_2);
				// printf("%d\n",variation);
				if (variation == 1)
				{
					char* old_hap = malloc(sizeof(char)*(strlen(vcf->hap)+1));
					old_hap[strlen(vcf->hap)] = '\0';
					strcpy(old_hap,vcf->hap); // keep record of the old haplotype before updating VCF line
					int min_pos = chr_pos; // minimum position after which we can have variations
					if (vcf->dif <= 0) // in case of deletion/SNP, we have to make sure the next var does not concern deleted/modified nucleotides
					{
						min_pos = chr_pos + strlen(vcf->ref);
					}
					int k = 0;
					while (k<strlen(vcf->ref))
					{
						pos++;
						k++;
						if (buffer[pos] == '\n') { pos++; }
						chr_pos++;
					}
					end_vcf = VCF_update(f_vcf,vcf);
					while (vcf->pos < min_pos && end_vcf == 0 && \
((vcf->hap[0] == 1 && old_hap[0] == 1) || (vcf->hap[2] == 1 && old_hap[2] == 1) || phased == 0))
					// we don't want to overlap the variations, but it must be possible if it's not on the same haplotype
					{
						end_vcf = VCF_update(f_vcf,vcf);
					}
					free(old_hap);
					old_hap = NULL;
				}
				else // the variation wasn't applied ==> we have to write the reference genome instead
				{
					fputc(buffer[pos],f_new);
					if (phased == 1) { fputc(buffer[pos],f_new2); }
					pos++;
					if (buffer[pos] == '\n')
					{
						fputc(buffer[pos],f_new);
                                        	if (phased == 1) { fputc(buffer[pos],f_new2); }
						pos++;
					}
					chr_pos++;
				}
				while (vcf->pos < chr_pos && end_vcf == 0) // it means we only apply one variation in the same position
				{
					end_vcf = VCF_update(f_vcf,vcf);
				}
			}
			while (pos < chr_size[i]) // there is no more variation for this chromosome, we write everything left
			{
				fputc(buffer[pos],f_new);
				if (phased == 1) { fputc(buffer[pos],f_new2); }
				pos++;
			}
		}
		free(buffer);
		free(chr);
		buffer = NULL;
		chr=NULL;
	}
	fclose(f_ref);
	fclose(f_vcf);
	fclose(f_new);
	fclose(f_new2);
	if (phased == 0) { remove(newName2); }

	end_time = now();
	total_time = end_time - start_time;
	if (phased == 1) {fprintf(stdout,"[crisflash] Finished creating phased genomes %s (%d variants) and %s (%d variants). (%f seconds)\n", newName, nb_of_var_1, newName2, nb_of_var_2, total_time); }
	else {fprintf(stdout,"[crisflash] A modified genome incorporating variant data (%d variants) has been saved in %s. (%f seconds)\n", nb_of_var_1, newName, total_time);}

	for (int i=0; i<n; i++) { free(chr_names[i]); }
	for (int i=0; i<m; i++) { free(chr_names_vcf[i]); }
	VCF_destroy(vcf);

	return phased;
}

