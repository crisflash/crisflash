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

#ifndef DEF_VCF
#define DEF_VCF

#ifndef CRISFLASH_VERSION
#define CRISFLASH_VERSION "1.1.0"
#endif

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <libgen.h>

#include "nary_tree.h"

int find_chromosome_positions_ref(char* gName, int* chr_pos, int* chr_size, char** chr_names);
int find_chromosome_positions_vcf(char* vcfName, int* chr_pos, char** chr_names);
int write_variation(VCF* vcf, char* buffer, int length_buffer, int pos, FILE* f_new, FILE* f_new2, int phased, int* var1, int* var2);
int VCF_to_genome(char* gName, char* vcfName, char* newName, char* newName2);

#endif
