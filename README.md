# Project Title

Crisflash is a tool for rapid design and validation of CRISPR/Cas9 guide RNAs. Crisflash is designed for validation speed
and for improved validation accuracy thorugh an option of incorporating user-supplied variant data in addition to
reference genome sequence.

Crisflash is currently supporting only 'NGG' PAM sequence.

## Getting Started

These instructions will get you a copy of crisflash and provide help for running the tool in your local machine.

Crisflash is written in C and dependent on GNU Compiler Collection (gcc) and 'make' for compilation. Crisflash
should compile under most linux distributsions and recent versions of macOS with Xcode installed.

### Installing

Crisflash sourcecode is available from github web page https://github.com/crisflash/crisflash/

To compile the source code, type:

```
cd src
make
```
crisflash binary will be available in bin/crisflash

### How to run

## Example 1.
```
crisflash -g genome.fa -s candidate_sequence.fa -o validated_gRNAs.bed -m 5
```
Validates all gRNAs identified in candidate_sequence.fa for off-targets with up to 5 mismatches genome wide for genome sequence in genome.fa.
Results are saved in validated_gRNAs.bed where gRNA sequence is preserved in 'name' field and off-target validation score in score field.

## Example 2.
```
crisflash -g genome.fa -@ candidate_gRNAs.bed -o validated_gRNAs.bed
```
As Example 1 but here candidate gRNAs are provided in candidate_gRNAs.bed.

## Example 3.
```
crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o validated_gRNAs.bed
```
As Example 1 but here the reference genome sequence in genome.fa is adapted to reflect the variant data provided in phased_variants.vcf.

## Example 4.
```
crisflash -g genome.fa -a genome_gRNAs_dump.bed
```
Lists and saves all gRNAs found in genome.fa to genome_gRNAs_dump.bed. Note that gRNA positions and sequence reported in
output are not being validated for off-targets.

## Crisflash options

-g FILE       FASTA format reference genome.
-b FILE       BED file containing all gRNAs of a genome (the file can be created with -a flag).
-s FILE       FASTA file containing candidate sequence.
-@ FILE       BED file containing candidate gRNA coordinates and sequences on comment field.
-V FILE       phased VCF file.
-o FILE       output saved in BED format. 'sequence/exact match count/off-target count' in comment field, off-target score in score field.
-a FILE       output of all gRNAs in a genome in BED format. No validation of off-targets.
-A FILE       output of all gRNAs in a genome in BED format with off-target validation. This takes long time even in multicore mode.
-m INT        Number of mismatches allowed. Default: 2.
-t INT        Number of threads to use for off-target validation. Default: 1.
-l    	      Include low complexity areas to off-target validation. Default: use only upper case sequences in soft masked genome.
-h    	      Print help.
-v    	      Print version.

## Versioning

Crisflash is currently in release 1.0.0. We aim to follow versioning principles described in https://semver.org/.

## Authors

* **Adrien Jacquin** - *Author for most of the code.*
* **Margus Lukk** - *Idea, initial prototype and current maintainer of the code.*

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://github.com/crisflash/crisflash/blob/master/LICENSE) file for details.

## Acknowledgments

* The work for initial crisflash release (v1.0.0) was funded by Cancer Research UK Cambridge Institute.
