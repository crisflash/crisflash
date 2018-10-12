# Project Title

Crisflash is a tool for rapid design and validation of CRISPR/Cas9 guide RNAs. Crisflash is designed for validation speed
and for improved validation accuracy thorugh an option of incorporating user-supplied variant data in addition to
reference genome sequence.

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

#### Example 1.
```
crisflash -g genome.fa -s candidate_sequence.fa -o validated_gRNAs.bed -m 5
```
Validates all gRNAs identified in candidate_sequence.fa for off-targets with up to 5 mismatches genome wide for genome sequence in genome.fa.
Results are saved in validated_gRNAs.bed where gRNA sequence is preserved in 'name' field and off-target validation score in score field.

#### Example 3.
```
crisflash -g genome.fa -s candidate_sequence.fa -o validated_gRNAs.cas-offinder -C
```
As Example 1 but the output is saved in Cas-OFFinder format. Cas-OFFinder and Crisflash outputs (with option -C) are expected to be identical
after the outputs from both tools have been sorted. E.g. by linux sort: 'sort -k1,1 -k2,2 -k3,3 file > file.sorted'.

#### Example 3.
```
crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o validated_gRNAs.bed
```
As Example 1 but here the reference genome sequence in genome.fa is adapted to reflect the variant data provided in phased_variants.vcf.

#### Example 4.
```
crisflash -g genome.fa -p NNNRYAC -o all_gRNAs.bed
```
Returns all gRNA sequences in genome for PAM sequence 'NNNRYAC'.

## Crisflash options
```
-g FILE       FASTA format reference genome.
-s FILE       FASTA file containing candidate sequence.
-p PAM	      PAM sequence. Default: NGG.
-V FILE       phased VCF file.
-o FILE       output saved in BED format. 'sequence/exact match count/off-target count' in comment field, off-target score in score field.
-B 	      Saves output in BED format, with sequence provided on comment field and off-target score on score field. (Default)
-C            Saves output in Cas-OFFinder format.
-m INT        Number of mismatches allowed. Default: 2.
-t INT        Number of threads to use for off-target validation. Default: 1.
-u 	      Excludes low complexity / soft masked genomic sequences (lowercase) from off-target search.
-h    	      Print help.
-v    	      Print version.
```

## Inputs and outputs

Crisflash inputs for reference genome (option -g) and design target sequence (option -s), are in FASTA format. sgRNAs are designed and matched only against sequences containing standard nucleotide symbols (A,T,G,C,a,t,g,c) are considered. Crisflash memory requirement for human genome with option -u is around ~50GB and ~90GB otherwise.

The tool can be run with either -g or -s option only, resulting in quick identification and print-out of all gRNAs and their location in the sequence with no matching and scoring.

PAM sequence (option -p) may contain any capital IUPAC symbols for nucleotides, hence making the tool universal to any present and future PAM sequences. For example, PAMs such as NGG, NRG, NNNRYAC, NNAGAW are all valid.

Variant information (option -V) is accepted in VCF format with the requirement for the variant file being phased. Internally, phased variants and reference genome are combined in creation of two new improved 'haploid' genomes which are then used for genome-wide scan and indexing of all gRNA sites. In the building / indexing phase, gRNAs of the same locus with identical sequences for both 'haploid genomes' are recorded once. However, should gRNAs of a locus differ in any of the base pairs, two separate gRNAs are being indexed.
 
By default, Crisflash output is provided in BED format and consist of entries to all gRNA candidates identified in the target sequence (provided by option -s). DNA sequences of the candidates are reported on BED file comment field in format 'sequence/nr_of_exact_matches/nr_of_approximate_matches', off-target match scores (ranging 0 to 1, less to more unique) on score field and reported genomic coordinates correspond to the position of the candidate sequence in the design target sequence, not in the reference genome. Candidate gRNAs with no exact match are assigned score 0.

Information about genome-wide matching information of the candidates is available by specifying Crsiflash output in Cas-OFFinder format (option -C). Briefly, the first column in Cas-OFFinder format is for the candidate gRNA sequence with the base pairs part of the PAM being masked by N. The second and third column record the match position of the candidate in the reference genome while the fourth column contains the matched sequence with mismatched bases displayed in lower case. The final two columns are for the chromosome strand and for the number of observed mismatches.

## Versioning

Crisflash is currently in release 1.1.0. We aim to follow versioning principles described in https://semver.org/.

## Authors

* **Adrien Jacquin** - *Author for most of the code.*
* **Margus Lukk** - *Idea, initial prototype and current maintainer of the code.*

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://github.com/crisflash/crisflash/blob/master/LICENSE) file for details.

## Acknowledgments

* The work for initial crisflash release (v1.0.0) was funded by Cancer Research UK Cambridge Institute.
