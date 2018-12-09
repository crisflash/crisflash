# Project Title

Crisflash is a tool for rapid design and potential off-target discovery tool of CRISPR/Cas9 guide RNAs. Crisflash is designed for speed and for improved gRNA matching and scoring accuracy through an option of incorporating user-supplied variant data.


## Getting Started

These instructions will get you a copy of Crisflash and provide help for running the tool in your local machine. Before using Crisflash with variant data, users should be familiar with the behaviour of the tool explained further below.

Crisflash is written in C and dependent on GNU Compiler Collection (gcc) and 'make' for compilation. Crisflash should compile under most linux distributions and under recent versions of MacOS having Xcode installed.

### Installing

Crisflash source code is available in Github: https://github.com/crisflash/crisflash/

To compile the source code, type:

```
cd src
make
```
Crisflash binary will appear in bin/crisflash.

### How to run

#### Example 1.
```
crisflash -g genome.fa -s candidate_sequence.fa -o scored_gRNAs.bed -m 5
```
Validates all gRNAs identified in candidate_sequence.fa for off-targets containing up to 5 mismatches for genome provided in genome.fa.
Results are saved in scored_gRNAs.bed where each gRNA sequence is preserved in BAM file 'name' field. Off-target score is in score field.

#### Example 3.
```
crisflash -g genome.fa -s candidate_sequence.fa -o results_gRNAs.cas-offinder -C
```
As Example 1 but the output is saved in Cas-OFFinder format. Cas-OFFinder and Crisflash outputs (with option -C) are expected to be identical
after the outputs from both tools have been sorted. E.g. by linux sort: 'sort -k1,1 -k2,2 -k3,3 file > file.sorted'.

#### Example 3.
```
crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o results_gRNAs.bed -A
```
As Example 1 but here the reference genome sequence in genome.fa is adapted to reflect the variant data provided in phased_variants.vcf and output in Cas-Offinder format has additional column containing information about haplotype the match was found and if the match in genome had PAM, protospacer or both changed by the variant data.

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
-C          Saves output in Cas-OFFinder format.
-A	      Saves output in Cas-OFFinder format with additional column for haplotype level variant info.
-m INT        Number of mismatches allowed. Default: 2.
-t INT        Number of threads to use for off-target validation. Default: 1.
-u 	      Excludes low complexity / soft masked genomic sequences (lowercase) from off-target search.
-h    	      Print help.
-v    	      Print version.
```

## Inputs and outputs

Crisflash inputs for reference genome (option -g) and design target sequence (option -s), are in FASTA format. sgRNAs are designed and matched only against sequences containing standard nucleotide symbols (A,T,G,C,a,t,g,c). Crisflash memory requirement for human genome with option -u is around ~50GB and ~90GB otherwise.

The tool can be run with either -g or -s option only, resulting a quick identification and print-out of all gRNAs and their locations in the sequence with no scoring. Note that this feature can be used in combination with variant data: e.g. providing -g -V options. For such option combination we recommend option -A for output. Sould option -C (default) be chosen, the BED file comment field in the output will contain ?grnaSequence:HM? where H is a character specifying the haplotype the gRNA was found and M variant modification profile for the gRNA (see further below).

PAM sequence (option -p) may contain any capital IUPAC symbols for nucleotides, hence making the tool universal to any present and future PAM sequences. For example, PAMs such as NGG, NRG, NNNRYAC, NNAGAW are all valid.

Variant information (option -V) is accepted in VCF format. Crisflash will proceed if the variant data is not phased. However, we highly recommend phasing the variants first as using unphased variants may result reporting non-existent gRNAs as it is the case with using standard reference genome without variant data.

Crisflash supports both SNPs and INDELS and edits the reference sequence accordingly before using it for genome wide gRNA discovery. When phased variants are provided, Crisflash will create two improved 'haploid genome? sequences and will run gRNA discovery on both. However, gRNAs identical in both haplotypes for a locus are reported only once. Crisflash also records and reports if any bases in gRNA were changed in PAM (mutation type ?b?), protospacer (mutation type ?c?) or both (mutation type reported as ?f?). Applying INDEL data may change the chromosome length and gRNA coordinates; and make comparison of gRNA locations difficult. For this reason, all coordinates reported by Crisflash are for standard reference, even if gRNA location in variant adjusted genome may have been shifted.

Identifying gRNAs from phased variant data may result in following three scenarios:

1. gRNAs in both haplotypes of a locus are identical - we report only one gRNA, specifying the haplotype as 3.

2. gRNAs are different - we report both gRNAs and label each according to the haplotype: haplotype 1 for variants labelled as 1|0, haplotype 2 for variants labelled as 0|1. Note that when candidate matches to the locus, gRNAs of the locus are treated as separate matches, hence reducing the total score.

3. gRNA is present only in one haplotype - this may be due to deletion or creation of PAM in one of the haplotypes or that the gRNA in the other haplotype has more mutations in it than specified by option -m. We report only the gRNA for the haplotype.

We suggest users to look out for candidates which match gRNAS with no mismatches or a small number of mismatches where the haplotype value is not 3 as it indicates heterozygosity (Crisflash output from option -A). We leave it for the user to interpret the effect it may have to their experiment. Note that Crisflash output from option -A may need to be sorted by chromosomes and chromosomal coordinates as hits from both haplotypes may not appear one after another in the output. Output from option -A could be sorted as follows:

'sort -k1,1 -k2,2 -k3,3 file > file.sorted'

In building / indexing phase, gRNAs of the same locus with identical sequences for both 'haploid genomes' are recorded once. However, should gRNAs of a locus differ in any of the base pairs, two separate gRNAs are being indexed.
 
By default, Crisflash output is provided in BED format and consist of entries to all gRNA candidates identified in the target sequence (provided by option -s). DNA sequences of the candidates are reported on BED file comment field('sequence/nr_of_exact_matches/nr_of_approximate_matches'), off-target match score (ranges from 0 to 1, smaller score indicates higher uniqueness) is on score field and reported genomic coordinates correspond to the position of the candidate sequence in the design target sequence, not in the reference genome. Candidate gRNAs with no exact match are assigned score 0.

Information about genome-wide matching information of the candidates is available by specifying Crsiflash output in Cas-OFFinder format (options -C and -A). Briefly, the first column in Cas-OFFinder format is for the candidate gRNA sequence with the base pairs part of the PAM being masked by N. The second and third column record the match position of the candidate in the reference genome while the fourth column contains the matched sequence with mismatched bases displayed in lower case. The final two columns are for the chromosome strand and for the number of observed mismatches. Option -A adds extra column, not present in Cas-Offinder output, containing two-character string where the first character specifies the source haplotype of the gRNA (values 1,2 or 3) and the second shows whether the sequence for PAM area (value ?b?), protospacer (value ?c?) or both (value ?f?) were changed by the variant data. 

## Versioning

Crisflash is currently in release 1.2.0. We aim to follow versioning principles described in https://semver.org/.

## Authors

* **Adrien Jacquin** - *Author*
* **Margus Lukk** - *Author and current maintainer*

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://github.com/crisflash/crisflash/blob/master/LICENSE) file for details.

## Acknowledgments

* The work for initial crisflash release (v1.0.0) was funded by Cancer Research UK Cambridge Institute.

