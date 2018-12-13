# Project Title

Crisflash is a tool for rapid design and a potential off-target discovery tool for CRISPR/Cas9 guide RNAs. Crisflash is designed for speed, improved gRNA matching and scoring accuracy by providing the option to incorporate user-supplied variant data.


## Getting Started

This step-by-step guide will provide you with a copy of Crisflash and give instructions for running the tool in your local machine. Before using Crisflash with variant data, users should be familiar with the behaviour of the tool explained further below.

Crisflash is written in C and is dependent on GNU Compiler Collection (gcc) and ‘'make'’ for compilation. Crisflash should compile under most linux distributions as well as under recent versions of MacOS having Xcode installed.

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
after the outputs from both tools have been sorted. E.g. by linux sort: ‘'sort -k1,1 -k2,2 -k3,3 file > file.sorted'’.

#### Example 3.
```
crisflash -g genome.fa -V phased_variants.vcf -s candidate_area.fa -o results_gRNAs.bed -A
```
As Example 1 but here the reference genome sequence in genome.fa is adapted to reflect the variant data provided in phased_variants.vcf. The output is in modified Cas-Offinder format described further below.


.

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

## Crisflash inputs

Crisflash inputs for reference genome (option -g) and design target sequence (option -s), are expected to be in FASTA format (file name ending with .fa). sgRNAs are designed and matched only against sequences containing standard nucleotide symbols (A,T,G,C,a,t,g,c). Crisflash memory requirement for human genome with option -u is around ~50GB and ~90GB otherwise.

Running the tool with either -g or -s option, opposed to using both, results in the print-out of all gRNAs with no matching and scoring.

PAM sequence (option -p) may contain any capital IUPAC symbols for nucleotides, consequently making the tool universal to any present and future PAM sequences. For example, PAMs such as NGG, NRG, NNNRYAC, NNAGAW are all valid.

Variant information (option -V) is accepted in VCF format (.vcf suffix expected). Crisflash will proceed if the variant data is not phased. However, we highly recommend phasing the variants first. Using Crisflash with unphased variants may produce incorrect results.

Crisflash supports both SNPs and INDELs and edits the reference sequence accordingly before being used for gRNA discovery. When phased variants are provided, Crisflash will create two improved ‘’haploid genome’’ sequences and will run gRNA discovery for both. Identified gRNAs are marked for their haplotype with possible values being ‘’1’’, ‘‘2’’ or ‘‘3’’ if gRNA is identical in both haplotypes. Crisflash also records and will report when genomic variation has resulted in any changes to the PAM area (mutation type ‘b’),in protospacer (mutation type ‘’c’’) or in both (mutation type ‘’f’’). In the case of no change compared to reference, the mutation type is set to ‘’0’’. 

Applying INDEL data may change the chromosome length and gRNA coordinates; making comparison of gRNA locations difficult. For this reason, all coordinates reported by Crisflash are provided for original reference, even if the gRNA location in variant-adjusted genome may be shifted due to incorporation of INDEL data.

Identification of gRNAs in phased data may result in following three scenarios:

1. gRNAs are identical in both haplotypes – Crisflash reports gRNA with haplotype value ‘‘3’’.

2. gRNAs for the locus are different - we report both gRNAs and label each according to its source haplotype: ‘’1’’ for variants labelled 1|0 in vcf; and ‘’2’’ for variants labelled as 0|1. Note: genome-wide search for target and off-target matches for candidates is based on sequence similarity for up to a specified number of mismatches. If one of the gRNAs of the locus (for one of the haplotypes) has more than an expected number of mismatches (specified by option -m), it will not appear in the search results and will not be contributing for scoring. On the other hand, if both gRNAs are matched, the mismatch score will be affected and will be consequently lower. In result, match against the locus with gRNAs differing in haplotypes but those not being too different are penalised by the scoring system.

3. gRNA is present only in one haplotype - this may be due to the deletion or appearance of a new PAM compared with a reference genome in one of the haplotypes.

We suggest careful examination of results for candidates with exact or close to exact matches to only one haplotype.

## Crisflash outputs

### Output in BED file format

By default, Crisflash output is provided in BED file format. In the BED file, the content of the comment filed varies. This depends on whether the program was executed by -g or -s option, or using both. For either -g or -s, the comment filed contains the gRNA sequence followed by ‘’:’’ and two characters. The first of the two characters is a number (1-3) indicating the haplotype while the other  show the mutation type (values ‘’b’’,’’c’’, ‘‘f’’, ‘‘0’’).

Executions with -g and -s options result comment field with following information separated by ‘/’: gRNA sequence, the number of of exact matches and the number of approximate matches.

The off-target match score, ranging from 0 to 1, is reported on BED score field. Candidate gRNAs with higher uniqueness have a score closer to 0. Genomic coordinates reported in BED are for the start and end positions of the gRNA in a sequence provided by the option -s. Candidate gRNAs with no exact match have score 0.

### Output in Cas-OFFinder format

Detailed information about all matching reference sites is provided in Cas-OFFinder format by specifying options -C and -A. Briefly, the first column is for candidate gRNA sequence where the section of PAM sequence is indicated by a sequence of ‘’N’’. The second and third column records candidate-match positions in reference while the fourth column contains matched-sequence in reference genome with mismatched bases being in lower case. The final two columns are for the chromosome strand and for the number of observed mismatches. Option -A adds an extra column not present in Cas-Offinder tool output. This is for the two-character string, similar to the one already discussed for BED file output: the first of the two characters show the haplotypic origin of the gRNA (values ‘’1’’ or ‘’2’’; or ‘’3’’ if gRNA is identical in both haplotypes); and the second character indicates whether the sequence in PAM area (value ‘'b’'), protospacer area (value '’c'’), both (value ‘'f'’) or neither (value ‘'0’'), was changed by applied variant data.

Crisflash output from option -A needs to be sorted by chromosomes and chromosomal coordinates in order for the gRNAs for both haplotypes (if present and identified by the match) appear on consecutive lines. An example sort command is shown below:
‘’’
'sort -k1,1 -k2,2 -k3,3 file > file.sorted'
‘’’
## Versioning

Crisflash is currently in release 1.2.0. We aim to follow versioning principles described in https://semver.org/.

## Authors

* **Adrien Jacquin** - *Author*
* **Margus Lukk** - *Author and current maintainer of the code*

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](https://github.com/crisflash/crisflash/blob/master/LICENSE) file for details.

## Acknowledgments

* The work for initial crisflash release (v1.0.0) was funded by Cancer Research UK Cambridge Institute.

