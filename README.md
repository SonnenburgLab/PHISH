# PHISH -- PHages In SearcHsra  
Written by Bryan D. Merrill  

## Introduction:  
This file describes how to use SearchSRA (https://www.searchsra.org/) to map 100,000 reads from each of >100,000 different metagenomes onto a sequence of interest to determine where the query can be found in the environment. SearchSRA runs in two modes: nucleotide mode (Search-SRA), and protein mode (Translated-DNA-Search). This document will focus primarily on protein mode (Translated-DNA-Search).

## Software prerequisites:  
* Internet browser  
* Bedtools (https://bedtools.readthedocs.io/en/latest/)  

## Running SearchSRA in protein mode
In order to run SearchSRA in protein mode (as of August 2019), a FASTA-formatted file containing only a single amino acid sequence needed to be uploaded. Since most sequences of interest (e.g. a whole phage genome) contain many individual protein sequences, they all need to be combined into a single amino acid sequence in order to meet the input requirements.  

We can use bash commands to do this. For example, the contents of `amino_acid_seqs.faa` might look like this, where each sequence occupies two lines, (1) a header, and (2) the amino acid sequence.  
```
>GenomeName_gp1
MIITTEKETILGNGSKSKAFSITASPKVFKILSSDLYTNKIRAV
>GenomeName_gp2
MKSYKVNLELFDKAVHREYRIIQRFFDMGEAEEFKTRFKDIRDK
>GenomeName_gp3
MKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADH
```

The amino acid sequence can't be spread across multiple lines or the folloiwng command won't work. To concatenate all amino acid sequences, separated by "XXX", you can run the following bash command:  
```
echo ">GenomeName" > concat_amino_acid_seqs.faa;
grep -v ">" amino_acid_seqs.faa | sed 's/$/XXX/' | tr -d '\n' >> concat_amino_acid_seqs.faa
echo "" >> concat_amino_acid_seqs.faa
```

Which writes the following text to `concat_amino_acid_seqs.faa`:  
```
>GenomeName
MIITTEKETILGNGSKSKAFSITASPKVFKILSSDLYTNKIRAVXXXMKSYKVNLELFDKAVHREYRIIQRFFDMGEAEEFKTRFKDIRDKXXXMKFVKIDSSSVDMKKYKLQNNVRRSIKSSSMNYANVAIMTDADHXXX
```
After you make an account with SearchSRA, upload `concat_amino_acid_seqs.faa` and run in Translated-DNA-Search mode.  

There are several caveats about this approach, given that short sequencing reads will be translated and aligned against this single concatenated amino acid sequence using the `diamond blastx` algorithm. First, not all genes in the sequence are separated by 3 amino acids (9 bp, represented by "XXX"). You can change the number of "X" or skip adding any "X" characters by eliminating the `sed` step above.  Second, although most genes on bacteriophage genomes tend to be facing the same direction, concatenated amino acid sequences do not adequately deal with sequence complementarity. It is a reasonable assumption that some reads will map on the edge of a protein sequence whose gene is reverse-complemented relative to a subsequent gene (3') and will be forced to pick one or the other to align to. This may have small but measurable impacts on coverage estimations later in the pipeline, but cannot be avoided due to the current (August 2019) implementation of protien (Translated-DNA-Search) SearchSRA.  

## Parsing SearchSRA results
When it is finished, download your results from SearchSRA. They will come down as a zip file. Unzip it, and you'll see it contains the following directory structure:  
```
results/
 -> 1/
   -> DRR000836.m8
   -> DRR000980.m8
   ...
   -> Up to 5,000 m8 files.
 -> 2/
 ...
 -> 44/
```

Each folder contains up to 5,000 M8-formatted `diamond blastx` result files, where each file represents one SRR or ERR used to map reads against your concatenated amino acid sequence query.  The M8 (also `blast outfmt=6`) format specification is here (http://www.metagenomics.wiki/tools/blast/blastn-output-format-6).  

However, many of these files are empty, as one file is written for each metagenome queried regardless of whether there were any hits. To make things easier, let's just remove all empty files!  
```
for folder in {1..44}; do rm `find $folder -type f -size 0`; done
```

It is important to note now (and later) that an SRR/ERR represents a single run of a DNA sequencer on a single sample (described by an SRS identifier). It is possible to sequence one sample (identified by one SRS) several times (identified by several SRRs) and therefore each SRR file in the results does not necessarily represent one unique sample. Same-sample SRR data ought to be combined if they came from the same sample (SRS) to represent the total sequencing data for that sample.  

Columns of interest in each M8 file include:  
- #2 (query name)  
- #8 (the first bp the read aligned to)  
- #9 (the last bp the read aligned to)  

If we extract just those 3 columns from each M8 file, we get a very simple BED file (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). We can do this by running the following on each of the directories, numbered 1 to 44. This will place a `*.bed` file next to each M8 file:  
```
for num in {1..44}; do 
for file in `ls $num/*`; do 
out=${file%%.*}.bed; 
echo $out; 
cut -f 2,8,9 $file > $out; 
done; 
done;
```

Next, we make a simple BED file for our genome. You can call it `GenomeName.bed`. It looks like this:  
```
GenomeName	0	10022
```
Where 10022 is the number of amino acid residues (including the "X" characters) are in your concatenated sequence you uploaded (which is 10023) - 1 = 100222. This is the easiest way to count:  
```
grep -v ">" concat_amino_acid_seqs.faa | tr -d '\n' | wc -m
```

Now, we can use bedtools to estimate how many reads mapped to each position of our genome, with single-base resolution. Unfortunately, this results in some pretty huge files. At least we can write a bash `for` loop to do this for us:  
```
printf "SRR\treads_mapped\tpos_detected\ttotal_length\tperc_detected\n" > GenomeName_results/ALL_bp_stats.txt
touch GenomeName_results/ALL_bp_COVERAGE.txt
rm GenomeName_results/ALL_bp_COVERAGE.txt
for folder in {1..44}; do 
	for file in GenomeName_results/$folder/*.bed; do
		name=`basename $file .bed`
		echo GenomeName_results/$folder/$name
		cov=`bedtools coverage -a GenomeName.full.bed -b $file -d | cut -f 5 | tr "\n" "\t"`
		printf "$name\t$cov\n" >> GenomeName_results/ALL_bp_COVERAGE.txt
		stats=`bedtools coverage -a GenomeName.full.bed -b $file | cut -f 4,5,6,7`
		printf "$name\t$stats\n" >> GenomeName_results/ALL_bp_stats.txt
	done
done
```

You can now use `ALL_bp_COVERAGE.txt` for making a heatmap directly, or use `ALL_bp_stats.txt` to browse summary statistics for each SRR/ERR mapped onto your query.  

File formats:
- `ALL_bp_COVERAGE.txt` = this table has no headers. The first column is the `SRR` or `ERR` accession number, and columns 2 through the end each represent one base pair positoin in the concatenated amino sequence you provided to SearchSRA.  
- `ALL_bp_stats.txt` = table headers are SRR, reads_mapped, pos_detected, total_length, and perc_detected.  

## Pre-processing the tables  
You can use the `SRAdb` R package to gather sample names and metadata for each SRR/ERR you are analyzing, and even combine rows where more than one SRR/ERR come from the same sample, ending with one row representing that sample.  

You're now ready to make a heatmap of your data!  
