# HapGene
An annotation pipeline for haplotype-resolved gene annotation in heterozygous diploid genomes.

## Dependencies
1. [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)

## Installation
### Install via github
```bash
git clone https://github.com/JiaoLab2021/HapGene.git
```

## Usage
```bash
HapGene.py [-h] -g genome1 genome2 [genome1 genome2 ...] -w workdir --protein protein -t threads -r rawdatadir -s species -p prefix1 prefix2 [prefix1 prefix2 ...] [--TE_anno] [--long] [--lib TE_library] [--threshold THRESHOLD] [--lencf length_cutoff] [--tpmcf TPM_Value_cutoff]
```


