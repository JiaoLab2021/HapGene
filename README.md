# HapGene
An annotation pipeline for haplotype-resolved gene annotation in heterozygous diploid genomes.

## Dependencies
1. [BRAKER3](https://github.com/Gaius-Augustus/BRAKER)
2. [Hisat2](https://github.com/DaehwanKimLab/hisat2)
3. [Samtools](https://github.com/samtools/samtools)
4. [Fastp](https://github.com/OpenGene/fastp)
5. [Miniprot](https://github.com/lh3/miniprot)
6. [IsoSeq](https://github.com/PacificBiosciences/IsoSeq)
7. [Minimap2](https://github.com/lh3/minimap2)
8. [Transdecoder](https://github.com/sghignone/TransDecoder)
9. [EvidenceModeler](https://github.com/EVidenceModeler)
10. [Blast](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/)
11. [LiftOn](https://github.com/Kuanhao-Chao/LiftOn)
12. [MMseq2](https://github.com/soedinglab/MMseqs2)
13. [PASApipeline](https://github.com/PASApipeline/PASApipeline)

## Installation
### Install via github
```bash
git clone https://github.com/JiaoLab2021/HapGene.git
```

## Usage
```bash
HapGene.py [-h] -g genome1 genome2 [genome1 genome2 ...] -w workdir --protein protein -t threads -r rawdatadir -s species -p prefix1 prefix2 [prefix1 prefix2 ...] [--TE_anno] [--long] [--lib TE_library] [--threshold THRESHOLD] [--lencf length_cutoff] [--tpmcf TPM_Value_cutoff]
```
### Input files
1. Reference Genome
2. Homologous Protein Sequences (.fasta)
3. Transcriptome Files
```bash
S1_1.fastq.gz S1_2.fastq.gz S2_1.fastq.gz S2_2.fastq.gz ...
```
Please note that transcriptome files must have the suffixes _1.fastq.gz and _2.fastq.gz.

### Running
Test data are provided in the test/ folder for users to try the pipeline.
Example command:
```bash
HapGene.py -g test/test.hap1.genome.fa test/test.hap2.genome.fa \
-w anno -t 30 -s species \
-r test \
--protein test/homo_protein.fa \
-p hap1 hap2 
```


