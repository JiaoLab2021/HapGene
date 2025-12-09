#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2025/7/24 11:02
# @Author  : jhuang
# 原始fastq文件必须更名为*_1.fastq.gz/ *_2.fastq.gz
import glob
import subprocess
import argparse
import os
import re
import multiprocessing
import sys
import time
import logging
from functools import partial


# 先激活环境braker3
def get_parser():
    parser = argparse.ArgumentParser(description="This script parallelizes the commands of the input file.")

    required = parser.add_argument_group('required arguments')
    required.add_argument("-g", "--genome", metavar="genome1 genome2", help="Input your genome", required=True, nargs='+')
    required.add_argument("-w", "--workdir", metavar="workdir", help="Input your workdir.", required=True)
    required.add_argument("--protein", metavar="protein", help="Input your protein sequence.", required=True)
    required.add_argument("-t", "--threads", metavar="threads", help="Input your threads.", required=True)
    required.add_argument("-r", "--rawdatadir", metavar="rawdatadir", help="Input rawdata dir.", required=True)
    required.add_argument("-s", "--species", metavar="species", help="Input species name.", required=True)

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--TE_anno', action='store_true', help='Run TE annotation step', required=False)
    optional.add_argument('--long', action='store_true', help='Pacbio long-reads mode', required=False)
    optional.add_argument('--lib', metavar="TE_library", help='Library file for RepeatMasker', required=False)
    required.add_argument("-p", "--prefix", metavar="prefix1 prefix2", help="Output prefix.", required=True, nargs='+')

    return parser


def cmd_linux(cmd):
    logger = setup_logger()
    logger.info(cmd)
    subprocess.run(cmd, shell=True, close_fds=True)


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        print("new path: {}".format(path))
    else:
        print("The path is exit !!!")


def setup_logger(log_filename=None):
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    if not logger.handlers:
        if not log_filename:
            log_filename = f"{os.path.splitext(os.path.basename(__file__))[0]}.log"
        fh = logging.FileHandler(log_filename, encoding='utf-8')
        sh = logging.StreamHandler()

        standard_format = logging.Formatter(
            "[%(asctime)s] [%(filename)s line:%(lineno)d][%(levelname)s]\t%(message)s")
        simple_format = logging.Formatter(
            "[%(levelname)s][%(asctime)s][%(filename)s:%(lineno)d] %(message)s")

        fh.setFormatter(standard_format)
        sh.setFormatter(simple_format)

        logger.addHandler(fh)
        logger.addHandler(sh)
    return logger


class GenomeAnnotation:
    def __init__(self, args):
        self.genome = args.genome
        self.protein = args.protein
        self.workdir = args.workdir.rstrip('/')
        self.threads = args.threads
        self.raw_data_dir = args.rawdatadir
        self.species = args.species
        self.TE_anno = args.TE_anno
        self.long = args.long
        self.lib = args.lib
        self.prefix = args.prefix
        self.intermediate_files_set = set()
        self.logger = setup_logger()
        self.hapgene_path = '/home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script'

        try:
            self.out_prefix = re.findall(r'(.*?)\.fa', os.path.basename(self.genome))[0]
        except Exception:
            self.logger.warning('Failed to extract prefix from genome filename. Using default "genome".')
            self.out_prefix = 'genome'

    def build_index(self):
        hisat2_index_dir = os.path.join(self.workdir, "hisat2_index")
        mkdir(hisat2_index_dir)
        prefix = os.path.basename(os.path.splitext(self.genome)[0])
        hisat2_index = os.path.join(hisat2_index_dir, prefix)
        cmd = f"hisat2-build -p 6 {self.genome} {hisat2_index}"
        return cmd, hisat2_index

    def hisat2_bam(self, hisat2_index, read1, read2, bamfile):
        x = os.path.basename(hisat2_index)
        cmd = f"hisat2 -p 7 --dta --summary-file {x}.summary -x {hisat2_index} -1 {read1} -2 {read2} | samtools view -Sb > {bamfile}"
        return cmd

    def run_braker_denovo(self):
        # build index
        index_cmd, hisat2_index = self.build_index()
        cmd_linux(index_cmd)

        bam_dir = os.path.join(self.workdir, "01_rna_bam")
        mkdir(bam_dir)
        raw_data = glob.glob(os.path.join(self.raw_data_dir, "*_1.fastq.gz"))
        if not raw_data:
            self.logger.error("Please rename your fastq data, such as 'A_1.fastq.gz' 'A_2.fastq.gz' !")
            sys.exit(1)

        jobs = []
        for file1 in raw_data:
            basen = os.path.basename(file1)
            prefix = re.findall(r"^(.*?)_1\.fastq\.gz", basen)[0]
            file2 = os.path.join(self.raw_data_dir, prefix + "_2.fastq.gz")
            bamfile = os.path.join(bam_dir, prefix + ".bam")
            hisat2_bam_cmd = self.hisat2_bam(hisat2_index, file1, file2, bamfile)
            jobs.append(hisat2_bam_cmd)

        pool = multiprocessing.Pool(5)
        for cmd in jobs:
            pool.apply_async(cmd_linux, (cmd,))
            time.sleep(0.2)
        pool.close()
        pool.join()

        # filter bam files
        files = glob.glob(os.path.join(bam_dir, "*.bam"))
        jobs = []
        for i in files:
            prefix = os.path.splitext(os.path.basename(i))[0]
            out = os.path.join(bam_dir, prefix + ".q30.bam")
            self.intermediate_files_set.add(i)
            cmd = f'samtools view -@ 15 -O BAM -q 30 -bh -o {out} {i}; rm {i}'
            jobs.append(cmd)
        pool = multiprocessing.Pool(7)
        for cmd in jobs:
            pool.apply_async(cmd_linux, (cmd,))
            time.sleep(0.2)
        pool.close()
        pool.join()

        # run braker.pl
        files = glob.glob(os.path.join(bam_dir, "*.bam"))
        files_str = ",".join(files)
        braker_dir = os.path.join(self.workdir, "02_braker")
        mkdir(braker_dir)
        cmd = f"braker.pl --genome={self.genome} --bam={files_str} --workingdir={braker_dir} --prot_seq={self.protein} --threads {self.threads} --gm_max_intergenic 10000 --skipOptimize --species={self.species} --useexisting"
        cmd_linux(cmd)

        # check and fix annotations
        if not os.path.exists(f'{braker_dir}/braker.gtf'):
            self.logger.error('Not find result of braker! Please check!')
            sys.exit(1)
        cmd_linux(f'mv {braker_dir}/braker.gtf {braker_dir}/braker.tmp.gtf')
        cmd_linux(f'python /home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script/braker_gtf.modify.py {braker_dir}/braker.tmp.gtf {braker_dir}/braker.gtf')
        cmd_linux(f'python /home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script/braker_gtf2gff.py {braker_dir}/braker.gtf {braker_dir}/braker.gff')
        cmd_linux(f'rm {braker_dir}/braker.tmp.gtf')
        os.chdir(braker_dir)
        cmd_linux(f'python /home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script/check_braker_CDS_3.py braker.gtf braker {self.genome}')
        cmd_linux(f'python /home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script/braker_gtf2gff.py braker.rmERROR.gtf braker.rmERROR.gff')

    def run_TE(self):
        TE_dir = os.path.join(self.workdir, "00_TE")
        mkdir(TE_dir)
        os.chdir(TE_dir)
        cmd1 = f'BuildDatabase -name {self.out_prefix}_TEdb {self.genome}'
        cmd_linux(cmd1)
        cmd2 = f'RepeatModeler -database {self.out_prefix}_TEdb -pa 10 -engine ncbi'
        cmd_linux(cmd2)

    def run_miniprot(self):
        protein_dir = os.path.join(self.workdir, "03_miniprot")
        mkdir(protein_dir)
        os.chdir(protein_dir)
        cmd1 = f"miniprot -t 30 --gff --gtf --aln {self.genome} {self.protein} | grep -v '#' | awk '$3 != \"gene\"' > {protein_dir}/miniprot.{self.out_prefix}.gtf"
        cmd_linux(cmd1)
        cmd2 = f'braker_gtf2gff.py {protein_dir}/miniprot.{self.out_prefix}.gtf {protein_dir}/miniprot.{self.out_prefix}.gff'
        cmd_linux(cmd2)

    def run_transdecoder(self):
        transdecoder_dir = os.path.join(self.workdir, "04_transdecoder")
        mkdir(transdecoder_dir)
        os.chdir(transdecoder_dir)
        stringtie_gtf = f'{os.path.join(self.workdir, "02_braker")}/GeneMark-ETP/rnaseq/stringtie/transcripts_merged.gff'
        cmd_linux(f'gtf_genome_to_cdna_fasta.pl {stringtie_gtf} {self.genome} > {transdecoder_dir}/transcripts.fasta')
        cmd_linux(f'gtf_to_alignment_gff3.pl {stringtie_gtf} > {transdecoder_dir}/transcripts.gff3')

        cmd_linux(f'TransDecoder.LongOrfs -t {transdecoder_dir}/transcripts.fasta -m 50')
        cmd_linux(f'diamond makedb --in /home/jhuang/Database/02_UniProt/uniprot_sprot.fasta --db {transdecoder_dir}/uniprot_sprot.fasta')
        cmd_linux(f'diamond blastp -d {transdecoder_dir}/uniprot_sprot.fasta -q {transdecoder_dir}/transcripts.fasta.transdecoder.longest_orfs.pep --evalue 1e-5 --max-target-seqs 1 > {transdecoder_dir}/blastp.outfmt6')
        cmd_linux(f'TransDecoder.Predict -t {transdecoder_dir}/transcripts.fasta --retain_blastp_hits {transdecoder_dir}/blastp.outfmt6')
        cmd_linux(f'cdna_alignment_orf_to_genome_orf.pl {transdecoder_dir}/transcripts.fasta.transdecoder.gff3 {transdecoder_dir}/transcripts.gff3 {transdecoder_dir}/transcripts.fasta > {transdecoder_dir}/transcripts.fasta.transdecoder.genome.gff3')
        cmd_linux(f'python /home/jhuang/script/pycharm/hj/transdecoder_augustus.py {transdecoder_dir}/transcripts.fasta.transdecoder.genome.gff3')
        cmd_linux(f'mv {transdecoder_dir}/longest.cds.gff {transdecoder_dir}/transdecoder.longest.gff')

    def run_EVM(self):
        self.evm_dir = os.path.join(self.workdir, "05_EVM")  
        mkdir(self.evm_dir)
        os.chdir(self.evm_dir)
        weights = f'{self.hapgene_path}/evm.weights.pb.txt' if self.long else f'{self.hapgene_path}/evm.weights.illumina.txt'
        cmd_linux(f'cp {weights} {self.evm_dir}')
        cmd_linux(
            f'gffread {os.path.join(self.workdir, "02_braker")}/braker.rmERROR.gff -T -o {self.evm_dir}/braker.gtf')
        cmd_linux(
            f'perl /home/jhuang/script/pycharm/01_Bio/braker_gtf_to_EVM_gff3.pl {self.evm_dir}/braker.gtf > {self.evm_dir}/braker_evm.gff3')
        cmd_linux(
            f'cat {self.evm_dir}/braker_evm.gff3 {os.path.join(self.workdir, "04_transdecoder")}/transdecoder.longest.gff > {self.evm_dir}/gene_prediction.gff')
        cmd_linux(
            f'bash /home/jhuang/script/pycharm/01_Bio/evm.merge.sh {self.genome} {self.evm_dir}/gene_prediction.gff {os.path.join(self.workdir, "03_miniprot")}/miniprot.{self.out_prefix}.gff {weights}')
        cmd = """
        cat EVM.all.gff3 | awk 'BEGIN {OFS="\t"} { if ($3 == "gene") $3 = "1"; else if ($3 == "mRNA") $3 = "2"; else if ($3 == "exon") $3 = "3"; else if ($3 == "CDS") $3 = "4"; print $0 }' |
        sort -k1,1 -k4,4n -k3,3n |
        awk 'BEGIN {OFS="\t"} { if ($3 == "1") $3 = "gene"; else if ($3 == "2") $3 = "mRNA"; else if ($3 == "3") $3 = "exon"; else if ($3 == "4") $3 = "CDS"; print $0 }' > EVM.all.gff
        """
        cmd_linux(cmd)

    def run_all(self):
        self.run_braker_denovo()
        if self.TE_anno:
            self.run_TE()
        self.run_miniprot()
        self.run_transdecoder()
        self.run_EVM()


def main():
    parser = get_parser()
    args = parser.parse_args()

    # if not args.TE_anno and not args.lib:
    #     parser.error("argument --lib is required when --TE_anno is not specified.")

    base_workdir = args.workdir.rstrip('/')

    for i, (genome_file, prefix_hap) in enumerate(zip(args.genome, args.prefix), 1):
        args_current = vars(args).copy()
        args_current['genome'] = genome_file
        args_current['prefix'] = prefix_hap
        print(genome_file)
        path_prefix = f"0{i}_{os.path.basename(genome_file)}"
        current_workdir = os.path.join(base_workdir, path_prefix)

        args_current['workdir'] = current_workdir
        args_obj = argparse.Namespace(**args_current)

        print(f"Running pipeline for genome: {genome_file} with workdir: {current_workdir}")
        mkdir(current_workdir)

        pipeline = GenomeAnnotation(args_obj)
        pipeline.run_all()

        evm_dir = pipeline.evm_dir
        os.chdir(evm_dir)
        cmd_linux(f'gff_id_rename.py EVM.all.gff hap{i}')
        cmd_linux(f'mv EVM.all.rename.gff {prefix_hap}.gff')


if __name__ == '__main__':
    main()