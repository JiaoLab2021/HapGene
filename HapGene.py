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
import pandas as pd


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
        # print(genome_file)
        path_prefix = f"0{i}_{prefix_hap}"
        current_workdir = os.path.join(base_workdir, path_prefix)

        args_current['workdir'] = current_workdir
        args_obj = argparse.Namespace(**args_current)

        evm_file = os.path.join(current_workdir, "05_EVM",
                                f"{prefix_hap}.EVM.all.gff")
        print(evm_file)
        # exit(1)
        if os.path.exists(evm_file):
            print(f"Found existing EVM result for {genome_file}, skipping pipeline.")
            continue

        print(f"Running pipeline for genome: {genome_file} with workdir: {current_workdir}")
        mkdir(current_workdir)

        pipeline = HapGenomeAnnotation(args_obj)
        pipeline.run_all()

        evm_dir = pipeline.evm_dir
        os.chdir(evm_dir)
        # cmd_linux(f'gff_id_rename.py {prefix}.EVM.all.gff hap{i}')
        # cmd_linux(f'mv {prefix}.EVM.all.rename.gff {prefix_hap}.gff')

    polish = AnnotationPolish(args)
    polish.run_all()

    filter = AnnotationFilter(args)
    filter.run_all()


# 先激活环境braker3
def get_parser():
    parser = argparse.ArgumentParser(description="This script parallelizes the commands of the input file.")

    required = parser.add_argument_group('required arguments')
    required.add_argument("-g", "--genome", metavar="genome1 genome2", help="Input your genome", required=True, nargs='+')
    required.add_argument("-w", "--workdir", metavar="workdir", help="Input your workdir.", required=True)
    required.add_argument("--protein", metavar="protein", help="Input your protein sequence.", required=True)
    required.add_argument("-t", "--threads", type=int, metavar="threads", help="Input your threads.", required=True)
    required.add_argument("-r", "--rawdatadir", metavar="rawdatadir", help="Input rawdata dir.", required=True)
    required.add_argument("-s", "--species", metavar="species", help="Input species name.", required=True)
    required.add_argument("-p", "--prefix", metavar="prefix1 prefix2", help="Output prefix.", required=True, nargs='+')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--TE_anno', action='store_true', help='Run TE annotation step', required=False)
    optional.add_argument('--long', action='store_true', help='Pacbio long-reads mode', required=False)
    optional.add_argument('--lib', metavar="TE_library", help='Library file for RepeatMasker', required=False)
    optional.add_argument('--threshold', type=float, default=0.9, help='Threshold for polish filtering alignments (default: 0.9)')
    optional.add_argument("--lencf", metavar="length_cutoff", type=float,  help="Mono-exonic genes length cutoff [default=1000]", required=False, default=1000)
    optional.add_argument("--tpmcf", metavar="TPM_Value_cutoff", type=float, help="Mono-exonic genes TPM value cutoff [default=1]", required=False, default=1)
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


class HapGenomeAnnotation:
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
        cmd = f"/home/jhuang/biosoft/braker_v3.0.8/scripts/braker.pl --genome={self.genome} --bam={files_str} --workingdir={braker_dir} --prot_seq={self.protein} --threads {self.threads} --gm_max_intergenic 10000 --skipOptimize --species={self.species} --useexisting"
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
        cmd1 = f'BuildDatabase -name {self.prefix}_TEdb {self.genome}'
        cmd_linux(cmd1)
        cmd2 = f'RepeatModeler -database {self.prefix}_TEdb -pa 10 -engine ncbi'
        cmd_linux(cmd2)

    def run_miniprot(self):
        protein_dir = os.path.join(self.workdir, "03_miniprot")
        mkdir(protein_dir)
        os.chdir(protein_dir)
        cmd1 = f"miniprot -t 30 --gff --gtf --aln {self.genome} {self.protein} | grep -v '#' | awk '$3 != \"gene\"' > {protein_dir}/miniprot.{self.prefix}.gtf"
        cmd_linux(cmd1)
        cmd2 = f'braker_gtf2gff.py {protein_dir}/miniprot.{self.prefix}.gtf {protein_dir}/miniprot.{self.prefix}.gff'
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
            f'bash /home/jhuang/script/pycharm/01_Bio/evm.merge.sh {self.genome} {self.evm_dir}/gene_prediction.gff {os.path.join(self.workdir, "03_miniprot")}/miniprot.{self.prefix}.gff {weights}')
        cmd = f"""
        cat EVM.all.gff3 | \
        awk 'BEGIN {{OFS="\\t"}} {{ if ($3 == "gene") $3 = "1"; else if ($3 == "mRNA") $3 = "2"; else if ($3 == "exon") $3 = "3"; else if ($3 == "CDS") $3 = "4"; print $0 }}' | \
        sort -k1,1 -k4,4n -k3,3n | \
        awk 'BEGIN {{OFS="\\t"}} {{ if ($3 == "1") $3 = "gene"; else if ($3 == "2") $3 = "mRNA"; else if ($3 == "3") $3 = "exon"; else if ($3 == "4") $3 = "CDS"; print $0 }}' > {self.prefix}.EVM.all.gff
        """
        cmd_linux(cmd)
        cmd = f'gff_id_rename.py {self.prefix}.EVM.all.gff {self.prefix}_evm'
        cmd_linux(cmd)

    def run_all(self):
        self.run_braker_denovo()
        if self.TE_anno:
            self.run_TE()
        self.run_miniprot()
        self.run_transdecoder()
        self.run_EVM()


class AnnotationPolish:
    def __init__(self, args):
        self.genome = args.genome
        self.workdir = args.workdir.rstrip('/')
        self.TE_anno = args.TE_anno
        self.long = args.long
        self.lib = args.lib
        self.prefix = args.prefix
        self.intermediate_files_set = set()
        self.logger = setup_logger()
        self.hapgene_path = '/home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script'
        self.threshold = args.threshold

    def genome_gff_modify(self):
        # ====================== Step1: 所有基因组和注释修改与整合 ======================
        modified_dir = f'{self.workdir}/03_polish/00_merge'
        polish_workdir = f'{self.workdir}/03_polish'
        mkdir(modified_dir)
        # 创建每个单倍型的polish目录
        polish_dirs = []
        gff_files = []
        genome_files = []
        os.chdir(modified_dir)

        for i, (genome_file, prefix_hap) in enumerate(zip(self.genome, self.prefix), 1):
            evm_file = os.path.join(self.workdir, f"0{i}_{prefix_hap}/05_EVM", f"{prefix_hap}.EVM.all.rename.gff")

            # 处理每个单倍型的GFF和基因组文件
            self.gff_chr_modi(genome_file, evm_file, prefix_hap, modified_dir)

            polish_dir = f'{polish_workdir}/0{i}_{prefix_hap}'
            mkdir(polish_dir)
            polish_dirs.append(polish_dir)
            gff_files.append(f'{modified_dir}/{prefix_hap}.gff')
            genome_files.append(f'{modified_dir}/{prefix_hap}.genome.fa')

        # 合并所有GFF和基因组文件
        cmd_linux(f'cat {" ".join(gff_files)} > {modified_dir}/all.gff')
        cmd_linux(f'cat {" ".join(genome_files)} > {modified_dir}/all.fa')

        cmd_linux(f'python {self.hapgene_path}/extract_gene_seq_header_modify.py all.gff all.gene all.fa')

        if not os.path.exists(f'{modified_dir}/mmseq_rep_seq.fasta'):
            cmd_linux(f'mmseqs easy-cluster all.gene.fa mmseq tmp_dir --min-seq-id 0.9 -c 0.9 --threads 36')

    def process_gene_data(self, gff, prefix, genome):
        cmd_linux(
            f'python {self.hapgene_path}/extract_gene_seq_header_modify.py {gff} {prefix}.gene {genome}')
        cmd_linux(
            f'cut -f 1,2,3,4 {prefix}.gene.bed > {prefix}.gene.bed1; mv {prefix}.gene.bed1 {prefix}.gene.bed')

    def mmseq_align(self):
        modified_dir = f'{self.workdir}/03_polish/00_merge'
        polish_workdir = f'{self.workdir}/03_polish'
        for i, (genome_file, prefix_hap) in enumerate(zip(self.genome, self.prefix), 1):
            # ====================== Step1: 对现有gene序列聚类，blastn比对到基因组 ======================
            os.chdir(modified_dir)
            # 处理每个单倍型的基因数据
            self.process_gene_data(f'{modified_dir}/{prefix_hap}.gff', prefix_hap, f'{modified_dir}/{prefix_hap}.genome.fa')

            os.chdir(f'{polish_workdir}/0{i}_{prefix_hap}')
            # 对每个单倍型进行mmseqs比对
            db_path = f'{prefix_hap}_db'
            aln_res = f'{prefix_hap}.alnRes.m8'
            if not os.path.exists(aln_res):
                cmd_linux(f'mmseqs createdb {modified_dir}/{prefix_hap}.genome.fa {db_path}')
                cmd_linux(f'mmseqs createindex {db_path} tmp --search-type 3')
                cmd_linux(f'mmseqs easy-search {modified_dir}/mmseq_rep_seq.fasta {db_path} {aln_res} tmp '
                          f'--cov-mode 2 --threads 36 -c 0.8 --search-type 3 '
                          f'--format-output "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qcov"')

            # 对每个单倍型进行polish和最终处理
            self.polish(prefix_hap, self.threshold, f'{modified_dir}/all.gff', f'{modified_dir}/{prefix_hap}.genome.fa',
                   f'{modified_dir}/all.fa')

            if os.path.exists('maybe_miss_gene.txt'):
                if os.path.getsize('maybe_miss_gene.txt') == 0:
                    print("文件为空")
                else:
                    self.add_final_gff(prefix_hap, f'{modified_dir}/{prefix_hap}.genome.fa', f'{modified_dir}/{prefix_hap}.gff')
            else:
                print("文件不存在")

            # 删除中间文件
            cmd_linux(f'rm *.fai *.mmi')

    def gff_chr_modi(self, genome, gff, prefix, modified_dir):
        gff_df = pd.read_csv(gff, sep='\t', comment='#', header=None,
                             names=['chr', 'source', 'feature', 'start', 'end', 'score',
                                    'strand', 'phase', 'attributes'])
        gff_df['chr'] = prefix + "_" + gff_df['chr'].astype(str)
        gff_df.to_csv(f'{modified_dir}/{prefix}.gff', sep='\t', index=False, header=False)

        # 处理 genome.fa
        output_lines = []
        with open(genome, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    name_id = re.findall(r'>(.*?)$', line)[0]
                    output_lines.append(f'>{prefix}_{name_id}\n')
                else:
                    output_lines.append(line + '\n')

        # 一次性写入，提高效率
        with open(f'{modified_dir}/{prefix}.genome.fa', 'w') as o:
            o.writelines(output_lines)

    def polish(self, target_prefix, threshold, query_gff, target_genome, query_genome):
        # ====================== Step2: 根据identity、coverage阈值,筛选blastn结果 ========================
        blastn_df = pd.read_csv(f'{target_prefix}.alnRes.m8', sep='\t', header=None)
        blastn_df.columns = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend',
                             'sstart',
                             'send', 'evalue', 'bitscore', 'qcovs']
        # print(blastn_df)
        gene_dict = {}
        with open(f'{self.workdir}/03_polish/00_merge/mmseq_rep_seq.fasta') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    gene_id = re.findall(r'>(.*?)$', line)[0]
                else:
                    length = len(line)
                    gene_dict[gene_id] = length
        # 使用map()方法来将基因ID映射到gene_dict中的基因长度
        blastn_df['gene_length'] = blastn_df['qseqid'].map(gene_dict)
        blastn_df['length_percent'] = (blastn_df['qend'] - blastn_df['qstart'] + 1) / blastn_df['gene_length'] * 100
        filtered_df = blastn_df[(blastn_df['pident'] > threshold) &
                                (blastn_df['qcovs'] > threshold) &
                                (blastn_df['length_percent'] > threshold)]

        def swap_columns(row):
            if row['sstart'] > row['send']:
                row['sstart'], row['send'] = row['send'], row['sstart']
            return row

        # 应用该函数到整个 DataFrame
        filtered_df = filtered_df.apply(swap_columns, axis=1)
        filtered_df = filtered_df[['sseqid', 'sstart', 'send', 'qseqid']]
        filtered_df.to_csv('blastn.filtered.bed', sep='\t', index=False, header=False)
        # print(filtered_df)
        # ===========================================================================================

        # ====================== Step2: bedtools查看overlap ==========================================
        cmd = f'bedtools intersect -a blastn.filtered.bed -b {self.workdir}/03_polish/00_merge/{target_prefix}.gene.bed -wao > bedtools.out'
        cmd_linux(cmd)
        # ===========================================================================================

        # ====================== Step3: 根据bedtools结果筛选miss gene =================================
        dup_set = set()
        single_set_tmp = set()  # 唯一没有overlap gene
        with open('bedtools.out', 'r') as f:
            for line in f:
                line = line.strip()
                line_li = line.split('\t')
                # print(line_li[3])
                try:
                    query = re.findall(r'(.*?)\.', line_li[3])[0]
                except IndexError:
                    query = line_li[3]
                target = line_li[7]

                if query in dup_set:
                    continue
                else:
                    if target == '.':
                        single_set_tmp.add(query)
                    else:
                        dup_set.add(query)
                        if query in single_set_tmp:
                            single_set_tmp.remove(query)
        # print(single_set_tmp)
        # print(dup_set)
        #
        with open('maybe_miss_gene.txt', 'w') as f3:
            for i in single_set_tmp:
                f3.write(i + '\n')

        if os.path.exists('maybe_miss_gene.txt'):
            if os.path.getsize('maybe_miss_gene.txt') == 0:
                print("文件为空")
            else:
                cmd_linux(
                    f'choose_gff_by_gene_id.py maybe_miss_gene.txt {query_gff} > ref_maybe_miss_gene.gff')  # ===========================================================================================

                # ====================== Step3: 根据bedtools结果筛选miss gene,进行LiftOn转移注释 =================================
                # LiftOn
                if os.path.getsize('ref_maybe_miss_gene.gff') != 0:
                    cmd = f'lifton -g ref_maybe_miss_gene.gff -o lifton_out.gff -copies {target_genome} {query_genome}'
                    cmd_linux(cmd)

                # 移除错误注释
                cmd = f'gffread lifton_out.gff -T -o lifton_out.gtf;python /home/jhuang/script/pycharm/01_Bio/check_braker_CDS_3.py lifton_out.gtf lifton_out_miss {target_genome}'
                cmd_linux(cmd)

                cmd = f'braker_gtf2gff.py lifton_out_miss.rmERROR.gtf lifton_out_miss.rmERROR.gff'
                cmd_linux(cmd)
                # ==================================================================================================
        else:
            print("文件不存在")

    def add_final_gff(self, target_prefix, target_genome, target_gff):
        # 先改名
        cmd_linux(
            f'gff_id_rename.py lifton_out_miss.rmERROR.gff lifton_{target_prefix}')
        # lifton_prefix = re.findall(r'(.*?)\.gff', lifton_rmERROR_gff)[0]
        cmd_linux(
            f'python /home/jhuang/script/pycharm/06_potato/extract_gene_seq_header_modify.py lifton_out_miss.rmERROR.rename.gff lifton_out_miss.gene {target_genome}')
        cmd_linux(
            f'cat lifton_out_miss.gene.bed | cut -f 1,2,3,4 > lifton_out_miss.gene.bed1; mv lifton_out_miss.gene.bed1 lifton_out_miss.gene.bed')
        cmd_linux(
            f'extract_gene_seq_header_modify.py {target_gff} {target_prefix}.gene {target_genome};cat {target_prefix}.gene.bed | cut -f 1,2,3,4 > {target_prefix}.gene.bed1; mv {target_prefix}.gene.bed1 {target_prefix}.gene.bed')
        cmd_linux(
            f'bedtools intersect -a lifton_out_miss.gene.bed -b {target_prefix}.gene.bed -wao > lifton_original.bedtools')

        bedtools_out_lifton_original_df = pd.read_csv('lifton_original.bedtools', sep='\t',
                                                      header=None)
        bedtools_out_lifton_original_df.columns = ['lchrom', 'lstart', 'lend', 'lid', 'tchrom', 'tstart', 'tend',
                                                   'tid',
                                                   'score']

        # 与原始注释无交叉
        df_no = bedtools_out_lifton_original_df[bedtools_out_lifton_original_df['tid'] == '.']
        # 与原始注释有交叉
        df_yes = bedtools_out_lifton_original_df[bedtools_out_lifton_original_df['tid'] != '.']
        # print(bedtools_out_lifton_transdecoder_df)
        # print(df_no)
        # print(df_yes)
        df_no_bed = df_no[['lchrom', 'lstart', 'lend', 'lid']]
        df_yes_bed = df_yes[['tchrom', 'tstart', 'tend', 'tid']]
        df_no_bed.to_csv('not_have_intersect.bed', sep='\t', header=False, index=False)
        df_yes_bed.to_csv('have_intersect.bed', sep='\t', header=False, index=False)

        # 与原注释交叉情况
        gene_dict_exon = {}
        with open('lifton_out_miss.rmERROR.rename.gff') as f:
            for line in f:
                line = line.strip()
                line_li = line.split('\t')
                if line_li[2] == 'gene':
                    gene_id = re.findall(r'ID=(.*?);', line_li[-1])[0]
                if line_li[2] == 'mRNA':
                    gene_dict_exon[gene_id] = 0
                if line_li[2] == 'exon':
                    gene_dict_exon[gene_id] += 1
        # for gene_id, num in gene_dict_exon.items():
        #     print(gene_id + '\t' + str(num))

        dfff = pd.read_csv('lifton_original.bedtools', sep='\t', header=None)
        dfff.columns = ['lchrom', 'lstart', 'lend', 'lid', 'tchrom', 'tstart', 'tend', 'tid', 'score']
        dfff['exon_num'] = dfff['lid'].map(gene_dict_exon)
        # dfff_no = dfff[(dfff['tid'] == '.') & (dfff['exon_num'] != 1)]  # 过滤单外显子基因
        dfff_no = dfff[(dfff['tid'] == '.')]
        dfff_no_id = dfff_no[['lid']]
        dfff_no_id.to_csv('lifton.id.txt', sep='\t', index=False, header=False)
        # print(dfff_no)
        unique_count1 = dfff_no_id['lid'].nunique()
        print('Final gene number in LiftOn: ' + str(unique_count1))

        cmd_linux(
            f'choose_gff_by_gene_id.py lifton.id.txt lifton_out_miss.rmERROR.rename.gff > lifton.id.gff')
        cmd_linux(
            f'cat lifton.id.gff {target_gff} > {target_prefix}_final_polish_find.gff;'
            f'sort_gff.py {target_prefix}_final_polish_find.gff;mv sorted.gff {target_prefix}_final_polish_find.gff3;rm {target_prefix}_final_polish_find.gff')

        # # 评估
        # cmd_linux(f'gffread {target_prefix}/{target_prefix}_final_polish_find.gff3 -g {target_genome} -y {target_prefix}/{target_prefix}.pep.fa')
        # cmd_linux(f'conda run -n braker3 busco -c 50 -m prot -l /home/jhuang/Database/busco/01_embryophyta_odb10 --offline '
        #           f'-i {target_prefix}/{target_prefix}.pep.fa -o {target_prefix}/{target_prefix}.pep.fa_busco_1614')
        # cmd_linux(f'conda run -n compleasm compleasm protein -p {target_prefix}/{target_prefix}.pep.fa '
        #           f'-l /home/jhuang/Database/busco/01_embryophyta_odb10 -t 50 -o {target_prefix}/{target_prefix}.pep.fa_compleasm_1614')

    def run_all(self):
        self.genome_gff_modify()
        self.mmseq_align()


class AnnotationFilter:
    def __init__(self, args):
        self.genome = args.genome
        self.workdir = args.workdir.rstrip('/')
        self.TE_anno = args.TE_anno
        self.long = args.long
        self.lib = args.lib
        self.prefix = args.prefix
        self.intermediate_files_set = set()
        self.logger = setup_logger()
        self.hapgene_path = '/home/jhuang/research/01_Citrus_sinensis/06_T2T_val/16_HapGene/script'
        self.threshold = args.threshold
        self.lencf = args.lencf
        self.tpmcf = args.tpmcf
        # required.add_argument("--gff", metavar="gff", help="Gff annotation", required=True)
        # required.add_argument("--go", metavar="GO", help="Interproscan result GO", required=True)
        # required.add_argument("--tpm", metavar="TPM_Value", help="Expression TPM", required=True)

    def gene_length_exon_num(self, input_gff, prefix):
        exon_num_dict = {}
        gene_length_dict = {}
        length_cut_set = set()
        gene_all_set = set()
        with open(input_gff, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    continue
                line_li = line.split('\t')
                feature = line_li[2]
                length1 = int(line_li[4]) - int(line_li[3]) + 1
                if feature == 'gene':
                    gene_id1 = re.findall(r'ID=(.*?);', line_li[-1])[0]
                    exon_num_dict[gene_id1] = {}
                    gene_length_dict[gene_id1] = length1
                if feature == 'mRNA':
                    trans_id = re.findall(r'ID=(.*?);', line_li[-1])[0]
                    gene_id = re.findall(r'Parent=(.*?);', line_li[-1])[0]
                    exon_num_dict[gene_id][trans_id] = 0
                if feature == 'exon':
                    exon_num_dict[gene_id][trans_id] += 1

        with open(f'{prefix}.gene.length.txt', 'w') as o:
            o.write('# gene_id\texon_num\tgene_length\ttrans_id\n')
            for gene_id, info in exon_num_dict.items():
                for trans_id, exon_num in info.items():
                    if exon_num == 1 and gene_length_dict[gene_id] <= self.lencf:
                        length_cut_set.add(gene_id)
                    o.write(
                        gene_id + '\t' + str(exon_num) + '\t' + str(gene_length_dict[gene_id]) + '\t' + trans_id + '\n')

        for gene_id in gene_length_dict.keys():
            gene_all_set.add(gene_id)

        return gene_all_set, length_cut_set

    def domain_go(self, go, gene_all_set):
        go_with_anno_set = set()
        with open(go, 'r') as go_file:
            for line in go_file:
                if line.startswith('#'):
                    continue
                line_li = line.strip().split('\t')
                go_trans_id = line_li[0]
                go_gene_id = re.findall(r'^(.*?)\.', go_trans_id)[0]
                go_with_anno_set.add(go_gene_id)

        go_without_anno_set = gene_all_set - go_with_anno_set
        return go_without_anno_set

    def expression_tpm(self, tpm):
        tpm_df = pd.read_csv(tpm, sep=',', header=0, index_col=0)
        # tpm_df.columns.values[0] = 'gene_id'
        df_filtered = tpm_df[(tpm_df <= self.tpmcf).all(axis=1)]
        expression_tpm_set = set(df_filtered.index)
        return expression_tpm_set

    def run_all(self):
        polish_workdir = f'{self.workdir}/03_polish'
        for i, (genome_file, prefix_hap) in enumerate(zip(self.genome, self.prefix), 1):
            os.chdir(f'{polish_workdir}/0{i}_{prefix_hap}')
            gff = f'{polish_workdir}/0{i}_{prefix_hap}/{prefix_hap}_final_polish_find.gff3'
            gff_df = pd.read_csv(gff, sep='\t', comment='#', header=None,
                                 names=['chr', 'source', 'feature', 'start', 'end', 'score',
                                        'strand', 'phase', 'attributes'])
            gff_df['chr'] = gff_df['chr'].str.replace(f'^{prefix_hap}_', '', regex=True)
            polish_gff = f'{polish_workdir}/0{i}_{prefix_hap}/{prefix_hap}.mo.polish.gff'
            gff_df.to_csv(polish_gff, sep='\t', index=False, header=False)
            cmd = f'gffread {polish_gff} -T -o {prefix_hap}.mo.polish.gtf'
            cmd = f'gffread {polish_gff} -g {genome_file} -y {prefix_hap}.mo.polish.pep.fa'
            cmd_linux(cmd)

            # run interproscan
            cmd1 = f'interproscan.sh -cpu 30 -i {prefix_hap}.mo.polish.pep.fa -b {prefix_hap}.GO -iprlookup -pa --goterms -f TSV --goterms'
            cmd_linux(cmd1)

            # expression
            cmd2 = f'featureCounts -T 50 -t exon -p -g gene_id -a {prefix_hap}.mo.polish.gtf -o {prefix_hap}.featureCounts.gene.txt {self.workdir}/0{i}_{prefix_hap}/01_rna_bam/*bam'
            cmd_linux(cmd2)
            cmd3 = f'Rscript {self.hapgene_path}/cal_tpm.R {prefix_hap}.featureCounts.gene.txt {prefix_hap}'
            cmd_linux(cmd3)

            gene_all_set, length_cut_set = self.gene_length_exon_num(polish_gff, prefix_hap)

            go = ''
            go_without_anno_set = self.domain_go(go, gene_all_set)

            tpm = ''
            expression_tpm_set = self.expression_tpm(tpm)

            false_positive_gene_set = length_cut_set & go_without_anno_set & expression_tpm_set
            false_positive_gene_list = sorted(list(false_positive_gene_set))

            filtered_gene_list = sorted(list(gene_all_set - false_positive_gene_set))

            with open(f'{prefix_hap}_false_positive_genes.txt', 'w') as o1:
                for i in false_positive_gene_list:
                    o1.write(i + '\n')

            with open(f'{prefix_hap}_filtered.txt', 'w') as o2:
                for j in filtered_gene_list:
                    o2.write(j + '\n')

            cmd = f'{self.hapgene_path}/choose_gff_by_gene_id.py {prefix_hap}_filtered.txt {polish_gff} > {prefix_hap}.filterd.gff'
            cmd_linux(cmd)


if __name__ == '__main__':
    main()