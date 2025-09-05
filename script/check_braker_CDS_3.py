#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/6/26 20:22
# @Author  : jhuang
# Identify genes with abnormal CDS lengths and transcripts containing incorrect start codons in the annotation
# modified Time: 2024/9/23 09:14
# Identify genes or transcripts with premature stop codons in the annotation results
# modified Time: 2024/12/19 15:11
# Remove genes or transcripts lacking stop codons
import re
from util import *
import sys


def main():
    if len(sys.argv) != 4:
        print("Usage: python check_braker_CDS_3.py <anno_gtf> <prefix> <genome.fa>")
        exit(1)
    file = sys.argv[1]
    prefix = sys.argv[2]
    output_removed_gtf = f'{prefix}.rmERROR.gtf'
    output_error = f'{prefix}_error.gtf'
    genome = sys.argv[3]
    dict = {}
    with open(file, 'r') as f1:
        for line in f1:
            line = line.strip()
            line_li = line.split('\t')
            chr = line_li[0]
            type = line_li[2]
            start = int(line_li[3])
            end = int(line_li[4])
            trans_id = re.findall(r'^transcript_id \"(.*?)\";', line_li[-1])[0]
            if type == "transcript":
                dict.update({trans_id: []})
            if type == "CDS":
                cds_len = end - start + 1
                dict[trans_id].append(cds_len)

    cds_error_trans = set()
    st_coden_error_trans = set()
    advance_stop_coden_error_trans = set()
    all_error_trans = set()
    for trans_id, len_li in dict.items():
        size = sum(len_li)
        if size % 3 != 0:
            cds_error_trans.add(trans_id)
            all_error_trans.add(trans_id)

    # Translate coding sequences to protein sequences
    cmd_linux(f'braker_gtf2gff.py {file} {prefix}.gff')
    cmd_linux(f'python /home/jhuang/script/pycharm/hj/gfftrans.py -gff {prefix}.gff cds -g {genome} -cdsp {prefix} -prop {prefix}')

    # Retrieve transcript IDs with annotation errors
    # cmd_linux(f'grep -B1 -v -e "^M" -e "^>" gffread.{prefix}.protein.fa | grep ">" | sed "s/>//g" > start_stop_coden_error.txt')
    with open(f'{prefix}.protein.fasta', 'r') as f2, open('advance_stop_coden_error.txt', 'w') as o2, open('start_stop_coden_error.txt', 'w') as o3:
        for line in f2:
            line = line.strip()
            if line.startswith('>'):
                current_id = re.findall(r'>(.*?)$', line)[0]
                # print(current_id)
                continue
            if not line.startswith('M') or not line.endswith('*'):
                o3.write(current_id + '\n')
            if '*' in line[:-1]:
                o2.write(current_id + '\n')

    with open("start_stop_coden_error.txt", 'r') as s1:
        for line in s1:
            line = line.strip()
            st_coden_error_trans.add(line)
            all_error_trans.add(line)
    with open("advance_stop_coden_error.txt", 'r') as s2:
        for line in s2:
            line = line.strip()
            advance_stop_coden_error_trans.add(line)
            all_error_trans.add(line)
    with open("cds_error_trans.txt", 'w') as q1, open("all_error_trans.txt", 'w') as q2:
        for i in cds_error_trans:
            q1.write(i + "\n")
        for j in all_error_trans:
            q2.write(j + "\n")

    with open(output_error, 'w') as f2, open(output_removed_gtf, 'w') as f3, open(file, "r") as f4:
        for line in f4:
            line = line.strip()
            if line.startswith('#'):
                continue
            line_li = line.split('\t')
            transcript_id = re.findall(r'transcript_id "(.*?)";', line_li[-1])[0]
            if transcript_id in all_error_trans:
                f2.write("\t".join(line_li) + "\n")
            else:
                f3.write("\t".join(line_li) + "\n")


if __name__ == '__main__':
    main()
