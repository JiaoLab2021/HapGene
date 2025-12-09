#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/11/8 15:26
# @Author  : jhuang
# 保留注释中最长cds对应转录本
from util import *
import re
import sys
import os
import random


def main():
    if len(sys.argv) != 2:
        print("Usage: python transdecoder_augustus.py <gff_file>")
        exit(1)
    gff_file = sys.argv[1]
    gene_cds_dict = {}
    with open(gff_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            line_li = line.split('\t')
            feature = line_li[2]
            start = int(line_li[3])
            end = int(line_li[4])
            if feature == 'gene':
                gene_id = re.search('ID=(.*?);', line_li[-1]).group(1)
                gene_cds_dict[gene_id] = {}
            if feature == 'mRNA':
                mRNA_id = re.search('ID=(.*?);', line_li[-1]).group(1)
                parent_id = re.search('Parent=(.*?);', line_li[-1]).group(1)
                gene_cds_dict[parent_id][mRNA_id] = 0
            if feature == 'CDS':
                cds_length = end - start + 1
                gene_cds_dict[parent_id][mRNA_id] += cds_length

    result = set()
    for gene_id, mRNA_info in gene_cds_dict.items():
        max_value = max(mRNA_info.values())
        max_keys = [k for k, v in mRNA_info.items() if v == max_value]
        result.add(random.choice(max_keys))

    with open(gff_file, 'r') as f, open('longest.cds.gff', 'w') as o:
        for line in f:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            line_li = line.split('\t')
            feature = line_li[2]
            if feature == 'gene':
                gene_line_li = line_li
                # o.write(line + '\n')
            elif feature == 'mRNA':
                mRNA_id = re.search('ID=(.*?);', line_li[-1]).group(1)
                if mRNA_id in result:
                    gene_line_li[3] = line_li[3]
                    gene_line_li[4] = line_li[4]
                    o.write('\t'.join(gene_line_li) + '\n')
                    o.write(line + '\n')
            else:
                try:
                    parent_id = re.findall('Parent=(.*?);', line_li[-1])[0]
                except IndexError:
                    parent_id = re.findall('Parent=(.*?)$', line_li[-1])[0]
                if parent_id in result:
                    o.write(line + '\n')


if __name__ == '__main__':
    main()
