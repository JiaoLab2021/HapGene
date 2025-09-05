#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/8/23 10:36
# @Author  : jhuang

import re
from util import *
import sys


def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_gene_seq_header_modify.py <gff_anno> <prefix> <genome>")
        exit(1)
    gff_anno = sys.argv[1]
    prefix = sys.argv[2]
    genome = sys.argv[3]
    cmd_linux(f'gfftransformer -gff {gff_anno} gln -p {prefix}')
    cmd_linux(f'seqkit subseq --bed {prefix}.bed -o {prefix}.fa {genome}')

    input_file_fasta = f'{prefix}.fa'
    with open(input_file_fasta, 'r') as f, open('rename00.fa', 'w') as o:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header_li = re.split(r'\s+', line)
                # print(header_li)
                header = header_li[1]
                o.write(f'>{header}\n')
            else:
                o.write(line + '\n')

    cmd_linux(f'rm {input_file_fasta}; mv rename00.fa {input_file_fasta}')


if __name__ == '__main__':
    main()
