#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/11/11 21:22
# @Author  : jhuang


from util import *
import re
import sys
import os


def main():
    if len(sys.argv) != 3:
        print("Usage: python choose_gff_by_gene_id <gene_id_file> <gff_file>")
        exit(1)
    gene_id_file = sys.argv[1]
    gff_file = sys.argv[2]
    gene_set = set()
    with open(gene_id_file, 'r') as f:
        for line in f:
            line = line.strip()
            try:
                gene_id = re.findall(r'(.*?)\.', line)[0]
            except IndexError:
                gene_id = line
            gene_set.add(gene_id)

    trans_set = set()
    with open(gff_file, 'r') as f2:
        for line in f2:
            line = line.strip()
            if line.startswith('#') or not line:
                continue
            line_li = line.split('\t')
            feature = line_li[2]
            if feature == "gene":
                gene_id = re.findall('ID=(.*?);', line_li[-1])[0]
                if gene_id in gene_set:
                    print(line)
            elif feature == "mRNA":
                try:
                    parent_id = re.findall('Parent=(.*?);', line_li[-1])[0]
                except IndexError:
                    parent_id = re.findall('Parent=(.*?)$', line_li[-1])[0]
                trans_id = re.findall('ID=(.*?);', line_li[-1])[0]
                if parent_id in gene_set:
                    print(line)
                    trans_set.add(trans_id)
            else:
                try:
                    parent_id = re.findall('Parent=(.*?);', line_li[-1])[0]
                except IndexError:
                    parent_id = re.findall('Parent=(.*?)$', line_li[-1])[0]
                if parent_id in trans_set:
                    print(line)


if __name__ == '__main__':
    main()
