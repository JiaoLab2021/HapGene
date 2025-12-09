#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/8/23 10:00
# @Author  : jhuang
# 修改gff注释文件id
# 格式为prefix_C*00010...

import sys
import re
import os
from collections import defaultdict


def main():
    if len(sys.argv) != 3:
        print("Usage: python gff_id_rename.py <input_gff> <prefix>")
        exit(1)
    input_gff = sys.argv[1]
    prefix = sys.argv[2]
    pre = re.findall(r'(.*?)\.gff.*$', os.path.basename(input_gff))[0]
    output_gff = pre + '.rename.gff'
    with open(input_gff) as f1, open(output_gff, 'w') as f2:
        gene_n = 10
        feature_counts = defaultdict(int)

        for line in f1:
            line = line.strip()
            if line.startswith('#'):
                f2.write(line + '\n')
                continue

            line_li = line.split('\t')
            feature_type = line_li[2]
            chr = line_li[0]

            if feature_type == 'gene':
                info_gene_last = line_li[-1]
                info_gene_last = re.sub(r'ID=[^;]+(?=;|$)', f'ID={prefix}{chr}_g{str(gene_n).zfill(6)}', info_gene_last)
                info_gene_last = re.sub(r'Name=[^;]+(?=;|$)', f'Name={prefix}{chr}_g{str(gene_n).zfill(6)}', info_gene_last)
                info_gene = '\t'.join(line_li[:8] + [info_gene_last])
                f2.write(info_gene + '\n')

                gene_n += 10
                trans_n = 1
                feature_counts.clear()

            elif feature_type == 'mRNA':
                info_mRNA_last = line_li[-1]
                info_mRNA_last = re.sub(r'ID=[^;]+(?=;|$)', f'ID={prefix}{chr}_g{str(gene_n - 10).zfill(6)}.{trans_n}', info_mRNA_last)
                info_mRNA_last = re.sub(r'Name=[^;]+(?=;|$)', f'Name={prefix}{chr}_g{str(gene_n - 10).zfill(6)}.{trans_n}', info_mRNA_last)
                info_mRNA_last = re.sub(r'Parent=[^;]+(?=;|$)', f'Parent={prefix}{chr}_g{str(gene_n - 10).zfill(6)}', info_mRNA_last)
                info_mRNA = '\t'.join(line_li[:-1] + [info_mRNA_last])
                f2.write(info_mRNA + '\n')

                trans_n += 1
                feature_counts.clear()

            else:
                # 处理其他feature
                feature_counts[feature_type] += 1
                feature_n = feature_counts[feature_type]
                info_other_last = (f'ID={prefix}{chr}_g{str(gene_n - 10).zfill(6)}.{trans_n - 1}.{line_li[2]}.{feature_n};'
                                   f'Parent={prefix}{chr}_g{str(gene_n - 10).zfill(6)}.{trans_n - 1};'
                                   f'Name={prefix}{chr}_g{str(gene_n - 10).zfill(6)}.{trans_n - 1}.{line_li[2]}.{feature_n}')
                info_other = '\t'.join(line_li[:-1] + [info_other_last])
                f2.write(info_other + '\n')


if __name__ == '__main__':
    main()