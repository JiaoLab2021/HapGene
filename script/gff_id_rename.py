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
        gene_n = 1
        feature_counts = defaultdict(int)  # 默认字典用于特征计数

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
                info_gene_last = re.sub(r'ID=[^;]+(?=;|$)', f'ID={prefix}_g{gene_n}', info_gene_last)
                info_gene_last = re.sub(r'Name=[^;]+(?=;|$)', f'Name={prefix}_g{gene_n}', info_gene_last)
                info_gene = '\t'.join(line_li[:-2] + ['.'] + [info_gene_last])
                f2.write(info_gene + '\n')

                gene_n += 1
                trans_n = 1
                feature_counts.clear()  # 清除所有特征的计数

            elif feature_type == 'mRNA':
                info_mRNA_last = line_li[-1]
                info_mRNA_last = re.sub(r'ID=[^;]+(?=;|$)', f'ID={prefix}_g{gene_n-1}.t{trans_n}', info_mRNA_last)
                info_mRNA_last = re.sub(r'Name=[^;]+(?=;|$)', f'Name={prefix}_g{gene_n-1}.t{trans_n}', info_mRNA_last)
                info_mRNA_last = re.sub(r'Parent=[^;]+(?=;|$)', f'Parent={prefix}_g{gene_n-1}', info_mRNA_last)
                info_mRNA = '\t'.join(line_li[:-1] + [info_mRNA_last])
                f2.write(info_mRNA + '\n')

                trans_n += 1
                feature_counts.clear()  # 重置特征计数以适应新的转录本

            else:
                # 处理其他feature
                feature_counts[feature_type] += 1  # 当前特征计数加1
                feature_n = feature_counts[feature_type]
                info_other_last = (f'ID={prefix}_g{gene_n-1}.t{trans_n-1}.{line_li[2]}.{feature_n};'
                                   f'Parent={prefix}_g{gene_n-1}.t{trans_n-1};'
                                   f'Name={prefix}_g{gene_n-1}.t{trans_n-1}.{line_li[2]}.{feature_n}')
                info_other = '\t'.join(line_li[:-1] + [info_other_last])
                f2.write(info_other + '\n')


if __name__ == '__main__':
    main()