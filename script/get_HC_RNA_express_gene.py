#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/7/24 21:07
# @Author  : jhuang
# Usage: python remove_Te_related_gene.py $annotation $tmp $output $overlap(%)
# 记得修改无表达基因是否删除！！！(47行)

from util import *
import pandas as pd
import re
import sys


def main():
    annotation = sys.argv[1]
    tmp = sys.argv[2]
    out_put = sys.argv[3]
    overlap = sys.argv[4]
    shred = float(sys.argv[5])

    df = pd.read_csv(tmp, header=0, index_col=0, low_memory=False)
    # 创建一个新的列 'Check'，并使用 apply 方法检查条件
    df['row_sum'] = df.sum(axis=1)
    df['Check'] = df.apply(lambda row: (row > shred).any(), axis=1)  # axis=1表示按行应用
    # 获取 'Check' 列为 True 的行的索引，并存储在一个列表中
    HC_gene = df[df['Check'] == True].index.tolist()
    HC_gene_set = set(HC_gene)

    # 提取 'sum' 列等于 0 的行
    no_expression_set = set(df[df['row_sum'] == 0].index.tolist())

    TErelated_gene = []
    with open("gene.TE.intersect.out", "r") as o1:
        for line in o1:
            line = line.strip()
            if line.startswith("#"):
                continue
            else:
                line_li = line.split("\t")
                if float(line_li[2]) >= float(overlap):  # overlap大于30时，认定与TE有关
                    gene = re.findall(r'(.*?)\.', line_li[1])[0]
                    TErelated_gene.append(gene)
                else:
                    continue

    TErelated_gene_set = set(TErelated_gene)
    remove_gene_set = TErelated_gene_set.difference(HC_gene_set)
    with open("remove.TE.gene.txt", "w") as oo:
        for i in remove_gene_set:
            oo.write(i + "\n")

    length1 = len(remove_gene_set)
    print(f"removed {length1} gene.")
    remove_trans_id_set = set()

    with open(out_put, "w") as o2:
        with open(annotation, "r") as o3:
            for line in o3:
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    line_li = line.split("\t")
                    if line_li[2] == "gene":
                        gene_id = re.findall(r'ID=(.*?);', line_li[-1])[0]
                        if gene_id not in remove_gene_set:
                            o2.write("\t".join(line_li) + "\n")
                    elif line_li[2] == "mRNA":
                        mRNA_id = re.findall(r'ID=(.*?);', line_li[-1])[0]
                        try:
                            parent_gene = re.findall(r'Parent=(.*?);', line_li[-1])[0]
                        except IndexError:
                            parent_gene = re.findall(r'Parent=(.*?)$', line_li[-1])[0]
                        if parent_gene in remove_gene_set:
                            remove_trans_id_set.add(mRNA_id)
                        else:
                            o2.write("\t".join(line_li) + "\n")
                    else:
                        try:
                            parent_mRNA = re.findall(r'Parent=(.*?);', line_li[-1])[0]
                        except IndexError:
                            parent_mRNA = re.findall(r'Parent=(.*?)$', line_li[-1])[0]
                        if parent_mRNA not in remove_trans_id_set:
                            o2.write("\t".join(line_li) + "\n")


if __name__ == '__main__':
    main()
