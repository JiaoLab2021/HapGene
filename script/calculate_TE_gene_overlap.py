#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/7/4 20:41
# @Author  : jhuang
# EDTA overlap
# gff文件exon最后1列必须包含ID=*;Parent=*

import os
import re
import sys


def main():
    TEbed = sys.argv[1]
    annotation = sys.argv[2]
    gene_bed = "gene.bed"
    overlap_file = "gene.TE.intersect.out"

    # ======================= gene.bed ===================
    len_dict = {}
    with open(annotation, 'r') as f1, open(gene_bed, 'w') as o1:
        for line in f1:
            line = line.strip()
            if line.startswith('#'):
                continue
            line_li = line.split("\t")
            chr = line_li[0]
            type_g = line_li[2]
            start = int(line_li[3])
            end = int(line_li[4])
            if type_g == "exon":
                info = line_li[-1].split(";")
                exon_id = re.findall(r"ID=(.*?)$", info[0])[0]
                gene_id = re.findall(r"Parent=(.*?)$", info[1])[0]
                if gene_id not in len_dict.keys():
                    len_dict.update({gene_id: end - start + 1})
                else:
                    len_dict[gene_id] += end - start + 1
                o1.write(chr + "\t" + str(start) + "\t" + str(end) + "\t" +
                         info[0] + "\t" + info[1] + "\n")
    # ===============================================================

    # intersectBed
    os.system(f"intersectBed -a {gene_bed} -b {TEbed} -wao | sort -k1,1 -k5,5 -k2,2n -k7,7n > intersectBed.txt")

    # overlap ca
    overlap = 0
    overlap_dict = {}
    g_id, scaf, g_start, g_end, inf, ex_id, ex_id2 = None, None, None, None, None, None, None
    te_start, te_end = None, None

    with open("intersectBed.txt", 'r') as f2, open(overlap_file, 'w') as o2:
        for line in f2:
            line = line.strip()
            line_li = line.split("\t")
            g_id2 = re.findall(r"Parent=(.*?)$", line_li[4])[0]
            ex_id2 = re.findall(r"ID=(.*?)$", line_li[3])[0]
            if not g_id:  # 初始行操作
                g_id = g_id2
                ex_id = ex_id2
                if int(line_li[9]) > 0:
                    overlap = int(line_li[9]) + 1
                scaf, g_start, g_end, inf = line_li[0], int(line_li[1]), int(line_li[2]), line_li[3]
                te_start, te_end = int(line_li[6]), int(line_li[7])
            elif g_id != g_id2:
                overlap = f"{100 * overlap / len_dict[g_id]:.1f}"
                # print(g_id + ":" + overlap + "\n")
                o2.write(scaf + "\t" + g_id + "\t" + str(overlap) + "\n")
                overlap_dict.update({g_id: overlap})
                # if g_id not in overlap_dict.keys():
                #     overlap_dict.update({g_id: overlap})
                # else:
                overlap = 0

                g_id = g_id2
                ex_id = ex_id2
                if int(line_li[9]) > 0:
                    overlap = int(line_li[9]) + 1
                scaf, g_start, g_end, inf = line_li[0], int(line_li[1]), int(line_li[2]), line_li[3]
                te_start, te_end = int(line_li[6]), int(line_li[7])
            elif ex_id != ex_id2:
                if int(line_li[9]) > 0:
                    overlap += int(line_li[9]) + 1
                ex_id = ex_id2
                te_start, te_end = int(line_li[6]), int(line_li[7])
            else:
                if int(line_li[6]) >= te_start and int(line_li[7]) <= te_end:
                    continue
                elif int(line_li[6]) >= te_start and int(line_li[6]) <= te_end:
                    if te_end+1 > int(line_li[2]):
                        continue
                    else:
                        if int(line_li[2]) <= int(line_li[7]):
                            overlap += int(line_li[2]) - (te_end+1) + 1
                        else:
                            overlap += int(line_li[7]) - (te_end+1) + 1
                    te_end = int(line_li[7])
                else:
                    overlap += int(line_li[9]) + 1
                    te_start, te_end = int(line_li[6]), int(line_li[7])

        overlap = f"{100 * overlap / len_dict[g_id]:.1f}"
        overlap_dict[g_id] = overlap
        o2.write(scaf + "\t" + g_id + "\t" + str(overlap) + "\n")


if __name__ == '__main__':
    main()
