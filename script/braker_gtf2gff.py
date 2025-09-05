#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/5/19 16:53
# @Author  : jhuang
# Transform the GTF annotation file from Braker into GFF format
# Usage: pytoonbraker_gtf2gff.py $gtf_file $output

import re
import sys


def main():
    if len(sys.argv) != 3:
        print("Usage: python braker_gtf2gff.py <gtf_file> <output_gff>")
        exit(1)
    gtf_file = sys.argv[1]
    gff_file = sys.argv[2]
    gene_num = set()

    with open(gtf_file, "r") as f1:
        with open(gff_file, "w") as f2:
            for line in f1:
                line = line.strip()
                if line.startswith("#"):
                    continue
                else:
                    line_list = line.split("\t")
                    symbol_type = line_list[2]
                    try:
                        transcript_id = re.findall(r'transcript_id "(.*?)";', line_list[-1])[0]
                    except IndexError:
                        transcript_id = re.findall(r'transcript_id "(.*?)"$', line_list[-1])[0]
                    try:
                        gene_id = re.findall(r'gene_id "(.*?)";', line_list[-1])[0]

                    except IndexError:
                        gene_id = re.findall(r'gene_id "(.*?)"$', line_list[-1])[0]
                    gene_information = line_list.copy()
                    gene_information[2] = "gene"
                    gene_information[-1] = f"ID={gene_id};Name={gene_id}"

                    if gene_id not in gene_num:
                        gene_num.add(gene_id)
                        f2.write("\t".join(gene_information) + "\n")
                    if symbol_type == "transcript":
                        other_symbol = {}
                        transcript_information = line_list.copy()
                        transcript_information[2] = "mRNA"
                        transcript_information[-1] = f"ID={transcript_id};Parent={gene_id};Name={transcript_id}"
                        f2.write("\t".join(transcript_information) + "\n")
                    else:
                        if symbol_type not in other_symbol:
                            other_symbol[symbol_type] = 1
                            symbol_name = f"{transcript_id}.{symbol_type}.{other_symbol[symbol_type]}"
                            f2.write("\t".join(line_list[:-1]) + "\t" + f"ID={symbol_name};Parent={transcript_id};Name={symbol_name}" + "\n")
                        else:
                            other_symbol[symbol_type] += 1
                            symbol_name = f"{transcript_id}.{symbol_type}.{other_symbol[symbol_type]}"
                            f2.write("\t".join(line_list[:-1]) + "\t" + f"ID={symbol_name};Parent={transcript_id};Name={symbol_name}" + "\n")


if __name__ == '__main__':
    main()
