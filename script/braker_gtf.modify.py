#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/3/29 10:47
# @Author  : jhuang
# 修改braker输出的gtf文件（格式不标准）
# Usage： python braker_gtf.modify.py $braker.gtf $modified.name
import re
import sys
from util import *


def main():
    before_gtf = sys.argv[1]
    modified_gtf = sys.argv[2]
    with open(modified_gtf, "w") as f:
        with open(before_gtf, "r") as f1:
            for line in f1:
                line = line.strip()
                line_li = line.split("\t")
                if line_li[2] == "gene":
                    continue
                if line_li[2] == "mRNA":
                    continue
                if line_li[2] == "transcript":
                    trans_id = line_li[-1]
                    gene_id = re.findall(r"^(.*?)\.", trans_id)[0]
                    info = f"transcript_id \"{trans_id}\"; gene_id \"{gene_id}\";"
                    f.write("\t".join(line_li[0:8]) + "\t" + info + "\n")
                else:
                    f.write(line + "\n")


if __name__ == '__main__':
    main()
