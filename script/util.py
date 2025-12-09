#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 2024/3/18 10:48
# @Author  : jhuang

import os
import subprocess
import logging


def cmd_linux(cmd):
    print(cmd)
    subprocess.run(cmd, shell=True, close_fds=True)


def mkdir(path):
    """ 创建文件夹 """
    folder = os.path.exists(path)
    if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径
        print("new path: {}".format(path))
    else:
        print("The path is exit !!!")


def log():
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>日志设置>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
    logger = logging.getLogger()  # 创建一个logging对象
    fh = logging.FileHandler("{}.log".format(os.path.join(os.getcwd(), os.path.basename(__file__))),
                             encoding='utf-8')  # 创建一个handler，用于写入日志文件
    sh = logging.StreamHandler()  # 再创建一个handler，用于输出到控制台
    # 配置显示格式
    standard_format = logging.Formatter(
        "[%(asctime)s] [%(filename)s line:%(lineno)d][%(levelname)s]\t%(message)s")
    simple_format = logging.Formatter(
        "[%(levelname)s][%(asctime)s][%(filename)s:%(lineno)d] %(message)s")
    fh.setFormatter(standard_format)
    sh.setFormatter(simple_format)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.setLevel(10)  # 日志水平
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
    return logger

