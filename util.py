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
    """ Create a directory """
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        print("new path: {}".format(path))
    else:
        print("The path is exit !!!")


def log():
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> log setting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> #
    logger = logging.getLogger()
    fh = logging.FileHandler("{}.log".format(os.path.join(os.getcwd(), os.path.basename(__file__))),
                             encoding='utf-8')
    sh = logging.StreamHandler()

    standard_format = logging.Formatter(
        "[%(asctime)s] [%(filename)s line:%(lineno)d][%(levelname)s]\t%(message)s")
    simple_format = logging.Formatter(
        "[%(levelname)s][%(asctime)s][%(filename)s:%(lineno)d] %(message)s")
    fh.setFormatter(standard_format)
    sh.setFormatter(simple_format)
    logger.addHandler(fh)
    logger.addHandler(sh)
    logger.setLevel(10)
    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< #
    return logger

