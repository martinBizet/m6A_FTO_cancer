# ------------------------------------------------------------
#
# Matthieu Defrance ULB 2016
# lib/base
#
# ------------------------------------------------------------

import os
import sys
import random
import time
from config import *

# ------------------------------------------------------------
# select
# ------------------------------------------------------------
def select(cond, ifTrue, ifFalse):
    if cond:
        return ifTrue
    else:
        return ifFalse

# ------------------------------------------------------------
# log / run
# ------------------------------------------------------------
def log(msg):
    sys.stderr.write('[%s] %s\n' % (time.strftime('%Y/%m/%d %H:%M'), msg))
    sys.stderr.flush()

def run(*args):
    cmd = ' '.join(['%s' % str(v) for v in args])
    log('%s' % cmd)
    if not DEBUG:
        os.system(cmd)

# ------------------------------------------------------------
# csv read / write
# ------------------------------------------------------------
def read_csv(filename, sep = ',', header = True):
    f = open(filename)
    lines = f.readlines()
    lines = lines[1:]
    return [line.strip().split(sep) for line in lines]

# ------------------------------------------------------------
# tmp
# ------------------------------------------------------------
def tmp_dir(default_dir = TMP_DIR):
    tmp = '%s/%s%d/' % (default_dir, time.strftime('%y%m%d%H%M'), random.randint(1, 1000))
    run('mkdir -p %s' % tmp)
    return tmp

def tmp_file(default_dir = TMP_DIR):
    tmp = '%s/%s%d.tmp' % (default_dir, time.strftime('%y%m%d%H%M'), random.randint(1, 1000))
    return tmp

def clean_tmp(tmp):
    run('rm -rf %s' % tmp)
