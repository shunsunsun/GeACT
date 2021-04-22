#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 16:15:09 2019

@author: lixiang
"""

import sys
import pandas as pd
import glob

IN = sys.argv[1]
OUT = sys.argv[2]
ct = sys.argv[3]
cells = sys.argv[4]
t = sys.argv[5]

data = pd.DataFrame()
f = open(cells,'r')
line = f.readline()
while(line):
    l = line.strip('\n').split('\t')
    Dir = '_'.join(l[0].split('_')[:-1])
    if t == 'narrow_split':
        filename = glob.glob('/'.join([IN,'*','02-mapping',Dir,l[0],
                                       '_'.join([l[0],ct,'NarrowPeaks_split.txt'])]))[0]
    elif t == 'narrow':
        filename = glob.glob('/'.join([IN,'*','02-mapping',Dir,l[0],
                                       '_'.join([l[0],ct,'NarrowPeaks.txt'])]))[0]
    elif t == 'TSS':
        filename = glob.glob('/'.join([IN,'*','02-mapping',Dir,l[0],
                                       '_'.join([l[0],ct,'TSS1kb.txt'])]))[0]
    tmp = pd.read_csv(filename,sep='\t',names=[l[1]])
    data = data.append(tmp.T)
    line = f.readline()
f.close()

#data = data.ix[:,(data>0).sum(axis=0) > 0.05*data.shape[0]]

if t == 'narrow_split':
    output_file = '_'.join(['insertcount',ct,'NarrowPeaks_split_filter.tsv.gz'])
elif t == 'narrow':
    output_file = '_'.join(['insertcount',ct,'NarrowPeaks_filter.tsv.gz'])
elif t == 'TSS':
    output_file = '_'.join(['insertcount',ct,'TSS1kb_filter.tsv.gz'])

data.to_csv(OUT+'/'+output_file,sep='\t',compression='gzip')