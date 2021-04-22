#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 16:23:26 2019

@author: lixiang
"""

import os
import sys

import networkx as nx

import pandas as pd

IN = sys.argv[1]
OUT = sys.argv[2]
sp = sys.argv[3]

Dirs = []
line = sys.stdin.readline()
while(line):
    Dirs.append(line.strip('\n'))
    line = sys.stdin.readline()

ids = []
for i in range(1,13):
    for j in range(1,9):
        ids.append('%d-%d' %(i,j))

frag_cell = dict()
cell_id = []
for Dir in Dirs:
    for i in ids:
        filename = '/'.join([IN,Dir,Dir+'_'+i,'_'.join([Dir,i,sp,'frag.bed'])])
        if os.path.exists(filename):
            cell_id.append(Dir+'_'+i)
            f = open(filename,'r')
            line = f.readline()
            while(line):
                l = line.strip('\n').split('\t')
                if '\t'.join(l[0:3]) not in frag_cell:
                    frag_cell['\t'.join(l[0:3])] = {Dir+'_'+i:(int(l[4]),set(l[3].split(',')))}
                else:
                    frag_cell['\t'.join(l[0:3])][Dir+'_'+i] = (int(l[4]),set(l[3].split(',')))
                line = f.readline()
            f.close()

print(sys.getsizeof(frag_cell))

cell_graph = nx.Graph()
cell_graph.add_nodes_from(cell_id)
for i,c1 in enumerate(cell_id[:-1]):
    l1 = c1.split('_')
    for c2 in cell_id[i+1:]:
        l2 = c2.split('_')
        if l1[-1] == l2[-1] or (l1[-2] == l2[-2] and \
             (l1[-1].split('-')[0]==l2[-1].split('-')[0] or \
              l1[-1].split('-')[1]==l2[-1].split('-')[1])):
            cell_graph.add_edge(c1,c2)

def argmax_dict(D,keys=None):
    if not keys:
        keys=D.keys()
    m = 0
    arg_m = []
    for i in keys:
        if D[i][0]>m:
            m = D[i][0]
            arg_m = [i]
        elif D[i][0]==m:
            arg_m.append(i)
    return(arg_m)

con_freq = pd.DataFrame(0,index=cell_id,columns=cell_id)

decon_file = dict()
con_file = dict()
for i in cell_id:
    Dir = '_'.join(i.split('_')[:-1])
    decon_file[i] = open('/'.join([IN,Dir,i,'_'.join([i,sp,'frag_decon.bed'])]),'w')
    con_file[i] = open('/'.join([IN,Dir,i,'_'.join([i,sp,'frag_con.bed'])]),'w')
        
for frag in frag_cell:
    if len(frag_cell[frag])==1:
        for i in frag_cell[frag]:
            decon_file[i].write('\t'.join([frag,','.join(sorted(list(frag_cell[frag][i][1]))),
                            str(frag_cell[frag][i][0])])+'\n')
            con_freq.loc[i,i] += 1
    else:
        H = cell_graph.subgraph(frag_cell[frag].keys()).copy()
        remove_edges = []
        for i,j in H.edges():
            if not (frag_cell[frag][i][1].issubset(frag_cell[frag][j][1]) or
                    frag_cell[frag][j][1].issubset(frag_cell[frag][i][1])):
                remove_edges.append((i,j))
        H.remove_edges_from(remove_edges)
        for cc in nx.connected_components(H):
            if len(cc) == 1:
                for i in cc:
                    decon_file[i].write('\t'.join([frag,','.join(sorted(list(frag_cell[frag][i][1]))),
                            str(frag_cell[frag][i][0])])+'\n')
                    con_freq.loc[i,i] += 1
            else:
                m = argmax_dict(frag_cell[frag],keys=cc)
                if len(m)==1:
                    for i in m:
                        decon_file[i].write('\t'.join([frag,','.join(sorted(list(frag_cell[frag][i][1]))),
                            str(frag_cell[frag][i][0])])+'\n')
                        con_freq.loc[i,i] += 1
                else:
                    for i in m:
                        con_file[i].write('\t'.join([frag,','.join(sorted(list(frag_cell[frag][i][1]))),
                            str(frag_cell[frag][i][0])])+'\n')
                for i in m:
                    for j in cc:
                        if j != i:
                            con_freq.loc[j,i] += 1

con_freq.to_csv(OUT+'/'+'con_freq_'+sp+'.csv.gz',compression='gzip')
for i in cell_id:
    decon_file[i].close()
    con_file[i].close()