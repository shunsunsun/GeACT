#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 23:46:06 2019

@author: lixiang
"""

import sys

l = []
line = sys.stdin.readline()
l.append(int(line.strip('\n').strip(' ').split(' ')[0]))
line = sys.stdin.readline()
l.append(int(line.strip('\n').strip(' ').split(' ')[0]))
line = sys.stdin.readline()
line = sys.stdin.readline()
l.append(int(line.strip('\n').strip(' ').split(' ')[0]))
line = sys.stdin.readline()
l.append(int(line.strip('\n').strip(' ').split(' ')[0])+l[2])

l.append(round(l[2]/l[0],4))
l.append(round(l[3]/l[0],4))

print('\t'.join([str(k) for k in l]))