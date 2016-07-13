# -*- coding: utf-8 -*-
"""
Created on Thu May 26 09:43:53 2016

@author: Zi
"""

import numpy as np
from collections import Counter
common=np.load('common.npz')
gm_only=np.load('gm-only.npz')
k_only=np.load('k-only.npz')
chrs=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
Sigs=['Cmyc','Cux1','Hcfc1','Max','Mxi1','Sp1','Spl1','Srf','Tbp','Yy1','Znf143','Znf384']
listc=[[],[],[],[],[],[],[],[],[],[],[],[]]
listg=[[],[],[],[],[],[],[],[],[],[],[],[]]
listk=[[],[],[],[],[],[],[],[],[],[],[],[]]
for c in chrs:
    for i in range(len(common[c])):
        for j in range(len(Sigs)):
            listc[j].append(tuple(common[c][i][j+1]))
    for i in range(len(k_only[c])):
        for j in range(len(Sigs)):
            listg[j].append(tuple(k_only[c][i][j+1]))
    for i in range(len(gm_only[c])):
        for j in range(len(Sigs)):
            listk[j].append(tuple(gm_only[c][i][j+1]))
for i in range(len(listc)):
    print Sigs[i],'\t',Counter(listc[i])[(0,0)],'\t',Counter(listc[i])[(0,1)],'\t',Counter(listc[i])[(1,0)],'\t',Counter(listc[i])[(1,1)],'\t',len(listc[1])
for i in range(len(listg)):
    print Sigs[i],'\t',Counter(listg[i])[(0,0)],'\t',Counter(listg[i])[(0,1)],'\t',Counter(listg[i])[(1,0)],'\t',Counter(listg[i])[(1,1)],'\t',len(listg[1])
for i in range(len(listk)):
    print Sigs[i],'\t',Counter(listk[i])[(0,0)],'\t',Counter(listk[i])[(0,1)],'\t',Counter(listk[i])[(1,0)],'\t',Counter(listk[i])[(1,1)],'\t',len(listk[1])
        
        