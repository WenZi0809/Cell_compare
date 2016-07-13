# -*- coding: utf-8 -*-
"""
Created on Wed May 25 18:59:01 2016

@author: Zi
"""

import numpy as np
import os
gmdata=np.load('/public/home/3dgenome/wenzi/co-binding/gm12878/ctcf-bind.npz')
kdata=np.load('/public/home/3dgenome/wenzi/co-binding/k562/ctcf-bind.npz')
gmfold='/public/home/3dgenome/wenzi/co-binding/gm12878/'
kfold='/public/home/3dgenome/wenzi/co-binding/k562'
Sigs=['Cmyc','Cux1','Hcfc1','Max','Mxi1','Sp1','Spl1','Srf','Tbp','Yy1','Znf143','Znf384']
chrs=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
wind=750
gk11={}
gk01={}
gk10={}
for c in chrs:   
    gk11[c]=list(set(kdata[c]) & set(gmdata[c]))
    gk10[c]=list(set(gmdata[c]) - set(kdata[c]))
    gk01[c]=list(set(kdata[c]) - set(gmdata[c]))

tf11={}
tf10={}
tf01={}    
for chrom in chrs:
    common=gk11[chrom]
    gm=gk10[chrom]
    k=gk01[chrom]
   
    commondatabase=np.zeros([len(common),len(Sigs)+1,2],int) ##common(1,1)
    gmdatabase=np.zeros([len(gm),len(Sigs)+1,2],int)##(1,0)
    kdatabase=np.zeros([len(k),len(Sigs)+1,2],int)##(0,1)
    
    for x in range(len(common)): 
        commondatabase[x][0][0]=1
        commondatabase[x][0][1]=1 
        
    for x in range(len(gm)): 
        gmdatabase[x][0][0]=1
        gmdatabase[x][0][1]=0 
        
    for x in range(len(k)): 
        kdatabase[x][0][0]=0
        kdatabase[x][0][1]=1 
        
    for n in range(len(Sigs)):
        sigchr=[] ##gm12878
        desFil = 'Gm12878'+Sigs[n]+'.narrowPeak'   
        filein = open(os.path.join(gmfold, desFil))
        for line in filein:
            lists=line.split('\t')
            if lists[0]==chrom:
                sigchr.append(lists[1:3])
        for i in sigchr:
            for j in range(len(common)):
                if common[j]>=((int(i[0])+int(i[1]))/2-wind) and common[j]<=((int(i[0])+int(i[1]))/2+wind):
                     commondatabase[j][n+1][0]=1
            for p in range(len(gm)):
                if gm[p]>=((int(i[0])+int(i[1]))/2-wind) and gm[p]<=((int(i[0])+int(i[1]))/2+wind):
                     gmdatabase[p][n+1][0]=1
            for l in range(len(k)):
                if k[l]>=((int(i[0])+int(i[1]))/2-wind) and k[l]<=((int(i[0])+int(i[1]))/2+wind):
                     kdatabase[l][n+1][0]=1

        ksigchr=[] ##k562
        kdesFil = 'K562'+Sigs[n]+'.narrowPeak'   
        kfilein = open(os.path.join(kfold, kdesFil))
        for line in kfilein:
            lists=line.split('\t')
            if lists[0]==chrom:
                ksigchr.append(lists[1:3])
        for i in ksigchr:
            for j in range(len(common)):
                if common[j]>=((int(i[0])+int(i[1]))/2-wind) and common[j]<=((int(i[0])+int(i[1]))/2+wind):
                     commondatabase[j][n+1][1]=1
            for p in range(len(gm)):
                if gm[p]>=((int(i[0])+int(i[1]))/2-wind) and gm[p]<=((int(i[0])+int(i[1]))/2+wind):
                     gmdatabase[p][n+1][1]=1
            for l in range(len(k)):
                if k[l]>=((int(i[0])+int(i[1]))/2-wind) and k[l]<=((int(i[0])+int(i[1]))/2+wind):
                     kdatabase[l][n+1][1]=1
    tf11[chrom]=commondatabase
    tf10[chrom]=gmdatabase
    tf01[chrom]=kdatabase
np.savez('common',**tf11)
np.savez('gm-only',**tf10)
np.savez('k-only',**tf01)
    
    
    
    
    
    
    
    
    