# -*- coding: utf-8 -*-
"""
Created on Fri May 13 09:21:41 2016

@author: Zi
"""
import os
import numpy as np
datafolder1 = 'C:\\Users\\Zi\\Desktop\\motif_alone'
datafolder2 = 'C:\\Users\\Zi\\Desktop'
D1={}
D2={}
chrs=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']
for chorm in chrs:
   file1=open(os.path.join(datafolder1, 'Ctcf-only.txt'))
   file2=open(os.path.join(datafolder2, 'co_Ctcf_narrowpeak.txt'))
   arr1 = []
   arr2 = []   
   for line in file1:
         line=line.strip('\n')
         list1=line.split('\t')
         if list1[1]==chorm :
             arr1.append(list1[2:4])
   nparr1=np.array(arr1,dtype=int)
   D1[chorm]=arr1
   
   for line in file2:
         list2=line.split('\t')
         if list2[0]==chorm :
             arr2.append(list2[1:3])
   nparr2=np.array(arr2,dtype=int)
   D2[chorm]=arr2
   file1.close()
   file2.close()

bind={}
unbind={}
for i in chrs:
    motif=D1[i]
    peak=D2[i]
    b=[]
    ub=[]
    ubset=[]
    usets=[]
    for x in peak:
        for y in motif:
            if int(y[0])>=int(x[0]) and int(y[1])<=int(x[1]):
                mid=(int(y[0])+int(y[1]))/2
                b.append(mid)
                
    for u in motif:
        mid=(int(u[0])+int(u[1]))/2
        ubset.append(mid)
    ub=list(set(ubset)-set(b))
    b=list(set(b))
    b.sort()
    ub.sort()
    for n in range(len(ub)/3):
        usets.append(ub[n*3])
    bind[i]=b
    unbind[i]=usets
    
np.savez('ctcf-Ctcf_bind.npz',**bind)           
#np.savez('ctcf-unbind.npz',**unbind) 


