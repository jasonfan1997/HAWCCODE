#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 23:01:36 2018

@author: jason
"""

from math import sqrt
import numpy as np
import math
from sklearn.metrics import confusion_matrix
import sys
if len(sys.argv) >= 2:
  seed = int(sys.argv[1])
else:
  seed = 0
np.random.seed(seed)
all_set =np.load("all.npy")



def q_factor(y_true, y_pred, threshold_shift=0):
    
    c = confusion_matrix(y_true,y_pred)
    if c.size!=1:
        hadrone=c[0,1]/float(c[0,0]+c[0,1])
        recall = c[1,1]/float(c[1,0]+c[1,1])
        re=recall/sqrt(hadrone)
        if re !=np.inf and not np.isnan(re):
            return re,recall,hadrone
        else:
            return 0,recall,0
    return 0,0,0

np.random.shuffle(all_set) #shuffle the whole dataset
training_set=all_set[0:int(math.floor(all_set.shape[0]*0.5))] #Split the dataset
prediction_set=all_set[int(math.floor(all_set.shape[0]*0.5)):]
new=prediction_set[:,9]/prediction_set[:,8]
fbin=np.append(np.insert(np.linspace(0.1, 0.9, num=9),0,0.05),1)
ebin=np.linspace(2.5, 5.5, num=13)
cut=np.loadtxt("standardcut.txt",delimiter=' ')
result=[]

for i in range(len(fbin)-1):
    
    
    for j in range(len(ebin)-1):
        temp=None
        true=None
        q=0
        row=None
        true=prediction_set[(new<fbin[i+1])&(new>fbin[i])]
        new2=(prediction_set[:,3])[(new<fbin[i+1])&(new>fbin[i])]

        
        true=true[(new2<ebin[j+1])&(new2>ebin[j])]   
        if i==0:
            result.append([fbin[i],ebin[j],0,0,0,0,0])
        row=cut[cut[:,0]==i]
        row=row[row[:,1]==j]
        if(len(true)!=0 and row.size != 0):

            temp=np.all([true[:,1]<row[0,3],true[:,2]<row[0,2]],axis=0).astype(int)
            t=q_factor(true[:,0],temp)
            #t[np.isnan(t)]=0
            #if t[1]>0.5:
            result.append([fbin[i],ebin[j],t[0],t[1],t[2],row[0,2],row[0,3]])
            #print(fbin[i],ebin[j],t[0],t[1],t[2],t[3])
            #else:
               # result.append([fbin[i],ebin[j],0.5])
            print(fbin[i],ebin[j],t)
        else:
            result.append([fbin[i],ebin[j],0,0,0,0,0])
            
np.savetxt("standardcutq.csv",np.array(result),delimiter=',',fmt='%f')
