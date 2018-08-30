#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 23:00:26 2018

@author: jason
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:20:13 2018

@author: user98
"""

from math import sqrt
import numpy as np
import pickle #pickle
from sklearn.metrics import roc_curve, auc
import sklearn
import pandas as pd
import math
from xgboost import XGBClassifier

seed = 0
np.random.seed(seed)


def qcut(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    q=tpr/np.sqrt(fpr)
    q[q == np.inf] = 0
    q[np.isnan(q)] = 0
    position=np.argmax(q)
    return np.array([q[position],thresholds[position],tpr[position],fpr[position]])

def qcut2(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    q=tpr/np.sqrt(fpr)
    q[q == np.inf] = 0
    q[np.isnan(q)] = 0
    q[tpr<0.2]=0
    #q[fpr>0.2]=0
    position=np.argmax(q)
    return np.array([q[position],thresholds[position],tpr[position],fpr[position]])

def qcut5(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    q=tpr/np.sqrt(fpr)
    q[q == np.inf] = 0
    q[np.isnan(q)] = 0
    q[tpr<0.5]=0
    #q[fpr>0.2]=0
    position=np.argmax(q)
    return np.array([q[position],thresholds[position],tpr[position],fpr[position]])

all_set =np.load("data.npy") #The data array store all the data
#The first column of data.npy will be 0(hadron) and 1(gamma)
#The rest is the input variables
np.random.shuffle(all_set)
training_set=all_set[0:int(math.floor(all_set.shape[0]*0.5))] #Split the dataset
prediction_set=all_set[int(math.floor(all_set.shape[0]*0.5)):]

clf=XGBClassifier(n_estimators=500, silent=False,learning_rate=0.1, max_depth=5, subsample=0.6,gamma=1,n_jobs=-1)
clf.fit(training_set[:,1:],training_set[:,0])
predicted_prob=clf.predict_proba(prediction_set[:,1:])
clf._Booster.dump_model("mcmodel.txt")
new=prediction_set[:,8]
fbin=np.array([0,0.067,0.105,0.162,0.247,0.356,0.485,0.618,0.74,0.84,1.01])
ebin=np.linspace(2.5, 5.5, num=13)
result=[]
for i in range(len(fbin)-1):
    
    
    for j in range(len(ebin)-1):
        temp=None
        true=None

        temp=(predicted_prob[:,1])[(new<fbin[i+1])&(new>fbin[i])]
        true=(prediction_set[:,0])[(new<fbin[i+1])&(new>fbin[i])]
        new2=(prediction_set[:,3])[(new<fbin[i+1])&(new>fbin[i])]
        temp=temp[(new2<ebin[j+1])&(new2>ebin[j])]
        true=true[(new2<ebin[j+1])&(new2>ebin[j])]   
    
        
        if(len(true)!=0):
            t=qcut(true,temp)
            t[np.isnan(t)]=0
            result.append([fbin[i],ebin[j],t[0],t[1],t[2],t[3]])
            print(fbin[i],ebin[j],t[0],t[1],t[2],t[3])
        else:
            result.append([fbin[i],ebin[j],0,0.5,0.5,0.5])

np.savetxt("optcut.csv",np.array(result),delimiter=',',fmt='%f')
result=np.array(result)
ENERGY_WIDTH = 0.25
ENERGY_OFFSET = 2.50
ename="rec.logNNEnergyV2"
#fbin = list(range(10))
ebin = list(range(12))
t=0
fbin=np.array([0,0.067,0.105,0.162,0.247,0.356,0.485,0.618,0.74,0.84,1.01])
for i in range(len(fbin)-1):
    for j in range(len(ebin)):
                if result[t,3]!=0.5 and result[t,3]<=1 and t>=12:
                    other_cut_string = str(t)+" \"{} * rec.nChAvail < rec.nHitSP20 && " \
    						   "rec.nHitSP20 <= {} * rec.nChAvail && "\
    						   "{} < {} && "             \
    						   "{} <= {} && "            \
    						   "rec.angleFitStatus == 0 && "          \
    						   "rec.nChAvail >= 700 && "              \
                            "rec.nChTot >= 800 && "              \
                            "rec.nChAvail>0.9*rec.nChTot && " \
                            "rec.zenithAngle < 0.785 && "              \
    						   "rec.coreFiduScale <= 100.0 &&"\
                            "{}<=rec.proba\"\n".format(fbin[i],fbin[i+1], ENERGY_WIDTH * ebin[j] + ENERGY_OFFSET, ename, ename, ENERGY_WIDTH * (ebin[j] + 1) + ENERGY_OFFSET,result[t,3])#,result[t,3],result[t,5])
                    with open("optcut","a") as myfile:
                        myfile.write(other_cut_string)
                t+=1      
