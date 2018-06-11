# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:20:13 2018

@author: user98
"""

from math import sqrt
import numpy as np
import pickle #pickle
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import sklearn
import sys
if len(sys.argv) >= 2:
  seed = int(sys.argv[1])
else:
  seed = 0
np.random.seed(seed)
from sklearn import preprocessing

import pandas as pd
import math
from xgboost import XGBClassifier
from sklearn import preprocessing

GAMMADIR = "./run.csv"
COLUMNS = ["rec.nHit/U/1","rec.CxPE40/F/0.01","rec.PINC/F/0.01","rec.logNNEnergy/F/0.01","rec.disMax/F/0.01","rec.LDFAge/F/0.01"
,"rec.LDFAmp/F/0.01","rec.nChAvail/U/1","rec.nHitSP20/U/1"]
allset = pd.read_csv(GAMMADIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
allset=allset[allset[:,1]!=0,:]
allset[:,0]= np.log10(np.divide(allset[:,1], allset[:,0],out=np.zeros_like(allset[:,0]), where=allset[:,0]!=0))
allset=np.delete(allset,1,1)


def q_factor(y_true, y_pred):
    
    precision = sklearn.metrics.precision_score(y_true,y_pred)
    nrecall = 1-sklearn.metrics.recall_score(y_true,y_pred)
    #return recall
    #beta_squared = beta ** 2
    return precision/sqrt(nrecall)

def qcut(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    q=tpr/np.sqrt(fpr)
    q[q == inf] = 0
    return np.max(q),thresholds[np.argmax(q)]
    


with open('clf.pickle', 'rb') as f:
    clf = pickle.load(f)
    
cut=pd.read_csv("cut.csv",skiprows=0).as_matrix()
predicted_prob=clf.predict_proba(allset)

'''
new=prediction_set[:,8]/prediction_set[:,7]
fbin=np.append(np.insert(np.linspace(0.1, 0.9, num=9),0,0.05),1)
ebin=np.linspace(2.5, 5.5, num=13)
result=[]
for i in range(len(fbin)):
	for j in range(len(ebin)):
        temp=None
        true=None
        if i == 0:
            temp=(predicted_prob[:,1])[np.where(new<fbin[i])]
            true=(prediction_set[:,0])[np.where(new<fbin[i])]
            new2=(prediction_set[:,3])[np.where(new<fbin[i])]
        else:
            #condlist = [predicted_prob<fbin[i], predicted_prob>fbin[i-1]]
            temp=(predicted_prob[:,1])[(new<fbin[i])&(new>fbin[i-1])]
            true=(prediction_set[:,0])[(new<fbin[i])&(new>fbin[i-1])]
            new2=(prediction_set[:,3])[(new<fbin[i])&(new>fbin[i-1])]

        if j == 0:
            temp=temp[np.where(new2<ebin[j])]
            true=true[np.where(new2<ebin[j])]
        else:
            temp=temp[(new2<ebin[j])&(new2>ebin[j-1])]
            true=true[(new2<ebin[j])&(new2>ebin[j-1])]   
            

        if(len(true)!=0):
            t=qcut(true,temp)
            result.append([fbin[i],ebin[j],t[0],t[1]])
            print(fbin[i],ebin[j],t[0],t[1])

'''