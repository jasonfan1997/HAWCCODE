# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:20:13 2018

@author: user98
"""

from math import sqrt
import numpy as np
import pickle #pickle模块
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
from keras.models import Model, Sequential
from keras.layers import Dense, GlobalAveragePooling2D,Dropout, Flatten,Input
from keras import backend as K
from keras.metrics import top_k_categorical_accuracy, sparse_top_k_categorical_accuracy
from keras import applications
from keras import regularizers
from sklearn import preprocessing
import keras
from keras.objectives import categorical_crossentropy
from sklearn.model_selection import cross_val_score
import pandas as pd
import math
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from sklearn.ensemble import BaggingClassifier
from sklearn import preprocessing
from sklearn.svm import SVC
import time
import datetime

starttime = datetime.datetime.now()
#GAMMADIR = "./gamma.csv"
#HADRONDIR="./hadron.csv"
#COLUMNS = ["rec.nHit/U/1","rec.CxPE40XnCh/U/1","rec.PINC/F/0.01","rec.disMax/F/0.01","rec.LDFAge/F/0.01"
#,"rec.LDFAmp/F/0.01"]

all_set =np.load("all.npy")
def fbeta(y_true, y_pred, threshold_shift=0):
    beta = 1

    # just in case of hipster activation at the final layer
    y_pred = K.clip(y_pred, 0, 1)

    # shifting the prediction threshold from .5 if needed
    y_pred_bin = K.round(y_pred + threshold_shift)

    tp = K.sum(K.round(y_true * y_pred_bin)) + K.epsilon()
    fp = K.sum(K.round(K.clip(y_pred_bin - y_true, 0, 1)))
    fn = K.sum(K.round(K.clip(y_true - y_pred, 0, 1)))

    precision = tp / (tp + fp)
    recall = tp / (tp + fn)

    beta_squared = beta ** 2
    return (beta_squared + 1) * (precision * recall) / (beta_squared * precision + recall + K.epsilon())


def qfactor(y_true, y_pred, threshold_shift=0):
    #beta = 1

    # just in case of hipster activation at the final layer
    #y_pred = K.clip(y_pred, 0, 1)
    gammatrue=K.sum(y_true)
    hadrontrue=K.sum( K.cast(K.equal(y_true,0), 'float32'))
    #gammapred=K.sum(K.cast(K.equal(K.cast(K.equal(K.round(y_pred),y_true),'float32'),y_true),'float32'))
    #hadronfalse=K.sum(K.cast(K.equal(K.cast(K.equal(y_true,0),'float32'),K.round(y_pred)),'float32'))
    # shifting the prediction threshold from .5 if needed
    gammapred=K.sum(K.cast(K.equal(K.round(y_pred),y_true), 'float32')*y_true)
    hadronfalse=K.sum(K.cast(K.equal(K.round(y_pred),1), 'float32')*K.cast(K.equal(y_true,0), 'float32'))
    precision = gammapred / (gammatrue+K.epsilon())
    recall = hadronfalse / (hadrontrue+K.epsilon())
    #return recall
    #beta_squared = beta ** 2
    return precision/K.sqrt(recall+K.epsilon())

def q_factor(y_true, y_pred, threshold_shift=0):
    
    precision = sklearn.metrics.precision_score(y_true,y_pred)
    nrecall = 1-sklearn.metrics.recall_score(y_true,y_pred)
    #return recall
    #beta_squared = beta ** 2
    return precision/sqrt(nrecall)

def qcut(y_true, y_pred):
    fpr, tpr, thresholds = roc_curve(y_true, y_pred)
    q=tpr/np.sqrt(fpr)
    q[q == inf] = 0
    position=np.argmax(q)
    return q[position],thresholds[position],tpr[position],fpr[position]
    

np.random.shuffle(all_set) #shuffle the whole dataset
training_set=all_set[0:math.floor(all_set.shape[0]*0.9)] #Split the dataset
prediction_set=all_set[math.floor(all_set.shape[0]*0.9):]

clf=XGBClassifier(n_estimators=4200, silent=False,learning_rate=0.0019, max_depth=5, subsample=0.6,n_jobs=-1)
clf.fit(training_set[:,1:],training_set[:,0])
predicted_prob=clf.predict_proba(prediction_set[:,1:])
predicted_label=(predicted_prob[:,1]>=0.5).astype(int)
print(q_factor(prediction_set[:,0],predicted_label))
print(sklearn.metrics.f1_score(prediction_set[:,0],predicted_label))
with open('clf.pickle', 'wb') as f:
    pickle.dump(clf, f)

endtime = datetime.datetime.now()
print ((endtime - starttime).seconds)
'''
with open('clf.pickle', 'rb') as f:
    clf = pickle.load(f)
    
predicted_prob=clf.predict_proba(prediction_set[:,1:])
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
            result.append([fbin[i],ebin[j],t[0],t[1],t[2],t[3]])
            print(fbin[i],ebin[j],t[0],t[1],t[2],t[3])

'''