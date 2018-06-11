#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:16:47 2018

@author: jason
"""

from math import sqrt
import numpy as np
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn import linear_model
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
import matplotlib.pyplot as plt
from sklearn.ensemble import BaggingClassifier
from sklearn import preprocessing
from sklearn.svm import SVC
GAMMADIR = "./gamma.csv"
HADRONDIR="./hadron.csv"
COLUMNS = ["rec.nHit/U/1","rec.CxPE40XnCh/U/1","rec.PINC/F/0.01","rec.logNNEnergy","rec.disMax/F/0.01","rec.LDFAge/F/0.01"
,"rec.LDFAmp/F/0.01"]
allset = pd.read_csv(GAMMADIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS).as_matrix()
hadron = pd.read_csv(HADRONDIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS).as_matrix()
allset=np.hstack((np.ones((len(allset),1)),allset))
hadron = np.hstack((np.zeros((len(hadron),1)),hadron))
allset=np.vstack((allset,hadron))
hadron=None
np.random.shuffle(allset)
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
    return precision/K.sqrt(recall)


model=Sequential() #Define the model

adam=keras.optimizers.Adam(lr=0.002)
model.add(Dense(7,input_shape=(allset.shape[1]-1,),activation='sigmoid')) #Add kernel_regularizer 
model.add(Dense(5,activation='sigmoid'))
model.add(Dense(1, activation='sigmoid'))
class_weight = {0 : 1.,1: len(allset)/np.count_nonzero(allset[:,0]==1),}  #calculate the class_weight
model.compile(loss='binary_crossentropy', optimizer=adam, metrics=[fbeta,qfactor,'accuracy'])
model.fit(allset[:,1:],allset[:,0],epochs=50,batch_size=1000,class_weight=class_weight,verbose=1,shuffle=False,validation_split=0.1) #turn off the shuffle for reproducible result

