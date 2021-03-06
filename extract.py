#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:16:47 2018

@author: jason
false_positive_rate, true_positive_rate, thresholds = roc_curve(training_set[:,0], result[:,2])
roc_auc = auc(false_positive_rate, true_positive_rate)
best=thresholds[np.argmin(false_positive_rate-true_positive_rate)]
"""

import numpy as np
if len(sys.argv) >= 2:
  seed = int(sys.argv[1])
else:
  seed = 0
np.random.seed(seed)
import pandas as pd

GAMMADIR = "./gamma.csv"
HADRONDIR="./hadron.csv"
COLUMNS = ["rec.nHit/U/1","rec.CxPE40/F/0.01","rec.PINC/F/0.01","rec.logNNEnergy/F/0.01","rec.disMax/F/0.01","rec.LDFAge/F/0.01"
,"rec.LDFAmp/F/0.01","rec.nChAvail/U/1","rec.nHitSP20/U/1"]
allset = pd.read_csv(GAMMADIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
hadron = pd.read_csv(HADRONDIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
#allset[:,0]=allset[:,0]/allset[:,1]
allset=allset[allset[:,1]!=0,:]
hadron=hadron[hadron[:,1]!=0,:]
allset[:,0]= np.log10(np.divide(allset[:,1], allset[:,0],out=np.zeros_like(allset[:,0]), where=allset[:,0]!=0))
allset=np.delete(allset,1,1)
allset=np.hstack((np.ones((len(allset),1)),allset))
hadron[:,0]= np.log10(np.divide(hadron[:,1],hadron[:,0], out=np.zeros_like(hadron[:,0]), where=hadron[:,0]!=0))
hadron=np.delete(hadron,1,1)
hadron = np.hstack((np.zeros((len(hadron),1)),hadron))
allset=np.vstack((allset,hadron))
hadron=None
np.random.shuffle(allset)

np.save("all.npy",allset)
