#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  9 15:16:47 2018

@author: jason
"""
import sys
import numpy as np
if len(sys.argv) >= 2:
  seed = int(sys.argv[1])
else:
  seed = 0
np.random.seed(seed)
import pandas as pd

GAMMADIR = "./Signal.csv"
HADRONDIR="./Background.csv"
COLUMNS = ["rec.nHit/U/1","rec.CxPE40/F/0.01","rec.PINC/F/0.01","rec.logNNEnergyV2/F/0","rec.disMax/F/0.01","rec.LDFAmp/F/0.01","rec.LDFChi2/F/0.01","rec.nChAvail/U/1","rec.nHitSP20/U/1"]
allset = pd.read_csv(GAMMADIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
hadron = pd.read_csv(HADRONDIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
#allset[:,0]=allset[:,0]/allset[:,1]
allset=allset[allset[:,1]!=0,:]
hadron=hadron[hadron[:,1]!=0,:]
allset[:,0]= np.log10(np.divide(allset[:,1], allset[:,0],out=np.zeros_like(allset[:,0]), where=allset[:,0]!=0))
allset=np.delete(allset,1,1)

allset[:,6]=np.divide(allset[:,7], allset[:,6],out=np.zeros_like(allset[:,6]), where=allset[:,6]!=0)
allset=np.delete(allset,7,1)
allset=np.hstack((np.ones((len(allset),1)),allset))
#allset=allset[:len(hadron),:] #balancing class
hadron[:,0]= np.log10(np.divide(hadron[:,1],hadron[:,0], out=np.zeros_like(hadron[:,0]), where=hadron[:,0]!=0))
hadron=np.delete(hadron,1,1)
hadron[:,6]=np.divide(hadron[:,7], hadron[:,6],out=np.zeros_like(hadron[:,6]), where=hadron[:,6]!=0)
hadron=np.delete(hadron,7,1)
hadron = np.hstack((np.zeros((len(hadron),1)),hadron))
hadron=hadron[:len(allset),:]
allset=np.vstack((allset,hadron))
'''
HADRONDIR="./Background.csv"
hadron = pd.read_csv(HADRONDIR, skipinitialspace=True,
                             skiprows=0, usecols=COLUMNS)[COLUMNS].as_matrix()
hadron=hadron[hadron[:,1]!=0,:]
hadron[:,0]= np.log10(np.divide(hadron[:,1],hadron[:,0], out=np.zeros_like(hadron[:,0]), where=hadron[:,0]!=0))
hadron=np.delete(hadron,1,1)
hadron = np.hstack((np.zeros((len(hadron),1)),hadron))
allset=np.vstack((allset,hadron))
'''
np.random.shuffle(allset)

np.save("data.npy",allset)
