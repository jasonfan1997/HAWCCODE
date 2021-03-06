#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  5 00:48:41 2018

@author: jason
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 10 16:20:13 2018

@author: user98
"""

import numpy as np
import argparse
import pickle #pickle
from xcdf import XCDFFile
import ntpath
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

import os
for prefix in ('OMP', 'MKL', 'NUMEXPR'):
    os.environ['%s_NUM_THREADS' %prefix] = '1'
#deg = np.pi/180.
parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+', help='path to the file')
parser.add_argument('out', nargs='+', help='path to the outfile')
args_namespace = parser.parse_args()
args = vars(args_namespace)['file'][0]
outdir= vars(args_namespace)['out'][0]
xf = XCDFFile(args)
name=path_leaf(args)
#xf = XCDFFile("reco_run006657_00001.xcd")
print args[0]
allset = []
#for record in xf.fields("rec.nHit,rec.CxPE40,rec.PINC,rec.logNNEnergy,rec.disMax,rec.LDFAge,rec.LDFAmp,rec.LDFChi2,rec.nChAvail,rec.nHitSP20"):
#for record in xf.fields("rec.nHit,rec.CxPE40,rec.PINC,rec.disMax,rec.LDFAge,rec.LDFAmp,rec.LDFChi2"):
for record in xf.fields("rec.nHit,rec.CxPE40,rec.PINC,rec.logNNEnergyV2,rec.disMax,rec.LDFAmp,rec.LDFChi2,rec.nChAvail,rec.nHitSP20"):
    allset.append(record)

allset=np.array(allset)
#allset=allset[allset[:,1]!=0,:]
allset[:,0]= np.log10(np.divide(allset[:,1], allset[:,0],out=np.zeros_like(allset[:,0]), where=allset[:,0]!=0))
allset[allset==-np.inf]=0
allset=np.delete(allset,1,1)
allset[:,6]=allset[:,7]/allset[:,6]
allset[allset==np.inf]=0
allset=np.delete(allset,7,1)



with open('clfMCv7v3.pickle', 'rb') as f:
    clf = pickle.load(f)

clf.n_jobs=1
args=path_leaf(args)
name=args.replace("xcd","csv")
predicted_prob=clf.predict_proba(allset)
np.savetxt(outdir+name,predicted_prob[:,1],delimiter=',', fmt='%f')
'''
xf = XCDFFile("reco_run006657_00001.xcd")

for record in xf.fields("rec.nHit,rec.CxPE40,rec.PINC,rec.logNNEnergy,rec.disMax,rec.LDFAge,rec.LDFAmp,rec.nChAvail,rec.nHitSP20"):
    allset = np.vstack((allset, np.array(record)))
'''
