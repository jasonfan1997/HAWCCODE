#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  6 22:19:33 2018

@author: jason
"""

import argparse
import os
import subprocess 
import sys
import numpy as np

ANGRE_WIDTH = 0.05
FHIT_WIDTH = 0.1
ENERGY_WIDTH = 0.25
ENERGY_OFFSET = 2.50
ename="rec.logNNEnergy"
parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+', help='path to the file')
parser.add_argument('out', nargs='+', help='path to the outfile')
args_namespace = parser.parse_args()
args = vars(args_namespace)['file'][0]
out= vars(args_namespace)['out'][0]
result=np.loadtxt(args,delimiter=',')
#os.remove(out)
#fbin=np.append(np.insert(np.linspace(0.1, 0.9, num=9),0,0.05),1)
#ebin=np.linspace(2.5, 5.5, num=13)
fbin = list(range(10))
ebin = list(range(12))
t=0

for i in range(len(fbin)):
    for j in range(len(ebin)):
                #bins+=1
                other_cut_string = str(t)+" \"{} * rec.nChAvail < rec.nHitSP20 && " \
						   "rec.nHitSP20 <= {} * rec.nChAvail && "\
						   "{} < {} && "             \
						   "{} <= {} && "            \
						   "rec.angleFitStatus == 0 && "          \
                           "rec.angleFitStatus == 0 && "          \
                           "rec.dec <= 0.41887902&& "          \
                           "0.34906585 <=rec.dec  && "          \
                           "rec.ra <= 1.49225651 && "          \
						   "1.42244334 <= rec.ra && "           \
						   "rec.nChAvail >= 700 && "              \
                           "rec.zenithAngle < 0.785 && "              \
						   "rec.coreFiduScale <= 100.0&&"   \
                           "{}<=rec.proba\"\n".format( FHIT_WIDTH *          fbin[i], FHIT_WIDTH * (fbin[i] + 1), ENERGY_WIDTH * ebin[i] + ENERGY_OFFSET, ename, ename, ENERGY_WIDTH * (ebin[i] + 1) + ENERGY_OFFSET,result[t,2])
                t+=1           
                with open(out,"a") as myfile:
                    myfile.write(other_cut_string)
    
        
