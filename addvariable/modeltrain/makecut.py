#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 00:05:22 2018

@author: jason
"""


out="mcrdthres"
name="mcrdthres.csv"
result=np.loadtxt(name,delimiter=',')
ENERGY_WIDTH = 0.25
ENERGY_OFFSET = 2.50
ename="rec.logNNEnergyV2"
#fbin = list(range(10))
ebin = list(range(12))
t=0
fbin=np.array([0,0.067,0.105,0.162,0.247,0.356,0.485,0.618,0.74,0.84,1.01])
for i in range(len(fbin)-1):
    for j in range(len(ebin)):
                if result[t,3]!=0.5:
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
                    with open(out,"a") as myfile:
                        myfile.write(other_cut_string)
                t+=1           
