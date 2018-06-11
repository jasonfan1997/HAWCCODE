#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  9 23:45:10 2018

@author: jason
"""

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

result=np.loadtxt("standardcutq.csv",delimiter=',')
out="scut"
'''
new=np.zeros((120,4))
for i in range(len(result)):
    new[int((result[i,0]-1)*10+result[i,1]),:]=result[i,:4]
'''
#new[:,3]=np.power(new[:,3],10)
                       
#"rec.dec <= 0.41887902&& "          \
#"0.34906585 <=rec.dec  && "          \
#"rec.ra <= 1.49225651 && "          \
# "1.42244334 <= rec.ra && "           \
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
						   "rec.nChAvail >= 700 && "              \
                           "rec.zenithAngle < 0.785 && "              \
						   "rec.coreFiduScale <= 100.0&&"   \
                           "rec.PINC<={} &&" \
                           "log10(rec.CxPE40 / rec.nHitSP20)<={}\"\n".format( FHIT_WIDTH *          fbin[i], FHIT_WIDTH * (fbin[i] + 1), ENERGY_WIDTH * ebin[j] + ENERGY_OFFSET, ename, ename, ENERGY_WIDTH * (ebin[j] + 1) + ENERGY_OFFSET,result[t,5],result[t,6])
                t+=1           
                with open(out,"a") as myfile:
                    myfile.write(other_cut_string)
    
        
