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


import os
for prefix in ('OMP', 'MKL', 'NUMEXPR'):
    os.environ['%s_NUM_THREADS' %prefix] = '1'
#deg = np.pi/180.
"""
parser = argparse.ArgumentParser()
parser.add_argument('file', nargs='+', help='path to the file')
parser.add_argument('num', nargs='+', help='path to the file')
args_namespace = parser.parse_args()
args = vars(args_namespace)['file'][0]
number = vars(args_namespace)['num'][0]
"""
args="temp.txt"
number=3
result=np.loadtxt(args,delimiter=" ",dtype="str")
np.random.shuffle(result)
goal=np.mean(result[:,3].astype(np.float))*number
temp=0
sumrow=0
last=0
for i in range(len(result)):
    sumrow+=float(result[i,3])
    if sumrow >=goal:
        np.savetxt("node"+str(temp),result[last:i,:3],fmt="%s") 
        sumrow=0
        temp+=1
        last=i
    