# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 09:58:57 2018

@author: kwoklung
"""

import numpy as np
import matplotlib.pyplot as plt
q=np.loadtxt("r.csv",delimiter=',')
j=0
for i in range(10):
    plt.title('fbin:'+str(q[i*12,0]))
    plt.plot(q[i*12:i*12+12,2],color='r')
    plt.plot(q[i*12:i*12+12,3],color='b')
    plt.plot(q[i*12:i*12+12,4],color='g')
    plt.savefig('plot/fbin'+str(q[i*12,0]*10)+".png")
    plt.close()
    