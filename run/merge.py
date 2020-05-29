#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import math

sume = np.zeros(20,dtype=np.float32)
sumde2 = np.zeros(20,dtype=np.float32)
sumde = np.zeros(20,dtype=np.float32)

for fid in range(0,200):
   fname=str(fid)+".dat"
   print(fname)
   with open(fname) as f:
       while(True):
         line = f.readline()
         if(not line):
 	    break
         (tmp,i,e,de) = line.split(' ')
         print(tmp,i,e,de)
         idx = int(i)-1
         print('idx',idx,float(de))
         sume[idx] = sume[idx] + float(e)
         sumde2[idx] = sumde2[idx] + float(de)*float(de)

file = open('fout.dat','w') 
for j in range(0,20):
   sumde[j] = math.sqrt(sumde2[j])
#   print("Detector",j+1,sume[j],sumde[j])
   ss="Detector "+str(j+1)+" "+str(sume[j])+" "+str(sumde[j])+"\n"
   file.write(ss)
file.close()
