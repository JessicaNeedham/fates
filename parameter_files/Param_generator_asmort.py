#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:08:44 2019

@author: JFNeedham
"""
import numpy as np
#import matplotlib.pyplot as plt
from tempfile import TemporaryFile
import modp
import os
from scipy.io import netcdf as nc
import argparse
import shutil
import tempfile
import sys
import datetime
import time


### 
n = 0  

# matrix for holding pairs that fall along the growth survival trade-off
param_mat1 = []
param_mat2 = []

while (n < 101):
    canopy_mort = np.random.uniform(0.005, 0.06)
    canopy_npp = np.random.uniform(0.4, 2.0)
    
    canopy_mort2 = np.random.uniform(0.005, 0.06)
    canopy_npp2 = np.random.uniform(0.4, 2.0)
    
    
    if canopy_mort < canopy_mort2 :
         if canopy_npp < canopy_npp2 :
            
                        p1 = [canopy_mort,canopy_npp]
                        p2 = [canopy_mort2,canopy_npp2]
                        
                        param_mat1.append(p1)
                        param_mat2.append(p2)
                        
                        n = n + 1
                        
                        
# convert lists to arrays                        
param_mat1 = np.array(param_mat1)
param_mat2 = np.array(param_mat2)

# plot them - make sure they look good
#fig, ax = plt.subplots()
#ax.scatter(param_mat1[:,1], param_mat1[:,3], label = 'PFT 1')
#ax.scatter(param_mat2[:,1], param_mat2[:,3], color = 'red', label = 'PFT 2')
#ax.set_xlabel('Canopy mortality', fontsize=15)
#ax.set_ylabel('Canopy NPP', fontsize=15)
#ax.set_title('Parameter space')
#ax.plot([param_mat1[:,1], param_mat2[:,1]], [param_mat1[:,3], param_mat2[:,3]], 'lightgrey')
#ax.legend()
#plt.savefig("Param_space.png")


# save them 

np.save('asmort_paramMatPFT1', param_mat1)
np.save('asmort_paramMatPFT2', param_mat2)



for i in range(1,101) :
    
    fileout = os.path.join('twopfts', 'fates_params_asmort_2pfts_%d.nc'  % (i))
    
    # over write pft 1 variables
    pft = 1
    fin = fileout
    
    
    # The first call to main we generate a new file
    var = 'fates_prescribed_mortality_canopy'
    modp.main(var = var, pft = pft, fin = 'fates_params_asmort_2pfts.nc', 
              val = param_mat1[i,0],  fout = fileout, O = 0)  
    # subsequent calls we overwrite that file
    var = 'fates_prescribed_npp_canopy'
    modp.main(var = var, pft = pft, fin = fin,
              val = param_mat1[i,1], fout = fileout, O = 1)
    
    # overwrite pft2 variables
    pft = 2
    var = 'fates_prescribed_mortality_canopy'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,0], fout = fileout, O = 1)  
      
    var = 'fates_prescribed_npp_canopy'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,1], fout = fileout, O = 1)
    
  
    
    
