#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:10:24 2019

@author: jessicaneedham
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 16:08:44 2019

@author: JFNeedham
"""
import numpy as np
import matplotlib.pyplot as plt
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

while (n < 100):
    canopy_mort = np.random.uniform(0.005, 0.06)
    under_mort = np.random.uniform(0.01, 0.1)
    canopy_npp = np.random.uniform(0.5, 2.0)
    under_npp = np.random.uniform(0.05, 0.2)
    recrt = np.random.uniform(0.125, 0.6)
    ip = np.random.uniform(100,300)
    
    canopy_mort2 = np.random.uniform(0.005, 0.06)
    under_mort2 = np.random.uniform(0.01, 0.1)
    canopy_npp2 = np.random.uniform(0.5, 2.0)
    under_npp2 = np.random.uniform(0.05, 0.2)
    recrt2 = np.random.uniform(0.125, 0.6)
    ip2 = np.random.uniform(100,300)
    
    if canopy_mort < canopy_mort2 :
        if under_mort < under_mort2 :
            if canopy_npp < canopy_npp2 :
                if under_npp < under_npp2 :
                    if recrt < recrt2 : 
                        
                        p1 = [canopy_mort,under_mort,canopy_npp,under_npp,recrt,ip]
                        p2 = [canopy_mort2,under_mort2,canopy_npp2,under_npp2,recrt2,ip2]
                        
                        param_mat1.append(p1)
                        param_mat2.append(p2)
                        
                        n = n + 1
                        
                        
# convert lists to arrays                        
param_mat1 = np.array(param_mat1)
param_mat2 = np.array(param_mat2)

# plot them - make sure they look good
fig, ax = plt.subplots()
ax.scatter(param_mat1[:,1], param_mat1[:,3], label = 'PFT 1')
ax.scatter(param_mat2[:,1], param_mat2[:,3], color = 'red', label = 'PFT 2')
ax.set_xlabel('Canopy mortality', fontsize=15)
ax.set_ylabel('Canopy NPP', fontsize=15)
ax.set_title('Parameter space')
ax.plot([param_mat1[:,1], param_mat2[:,1]], [param_mat1[:,3], param_mat2[:,3]], 'lightgrey')
ax.legend()
plt.savefig("Param_space_smort.png")


# save them 
paramMatPFT1 = TemporaryFile()
paramMatPFT2 = TemporaryFile()
np.save(paramMatPFT1, param_mat1)
np.save(paramMatPFT2, param_mat2)

print(param_mat1.shape[0])

for i in range(0,param_mat1.shape[0]) :
    
    fileout = os.path.join('twopfts', 'fates_params_smort_2pfts_%d.nc'  % (i))
    
    # over write pft 1 variables
    pft = 1
    fin = fileout
    
    
    # The first call to main we generate a new file
    var = 'fates_prescribed_mortality_canopy'
    modp.main(var = var, pft = pft, fin = 'fates_params_smort_2pfts.nc', 
              val = param_mat1[i,0],  fout = fileout, O = 0)  
    # subsequent calls we overwrite that file
    var = 'fates_prescribed_mortality_understory'
    modp.main(var = var, pft = pft, fin = fin,
              val = param_mat1[i,1], fout = fileout, O = 1)  
    
    var = 'fates_prescribed_npp_canopy'
    modp.main(var = var, pft = pft, fin = fin,
              val = param_mat1[i,2], fout = fileout, O = 1)
    
    var = 'fates_prescribed_npp_understory'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat1[i,3], fout = fileout, O = 1) 
    
    var = 'fates_prescribed_recruitment'
    modp.main(var = var, pft = pft, fin = fin,  
              val = param_mat1[i,4], fout = fileout, O = 1)  
   
    var = 'fates_mort_ip_senescence'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat1[i,5], fout = fileout, O = 1)
    
    # overwrite pft2 variables
    pft = 2
    var = 'fates_prescribed_mortality_canopy'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,0], fout = fileout, O = 1)  
    
    var = 'fates_prescribed_mortality_understory'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,1], fout = fileout, O = 1)  
    
    var = 'fates_prescribed_npp_canopy'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,2], fout = fileout, O = 1)
    
    var = 'fates_prescribed_npp_understory'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,3], fout = fileout, O = 1) 
    
    var = 'fates_prescribed_recruitment'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,4], fout = fileout, O = 1)  
   
    var = 'fates_mort_ip_senescence'
    modp.main(var = var, pft = pft, fin = fin, 
              val = param_mat2[i,5], fout = fileout, O = 1)
  
    
    
