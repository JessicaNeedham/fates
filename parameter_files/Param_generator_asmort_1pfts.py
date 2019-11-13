#!/usr/bin/env python3                                                                                                                                                                      
# -*- coding: utf-8 -*-                                                                                                                                                                     
"""                                                                                                                                                                                         
Created on Wed Sep 11 16:08:44 2019                                                                                                                                                         
                                                                                                                                                                                            
@author: JFNeedham                                                                                                                                                                          
"""
import numpy as np
#import matplotlib.pyplot as plt                                                                                                                                                            
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
n = 1

canopy_mort = np.linspace(0.005, 0.05, num=10)
ip = np.linspace(180,360, num=10)
print(canopy_mort)
print(ip)


for i in range(0,10) :
    
    for j in range(0,10) :

        fileout = os.path.join('onepft_ensembles', 'fates_params_asmort_1pft_%d.nc'  % (n))                                                                                                                                                          
        pft = 1
        fin = fileout


        # The first call to main we generate a new file                                                                                                                                         
        var = 'fates_prescribed_mortality_canopy'
        modp.main(var = var, pft = pft, fin = 'onepft_ensembles/fates_params_asmort_1pft_1.nc',
              val = canopy_mort[i],  fout = fileout, O = 0)
   
        var = 'fates_mort_ip_age_senescence'
        modp.main(var = var, pft = pft, fin = fin,
              val = ip[j], fout = fileout, O = 1)
    
        n = n+1
        print(n)

     
