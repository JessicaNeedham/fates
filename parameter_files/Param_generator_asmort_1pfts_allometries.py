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

d2h = np.linspace(0.05, 0.5, num=10)
d2bl2 = np.linspace(1.1, 1.4, num=10)


for i in range(0,10) :
    
    for j in range(0,10) :

        fileout = os.path.join('onepft_ensembles', 'fates_params_asmort_1pft_allom_%d.nc'  % (n))                                                                                                                                                          
        pft = 1
        fin = fileout


        # The first call to main we generate a new file                                                                                                                                         
        var = 'fates_allom_d2h2'
        modp.main(var = var, pft = pft, fin = 'onepft_ensembles/fates_params_asmort_1pft_allom_1.nc',
              val = d2h[i],  fout = fileout, O = 0)
   
        var = 'fates_allom_d2bl2'
        modp.main(var = var, pft = pft, fin = fin,
              val = d2bl2[j], fout = fileout, O = 1)
    
        n = n+1
        print(n)

     
