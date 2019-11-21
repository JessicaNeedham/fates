#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  5 08:56:41 2019

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


param_mat1 = []
param_mat2 = []

n = 0


for i in range(1,101) :
    filename_in = os.path.join('twopfts', 'fates_params_smort_2pfts_%d.nc'  % (i))
    fin = nc.netcdf_file(filename_in)

    canopy_mort = fin.variables['fates_prescribed_mortality_canopy'][0]
    canopy_npp = fin.variables['fates_prescribed_npp_canopy'][0]
    
    ps = [canopy_mort, canopy_npp]
    
    param_mat1.append(ps)
    
    canopy_mort2 = fin.variables['fates_prescribed_mortality_canopy'][1]
    canopy_npp2 = fin.variables['fates_prescribed_npp_canopy'][1]

    ps2 = [canopy_mort2, canopy_npp2]
    
    param_mat2.append(ps2)
    
    
    
param_mat1 = np.array(param_mat1)
param_mat2 = np.array(param_mat2)

np.save('paramMatPFT1.npy', param_mat1)
np.save('paramMatPFT2.npy', param_mat2)
