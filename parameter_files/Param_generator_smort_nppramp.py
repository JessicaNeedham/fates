#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 14:18:17 2019

@author: jessicaneedham
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



for i in range(1,101) :
    
    fin = os.path.join('twopfts', 'fates_params_smort_2pfts_%d.nc' % (i))
    fileout = os.path.join('twopfts', 'fates_params_nppramp_smort_2pfts_%d.nc'  % (i))
    

    # The first call to main we generate a new file
    var = 'fates_prescribed_npp_max'
    val = 0.25
    pft = 1
    modp.main(var = var, pft = pft, fin = fin, 
              val = val,  fout = fileout, O = 0)  
   
    # The now overwrite it
    pft =2
    modp.main(var = var, pft = pft, fin = fileout, 
              val = val,  fout = fileout, O = 1)
    

