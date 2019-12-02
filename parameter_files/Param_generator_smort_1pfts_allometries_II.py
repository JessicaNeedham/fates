#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 15:55:23 2019

@author: JFNeedham
"""

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

d2h2 = np.linspace(0.55, 2.0, num=10)
d2bl_p2 = np.linspace(1.1,1.4, num=10)


for i in range(0,10) :
    
    for j in range(0,10) :

        fileout = os.path.join('onepft_ensembles', 'fates_params_smort_1pft_allom_v2_%d.nc'  % (n))                                                                                                                                                          
        pft = 1
        fin = fileout


        # The first call to main we generate a new file                                                                                                                                         
        var = 'fates_allom_d2h2'
        modp.main(var = var, pft = pft, fin = 'onepft_ensembles/fates_params_smort_1pft_allom_v2_1.nc',
              val = d2h2[i],  fout = fileout, O = 0)
   
        var = 'fates_allom_d2bl2'
        modp.main(var = var, pft = pft, fin = fin,
              val = d2bl_p2[j], fout = fileout, O = 1)
    
        n = n+1
        print(n)

     

