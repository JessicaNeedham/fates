            #!/usr/bin/env python
            
            #### this script modifies a FATES parameter file. It accepts the following flags
            # --var or --variable: variable.
            # --pft or --PFT: PFT number. If this is missing, script will assume that its a global variable that is being modified.
            # --input or --fin: input filename.
            # --output or --fout: output filename.  If missing, will assume its directly modifying the input file, and will prompt unless -O is specified
            # --O or --overwrite: overwrite output file without asking.
            # --value or --val: value to put in variable
            # --s or --silent: don't write anything on successful execution.
            ####
            #
            # Written by C. Koven, 2018
            #
            
            # =======================================================================================
            # =======================================================================================
        
import os
from scipy.io import netcdf as nc
import shutil
import tempfile
import numpy as np
            
            # ========================================================================================
            # ========================================================================================
            #                                        Main
            # ========================================================================================
            # ========================================================================================
            
def main(var, pft, fin, val, fout, O):
            
    # work with the file in some random temporary place so that if something goes wrong, then nothing happens to original file and it doesn't make a persistent output file
    tempdir = tempfile.mkdtemp()
    tempfilename = os.path.join(tempdir, 'temp_fates_param_file.nc')
    ncfile_old = None
    outputfname = fout
       
    try:
        outputval = float(val)
        
    except:
        try:
            print('output variable not interpretable as real. trying array')
            outputval = np.fromstring(val, sep=',', dtype=np.float32)
            if len(outputval) == 0:
                raise RuntimeError('output variable needs to have size greater than zero')
        except:
            raise RuntimeError('output variable not interpretable as real or array')
    #
    #
     
    try:
        shutil.copyfile(fin, tempfilename)
                    #
        ncfile = nc.netcdf_file(tempfilename, 'a')
                    #
                  
        var = ncfile.variables[var]
        #print(var)
                   #
        ### check to make sure that, if a PFT is specified, the variable has a PFT dimension, and if not, then it doesn't. and also that shape is reasonable.
        ndim_file = len(var.dimensions)
        ispftvar = False
        # for purposes of current state of this script, assume 1D 
        if ndim_file > 2:
            raise ValueError('variable dimensionality is too high for this script')
        for i in range(ndim_file):
            if var.dimensions[i] == 'fates_pft':
                ispftvar = True
                npft_file = var.shape[i]
                pftdim = i
                otherdimpresent = False
            elif var.dimensions[i] in ['fates_history_age_bins','fates_history_size_bins','fates_history_height_bins','fates_NCWD','fates_litterclass','fates_leafage_class','fates_prt_organs','fates_hydr_organs','fates_variants']:
                otherdimpresent = True
                otherdimname = var.dimensions[i]
                otherdimlength = var.shape[i]
            else:
                raise ValueError('variable is not on either the PFT or scalar dimension')
        #
        
        if pft != None and ispftvar:
            if pft > npft_file:
                raise ValueError('PFT specified ('+str(pft)+') is larger than the number of PFTs in the file ('+str(npft_file)+').')
            if pftdim == 0:
                var[pft-1] = outputval
            if pftdim == 1:
                var[:,pft-1] = outputval
        elif pft == None and not ispftvar and ndim_file > 0:
            if not otherdimpresent:
                var[:] = outputval
            else:
                var[:] = outputval
        elif ndim_file < 1:
            var.assignValue(outputval)
        else:
            raise ValueError('Nothing happened somehow.')
        #           
        ncfile.close()
        if type(ncfile_old) != type(None):
            ncfile_old.close()
        #
        #
        # now move file from temporary location to final location
        if O == 1:   
            os.remove(outputfname)
        shutil.move(tempfilename, outputfname)
        shutil.rmtree(tempdir, ignore_errors=True)
    except:
        shutil.rmtree(tempdir, ignore_errors=True)
        raise
            
   # =======================================================================================
   # This is the actual call to main
               
if __name__ == "__main__":
    main()