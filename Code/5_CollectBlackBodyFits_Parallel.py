import pandas as pd 
import numpy as np 
import pandas as pd
import numpy as np 
import glob
import re
from inlist import *
import os
from multiprocessing.dummy import Pool as ThreadPool
from itertools import permutations
import time



def collect_blackbodyfits(i):
    """
        collect the blackbody fitting results.
    """
    temp = pd.DataFrame(index=[0], columns=columns)
    starID = os.path.basename(files[i]).split(".")[0]
    data = pd.read_csv(files[i], dtype={"starID":str})
    temp.loc[0, "starID"]  = starID
    if i%100==0: print("#{:d}th target:  {:s}".format(i, starID))
    temp.iloc[0, 1:] = data.iloc[0, 1:].copy()
    return temp.values[0]


def combine_parallel(numb=0, columns=None, func=None):
    """
        Combine the collected results and convert to a pandas DataFrame.
    """
    pool = ThreadPool(10)
    t0 = time.time()
    results = pool.map(func, range(numb))
    t1 = time.time()
    print("Collecting data for #{:} target, {:.2f} seconds used.".format(numb, t1-t0))
    target = pd.DataFrame(columns=columns, data=results)
    return target


# delete txt files that generated when cleaning old batch files.
files_txt = glob.glob(save_blackbody_flux_path_single+"*.txt")
if len(files_txt)>0: os.system("rm "+glob.glob(save_blackbody_flux_path_single+"*.txt")[0])
files_txt = glob.glob(save_blackbody_flux_path_batch+"*.txt")
if len(files_txt)>0: os.system("rm "+glob.glob(save_blackbody_flux_path_batch+"*.txt")[0])



# combine all the fitting results.
if Lbatch:
    # combine the fitting results for all the stars.
    files = np.array(glob.glob(save_blackbody_flux_path_batch+"*.csv"))
    columns = ["starID","Qflg","Av","ScaleF","ScaleFerr","Teff","Tefferr","RedChi"]
    target = pd.DataFrame()
    for sub, file in enumerate(files):
        file_tmp = pd.read_csv(file, dtype={"starID":str})
        target = pd.concat([target, file_tmp], axis=0 , ignore_index=True)
    # save bands used for each star.
    target.drop_duplicates(subset=['starID'], inplace=True)
    target = target[columns]
    target = target.sort_values(by=["starID"]).reset_index(drop=True)
    
else:
    # combine the fitting results for all the stars.
    files = np.array(glob.glob(save_blackbody_flux_path_single+"*.csv"))
    columns = ["starID","Qflg","Av","ScaleF","ScaleFerr","Teff","Tefferr","RedChi"]    
    target = combine_parallel(numb=len(files), columns=columns, func=collect_blackbodyfits)

# set digits.
target['Av'] = target['Av'].map(lambda x: '{0:.2f}'.format(x)) 
target['ScaleF'] = target['ScaleF'].map(lambda x: '{0:.3e}'.format(x)) 
target['ScaleFerr'] = target['ScaleFerr'].map(lambda x: '{0:.3e}'.format(x)) 
target['Teff'] = target['Teff'].map(lambda x: '{0:.3e}'.format(x)) 
target['Tefferr'] = target['Tefferr'].map(lambda x: '{0:.3e}'.format(x)) 
target['RedChi'] = target['RedChi'].map(lambda x: '{0:.3e}'.format(x)) 

# ensure that very same line of the tables of stellar parameters, input photometry, and this blackbody fits are consistent.
input = pd.read_csv(file_input_params, dtype={"starID":str})
table = pd.read_csv(file_input_photometry, dtype={"starID":str})
if (input.starID==table.starID).all(): 
    print("The tables of stellar parameters and input photometry are consistent")
    target = pd.merge(table[["starID"]], target, on="starID").reset_index(drop=True).copy()
    if len(input)==len(target): 
        if (input.starID==target.starID).all():
            print("The tables of stellar parameters and blackbody fits are consistent")
        else:
            print("The tables of stellar parameters and blackbody fits are inconsistent. Check!!!")
            os.system("exit")    
    else:
        print("the number of stars with valid blackbody fits does not equal the number of stars in the input list")               
else:
    os.system("exit")    

# output
target.to_csv(Combined_BlackBody_Fits_path+"Output_BlackBody_Fits.csv", index=False,float_format="%.4e") 


