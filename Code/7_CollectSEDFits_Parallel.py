import pandas as pd
import numpy as np 
import glob
import re
from inlist import *
import SedFitFunctions
from uncertainties import unumpy 
import os 
from multiprocessing.dummy import Pool as ThreadPool
from itertools import permutations
import time



def collect_sed(i):
    """
        collect the SED fitting results obtained in the interpolation mode.
    """
    starID = os.path.basename(files[i]).split(".")[0]
    data = pd.read_csv(files[i], dtype={"starID":str})
    if i%100==0: print("#{:d}th target (interpolation):  {:s}".format(i, starID))
    return data.values[0]


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



################################################################################################
modelParams = pd.read_csv(grid_list_path+"PassbandFlux_"+grid+".csv").rename(columns={"modelindex":"ModelIndex"}) 
filters = pd.read_csv(filter_zero_point_path)
if len(glob.glob(single_data_sedfit_path+"*.txt"))>0: os.system("rm "+glob.glob(single_data_sedfit_path+"*.txt")[0])
if len(glob.glob(batch_data_sedfit_path+"*.txt"))>0: os.system("rm "+glob.glob(batch_data_sedfit_path+"*.txt")[0])


################################################################################################
# collect results from the SED fitting without the model interpolation.
print("Collecting data from the fitting without interpolation")
# files = glob.glob(save_sedfit_flux_path+"*.BestfitModels.csv")
if Lbatch:
    files = glob.glob(batch_data_sedfit_path+"*.BestfitModels.csv")
    target = pd.DataFrame()
    for sub, file in enumerate(files):
        try: 
            file_tmp = pd.read_csv(file, dtype={"starID":str})
            target = pd.concat([target, file_tmp], axis=0, ignore_index=True)
        except:
            print("stars are skipped due to input stellar parameters are off the grid")      
else:
    files = glob.glob(single_data_sedfit_path+"*.BestfitModels.csv")
    # parallel the loop for colloecting the SED fitting result (no intermpolation)
    columns = pd.read_csv(files[0]).columns
    target = combine_parallel(numb=len(files), columns=columns, func=collect_sed)


# select model stellar parameters.
target = pd.merge(target, modelParams[["teff", "logg", "feh", "ModelIndex"]], on="ModelIndex", how="left")
target.rename(columns={"teff":"modelteff", "logg":"modellogg", "feh":"modelfeh"}, inplace=True)
columns = ["starID", "Av", "Averr", "ScaleF", "ScaleFerr", "Redchi", "Fbol", "ModelIndex", "RelDiff", "npoint", "JumpBy", "modelteff", "modellogg", "modelfeh"]
target = target[columns]
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.drop_duplicates(subset=['starID'], inplace=True)

# load gaia distance.        
params = pd.read_csv(file_input_params, dtype={"starID":str})[["starID", "d", "derr", "teff", "tefferr", "logg", "loggerr", "feh", "feherr"]]
target = pd.merge(target, params, on="starID", how="left")
# calculate bolometric flux uncertainty.
target.loc[:, "Fbolerr"] = np.log(10)*target.ScaleFerr.copy()*target.Fbol.copy()
# calculate angular radius and its unceritainty in units of mas
unity2mas = (3.6*1.8/np.pi)*1e8 # convert 1 radian to 1 mas
target.loc[:, "angRad"] = 10**(0.5*target.ScaleF.copy()) * unity2mas # demensionless unit.
target.loc[:, "angRaderr"] = target.loc[:, "angRad"] * np.log(10) * 0.5*target.ScaleFerr.copy() 
# calculate log luminosity from bolometric flux and distance
target.loc[:, "lumi"], target.loc[:, "lumierr"] = \
SedFitFunctions.CalcLumi(dist=target.d, disterr=target.derr, Fbolo=target.Fbol, Fboloerr=target.Fbolerr)
# calculate radius and its uncertainty. Note that in this case, radius is derived from angular radius and diatance. 
target.loc[:, "radius"], target.loc[:, "radiuserr"] = \
SedFitFunctions.radius_From_angRadius(angRadius=target.angRad, angRadiuserr=target.angRaderr, dist=target.d, disterr=target.derr)
# sort the columns.
columns = ['starID', 'modelteff', 'modellogg', 'modelfeh', 'teff', 'tefferr', 'logg', 'loggerr', 'feh',
'feherr', 'd', 'derr', 'Av', 'Averr', 'ScaleF', 'ScaleFerr', 'Fbol', 'Fbolerr', 'angRad', 'angRaderr', 
'lumi', 'lumierr', 'radius', 'radiuserr', 'ModelIndex', 'Redchi', 'RelDiff', 'npoint', 'JumpBy']
target = target[columns]
# specify the format of each column.
formats = {'starID': '{:s}', 'modelteff': '{:.3f}', 'modellogg': '{:.3f}', 'modelfeh': '{:.3f}',
'teff': '{:.3f}', 'tefferr': '{:.3f}', 'logg': '{:.3f}', 'loggerr': '{:.3f}', 'feh': '{:.3f}',
'feherr': '{:.3f}', 'd': '{:.5f}', 'derr':'{:.5f}', 'Av':'{:.5f}', 'Averr':'{:.5f}', 'ScaleF':'{:.5f}', 
'ScaleFerr':'{:.5f}', 'Fbol':'{:.5e}', 'Fbolerr':'{:.5e}', 'angRad':'{:.5f}', 'angRaderr':'{:.5f}', 
'lumi':'{:.5f}', 'lumierr':'{:.5f}', 'radius':'{:.5f}', 'radiuserr':'{:.5f}', 'ModelIndex':'{:d}', 
'Redchi':'{:.5f}', 'RelDiff':'{:.5f}', 'npoint':'{:d}', 'JumpBy':'{:s}'}
for col, f in formats.items(): target[col] = target[col].map(lambda x: f.format(x))
# output the final table
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.to_csv(combined_data_sedfit_path+"Output_SED_Fits_No_Interp.csv", index=False) 




################################################################################################
# collect results from the SED fitting with the model interpolation.
print("Collecting data from the fitting with interpolation")
if Lbatch:
    files = glob.glob(batch_data_sedfit_path+"*.interpBestModel.csv")
    target = pd.DataFrame()
    for sub, file in enumerate(files):
        try: 
            file_tmp = pd.read_csv(file, dtype={"starID":str})
            target = pd.concat([target, file_tmp], axis=0, ignore_index=True)
        except:
            print("stars are skipped due to input stellar parameters are off the grid")    
else:
    # parallel the loop for colloecting the SED fitting result (no intermpolation)
    files = glob.glob(single_data_sedfit_path+"*.interpBestModel.csv")
    columns = pd.read_csv(files[0]).columns
    target = combine_parallel(numb=len(files), columns=columns, func=collect_sed)
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.drop_duplicates(subset=['starID'], inplace=True)

###################################################################
######### (1) for mode estimates ##################################
###################################################################
target_all = target.copy()
columns = ['starID', 'Redchi_mod', 'teff_mod', 'logg_mod', 'Av_mod', 'Averr_mod', 'ScaleF_mod', 'ScaleFerr_mod', 
'fbol_mod', 'fbolerr_mod']
target = target_all[columns].copy()
# columns = target[columns].columns.str.split("_mod").str[0]
target.columns = target[columns].columns.str.split("_mod").str[0]
target.rename(columns={"teff":"modelteff", "logg":"modellogg"}, inplace=True)    
# load gaia distance.        
columns = ["starID", "d", "derr", "teff", "tefferr", "logg", "loggerr", "feh", "feherr"]
params = pd.read_csv(file_input_params, usecols=columns, dtype={"starID":str})
target = pd.merge(target, params, on="starID", how="left")
# calculate angular radius and its unceritainty in units of mas
unity2mas = (3.6*1.8/np.pi)*1e8
target.loc[:, "angRad"] = 10**(0.5*target.ScaleF.copy()) * unity2mas # demensionless unit.
target.loc[:, "angRaderr"] = target.loc[:, "angRad"] * np.log(10) * 0.5*target.ScaleFerr.copy() 
# calculate radius and its uncertainty. Note that in this case, radius is derived from angular radius and diatance. 
target.loc[:, "radius"], target.loc[:, "radiuserr"] = \
SedFitFunctions.radius_From_angRadius(angRadius=target.angRad, angRadiuserr=target.angRaderr, dist=target.d, disterr=target.derr)
# calculate log luminosity from teff and radius. 
x = unumpy.uarray(target.teff, target.tefferr)
y = unumpy.uarray(target.radius, target.radiuserr)
z = unumpy.log10((x/5772)**4 * y**2)
target.loc[:, "lumi"] = unumpy.nominal_values(z)
target.loc[:, "lumierr"] = unumpy.std_devs(z)
# convert bolometric fluxed from the log scale to the linear scale.
target.loc[:, "fbol"] = 10**target.loc[:, "fbol"]
target.loc[:, "fbolerr"] = target.loc[:, "fbol"] * np.log(10) * target.loc[:, "fbolerr"]
target.rename(columns={"fbol":"Fbol", "fbolerr":"Fbolerr"}, inplace=True)
# calculate log luminosity from bolometric flux and distance. 
target.loc[:, "lumibol"], target.loc[:, "lumibolerr"] = \
SedFitFunctions.CalcLumi(Fbolo=target.Fbol, Fboloerr=target.Fbolerr, dist=target.d, disterr=target.derr)
# calculate angular radius and radius from bolometric flux and temperature.
target.loc[:, "angRadbol"] =  np.sqrt(target.Fbol/5.670374e-5/target.teff**4) * unity2mas 
target.loc[:, "angRadbolerr"] = np.sqrt(0.25*(target.Fbolerr/target.Fbol)**2+4*(target.tefferr/target.teff)**2) * target.loc[:, "angRadbol"]
target.loc[:, "radiusbol"], target.loc[:, "radiusbolerr"] = \
SedFitFunctions.radius_From_angRadius(angRadius=target.angRadbol, angRadiuserr=target.angRadbolerr, dist=target.d, disterr=target.derr)
# sort the columns.
columns = ['starID', 'modelteff', 'modellogg', 'teff', 'tefferr', 'logg', 'loggerr', 'feh',
'feherr', 'd', 'derr', 'Av', 'Averr', 'ScaleF', 'ScaleFerr', 'Fbol', 'Fbolerr', 'angRad', 'angRaderr', 
'lumi', 'lumierr', 'lumibol', 'lumibolerr', 'radius', 'radiuserr', 'angRadbol', 'angRadbolerr', 'radiusbol', 'radiusbolerr', 'Redchi']
target = target[columns]
# specify the format of each column.
formats = {'starID': '{:s}', 'modelteff': '{:.3f}', 'modellogg': '{:.3f}', 
'teff': '{:.3f}', 'tefferr': '{:.3f}', 'logg': '{:.3f}', 'loggerr': '{:.3f}', 'feh': '{:.3f}',
'feherr': '{:.3f}', 'd': '{:.5f}', 'derr':'{:.5f}', 'Av':'{:.5f}', 'Averr':'{:.5f}', 'ScaleF':'{:.5f}', 
'ScaleFerr':'{:.5f}', 'Fbol':'{:.5e}', 'Fbolerr':'{:.5e}', 'angRad':'{:.5f}', 'angRaderr':'{:.5f}', 
'lumi':'{:.5f}', 'lumierr':'{:.5f}', 'lumibol':'{:.5f}', 'lumibolerr':'{:.5f}', 'radius':'{:.5f}', 
'radiuserr':'{:.5f}', 'angRadbol':'{:.5f}', 'angRadbolerr':'{:.5f}', 'radiusbol':'{:.5f}', 'radiusbolerr':'{:.5f}', 'Redchi':'{:.5f}'}
for col, f in formats.items(): target[col] = target[col].map(lambda x: f.format(x))
# output the final table
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.to_csv(combined_data_sedfit_path+"Output_SED_Fits_Interp_MOD.csv", index=False) 


##################################################################
############# (2) Gaussian fit ###################################
##################################################################
columns = ['starID', 'teff', 'tefferr', 'logg', 'loggerr', 'Av', 'Averr', 'ScaleF', 'ScaleFerr', 'fbol', 'fbolerr']
target = target_all[columns].copy()
target.rename(columns={"teff":"modelteff", "tefferr":"modeltefferr", "logg":"modellogg", "loggerr":"modelloggerr"}, inplace=True)    
# load gaia distance.        
columns = ["starID", "d", "derr", "teff", "tefferr", "logg", "loggerr", "feh", "feherr"]
params = pd.read_csv(file_input_params, usecols=columns, dtype={"starID":str})
target = pd.merge(target, params, on="starID", how="left")
# calculate angular radius and its unceritainty in units of mas
unity2mas = (3.6*1.8/np.pi)*1e8
target.loc[:, "angRad"] = 10**(0.5*target.ScaleF.copy()) * unity2mas # demensionless unit.
target.loc[:, "angRaderr"] = target.loc[:, "angRad"] * np.log(10) * 0.5*target.ScaleFerr.copy() 
# calculate radius and its uncertainty. Note that in this case, radius is derived from angular radius and diatance. 
target.loc[:, "radius"], target.loc[:, "radiuserr"] = \
SedFitFunctions.radius_From_angRadius(angRadius=target.angRad, angRadiuserr=target.angRaderr, dist=target.d, disterr=target.derr)
# calculate log luminosity from teff and radius. 
x = unumpy.uarray(target.teff, target.tefferr)
y = unumpy.uarray(target.radius, target.radiuserr)
z = unumpy.log10((x/5772)**4 * y**2)
target.loc[:, "lumi"] = unumpy.nominal_values(z)
target.loc[:, "lumierr"] = unumpy.std_devs(z)
# convert bolometric fluxed from the log scale to the linear scale.
target.loc[:, "fbol"] = 10**target.loc[:, "fbol"]
target.loc[:, "fbolerr"] = target.loc[:, "fbol"] * np.log(10) * target.loc[:, "fbolerr"]
target.rename(columns={"fbol":"Fbol", "fbolerr":"Fbolerr"}, inplace=True)
# calculate log luminosity from bolometric flux and distance. 
target.loc[:, "lumibol"], target.loc[:, "lumibolerr"] = \
SedFitFunctions.CalcLumi(Fbolo=target.Fbol, Fboloerr=target.Fbolerr, dist=target.d, disterr=target.derr)
# calculate angular radius and radius from bolometric flux and temperature.
target.loc[:, "angRadbol"] =  np.sqrt(target.Fbol/5.670374e-5/target.teff**4) * unity2mas
target.loc[:, "angRadbolerr"] = np.sqrt(0.25*(target.Fbolerr/target.Fbol)**2+4*(target.tefferr/target.teff)**2) * target.loc[:, "angRadbol"]
target.loc[:, "radiusbol"], target.loc[:, "radiusbolerr"] = \
SedFitFunctions.radius_From_angRadius(angRadius=target.angRadbol, angRadiuserr=target.angRadbolerr, dist=target.d, disterr=target.derr)
# sort the columns.
columns = ['starID', 'modelteff', 'modeltefferr', 'modellogg', 'modelloggerr', 'teff', 'tefferr', 'logg', 'loggerr', 'feh',
'feherr', 'd', 'derr', 'Av', 'Averr', 'ScaleF', 'ScaleFerr', 'Fbol', 'Fbolerr', 'angRad', 'angRaderr', 
'lumi', 'lumierr', 'lumibol', 'lumibolerr', 'radius', 'radiuserr', 'angRadbol', 'angRadbolerr', 'radiusbol', 'radiusbolerr']
target = target[columns]
# specify the format of each column.
formats = {'starID': '{:s}', 'modelteff': '{:.3f}', 'modeltefferr': '{:.3f}', 'modellogg': '{:.3f}', 'modelloggerr': '{:.3f}',
'teff': '{:.3f}', 'tefferr': '{:.3f}', 'logg': '{:.3f}', 'loggerr': '{:.3f}', 'feh': '{:.3f}',
'feherr': '{:.3f}', 'd': '{:.5f}', 'derr':'{:.5f}', 'Av':'{:.5f}', 'Averr':'{:.5f}', 'ScaleF':'{:.5f}', 
'ScaleFerr':'{:.5f}', 'Fbol':'{:.5e}', 'Fbolerr':'{:.5e}', 'angRad':'{:.5f}', 'angRaderr':'{:.5f}', 
'lumi':'{:.5f}', 'lumierr':'{:.5f}', 'lumibol':'{:.5f}', 'lumibolerr':'{:.5f}', 'radius':'{:.5f}', 
'radiuserr':'{:.5f}', 'angRadbol':'{:.5f}', 'angRadbolerr':'{:.5f}', 'radiusbol':'{:.5f}', 'radiusbolerr':'{:.5f}', }
for col, f in formats.items(): target[col] = target[col].map(lambda x: f.format(x))
# output the final table
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.to_csv(combined_data_sedfit_path+"Output_SED_Fits_Interp.csv", index=False) 



################################################################################################
# collect what passbands used for the SED fitting.
print("Collect bands used for the SED fitting for indiviudal stars")
if Lbatch:
    files = glob.glob(batch_data_sedfit_path+"*.flux.csv")
    target = pd.DataFrame()
    for sub, file in enumerate(files):
        try: 
            file_tmp = pd.read_csv(file, dtype={"starID":str})
            target = pd.concat([target, file_tmp], axis=0, ignore_index=True)
        except:
            print("stars are skipped due to input stellar parameters are off the grid")  
else:
    files = glob.glob(single_data_sedfit_path+"*.flux.csv")
    columns = pd.read_csv(files[0]).columns
    target = combine_parallel(numb=len(files), columns=columns, func=collect_sed)


target["band"] = target["flux"].apply(lambda x: "|".join(x.split("|")[0::3]))   
target["flux_obs"] = target["flux"].apply(lambda x: "|".join(x.split("|")[1::3]))  
target["flux_pred"] = target["flux"].apply(lambda x: "|".join(x.split("|")[2::3]))  
target.drop(columns=["flux"], inplace=True)
target.drop_duplicates(subset=['starID'], inplace=True)
target = target.sort_values(by=["starID"]).reset_index(drop=True)
target.to_csv(combined_data_sedfit_path+"Output_SED_Fits_bandsUsed.csv", index=False)         




#################################################################################################
# Rejecting low-quality results.
Lcut_npoint = True
Lcut_noptic = False
npoint_min = 5
n_band_optic_min = 4
filename_interp = combined_data_sedfit_path+"/Output_SED_Fits_Interp.csv"
filename_no_interp = combined_data_sedfit_path+"/Output_SED_Fits_No_Interp.csv"
filename_bands = combined_data_sedfit_path+"/Output_SED_Fits_bandsUsed.csv"

# load fitting result.
target = pd.read_csv(filename_interp, dtype={"starID":str})
# cut with number of photometric measurements used for the Fitting.
if Lcut_npoint:
    SedFits = pd.read_csv(filename_no_interp, dtype={"starID":str}, usecols=["starID", "JumpBy", "npoint", "Redchi"])
    # SedFits = SedFits[["starID", "JumpBy", "npoint", "Redchi"]]
    SedFits = SedFits[SedFits.JumpBy=="DiffMet"].reset_index(drop=True)  # clip outliers
    SedFits = SedFits[SedFits.npoint>=npoint_min].reset_index(drop=True)  # clip outliers
    target = pd.merge(target, SedFits, on="starID").reset_index(drop=True)
# cut with number of photometric measuremtns in the optical
if Lcut_noptic:
    sed = pd.read_csv(filename_bands, dtype={"starID":str})
    sed.loc[:, "band"] = sed.band.str.replace(r"(\|twomassj)|(\|twomassh)|(\|twomassk)|(\|wisew1)|(\|wisew2)", "", regex=True)
    sed.loc[:, "n_band_optic"] = sed.band.str.split("|").str.len()
    sed = sed[sed.n_band_optic>=n_band_optic_min].reset_index(drop=True)
    target = pd.merge(target, sed, on="starID").reset_index(drop=True)
# ensure angular radius and bolometric flux are positive.
target = target[target.angRad>0].reset_index(drop=True)
target = target[target.Fbol>0].reset_index(drop=True)
target = target[np.isfinite(target.Fbol)].reset_index(drop=True)

# specify the format of each column.
formats = {'starID': '{:s}', 'modelteff': '{:.3f}', 'modeltefferr': '{:.3f}', 'modellogg': '{:.3f}', 'modelloggerr': '{:.3f}',
'teff': '{:.3f}', 'tefferr': '{:.3f}', 'logg': '{:.3f}', 'loggerr': '{:.3f}', 'feh': '{:.3f}', 'feherr': '{:.3f}', 
'd': '{:.5f}', 'derr':'{:.5f}', 'Av':'{:.5f}', 'Averr':'{:.5f}', 'ScaleF':'{:.5f}', 'ScaleFerr':'{:.5f}', 
'angRad':'{:.5f}', 'angRaderr':'{:.5f}', 'radius':'{:.5f}', 'radiuserr':'{:.5f}', 'lumi':'{:.5f}', 'lumierr':'{:.5f}', 
'Fbol':'{:.5e}', 'Fbolerr':'{:.5e}', 'lumibol':'{:.5f}', 'lumibolerr':'{:.5f}', 'angRadbol':'{:.5f}', 'angRadbolerr':'{:.5f}', 
'radiusbol':'{:.5f}', 'radiusbolerr':'{:.5f}', 'Redchi':'{:.2f}', 'npoint':'{:d}', 'JumpBy':'{:s}'}
for col, f in formats.items(): target[col] = target[col].map(lambda x: f.format(x))

# sort the columns and save the result.
columns = ['starID', 'modelteff', 'modeltefferr', 'modellogg', 'modelloggerr', 'teff', 'tefferr', 'logg', 'loggerr', 'feh', 'feherr', 
'd', 'derr', 'Av', 'Averr', 'ScaleF', 'ScaleFerr', 'angRad', 'angRaderr', 'radius', 'radiuserr', 'lumi', 'lumierr', 
'Fbol', 'Fbolerr', 'lumibol', 'lumibolerr', 'angRadbol', 'angRadbolerr', 'radiusbol', 'radiusbolerr', 'Redchi', 'npoint', 'JumpBy'] 
target[columns].to_csv(combined_data_sedfit_path+"Output_SED_Fits_Final.csv", index=False)