import glob
import pandas as pd
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from lmfit import Minimizer, Parameters, report_fit
import time
import re
import copy
import warnings
from collections import deque
from SedFitFunctions import *
import time
import sys
from scipy.interpolate import griddata
from inlist import *
import inlist
from phomKit import ExtLaws
from matplotlib.backends.backend_agg import FigureCanvas
import h5py
warnings.filterwarnings("ignore")

# run in the batch mode
if Lbatch==True:
    matplotlib.use('agg')
    index_start, index_end, job_index = int(sys.argv[1]), int(sys.argv[2])+1, int(sys.argv[3])
    if L_fig_h5py: hf = h5py.File(str(job_index), 'w')
else:
    matplotlib.use('Qt5Agg')

# ensure the entries at the same row of the three dataframes correspond to the same star.
params = pd.read_csv(file_input_params, dtype={"starID":str})
photo = pd.read_csv(file_input_photometry, dtype={"starID":str})
blackbody = pd.read_csv(Combined_BlackBody_Fits_path+"Output_BlackBody_Fits.csv",dtype={"starID":str})

# ensure to remove stars from the input list that have no valid blackbody fits.
params = pd.merge(blackbody[["starID"]], params, on="starID").reset_index(drop=True)
photo = pd.merge(blackbody[["starID"]], photo, on="starID").reset_index(drop=True)
if not ((params.starID==photo.starID).all() & (photo.starID==blackbody.starID).all()): 
    print("The tables of input parameters, input photometry, and blackbody fits are not consistent, thus exit.")
    exit()
# update index_start and index_end
if index_start>len(params)-1:
    exit()
else:
    if index_end > len(params): index_end = len(params)

# ensure in the input photometry table that, a magnitude error column of a given bandpass immediately follows its magnitude column.
subs = np.where(photo.columns!="starID")[0]
photo = photo[np.append("starID", photo.columns[subs].sort_values())].copy()

# Load pre-computed model flux densities, and stellar model parameters.
modelFluxDensity_ori = pd.read_csv(grid_list_path+"PassbandFlux_"+grid+".csv") 

# Stellar parameters.
# obsParams = pd.DataFrame(columns=["starID", "teff", "tefferr", "logg", "loggerr", "feh", "mass", "tv", "Av", "Averr", "Fbol", "Ref"])
obsParams = pd.DataFrame(columns=["starID", "teff", "tefferr", "logg", "loggerr", "feh", "mass", "tv", "Av", "Averr", "Fbol", "Ref"])
columns = ["starID", "teff", "tefferr", "logg", "loggerr", "Av", "Averr"]
if ~(["Av"] in params.columns.values): params.loc[:, "Av"] = np.nan
if ~(["Averr"] in params.columns.values): params.loc[:, "Averr"] = np.nan
obsParams[columns] = params[columns]

# metallicity
obsParams.loc[:, "feh"] = params[["feh"]] if LvaryFeh else np.zeros(len(obsParams)) 
obsParams.loc[:, "mass"] = np.ones(len(obsParams))      # assume 1.0 Msun, MARCS models most abundent
obsParams.loc[:, "tv"] = np.ones(len(obsParams))*2.0    # assume 2.0 km/s, MARCS models most abundent

# update A_lambda/A_V with a desired extinction law and R coefficient.
filters = pd.read_csv(filter_zero_point_path).set_index('band') 
filters.loc[:, "Alambda_Over_Av"] = ExtLaws(filters.loc[:, "lambda"].reset_index(), law=ExtinctionLaw, Rv=Rv)

# uses to save files for saving flux densities, result without interpolation and with interpolation.
FileOut, FileFinal, FileOptimal = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()


index_end=np.min([index_end, len(photo)])
#####################################################################################################
# # Part 1:  carry SED fitting using the original grid of model spectra, i.e., without interpolating 
for ii in range(index_start, index_end):  
    #############################################################################
    # (1) load input photometry
    t1=time.time() 
    starID, teff, logg = obsParams.loc[ii, ["starID", "teff", "logg"]]
    if obsParams.loc[ii, "tefferr"]/obsParams.loc[ii, "teff"]<0.01: obsParams.loc[ii, "tefferr"]=obsParams.loc[ii, "teff"]*0.01
    if obsParams.loc[ii, "tefferr"]/obsParams.loc[ii, "teff"]>0.04: obsParams.loc[ii, "tefferr"]=obsParams.loc[ii, "teff"]*0.04
    if obsParams.loc[ii, "loggerr"]<0.01: obsParams.loc[ii, "loggerr"]=0.01
    if obsParams.loc[ii, "loggerr"]>0.05: obsParams.loc[ii, "loggerr"]=0.05
    bands_dropped = ""

    # chech if the star is outsize of the grid.
    if grid=="BOSZ": 
        if teff<3500: 
            print("This star located outside the teff space is skipped for the SED fitting.")
            continue
    if grid=="MARCS":    
        if ((teff>8000) | (teff<2500)): 
            print("This star located outside the teff space is skipped for the SED fitting.")
            continue

    # Fetch observed SEDs. 
    photometry = pd.DataFrame(columns=["band", "mag", "magerr"])
    photometry["band"] = photo.iloc[ii, 1::2].index
    photometry["mag"] = photo.iloc[ii, 1::2].values
    photometry["magerr"] = photo.iloc[ii, 2::2].values
    photometry[["mag", "magerr"]] = photometry[["mag", "magerr"]].apply(pd.to_numeric, axis=1, errors='coerce')

    # read in the filter data, and sort bands by effective wavelength            
    photometry = photometry.set_index('band')
    photometry = pd.merge(photometry, filters, on="band").sort_values(by=["lambda"])  

    # original input photometry.
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input photometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print(photometry)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% input photometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

    # modify the offset of observed magntidues using prepared values from the table of the filter file.
    photometry.loc[:, "mag"] = photometry.loc[:, "mag"] + photometry.loc[:, "m0"]

    # don't use the WISE W3-4 bands. Blackbody.loc[ii, "Qflg"] sorted in wavelength,the same as photometry
    subs = np.where(np.isin(photometry.index.values, bandclip))[0]
    if len(subs)>0: 
        Qflg = list(blackbody.loc[ii, "Qflg"])
        for indx in subs: Qflg[indx] = "N"
        blackbody.loc[ii, "Qflg"] = "".join(Qflg)

    # save original photometry for a later comparison.
    # photometry_ori = photometry.copy()
    photometry_valid = photometry[~np.isnan(photometry.mag.values)].copy()
    photometry_ori = photometry_valid.copy()

    # Drop nan values and outliers using the blackbody fits.  
    if Lblackbody==True:
        # photometry = photometry_ori[np.array(list(blackbody.loc[ii, "Qflg"]))=="Y"]
        photometry = photometry[np.array(list(blackbody.loc[ii, "Qflg"]))=="Y"]
    else:
        photometry = photometry[~np.isnan(photometry.mag)]
    # up here, valid input photometry has been loaded.
    #############################################################################






    #############################################################################
    # (2) Select models whose stellar parametes are closest to the observed values. 
    modelParams = modelFluxDensity_ori[['teff', 'logg', 'feh', 'tv', 'mass', 'modelindex', 'fbol']]
    if np.isfinite(obsParams.loc[ii, "feh"]):
        feh_closest = modelParams.loc[np.argmin(abs(modelParams.feh.values-obsParams.loc[ii, "feh"])), "feh"]
    else:
        feh_closest = 0    
    modelParams = modelParams[modelParams.feh==feh_closest]     # Fix FeH=0.0

    # select models with teff close to the observed teff.
    if Lteffdiff==True: 
        # select a subsample of model spectra of which their Teffs close to the Teff of the target stars.
        if np.isfinite(obsParams.loc[ii, "teff"]):
            modelParams = modelParams[abs(modelParams.teff-obsParams.loc[ii, "teff"])<Teffdiff] # Model Teffs differ <1000K
        else:
            modelParams = modelParams[modelParams.teff<MaxTeff] # Model Teffs differ <6000K    
    else:         
        # select models with two closest logg values with respect to that of the target star.
        if np.isfinite(obsParams.loc[ii, "teff"]):
            subs = np.argsort(abs(np.unique(modelParams.teff)-obsParams.loc[ii, "teff"]))[:Nteff]
            teffgrid = np.unique(modelParams.teff)[subs]
        else:
            teffgrid = np.unique(modelParams[modelParams.teff<=MaxTeff])
    modelParams = modelParams[np.isin(modelParams.teff, teffgrid)].reset_index(drop=True)

    # select models with two closest logg values with respect to that of the target star.
    if np.isfinite(obsParams.loc[ii, "logg"]):
        subs = np.argsort(abs(np.unique(modelParams.logg)-obsParams.loc[ii, "logg"]))[:Nlogg]
        logggrid = np.unique(modelParams.logg)[subs]
    else:
        logggrid = np.unique(modelParams[modelParams.logg<=Maxlog].logg)
    # there is no models, given teff, select model with logg cloest to observed one.
    if len(modelParams[np.isin(modelParams.logg, logggrid)])==0: 
        bestlogg = np.unique(modelParams.logg)[np.argmin(abs(np.unique(modelParams.logg)-obsParams.loc[ii, "logg"]))]
        modelParams = modelParams[modelParams.logg==bestlogg].reset_index(drop=True).copy()
    else:
        modelParams = modelParams[np.isin(modelParams.logg, logggrid)].reset_index(drop=True).copy()

    # sort models with mass and tv
    modelParams = modelParams.sort_values(by=["mass", "tv"]).reset_index(drop=True).copy()
    if len(modelParams)>1: 
        # selec the model with one solar mass. This means spherical-geometry models are perfered over plane-parallel models
        if len(np.where(modelParams.mass==1.0)[0])>0:
            modelParams = modelParams[modelParams.mass==1.0].reset_index(drop=True).iloc[[0]]
        else:
            modelParams = modelParams.iloc[[0]]
    # up here, a few models with stellar parameters closest to observation are selected 
    #############################################################################


    #############################################################################
    # Jump for excluding observed photometric data that are too different from model predictions.
    while True: 
        minnpoint = copy.copy(minnpoint_ori)
        maxnpoint = copy.copy(maxnpoint_ori)

        # Convert magnitudes to flux densities. Below Lfluxnu has to be False, namely calculating Flux_lambda (flux
        # per unit wavelength), which are compared to model flux densities, also Flux_lambda. 
        point_lambda_obs = mag2fluxdensity(photometry, Av=0, Lfluxnu=False) 

        # DataFrame loop_model_rslt: for saving fitting results for each model at given npoint.
        columns = ["Av", "Averr", "ScaleF", "ScaleFerr", "Redchi", "Fbol", "ModelIndex", "RelDiff", "npoint"]   
        loop_model_rslt = pd.DataFrame(index=range(len(modelParams)), columns=columns)   
        if Lpeakpoints==True: maxnpoint = np.min([maxnpoint, len(point_lambda_obs)])

        # select the data points whose flux is not less tha FluxDropRatio times the max flux if the flag is set.
        if LFluxDropRatio==True:
            temp_n = len(point_lambda_obs[point_lambda_obs.flux>point_lambda_obs.flux.max()*FluxDropRatio])
            maxnpoint = np.min([maxnpoint, temp_n])

        # DataFrane final: for saving results for selected number of photometric point.    
        final =  pd.DataFrame(index=range(1), columns=columns)   
        npoint = copy.copy(maxnpoint)

        ###############################################################################################################
        # (3) given fixed number of bands of magnitudes, fit each model spectrum with the result is 
        # saved in "loop_model_rslt" and the best-fitting result is saved in final.
        for kk in range(len(modelParams)):
            # load model flux densities that were pre-calculated. Note that the flux densities are in the log scale. 
            modelFluxDensity = modelFluxDensity_ori[modelFluxDensity_ori.modelindex==modelParams.loc[kk, "modelindex"]].copy()
            columns = modelFluxDensity.columns.values[np.isin(modelFluxDensity.columns, photometry.index)]
            
            # note that band order in modelFluxDensity is probably different from the order in photometry DataFrame.
            fluxdensity = modelFluxDensity.iloc[0][columns]
            fluxdensity = pd.DataFrame({"band":fluxdensity.index, "flux":fluxdensity.values})
            fluxdensity = pd.merge(photometry[["lambda"]].reset_index(), fluxdensity, on="band").copy()
            point_lambda = fluxdensity.drop(columns=["band"]).rename(columns={"lambda":"wavelength"}).copy()


            # Make an initial guess of the scaling factor, defined as square angular diameter.   
            if Lblackbody==True: 
                scalingfactor = blackbody.loc[ii, "ScaleF"]*1.0 
            else: 
                if "gaiag" in point_lambda_obs.index: 
                    argminG = np.argmin(abs(point_lambda.wavelength-0.58363))
                    scalingfactor = point_lambda.loc[argminG, "flux"]/point_lambda_obs.loc["gaiag", "flux"]
                elif "gaiabp" in point_lambda_obs.index: 
                    argminBP = np.argmin(abs(point_lambda.wavelength-0.50209))
                    scalingfactor = point_lambda.loc[argminBP, "flux"]/point_lambda_obs.loc["gaiabp", "flux"]
                elif "gaiarp" in point_lambda_obs.index:              
                    argminRP = np.argmin(abs(point_lambda.wavelength-0.75888))
                    scalingfactor = point_lambda.loc[argminRP, "flux"]/point_lambda_obs.loc["gaiarp", "flux"]
                elif "twomassj" in point_lambda_obs.index: 
                    argminJ = np.argmin(abs(point_lambda.wavelength-1.235))
                    scalingfactor = point_lambda.loc[argminJ, "flux"]/point_lambda_obs.loc["twomassj", "flux"]
                elif "twomassh" in point_lambda_obs.index: 
                    argminH = np.argmin(abs(point_lambda.wavelength-1.662))
                    scalingfactor = point_lambda.loc[argminH, "flux"]/point_lambda_obs.loc["twomassh", "flux"]
                elif "twomassk" in point_lambda_obs.index: 
                    argminK = np.argmin(abs(point_lambda.wavelength-2.159))
                    scalingfactor = point_lambda.loc[argminK, "flux"]/point_lambda_obs.loc["twomassk", "flux"]                  
                else:
                    print("No GAIA and 2MASS photometry available. SED fitting ended here.")    
                    scalingfactor = 1e5

            # initial guess of the scaling factor is ready
            ScaleF = np.log10(scalingfactor)

            # select only top npoint photometric data points. 
            subs = np.argsort(np.array(point_lambda_obs["flux"], dtype=float))[::-1][:npoint] 

            # precalculated flux densities are in the log scale.
            XScaleF = np.array(point_lambda["flux"], dtype=float)[subs]
            YScaleF = np.log10(np.array(point_lambda_obs["flux"], dtype=float))[subs]
            YScaleFerr = np.array(point_lambda_obs.fluxerr/point_lambda_obs.flux/np.log(10), dtype=float)[subs]
            ExtinctCoeff =  np.array(point_lambda_obs["Alambda_Over_Av"][subs].copy())  

            # Start the fitting.
            if Lprint==True: 
                print("######################## loop over  models ##############################")
                print("Fitting using {:d} data points".format(npoint))
                print(photometry.iloc[subs,:])
                print(point_lambda_obs.iloc[subs,:])
                print("######################## loop over  models ##############################")
                print("\n\n\n")

            if LFixedAv==True:
                # Method 1: Av is known
                Av = obsParams.loc[ii, "Av"].copy()
                Averr = np.nan  
                # optimize the scaling factor
                ScaleF, ScaleFerr, Redchi = Fit_ScaleFactor(ScaleF, XScaleF-0.4*ExtinctCoeff*Av, YScaleF, err=YScaleFerr)
                # FluxPred = (XScaleF-0.4*ExtinctCoeff*Av) + ScaleF
                # RelDiff = np.sum(abs((FluxPred-YScaleF)/YScaleF))/len(YScaleF)  
            else:    
                # Method 2: Av is unknow, optimize extinction and the scaling factor simultaneously.
                Av, Averr, ScaleF, ScaleFerr, Redchi = Fit_Extinction_ScaleFactor(ExtinctCoeff, YScaleF-XScaleF, yerr=YScaleFerr)

            # calculate the flux density difference between observation and prediction.
            FluxPred = (XScaleF-0.4*ExtinctCoeff*Av) + ScaleF
            RelDiff = np.sum(abs((FluxPred-YScaleF)/YScaleF))/len(YScaleF)
            FbolatopEarth = np.nan

            # save 
            loop_model_rslt.loc[kk,:]=Av, Averr,ScaleF,ScaleFerr,Redchi,FbolatopEarth,modelParams.loc[kk,"modelindex"], RelDiff, npoint

        # save the result from the best fitting  spectrum to the DataFrame "final". 
        # Here, the number of bands used for fitting is fixed. 
        subs = np.argmin(np.array(loop_model_rslt.Redchi))
        final.iloc[0, :] = loop_model_rslt.loc[subs, :].copy()
        # up here, given the current number of passbands, the best fitting  model spectrum has been found out of all
        # the preselected models spectra,    
        ######################################################################################################


        ######################################################################################################
        # (4) Remove photometric outliers by calculating the difference between observed and best-fitting model fluxes.
        # this part does not do any fitting. Note that the best-fitting model spectrum has been found just above. 
        Av, Averr, ScaleF, ScaleFerr, Redchi, FbolatopEarth, ModelIndex, RelDiff, npoint = final.iloc[0, :]
        point_lambda_obs = mag2fluxdensity(photometry, Av=Av, Lfluxnu=False)

        # load model flux densities that were pre-calculated. Note that the flux densities are in the log scale. 
        # modelFluxDensity = modelFluxDensity_ori[modelFluxDensity_ori.modelindex==modelParams.loc[kk, "modelindex"]].copy()
        modelFluxDensity = modelFluxDensity_ori[modelFluxDensity_ori.modelindex==ModelIndex].copy()
        filename, modelteff, modellogg = modelFluxDensity.reset_index(drop=True).loc[0, ["filename", "teff", "logg"]]
        columns = modelFluxDensity.columns.values[np.isin(modelFluxDensity.columns, photometry.index)]

        # note that the band order in modelFluxDensity is probably different from the order in photometry DataFrame.
        fluxdensity = modelFluxDensity.iloc[0][columns]
        fluxdensity = pd.DataFrame({"band":fluxdensity.index, "flux":fluxdensity.values})
        fluxdensity = pd.merge(photometry[["lambda"]].reset_index(), fluxdensity, on="band").copy()
        point_lambda = fluxdensity.drop(columns=["band"]).rename(columns={"lambda":"wavelength"}).copy()

        # here model flux densities are in log scale, so is converted to the same scale as the observed conoutpart
        point_lambda.loc[:, "flux"]  = 10**point_lambda.loc[:, "flux"].copy() 

        # Scale the best fitting grid model and convolved/predicted fluxes.  
        scalingfactor = 10**ScaleF 
        point_lambda.loc[:, "flux"] *= scalingfactor
        FbolatopEarth=np.nan

        # for debuging
        if Lprint==True: 
            print("################################### loop over photometric bands #####################################")
            subs = np.argsort(np.array(point_lambda_obs["flux"], dtype=float))[::-1][:npoint] 
            print("Fitting using {:d} data points".format(npoint))
            print(photometry.iloc[subs,:])
            print(point_lambda_obs.iloc[subs,:])
            print("\n")
            print("Best Model: {:}".format(filename))
            print("log(ScalingFactor)=log(square angular diameter)={:.4f}+/-{:.4f}, Av={:.4f}".format(ScaleF, ScaleFerr, Av))
            print("Best fitting result with this number of bands")
            print(final)
            print("################################### loop over photometric bands #####################################")
            print("\n\n\n")
            

            
        # Drop data points which are far from predictions, here ignoring W3-4 bands.
        FluxWavelength = np.array(point_lambda_obs.loc[:, "wavelength"].copy())
        FluxDensityObs = np.array(point_lambda_obs.loc[:, "flux"].copy(), dtype=float)
        FluxDensitypred = np.array(point_lambda.loc[:, "flux"].copy(), dtype=float)
        MaxDiff = np.max(abs(FluxDensityObs-FluxDensitypred)/FluxDensitypred)
        # use the max flux densities as reference
        # MaxDiff = np.max(abs(FluxDensityObs-FluxDensitypred)/np.max(FluxDensityObs))

        if MaxDiff<FluxMaxDiff: 
            # print("=============================== result summary ================================")
            print("=============================== result summary ================================")
            print("#{:d} SED fitting for the target {:s}".format(ii, starID))
            print("Model spectrum used for the fitting (uninterpolated):")
            print(modelParams, "\n")
            print("Photometric data points dropped:")
            print(bands_dropped, "\n")
            print("Jumped out of the loop to optimise the number of bands for fitting: ")
            print("diff. between obs. and pred. fluxes is less than the requirement", "\n")
            subs = np.argsort(np.array(point_lambda_obs["flux"], dtype=float))[::-1][:npoint] 
            print("Final photometry input for the fitting using {:d} data points".format(npoint))
            Jump="DiffMet"
            # print(photometry.iloc[subs,:])
            print(point_lambda_obs.iloc[subs,:], "\n")
            break
        else:
            if len(FluxDensityObs)<=minnpoint: 
                print("grid model atmospheres used for the fitting:")
                print(modelParams, "\n")
                print("Jumped out of the loop to optimise the number of bands for fitting: ")
                print("condition of the min number of bands satisfied.", "\n")
                subs = np.argsort(np.array(point_lambda_obs["flux"], dtype=float))[::-1][:npoint] 
                print("Final photometry input for the fitting using {:d} data points".format(npoint))
                Jump="NumbMet"
                # print(photometry.iloc[subs,:])
                print(point_lambda_obs.iloc[subs,:], "\n")
                break  
            SubsDiff = np.argmax(abs(FluxDensityObs-FluxDensitypred)/FluxDensitypred)
            # SubsDiff = np.argmax(abs(FluxDensityObs-FluxDensitypred)/np.max(FluxDensitypred))
            # SubsDiff = np.argmax(abs(FluxDensityObs-FluxDensitypred)/np.max(FluxDensityObs))
            if Lprint==True: 
                print("=============================== outlier photometry rejection ================================")
                print(" {:} band dropped as flux density diff. is {:.3f}".format(photometry.index[SubsDiff], MaxDiff))
                print("=============================== outlier photometry rejection ================================")
                print("\n\n\n")
            bands_dropped += photometry.index[SubsDiff] + "|"
            #Note that although the length of photometry and FluxDensityObs is different, but the heading elements in both are the same
            photometry = photometry.drop([photometry.index[SubsDiff]]).copy()

    ###################################################################################################################################
    # up here, photometry outliers have been clipped.
    ###################################################################################################################################


    ###################################################################################################################################
    # (5) For the MARCS grid, extrapolate model spectra with wavelength up to 30 micron, 10 micron longer exteneded beyong the red end.
    # Since the maximum wavelength of marcs spectra is 20 micron and the W4 bandpass is nearly entirely beyond this wavelength, MARCS
    # model spectra need to be extrapolated. Extrapolate the spectrum with wavelength longer than 20 micron by carrying out a 3rd-order 
    # polynomial fitting to log10(flux) as a function of wavelength.
    spec_ori = readspec(grid_spectrum_path, filename, grid=grid)
    if grid=="MARCS": spec_ori = flux_extrap(spec_ori.copy(), Minwavelength=8, Maxwavelength=30, order=3)

    ###################################################################################################################################
    # (6) save the best SED-fitting parameters, and predicted and extinction-corrected observed flux densities.
    # Calculate the model fluxes.
    spec_lambda = resolution(spec_ori, res=100)
    fluxdensity= FluxDensity(spec_ori, bandpass=np.array(photometry.index), path=filter_transmission_path)
    wavelength =  np.array(photometry["lambda"])
    point_lambda = pd.DataFrame({"wavelength":wavelength, "flux":fluxdensity})

    # scale the best-fitting model spectrum    
    scalingfactor = 10**ScaleF 
    spec_lambda.loc[:, "flux"] *= scalingfactor
    spec_lambda.loc[:, "flux_sm"] *= scalingfactor
    point_lambda.loc[:, "flux"] *= scalingfactor

    # calculate the bolometric flux
    sed = spec_lambda.copy()
    FbolatopEarth=np.trapz(sed.flux, x=sed.wavelength*1e4)
    final.loc[0, "Fbol"] = FbolatopEarth

    # save predicted and deextincted, observed fluxes.
    FluxBandpass = np.array(point_lambda_obs.index.copy())
    FluxDensityObs = np.array(point_lambda_obs.loc[:, "flux"].copy(), dtype=float)
    FluxDensitypred = np.array(point_lambda.loc[:, "flux"].copy(), dtype=float)
    out = pd.DataFrame({"bandpass":FluxBandpass, "FluxObs":FluxDensityObs, "Fluxpred":FluxDensitypred})

    # define output format.
    out.loc[:, "FluxObs"] = np.log10(out.loc[:, "FluxObs"]).round(6)
    out.loc[:, "Fluxpred"] = np.log10(out.loc[:, "Fluxpred"]).round(6)
    out_tmp = pd.DataFrame()
    out_tmp.loc[0, "flux"] = '|'.join([''.join(i) for i in out.values.astype(str).flatten()])
    out_tmp.loc[0, "starID"] = str(starID)

    if Lbatch:
        FileOut = pd.concat([FileOut, out_tmp], axis=0, ignore_index=True)
    else:
        if Lsave: out_tmp.to_csv(single_data_sedfit_path+str(starID)+".flux.csv",  float_format="%.5e", index=False)

    # save best-fitting parameters.
    final = final.dropna(subset=["ScaleF"])
    final.loc[0,"JumpBy"]=Jump
    final.loc[0, "starID"] = str(starID)
    if Lbatch:
        # final.loc[0, "starID"] = str(starID)
        FileFinal = pd.concat([FileFinal, final], axis=0, ignore_index=True)
    else:    
        if Lsave: final.to_csv(single_data_sedfit_path+str(starID)+".BestfitModels.csv",  float_format="%.5e", index=False)

    # Show all the original photometric data, namely without discarding bad observed data, later used for plotting.
    point_lambda_obs_ori = mag2fluxdensity(photometry_ori, Av=Av, Lfluxnu=False)
    fluxdensity_ori= FluxDensity(spec_ori, bandpass=np.array(photometry_ori.index), path=filter_transmission_path)
    wavelength_ori =  np.array(photometry_ori["lambda"])
    point_lambda_ori = pd.DataFrame({"wavelength":wavelength_ori, "flux":fluxdensity_ori})   
    point_lambda_ori.loc[:, "flux"] *= scalingfactor
    ##################################################################################################################################
    # up here, all the SED fitting for given  models is done.
    ##################################################################################################################################




    ###################################################################################################################################
    # Part 2: Interpolate model flux densities.
    if Linterpolate==True:
        while True:
            # (1) select 4 closest model logg values and 4 model teff values, leading to 16 points in the teff-logg space. 
            # The model flux densities are in the log scale.
            sedmodels = modelFluxDensity_ori.copy()
            sedmodels = sedmodels[sedmodels.feh==feh_closest].copy()

            # select models with four closest mldel Teff values with respect to that of the target star.
            if np.isfinite(obsParams.loc[ii, "teff"]):
                subs = np.argsort(abs(np.unique(sedmodels.teff)-obsParams.loc[ii, "teff"]))[:4]
                teffgrid = np.unique(sedmodels.teff)[subs]

            # select models with four closest model logg values with respect to that of the target star.
            logg_ori=np.nan
            if not np.isfinite(float(obsParams.loc[ii, "logg"])):
                logg_ori = obsParams.loc[ii, "logg"].copy()
                obsParams.loc[ii, "logg"] = float(modellogg)
            subs = np.argsort(abs(np.unique(sedmodels.logg)-obsParams.loc[ii, "logg"]))[:4]
            logggrid = np.unique(sedmodels.logg)[subs]

            # check if the star is located outside the the  parameter space, the interpolation part was skipped.
            LoffGrid=False
            if (obsParams.loc[ii, "logg"]<np.min(logggrid)) | (obsParams.loc[ii, "logg"]>np.max(logggrid)): 
                print("the star is located outside the  logg space, the interpolation part was skipped.")
                LoffGrid=True
                break
            if (obsParams.loc[ii, "teff"]<np.min(teffgrid)) | (obsParams.loc[ii, "teff"]>np.max(teffgrid)): 
                print("the star is located outside the teff space, the interpolation part was skipped.")
                LoffGrid=True
                break
            sedmodels = sedmodels[np.isin(sedmodels.teff, teffgrid)].reset_index(drop=True).copy()
            sedmodels = sedmodels[np.isin(sedmodels.logg, logggrid)].reset_index(drop=True).copy()
            
            if (sedmodels.teff.unique().shape[0]<2) | (sedmodels.logg.unique().shape[0]<2): 
                print("There are less than four models with diff. teff & logg combinations.")
                LoffGrid=True
                break

            # if there are too many model spectra, pick relevant.
            sedmodels = sedmodels.sort_values(by=['teff', 'logg', 'feh', 'mass', 'tv']).reset_index(drop=True).copy()
            if grid=="MARCS": sedmodels = sedmodels.drop_duplicates(subset=['teff', 'logg', 'feh'], keep='first')
            if grid=="BOSZ": sedmodels = sedmodels[sedmodels.alpha==0.0].reset_index(drop=True).copy()
            ##################################################################################################################################
            # Up here, up to 16 model bracketing the input teff and logg values are selected.
            ##################################################################################################################################

            # (2) interpolate model fluxes to a grid of Teff and logg.
            # to perform irregular grid interpolation, follow this example:
            # https://scipython.com/book/chapter-8-scipy/examples/two-dimensional-interpolation-with-scipyinterpolategriddata/
            teffs = np.arange(min(sedmodels.teff), max(sedmodels.teff)+teffStep, teffStep)
            loggs = np.arange(min(sedmodels.logg), max(sedmodels.logg)+loggStep, loggStep)
            teff_input = sedmodels.teff.values.copy()
            logg_input = sedmodels.logg.values.copy()

            # Now ensure that sedmodels have the same columns order and content as photometry dataframe, except for the last column for fbol
            sedmodels = sedmodels[np.append(photometry.index.values, "fbol")].reset_index(drop=True).copy()
            flux_input = sedmodels.values.copy()

            # note that the shape of TEFF np.shape(TEFFs) is not (len(teffs), len(loggs)), it is (len(loggs), (len(teffs)). This 
            # is accutally consistent, since the first dimension of array in Python is row (y axis), followed by column (x axis).
            TEFFs, LOGGs = np.meshgrid(teffs, loggs)
            
            # note the shape, which is correct, being consistent with those of TEFFs and LOGGs
            flux_interp = np.zeros((len(loggs), len(teffs), len(photometry)+1)) # add column for fbol 


            # interpolate model flux densities and bolometric flux. 
            for indm in range(len(photometry)+1):
                # note that the shape of griddata output is the same as that of TEFFs, or LOGGs.
                # for flux_interp, the first demension (i.e., left) corresponds to logg, the second to teff, and the third to passband.
                flux_interp[:, :, indm] = griddata((teff_input, logg_input), flux_input[:, indm], (TEFFs, LOGGs), method="linear")

            ###################################################################################################################################
            # up here the flux interpolation is done.    
            #############################################################################################################################


            ###################################################################################################################################
            # (3) optimize Av and scaling factor for each node of the grid.
            # grid Teff values stored in teffs, grid logg values stored in loggs, and interpolated flux stored flux_interp. Note that for 
            # flux_interp (in the log scale), the first demension (i.e., left) corresponds to logg, the second to teff, and the third to passband.
            # interp_FitResult: save SED fitting results for each grid node.
            columns = ["logg", "teff", "Av", "Averr", "ScaleF", "ScaleFerr", "Redchi", "fbol_emergent"]
            interp_FitResult = pd.DataFrame(index=range(np.product(np.shape(TEFFs))), columns=columns)    

            # prepare data for the SED fitting.
            interp_point_lambda_obs = mag2fluxdensity(photometry, Av=0, Lfluxnu=False)
            subs_nopint = np.argsort(np.array(interp_point_lambda_obs["flux"], dtype=float))[::-1][:npoint] 
            YScaleF = np.log10(np.array(interp_point_lambda_obs["flux"], dtype=float))[subs_nopint]
            YScaleFerr = np.array(interp_point_lambda_obs.fluxerr/interp_point_lambda_obs.flux/np.log(10), dtype=float)[subs_nopint]
            ExtinctCoeff =  np.array(interp_point_lambda_obs["Alambda_Over_Av"][subs_nopint].copy())

            # start the SED-fitting for each node.
            interp_subs = 0 
            modelMiss=False
            for indx in range(len(loggs)):
                for indy in range(len(teffs)):        
                    # note that model_flux_grid (band flux) has been sorted in the order of wavelength, in line with the order of the variable "photometry".
                    # Note that the last column of flux_interp for a given logg and teff is fbol
                    XScaleF = flux_interp[indx, indy, :-1][subs_nopint]

                    # check if an interpolated model is available or not. The unavailability can exist if original models are not available, for instance, near 
                    # the edge of the 4 x 4 grid (teff vs logg).
                    if not (np.isnan(XScaleF).any()):   
                        # Start the fitting.
                        if Lprint==True: 
                            # print("######################## loop over interpolated models ##############################")
                            print("#{:d}th of logg list, #{:d}th of Teff list".format(indx, indy))
                            # print("######################## loop over MARCS models ##############################")
                        if LFixedAv==True:
                            # Method 1: Av is known
                            interp_Av = obsParams.loc[ii, "Av"].copy()
                            # optimize the scaling factor
                            interp_ScaleF, interp_ScaleFerr, interp_Redchi = Fit_ScaleFactor(ScaleF, XScaleF-0.4*ExtinctCoeff*Av, YScaleF, err=YScaleFerr)      
                            interp_FitResult.loc[interp_subs, :] = loggs[indx], teffs[indy], interp_Av, 0.0, interp_ScaleF, interp_ScaleFerr, interp_Redchi, flux_interp[indx, indy, -1]        
                        else:    
                            # Method 2: Av is unknow, optimise Av and the scaling factor simultaneously.
                            Av, Averr, ScaleF, ScaleFerr, Redchi = Fit_Extinction_ScaleFactor(ExtinctCoeff, YScaleF-XScaleF, yerr=YScaleFerr)
                            interp_FitResult.loc[interp_subs, :] = loggs[indx], teffs[indy], Av, Averr, ScaleF, ScaleFerr, Redchi, flux_interp[indx, indy, -1]
                    else:
                        interp_FitResult.loc[interp_subs, :] = np.append([loggs[indx], teffs[indy]], np.full([6], np.nan))
                        if not modelMiss: print("models missing")
                        modelMiss=True
                    interp_subs += 1
                    

            # compute bolometric fluxed received on Earth from that on the stellar surface.       
            interp_FitResult.loc[:, "fbol"] = interp_FitResult.fbol_emergent + interp_FitResult.ScaleF    
            interp_FitResult.loc[:, "fbolerr"] = interp_FitResult.ScaleFerr.copy()
            interp_FitResult.drop(columns=["fbol_emergent"], inplace=True)
            ###################################################################################################################################
            # up here the SED-fitting has done for all the interpolated models, and the fitting result for each node is saved in interp_FitResult   
            #############################################################################################################################          

            # (4) find the best fitting model out of all the interpolated and original models using a Bayesian method. 
            # Using Gaussian priors for teff.
            if ((np.isnan(obsParams.loc[ii, "tefferr"])) | (obsParams.loc[ii, "tefferr"]<=0)): 
                priors = gaussian_func(TEFFs, obsParams.loc[ii, "teff"], 400) 
            else: 
                priors = gaussian_func(TEFFs, obsParams.loc[ii, "teff"], np.sqrt(obsParams.loc[ii, "tefferr"]**2 +  tefferr_sys**2)) 
            priors /= np.max(priors)    
            # Using Gaussian priors for teff.
            if ((np.isnan(obsParams.loc[ii, "loggerr"])) | (obsParams.loc[ii, "loggerr"]<=0)):
                priors *= gaussian_func(LOGGs, obsParams.loc[ii, "logg"], 0.2)
            else:
                priors *= gaussian_func(LOGGs, obsParams.loc[ii, "logg"], np.sqrt(obsParams.loc[ii, "loggerr"]**2 +  loggerr_sys**2))  
            priors /= np.max(priors)



            # The likelihood function is computed by taking the exponential of the reduced Chi sqare, and is multiplied by the prior to estimate a posterior.
            posteriors = priors*np.exp(-interp_FitResult.Redchi.astype(float).values).reshape(len(loggs), len(teffs))
            if (np.nanmax(posteriors)==0) | (np.isnan(np.nanmax(posteriors))): 
                LoffGrid=True
                print("ReChi2 of interpolated models too large, thus this star is skipped")
                print("=============================== result summary ================================")
                break

            # Using Gaussian priors for Av.
            if LAvprior:
                Av_mesh = interp_FitResult.Av.astype(float).values.reshape(len(loggs), len(teffs))
                if ((np.isnan(obsParams.loc[ii, "Averr"])) | (obsParams.loc[ii, "Averr"]<=0)):
                    posteriors *= gaussian_func(Av_mesh, obsParams.loc[ii, "Av"], 0.05)
                else:
                    posteriors *= gaussian_func(Av_mesh, obsParams.loc[ii, "Av"], np.sqrt(obsParams.loc[ii, "Averr"]**2 +  averr_sys**2))  
            posteriors /= np.nanmax(posteriors)

            

            # now calculate the estimates of fitting parameters and their errors by sampling the posterior.
            interp_FitResult.loc[:, "posterior"] = posteriors.flatten()/np.nansum( posteriors.flatten())
            samples = np.random.choice(interp_FitResult.index, size=int(2e5), p=interp_FitResult.posterior.fillna(0))
            interp_FitResult_samples = interp_FitResult.loc[np.array(samples, dtype=int), :]
            interp_FitResult[interp_FitResult.columns] = interp_FitResult[interp_FitResult.columns].apply(pd.to_numeric, errors='coerce')

            # drop nan values.
            interp_FitResult_dropna = interp_FitResult.dropna()

            # use this star to understand why the posterior is a mixture of multiple Gaussians that share same range of the independent variable (i.e., x).
            # The reason is that: although the joint prior of Teff and logg is Gaussian as a continous function, but can be asymmetric on a grid of Teff
            # and logg. 
            LcheckPost=False
            if ((LcheckPost) & (obsParams.loc[ii, "starID"]=='1432587')):
                interp_FitResult_sub = interp_FitResult[(interp_FitResult.teff<4270)&(interp_FitResult.teff>4170)&
                (interp_FitResult.posterior>0.0045)&(interp_FitResult.posterior<0.01)]
                interp_FitResult_sub2 = interp_FitResult[(interp_FitResult.teff<4270)&(interp_FitResult.teff>4170)&
                (interp_FitResult.posterior>0.000)&(interp_FitResult.posterior<0.003)]
                plt.contourf(TEFFs, LOGGs, priors, levels=14, cmap="RdBu_r")
                plt.scatter(interp_FitResult.loc[:, "teff"], interp_FitResult.posterior, c="g")
                plt.scatter(interp_FitResult_sub2.loc[:, "teff"], interp_FitResult_sub2.posterior, c="k")
                plt.scatter(interp_FitResult_sub2.loc[:, "teff"], interp_FitResult_sub2.loc[:, "logg"], c="k")
                plt.scatter(interp_FitResult_sub.loc[:, "teff"], interp_FitResult_sub.posterior, c="r")
                plt.scatter(interp_FitResult_sub.loc[:, "teff"], interp_FitResult_sub.loc[:, "logg"], c="r")
                plt.scatter(obsParams.loc[ii, "teff"], obsParams.loc[ii, "logg"])


            # (a) Search for the interpolated model whose posterior is the highest (mode estimator). 
            subs = np.where(posteriors==np.nanmax(posteriors))
            optm_teff, optm_logg = TEFFs[subs][0], LOGGs[subs][0]
            optimal_mod = interp_FitResult[(interp_FitResult.teff==optm_teff)&(interp_FitResult.logg==optm_logg)].reset_index(drop=True).copy()
            optimal_mod.drop(columns=["posterior"], inplace=True)
            optimal_mod.columns += "_mod"
            print("Mode estimates (interpolated):")
            print(optimal_mod, "\n")

            # (b) compute weighted mean and std using posteriors.
            optimal = pd.DataFrame(index=range(1))
            optimal.loc[:, "logg"], optimal.loc[:, "loggerr"] = weighted_mean_std(interp_FitResult_dropna.logg, prob=interp_FitResult_dropna.posterior)
            optimal.loc[:, "teff"], optimal.loc[:, "tefferr"] = weighted_mean_std(interp_FitResult_dropna.teff, prob=interp_FitResult_dropna.posterior)
            optimal.loc[:, "Av"], optimal.loc[:, "Averr"] = weighted_mean_std(interp_FitResult_dropna.Av, prob=interp_FitResult_dropna.posterior)
            optimal.loc[:, "ScaleF"], optimal.loc[:, "ScaleFerr"] = weighted_mean_std(interp_FitResult_dropna.ScaleF, prob=interp_FitResult_dropna.posterior)
            optimal.loc[:, "fbol"], optimal.loc[:, "fbolerr"] = weighted_mean_std(interp_FitResult_dropna.fbol, prob=interp_FitResult_dropna.posterior)
            if LweightedMean: optimal = pd.concat([optimal, optimal_mod], axis=1)


            # Fit gassians to marginalized posteriors.
            # use weighted average to approximate logg, while use Gaussian fits to approximate the the rest.
            optimal.loc[:, "logg"], optimal.loc[:, "loggerr"] = weighted_mean_std(interp_FitResult_dropna.logg, prob=interp_FitResult_dropna.posterior)  

            # for the ["teff", "Av", "ScaleF", "fbol"] fit a Gaussian to each.
            interp_FitResult_samples_fit = interp_FitResult_samples[["teff", "Av", "ScaleF", "fbol"]].copy()
            fit_hist = pd.DataFrame(columns=["teff_bin", "teff_prob", "Av_bin", "Av_prob", "ScaleF_bin", "ScaleF_prob", "fbol_bin", "fbol_prob"])
            fit_hist_params = pd.DataFrame(columns=["H", "Herr", "P", "Perr", "W", "Werr", "N", "Nerr"], index=["teff", "Av", "ScaleF", "fbol"])
            for col in interp_FitResult_samples_fit.columns:
                # show the histograms. Note that spikes generally occur in high resolution histograms (large samples), due to the quantisation noise. 
                # Recall that the quantisation noise is the effect of representing an analog continuous signal with a discrete number. To suppress the 
                # spikes, a Gaussian random noise (mu=0, sigma=5/6*bin size) is added to each the count number in each bin. The sigma is choosen smaller 
                # enough to not significantly change the statistics of the distribution.
                mu = 0
                sigma = (interp_FitResult_samples.loc[:, col].max()-interp_FitResult_samples.loc[:, col].min())/120 # mean and standard deviation
                if sigma==0: 
                    LoffGrid=True
                    print("ReChi2 of interpolated models too large, thus this star is skipped")
                    print("=============================== result summary ================================")
                    break
                interp_FitResult_samples.loc[:, col] += np.random.normal(mu, sigma, len(interp_FitResult_samples_fit))
                
                # calculate histograms.
                bins=np.linspace(interp_FitResult_samples.loc[:, col].min(), interp_FitResult_samples.loc[:, col].max(), 100)
                density, bin_edges = np.histogram(interp_FitResult_samples.loc[:, col], bins=bins, density=True)
                bin_center = bin_edges[:-1] + (bin_edges[1]-bin_edges[0])/2
                
                # save histograms.
                fit_hist.loc[:, col+"_bin"] = bin_center
                fit_hist.loc[:, col+"_prob"] = density

                # fit a Gaussian to the histogram, namely a marginalized posterior.
                param_guess = [np.nanmax(density), optimal_mod.loc[0, col+"_mod"], optimal.loc[0, col+"err"], 0] 
                H, Herr, P, Perr, W, Werr, N, Nerr = Fit_Gaussian(bin_center, density, param_guess)
                optimal.loc[:, col], optimal.loc[:, col+"err"] = P, W

                #store fitting results.
                fit_hist_params.loc[col, :] =  H, Herr, P, Perr, W, Werr, N, Nerr 
            if LGassianFit: optimal = pd.concat([optimal, optimal_mod], axis=1)

            # for output
            optimal.loc[0, "starID"] = str(starID)
            if Lbatch:
                optimal.loc[0, "starID"] = str(starID)
                FileOptimal = pd.concat([FileOptimal, optimal], axis=0, ignore_index=True)
            else:
                if Lsave: optimal.to_csv(single_data_sedfit_path+str(starID)+".interpBestModel.csv",  float_format="%.5e", index=False)

            # Since this star has no logg in the original input, it's set to modellogg of the best-fitting MARCS model. Now set it back to NaN.   
            if np.isfinite(logg_ori): obsParams.loc[ii, "logg"] = np.nan
        
            # now the interpolation process is done.    
            break
            ###################################################################################################################################
            # up here the SED fitting over interpolated fluxes is done.    
            #############################################################################################################################



    ######################################################################################################################################################
    # Save figures.
    # Panel 1: show interpolation grid.
    if LoffGrid==False:
        bandpass_show="gaiag"
        if len(np.where(photometry.index.values==bandpass_show)[0])==0: bandpass_show=photometry.index.values[0]
        subs = np.where(photometry.index.values==bandpass_show)[0][0]
        
        fig, ax = plt.subplots(3,3, figsize=(15,12))
        ax = ax.flat
        cntr1 = ax[0].contourf(TEFFs, LOGGs, flux_interp[:,:, subs], levels=14, cmap="RdBu_r")
        cbar = fig.colorbar(cntr1, ax=ax[0])
        cbar.ax.set_ylabel("log(Flux)")
        ax[0].scatter(TEFFs, LOGGs,  edgecolor="k", facecolor="none", lw=0.1)
        ax[0].scatter(teff_input, logg_input, s=100, marker="D", edgecolor="k", facecolor="none",lw=3)
        ax[0].scatter(obsParams.loc[ii, "teff"], obsParams.loc[ii, "logg"], s=100, marker="s", c="r")
        ax[0].scatter(modelParams.loc[:, "teff"], modelParams.loc[:, "logg"], s=100, marker="s", c="b")
        ax[0].set_xlim(np.flip(ax[0].get_xlim()))
        ax[0].set_ylim(np.flip(ax[0].get_ylim()))
        ax[0].set_xlabel("Effective temperature (K)", fontsize=14)
        ax[0].set_ylabel("Surface gravity", fontsize=14)

        cntr1 = ax[1].contourf(TEFFs, LOGGs, interp_FitResult.Redchi.values.reshape(len(loggs), len(teffs)), levels=14, cmap="RdBu_r")
        cbar = fig.colorbar(cntr1, ax=ax[1])
        cbar.ax.set_ylabel("Chi2")
        subs = interp_FitResult.Redchi.astype(float).idxmin()
        ax[1].scatter(interp_FitResult.loc[subs,"teff"], interp_FitResult.loc[subs,"logg"], s=150, marker="*", color="r")  
        ax[1].scatter(teff_input, logg_input, s=100, marker="D", edgecolor="k", facecolor="none",lw=3)
        ax[1].scatter(obsParams.loc[ii, "teff"], obsParams.loc[ii, "logg"], s=100, marker="s", c="r")
        ax[1].scatter(modelParams.loc[:, "teff"], modelParams.loc[:, "logg"], s=100, marker="s", c="b")
        ax[1].set_xlim(np.flip(ax[1].get_xlim()))
        ax[1].set_ylim(np.flip(ax[1].get_ylim()))
        ax[1].set_xlabel("Effective temperature (K)", fontsize=14)
        ax[1].set_ylabel("Surface gravity", fontsize=14)

        cntr1 = ax[2].contourf(TEFFs, LOGGs, interp_FitResult.Av.values.reshape(len(loggs), len(teffs)), levels=14, cmap="RdBu_r")
        cbar = fig.colorbar(cntr1, ax=ax[2])
        cbar.ax.set_ylabel("Extinction (mag)")
        ax[2].scatter(teff_input, logg_input, s=100, marker="D", edgecolor="k", facecolor="none",lw=3)
        ax[2].scatter(obsParams.loc[ii, "teff"], obsParams.loc[ii, "logg"], s=100, marker="s", c="r")
        ax[2].scatter(modelParams.loc[:, "teff"], modelParams.loc[:, "logg"], s=100, marker="s", c="b")
        subs = interp_FitResult.Redchi.astype(float).idxmin()
        ax[2].scatter(interp_FitResult.loc[subs,"teff"], interp_FitResult.loc[subs,"logg"], s=150, marker="*", color="r")  
        subs = np.where(posteriors==np.nanmax(posteriors))
        ax[2].scatter(TEFFs[subs], LOGGs[subs], s=150, marker="*", color="lime")  
        ax[2].set_xlim(np.flip(ax[2].get_xlim()))
        ax[2].set_ylim(np.flip(ax[2].get_ylim()))
        ax[2].set_xlabel("Effective temperature (K)", fontsize=14)
        ax[2].set_ylabel("Surface gravity", fontsize=14)

        cntr1 = ax[3].contourf(TEFFs, LOGGs, posteriors, levels=14, cmap="RdBu_r")
        cbar = fig.colorbar(cntr1, ax=ax[3])
        cbar.ax.set_ylabel("Posterior")
        subs = interp_FitResult.Redchi.astype(float).idxmin()
        ax[3].scatter(interp_FitResult.loc[subs,"teff"], interp_FitResult.loc[subs,"logg"], s=150, marker="*", color="r")  
        ax[3].scatter(teff_input, logg_input, s=100, marker="D", edgecolor="k", facecolor="none",lw=3)
        ax[3].scatter(obsParams.loc[ii, "teff"], obsParams.loc[ii, "logg"], s=100, marker="s", c="r")
        ax[3].scatter(modelParams.loc[:, "teff"], modelParams.loc[:, "logg"], s=100, marker="s", c="b")
        subs = np.where(posteriors==np.nanmax(posteriors))
        ax[3].scatter(TEFFs[subs], LOGGs[subs], s=150, marker="*", color="lime")  
        ax[3].set_xlim(np.flip(ax[3].get_xlim()))
        ax[3].set_ylim(np.flip(ax[3].get_ylim()))
        ax[3].set_xlabel("Effective temperature (K)", fontsize=14)
        ax[3].set_ylabel("Surface gravity", fontsize=14)

        ax[4].grid(b=True, which='minor', zorder=1)
        ax[4].grid(b=True, which='major', zorder=2)
        ax[4].plot(spec_lambda.wavelength, spec_lambda.flux_sm, c="k", label="smoothed model spectrum",zorder=3)
        ax[4].scatter(point_lambda_ori.wavelength, point_lambda_ori.flux, marker="D", c="r", label="model flux", s=25, zorder=4)
        xerr, yerr =point_lambda_obs_ori.bandwidth/2,point_lambda_obs_ori.fluxerr
        ax[4].errorbar(point_lambda_obs_ori.wavelength, point_lambda_obs_ori.flux,xerr=xerr,yerr=yerr,fmt="o", c="b", ms=4, label="observed flux", zorder=5)
        ax[4].scatter(point_lambda_obs.wavelength, point_lambda_obs.flux, \
        marker="o", facecolors='none', edgecolors='purple', s=90, linewidth=2,label="observed flux for fitting", zorder=6)
        ax[4].set_xlim(0.1,40)
        ax[4].set_xscale("log")
        ax[4].legend(fontsize=8)
        ax[4].set_xlabel(r"$\rm Wavelength\ (\mu m)$", fontsize=14)
        ax[4].set_ylabel(r"$\rm Flux\ (erg/cm^2/s/\AA)$", fontsize=14)


        if LGassianFit: 
            # Teff
            width = fit_hist.loc[1, "teff_bin"]-fit_hist.loc[0, "teff_bin"]
            ax[5].bar(fit_hist.loc[:, "teff_bin"], fit_hist.loc[:, "teff_prob"], edgecolor="k", alpha=0.5, lw=0.5, width =width)
            ylim = ax[5].get_ylim()
            ylim = ax[5].get_ylim()
            ax[5].vlines([optimal.teff, optimal.teff-optimal.tefferr, optimal.teff+optimal.tefferr], ylim[0], ylim[1], lw=2.5, color="r", ls="--")
            ax[5].set_xlabel("Effective temperature (K)", fontsize=14)
            ax[5].set_ylabel("Probability density", fontsize=14)

            # Av
            width = fit_hist.loc[1, "Av_bin"]-fit_hist.loc[0, "Av_bin"]
            ax[6].bar(fit_hist.loc[:, "Av_bin"], fit_hist.loc[:, "Av_prob"], edgecolor="k", alpha=0.5, lw=0.5, width =width)
            ylim = ax[6].get_ylim()
            ax[6].vlines([optimal.Av, optimal.Av-optimal.Averr, optimal.Av+optimal.Averr], ylim[0], ylim[1], lw=2.5, color="r", ls="--")    
            ax[6].set_xlabel("Av (mag)", fontsize=14)
            ax[6].set_ylabel("Probability density", fontsize=14)

            # ScaleF
            width = fit_hist.loc[1, "ScaleF_bin"]-fit_hist.loc[0, "ScaleF_bin"]
            ax[7].bar(fit_hist.loc[:, "ScaleF_bin"], fit_hist.loc[:, "ScaleF_prob"], edgecolor="k", alpha=0.5, lw=0.5, width =width)
            ylim = ax[7].get_ylim()
            ax[7].vlines([optimal.ScaleF, optimal.ScaleF-optimal.ScaleFerr, optimal.ScaleF+optimal.ScaleFerr], ylim[0], ylim[1], lw=2.5, color="r", ls="--")        
            ax[7].set_xlabel(r"$\rm log(\theta^2)}$", fontsize=14)
            ax[7].set_ylabel("Probability density", fontsize=14)      

            # Fbol.
            width = fit_hist.loc[1, "fbol_bin"]-fit_hist.loc[0, "fbol_bin"]
            ax[8].bar(fit_hist.loc[:, "fbol_bin"], fit_hist.loc[:, "fbol_prob"], edgecolor="k", alpha=0.5, lw=0.5, width =width)
            ylim = ax[8].get_ylim()
            ax[8].vlines([optimal.fbol, optimal.fbol-optimal.fbolerr, optimal.fbol+optimal.fbolerr], ylim[0], ylim[1], lw=2.5, color="r", ls="--")        
            ax[8].set_xlabel("log (bolometric flux)", fontsize=14)
            ax[8].set_ylabel("Probability density", fontsize=14)

        plt.tight_layout()
        if Lsave: plt.savefig(fig_sedfit_path+str(obsParams.loc[ii, "starID"])+"."+figext)
        if Lshow: plt.show() 
        plt.close("all")

        # # if in the batch mode, combine all the figures into a hdf5 file.
        if Lbatch and L_fig_h5py:
            # Convert Matplotlib to NumpyArray        
            canvas = FigureCanvas(fig)
            canvas.draw()
            FigArray = np.array(canvas.renderer.buffer_rgba())
            # note that all the fingures will be saved into a single file. This file is generally very large.
            hf.create_dataset(str(obsParams.loc[ii, "starID"]), data=FigArray, dtype='uint8', compression="gzip", compression_opts=4)

        t2=time.time()
        print("Input stellar parameters: #{:d}-{:}, Teff={:.1f}K, logg={:.2f}, Av={:.2f} mag".format(ii, starID, teff, logg, obsParams.loc[ii, "Av"]))
        print("Best-fitting result: log(ScaleF)={:.3f}+/-{:.3f}, Av={:.3f} mag, log[Fbol (erg/cm^2/s)]={:.3f}".format(
            optimal.ScaleF.values[0], optimal.ScaleFerr.values[0],optimal.Av.values[0], optimal.fbol.values[0]))     
        print("=============================== result summary ================================\n")

######################################################################################################################################################
# save files.
if Lbatch:
    # assign the job index to the filename
    # Best-fitting result without interpolation
    filename = batch_data_sedfit_path+str(job_index)+".BestfitModels.csv"
    if Lsave: FileFinal.to_csv(filename,  float_format="%.5e", index=False)

    # Best-fitting result with interpolation
    filename = batch_data_sedfit_path+str(job_index)+".interpBestModel.csv"
    if Lsave: FileOptimal.to_csv(filename,  float_format="%.5e", index=False)

    filename = batch_data_sedfit_path+str(job_index)+".flux.csv"
    if Lsave: FileOut.to_csv(filename,  float_format="%.5e", index=False)

    # close hf
    if L_fig_h5py: hf.close()