import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib
from lmfit import Minimizer, Parameters, report_fit
import time
import sys    
import copy
import SedFitFunctions
from inlist import *
import inlist
import os 
from phomKit import ExtLaws, Saturation_Reject

# Some general remarks.
# This script is used to exclude outlier photometry by comparing observed and predicted flux densites. 
# The latter is derived by fitting blackbody to the observed SED. Note that in this circumstance, extinction 
# can be simply set to zero, because it would not affect rejecting photometric outliers. Effective temperature  
# and the scaling factor optimized from this fitting are of course biased accordingly.  

# The steps of this script is as follows:
# (1) load input photometry and convert magnitudes to fluxes.
# (2) exclude invalid data points by fitting a planck function.
# (3) determine Av. This can be done through a blackbody fitting. Pre-known can be loaded, or simply set to be 0.
# (4) with Av determined above, optimize Teff    
# (5) save predicted flux densities to files, and save figures.



if Lbatch==True:
    matplotlib.use('agg')
    index_start, index_end, job_index = int(sys.argv[1]), int(sys.argv[2])+1, int(sys.argv[3])
else:
    matplotlib.use('Qt5Agg')



# Only used to read in wavelength. Read in the filter data.
starfiles = np.array(["p5000_g+3.5_m0.0_t01_st_z-0.25_a+0.10_c+0.00_n+0.00_o+0.10_r+0.00_s+0.00.flx"])
filters = pd.read_csv(filter_zero_point_path).set_index('band')
filters.loc[:, "Alambda_Over_Av"] = ExtLaws(filters.loc[:, "lambda"].reset_index(), law=ExtinctionLaw, Rv=Rv)
photo = pd.read_csv(file_input_photometry, dtype={"starID":str})
params = pd.read_csv(file_input_params, dtype={"starID":str})
if not (params.starID==photo.starID).all(): 
    print("The tables of input parameters, input photometry, and blackbody fits are not consistent, thus exit.")
    exit()

# ensure in the dataframe table, the photometric error column for a given band immediately follows the magnitude column.
subs = np.where(photo.columns!="starID")[0]
photo = photo[np.append("starID", photo.columns[subs].sort_values())].copy()

# Save fitting results
columns =  ["starID", "Qflg", "Av", "ScaleF", "ScaleFerr", "Teff", "Tefferr", "RedChi"]
result = pd.DataFrame(index=range(len(photo.starID)), columns=columns)
result["Qflg"] = np.repeat("Y"*((len(photo.columns)-1)//2), len(photo))

# used to store fitting result for a batch.
FileBlkBdFits = pd.DataFrame()
index_end=np.min([index_end, len(photo)])

# (3) Fit Planck functions
for ii in range(index_start, index_end):     
    ##############################################################################################################################
    # load input photometry and convert magnitudes to fluxes.
    t0 = time.time()
    result.loc[ii, "starID"] = photo.loc[ii, "starID"]
    
    # Fetch observed SEDs. 
    photometry = pd.DataFrame(columns=["band", "mag", "magerr"])
    photometry["band"] = photo.iloc[ii, 1::2].index
    photometry["mag"] = photo.iloc[ii, 1::2].values
    photometry["magerr"] = photo.iloc[ii, 2::2].values
    photometry[["mag", "magerr"]] = photometry[["mag", "magerr"]].apply(pd.to_numeric, errors='coerce', axis=1)
    
    # select filters whenever available.
    photometry = photometry.set_index('band')    
    photometry = pd.merge(photometry, filters, on="band").sort_values(by=["lambda"])  

    # read in the wavelength and flux of MARCS models.
    # spec_lambda = SedFitFunctions.readspec(marcs_spectrum_path, starfiles[0])
    modelFluxDensity = pd.read_csv(grid_list_path+"PassbandFlux_"+grid+".csv") 
    spec_lambda = SedFitFunctions.readspec(grid_spectrum_path, modelFluxDensity.loc[0, "filename"], grid=grid)

    # convert magnitudes to flux densities.
    point_lambda_obs_ori = SedFitFunctions.mag2fluxdensity(photometry.copy(), Av=0, Lfluxnu=False)
    point_lambda_obs_ori = point_lambda_obs_ori.apply(pd.to_numeric, errors='coerce', axis=1)
    point_lambda_obs_ori.loc[:, "logFlux"] = np.log10(point_lambda_obs_ori.flux)
    point_lambda_obs_ori.loc[:, "logFluxerr"] = point_lambda_obs_ori.fluxerr/point_lambda_obs_ori.flux/np.log(10)
    # print("Input fluxes:")
    # print(point_lambda_obs_ori)


    # check if magnitudes are saturated.
    photometry = Saturation_Reject(photometry)
    photometry_ori = photometry.copy()


    # flag invalid data points. if a photometric measurement is invalid, set Qflg to be "N", otherwise "Y".
    for nn in range(len(photometry_ori.mag)):
        if ~np.isfinite(photometry_ori.loc[photometry_ori.index[nn], "mag"]):
            Qflg = list(result.loc[ii, "Qflg"])
            Qflg[nn] = "N"
            result.loc[ii, "Qflg"] = "".join(Qflg)

    # remove invalid data points, e.g. NaN, prepared for the fitting.
    photometry = photometry[~np.isnan(photometry.mag)]
    band_all = photometry_ori.index.values
    bands_dropped = ""
    if len(photometry)<5: 
        print("Only {:d} photometric measurements. As < 4 (few measurements or more available but saratued), skip this star".format(len(photometry)))
        continue

    ##############################################################################################################################
    # exclude invalid data points by fitting a planck function.
    while True: 
        # Convert magnitude to flux density.
        point_lambda_obs = SedFitFunctions.mag2fluxdensity(photometry.copy(), Av=0, Lfluxnu=False)

        # Do the fit 
        WaveObs = np.array(point_lambda_obs.wavelength/1e6, dtype="float") # wavelength in units of meter
        FluxObs = np.array(point_lambda_obs.flux, dtype="float")
        
        # After checking results, I find giving same weights to all the photometric measurements works better to removed photometric outliers.
        # FluxErrObs = np.array(point_lambda_obs.fluxerr, dtype="float")
        FluxErrObs = np.ones_like(point_lambda_obs.flux.astype(float))
        BlackbodyScaleF, BlackbodyScaleFerr, BlackbodyTeff, BlackbodyTefferr, RedChi = \
            SedFitFunctions.Fit_BlackBody(WaveObs, FluxObs, Teffmax=Teffmax, Teffmin=Teffmin, err=FluxErrObs)
        
        # Generate a black-body spectrum given wavelength and effective temperature. Here extinction is ingored.
        BlackbodyFlux = SedFitFunctions.Flux_Blackbody(spec_lambda.wavelength/1e6, BlackbodyTeff)*BlackbodyScaleF
        Blackbody = pd.DataFrame({"wavelength":spec_lambda.wavelength, "flux":BlackbodyFlux})
            
        # Now check the fitting and reject bad data point.
        fluxobs = np.array(point_lambda_obs.flux, dtype=float)
        fluxpred = SedFitFunctions.Flux_Blackbody(np.array(point_lambda_obs.wavelength, dtype=float)/1e6, BlackbodyTeff)*BlackbodyScaleF
        
        # Show
        LAvShow=False
        if LAvShow==True:
            fig, ax = plt.subplots(1,1)
            ax.plot(WaveObs, fluxobs, marker="D", c="r", label="model flux density")
            ax.plot(WaveObs, fluxpred, marker="o", c="b",label="observed flux density")
            ax.legend()
            ax.set_xscale("log")
            plt.grid(b=True, which='minor')
            plt.grid(b=True, which='major')
            plt.tight_layout()
            plt.show()

        # if relative difference > 0.25, discard this data point.
        if max(abs(fluxobs-fluxpred)/BlackbodyFlux.max())>0.25:
            ind = np.argmax(abs(fluxobs-fluxpred)/BlackbodyFlux.max())   
            if photometry.loc[photometry.index[ind], "lambda"]<1.0: 
                subs = np.where(photometry_ori.index==photometry.index[ind])[0][0]
                Qflg = list(result.loc[ii, "Qflg"])
                Qflg[subs] = "N"
                result.loc[ii, "Qflg"] = "".join(Qflg)

                # record drop bands
                bands_dropped += photometry.index[ind] + "|"
                photometry = photometry.drop([photometry.index[ind]])
                continue
        break    
    
    # update flux densities using valid photometric data.
    point_lambda_obs = SedFitFunctions.mag2fluxdensity(photometry.copy(), Av=0, Lfluxnu=False)
    Wave = np.array(point_lambda_obs.wavelength/1e6, dtype="float") # wavelength in units of meter
    ExtLaw = np.array(photometry["Alambda_Over_Av"])
    print("Photometric data points dropped: ", bands_dropped)

    ##############################################################################################################################
    # determine Av.
    if Optimize_Av==True:
        Av = np.arange(0, 1.0, 0.01)
        FitBD = pd.DataFrame(index=range(len(Av)), columns=["Av", "ScaleF", "ScaleFerr", "Teff", "Tefferr", "RedChi", "RltDiff"])
        FitBD.loc[:, "Av"] = Av
        DiffMin = np.zeros_like(Av)

        for mm in range(len(Av)):
            Flux = np.array(point_lambda_obs.flux, dtype="float")*10**(0.4*ExtLaw*Av[mm])
            Fluxerr = np.array(point_lambda_obs.fluxerr, dtype="float")*10**(0.4*ExtLaw*Av[mm])
            FitBD.iloc[mm:, 1:-1] = SedFitFunctions.Fit_BlackBody(Wave, Flux, Teffmax=Teffmax, Teffmin=Teffmin, err=Fluxerr)

            # Generate the Planck function with optimal parameters.
            PlanckFlux = SedFitFunctions.Flux_Blackbody(spec_lambda.loc[:, "wavelength"]/1e6, FitBD.loc[mm, "Teff"])*FitBD.loc[mm, "ScaleF"]
            Planck = pd.DataFrame({"wavelength":spec_lambda.loc[:, "wavelength"], "flux":PlanckFlux})

            # Calculate the maximum relative flux difference.
            fluxobs = np.array(point_lambda_obs.flux, dtype=float)
            waveobs = np.array(point_lambda_obs.wavelength, dtype=float)/1e6
            fluxpred = SedFitFunctions.Flux_Blackbody(waveobs, FitBD.loc[mm, "Teff"])*FitBD.loc[mm, "ScaleF"]
            FitBD.loc[mm, "RltDiff"] = np.sum(abs(fluxpred-fluxobs)/fluxobs)

            LAvShow=False
            if LAvShow==True:
                if mm%5==0:
                    fig, ax = plt.subplots(1,1)
                    ax.plot(Planck.wavelength, Planck.flux)
                    ax.scatter(point_lambda_obs.wavelength, Flux, marker="D", c="k")
                    ax.set_xscale("log")
                    plt.grid(b=True, which='minor')
                    plt.grid(b=True, which='major')
                    plt.tight_layout()

        # Do the final fitting.
        subs = np.argmin(np.array(FitBD["RltDiff"]))  
        OptimalAv = FitBD.loc[subs,"Av"]
    
    elif LFixedAv:
        # load known extinctions
        OptimalAv = params[photo.loc[ii, "starID"]==params.starID].Av.values[0]
    else:
        OptimalAv = 0

    ##############################################################################################################################
    # with Av determined, now optimize Teff    
    Flux = np.array(point_lambda_obs.flux, dtype="float")*10**(0.4*ExtLaw*OptimalAv)
    PlankScaleF, PlankScaleFerr, PlankTeff, PlankTefferr, RedChi  = SedFitFunctions.Fit_BlackBody(Wave, Flux, Teffmax=Teffmax, Teffmin=Teffmin, err=FluxErrObs)
    PlanckFlux = SedFitFunctions.Flux_Blackbody(spec_lambda.loc[:, "wavelength"]/1e6, PlankTeff)*PlankScaleF
    Planck = pd.DataFrame({"wavelength":spec_lambda.loc[:, "wavelength"], "flux":PlanckFlux})
    fluxobs = np.array(point_lambda_obs.flux, dtype=float)
    waveobs = np.array(point_lambda_obs.wavelength, dtype=float)/1e6
    fluxpred = SedFitFunctions.Flux_Blackbody(waveobs, PlankTeff)*PlankScaleF

    # W3 and W4 photometry can be rejected for fitting.
    LexcludeW3W4=False
    if LexcludeW3W4==True:
        photometry_final = photometry[photometry["lambda"]<10].reset_index(drop=True)
        point_lambda_obs_final = SedFitFunctions.mag2fluxdensity(photometry_final, Av=0, Lfluxnu=False)
        # Do the fit 
        WaveObs = np.array(point_lambda_obs_final.wavelength/1e6, dtype="float") # wavelength in units of meter
        FluxObs = np.array(point_lambda_obs_final.flux, dtype="float")
        BlackbodyScaleF, BlackbodyScaleFerr, BlackbodyTeff, BlackbodyTefferr, RedChi = SedFitFunctions.Fit_BlackBody(WaveObs, FluxObs, Teffmax=Teffmax, Teffmin=Teffmin)
        # Generate a black-body spectrum given wavelength and effective temperature. Here extinction is ingored.
        BlackbodyFlux = SedFitFunctions.Flux_Blackbody(spec_lambda.wavelength/1e6, BlackbodyTeff)*BlackbodyScaleF
        Blackbody_final = pd.DataFrame({"wavelength":spec_lambda.wavelength, "flux":BlackbodyFlux})
        # Now check the fitting and reject bad data point.
        fluxobs = np.array(point_lambda_obs.flux, dtype=float)
        fluxpred = SedFitFunctions.Flux_Blackbody(np.array(point_lambda_obs.wavelength, dtype=float)/1e6, BlackbodyTeff)*BlackbodyScaleF
        # save predicted flux densities.
        point_lambda_obs.loc[:, "fluxobs"] = fluxobs
        point_lambda_obs.loc[:, "fluxpred"] = fluxpred
        if Lsave: point_lambda_obs[["fluxobs", "fluxpred"]].to_csv(save_blackbody_flux_path_single+str(photo.loc[ii, "starID"])+".cut.csv", float_format="%.4e")

    # ##############################################################################################################################
    # save predicted flux densities to files and save figures.
    if LOutput_preditcted_flux==True:
        subs = [jj=="Y" for jj in result.loc[ii, "Qflg"]]
        bands = photometry_ori.index[subs]
        out = pd.DataFrame(index=range(len(fluxobs)), columns=["band", "fluxobs", "fluxpred"])
        out.loc[:, "band"] = bands
        out.loc[:, "fluxobs"] = fluxobs
        out.loc[:, "fluxpred"] = fluxpred
        if Lsave: out.to_csv(save_blackbody_flux_path_single+str(photo.loc[ii, "starID"])+".flux.csv", index=False, float_format="%.4e")
    if LOutput_optimized_Teff_Av_ScaleFactor==True:
        # save the fitting results.
        result.loc[ii, "Av"] = OptimalAv
        result.loc[ii, "ScaleF"] = PlankScaleF
        result.loc[ii, "ScaleFerr"] = PlankScaleFerr
        result.loc[ii, "Teff"] = PlankTeff
        result.loc[ii, "Tefferr"] = PlankTefferr
        result.loc[ii, "RedChi"] = RedChi
        # save
        temp = pd.DataFrame(columns=result.columns)
        temp.loc[0]=result.iloc[ii,:]
        # save result from running in a batch model.
        if Lbatch:
            FileBlkBdFits = pd.concat([FileBlkBdFits, temp], axis=0, ignore_index=True)
        else:    
            if Lsave: temp.to_csv(save_blackbody_flux_path_single+str(result.loc[ii, "starID"])+".csv",index=False,float_format="%.4e")

    # plot the fitting
    fontsize=16
    fig, ax = plt.subplots(1,1, figsize=(8,6))
    # orginal data
    ax.errorbar(point_lambda_obs_ori.wavelength, point_lambda_obs_ori.flux, \
        xerr=point_lambda_obs_ori.bandwidth/2, yerr=point_lambda_obs_ori.fluxerr, \
        fmt="o", c="k",  zorder=1, label="observed  flux density, outliers")
    # without considering extinctions.
    ax.errorbar(point_lambda_obs.wavelength, point_lambda_obs.flux, xerr=point_lambda_obs.bandwidth/2, \
        yerr=point_lambda_obs.fluxerr,fmt="o", c="r",  zorder=2, label="observed  flux density")
    ax.plot(Blackbody.wavelength, Blackbody.flux, c="r", lw=2, zorder=3, label="Best-fitting planck, undeextincted")
    # considering extinctions.
    ax.errorbar(Wave*1e6, Flux, xerr=point_lambda_obs.bandwidth/2, yerr=point_lambda_obs.fluxerr,fmt="o", c="b", \
        zorder=4, label="deextincted  flux density")
    ax.plot(Planck.wavelength, Planck.flux, c="b", lw=2, zorder=5, label="Best-fitting planck, deextincted")
    # layout adjustment
    ax.set_xlim(0.1,30)
    ax.set_xscale("log")
    ax.set_xlabel("Wavelength (micron)", fontsize=fontsize)
    ax.set_ylabel("Flux (erg/cm2/s/A)", fontsize=fontsize)
    # reset items order in legend.
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[2], handles[0], handles[3], handles[1], handles[4]]
    labels = [labels[2], labels[0], labels[3], labels[1], labels[4]]
    ax.legend(handles,labels,fontsize=fontsize-6)
    plt.grid(b=True, which='minor')
    plt.grid(b=True, which='major')
    plt.title(str(photo.loc[ii, "starID"]))
    plt.tight_layout()
    if Lsave: plt.savefig(save_blackbody_plot_path+str(photo.loc[ii, "starID"])+"."+figext)  
    plt.show() if Lshow==True else plt.close("all") 
    t1 = time.time()
    print("#{:d}-{:}, Av={:.2f}, Teff={:.1f}, {:.3f}s consumed".format(ii, photo.loc[ii, "starID"], OptimalAv, PlankTeff, t1-t0)) 



if Lbatch:
    # assign the job index to the filename
    filename = save_blackbody_flux_path_batch+str(job_index)+".csv"
    if Lsave: FileBlkBdFits.to_csv(filename,  float_format="%.5e", index=False)
