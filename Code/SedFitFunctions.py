import glob
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt  
import matplotlib
from lmfit import Minimizer, Parameters, report_fit
import time
import re
import math
import warnings
from collections import deque
from astropy.io import ascii
from astropy.table import Table
import sys
import pylab
import json
import requests
from urllib.parse import quote as urlencode
from urllib.request import urlretrieve
import http.client as httplib 
import astropy.coordinates as coords
from astroquery.vizier import Vizier
import astropy.units as units
import astropy.units as u
from uncertainties import unumpy 



def readspec(marcs_path, filename, grid="MARCS"):
    if grid=="MARCS":
        # read in the wavelength (in units of micron) and flux (in units of erg/cm^2/s/A) of MARCS models.
        wave = pd.read_csv(marcs_path+"flx_wavelengths.vac", header=None, names=["wavelength"])
        flux = pd.read_csv(marcs_path+filename, header=None, names=["flux"])
        spec = pd.DataFrame({"wavelength":np.array(wave.wavelength/1e4), "flux":np.array(flux.flux)})
    elif grid=="BOSZ":
        spec = Table.read(marcs_path+filename).to_pandas()
        spec = spec.drop(columns=["Continuum"]).rename(columns={"Wavelength":"wavelength", "SpecificIntensity":"flux"})
        spec.loc[:, "wavelength"] /= 1e4
        spec.loc[:, "flux"] *= np.pi
    else:
        print("The grid is not resolved. CHECK the grid")
    return spec

def resolution(atmos, res=200):
    "smooth original MARCS spectra"
    model = atmos.copy()
    sampfactor = int(20000/res)
    model["flux_sm"] = model["flux"].rolling(sampfactor, win_type='boxcar').mean()
    return model

def flux_extrap(spec, Minwavelength=8, Maxwavelength=30, order=3):
    """extrapolate flux from 20 to 30 micron for MARCS spectra whose max wavelength is 20 micon
       spec: it's an input spectrum, which is a pandas DataFrame with two columns: wavelength and flux
       Minwavelength: min wavelength of a subset of the input spectrum for extrapolation.
       Maxwavelength: max wavelength of a subset of the input spectrum for extrapolation.
       order: the order of the polynomial used for extrapolation.
       Output: original spectrum with the extrapolated spectrum concatenated.
    """
    extrap = spec[spec.wavelength>Minwavelength].reset_index()
    extrap_wave_step = np.median(np.array(extrap.loc[1:, "wavelength"])-np.array(extrap.loc[0:len(extrap)-2, "wavelength"])) 
    extrap_wave = np.arange(np.min(extrap.wavelength), Maxwavelength, extrap_wave_step) 
    extrap_func = np.poly1d(np.polyfit(extrap.wavelength.values, np.log10(extrap.flux.values), order))
    extrap = pd.DataFrame({"wavelength":extrap_wave, "flux":10**(extrap_func(extrap_wave))})
    spec = pd.concat([spec, extrap[extrap.wavelength>np.max(spec.wavelength)].reset_index(drop=True)]).reset_index(drop=True)
    return spec

def Flambda2Fnu(atmos):
    """Change flux in erg/cm^2/s/anstrogm to flux in Jy, where 1 Jy = 10^(-23) erg/cm^2/s/Hz. Note that
       wavelength has the units of microns, instead of anstrogms"""
    model = atmos.copy()
    model["flux"] = (atmos["wavelength"]*1e4)**2 / 2.998 * 10**5 * atmos["flux"]
    return model

def Fnu2Flambda(atmos):
    """Change flux in Jy to that in erg/cm^2/s/anstrogm, where 1 Jy = 10^(-23) erg/cm^2/s/Hz"""
    model = atmos.copy()
    model["flux"] = atmos["flux"] / ((atmos["wavelength"]*1e4)**2 / 2.998 * 10**5) 
    return model


def pivotwave(bandpass=None, path=None):
    # collect response functions for individual filters.
    filterpaths = np.sort(glob.glob(path+"*.dat"))
    filternames = [os.path.basename(i)[:-4] for i in filterpaths]
    band = pd.DataFrame({"band":filternames, "path":filterpaths})
    filters = pd.read_csv("../Data/MetaData/filters.dat")[["band", "device"]]
    filters = pd.merge(filters, band, on="band")
    wave = np.zeros(len(bandpass))
    # calculate pivot wavelength using equation given on the page 396 (immediately above equa. 7) of Casagrande+2014
    for i in range(len(bandpass)):
        if np.isin(bandpass[i], filters.band):
            filename = filters[filters.band==bandpass[i]].path.iloc[0]
            device = filters[filters.band==bandpass[i]].device.iloc[0]            
            responseFunc = pd.read_csv(filename)
            if device=="energy":
                # for energy-counting devices, filter transmission is not mupltiplied by wavelength             
                numerator = np.trapz(responseFunc.response, x=responseFunc.wavelength)
                denominator = np.trapz(responseFunc.response/responseFunc.wavelength**2, x=responseFunc.wavelength)                
            else:    
                # for photon-counting devices, filter transmission must be mupltiplied by wavelength 
                numerator = np.trapz(responseFunc.response*responseFunc.wavelength, x=responseFunc.wavelength)
                denominator = np.trapz(responseFunc.response/responseFunc.wavelength, x=responseFunc.wavelength)    
            wave[i] = (numerator/denominator)**0.5
            # For extra broad filters, the diffence can be non-negligeble. 
        else:
            print("The bandpass {:s} has not filter profiles prepared, thus the computation can be performed.".format(bandpass[i]))   
            wave[i] = None
    return wave

def FluxDensity(atmos, bandpass=None, path=None):
    """the flux of the MARCS model atmosphere is in units of erg/cm^2/s/anstrogm, and 
       the wavelength in the units of anstrogm.
    """
    filterpaths = np.sort(glob.glob(path+"*.dat"))
    filternames = [os.path.basename(i)[:-4] for i in filterpaths]
    band = pd.DataFrame({"band":filternames, "path":filterpaths})
    filters = pd.read_csv("../Data/MetaData/filters.dat")[["band", "device"]]
    filters = pd.merge(filters, band, on="band")
    fluxes = np.zeros(len(bandpass))
    for i in range(len(bandpass)):
        if np.isin(bandpass[i], filters.band):
            filename = filters[filters.band==bandpass[i]].path.iloc[0]
            device = filters[filters.band==bandpass[i]].device.iloc[0]            
            responseFunc = pd.read_csv(filename)
            # force a filter transmission curve to have the same wavelength resolution as the spectrum. 
            # This is recommmended by Bessell2012.
            # flux_interp = np.interp(responseFunc.wavelength, atmos.wavelength, atmos.flux)
            response_interp = np.interp(atmos.wavelength, responseFunc.wavelength, responseFunc.response)
            atmos.loc[:, "response_interp"] = response_interp
            atmos_trunc = atmos[(atmos.wavelength>=responseFunc.wavelength.min()) & (atmos.wavelength<=responseFunc.wavelength.max())]  

            if device=="energy":
                # for energy-counting devices, filter transmission is not mupltiplied by wavelength             
                # numerator = np.trapz(responseFunc.response*flux_interp, x=responseFunc.wavelength)
                # denominator = np.trapz(responseFunc.response, x=responseFunc.wavelength)  
                numerator = np.trapz(atmos_trunc.response_interp*atmos_trunc.flux, x=atmos_trunc.wavelength)
                denominator = np.trapz(atmos_trunc.response_interp, x=atmos_trunc.wavelength)                          
            else:    
                # for photon-counting devices, filter transmission must be mupltiplied by wavelength 
                # numerator = np.trapz(responseFunc.response*flux_interp*responseFunc.wavelength, x=responseFunc.wavelength)
                # denominator = np.trapz(responseFunc.response*responseFunc.wavelength, x=responseFunc.wavelength)
                numerator = np.trapz(atmos_trunc.response_interp*atmos_trunc.wavelength*atmos_trunc.flux, x=atmos_trunc.wavelength)
                denominator = np.trapz(atmos_trunc.response_interp*atmos_trunc.wavelength, x=atmos_trunc.wavelength)
            fluxes[i] = numerator/denominator
            # For extra broad filters, the diffence can be non-negligeble. 
        else:
            print("The bandpass you inputed has not filter profiles prepared, thus the computation can be performed.")     
            fluxes[i] = None
    return fluxes   

def GaiaMag(atmos, bandpass=None, path=None):
    """
    Calculate Gaia G, BP, RP magnitudes from an spectrum for which the units of the flux is erg/cm^2/s/A and 
    the units of wavelength is micron. 
    """
    PA = 0.7278 #m2
    c = 2.997924581e8 #m/s
    h = 6.62607004e-34# m2 kg / s
    const = PA / 1e9 / h / c 
    # load filter profiles
    filterpaths = np.sort(glob.glob(path+"*.dat"))
    filternames = [os.path.basename(i)[:-4] for i in filterpaths]
    band = pd.DataFrame({"band":filternames, "path":filterpaths})
    filters = pd.read_csv("../Data/MetaData/filters.dat")[["band", "device"]]
    filters = pd.merge(filters, band, on="band")
    mag = np.ones(len(bandpass))
    for i in range(len(bandpass)):
        if np.isin(bandpass[i], filters.band):
            filename = filters[filters.band==bandpass[i]].path.iloc[0]
            device = filters[filters.band==bandpass[i]].device.iloc[0]            
            responseFunc = pd.read_csv(filename)
            # force a filter transmission curve to have the same wavelength resolution as the spectrum. This is recommmended by Bessell2012.
            response_interp = np.interp(atmos.wavelength, responseFunc.wavelength, responseFunc.response)
            atmos.loc[:, "response_interp"] = response_interp
            atmos_trunc = atmos[(atmos.wavelength>=responseFunc.wavelength.min()) & (atmos.wavelength<=responseFunc.wavelength.max())]  
            flux_mean = np.trapz(atmos_trunc.response_interp*atmos_trunc.wavelength*1e3*atmos_trunc.flux*1e-2, x=atmos_trunc.wavelength*1e3)
            if bandpass[i]=="gaiag": mag[i] = -2.5*np.log10(const * flux_mean) + 25.6874
            if bandpass[i]=="gaiabp": mag[i] = -2.5*np.log10(const * flux_mean) + 25.3385
            if bandpass[i]=="gaiarp": mag[i] = -2.5*np.log10(const * flux_mean) + 24.7479
    return mag

def mag2fluxdensity(photometry, Lfluxnu=True, Av=0.0):    
    fluxdensity = pd.DataFrame(index=range(len(photometry.index)), columns=["band", "wavelength", "flux", "fluxerr"])
    fluxdensity["band"] = photometry.index
    fluxdensity= fluxdensity.set_index('band')
    for i in photometry.index: 
        # convert magnitudes to flux densities per frequency unit
        if photometry.loc[i, "magsystem"]=="Vega":
            fluxdensity.loc[i, "flux"] = 10**(photometry.loc[i, "mag"]/(-2.5))*photometry.loc[i, "fnu0"]
        if photometry.loc[i, "magsystem"]=="AB":
            fluxdensity.loc[i, "flux"] = 10**(photometry.loc[i, "mag"]/(-2.5))*3631
        # convert flux densities from per frequency unit to per wavelength unit, namely flux_nu to flux_lambda. 
        # Noth that the ref wavelength is in units of micron.
        if Lfluxnu==False: 
            fluxdensity.loc[i, "flux"] = fluxdensity.loc[i, "flux"]*2.998/1e13/photometry.loc[i, "lambda"]**2
        # collect the wavelength and bandwidth
        fluxdensity.loc[i, "wavelength"] = photometry.loc[i, "lambda"].copy()
        fluxdensity.loc[i, "bandwidth"] = photometry.loc[i, "bandwidth"].copy()
        fluxdensity.loc[i, "Alambda_Over_Av"] = photometry.loc[i, "Alambda_Over_Av"].copy()
        # correct extinction, assuming the extinction is known.
        fluxdensity.loc[i, "flux"] *= 10**(0.4*photometry.loc[i, "Alambda_Over_Av"]*Av)
        # evaluate flux density uncertainty with error propogation. Note that flux density 
        # uncertainty can be calculated via the same equation, no matter what magnitude system is.
        # |df| = ln10/2.5*f*dm
        if photometry.loc[i, "magerr"]  is not np.ma.masked:
            fluxdensity.loc[i, "fluxerr"] = fluxdensity.loc[i, "flux"]*np.log(10)/2.5*photometry.loc[i, "magerr"]
        else:
            fluxdensity.loc[i, "fluxerr"] = np.nan
        # flux densities in units of erg/cm^2/s/A    
    return fluxdensity

def ClipPhotometry(band_clip, band_all, Qflgs_all, photometry):
    if band_clip in photometry.index: 
        # set Qflg to be "N"
        subs = np.where(band_all==band_clip)[0][0]
        if subs==0:
            Qflgs_all = "N"+Qflgs_all[subs+1:]
        elif subs==len(band_all)-1:
            Qflgs_all = Qflgs_all[:subs]+"N"
        else:
            Qflgs_all = Qflgs_all[:subs]+"N"+Qflgs_all[subs+1:] 
        photometry = photometry.drop(band_clip, axis=0)
    return photometry, Qflgs_all

# Fit black-body spectra to multi-wavelength photometry to estimate Teff.
def Flux_Blackbody(wave, T):
    """Calculate the flux of a black-body radiation. Note that the intensity and
    the flux differ by a factor of pi, in the sense of I_lambda * pi = Flux_lambda
    (see Page 123 of the textbook Stellar Photospheres by David F. Gary)"""
    # define constant
    h = 6.626070e-34  # in units of J/s
    c = 2.997925e+8     # in units of m/s
    k = 1.380649e-23   # in units of J/K
    # intensity has the units of erg/s/cm^2/cm
    a = 2.0*h*c**2
    b = h*c/(wave*k*T)
    intensity = a/ ( (wave**5) * (np.exp(b) - 1.0) ) * 10
    # flux has the units of erg/s/cm^2/A
    flux = np.pi * intensity /1e8
    return flux

def Resi_BlackBody(params, wave, f, err):
    alpha = params["alpha"]
    beta = params["beta"]    
    # model = alpha * 2.0*h*c**2/ ( (wave**5) * (np.exp(h*c/(wave*k*beta)) - 1.0) ) * 10 * np.pi /1e8  
    model = alpha * Flux_Blackbody(wave, beta)
    return model-f if err is None else (model-f)/err

def Fit_BlackBody(x, f, Teffmax=0, Teffmin=0, err=None):
    params = Parameters()
    params.add("alpha", value=10**(-17))
    params.add("beta", value=4000, min=Teffmin, max=Teffmax)
    minner = Minimizer(Resi_BlackBody, params, fcn_args=(x, f, err))
    out = minner.minimize()
    alpha = out.params["alpha"].value
    beta = out.params["beta"].value
    alphaerr = out.params["alpha"].stderr
    betaerr = out.params["beta"].stderr
    redchi = out.redchi
    return alpha, alphaerr, beta, betaerr, redchi

def Resi_ScaleFactor(params, x, y, err):
    alpha = params["alpha"]
    model = x+alpha
    return model-y if err is None else (model-y)/err

def Fit_ScaleFactor(scalef, x, y, err=None):
    params = Parameters()
    params.add("alpha", value=scalef)
    minner = Minimizer(Resi_ScaleFactor, params, fcn_args=(x, y, err))
    out = minner.minimize()
    alpha = out.params["alpha"].value
    alphaerr = out.params["alpha"].stderr
    redchi = out.redchi
    return alpha, alphaerr, redchi

def Resi_Fit_Extinction_ScaleFactor(params, x, y, yerr):
    alpha = params["alpha"]
    beta = params["beta"]
    model = -0.4*alpha*x + beta
    return model-y if yerr is None else (model-y)/yerr

def Fit_Extinction_ScaleFactor(x,y,yerr=None):
    params = Parameters()
    params.add("alpha", value=1.0, min=-5, max=20)
    params.add("beta", value=1.0)
    minner = Minimizer(Resi_Fit_Extinction_ScaleFactor, params, fcn_args=(x, y, yerr))
    out = minner.minimize()
    alpha = out.params["alpha"].value
    alphaerr = out.params["alpha"].stderr
    beta = out.params["beta"].value
    betaerr = out.params["beta"].stderr    
    redchi = out.redchi
    return alpha, alphaerr, beta, betaerr, redchi


def Gaussian(x, H, P, W, N):
    return N + H * np.exp(-((x-P)/W)**2/2)

def Resi_Gaussian(params, x, y, err):
    H = params["H"]
    P = params["P"]
    W = params["W"]
    N = params["N"]
    model = Gaussian(x, H, P, W, N)
    return model-y if err is None else (model-y)/err

def Fit_Gaussian(x, y, param_guess, err=None):
    params = Parameters()
    params.add("H", value=param_guess[0], min=0)
    params.add("P", value=param_guess[1])
    params.add("W", value=param_guess[2], min=0)
    params.add("N", value=param_guess[3], min=0)
    minner = Minimizer(Resi_Gaussian, params, fcn_args=(x, y, err))
    out = minner.minimize()
    H = out.params["H"].value
    P = out.params["P"].value
    W = out.params["W"].value
    N = out.params["N"].value
    Herr = out.params["H"].stderr
    Perr = out.params["P"].stderr
    Werr = out.params["W"].stderr    
    Nerr = out.params["N"].stderr  
    return H, Herr, P, Perr, W, Werr, N, Nerr

def CalcLumi(dist=None, disterr=None, Fbolo=None, Fboloerr=None):
    """Calculate luminosity with distance and bolometric flux L=4*pi*d^2 * Fbol / Lsun
       distance: in the units of parsec
       bolometric flux: in the units of erg/s/cm^2 
       output: log(L/Lsun)
    """
    parsec2cm = 3.086e18
    solarlumi = 3.85e33  #https://link.springer.com/referenceworkentry/10.1007%2F1-4020-4520-4_374     
    x = unumpy.uarray(dist, disterr)
    y = unumpy.uarray(Fbolo, Fboloerr)
    z = unumpy.log10(4*np.pi*(x*parsec2cm)**2  * y/solarlumi)
    logL = unumpy.nominal_values(z)
    logLerr = unumpy.std_devs(z)
    return logL, logLerr

def radius_From_angRadius(angRadius=None, angRadiuserr=None, dist=None, disterr=None):
    """calculate radius from angular radius and distance"""
    # Input:
    # angRadius: milli arcsec
    # distance: in the units of parsec
    radiusSun = 6.957*1e10 # in units of cm, see Table 2 of Choi et al. 2016.
    parsec2cm = 3.086e18
    mas2unity = (np.pi/3.6/1.8)*1e-8
    x = unumpy.uarray(angRadius, angRadiuserr)
    y = unumpy.uarray(dist, disterr)
    z = ((x*mas2unity) * (y*parsec2cm))/radiusSun
    R = unumpy.nominal_values(z)
    Rerr = unumpy.std_devs(z)
    return R, Rerr


def gaussian_func(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))



def weighted_mean_std(x, prob=0):
    """Compute weighted mean and std given probability"""
    average = np.average(x, weights=prob)
    variance = np.average((x-average)**2, weights=prob)
    return (average, math.sqrt(variance))