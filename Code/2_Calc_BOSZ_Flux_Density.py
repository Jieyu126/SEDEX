import pandas as pd 
import glob
from astropy.table import Table
import matplotlib.pyplot as plt 
import re 
from SedFitFunctions import *
import numpy as np
from inlist import *
from joblib import Parallel, delayed


def fetch_filename(filename): return filename.split("/")[-1]

def fetch_params(filename):
    feh = -1*float(filename[2:5][1:])/10 if filename[2:5][0]=="m" else float(filename[2:5][1:])/10
    carbon = -1*float(filename[6:9][1:])/10 if filename[6:9][0]=="m" else float(filename[6:9][1:])/10
    alpha = -1*float(filename[10:13][1:])/10 if filename[10:13][0]=="m" else float(filename[10:13][1:])/10
    teff = int(re.search("t\d+(?=g)", filename).group(0)[1:])
    logg = float(re.search("g\d+(?=v)", filename).group(0)[1:])/10
    tv = float(re.search("v\d+(?=modrt)", filename).group(0)[1:])/10
    rot = float(re.search("modrt\d+(?=b)", filename).group(0)[5:])/10
    res = int(re.search(r"\d+(?=rs.fits)", filename).group(0))
    return [filename, teff, logg, feh, carbon, alpha, tv, rot, res]

def grid_params(filename, Lextract=False):
    """if extract=True, extract stellar parameters from the filename of each model spectrum.
       filename: the filename of bosz_wget file downloaded from mast:
       https://archive.stsci.edu/hlsp/bosz/search.php. Note that only [alpha/M]=0
       and [C/M]=0 models are downloaded.  
    """
    if Lextract:
        # load all the BOSZ models with c/m=0 and alpha/m=0.
        grid = pd.read_csv(filename, header=None, names=["scheme", "arg", "filename"], sep="\s+")
        grid.loc[:, "filename"] = grid["filename"].map(fetch_filename)
        columns = ["filename", "teff", "logg", "feh", "carbon", "alpha", "tv", "rot", "res"]
        params = pd.DataFrame(index=range(len(grid)), columns=columns)
        params.loc[:, :] = np.vstack(grid["filename"].map(fetch_params).values)
        params[columns[1:]] = params[columns[1:]].apply(pd.to_numeric, axis=1)
        # now correct rounded [M/H] values
        params.loc[:, 'feh'] = params['feh'].replace([-2.3, -1.8, -1.3, -0.8, -0.3, 0.3], [-2.25, -1.75, -1.25, -0.75, -0.25, 0.25]) 
        # now correct rounded [alpha/M] values
        params.loc[:, 'alpha'] = params['alpha'].replace([-0.3, 0.3], [-2.25, 0.25])         
    else:
        # directly load atmospheric parameters from the official table downloaded from MAST.
        params = pd.read_csv(filename)
        params.loc[:, "tv"] /= 10

    # add model index
    params.insert(1, "modelindex", np.arange(len(params)))
    # add geometry of the atmosphere, i.e., all plane parallel.
    params.loc[:, "geom"] = "p"
    # add a mass column, but set it to 0 to be consistent with the MARCS grid, since plane-parallel calculations does not rely on mass. 
    params.loc[:, "mass"] = 0
    params = params[['filename', 'modelindex', 'teff', 'logg', 'feh', 'carbon', 'alpha', 'tv', 'rot', 'res', 'geom', 'mass']]
    return params

def calFluxDensity(grid_path, filename, filter_profile_path, filters, grid="MARCS"):
    print("calculate flux densitis for the model: {:s}".format(filename))
    # read in the wavelength and flux of MARCS models.
    spec_ori = readspec(grid_path, filename, grid=grid)
    # convolve the model atmosphere with the filter transmission for each band chosen. The model spectrum has flux in 
    # units of erg/cm^2/s/Anstrogm. The resulting flux density for each band chosen are later used for minimisation.
    return np.append(filename, np.log10(FluxDensity(spec_ori, bandpass=filters.band.values.copy(), path=filter_profile_path)))

def calFbol(grid_path, filename, grid="MARCS"):
    print("calculate integrated emergent flux for the model: {:s}".format(filename))
    # read in the wavelength and flux of MARCS models.
    spec_ori = readspec(grid_path, filename, grid=grid)
    # extrapolate spectra from 20 to 30 micron for MARCS model spectra.
    if grid=="MARCS": spec_ori = flux_extrap(spec_ori.copy(), Minwavelength=8, Maxwavelength=30, order=3)
    # calculate bolometric flux on the stellar surface, in log10 scale.
    # if grid=="BOSZ":  spec_ori = spec_ori[spec_ori.wavelength>0.13].reset_index(drop=True)
    Fbol = np.log10(np.trapz(spec_ori.flux, x=spec_ori.wavelength*1e4))
    return filename, Fbol


##########################################################################################################
# load the parameters of the grid.
Ltest=False
filters = pd.read_csv("/home/yujie/Documents/SoftLinks/Scratch/SEDFit/Version/Data/MetaData/filters.dat")
grid_path = "/home/yujie/Documents/SoftLinks/Scratch/SpectrumGrid/BOSZ/Grid/"

# collect atmospheric parameters from the filenames of the spectrum models or directly load the table downloaded
# from MAST. Both methods yield the identical results.
models = grid_params("~/Documents/SoftLinks/Scratch/SpectrumGrid/BOSZ/Grid/bosz_wget.sh", Lextract=True)
# models = grid_params("~/Documents/SoftLinks/Scratch/SpectrumGrid/BOSZ/bosz.txt", Lextract=False)

# calculate flux densities for the entire grid
FluxDensity = Parallel(n_jobs=10)(delayed(calFluxDensity)(grid_path, \
filename, filter_transmission_path, filters, grid="BOSZ") for filename in  models.loc[:, "filename"])
FluxDensity = pd.DataFrame(data = np.vstack(FluxDensity))
FluxDensity.columns = np.append("filename", filters.band.values)
columns = FluxDensity.columns.values[1:]
FluxDensity[columns] = FluxDensity[columns].apply(pd.to_numeric)

# calculate integrated emergent flux for the entire grid
FbolFlux = Parallel(n_jobs=10)(delayed(calFbol)(grid_path, filename, grid="BOSZ") for filename in  models.loc[:, "filename"])
FbolFlux = pd.DataFrame(data = np.vstack(FbolFlux))
FbolFlux.columns = ["filename", "fbol"]
FbolFlux["fblo"] = FbolFlux["fbol"].apply(pd.to_numeric)

# combine atmospheric parameters and flux densities.
result = pd.merge(FbolFlux, FluxDensity, on="filename").reset_index(drop=True)
result = pd.merge(models, result, on="filename").reset_index(drop=True)
result.to_csv(grid_list_path+"PassbandFlux_BOSZ.csv", index=False, float_format="%.8f")

# table = pd.read_csv(grid_list_path+"PassbandFlux_BOSZ.csv")
# table = table.drop(columns=["fbol"])
# table = pd.merge(table, FbolFlux, on="filename")
# table.to_csv(grid_list_path+"PassbandFlux_BOSZ.csv", index=False, float_format="%.8f")

##########################################################################################################
# Compare MARCS and BOSZ model spectra.
##########################################################################################################
if Ltest:
    # (1) display the overall parameter grid.
    fig, ax = plt.subplots(1,1,figsize=(11,6))
    ax.scatter(models.teff, models.logg)
    ax.set_xlim(30500, 3000)
    ax.set_ylim(5.5, -0.5)
    ax.set_xlabel("Effective temperature (K)", fontsize=16)
    ax.set_ylabel("log g", fontsize=16)
    plt.tight_layout()
    plt.show()

    # (2) parameter grid for each metallicity.
    for i in np.unique(models.loc[:, 'feh']):
        print("for [M/H]={:.2f}, {} models available".format(i, len(models[models.feh==i])))
        fig, ax = plt.subplots(1,1,figsize=(11,6))
        ax.scatter(models[models.feh==-2.50].teff, models[models.feh==-2.50].logg, c="r")
        ax.scatter(models[models.feh==i].teff, models[models.feh==i].logg, c="k")
        ax.set_xlabel("Effective temperature (K)", fontsize=16)
        ax.set_ylabel("log g", fontsize=16)
        ax.set_xlim(30500, 3000)
        ax.set_ylim(5.5, -0.5)
        plt.tight_layout()
        plt.show()


    # (3) compare MARCS and BOSZ model spectra
    # load a MARCS model spectrum.
    # marcs_spectrum_path = "/home/yujie/Documents/SoftLinks/Scratch/SpectrumGrid/MARCS/Version3/Models/"
    marcs = pd.read_csv(grid_list_path+"PassbandFlux.csv")    
    marcs = marcs[(marcs.teff==6000)&(marcs.logg==3.5)&(marcs.feh==0)&(marcs.geom=="p")&(marcs.tv==2)].reset_index(drop=True)
    spec_marcs = readspec(marcs_spectrum_path, marcs.filename.values[0])
    # load a BOSZ model spectrum.
    bosz = models[(models.teff==6000)&(models.logg==3.5)&(models.feh==0)&(models.tv==2)].reset_index(drop=True)
    spec_bosz = readspec(grid_path, bosz.filename.values[2], grid="BOSZ")
    spec_bosz.loc[:, "flux_interp"] = np.interp(spec_bosz.wavelength, spec_marcs.wavelength, spec_marcs.flux)
    # compare integrated flux densities.
    Fbol_marcs=np.trapz(spec_marcs.flux, x=spec_marcs.wavelength)
    Fbol_bosz=np.trapz(spec_bosz.flux, x=spec_bosz.wavelength)
    print("ratio of integrated flux densities (MARCS/BOSZ): {:.6f}".format(Fbol_marcs/Fbol_bosz))
    # compare the two spectra.
    fig, ax = plt.subplots(2, 1, sharex=True, figsize=(10,8))
    ax[0].plot(spec_bosz.wavelength, spec_bosz.flux)
    ax[0].plot(spec_marcs.wavelength, spec_marcs.flux, alpha=0.4)
    ax[0].set_xscale("log")
    subs = np.where((spec_bosz.wavelength>=spec_marcs.wavelength.min())&(spec_bosz.wavelength<=spec_marcs.wavelength.max()))[0]
    ax[1].plot(spec_bosz.loc[subs, "wavelength"], spec_bosz.loc[subs, "flux"]/spec_bosz.loc[subs, "flux_interp"])
    ax[1].axhline(1.0, lw=2.0, ls="--", c="r")
    ax[1].set_xscale("log")
    plt.show()

    # (4) compare flux densities and integrated emergent flux calculated with the MARCS and BOSZ grids.
    # BOSZ
    filename = bosz.filename.values[2]
    FluxDensity_bosz = calFluxDensity(grid_path, filename, filter_transmission_path, filters, grid="BOSZ")
    FluxDensity_bosz = np.float64(FluxDensity_bosz[1:])
    Fbol_bosz = calFbol(grid_path, filename, grid="BOSZ")
    # MARCS
    filename = marcs.filename.values[0]
    FluxDensity_marcs = calFluxDensity(marcs_spectrum_path, filename, filter_transmission_path, filters, grid="MARCS")
    FluxDensity_marcs = np.float64(FluxDensity_marcs[1:])
    Fbol_marcs = calFbol(marcs_spectrum_path, filename, grid="MARCS")
    # visualize.
    plt.scatter(FluxDensity_bosz, FluxDensity_marcs/FluxDensity_bosz)
    plt.plot([0,8], [1, 1])
    plt.show()





