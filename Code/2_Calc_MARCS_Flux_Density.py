import pandas as pd 
import glob 
import numpy as np 
from SedFitFunctions import *
from inlist import *
from multiprocessing.dummy import Pool as ThreadPool


def extract_params(chemcomp="st"):
    # each filename has the exactly same length.
    models = pd.read_csv(marcs_list_path+"AllModels.list", names=["filename"])
    # select particular models
    if chemcomp=="st": models = models[models.filename.str.contains("st")]
    # if Lspheric==True: models = models[~models.filename.str.contains("p")]
    # if Lparallel==True: models = models[models.filename.str.contains("p")]
    teff = np.array(list(map(lambda x: x[1:5], models.filename)))
    logg = np.array(list(map(lambda x: x[7:11], models.filename)))
    feh = np.array(list(map(lambda x: x[25:30], models.filename)))
    mass = np.array(list(map(lambda x: x[13:16], models.filename)))
    geom = np.array(list(map(lambda x: x[0:1], models.filename)))
    tv = np.array(list(map(lambda x: x[18:20], models.filename)))
    chem = np.array(list(map(lambda x: x[21:23], models.filename)))
    # save for the final results.
    data = pd.DataFrame({"filename":models.filename,"modelindex":models.index,
                        "teff":teff, "logg":logg, "feh":feh, 
                        "mass":mass, "geom":geom, "tv":tv, "chem":chem})
    # note that the following modification was for making the values exactly equal to 
    # the input required for downloading MARCS model atmosphere.                     
    # remove +, such as for +2.0.
    data.loc[:, "logg"].replace(regex=True, inplace=True, to_replace="\+", value="")
    data.loc[:, "feh"].replace(regex=True, inplace=True, to_replace="\+", value="")
    # remove trailing .0, such as 2.0.
    data.loc[:, "logg"].replace(regex=True, inplace=True, to_replace="\.0$", value="")
    # remove trailing .0, such as 2.00.
    data.loc[:, "feh"].replace(regex=True, inplace=True, to_replace="\.00$", value="")
    # remove trailing .0, such as 2.20, 0.00.
    data.loc[:, "feh"].replace(regex=True, inplace=True, to_replace="0$", value="")
    # remove trailing ., such as for mass, 15.
    data.loc[:, "mass"].replace(regex=True, inplace=True, to_replace="\.$", value="")
    # recover 0
    data.loc[:, "feh"].replace(regex=True, inplace=True, to_replace="", value="0")
    data.loc[:, "mass"].replace(regex=True, inplace=True, to_replace="\.0$", value="")
    # remove heading 0 for turbulence velocity.
    data.loc[:, "tv"].replace(regex=True, inplace=True, to_replace="^0", value="")
    # output
    # data.to_csv("../Data/MetaData/AllModelsParams.csv",index=False)
    return data.reset_index(drop=True).copy()


def loop_calFluxDensity(ind):
    print("#{:d} calculate flux densityeis for the model: {:s}".format(ind, Fluxes.loc[ind, "filename"]))
    # read in the wavelength and flux of MARCS models.
    spec_ori = readspec(marcs_spectrum_path, Fluxes.loc[ind, "filename"])
    # convolve the model atmosphere with the filter transmission for each band chosen. The model spectrum has flux in 
    # units of erg/cm^2/s/Anstrogm. The resulting flux density for each band chosen are later used for minimisation.
    return Fluxes.loc[ind, "filename"], np.log10(FluxDensity(spec_ori, bandpass=filters.band.values.copy(), path=filter_transmission_path))


def loop_calFbol(ind):
    print("#{:d} calculate bolemetric fluxes for the model: {:s}".format(ind, models.loc[ind, "filename"]))
    # read in the wavelength and flux of MARCS models.
    spec_ori = readspec(marcs_spectrum_path, models.loc[ind, "filename"])
    # extrapolate spectra from 20 to 30 micron
    spec_ori = flux_extrap(spec_ori.copy(), Minwavelength=8, Maxwavelength=30, order=3)
    # calculate bolometric flux on the stellar surface, in log10 scale.
    Fbol = np.log10(np.trapz(spec_ori.flux, x=spec_ori.wavelength*1e4))
    return models.loc[ind, "filename"], Fbol


# The Spherical-geometry models come with 4 different masses: 0.5, 1.0, 2.0, and 5.0 MSun. 
# Only the 1.0 MSun subgrids have all the combinations of Teff and log g. The plane-parallel 
# models are independent of mass, indicated by M=0.0 MSun, so we manually re-label it to M=1.0 Msun.  

##########################################################################################################
# load the model filenames and extract stellar parameters from the filenames.
models = extract_params(chemcomp="st")
models.loc[:, "filename"] += ".flx"
model_avail = pd.DataFrame({"filename": map(os.path.basename, np.sort(glob.glob(marcs_spectrum_path+"*st*.flx")))})
# Check if the models in the model list are available.
for i in range(len(models)):
    if not (models.loc[i, "filename"] in model_avail.filename.values): 
        print("Note that this model name is in the model list, but the spectrum is not available at MARCS webpage.")
        print(models.loc[i, "filename"])
models = pd.merge(models, model_avail, on="filename").reset_index(drop=True).copy() 
# calculate bolometric flux on the stellar surface.
pool = ThreadPool(10)
t0 = time.time()
results = pool.map(loop_calFbol, range(len(models)))
t1 = time.time()
print("Calculate flux densities using #{:} MARCS models, {:.2f} seconds used.".format(len(models), t1-t0))
Fbol = pd.DataFrame(columns=["filename", "fbol"], data=results)
models = pd.merge(models, Fbol, on="filename")
# save model parameters.
models[["teff","logg","feh","tv","mass","modelindex", "fbol"]].to_csv(marcs_list_path+"StandardComposition_SphericalPlane.txt", index=False)


##########################################################################################################
# load the passbands against which model fluxes will be evaluated.
filters = pd.read_csv("../Data/MetaData/filters.dat")
Fluxes = pd.DataFrame(index=range(len(models)), columns=["filename"]+filters.band.values.copy().tolist())
Fluxes.loc[:, "filename"] = models.filename.copy()

# parallel the loop.
pool = ThreadPool(10)
t0 = time.time()
results = pool.map(loop_calFluxDensity, range(len(Fluxes)))
t1 = time.time()
print("Calculate flux densities using #{:} MARCS models, {:.2f} seconds used.".format(len(Fluxes), t1-t0))
for i in range(len(results)): Fluxes.iloc[i, 0], Fluxes.iloc[i, 1:] = results[i][0], results[i][1]
for i in Fluxes.columns[1:]: Fluxes.loc[:, i] = pd.to_numeric(Fluxes.loc[:, i], downcast='float')
Fluxes = pd.merge(models, Fluxes, on="filename").copy()
Fluxes.to_csv(marcs_list_path+"PassbandFlux_MARCS.csv",index=False, float_format="%.8f")

