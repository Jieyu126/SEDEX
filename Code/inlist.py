import os
import numpy as np

####################################################################################################################################### 
Lbatch = False           # run code in the batch mode or in the single-star mode, the latter of which is helpful for debugging.
Lshow = False            # show figures interactively, stwich this off when run the pipeline in the batch mode.
Lsave = True             # whether to save fitting results, no mather in the batch mode or single star mode. Useful to unset it for debugging..
index_start = 0          # index starting from 0. Only works when Lbatch=False
index_end = 29           # if you have 100 stars, for instance, then index_end should be 100. Only works when Lbatch=False
sample = "TOIs"          # if one wants to run different samples individually: "TOIs", "Seismic", "APOGEE", "GALAH", "RAVE"

####################################################################################################################################### 
file_target_list = "../Data/Input_Fits/"+sample+"/UserInputData.csv"    # taret list
Xpath = "../Data/Input_Fits/"+sample+"/CrossMatchTables/"               # to save crossmatch tables fetched from the Gaia Archive
bandclip = np.array(["wisew3", "wisew4"])                               # exclude certain bands for fitting


####################################################################################################################################### 
# Select a extinction law
ExtinctionLaw = "F19"    # CCM89, O94, F99, F04, M14, G16, F19. For more detains see: https://dust-extinction.readthedocs.io/en/stable/index.html
Rv = 3.1                 # R(V) = A(V)/E(B-V)
grid = "MARCS"           # MARCS or BOSZ. BOSZ for Teff>8000K.
LvaryFeh = True          # Use known metallicity rather than assume solar metallicity. Techiniqly, this is [M/H] instead of [Fe/H]
LFixedAv = False         # use known Av
LAvprior = False         # use known Av as a prior.

####################################################################################################################################### 
# define some paths for saving results from Blackbody fitting and SED fitting.
fontsize = 15                                     # fontsize of the main figure 
figext = "png"                                    # extension of the filename of figures: "png", "pdf"
L_fig_h5py = False                                # whether to combine all the figures generated in a single batch mode. Works when Lbatch is set.
grid_list_path = "../Data/MetaData/"
filter_zero_point_path = "../Data/MetaData/filters.dat"
filter_transmission_path = "../Data/MetaData/SpectralWindows/"
if grid == "MARCS": grid_spectrum_path = "/home/yujie/Documents/SoftLinks/Scratch/SpectrumGrid/MARCS/Version3/Models/"
if grid == "BOSZ": grid_spectrum_path = "/home/yujie/Documents/SoftLinks/Scratch/SpectrumGrid/BOSZ/Grid/"
file_input_params = "../Data/Input_Fits/"+sample+"/Input_params.csv"
file_input_photometry = "../Data/Input_Fits/"+sample+"/Input_Photometry.csv"

######################################################################################################################################
# For 4_BlackBodyFit.py
Optimize_Av = False                           # if True, optimize Av
LOutput_preditcted_flux = False               # If True, save best-fitting fluxes to a file
LOutput_optimized_Teff_Av_ScaleFactor = True  # If True, save optimised Teff and Av to a file
Teffmin = 1000                                # mim model Teff
Teffmax = 30000                               # max model Teff
save_blackbody_flux_path_single = "../Data/Output_Fits/"+sample+"/BlackBodyFits/SingleMode/"
save_blackbody_flux_path_batch = "../Data/Output_Fits/"+sample+"/BlackBodyFits/BatchMode/"
save_blackbody_plot_path = "../Figures/"+sample+"/BlackBodyFits/"
Combined_BlackBody_Fits_path = "../Data/Output_Fits/"+sample+"/BlackBodyFits/"


######################################################################################################################################
# mainly for 6_SedFit.py
# (1) Meta parameters to control the code.
LPhotometryDataPreapred = True  # use prepared data.
Lprint = False                  # print detailed information results for debugging.
Lblackbody = True               # use blackbody fits to exclude bad data points and to set initial guesses for for the angular radius.
Linterpolate = True             # if true, interpolate marcs model fluxes in a grid space of teff and logg.
# (2) for interpolation
tefferr_sys = 0                 # add in quadruture to the teff formal error
loggerr_sys = 0                 # add in quadruture to the logg formal error
averr_sys = 0                   # add in quadruture to the av formal error
teffStep = 5                    # teff step size of the teff-logg grid for model flux interpolation
loggStep = 0.25                 # logg step size of the teff-logg grid for model flux interpolation
# (3) pick data points used for the fitting.
Lpeakpoints = True              # using given number of data points for the fitting.
minnpoint_ori = 4               # min number of data points used to optimize scaling factor. This flag works even when the data-point droppin method is used. 
maxnpoint_ori = 50              # max number of data points used to optimize scaling factor. if maxnpoint = -1, all data points used.
LFluxDropRatio = False          # select the datapoints for fitting based on flux, flux<max(flux)*FluxDropRatio
FluxDropRatio = 0.01            # using data points if their fluxes are not less than FluxDropRatio times the maximum flux.
FluxMaxDiff = 0.10              # drop observed photometric data points if the difference between the observed and predicted is larger than this number.
# (4)  Picking models.  
Lteffdiff = False               # select models use teff diff.
Teffdiff = 300                  # Search MARCS models whose Teff meets the condition |Teff-Teff,model|<Teffdiff
MaxTeff = 10000                 # select models with Teff<=MaxTeff when Teff is not available.
Nteff = 1                       # select the nearest model given teff and logg   
Nlogg = 1                       # select the nearest model given teff and logg   
Maxlog = 6.5                    # select models with logg<=Maxlog when logg is not available.
extrap_wave_max = 30            # maximum wavelength for extrapolation, only for the MARCS grid.
LweightedMean = False           # use posterior-weighted mean to approximate logg, Teff, Av, ScaleF, and fbol. Can be biased if pdf is truncated.   
LGassianFit = True              # Fit Gaussians to posteriors to estimate Teff, Av, ScaleF, and fbol. This works properly when posteriors are truncated, 
                                # which happens when input parameters are near the edge of spectral models, or are significantly off the ground truth. 
# (5) paths for input and output data and figures.
fig_sedfit_path = "../Figures/"+sample+"/SEDFits/"+grid+"/"
single_data_sedfit_path = "../Data/Output_Fits/"+sample+"/SEDFits/"+grid+"/SingleMode/"
batch_data_sedfit_path =  "../Data/Output_Fits/"+sample+"/SEDFits/"+grid+"/BatchMode/"
combined_data_sedfit_path = "../Data/Output_Fits/"+sample+"/SEDFits/"+grid+"/"


# creat some folders.
if __name__ == '__main__':
    # Execute when the module is not initialized from an import statement.
    # to save downloaded input photometry
    if not os.path.isdir(Xpath): os.system("mkdir -p " + Xpath) 
    # to save results returned from the blackbody fitting module.
    if not os.path.isdir(save_blackbody_flux_path_single): os.system("mkdir -p " + save_blackbody_flux_path_single)
    if not os.path.isdir(save_blackbody_flux_path_batch): os.system("mkdir -p " + save_blackbody_flux_path_batch)
    if not os.path.isdir(save_blackbody_plot_path): os.system("mkdir -p " + save_blackbody_plot_path)
    if not os.path.isdir(Combined_BlackBody_Fits_path): os.system("mkdir -p " + Combined_BlackBody_Fits_path)
    # to save results returned from the SED fitting module
    if not os.path.isdir(fig_sedfit_path): os.system("mkdir -p " + fig_sedfit_path)
    if not os.path.isdir(single_data_sedfit_path): os.system("mkdir -p " + single_data_sedfit_path)
    if not os.path.isdir(batch_data_sedfit_path): os.system("mkdir -p " + batch_data_sedfit_path)
    if not os.path.isdir(combined_data_sedfit_path): os.system("mkdir -p " + combined_data_sedfit_path)
