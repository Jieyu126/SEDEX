import pandas as pd 
import os 
import numpy as np 
from inlist import * 
import matplotlib.pyplot as plt
import matplotlib.image as mpimg



starID = input("Input the target ID for visual inspection: ")

# get the index of the target star.
params = pd.read_csv(file_input_params, dtype={"starID":str})
photo = pd.read_csv(file_input_photometry, dtype={"starID":str})

# check the index
subs = np.where(photo.starID==str(starID))[0]
print("")
print("(0) ===================================================================================================")
print("index of this target in the whole target list: {:d}".format(subs[0]))


# load filter list.
filters = pd.read_csv("../Data/MetaData/filters.dat")
print("")
print("(1) ===================================================================================================")
print("All Filters collected for fitting:")
print(filters)

# load photometry
pho_ori = pd.read_csv(file_input_photometry, dtype={"starID":str})
pho_ori = pho_ori[pho_ori.starID==str(starID)].reset_index(drop=True)
photometry = pd.DataFrame(columns=["band", "mag", "magerr"])
photometry["band"] = pho_ori.iloc[0, 1::2].index
photometry["mag"] = pho_ori.iloc[0, 1::2].values
photometry["magerr"] = pho_ori.iloc[0, 2::2].values
photometry[["mag", "magerr"]] = photometry[["mag", "magerr"]].apply(pd.to_numeric, errors='coerce', axis=1)
photometry = pd.merge(photometry, filters, on="band").sort_values(by=["lambda"])[["band", "mag", "magerr"]]
print("\n\n\n")
print("(2) ===================================================================================================")
print("original photometry:")
print(photometry)

# load Black-body fitting output
pho_blkbd = pd.read_csv(Combined_BlackBody_Fits_path+"Output_BlackBody_Fits.csv", dtype={"starID":str})
print("\n\n\n")
print("(3) ===================================================================================================")
print("Black-body fitting output:")
print(pho_blkbd[pho_blkbd.starID==str(starID)])

# # load filters used for sed fitting.
# sed_filters = pd.read_csv(save_sedfit_flux_path + starID + ".flux.csv", dtype={"starID":str})
# print("\n\n\n")
# print("(4) ===================================================================================================")
# print("filters SED fitting used")
# print(sed_filters)

# # load best fitting result

# sedmodel = pd.read_csv(save_sedfit_flux_path + starID +".BestfitModels.csv", dtype={"starID":str})
# print("\n\n\n")
# print("(5) ===================================================================================================")
# print("SED fitting output:")
# print(sedmodel)


# load SED fit
fig, ax = plt.subplots(1,1)
img=mpimg.imread(fig_sedfit_path + starID +".png")
ax.imshow(img)


# load blackbody fit
fig, ax = plt.subplots(1,1)
img=mpimg.imread(save_blackbody_plot_path + starID +".png")
ax.imshow(img)
plt.show()


