import pandas as pd
import numpy as np 
from astroquery.gaia import Gaia
from dustmaps.bayestar import BayestarQuery
import astropy.units as units
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy.io import votable
import glob
import re 
import time
import os 
import phomKit
from inlist import *
import copy as cp



LgaiaPhot=True
Lgaiaedr3Distance=True
Lapass=True
Lsdss=True
Lhipparcos=True
Ltycho2=True
Lps1=True
L2mass=True
Lallwise=True
Lskymapper=True

########################################################################################################################################################
# (1) Input stellar parameters, namely Teff, logg, and [Fe/H] and their uncertainties.
upload_filename  = file_target_list
user_filename = "SEDEX"

# (1.1) Gaia photometry.
if LgaiaPhot:   
    # download Gaia photometry. 
    out_filename = Xpath+"Photometry_GAIA_DR3.csv"   
    query = ("SELECT targets.starID, gaia.* \
    FROM user_jyu01.{:s} AS targets \
    INNER JOIN gaiadr3.gaia_source AS gaia \
    ON gaia.source_id = targets.starID".format(user_filename))
    phomKit.access_crossmatch_from_Gaia(query, upload_filename=upload_filename, user_filename=user_filename, format="csv", out_filename=out_filename) 
    # calculate magnitude errors in the G, Bp, Rp bands
    result = pd.read_csv(out_filename, dtype={"starid":str, "source_id":str})
    gflux = result.loc[:, "phot_g_mean_flux"]
    gfluxerr = result.loc[:, "phot_g_mean_flux_error"]
    bpflux = result.loc[:, "phot_bp_mean_flux"]
    bpfluxerr = result.loc[:, "phot_bp_mean_flux_error"]
    rpflux = result.loc[:, "phot_rp_mean_flux"]
    rpfluxerr = result.loc[:, "phot_rp_mean_flux_error"]
    gmagerr, bpmagerr, rpmagerr = phomKit.Gaia_mag_errors(gflux=gflux, gfluxerr=gfluxerr, bpflux=bpflux, bpfluxerr=bpfluxerr, rpflux=rpflux, rpfluxerr=rpfluxerr)
    result.loc[:, "phot_g_mean_mag_error"] = gmagerr
    result.loc[:, "phot_bp_mean_mag_error"] = bpmagerr
    result.loc[:, "phot_rp_mean_mag_error"] = rpmagerr
    # select needed columns
    columns = ["starid", "source_id","phot_g_mean_mag","phot_g_mean_mag_error","phot_bp_mean_mag",\
        "phot_bp_mean_mag_error","phot_rp_mean_mag","phot_rp_mean_mag_error", "parallax", "parallax_error"]
    result = result[columns]
    result.to_csv(Xpath+"Photometry_GAIA_DR3_Subset"+".csv", index=False, float_format="%.4f")

# (1.2) download Gaia distance
if Lgaiaedr3Distance:
    out_filename = Xpath+"Distance_GAIA_DR3.csv"   
    query = ("SELECT targets.starID, dist.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN external.gaiaedr3_distance AS dist \
            ON dist.source_id = targets.starid ".format(user_filename))
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename) 
    # refine distance.
    dist = pd.read_csv(Xpath+"Distance_GAIA_DR3.csv", dtype={"starid":str, "source_id":str})
    dist.loc[:, "r_phogeo"] = dist.r_med_photogeo.copy()
    dist.loc[:, "r_phogeo_err"] =  (dist.r_hi_photogeo.copy()-dist.r_lo_photogeo.copy())/2
    dist.loc[:, "r_geo"] = dist.r_med_geo.copy()
    dist.loc[:, "r_geo_err"] =  (dist.r_hi_geo.copy()-dist.r_lo_geo.copy())/2
    dist.loc[:,"d"] = dist.loc[:, "r_phogeo"].copy()
    dist.loc[:, "derr"] = dist.loc[:, "r_phogeo_err"].copy()
    subs = np.where(pd.isna(dist.loc[:, "r_phogeo"]))[0]
    dist.loc[subs, "d"] = dist.loc[subs, "r_geo"].copy()
    dist.loc[subs, "derr"] = dist.loc[subs, "r_geo_err"].copy() 
    dist = dist[["starid", "source_id", "d", "derr"]].copy()
    dist.to_csv(Xpath+"Distance_GAIA_DR3_Subset.csv", index=False, float_format=("%.4f"))

# (1.3) access APASS crossmatch from Gaia archive.
if Lapass:
    out_filename = Xpath+"Photometry_APASS_DR9.csv"   
    query = ("SELECT targets.starID, apass.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.apassdr9_best_neighbour AS apass \
            ON apass.source_id = targets.starID \
            INNER JOIN external.apassdr9 AS catalog \
            ON catalog.recno=apass.original_ext_source_id ".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

# (1.4) access SDSS crossmatch from Gaia archive.
if Lsdss:
    out_filename = Xpath+"Photometry_SDSS_DR13.csv" 
    query = ("SELECT targets.starID, sdss.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.sdssdr13_best_neighbour AS sdss \
            ON sdss.source_id = targets.starID \
            INNER JOIN external.sdssdr13_photoprimary AS catalog \
            ON catalog.objid=sdss.original_ext_source_id ".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)
            
# (1.5) access Hipparcos crossmatch from Gaia archive.
if Lhipparcos:
    out_filename = Xpath+"Photometry_Hipparcos2.csv"     
    query = ("SELECT targets.starID, hp.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.hipparcos2_best_neighbour AS hp \
            ON hp.source_id = targets.starID \
            INNER JOIN public.hipparcos_newreduction AS catalog \
            ON catalog.hip=hp.original_ext_source_id ".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

# (1.6) access Tycho2 crossmatch from Gaia archive.
if Ltycho2:
    out_filename = Xpath+"Photometry_TYCHO2.csv"     
    query = ("SELECT targets.starID, tycho2.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.tycho2tdsc_merge_best_neighbour AS tycho2 \
            ON tycho2.source_id = targets.starID \
            INNER JOIN public.tycho2 AS catalog \
            ON catalog.id=tycho2.original_ext_source_id ".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

##########################################################################################################
# (1.7) access PAN-STARRS DR1 crossmatch from Gaia archive.
if Lps1:
    out_filename = Xpath+"Photometry_PS1.csv"     
    query = ("SELECT targets.starID, ps1.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.panstarrs1_best_neighbour AS ps1 \
            ON ps1.source_id = targets.starID \
            INNER JOIN gaiadr2.panstarrs1_original_valid AS catalog \
            on catalog.obj_id=ps1.original_ext_source_id".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

##########################################################################################################
# (1.8) access 2MASS crossmatch from Gaia archive.
if L2mass:
    out_filename = Xpath+"Photometry_2MASS.csv"     
    query = ("SELECT targets.starID, twomass.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS twomass \
            ON twomass.source_id = targets.starID \
            INNER JOIN gaiadr1.tmass_original_valid as catalog \
            on catalog.designation=twomass.original_ext_source_id".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

##########################################################################################################
# (1.9) access ALLWISE crossmatch from Gaia archive.
if Lallwise:
    out_filename = Xpath+"Photometry_ALLWISE.csv"     
    query = ("SELECT targets.starID, allwise.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.allwise_best_neighbour AS allwise \
            ON allwise.source_id = targets.starID \
            INNER JOIN gaiadr1.allwise_original_valid as catalog \
            ON catalog.designation=allwise.original_ext_source_id ".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)

##########################################################################################################
# (1.10) access SKymapper DR2 crossmatch from Gaia archive.
if Lskymapper:
    out_filename = Xpath+"Photometry_SKYMAPPER_DR2.csv"     
    query = ("SELECT targets.starID, skymapper.*, catalog.* \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.skymapperdr2_best_neighbour AS skymapper \
            ON skymapper.source_id = targets.starID \
            INNER JOIN external.skymapperdr2_master as catalog \
            ON catalog.object_id=skymapper.original_ext_source_id".format(user_filename))   
    phomKit.access_crossmatch_from_Gaia(query, user_filename=user_filename, out_filename=out_filename)






########################################################################################################################################################
# (2) combine photometry from all of the catalogues.
# Gaia DR3 photometry.
columns = ["starid","phot_g_mean_mag","phot_g_mean_mag_error","phot_bp_mean_mag","phot_bp_mean_mag_error","phot_rp_mean_mag","phot_rp_mean_mag_error"]
gaia = pd.read_csv(Xpath+"Photometry_GAIA_DR3_Subset.csv", usecols=columns, dtype={"starid":str})
gaia = gaia.rename(columns={
    "phot_g_mean_mag":"gaiadr3g",
    "phot_g_mean_mag_error":"gaiadr3g_e",
    "phot_bp_mean_mag":"gaiadr3bp",
    "phot_bp_mean_mag_error":"gaiadr3bp_e",
    "phot_rp_mean_mag":"gaiadr3rp",    
    "phot_rp_mean_mag_error":"gaiadr3rp_e"})

# Hipparcos2
columns = ['starid', 'hp_mag', 'e_hp_mag']
hipparcos2 = pd.read_csv(Xpath+"Photometry_Hipparcos2.csv", usecols=columns, dtype={"starid":str})
hipparcos2 = hipparcos2.rename(columns={
    "hp_mag": "hipparcos2hp",
    'e_hp_mag': "hipparcos2hp_e"})



# Tycho2
columns = ['starid', 'bt_mag', 'e_bt_mag', 'vt_mag', 'e_vt_mag']
tycho2 = pd.read_csv(Xpath+"Photometry_TYCHO2.csv", usecols=columns, dtype={"starid":str})
tycho2 = tycho2.rename(columns={
    "bt_mag": "tycho2bt",
    'e_bt_mag': "tycho2bt_e",
    'vt_mag': "tycho2vt",
    'e_vt_mag': "tycho2vt_e"})

# APASS
columns = ["starid", 'vmag','e_vmag','bmag','e_bmag','g_mag','e_g_mag','r_mag','e_r_mag','i_mag','e_i_mag']
apass9 = pd.read_csv(Xpath+"Photometry_APASS_DR9.csv", usecols=columns, dtype={"starid":str})
apass9 = apass9.rename(columns={
    "vmag":"apassv",
    "e_vmag":"apassv_e",
    "bmag":"apassb",    
    "e_bmag":"apassb_e",
    "g_mag":"apassg",    
    "e_g_mag":"apassg_e",  
    "r_mag":"apassr",
    "e_r_mag":"apassr_e",
    "i_mag":"apassi",
    "e_i_mag":"apassi_e"})

# SDSS photometry
columns = ['starid','petromag_u','petromagerr_u','petromag_g','petromagerr_g','petromag_r','petromagerr_r','petromag_i','petromagerr_i','petromag_z','petromagerr_z']
sdss = pd.read_csv(Xpath+"Photometry_SDSS_DR13.csv", usecols=columns, dtype={"starid":str})
sdss = sdss.rename(columns={
    "petromag_u":"sdssu",
    "petromagerr_u":"sdssu_e",
    "petromag_g":"sdssg",
    "petromagerr_g":"sdssg_e",
    "petromag_r":"sdssr",
    "petromagerr_r":"sdssr_e",
    "petromag_i":"sdssi",
    "petromagerr_i":"sdssi_e",
    "petromag_z":"sdssz",
    "petromagerr_z":"sdssz_e"})

# Pan-Starrs photometry.
columns = ['starid', 'g_mean_psf_mag', 'g_mean_psf_mag_error',
        'r_mean_psf_mag', 'r_mean_psf_mag_error','i_mean_psf_mag', 'i_mean_psf_mag_error', 
        'z_mean_psf_mag', 'z_mean_psf_mag_error','y_mean_psf_mag', 'y_mean_psf_mag_error']
ps1 = pd.read_csv(Xpath+"Photometry_PS1.csv", usecols=columns, dtype={"starid":str})
ps1 = ps1.rename(columns={
    "g_mean_psf_mag":"panstarrsps1g",
    "g_mean_psf_mag_error":"panstarrsps1g_e",
    "r_mean_psf_mag":"panstarrsps1r",  
    "r_mean_psf_mag_error":"panstarrsps1r_e",
    "i_mean_psf_mag":"panstarrsps1i",   
    "i_mean_psf_mag_error":"panstarrsps1i_e",
    "z_mean_psf_mag":"panstarrsps1z",  
    "z_mean_psf_mag_error":"panstarrsps1z_e",
    "y_mean_psf_mag":"panstarrsps1y",         
    "y_mean_psf_mag_error":"panstarrsps1y_e"})


# SkyMapper photometry.
columns = ["starid", "u_psf", "e_u_psf", "v_psf", "e_v_psf", "g_psf", "e_g_psf",
                    "r_psf", "e_r_psf", "i_psf", "e_i_psf", "z_psf", "e_z_psf"]
skymapper = pd.read_csv(Xpath+"Photometry_SKYMAPPER_DR2.csv", usecols=columns, dtype={"starid":str})
skymapper = skymapper.rename(columns={
    "u_psf":"skymapperu",
    "e_u_psf":"skymapperu_e",
    "v_psf":"skymapperv",  
    "e_v_psf":"skymapperv_e",
    "g_psf":"skymapperg",   
    "e_g_psf":"skymapperg_e",
    "r_psf":"skymapperr",  
    "e_r_psf":"skymapperr_e",
    "i_psf":"skymapperi",         
    "e_i_psf":"skymapperi_e",
    "z_psf":"skymapperz",         
    "e_z_psf":"skymapperz_e"})                       


# ALLWISE photometry.
columns = ["starid", 'w1mpro', 'w1mpro_error', 'w2mpro', 'w2mpro_error', 'w3mpro', 'w3mpro_error', 'w4mpro', 'w4mpro_error']
allwise = pd.read_csv(Xpath+"Photometry_ALLWISE.csv", usecols=columns, dtype={"starid":str})
allwise = allwise.rename(columns={
    "w1mpro":"wisew1",
    "w1mpro_error":"wisew1_e",
    "w2mpro":"wisew2",
    "w2mpro_error":"wisew2_e",
    "w3mpro":"wisew3",
    "w3mpro_error":"wisew3_e",
    "w4mpro":"wisew4",
    "w4mpro_error":"wisew4_e"})

# 2MASS photometry.
columns = ["starid", 'j_m', 'j_msigcom', 'h_m', 'h_msigcom', 'ks_m', 'ks_msigcom']
twomass = pd.read_csv(Xpath+"Photometry_2MASS.csv", usecols=columns, dtype={"starid":str})
twomass = twomass.rename(columns={
    "j_m":"twomassj", 
    "j_msigcom":"twomassj_e",
    "h_m":"twomassh",
    "h_msigcom":"twomassh_e",   
    "ks_m":"twomassk", 
    "ks_msigcom":"twomassk_e"})


# combine photometry results.
photo = pd.merge(gaia, tycho2, on="starid", how="left")
photo = pd.merge(photo, hipparcos2, on="starid", how="left")
photo = pd.merge(photo, apass9, on="starid", how="left")
photo = pd.merge(photo, sdss, on="starid", how="left")
photo = pd.merge(photo, allwise, on="starid", how="left")
photo = pd.merge(photo, twomass, on="starid", how="left")
photo = pd.merge(photo, ps1, on="starid", how="left")
photo = pd.merge(photo, skymapper, on="starid", how="left")
photo = photo.rename(columns={
"starid":"starID",    
"gaiadr3g":"gaiag", 
"gaiadr3g_e":"gaiag_e", 
"gaiadr3bp":"gaiabp", 
"gaiadr3bp_e":"gaiabp_e", 
"gaiadr3rp": "gaiarp", 
"gaiadr3rp_e":"gaiarp_e"})


# prepare the input photometry table for implementing the SED fitting. 
tycho2_duplicatesID = tycho2[tycho2[["starid"]].duplicated()].starid.values
photo = photo[~np.isin(photo.starID, tycho2_duplicatesID)].reset_index(drop=True)
# replace negative magnitudes and errors with NaN.
photo.iloc[:, 1:] = photo.iloc[:, 1:].where(photo.iloc[:, 1:]>=0, np.nan)
# change nan from object to np.nan
photo.iloc[:, 1:] = photo.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
# inflat errors.
photo = phomKit.InfltError(InputPhotoMetry=photo.copy())
# Keep photo that have at least 5 bands photometry.
photo = photo[((len(photo.columns)-1)-photo.isnull().sum(axis=1))/2>4].reset_index(drop=True)






##########################################################################################################
# (3) combine stellar parameters and photometry.
params = pd.read_csv(file_target_list, dtype={"starID":str})
params = params[(~params.teff.isna()) & (~params.logg.isna())].reset_index(drop=True)
dist = pd.read_csv(Xpath+"Distance_GAIA_DR3_Subset.csv", usecols=["starid", "d", "derr"], dtype={"starid":str})
dist = dist.rename(columns={"starid":"starID"})
params = pd.merge(params, dist, on="starID").reset_index(drop=True)
# make parameter and photometry dataframes consistent.
params_columns = cp.copy(params.columns)
photom_columns = cp.copy(photo.columns)
table = pd.merge(params, photo, on="starID").reset_index(drop=True)
params = table[params_columns]
photo = table[photom_columns]
params.to_csv(file_input_params, index=False, float_format="%.6f")
photo.to_csv(file_input_photometry, index=False, float_format="%.4f")
