import pandas as pd

# reformat the 2MASS filter transmission data
twomassj = pd.read_csv("original_files/2massj.dat", sep="\s+", header=None,names=["wavelength", "response"])
twomassj.to_csv("twomassj.dat",index=False, float_format="%.8f")
twomassh = pd.read_csv("original_files/2massh.dat", sep="\s+", header=None,names=["wavelength", "response"])
twomassh.to_csv("twomassh.dat",index=False, float_format="%.8f")
twomassk = pd.read_csv("original_files/2massk.dat", sep="\s+", header=None,names=["wavelength", "response"])
twomassk.to_csv("twomassk.dat",index=False, float_format="%.8f")

# reformat the WISE filter transmission data
# wisew1 = pd.read_csv("original_files/wisew1.dat", sep="\s+", usecols=(0,1))
# wisew1.to_csv("wisew1.dat", index=False, float_format="%.8f")
# wisew2 = pd.read_csv("original_files/wisew2.dat", sep="\s+", usecols=(0,1))
# wisew2.to_csv("wisew2.dat", index=False, float_format="%.8f")
# wisew3 = pd.read_csv("original_files/wisew3.dat", sep="\s+", usecols=(0,1))
# wisew3.to_csv("wisew3.dat", index=False, float_format="%.8f")
# wisew4 = pd.read_csv("original_files/wisew4.dat", sep="\s+", usecols=(0,1))
# wisew4.to_csv("wisew4.dat", index=False, float_format="%.8f")
wisew1 = pd.read_csv("original_files/wisew1.dat", skiprows=2, header=None, names=["wavelength", "response"], sep="\s+", usecols=(0,1))
wisew1.to_csv("wisew1.dat", index=False, float_format="%.8f")
wisew2 = pd.read_csv("original_files/wisew2.dat", skiprows=2, header=None, names=["wavelength", "response"], sep="\s+", usecols=(0,1))
wisew2.to_csv("wisew2.dat", index=False, float_format="%.8f")
wisew3 = pd.read_csv("original_files/wisew3.dat", skiprows=2, header=None, names=["wavelength", "response"], sep="\s+", usecols=(0,1))
wisew3.to_csv("wisew3.dat", index=False, float_format="%.8f")
wisew4 = pd.read_csv("original_files/wisew4.dat", skiprows=2, header=None, names=["wavelength", "response"], sep="\s+", usecols=(0,1))
wisew4.to_csv("wisew4.dat", index=False, float_format="%.8f")


# reformat the SDSS filter transmission data
sdssu = pd.read_csv("original_files/sdssu.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
sdssu["wavelength"] =  sdssu.wavelength/10000
sdssu.to_csv("sdssu.dat", index=False,float_format="%.8f")
sdssg = pd.read_csv("original_files/sdssg.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
sdssg["wavelength"] =  sdssg.wavelength/10000
sdssg.to_csv("sdssg.dat", index=False,float_format="%.8f")
sdssr = pd.read_csv("original_files/sdssr.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
sdssr["wavelength"] =  sdssr.wavelength/10000
sdssr.to_csv("sdssr.dat", index=False,float_format="%.8f")
sdssi = pd.read_csv("original_files/sdssi.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
sdssi["wavelength"] =  sdssi.wavelength/10000
sdssi.to_csv("sdssi.dat", index=False,float_format="%.8f")
sdssz = pd.read_csv("original_files/sdssz.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
sdssz["wavelength"] =  sdssz.wavelength/10000
sdssz.to_csv("sdssz.dat", index=False,float_format="%.8f")

# reformat the Gaia filter transmission data for Gaia DR2 revised.
# Gaia spectral response data
# gaia = pd.read_csv("original_files/GaiaDR2_RevisedPassbands.dat", sep="\s+")
# gaiag = gaia[["wavelength", "Gresponse"]]
# gaiabp = gaia[["wavelength", "BPresponse"]]
# gaiarp = gaia[["wavelength", "RPresponse"]]
# gaiag = gaiag.rename(columns={"Gresponse":"response"})
# gaiabp = gaiabp.rename(columns={"BPresponse":"response"})
# gaiarp = gaiarp.rename(columns={"RPresponse":"response"})
# gaiag["wavelength"]  = gaiag.wavelength /1000
# gaiabp["wavelength"] = gaiabp.wavelength/1000
# gaiarp["wavelength"] = gaiarp.wavelength/1000
# gaiag = gaiag[gaiag.response!=99.99].reset_index(drop=True)
# gaiabp = gaiabp[gaiabp.response!=99.99].reset_index(drop=True)
# gaiarp = gaiarp[gaiarp.response!=99.99].reset_index(drop=True)
# gaiag.to_csv("gaiag.dat",index=False, float_format="%.8f")
# gaiabp.to_csv("gaiabp.dat",index=False, float_format="%.8f")
# gaiarp.to_csv("gaiarp.dat",index=False, float_format="%.8f")
# reformat the Gaia filter transmission data for Gaia EDR3 revised.
# gaiag = pd.read_csv("original_files/GAIA_GAIA3.G.dat", sep="\s+")
# gaiabp = pd.read_csv("original_files/GAIA_GAIA3.Gbp.dat", sep="\s+")
# gaiarp = pd.read_csv("original_files/GAIA_GAIA3.Grp.dat", sep="\s+")
# gaiag["wavelength"]  = gaiag.wavelength /10000
# gaiabp["wavelength"] = gaiabp.wavelength/10000
# gaiarp["wavelength"] = gaiarp.wavelength/10000
# gaiag.to_csv("gaiag.dat",index=False, float_format="%.8f")
# gaiabp.to_csv("gaiabp.dat",index=False, float_format="%.8f")
# gaiarp.to_csv("gaiarp.dat",index=False, float_format="%.8f")
# filter transmission curves from Gaia website: https://www.cosmos.esa.int/web/gaia/edr3-passbands
columns = ["wavelength", "g", 'bp', "rp"]
gaia = pd.read_csv("original_files/GAIAEDR3.dat", sep="\s+", header=None, names=columns, usecols=(0,1,3,5))
gaiag = gaia[["wavelength", "g"]].copy()
gaiag = gaiag.rename(columns={"g":"response"})
gaiag = gaiag[gaiag.response<99].reset_index(drop=True)
gaiag["wavelength"]  = gaiag.wavelength /1000
gaiag.to_csv("gaiag.dat",index=False, float_format="%.8f")
gaiabp = gaia[["wavelength", "bp"]].copy()
gaiabp = gaiabp.rename(columns={"bp":"response"})
gaiabp = gaiabp[gaiabp.response<99].reset_index(drop=True)
gaiabp["wavelength"] = gaiabp.wavelength/1000
gaiabp.to_csv("gaiabp.dat",index=False, float_format="%.8f")
gaiarp = gaia[["wavelength", "rp"]].copy()
gaiarp = gaiarp.rename(columns={"rp":"response"})
gaiarp = gaiarp[gaiarp.response<99].reset_index(drop=True)
gaiarp["wavelength"] = gaiarp.wavelength/1000
gaiarp.to_csv("gaiarp.dat",index=False, float_format="%.8f")

#Tycho2
# from Mann & von Braun 2015
tycho2bt = pd.read_csv("original_files/TYCHO_BT_Mann2015.dat", sep="\s+")
tycho2vt = pd.read_csv("original_files/TYCHO_VT_Mann2015.dat", sep="\s+")
tycho2bt["wavelength"] = tycho2bt.wavelength /10000
tycho2vt["wavelength"] = tycho2vt.wavelength/10000
tycho2bt.to_csv("tycho2bt.dat",index=False, float_format="%.8f")
tycho2vt.to_csv("tycho2vt.dat",index=False, float_format="%.8f")

# Hipparcos
# from Mann & von Braun 2015
hipparcoshp = pd.read_csv("original_files/Hipparcos_Mann2015.dat", sep="\s+")
hipparcoshp["wavelength"] = hipparcoshp.wavelength /10000
hipparcoshp.to_csv("hipparcos2hp.dat",index=False, float_format="%.8f")

# SPIZTER/IRAC spectral response data
irac1 = pd.read_csv("original_files/irac1.dat", sep="\s+")
irac1.to_csv("iracl1.dat", index=False, float_format="%.8f")
irac2 = pd.read_csv("original_files/irac2.dat", sep="\s+")
irac2.to_csv("iracl2.dat", index=False, float_format="%.8f")
irac3 = pd.read_csv("original_files/irac3.dat", sep="\s+")
irac3.to_csv("iracl3.dat", index=False, float_format="%.8f")
irac4 = pd.read_csv("original_files/irac4.dat", sep="\s+")
irac4.to_csv("iracl4.dat", index=False, float_format="%.8f")

# IRAS spectral response data
iras1 = pd.read_csv("original_files/iras1.dat", sep="\s+")
iras1.to_csv("iras1.dat", index=False, float_format="%.8f")
iras2 = pd.read_csv("original_files/iras2.dat", sep="\s+")
iras2.to_csv("iras2.dat", index=False, float_format="%.8f")
iras3 = pd.read_csv("original_files/iras3.dat", sep="\s+")
iras3.to_csv("iras3.dat", index=False, float_format="%.8f")
iras4 = pd.read_csv("original_files/iras4.dat", sep="\s+")
iras4.to_csv("iras4.dat", index=False, float_format="%.8f")

# reformat the APASS filter transmission data
apassb = pd.read_csv("original_files/johnsoncousinsb.dat")
apassb.loc[:, "wavelength"] = apassb["wavelength"]/10000
apassb.to_csv("apassb.dat", index=False, float_format="%.8f")
apassv = pd.read_csv("original_files/johnsoncousinsv.dat")
apassv.loc[:, "wavelength"] = apassv["wavelength"]/10000
apassv.to_csv("apassv.dat", index=False, float_format="%.8f")
apassg = pd.read_csv("original_files/sdssg.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassg["wavelength"] =  apassg.wavelength/10000
apassg.to_csv("apassg.dat", index=False,float_format="%.8f")
apassr = pd.read_csv("original_files/sdssr.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassr["wavelength"] =  apassr.wavelength/10000
apassr.to_csv("apassr.dat", index=False,float_format="%.8f")
apassi = pd.read_csv("original_files/sdssi.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassi["wavelength"] =  apassi.wavelength/10000
apassi.to_csv("apassi.dat", index=False,float_format="%.8f")
apassg = pd.read_csv("original_files/sdssg.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassg["wavelength"] =  apassg.wavelength/10000
apassg.to_csv("apassg.dat", index=False,float_format="%.8f")
apassr = pd.read_csv("original_files/sdssr.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassr["wavelength"] =  apassr.wavelength/10000
apassr.to_csv("apassr.dat", index=False,float_format="%.8f")
apassi = pd.read_csv("original_files/sdssi.dat", sep="\s+", usecols=(0,1), header=None, skiprows=6,names=["wavelength", "response"])
apassi["wavelength"] =  apassi.wavelength/10000
apassi.to_csv("apassi.dat", index=False,float_format="%.8f")




# reformat Johnson-Cousins UBVRI filter transmission data
johnsoncousinsu = pd.read_csv("original_files/johnsoncousinsu.dat")
johnsoncousinsu.loc[:, "wavelength"] = johnsoncousinsu["wavelength"]/10000
johnsoncousinsu.to_csv("johnsoncousinsu.dat", index=False, float_format="%.8f")
johnsoncousinsb = pd.read_csv("original_files/johnsoncousinsb.dat")
johnsoncousinsb.loc[:, "wavelength"] = johnsoncousinsb["wavelength"]/10000
johnsoncousinsb.to_csv("johnsoncousinsb.dat", index=False, float_format="%.8f")
johnsoncousinsv = pd.read_csv("original_files/johnsoncousinsv.dat")
johnsoncousinsv.loc[:, "wavelength"] = johnsoncousinsv["wavelength"]/10000
johnsoncousinsv.to_csv("johnsoncousinsv.dat", index=False, float_format="%.8f")
johnsoncousinsr = pd.read_csv("original_files/johnsoncousinsr.dat")
johnsoncousinsr.loc[:, "wavelength"] = johnsoncousinsr["wavelength"]/10000
johnsoncousinsr.to_csv("johnsoncousinsr.dat", index=False, float_format="%.8f")
johnsoncousinsi = pd.read_csv("original_files/johnsoncousinsi.dat")
johnsoncousinsi.loc[:, "wavelength"] = johnsoncousinsi["wavelength"]/10000
johnsoncousinsi.to_csv("johnsoncousinsi.dat", index=False, float_format="%.8f")

# Pan-STARRs
# take filter profiles from: https://ipp.ifa.hawaii.edu/ps1.filters/
columns = ["wavelength", "g", "r", "i", "z", "y"]
ps1 = pd.read_csv("original_files/panstarrsps1.dat", sep="\s+", header=None, skiprows=26, usecols=(0,2,3,4,5,6), names=columns)
panstarrsps1g = ps1[["wavelength", "g"]].rename(columns={"g":"response"})
panstarrsps1g["wavelength"] = panstarrsps1g["wavelength"]/1000
panstarrsps1g = panstarrsps1g.sort_values(by=["wavelength"]).reset_index(drop=True)
panstarrsps1g.to_csv("panstarrsps1g.dat", index=False, float_format="%.8f") 
panstarrsps1r = ps1[["wavelength", "r"]].rename(columns={"r":"response"})
panstarrsps1r["wavelength"] = panstarrsps1r["wavelength"]/1000
panstarrsps1r = panstarrsps1r.sort_values(by=["wavelength"]).reset_index(drop=True)
panstarrsps1r.to_csv("panstarrsps1r.dat", index=False, float_format="%.8f") 
panstarrsps1i = ps1[["wavelength", "i"]].rename(columns={"i":"response"})
panstarrsps1i["wavelength"] = panstarrsps1i["wavelength"]/1000
panstarrsps1i = panstarrsps1i.sort_values(by=["wavelength"]).reset_index(drop=True)
panstarrsps1i.to_csv("panstarrsps1i.dat", index=False, float_format="%.8f") 
panstarrsps1z = ps1[["wavelength", "z"]].rename(columns={"z":"response"})
panstarrsps1z["wavelength"] = panstarrsps1z["wavelength"]/1000
panstarrsps1z = panstarrsps1z.sort_values(by=["wavelength"]).reset_index(drop=True)
panstarrsps1z.to_csv("panstarrsps1z.dat", index=False, float_format="%.8f") 
panstarrsps1y = ps1[["wavelength", "y"]].rename(columns={"y":"response"})
panstarrsps1y["wavelength"] = panstarrsps1y["wavelength"]/1000
panstarrsps1y = panstarrsps1y.sort_values(by=["wavelength"]).reset_index(drop=True)
panstarrsps1y.to_csv("panstarrsps1y.dat", index=False, float_format="%.8f") 

# The filter transmission curves of the SkyMapper system is take from Bessell+2011 (arxive source file).
f=open("original_files/skymapper_bessell_2011.dat","r")
lines=f.readlines()
# define pandas dataframe for each file.
skymapperu = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
skymapperv = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
skymapperg = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
skymapperr = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
skymapperi = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
skymapperz = pd.DataFrame(index=range(len(lines)-1), columns=["wavelength", "response"])
# load filter profiles from the raw file.
for i, x in enumerate(lines[1:]):
    skymapperu.loc[i, "wavelength"] = x[0:6]
    skymapperu.loc[i, "response"] = x[6:12]
    skymapperv.loc[i, "wavelength"] = x[20:26]
    skymapperv.loc[i, "response"] = x[26:33]
    skymapperg.loc[i, "wavelength"] = x[40:46]
    skymapperg.loc[i, "response"] = x[46:53]
    skymapperr.loc[i, "wavelength"] = x[60:66]
    skymapperr.loc[i, "response"] = x[66:73]
    skymapperi.loc[i, "wavelength"] = x[80:86]
    skymapperi.loc[i, "response"] = x[86:93]
    skymapperz.loc[i, "wavelength"] = x[100:106]
    skymapperz.loc[i, "response"] = x[106:113]                    
f.close()
# u band
skymapperu[["wavelength", "response"]] = skymapperu[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperu = skymapperu.dropna()
skymapperu.loc[:, "wavelength"] = skymapperu.loc[:, "wavelength"]/1e4
skymapperu = skymapperu.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperu.to_csv("skymapperu.dat", index=False, float_format="%.8f") 
# v band
skymapperv[["wavelength", "response"]] = skymapperv[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperv = skymapperv.dropna()
skymapperv.loc[:, "wavelength"] = skymapperv.loc[:, "wavelength"]/1e4
skymapperv = skymapperv.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperv.to_csv("skymapperv.dat", index=False, float_format="%.8f") 
# g band
skymapperg[["wavelength", "response"]] = skymapperg[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperg = skymapperg.dropna()
skymapperg.loc[:, "wavelength"] = skymapperg.loc[:, "wavelength"]/1e4
skymapperg = skymapperg.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperg.to_csv("skymapperg.dat", index=False, float_format="%.8f") 
# r band
skymapperr[["wavelength", "response"]] = skymapperr[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperr = skymapperr.dropna()
skymapperr.loc[:, "wavelength"] = skymapperr.loc[:, "wavelength"]/1e4
skymapperr = skymapperr.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperr.to_csv("skymapperr.dat", index=False, float_format="%.8f") 
# i band
skymapperi[["wavelength", "response"]] = skymapperi[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperi = skymapperi.dropna()
skymapperi.loc[:, "wavelength"] = skymapperi.loc[:, "wavelength"]/1e4
skymapperi = skymapperi.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperi.to_csv("skymapperi.dat", index=False, float_format="%.8f") 
# z band
skymapperz[["wavelength", "response"]] = skymapperz[["wavelength", "response"]].apply(pd.to_numeric, errors='coerce')
skymapperz = skymapperz.dropna()
skymapperz.loc[:, "wavelength"] = skymapperz.loc[:, "wavelength"]/1e4
skymapperz = skymapperz.sort_values(by=["wavelength"]).reset_index(drop=True)
skymapperz.to_csv("skymapperz.dat", index=False, float_format="%.8f") 

