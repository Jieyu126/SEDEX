import http.client as httplib
from urllib.request import urlretrieve
from urllib.parse import quote as urlencode
import json
import sys
from astroquery.mast import Catalogs
import pandas as pd 
import numpy as np
from astroquery.simbad import Simbad
import lightkurve as lk
import ipyaladin.aladin_widget as ipyal
import astropy.units as u
from astropy.table import QTable
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import astropy.units as u
import extinction
import dust_extinction.parameter_averages as extLaw






class StaridResolver():
    """
    Aim: resolve star ID to coordinates. For more details, refer to https://ps1images.stsci.edu/ps1_dr2_api.html
    Example: 
        obj = StaridResolver(starID="12 Oph")
        ra, dec = obj.resolve()
    """
    def __init__(self, starID=""):
        self.name=starID
    
    def mastQuery(self, request):
        """Perform a MAST query.
        Parameters
        ----------
        request (dictionary): The MAST request json object
        Returns head,content where head is the response HTTP headers, and content is the returned data
        """
        server='mast.stsci.edu'
        # Grab Python Version 
        version = ".".join(map(str, sys.version_info[:3]))
        # Create Http Header Variables
        headers = {"Content-type": "application/x-www-form-urlencoded",
                "Accept": "text/plain",
                "User-agent":"python-requests/"+version}
        # Encoding the request as a json string
        requestString = json.dumps(request)
        requestString = urlencode(requestString)
        # opening the https connection
        conn = httplib.HTTPSConnection(server)
        # Making the query
        conn.request("POST", "/api/v0/invoke", "request="+requestString, headers)
        # Getting the response
        resp = conn.getresponse()
        head = resp.getheaders()
        content = resp.read().decode('utf-8')
        # Close the https connection
        conn.close()
        return head,content

    def resolve(self):
        """Get the RA and Dec for an object using the MAST name resolver
        Parameters
        ----------
        name (str): Name of object
        Returns RA, Dec tuple with position"""
        # Since this routine does not resolve K2 stars. The following 
        # part is added by Jie Yu, October 02, 2021.
        if "EPIC" in self.name: 
            search_result = lk.search_lightcurve(self.name, author='k2')
            ra, dec = search_result.table.columns["s_ra"][0], search_result.table.columns["s_dec"][0]
            return ra, dec
        else: 
            resolverRequest = {'service':'Mast.Name.Lookup',
                            'params':{'input':self.name,
                                        'format':'json'
                                        },
                            }
            headers,resolvedObjectString = self.mastQuery(resolverRequest)
            resolvedObject = json.loads(resolvedObjectString)
            # The resolver returns a variety of information about the resolved object, 
            # however for our purposes all we need are the RA and Dec
            try:
                objRa = resolvedObject['resolvedCoordinate'][0]['ra']
                objDec = resolvedObject['resolvedCoordinate'][0]['decl']
            except IndexError as e:
                raise ValueError("Unknown object '{}'".format(self.name))
            return (objRa, objDec)     





class photoNeighbours():
    def __init__(self, starID=None, ra=None, dec=None, fov=5/60, gaia_search_radius=1.0):
        self.starID = starID
        self.fov = fov
        self.gaia_search_radius = gaia_search_radius
        self.ra = ra
        self.dec = dec
        
    def download_gaia_data(self):
        if self.starID is not None:
            obj = StaridResolver(starID=self.starID)
            self.ra, self.dec = obj.resolve()  
        # check Gaia stars near the target and fetch Gaia data, 
        # use astroquery: https://astroquery.readthedocs.io/en/latest/gaia/gaia.html
        coord = SkyCoord(ra=self.ra, dec=self.dec, unit=(u.degree, u.degree), frame='icrs')
        radius = u.Quantity(self.gaia_search_radius, u.deg)
        j = Gaia.cone_search_async(coord, radius)
        r = j.get_results()
        self.r=r.to_pandas() 

    def aladin_show(self):
        # survey: from which choose an image to display. For more, check https://aladin.u-strasbg.fr/hips/list
        # fov: field of view in decimal degrees
        # target: you can provide either star name directly or its ra and dec
        # reticle_color: color of the text and marker.
        # for more arguments, see: https://aladin.u-strasbg.fr/AladinLite/doc/API/
        aladin= ipyal.Aladin(survey= 'P/2MASS/color', 
                             fov= self.fov, target="{:} {:}".format(self.ra, self.dec), 
                             reticle_color='#ff89ff')
        self.aladin = aladin
        return self.aladin
  
    def addons(self):   
        # add any columns from GAIA EDR3 that you'd like to inspect on the image above.
        # After running this cell, and you can click the boxes to view the data.
        # If the color of boxes are not clearly visible, run this cell several times until you're satisfied with the color
        # run r.columns to check column names.
        t = QTable([self.r.ra, self.r.dec, self.r.source_id, self.r.parallax, self.r.parallax_error, self.r.dist*3600, self.r.phot_g_mean_mag], 
                   names=('ra', 'dec', 'source_id', 'parallax', 'parallax_error', 'angdist', 'phot_g_mean_mag'))
        self.aladin.add_table(t)

        # add circles with radii equal to integer multiples of 4''
        self.aladin.add_overlay_from_stcs('Circle ICRS {:} {:} {:}'.format(self.ra, self.dec, 4/3600), {'color': 'red'})
        self.aladin.add_overlay_from_stcs('Circle ICRS {:} {:} {:}'.format(self.ra, self.dec, 8/3600), {'color': 'red'})
        self.aladin.add_overlay_from_stcs('Circle ICRS {:} {:} {:}'.format(self.ra, self.dec, 12/3600), {'color': 'red'})

        # choose coordinate systems: icrs, galactic
        # aladin.coo_frame= 'icrs'
        # aladin.coo_frame= 'galactic'  




def InfltError(InputPhotoMetry=None):
    InfltErrorTable = pd.DataFrame(columns=["band", "sys", "infltNumb"])
    InfltErrorTable = InfltErrorTable.set_index('band') 
    InfltErrorTable.loc["gaiag", :] = 0.02, 3
    InfltErrorTable.loc["gaiabp", :] = 0.02, 3
    InfltErrorTable.loc["gaiarp", :] = 0.02, 3
    InfltErrorTable.loc["hipparcos2hp", :] = 0.02, 3
    InfltErrorTable.loc["tycho2bt", :] = 0.02, 3
    InfltErrorTable.loc["tycho2vt", :] = 0.02, 3
    InfltErrorTable.loc["twomassj", :] = 0.02, 3
    InfltErrorTable.loc["twomassh", :] = 0.02, 3
    InfltErrorTable.loc["twomassk", :] = 0.02, 3
    InfltErrorTable.loc["wisew1", :] = 0.03, 3
    InfltErrorTable.loc["wisew2", :] = 0.03, 3
    InfltErrorTable.loc["wisew3", :] = 0.03, 3
    InfltErrorTable.loc["wisew4", :] = 0.03, 3
    InfltErrorTable.loc["johnsoncousinsu", :] = 0.06, 3
    InfltErrorTable.loc["johnsoncousinsb", :] = 0.06, 3
    InfltErrorTable.loc["johnsoncousinsv", :] = 0.06, 3
    InfltErrorTable.loc["johnsoncousinsr", :] = 0.06, 3
    InfltErrorTable.loc["johnsoncousinsi", :] = 0.06, 3
    InfltErrorTable.loc["sdssu", :] = 0.06, 3
    InfltErrorTable.loc["sdssg", :] = 0.06, 3
    InfltErrorTable.loc["sdssr", :] = 0.06, 3
    InfltErrorTable.loc["sdssi", :] = 0.06, 3
    InfltErrorTable.loc["sdssz", :] = 0.06, 3
    InfltErrorTable.loc["panstarrsps1g", :] = 0.06, 3
    InfltErrorTable.loc["panstarrsps1r", :] = 0.06, 3
    InfltErrorTable.loc["panstarrsps1i", :] = 0.06, 3
    InfltErrorTable.loc["panstarrsps1z", :] = 0.06, 3
    InfltErrorTable.loc["panstarrsps1y", :] = 0.06, 3
    InfltErrorTable.loc["apassb", :] = 0.06, 3
    InfltErrorTable.loc["apassv", :] = 0.06, 3
    InfltErrorTable.loc["apassg", :] = 0.06, 3
    InfltErrorTable.loc["apassr", :] = 0.06, 3
    InfltErrorTable.loc["apassi", :] = 0.06, 3   
    InfltErrorTable.loc["skymapperu", :] = 0.06, 3
    InfltErrorTable.loc["skymapperv", :] = 0.06, 3
    InfltErrorTable.loc["skymapperg", :] = 0.06, 3
    InfltErrorTable.loc["skymapperr", :] = 0.06, 3
    InfltErrorTable.loc["skymapperi", :] = 0.06, 3    
    InfltErrorTable.loc["skymapperz", :] = 0.06, 3    

    #update photometric errors
    for band in InputPhotoMetry.columns:
        if band in InfltErrorTable.index:
            # in the case where errors are available
            InputPhotoMetry.loc[:, band+"_e"] = (InfltErrorTable.loc[band, "sys"]**2 + InputPhotoMetry.loc[:, band+"_e"]**2)**0.5
            # in the case where errors are unavailable
            rep_val = InfltErrorTable.loc[band, "sys"] * InfltErrorTable.loc[band, "infltNumb"]
            InputPhotoMetry.loc[:, band+"_e"] = InputPhotoMetry.loc[:, band+"_e"].fillna(rep_val).copy()
            # if magnitudes in the current band are NaN, do not update their errors with numerical values but keep them as NaN
            subs = ~np.isfinite(np.array(InputPhotoMetry.loc[:, band], dtype=float)) # change nan from object to float
            InputPhotoMetry.loc[subs, band+"_e"] = np.nan
    return InputPhotoMetry



def Saturation_Reject(InputPhotoMetry):
    # references:
    # SDSS: https://classic.sdss.org/dr7/instruments/technicalPaper/index.html
    saturateTable = pd.DataFrame(columns=["band", "cutoffmag"])
    saturateTable = saturateTable.set_index('band') 
    saturateTable.loc["sdssu", "cutoffmag"] = 13
    saturateTable.loc["sdssg", "cutoffmag"] = 14
    saturateTable.loc["sdssr", "cutoffmag"] = 14
    saturateTable.loc["sdssi", "cutoffmag"] = 14
    saturateTable.loc["sdssz", "cutoffmag"] = 12
    # saturateTable.loc["gaiag", :] = 0.02, 3
    # saturateTable.loc["gaiabp", :] = 0.02, 3
    # saturateTable.loc["gaiarp", :] = 0.02, 3
    # saturateTable.loc["hipparcos2hp", :] = 0.02, 3
    # saturateTable.loc["tycho2bt", :] = 0.02, 3
    # saturateTable.loc["tychovbt", :] = 0.02, 3
    # saturateTable.loc["twomassj", :] = 0.02, 3
    # saturateTable.loc["twomassh", :] = 0.02, 3
    # saturateTable.loc["twomassk", :] = 0.02, 3
    # saturateTable.loc["wisew1", :] = 0.03, 3
    # saturateTable.loc["wisew2", :] = 0.03, 3
    # saturateTable.loc["wisew3", :] = 0.03, 3
    # saturateTable.loc["wisew4", :] = 0.03, 3
    # saturateTable.loc["johnsoncousinsu", :] = 0.06, 3
    # saturateTable.loc["johnsoncousinsb", :] = 0.06, 3
    # saturateTable.loc["johnsoncousinsv", :] = 0.06, 3
    # saturateTable.loc["johnsoncousinsr", :] = 0.06, 3
    # saturateTable.loc["johnsoncousinsi", :] = 0.06, 3
    # saturateTable.loc["panstarrsps1g", :] = 0.06, 3
    # saturateTable.loc["panstarrsps1r", :] = 0.06, 3
    # saturateTable.loc["panstarrsps1i", :] = 0.06, 3
    # saturateTable.loc["panstarrsps1z", :] = 0.06, 3
    # saturateTable.loc["panstarrsps1y", :] = 0.06, 3
    # saturateTable.loc["apassb", :] = 0.06, 3
    # saturateTable.loc["apassv", :] = 0.06, 3
    # saturateTable.loc["apassg", :] = 0.06, 3
    # saturateTable.loc["apassr", :] = 0.06, 3
    # saturateTable.loc["apassi", :] = 0.06, 3   
    # saturateTable.loc["skymapperu", :] = 0.06, 3
    # saturateTable.loc["skymapperv", :] = 0.06, 3
    # saturateTable.loc["skymapperg", :] = 0.06, 3
    # saturateTable.loc["skymapperr", :] = 0.06, 3
    # saturateTable.loc["skymapperi", :] = 0.06, 3    
    # saturateTable.loc["skymapperz", :] = 0.06, 3    
    # use a consevative threshold for each band by lowering the standard reference value by 2 mag.
    saturateTable.loc[:, "cutoffmag"] -= 2 
    
    # check if photometry in certain bands are saturations, if so, set the magntidues to NaN.
    subs = np.where(np.isin(InputPhotoMetry.index, saturateTable.index))[0]
    columns = InputPhotoMetry.index[subs].values
    for i in columns:
        if np.isfinite(InputPhotoMetry.loc[i, "mag"]):
            if InputPhotoMetry.loc[i, "mag"]<saturateTable.loc[i, "cutoffmag"]:
                InputPhotoMetry.loc[i, "mag"] = np.nan
                InputPhotoMetry.loc[i, "magerr"] = np.nan
    return InputPhotoMetry    



def target_data(id=None, radius=5.0, catalog="TIC", Lprint=False):
    """downloading target properties of TIC /and Gaia, as the below example  
       table = target_data(id = "KIC 1435467", catalog="TIC")
       table = target_data(id = "KIC 1435467", catalog="TIC", Lprint=True) # print out properties for the enquiry target
       table = target_data(id = "KIC 1435467", catalog="Gaia") # search for Gaia data including parallaxes, the defaul version is the highest.
    """
    # convert search radius in arcsec to degree since the default is degree
    radius = radius/60/60 
    # Query the TESS Input Catalog centered on HD 209458 with a 0.2 degree radius.
    catalogTIC = Catalogs.query_object(id, radius=radius, catalog=catalog).to_pandas()
    if len(catalogTIC)>0:
        starid = catalogTIC.loc[catalogTIC.dstArcSec.idxmin(), "ID"]
        catalogTIC=catalogTIC[catalogTIC.ID==starid].reset_index(drop=True).copy()
        # print out quiry results.
        if Lprint==True:
            # Print out the number of returned rows. 
            print("Number of TIC objects within %f deg of %s: %u" % (radius, starid, len(catalogTIC)))
            # print out target parameters.
            for i in catalogTIC.columns: print('{:s}: {:}'.format(i, catalogTIC.loc[0, i]))
        return catalogTIC
    else:
        print("No TIC objects within %f deg of %s" % (radius, id))


def Simbad_IDs(starID="", DR="DR2"):
    """fetching a target idenfifiers over various catalogs using SIMBAD.
       Example1: Simbad_IDs(id="HD 3457") 
       Example2: Simbad_IDs(id="2MASS J00373050+0308070")
       Example3: Simbad_IDs(id="Gaia DR2 2550687966499222144"] 
       Example4: to know the naming method, just run Simbad_IDs(id="HD 3457") and check out the output.
    """
    result = Simbad.query_objectids(starID).to_pandas()
    if len(result)>0:
        result.loc[:, "ID"] = result.ID.str.decode("utf8")
        if ~result.loc[:, "ID"].str.contains("Gaia").any():
            # Ocassionally, Simbad does not provide Gaia source_id with a Given star ID. 
            # This may change if you use a different star ID for the same star
            for i in range(len(result)):
                    result = Simbad.query_objectids(result.loc[0, "ID"]).to_pandas()
                    result.loc[:, "ID"] = result.ID.str.decode("utf8")
                    if result.loc[:, "ID"].str.contains("Gaia").any(): break  
        if DR=="DR3":          
            if ~result.loc[:, "ID"].str.contains("Gaia3").any():
                # in case Gaia DR3 does not included but Gaia DR2 ID do, search again using Gaia DR2 ID.
                result = result[result.loc[:, "ID"].str.contains("Gaia")].reset_index(drop=True)
                if len(result)>0: 
                    result = Simbad.query_objectids(result.loc[len(result)-1, "ID"]).to_pandas()
                    result.loc[:, "ID"] = result.ID.str.decode("utf8")
    return result


def ExtLaws(filters, law="F19", Rv=3.1, fA=1.0, Lprint=False):
    """
    Input:
        wavelength: in units of micron, in the range of 0.3 <= x <= 8.7, x is 1/wavelength
        law: "CCM89", "O94", "F99", "F04", "M14", "G16", "F19"
        fA: only for G16, the default is 1.0.
        For more details, refer to https://dust-extinction.readthedocs.io/en/stable/
    Output: A_lambda / A_V
    Note: the typic wavelength range in dust_extinction is 0.3 <= x <= 8.7, x has units 1/micron. 
        This means WISE W1 and W2 bands are not covered, but these two are commonly used for the 
        SED fiting. As thus, CCM89 is used and implemented using another python package extinction 
          https://extinction.readthedocs.io/en/latest/. In other words, wavelength longward of 1/0.3 
          is calculated with "extinction". In fact, CCM89 is very close to F99.   
    """
    # ext = pd.DataFrame(columns=["band", "wave", "A"])
    # ext.loc[:, "band"] = filters.band
    # ext.loc[:, "wave"] = filters.
    # use dust_extinction
    filters.loc[:, "A"] = np.nan
    subs = np.where((1/8.7 <= filters.loc[:, "lambda"]) & (filters.loc[:, "lambda"] <= 1/0.3))[0]

    if len(subs)>0:
        if law=="CCM89": filters.loc[subs, "A"] = extLaw.CCM89(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "O94": filters.loc[subs, "A"] =  extLaw.O94(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "F99": filters.loc[subs, "A"] =  extLaw.F99(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "F04": filters.loc[subs, "A"] =  extLaw.F04(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "M14": filters.loc[subs, "A"] =  extLaw.M14(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "G16": filters.loc[subs, "A"] =  extLaw.G16(RvA = Rv, fA = fA)(filters.loc[subs, "lambda"].values * u.micron)
        if law == "F19": filters.loc[subs, "A"] =  extLaw.F19(Rv = Rv)(filters.loc[subs, "lambda"].values * u.micron)
    # use extinction
    subs = np.where(((1/8.7 > filters.loc[:, "lambda"]) | (filters.loc[:, "lambda"] > 1/0.3)))[0]
    if len(subs)>0: 
        # wavelength in units of Angstroms
        if law=="O94": 
            filters.loc[subs, "A"] = extinction.odonnell94(filters.loc[subs, "lambda"].values*1e4, 1.0, Rv)
        elif law=="F99":
            filters.loc[subs, "A"] = extinction.fitzpatrick99(filters.loc[subs, "lambda"].values*1e4, 1.0, Rv)
        else:    
            if Lprint:
                print("Note that you're using CCM89 extinction law for these bands:")
                print(filters.iloc[subs, :-1])
            filters.loc[subs, "A"] = extinction.ccm89(filters.loc[subs, "lambda"].values*1e4, 1.0, Rv) 
    # combine    
    return filters.A.values







def access_crossmatch_from_Gaia(query, angDist=3, upload_filename="", user_filename="", format="csv", out_filename=""):
    """In the current version of Gaia archive (as of 2021 July), star IDs for databases of Gaia ERD3,
    Tycho2, APASS DR9, SDSS DR13, Skymapper DR2, Pan-starrs DR1, ALLWISE, TWOMASS are available. But the photometry
    is only available from Gaia archive for the first five but not for the last three. In this unvailable case, star IDs
    are obtained from GAIA crossmatch, and then used to select targets from the cross-matches access from CDS/XMATCh using 
    astroquery.XMATCH implemented in Python.  
    """
    # login
    Gaia.login(user='jyu01', password='Yping:126:')
    if upload_filename != "": 
        try:
            job = Gaia.delete_user_table(user_filename)
            print("Table already uploaded to the Gaia Archive. So delete it and upload a newer version.") 
        except:
            print("No table already uploaded to the Gaia Archive. So don't need to delete.")   
        job = Gaia.upload_table(upload_resource=upload_filename, table_name=user_filename, format=format)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Access TWOMASS IDs.
    job = Gaia.launch_job_async(query=query)
    results = job.get_results().to_pandas()
    if ("angular_distance" in results.columns) & ("starid" in results.columns): 
        results = results.loc[results.groupby("starid")["angular_distance"].idxmin()].reset_index(drop=True)
        results = results[results.angular_distance<=angDist].reset_index(drop=True)
        results.loc[:, "starid"] = results.loc[:, "starid"].astype("str")  
    results.to_csv(out_filename,index=False)


def GaiaID_convertor(upload_filename="", user_filename="", out_filename="", input="GaiaDR2", output="GaiaEDR3"):
    """ convert Gaia DR2 IDs to EDR3 IDs
        upload_filename: filename of the file to be uploaded to Gaia Archive containing Gaia IDs.
        user_filename: the name of the table to appear in the dropdown of User tables on Gaia Archive
        out_filename: the filename of the file after cross matching.
        input: the input Gaia IDs, such as "GaiaDR2"
        output: the output Gaia IDs, such as "GaiaEDR3"
        Note that the column name for the star IDs mush be starID.
    """
    if ((input=="GaiaDR2") & (output=="GaiaEDR3")):    
        query = ("SELECT targets.*, gaia.dr3_source_id, gaia.angular_distance \
            FROM user_jyu01.{:s} AS targets \
            INNER JOIN gaiadr3.dr2_neighbourhood AS gaia \
            ON gaia.dr2_source_id = targets.starid ".format(user_filename))   
        access_crossmatch_from_Gaia(query, upload_filename=upload_filename, user_filename=user_filename, out_filename=out_filename)


def clean_Gaia_distance(in_filename="", out_filename=""):
    # Where photogeometric distances are not available, they replaced by geometric distances.
    dist = pd.read_csv(in_filename, dtype={"source_id":str})
    dist.loc[:, "r_phogeo"] = dist.r_med_photogeo.copy()
    dist.loc[:, "r_phogeo_err"] =  (dist.r_hi_photogeo.copy()-dist.r_lo_photogeo.copy())/2
    dist.loc[:, "r_geo"] = dist.r_med_geo.copy()
    dist.loc[:, "r_geo_err"] =  (dist.r_hi_geo.copy()-dist.r_lo_geo.copy())/2
    dist.loc[:,"d"] = dist.loc[:, "r_phogeo"].copy()
    dist.loc[:, "derr"] = dist.loc[:, "r_phogeo_err"].copy()
    subs = np.where(pd.isna(dist.loc[:, "r_phogeo"]))[0]
    dist.loc[subs, "d"] = dist.loc[subs, "r_geo"].copy()
    dist.loc[subs, "derr"] = dist.loc[subs, "r_geo_err"].copy() 
    dist = dist[["source_id", "d", "derr"]].copy()
    dist.to_csv(out_filename, index=False, float_format=("%.4f"))

    
    
def clean_Gaia_photometry(in_filename="", out_filename=""):
    # correct G-band photometry and flux.
    gaiaedr3 = pd.read_csv(in_filename, dtype={"starid":str, "source_id":str})
    gmag_corr, gflux_corr = kepler.correct_gband(gaiaedr3.bp_rp, gaiaedr3.astrometric_params_solved, gaiaedr3.phot_g_mean_mag, gaiaedr3.phot_g_mean_flux)
    gaiaedr3.loc[:, "phot_g_mean_mag"] = gmag_corr
    gaiaedr3.loc[:, "phot_g_mean_flux"] = gflux_corr
    # calculate magnitude errors in the G, Bp, Rp bands
    gflux = gaiaedr3.loc[:, "phot_g_mean_flux"]
    gfluxerr = gaiaedr3.loc[:, "phot_g_mean_flux_error"]
    bpflux = gaiaedr3.loc[:, "phot_bp_mean_flux"]
    bpfluxerr = gaiaedr3.loc[:, "phot_bp_mean_flux_error"]
    rpflux = gaiaedr3.loc[:, "phot_rp_mean_flux"]
    rpfluxerr = gaiaedr3.loc[:, "phot_rp_mean_flux_error"]
    gmagerr, bpmagerr, rpmagerr = kepler.Gaia_mag_errors(gflux=gflux, gfluxerr=gfluxerr, bpflux=bpflux, bpfluxerr=bpfluxerr, rpflux=rpflux, rpfluxerr=rpfluxerr)
    gaiaedr3.loc[:, "phot_g_mean_mag_error"] = gmagerr
    gaiaedr3.loc[:, "phot_bp_mean_mag_error"] = bpmagerr
    gaiaedr3.loc[:, "phot_rp_mean_mag_error"] = rpmagerr
    # select needed columns
    columns = ["starid", "source_id","phot_g_mean_mag","phot_g_mean_mag_error","phot_bp_mean_mag",\
        "phot_bp_mean_mag_error","phot_rp_mean_mag","phot_rp_mean_mag_error", "parallax", "parallax_error"]
    gaiaedr3 = gaiaedr3[columns]
    gaiaedr3.to_csv(out_filename, index=False, float_format="%.4f")


def access_polluters_from_Gaia(query, upload_filename="", user_filename="", format="csv", out_filename=""):
    """In the current version of Gaia archive (as of 2021 July), star IDs for databases of Gaia ERD3,
        Tycho2, APASS DR9, SDSS DR13, Skymapper DR2, Pan-starrs DR1, ALLWISE, TWOMASS are available. But the photometry
        is only available from Gaia archive for the first five but not for the last three. In this unvailable case, star IDs
        are obtained from GAIA crossmatch, and then used to select targets from the cross-matches access from CDS/XMATCh using 
        astroquery.XMATCH implemented in Python.  
    """
    # login
    Gaia.login(user='jyu01', password='Yping:126:')
    if upload_filename != "": job = Gaia.upload_table(upload_resource=upload_filename, table_name=user_filename, format=format)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    query = "SELECT targets.starID, gaia.ra, gaia.dec, distance(POINT('ICRS', gaia.ra, gaia.dec),POINT('ICRS', source.ra, source.dec)) AS AngDist, source.*\
            FROM user_jyu01.seismic AS targets \
            INNER JOIN gaiaedr3.dr2_neighbourhood AS Bneigb\
            ON Bneigb.dr2_source_id = targets.starID\
            INNER JOIN gaiaedr3.gaia_source AS gaia\
            ON gaia.source_id = Bneigb.dr3_source_id \
            JOIN gaiaedr3.gaia_source AS source\
            ON 1=CONTAINS(POINT('ICRS', gaia.ra, gaia.dec),CIRCLE('ICRS', source.ra, source.dec, 0.002222222222224))"
    job = Gaia.launch_job_async(query=query)
    results = job.get_results().to_pandas()
    results.to_csv(out_filename,index=False)




def download_twomass_photometry(targetlist=None, maxrec=500000, Lsplit=True, nchunks=10, filename=None):
    """
    Aim: download 2MASS photometry given 2MASS source ID that can be obtained from Gaia archive.
    Input: 
        targetlist: this is a pandas DataFrame that must include one column called original_ext_source_id that can be fetched from Gaia archive.
        maxrec: maximum number of records downloaded for each query. Note that if the target list is not split, then maxrec is specific for the 
                entire target list. This number has to be reasonably large in order not to miss any matched records. This argument can be twice 
                larger than the size of each chunk.
        Lsplit: whether to split the target list. This is necessary if the target list contains too many targets. Suggest to split it when the 
                size is great than 500000. 
        nchunks: number of chunks to split the target list into. It takes 20s if each chunk contains 200000 entries.                 
    Output:
        a file of twomass photometry. If the target list is splitted, the file is also combined.
    Example:
        download_twomass_photometry(targetlist=twomass.copy(), maxrec=500000, Lsplit=True, nchunks=10)
    """    
    t0 = time.time()    
    # step 1: download twomass photometry.    
    index = np.array(np.linspace(0, len(targetlist)-1, nchunks+1), dtype=np.int64)
    query="select *\
            from twomass \
            inner join tap_upload.t1 as mine \
            on twomass.mainid=mine.original_ext_source_id"
    if Lsplit==False:
        result = pyvo_query(tapurl="https://gaia.ari.uni-heidelberg.de/tap", query=query, table=targetlist, maxrec=maxrec)   
        result.to_csv(filename,index=False)             
    else:
        # split the input target list into chunks    
        for i in np.arange(len(index)-1):
            print("Downloading chunk {:d} out of {:d} chunks".format(i+1, len(index)-1))    
            if i==len(index)-2:
                chunk = targetlist.iloc[index[i]:, :].reset_index(drop=True).copy() 
            else:
                chunk = targetlist.iloc[index[i]:index[i+1], :].reset_index(drop=True).copy() 
            result = pyvo_query(tapurl="https://gaia.ari.uni-heidelberg.de/tap", query=query, table=chunk.copy(), maxrec=maxrec)   
            result.to_csv("chunk"+str(i)+".csv",index=False) 

        # step 2: combine chunks and output    
        filenames = np.array(glob.glob("chunk*.csv"))
        ind = np.argsort(np.array([re.search(r"\d+", i).group(0) for i in filenames], dtype=int))
        filenames = filenames[ind]
        result = pd.DataFrame()
        for ind, i in enumerate(filenames):
                print("Combine individual chunks: the {:d}th of {:d} chunks".format(ind+1, len(filenames)))    
                temp = pd.read_csv(i)
                result = pd.concat([result, temp])
        print("Saving twomass photometry to a file")
        result.to_csv(filename,index=False) 
        
        # step 3: remove chunks    
        for i in filenames: 
                print("Remove {:s}".format(i))
                os.system("rm "+i) 
        t1 = time.time()
        print("Time consumed: {:.1f}s".format(t1-t0))

def correct_gband(bp_rp, astrometric_params_solved, phot_g_mean_mag, phot_g_mean_flux):
    """
    Correct the G-band fluxes and magnitudes for the input list of Gaia EDR3 data.
    Note that this is copied from the public resource at https://github.com/agabrown/gaiaedr3-6p-gband-correction
    Parameters
    ----------
    
    bp_rp: float, numpy.ndarray
        The (BP-RP) colour listed in the Gaia EDR3 archive.
    astrometric_params_solved: int, numpy.ndarray
        The astrometric solution type listed in the Gaia EDR3 archive.
    phot_g_mean_mag: float, numpy.ndarray
        The G-band magnitude as listed in the Gaia EDR3 archive.
    phot_g_mean_flux: float, numpy.ndarray
        The G-band flux as listed in the Gaia EDR3 archive.
        
    Returns
    -------
    
    The corrected G-band magnitudes and fluxes. The corrections are only applied to
    sources with a 2-paramater or 6-parameter astrometric solution fainter than G=13, 
    for which a (BP-RP) colour is available.
    
    Example
    -------
    
    gmag_corr, gflux_corr = correct_gband(bp_rp, astrometric_params_solved, phot_g_mean_mag, phot_g_mean_flux)
    """

    if np.isscalar(bp_rp) or np.isscalar(astrometric_params_solved) or np.isscalar(phot_g_mean_mag) \
                    or np.isscalar(phot_g_mean_flux):
        bp_rp = np.float64(bp_rp)
        astrometric_params_solved = np.int64(astrometric_params_solved)
        phot_g_mean_mag = np.float64(phot_g_mean_mag)
        phot_g_mean_flux = np.float64(phot_g_mean_flux)
    
    if not (bp_rp.shape == astrometric_params_solved.shape == phot_g_mean_mag.shape == phot_g_mean_flux.shape):
        raise ValueError('Function parameters must be of the same shape!')
    
    do_not_correct = np.isnan(bp_rp) | (phot_g_mean_mag<13) | (astrometric_params_solved == 31)
    bright_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>=13) & (phot_g_mean_mag<=16)
    faint_correct = np.logical_not(do_not_correct) & (phot_g_mean_mag>16)
    bp_rp_c = np.clip(bp_rp, 0.25, 3.0)
    
    correction_factor = np.ones_like(phot_g_mean_mag)
    correction_factor[faint_correct] = 1.00525 - 0.02323*bp_rp_c[faint_correct] + \
        0.01740*np.power(bp_rp_c[faint_correct],2) - 0.00253*np.power(bp_rp_c[faint_correct],3)
    correction_factor[bright_correct] = 1.00876 - 0.02540*bp_rp_c[bright_correct] + \
        0.01747*np.power(bp_rp_c[bright_correct],2) - 0.00277*np.power(bp_rp_c[bright_correct],3)
    
    gmag_corrected = phot_g_mean_mag - 2.5*np.log10(correction_factor)
    gflux_corrected = phot_g_mean_flux * correction_factor
    
    return gmag_corrected, gflux_corrected

def Gaia_mag_errors(gflux=None, gfluxerr=None, bpflux=None, bpfluxerr=None, rpflux=None, rpfluxerr=None, version="EDR3"): 
    """
    Aim: calculate magnitude errors in Gaia EDR3 G, Bp, Rp bands.
    Input: 
        gflux: phot_g_mean_flux of Gaia EDR3 
        gfluxerr: phot_g_mean_flux_error of Gaia EDR3 
        bpflux: phot_bp_mean_flux of Gaia EDR3 
        bpfluxerr: phot_bp_mean_flux_error of Gaia EDR3 
        rpflux: phot_rp_mean_flux of Gaia EDR3 
        rpfluxerr: phot_rp_mean_flux_error of Gaia EDR3 
    Output:
        mag errors in the G, Bp, and Rp bands
    """    
    if version=="EDR3":
        gmagerr = np.sqrt((2.5*gfluxerr/gflux/np.log(10))**2 + (0.0027553202 )**2)
        bpmagerr = np.sqrt((2.5*bpfluxerr/bpflux/np.log(10))**2 + (0.0027901700)**2)
        rpmagerr = np.sqrt((2.5*rpfluxerr/rpflux/np.log(10))**2 + (0.0037793818)**2)
        return gmagerr, bpmagerr, rpmagerr

def Gaia_XP_RVS_Spectra(query, upload_filename="", user_filename="", format_upload_table="csv", format_download_table="fits", 
    out_filename="", retrieval_type="", data_structure="", data_release=""):
    """
    Aim: 
        download six different products: epoch photometry, medium- and low-resolution spectra, and probability 
        density distributions for the different astrophysical parameters.
    arguments: 
        retrieval_type: 'EPOCH_PHOTOMETRY', 'MCMC_GSPPHOT', 'MCMC_MSC', 'XP_SAMPLED', 'XP_CONTINUOUS', 'RVS', 'ALL'
        data_structure: 'INDIVIDUAL', 'COMBINED', 'RAW'
        data_release: 'Gaia DR3' (default), 'Gaia DR2'
        format: ‘csv’, ‘ecsv’,’json’,’votable_plain’ and ‘fits’
    Reference: 
        https://www.cosmos.esa.int/web/gaia-users/archive/datalink-products
    """
    Gaia.login(user='jyu01', password='Yping:126:')
    if upload_filename != "": 
        try:
            job = Gaia.delete_user_table(user_filename)
            print("Table already uploaded to the Gaia Archive. So delete it and upload a newer version.") 
        except:
            print("No table already uploaded to the Gaia Archive. So don't need to delete.")   
        job = Gaia.upload_table(upload_resource=upload_filename, table_name=user_filename, format=format_upload_table)


    # fetch spectra.
    job = Gaia.launch_job_async(query=query)
    results = job.get_results()
    datalink = Gaia.load_data(ids=results['source_id'], data_release=data_release, 
                            retrieval_type=retrieval_type, data_structure=data_structure, 
                            verbose=False, output_file=out_filename, format=format_download_table)
                            
