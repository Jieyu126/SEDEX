#!/bin/bash
# (0) set environment, skip this step if it has been set up.
# . ~/.bashrc_yujie  

# (1) configure inlist. Recommond to use default values.
echo "Step 1: configure inlist"
python inlist.py

# (2) download input photometry
echo "Step 2: fetch input photometry via Gaia Archive, given Gaia DR3 source IDs"
python 3_Prepare_Input_Params_Photometry.py

# (3) photometric outlier clipping by carrying out blackbody fitting.
echo "Step 3: clip photometric outliers (blackbody fitting)"
python 4_BlackBodyFit.py

# (4) combine the blackbody fitting results (good bandpasses).
echo "Step 4: combine blackbody-fitting results"
python 5_CollectBlackBodyFits_Parallel.py 

# (5) perform SED fitting
echo "Step 5: perform SED fitting"
python 6_SedFit.py 

# (6) Combine the SED-fitting results. A table is generated with all the derived parameter estimates and input parameters.
echo "Step 6: collecte the SED-fitting results."
python 7_CollectSEDFits_Parallel.py