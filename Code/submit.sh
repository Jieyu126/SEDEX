#!/bin/bash
# (1) set environment, skip this step if it has been set up.
# . ~/.bashrc_yujie  

# (2) configure inlist. Recommond to use the default values.
echo "Configuring inlist"
python inlist.py

# (3) download input photometry
echo "Fetching input photometry via Gaia Archive, given Gaia DR3 source IDs"
python 3_Prepare_Input_Params_Photometry.py

# (4) photometric outlier clipping by carrying out blackbody fitting.
echo "Photometric outlier clipping (blackbody fitting)"
python 4_BlackBodyFit.py

# (5) combine the blackbody fitting results (good bandpasses).
echo "Combine blackbody-fitting results"
python 5_CollectBlackBodyFits_Parallel.py 

# (6) perform SED fitting
echo "Perform SED fitting"
python 6_SedFit.py 

# (7) Combine the SED-fitting results. A table is generated with all the derived parameter estimates and input parameters.
echo "Collecting the SED-fitting results."
python 7_CollectSEDFits_Parallel.py