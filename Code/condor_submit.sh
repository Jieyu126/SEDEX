# number of jobs to be submitted.
njobs=160

# Spicify the sample
# sample=Seismic     
# nstars=6957
# sample=TOIs     
# nstars=200
# sample=Kepler     
# nstars=177652
# sample=APOGEE     
# nstars=610602
# sample=GALAH     
# nstars=572702
# sample=RAVE     
# nstars=403622
# sample=Queiroz2020   
# nstars=668478
sample=RVS
nstars=5529249


# Specify for blackbody fits or SED fits.
# script=5_BlackBodyFit.py
# fits=BlackBodyFits
script=7_SedFit.py
fits=SEDFits_MARCS
# script=7_SedFit.py
# fits=SEDFits_BOSZ

# include : ./test.sh |
# environment="newvar=$(var)"

# fits=SEDFits
filePath=/home/yujie/Documents/SoftLinks/Scratch/Projects/3DustMaps/SEDEX/BigData/$(sample)/Tables/$(fits)
logPath=/home/yujie/Documents/SoftLinks/Scratch/Projects/3DustMaps/SEDEX/BigData/$(sample)/Batch/$(fits)

# scripts of jobs
universe=vanilla
getenv=True
executable=/bin/bash
arguments=wrapper.sh $(process) $(script) $(filePath) $(logPath) $(njobs) $(nstars)

# log files
output=$(logPath)/slot$(process).out
error=$(logPath)/slot$(process).err
log=$(logPath)/slot$(process).log


# to use or to avoid certain machines
requirements=(machine!="seismo22.mps.mpg.de")&&(machine!="seismo26.mps.mpg.de")&&(machine!="seismo27.mps.mpg.de")
# requirements = (machine != "seismo22.mps.mpg.de")

# image_size in KB
image_size=300000 

# number of CPUs requested. If you submit you jobs through seismo16, then you can only use seismo17-31, 
# alternatively, if you submit your jobs through seismo1, then you can only use seismo2-15. This is because 
# seismo16-31 use Intel hardware, while seismo16-31 use AMD hardware.
# queue 150
queue $(njobs)



# how to run in the batch mode
# (0) set up inlist, and in ipython terminal type import inlist to creat some directories.

# (1) check the availability of servers 
# http://vlx04.mps.mpg.de/ganglia/?m=&r=hour&s=by%20name&hc=6

# (2) ssh seismo16 (submit jobs to seismo18-31 through seismo16)
# Currently, I am allowded to use the cluster seismo.
# To submit a job, log in to either seismo1 or seismo16, as these two 
# servers permit job submissions. Refer to this page for more details.
# https://www.mps.mpg.de/de/services/intranet/rechenzentrum/hardware/ServersStarsInterior
# ssh seismo1 (submit jobs to seismo2-15 through seismo1)

# (3) condor_submit submit.sh
# Two variables have to be set up. (1) image_size. It is necessary to estimate the total size of 
# figures to be generated and stored. (2) queue. This refers to the number of servers used to executing
# your job. NOTE that once this vaviable is set up, you have to assign the same value to the variable 
# in the script queue_numb. Dont forget to assign the variable "loop_size" in the script wrapper.sh with 
# the number of stars. # If you submit you jobs through seismo16, then you can only use seismo17-31, 
# alternatively, if you submit your jobs through seismo1, you can only use seismo2-15. This is because 
# seismo16-31 use Intel hardware, while seismo16-31 use AMD hardware.

# (4) some widely used commands to check the status of jobs submitted.
# reference: https://htcondor.readthedocs.io/en/latest/users-manual/managing-a-job.html
# condor_q
# condor_q -nobatch
# condor_status 
# condor_status -run 
# condor_status -submitters
# condor_status -constraint 'RemoteUser == "yujie@mps.mpg.de"'
# condor_rm <user>