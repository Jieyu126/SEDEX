#!/bin/bash
#SBATCH --partition=helio
#SBATCH --qos=helio_default
#SBATCH --account=helio
#SBATCH --mem=5G
##SBATCH --time=00:20:00
#SBATCH --time=2-00:00:00
#SBATCH --mail-type=FAIL,REQUEUE,STAGE_OUT
#SBATCH --mail-user=yujie@mps.mpg.de
#SBATCH --cpus-per-task=1
#SBATCH --array=0-129
#SBATCH --job-name=SED_Entire
#SBATCH --error=/home/yujie/Documents/SoftLinks/Scratch/Projects/SEDFit/SEDEX/Data/Output_Fits/Seismic/SEDFits/MARCS/BatchMode/Log/%x_slurm%A_%a.err
#SBATCH --output=/home/yujie/Documents/SoftLinks/Scratch/Projects/SEDFit/SEDEX/Data/Output_Fits/Seismic/SEDFits/MARCS/BatchMode/Log/%x_slurm%A_%a.out
nstars=6957  # Number of stars






############################################################################################
# The following part is not required to be changed
############################################################################################
# cpus-per-task.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK 
# set environment.
. ~/.bashrc_yujie        
# to creat some folders.
srun python ./inlist.py  


###########################################################################################
# Part I (only this part needs to be reset): configure the following parameters.
ntasks=$SLURM_ARRAY_TASK_COUNT  # Number of tasks of this job
index=$SLURM_ARRAY_TASK_ID      # task index  
# specify directories for saving fiels.
if [ $fits = "BlackBodyFits" ]
then
    script="4_BlackBodyFit.py"
else
    script="6_SedFit.py"
fi

# extract information from the filename of the standard slurm output stated above.
input=`basename "$0"`
while IFS= read -r line
    do
        # search for the path for saving standard output
        b=${line:0:16}
        if [[ "$b" == "#SBATCH --output" ]]; then
            # determine the paths for saving log files and SED-fitting output
            logFilePattern=`echo "$line" | rev | cut -d "/" -f1 | rev`
            size=${#logFilePattern}
            path_log_files=${line:17:-$size-1} 
            name=`echo "$path_log_files" | rev | cut -d "/" -f1 | rev`
            size=${#name}
            path_output_files=${path_log_files:0:-$size-1} 
            # search for the sample, grid, and the module to run.
            sample=`echo "$path_log_files" | rev | cut -d "/" -f5 | rev`
            fits=`echo "$path_log_files" | rev | cut -d "/" -f4 | rev`
            grid=`echo "$path_log_files" | rev | cut -d "/" -f3 | rev`
        fi
    done < "$input"

# Creat the folder for your log files
mkdir -p $path_log_files

# ###########################################################################################
# Part II: delete all the files generated from previous batch jobs
rm -rf `find $path_log_files -name "*.*" | grep  -v "$SLURM_ARRAY_JOB_ID"`
rm -rf `find $path_output_files -name "*.csv" | grep  -v "$SLURM_ARRAY_JOB_ID"`



###########################################################################################
# Part III: run the python code
# staring index of the target list
N0=0                            
# number of stars to be analyzed for each task	
task_size=$((($nstars+$ntasks-1)/$ntasks))
# set lower and upper indices
index_start=$(($N0+$index*$task_size))
index_end=$(($N0+($index+1)*$task_size-1))
# update the upper index if its distance to the max index is less than the step size.
if [ $index_start -le $(($nstars-1)) ]; then
    # update the upper index if its distance to the max index is less than the step size.
    if [ $index_end -gt $(($nstars-1)) ]; then
        index_end=$(($nstars-1))
    fi
    # python $script $index_start $index_end $index
    srun python $script $index_start $index_end $index
fi



###########################################################################################
# Part IV: some instructions
# ref: https://slurm.schedmd.com/job_array.html
# ref: https://slurm.schedmd.com/sbatch.html#lbAH
# %x: --job-name
# %A: Job array's master job allocation number.
# %a: Job array ID (index) number.
# %j and %J give the same job id? (%j=%A+%a i guess)

# (1.0) submit a job
# sbatch submit_slurm.sh

# (1.1) find the number of your running tasks 
# squeue -u yujie -h -t running -r | wc -l

# (1.2) To cancel all jobs submitted to a partition:
# scancel --user=$USER --partition=<partition_name>  
# scancel --user=yujie --partition=helio

# (1.3) view information about jobs located in the Slurm scheduling queue.
# squeue --user=$USER --partition=<partition_name> 
# squeue --user=yujie --partition=helio

# (1.4) View information about Slurm nodes and partitions. 
# sinfo -p helio

# (1.5) One or more Jobs submitted to a partition can be canceled from the queue or while running:
# scancel jobID_1,...,jobID_N   e.g.: scancel 1499,1500

# (1.6) Cancel a specific task of a job:
# scancel jobID.task_K   e.g.: scancel 1701.5 to cancel task 5 of job 1701  

# (1.7) find max number of tasks that I can run
# sinfo -o%C #the last number in the output will be the total number of CPUs.

# (1.8) find the batch script of a task given JOBID
# scontrol show job=<JOBID>
# scontrol show job=3696080_81