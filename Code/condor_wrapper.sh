#!/bin/bash
# to load anaconda installed locally
. ~/.bashrc_yujie  

# server index       				
index=$1   
script=$2
pathFiles=$3
logFiles=$4
njobs=$5  
loop_size=$6

# initial index
N0=0





numb=$(ls "$pathFiles"/*.csv 2> /dev/null)
if [ ! -z "$numb" ]; then
    # in this case, there are files in there.
        for file in "$pathFiles"/*.csv; do
            # extract filename 
            filename="$(basename -- $file)"
            # get the index number from the filename
            subs=$(echo "$filename" | grep -o -E '[0-9]+')
            # if this number is larger than the number of jobs, delete the file.
            if [[ "$subs" -ge "$njobs" ]]; then
                # delete files
                rm -f $file
            fi    
        done
fi


numb=$(ls "$logFiles"/*.err 2> /dev/null)
if [ ! -z "$numb" ]; then
    echo "1"
    # in this case, there are files in there.
        for file in "$logFiles"/*.err; do
            # extract filename 
            filename="$(basename -- $file)"
            # get the index number from the filename
            subs=$(echo "$filename" | grep -o -E '[0-9]+')
            # if this number is larger than the number of jobs, delete the file.
            if [[ "$subs" -ge "$njobs" ]]; then
                # delete files
                rm -f $file
            fi    
        done
fi

numb=$(ls "$logFiles"/*.log 2> /dev/null)
if [ ! -z "$numb" ]; then
    # in this case, there are files in there.
        for file in "$logFiles"/*.log; do
            # extract filename 
            filename="$(basename -- $file)"
            # get the index number from the filename
            subs=$(echo "$filename" | grep -o -E '[0-9]+')
            # if this number is larger than the number of jobs, delete the file.
            if [[ "$subs" -ge "$njobs" ]]; then
                # delete files
                rm -f $file
            fi    
        done
fi

numb=$(ls "$logFiles"/*.out 2> /dev/null)
if [ ! -z "$numb" ]; then
    # in this case, there are files in there.
        for file in "$logFiles"/*.out; do
            # extract filename 
            filename="$(basename -- $file)"
            # get the index number from the filename
            subs=$(echo "$filename" | grep -o -E '[0-9]+')
            # if this number is larger than the number of jobs, delete the file.
            if [[ "$subs" -ge "$njobs" ]]; then
                # delete files
                rm -f $file
            fi    
        done
fi






# number of stars to be analyzed for each job	
job_size=$((($loop_size+$njobs-1)/$njobs))

# set lower and upper indices
index_start=$(($N0+$index*$job_size))
index_end=$(($N0+($index+1)*$job_size-1))

# update the upper index if its distance to the max index is less than the step size.
if [ $index_start -le $(($loop_size-1)) ]; then
    # update the upper index if its distance to the max index is less than the step size.
    if [ $index_end -gt $(($loop_size-1)) ]; then
        index_end=$(($loop_size-1))
    fi
    python $script $index_start $index_end $index
fi


