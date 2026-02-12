#!/bin/sh
#SBATCH --job-name=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=12:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue

module load "MATLAB/2023a"

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, exiting."
    exit 2
fi

echo "Starting job on $HOSTNAME"


cd /home/$USER/project/cluster_code/; matlab -nodisplay -r "pre_hemo "$1";exit";

if [ "$?" -ne 0 ]; then
    echo "Error in detrending script, exiting."
    exit 2
fi

cd /home/$USER/project/cluster_code/ ; sbatch hemo_correct.sh $1 $2
