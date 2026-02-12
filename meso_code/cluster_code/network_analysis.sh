#!/bin/sh
#SBATCH --job-name=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=12:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue

module load MATLAB/2023a


echo "Starting job on $HOSTNAME"


cd /home/$USER/project/cluster_code/; matlab -nodisplay -r "cluster_network_analysis "$1"; exit";


