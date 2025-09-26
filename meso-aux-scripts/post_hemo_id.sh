#!/bin/sh
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=2:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --output="./outs/post_hemo/%j.out"

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters, exiting."
    exit 1
fi

module load MATLAB/2023a


matlab -nodisplay -r "post_hemo_id "$1";exit";
