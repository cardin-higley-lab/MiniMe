#!/bin/sh
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --time=00:10:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --array=1-256
#SBATCH --output=./outs/detrend/%j.out

module load MATLAB/2022a

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters, exiting."
    exit 1
fi

t=$SLURM_ARRAY_TASK_ID
echo $t

matlab -nodisplay -r "detrend_id "$1" "$t" ;exit";