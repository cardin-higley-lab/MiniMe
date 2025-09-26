#!/bin/sh
#SBATCH --job-name=1-1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=1:00:00
#SBATCH --partition="scavenge"
#SBATCH --array=1-64
#SBATCH --output="./outs/hemo_correct/%j.out"
#SBATCH --requeue
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters, exiting."
    exit 1
fi

module load MATLAB/2022a
echo "Starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME"
t=$SLURM_ARRAY_TASK_ID
echo $t

matlab -nodisplay -r "hemo_correct_id "$1" "$t" ;exit";