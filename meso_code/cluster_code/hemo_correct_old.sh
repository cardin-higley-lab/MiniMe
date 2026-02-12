#!/bin/sh
#SBATCH --job-name=1-1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=40G
#SBATCH --time=0:20:00
#SBATCH --partition="scavenge"
#SBATCH --array=1-64

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, exiting."
    exit 2
fi

module load MATLAB/2023a
echo "Starting job $SLURM_ARRAY_TASK_ID on $HOSTNAME"
t=$SLURM_ARRAY_TASK_ID
echo $t

echo $SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_JOB_ID

matlab -nodisplay -r "hemo_correct "$1" "$t" ;exit";

if [ "$t" -eq 64 ]; then 
    cd /home/$USER/project/cluster_code/ ; sbatch --dependency=afterok:$SLURM_ARRAY_JOB_ID post_hemo.sh $1 $2
fi


