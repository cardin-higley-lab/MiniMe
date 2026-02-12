#!/bin/sh
#SBATCH --job-name=1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=64G
#SBATCH --time=12:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --cpus-per-task 4
#SBATCH --output="./outs/post_hemo/%j.out"

module load MATLAB/2023a


if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, exiting."
    exit 2
fi

echo "Starting job on $HOSTNAME"

echo $1
echo $2
cd /home/$USER/project/cluster_code/; matlab -nodisplay -r "post_hemo('"$1"','"$2"'); exit";


