#!/bin/sh
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=1:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --output="./outs/separate_channel/%j.out"


module load miniconda

# Initialize conda for batch shell
source $(conda info --base)/etc/profile.d/conda.sh

# Check if environment exists
if conda env list | grep -q "analysis"; then
    echo "Environment analysis already exists."
else
    echo "Creating analysis environment..."
    conda env create -f analysis.yml
fi

conda activate analysis

DATA_FOLDER=/vast/palmer/scratch/higley/hd362/HD_Mouse_Training

python 	separate_channel_Hao_v14_auto_limit_64col.py $1 $2 $DATA_FOLDER
