#!/bin/sh
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=1:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --output="./outs/separate_channel/%j.out"


module load miniconda
conda activate DEEPLABCUT
python 	separate_channel_Hao_v14_auto_limit_64col.py $1 $2
