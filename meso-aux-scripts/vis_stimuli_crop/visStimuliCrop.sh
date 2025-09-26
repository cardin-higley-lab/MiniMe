#!/bin/sh
#SBATCH --ntasks=1 --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --partition="scavenge"
#SBATCH --requeue
#SBATCH --output="./outs/visCrop/%j.out"


module load GCCcore/12.2.0 cuDNN/8.7.0.84-CUDA-11.8.0 miniconda
conda activate DEEPLABCUT
python vis_stimuli_on_crop_06242025.py
