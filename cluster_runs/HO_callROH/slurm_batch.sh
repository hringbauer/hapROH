#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="HO_HAPSBURG"
#SBATCH --time=2:00:00
#SBATCH --mem=8G
# #SBATCH --mail-user=hringbauer@uchicago.edu
# #SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=broadwl
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=500-598 #598
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

# Execute the following tasks
module load python/3.7.0
python3 run_HO_longROH.py $SLURM_ARRAY_TASK_ID 
