#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="HAPSBURG"
#SBATCH --time=1:00:00
#SBATCH --mem=12G
# #SBATCH --mail-user=hringbauer@uchicago.edu
# #SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=broadwl
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=500-999   # max 1094
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

# Execute the following tasks
module load python/3.7.0
python3 run_individual.py $SLURM_ARRAY_TASK_ID 
