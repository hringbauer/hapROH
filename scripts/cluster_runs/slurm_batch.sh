#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="HO_HAPSBURG"
#SBATCH --time=1:20:00
#SBATCH --mem=7G
# #SBATCH --mail-user=hringbauer@uchicago.edu
# #SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=broadwl
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=131-599
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

# Execute the following tasks
module load python/3.7.0
python3 run_individual_HO.py $SLURM_ARRAY_TASK_ID 
