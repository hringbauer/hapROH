#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --job-name="PREP_REF_HDF5"
#SBATCH --time=6:00:00
#SBATCH --mem=40G
# #SBATCH --mail-user=hringbauer@uchicago.edu
# #SBATCH --mail-type=ALL
#SBATCH --no-requeue
#SBATCH --partition=broadwl
#SBATCH --export=NONE
#SBATCH --output=./logs/%A_%a.out
#SBATCH --error=./logs/%A_%a.err
#SBATCH --array=1-20
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

# Execute the following tasks
module load python
module load bcftools

python3 prep_ch.py $SLURM_ARRAY_TASK_ID 