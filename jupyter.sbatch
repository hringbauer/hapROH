#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --partition=jnovembre #jnovembre  #broadwl
#SBATCH --account=pi-jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=2g
#SBATCH --job-name=HAPSBURG_Parallel
#SBATCH --output=nb-log-%J.out
#SBATCH --error=nb-log-%J.err

# module load python/3.7.0
module load python
module load plink
module load bcftools

ip=$(/sbin/ip route get 8.8.8.8 | awk '{print $NF;exit}')

echo "http://"$ip":"$port"/"
jupyter lab --no-browser --ip=$ip --port=9279 --log-level='ERROR'
