#!/bin/bash

#SBATCH --time=20:00:00
#SBATCH --partition=jnovembre
#SBATCH --account=pi-jnovembre
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=3g
#SBATCH --job-name=HAPSBURG_Parallel
#SBATCH --output=nb-log-%J.out
#SBATCH --error=nb-log-%J.err

module load python/3.7.0

ip=$(/sbin/ip route get 8.8.8.8 | awk '{print $NF;exit}')
port=$((10000+ $RANDOM % 20000))

echo "http://"$ip":"$port"/"
jupyter lab --no-browser --ip=$ip --port=$port --log-level='ERROR'
