#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N decodeAuto #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
#$ -m e #send an email at the end of the job
#$ -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=10G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:72:1

id=$SGE_TASK_ID
id=$(($id-1))

cons=(0.0 0.05 0.1)
con_index=$(($id/24))
con=${cons[$con_index]}
id=$(($id-$con_index*24))

blocks=(1 2 4 8)
block_index=$(($id/6))
nblocks=${blocks[$block_index]}

cov_index=$(($id - 6*$block_index))
coverages=(0.05 0.1 0.5 1.0 2.0 5.0)
cov=${coverages[$cov_index]}

echo nblocks$nblocks
echo coverage$cov
echo contamination$con

python3 decodeAuto.py --cov $cov --con $con --nblock $nblocks
