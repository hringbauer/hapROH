#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N HapROH #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -pe smp 4 #needs 8 CPU cores
#$ -l h_vmem=100G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID
# -tc 20

for SAMPLE in 'NG10'; do
    /home/xiaowen_jia/anaconda3/bin/python3 haproh.py $SAMPLE
done