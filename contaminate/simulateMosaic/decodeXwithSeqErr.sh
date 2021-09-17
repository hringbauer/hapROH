#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapCon #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
#$ -m e #send an email at the end of the job
#$ -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=25G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 9:10:1

python3 decodeXwithSeqErr.py --err $SGE_TASK_ID --cov 0.1
python3 decodeXwithSeqErr.py --err $SGE_TASK_ID --cov 0.5