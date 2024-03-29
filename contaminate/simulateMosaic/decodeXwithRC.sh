#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapCon #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
#$ -m e #send an email at the end of the job
#$ -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=100G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 45:45:1

i=$SGE_TASK_ID
i=$(($i-1))
coverages=(0.05 0.1 0.5 1.0 2.0 5.0 10 20)
cov_index=$(($i/6))
cov=${coverages[$cov_index]}

cons=(0.0 0.05 0.1 0.15 0.2 0.25)
con_index=$(($i-6*$cov_index))
con=${cons[$con_index]}

echo cov$cov
echo con$con

python3 decodeXwithRC.py --cov $cov --con $con --err 1e-2 --eref 1e-3
