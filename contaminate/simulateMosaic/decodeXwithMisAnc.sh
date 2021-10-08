#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapCon_infer #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=25G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:54:1

i=$SGE_TASK_ID
i=$(($i-1))
coverages=(0.05 0.1 0.5 1.0 2.0 5.0)
cov_index=$(($i/9))
cov=${coverages[$cov_index]}
i=$(($i-9*$cov_index))

pop1=(CEU YRI CHB)
pop2=(CEU YRI CHB)
pop1_index=$(($i/3))
pop1=${pop1[$pop1_index]}
pop2_index=$(($i-3*$pop1_index))
pop2=${pop2[$pop2_index]}

echo coverage$cov
echo pop1$pop1
echo pop2$pop2

python3 decodeXwithMisAnc.py --cov $cov --con 0.1 --err 1e-2 --eref 1e-3 --conpop1 $pop1 --conpop2 $pop2
