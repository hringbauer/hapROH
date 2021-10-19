#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N misContam #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -q archgen.q #queue
# -m e #send an email at the end of the job
# -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=40G #request 4Gb of memory
#$ -V # load personal profile
#$ -o $JOB_NAME.o.$JOB_ID.$TASK_ID
#$ -t 1:24:1

id=$SGE_TASK_ID
id=$(($id-1))

conSlist=(0 0.05 0.1 0.15)
conSindex=$(($id/24))
conS=${conSlist[$conSindex]}
id=$(($id-24*$conSindex))

#conIlist=(0.0 0.05 0.1 0.15 0.2 0.25)
conIlist=(0.0 0.01 0.02 0.05 0.075 0.1)
conIindex=$(($id/4))
conI=${conIlist[$conIindex]}
id=$(($id-4*$conIindex))

coverage=(0.5 1.0 2.0 5.0)
cov=${coverage[$id]}

echo conS$conS
echo conI$conI
echo coverage$cov

python3 misContam.py --cov $cov --conI $conI --conS $conS --nblock 5

