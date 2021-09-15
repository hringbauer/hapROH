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
#$ -t 1:6:1

#iids=("ORC008" "SUC005" "SUA003" "S1252" "S1253" "SEC002" "ORC002" "SUC009" "ORC003" "MA110" "ISB001"
#"MA78" "MA73" "SUA001" "ORC007" "MA89" "ISC001" "SEC001" "SUA002" "SUC006" "SUC007" "SUC003" "ORC004"
#"S1250" "MA138" "MA81" "COR002" "MA100" "MA112" "S1249" "ORC006")


i=$SGE_TASK_ID
i=$(($i-1))
coverages=(0.05 0.1 0.5 1.0 2.0 5.0)
cov=${coverages[$i]}

python3 simXwithRC.py -n 100 --cov $cov --con 0.0 --err 1e-3 --eref 1e-3
python3 simXwithRC.py -n 100 --cov $cov --con 0.05 --err 1e-3 --eref 1e-3 
python3 simXwithRC.py -n 100 --cov $cov --con 0.1 --err 1e-3 --eref 1e-3 