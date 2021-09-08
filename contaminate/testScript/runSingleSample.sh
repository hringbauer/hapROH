#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapCon #Name of the command that will be listed in the queue
#$ -cwd #change to current directory
#$ -j y #join error and standard output in one file, no error file will be written
#$ -o output #standard output file or directory (joined with error because of -j y)
#$ -q archgen.q #queue
#$ -m e #send an email at the end of the job
#$ -M yilei_huang@eva.mpg.de #send email to this address
# -pe make 2 #needs 8 CPU cores
#$ -l h_vmem=20G #request 4Gb of memory
#$ -V # load personal profile

#iids=("ORC008" "SUC005" "SUA003" "S1252" "S1253" "SEC002" "ORC002" "SUC009" "ORC003" "MA110" "ISB001"
#"MA78" "MA73" "SUA001" "ORC007" "MA89" "ISC001" "SEC001" "SUA002" "SUC006" "SUC007" "SUC003" "ORC004"
#"S1250" "MA138" "MA81" "COR002" "MA100" "MA112" "S1249" "ORC006")


python3 runSingleSample.py -i SUC005
