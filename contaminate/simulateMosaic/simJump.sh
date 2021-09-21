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
#$ -t 1:10:1


i=$SGE_TASK_ID
i=$(($i-1))

jumps=(100 129 167 215 278 359 464 599 774 1000)
jump=${jumps[$i]}

python3 simXwithRC.py -n 100 --cov 0.5 --con 0.075 --err 1e-2 --eref 1e-3 --prefix jump$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXjumpErr/
python3 simXwithRC.py -n 100 --cov 0.1 --con 0.075 --err 1e-2 --eref 1e-3 --prefix jump$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXjumpErr/
