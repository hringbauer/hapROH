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

refErrs=(0.0001 0.00016681 0.00027826 0.00046416 0.00077426 0.00129155 0.00215443 0.00359381 0.00599484 0.01)
refErr=${refErrs[$i]}

python3 simXwithRC.py -n 100 --cov 1.0 --con 0.1 --err 1e-2 --eref $refErr --prefix ref$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXrefErr/ --hetero
python3 simXwithRC.py -n 100 --cov 0.5 --con 0.1 --err 1e-2 --eref $refErr --prefix ref$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXrefErr/ --hetero
python3 simXwithRC.py -n 100 --cov 0.1 --con 0.1 --err 1e-2 --eref $refErr --prefix ref$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXrefErr/ --hetero
python3 simXwithRC.py -n 100 --cov 0.05 --con 0.1 --err 1e-2 --eref $refErr --prefix ref$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXrefErr/ --hetero
