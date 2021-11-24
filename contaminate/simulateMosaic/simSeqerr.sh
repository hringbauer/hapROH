#!/bin/bash
#$ -S /bin/bash #defines bash as the shell for execution
#$ -N hapCon_seqerr #Name of the command that will be listed in the queue
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

errs=(0.001 0.0016681 0.00278256 0.00464159 0.00774264 0.0129155 0.02154435 0.03593814 0.05994843 0.1)
err=${errs[$i]}

#python3 simXwithRC.py -n 100 --cov 0.5 --con 0.1 --err $err --eref 1e-3 --prefix err$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXseqErr/ --hetero
#python3 simXwithRC.py -n 100 --cov 0.1 --con 0.1 --err $err --eref 1e-3 --prefix err$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXseqErr/ --hetero

python3 simXwithRC.py -n 100 --cov 0.05 --con 0.1 --err $err --eref 1e-3 --prefix err$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXseqErr/ --hetero
python3 simXwithRC.py -n 100 --cov 1.0 --con 0.1 --err $err --eref 1e-3 --prefix err$SGE_TASK_ID -b /mnt/archgen/users/yilei/tools/hapROH/simulated/1000G_Mosaic/TSI/maleXseqErr/ --hetero

