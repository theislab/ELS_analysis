#!/bin/bash
#SBATCH -J ELS_DIFFXPY
#SBATCH -o /home/icb/carlo.dedonno/projects/ELS/code/outs/ELS_DIFFXPY.out
#SBATCH -e /home/icb/carlo.dedonno/projects/ELS/code/outs/ELS_DIFFXPY.err
#SBATCH -p cpu_p
#SBATCH -c 32
#SBATCH -t 2-00:00:00
#SBATCH --mem=50G
#SBATCH --nice=1000

echo $SLURM_NODENAME
source /home/icb/carlo.dedonno/anaconda3/bin/activate sc
python /home/icb/carlo.dedonno/projects/ELS/code/diffxpy_standard.py