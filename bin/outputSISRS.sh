#!/bin/bash
#SBATCH --job-name="WGS_Pairwise"
#SBATCH --time=48:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node 
#SBATCH --mail-user="Robert.Literman@fda.hhs.goc"
#SBATCH --mail-type=END,FAIL
# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
cd $SLURM_SUBMIT_DIR

python Output_SISRS.py WGS_Pairwise
