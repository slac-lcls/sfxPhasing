#!/usr/bin/env bash
#SBATCH --image=docker:slaclcls/autosfx:latest
#SBATCH --nodes=20
#SBATCH --tasks-per-node=10
#SBATCH --constraint=haswell
#SBATCH --qos=realtime
#SBATCH --gres=craynetwork:4
#SBATCH --job-name=sfx
#SBATCH --account=lcls
#SBATCH --time=05:00:00
#run the application:
module load python
ulimit -u 8192
ulimit -n 60000
python MR_batch.py -rfl input.mtz -pdb input.pdb -seq input.fasta -n 4 --shifter 1
#python MR_batch.py -rfl lyso_gd_test2_xscale.mtz -pdb 1LSG_A.pdb -seq hewl.seq -n 2 --shifter 1
