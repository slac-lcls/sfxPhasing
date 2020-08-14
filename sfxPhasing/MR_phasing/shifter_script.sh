#!/usr/bin/env bash
#SBATCH --image=docker:registry.services.nersc.gov/annagian/sfx_6
#SBATCH --nodes=20
#SBATCH --tasks-per-node=10
#SBATCH --constraint=haswell
#SBATCH --qos=realtime
#SBATCH --gres=craynetwork:4
#SBATCH --job-name=sfx
#SBATCH --account=lcls
#SBATCH --mail-user=agiannakou@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --time=05:00:00
#run the application:
module load python
ulimit -u 8192
ulimit -n 60000
python MR_batch.py -rfl mt2_push1d2_xscale_ccp4if_Rfree.mtz -pdb 5T1A_A_sculpt.pdb -seq mt2_cryst.fasta -n 4 --shifter 1
#python MR_batch.py -rfl lyso_gd_test2_xscale.mtz -pdb 1LSG_A.pdb -seq hewl.seq -n 2 --shifter 1
