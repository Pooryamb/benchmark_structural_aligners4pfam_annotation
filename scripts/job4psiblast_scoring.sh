#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --output=tmp/logs/psi_blast_scoring.txt


module load CCEnv StdEnv python/3.13 scipy-stack arrow

source ~/pooryamb/str_align_benchmark/bin/activate

touch tmp/timestamps/split_pf_seq/rescoring_started

python scripts/calc_psiblast_raw_score.py

touch tmp/timestamps/split_pf_seq/rescoring_ended
