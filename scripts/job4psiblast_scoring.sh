#!/bin/bash
#SBATCH --time=00:40:00
#SBATCH --output=tmp/logs/psi_blast_scoring_%a.txt
#SBATCH --array=1-16

module load StdEnv python/3.13 scipy-stack arrow

source ~/pooryamb/str_align_benchmark/bin/activate

touch tmp/timestamps/split_pf_seq/rescoring_started_reseek_${SLURM_ARRAY_TASK_ID}

python scripts/calc_psiblast_raw_score.py --alis ./tmp/alis/split_pf_seq/reseek_sens_B${SLURM_ARRAY_TASK_ID}.tsv

touch tmp/timestamps/split_pf_seq/rescoring_ended_reseek_${SLURM_ARRAY_TASK_ID}
