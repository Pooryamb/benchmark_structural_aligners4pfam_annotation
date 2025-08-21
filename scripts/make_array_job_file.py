import os
import argparse


parser= argparse.ArgumentParser("""This script takes the path to a .sh file and also the time needed to run 
the command and generates a job file that can be submitted to a machine whose jobs are moderated with SLURM""")


parser.add_argument("--input_sh_path", type=str, required=True, help="Path to the .sh file. Each line should be a separate command")
parser.add_argument("--time", type=str, required=True, help="The time needed for the job to run. It must be in HH:MM:SS format")
parser.add_argument("--search_category", type=str, required=False, default="sample_pf", help="The time needed for the job to run. It must be in HH:MM:SS format")

args = parser.parse_args()

job_file_name = os.path.basename(args.input_sh_path)
search_category = args.search_category

job_name, _ = os.path.splitext(job_file_name)

commands_num = len([x for x in open(args.input_sh_path).readlines() if (x.strip())])

job_content = """#!/bin/bash
#SBATCH --time={time}
#SBATCH --output=tmp/logs/search/{search_category}/{job_name}_%a.out
#SBATCH --array=1-{batch_num}

module load StdEnv scipy-stack hmmer
source ~/iprs/bin/activate

touch tmp/timestamps/{search_category}/{job_name}_started_B${{SLURM_ARRAY_TASK_ID}}.txt

command=$(sed -n "$((SLURM_ARRAY_TASK_ID + 0))p" {commands_path})
eval "$command"

touch tmp/timestamps/{search_category}/{job_name}_ended_B${{SLURM_ARRAY_TASK_ID}}.txt
"""

with open(f"tmp/jobs/{job_name}_slurm_job.sh", 'w') as j_file:
     j_file.write(job_content.format(time=args.time, job_name=job_name, batch_num=commands_num, commands_path=args.input_sh_path, search_category=search_category))
