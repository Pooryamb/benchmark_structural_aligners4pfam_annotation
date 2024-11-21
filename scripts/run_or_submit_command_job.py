import os
import argparse



parser = argparse.ArgumentParser(description="""For this benchmarking, I had to run the commands on a computational
server moderated by SLURM. This script has been designed to make it possible to run the same commands on both servers
moderated by a job scheduler or servers that don't have one and work with the intact data (that has not been cut into chunks""")

parser.add_argument("--command", type=str, required=True, help="""The command to be executed.
If a path ends with #, such as "path/data#", this script changes them to "path/B${i}/data,
where i shows the batch number. If a path ends with ? mark, such as "path/file.txt?", it will
be changed to "path/file${i}.txt", where i shows the batch number.""")

parser.add_argument("--job_commands", type=str, help="""The commands that must be run on the server 
before executing the main commands. This is only for the cases that the command is going to be executed on
a server moderated by a job scheduler. This can be used to input loading commands, etc.""")

parser.add_argument("--how2run", required=True, help="It shows if the user wants to run commands using job_scheduler/gnu/none")

parser.add_argument("--time", type=str, help="The run time for the jobs. It must be in HH:MM:SS format. Alternatively, you can provide the server commands all together")

parser.add_argument("--add_default_header", action="store_true", help="""It shows if the default heading should be used.
If this option is used, it requires specifying --time as well""")

parser.add_argument("--batches", type=int, help="It shows the number of batches")

parser.add_argument("--job_name", type=str, help="""This will be used for naming the output file""")

args = parser.parse_args()

if args.how2run == "job_scheduler":
    if (not(args.job_commands) and not(args.add_default_header)):
        raise Exception("""If your command is going to be run on a server, you have to specify the server commands as well. 
Alternatively, you can use default header""")
    if args.add_default_header and (not(args.time) or not(args.job_name)):
        raise Exception("You must speicfy the run time and job name even when you use the default header")

default_header="""#!/bin/bash
#SBATCH --time={time}
#SBATCH --output=./tmp/{job_name}_log_%a.out
#SBATCH --array=1-{batches}

module load CCEnv StdEnv scipy-stack
source ~/iprs/bin/activate
i=${{SLURM_ARRAY_TASK_ID}}

"""

def convert2par(command):
    slurm_parts = []
    parts = command.split()

    for part in parts:
        if part.endswith("#") or part.endswith("?"):
            directory, base_name = os.path.split(part)
            if base_name.endswith("#"):
                bs_no_sign = base_name.replace("#", "")
                slurm_parts.append(os.path.join(directory, "B${i}", bs_no_sign))
            elif base_name.endswith("?"):
                bs_no_sign = base_name.replace("?", "")
                new_file_name = os.path.splitext(bs_no_sign)[0] + "_B${i}" + os.path.splitext(bs_no_sign)[1]
                slurm_parts.append(os.path.join(directory, new_file_name))
        else:
            slurm_parts.append(part)
    return " ".join(slurm_parts)


main_command = args.command

if args.how2run == "none":
    os.system(main_command.replace("#", "").replace("?", ""))

else:
    commands = args.command.split(";")
    parallel_command = ""
    for command in commands:
        parallel_command = parallel_command + convert2par(command) + ' ; '


    if args.how2run == "job_scheduler":
        with open(f"./tmp/tmp_{args.job_name}.sh", 'w') as job_file:
            if args.add_default_header:
                job_file.write(default_header.format(time=args.time, job_name=args.job_name, batches=args.batches))
            else:
                job_file.write(args.job_commands + '\n')

            job_file.write(parallel_command + "\n\n")
        os.system(f"sbatch ./tmp/tmp_{args.job_name}.sh")
    elif args.how2run == "gnu":
        os.system("parallel --dry-run 'for i in {}; do " + parallel_command +  "done' ::: {1.." + str(args.batches) + "}")
