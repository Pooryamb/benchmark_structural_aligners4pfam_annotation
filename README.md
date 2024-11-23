# Introduction

This repository contains the scripts used for comparing Foldseek and Reseek 
for the structure-based domain annotation task.
You require the following libraries for this benchmarking:
* Pandas (version 2.2.2)
* Polars (version 1.0.0)


To start, make sure that both reseek, foldseek, and mmseqs have been added to 
your $PATH variable. 


The tar file of Pfam structures has to be in data/pfam_cifs.tar file. 
All the databases will be stored in the data/dbs directory.
The following dbs should be present in data/dbs when the process starts:
* data/dbs/pfam_fs_cut/pfam :This database can be made by InterProSDB repository. Please note that for making 
this database, tsv files are converted to the Foldseek database files. The tsv files used for this operation must be inside the copied database.
* data/dbs/pfam_fl/pfam_fl : This includes the full-length database of the proteins used for making the Pfam database. This database also can be 
made using the InterProSDB repository. It can be found inside the tmp directory.



Using the following commands, make the database from cif files
```
mkdir -p data/dbs/pfam_cif_cut
foldseek createdb data/pfam_cifs.tar ./data/dbs/pfam_cif_cut/pfam
```

Cluster the Pfam database and select the cluster rep with the highest average plddt

```
foldseek convert2fasta ./data/dbs/pfam_cif_cut/pfam ./data/dbs/pfam_cif_cut/pfam.fasta
mkdir -p tmp/
mmseqs easy-cluster ./data/dbs/pfam_cif_cut/pfam.fasta ./data/dbs/pfam_cif_cut/pfam_clust tmp/pfam_clust --min-seq-id 0.4 -c 0.8 --cov-mode 0 -s 9
python scripts/extract_plddt.py --input ./data/pfam_cifs.tar --output ./data/pfam_plddt.json # Extracts plddt from structures
python scripts/calc_plddt_average_from_json.py # It converts the json of plddts to a tsv file containing the average plddts
python ./scripts/select_clusterrep_by_plddt_and_size.py #It selects the one with highest pLDDT as cluster representative. It also removes domains that are shorter than 10 amino acids
# First, we select the representatives from cif_cut database
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_clust_reps.tsv" SRC_DB="./data/dbs/pfam_cif_cut/pfam" OUT_DB="./data/dbs/pfam_cif_cut_clust/pfam"
# Then, we select the representatives from fs_cut database
python ./scripts/select_subtsv.py --input_basename "./data/dbs/pfam_fs_cut/pfam" --list2sel ./data/pfam_clust_reps.tsv --output_basename "./data/dbs/pfam_fs_cut_clust/pfam" --remove_cif_extension True
./scripts/convert_tsv2fsdb.sh ./data/dbs/pfam_fs_cut_clust/pfam ./tmp/make_fscut_db_log
```

Select a sample of Pfam to query against PfamSDB

```
# This will select a sample of cif_cut database
python scripts/select_random_sample_from_pfam.py --number=10000  # Selects 10000 random seeds to be queried against clustered Pfam
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_sample_reps.tsv" SRC_DB="./data/dbs/pfam_cif_cut_clust/pfam" \
OUT_DB="./data/dbs/pfam_cif_cut_sample/pfam"

# This will select a sample of the fs_cut database
python ./scripts/select_subtsv.py --input_basename "./data/dbs/pfam_fs_cut_clust/pfam" --list2sel ./data/pfam_sample_reps.tsv --output_basename "./data/dbs/pfam_fs_cut_sample/pfam" --remove_cif_extension True
./scripts/convert_tsv2fsdb.sh ./data/dbs/pfam_fs_cut_sample/pfam ./tmp/make_fscut_sample_db_log
```

All-against-all searching takes a long time. I had some limitations regarding the number of hours
that my job could be run on our computational server. Therefore, I cut each database into smaller
pieces and search each chunk separately. Please note that although I could search a chunk against
intact clustered version of Pfam using Foldseek 3Di+AA alignment, I couldn't do it using TM-align
. So, I will also cut the target database for TM-align.

```
CHUNK_NUM=16
total_lines=$(wc -l < ./data/pfam_sample_reps.tsv)
lines_per_batch=$(( (total_lines + CHUNK_NUM - 1) / CHUNK_NUM ))
split -l $lines_per_batch ./data/pfam_sample_reps.tsv ./tmp/sample_reps_
n=1  
for file in ./tmp/sample_reps_*; do  
    mv "$file" ./tmp/sample_reps_b${n}  
    n=$((n + 1))
done

for i in $(seq 1 $CHUNK_NUM); do
    make -f ./scripts/Makefile.make_subdb REP_INFO="./tmp/sample_reps_b${i}" SRC_DB="./data/dbs/pfam_cif_cut_sample/pfam" OUT_DB="./data/dbs/pfam_cif_cut_sample/B${i}/pfam"
done


for i in $(seq 1 $CHUNK_NUM); do
    python ./scripts/select_subtsv.py --input_basename "./data/dbs/pfam_fs_cut_sample/pfam" --list2sel ./tmp/sample_reps_b${i} --output_basename "./data/dbs/pfam_fs_cut_sample/B${i}/pfam" --remove_cif_extension True
    ./scripts/convert_tsv2fsdb.sh ./data/dbs/pfam_fs_cut_sample/B${i}/pfam ./tmp/make_fscut_sample_db_log_B${i}
done
```


The following snippet, converts each fsdb to cal format and then converts the cal to bca format:
```
for file in $(find "./data/dbs/" -iname "*.lookup"); do  
    python scripts/convert_fsdb2cal.py --input_basename ${file/.lookup/}
    reseek -convert ${file/.lookup/.cal} -bca ${file/.lookup/.bca}
done
```

The following part is for searching a sample of domains against the clustered Pfam.
First, it creates a file containing all of the commands. Then, user can run them
either regularly or submit them as a job to a machine that has SLURM job scheduler.
The instructions for running the scripts are found at the end of the code snippet.
I used the following run times for a cluster with 80 threads: 
reseek: 1:00:00, mmseqs, fs_cut and cif_cut: 00:15:00, tm_align: 3:00:00
```
mkdir -p tmp/start_end_time/
mkdir -p data/alis/sample_vs_pfam
mkdir -p ./tmp/alidbs/
mkdir -p ./tmp/job_logs/


sh_path=./tmp/rs_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "reseek -search ./data/dbs/pfam_cif_cut_sample/B${i}/pfam.bca -db ./data/dbs/pfam_cif_cut_clust/pfam.bca -output data/alis/sample_vs_pfam/reseek_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -verysensitive -evalue 1e99" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/fscut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/dbs/pfam_fs_cut_sample/B${i}/pfam ./data/dbs/pfam_fs_cut_clust/pfam data/alis/sample_vs_pfam/fs_cut_B${i}.tsv tmp/pfam_fs_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


sh_path=./tmp/cifcut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/dbs/pfam_cif_cut_sample/B${i}/pfam ./data/dbs/pfam_cif_cut_clust/pfam data/alis/sample_vs_pfam/cif_cut_B${i}.tsv tmp/pfam_cif_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/mm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "mmseqs easy-search --prefilter-mode 2 -e inf ./data/dbs/pfam_cif_cut_sample/B${i}/pfam.fasta ./data/dbs/pfam_cif_cut_clust/pfam.fasta data/alis/sample_vs_pfam/mm_B${i}.tsv tmp/pfam_mm_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/tm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do  
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/dbs/pfam_cif_cut_sample/B${i}/pfam ./data/dbs/pfam_cif_cut_clust/pfam data/alis/sample_vs_pfam/tm_B${i}.tsv tmp/pfam_tm_B${i} --alignment-type 1 --tmscore-threshold 0.0 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore" >> ${sh_path}
done  
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

#To run on a single machine:
#bash $sh_path
#To run on a machine operated by JOB scheduler:
#python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "3:00:00"; sbatch ./tmp/${sh_path/.sh/_slurm_job.sh}
#

```


Now, we find sensitivity by finding the number of TPs before the first FP.
We want to find labels both at family and clan level. So, we need to 
download the file containing the clan data, and make a fake clan id for 
those that don't belong to a clan. For such cases, we use the family id
as the clan id.

```
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.clans.tsv.gz -P ./tmp # At the time of working on this project, version 37 was downloaded. Feel free to change the version if needed
gunzip ./tmp/Pfam-A.clans.tsv.gz
python scripts/process_pf_clans_info.py # This will make a tsv file whose first column is Pfam id and second column is clan id
```

The next scripts will be used for finding the sensitivity based on the number of FPs before the first TP.
```
mkdir -p ./data/first_label_occ

file_paths=$(find ./data/alis/sample_vs_pfam/ -type f -name "*.tsv")
echo "$file_paths" | parallel "python scripts/find_nonred_labels.py --input {}"

```


```
# The next lines will make the clustered version of pfam_fl, which I skip for now
# python scripts/find_pfam_fl_clust_reps.py #This will make a file containing the name of proteins that were used for making PfamSDB
# make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_fl_clust_reps.tsv" SRC_DB="./data/dbs/pfam_fl/pfam_fl" OUT_DB="./data/dbs/pfam_fl_clust/pfam_fl"
```

## Searching a subset of Pfam_fls to search against clustered PfamSDB
```
shuf --random-source=./data/pfam_fl_clust_reps.tsv -n 1000 ./data/pfam_fl_clust_reps.tsv -o ./data/pfam_fl_sample1.tsv
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_fl_sample1.tsv" SRC_DB="./data/dbs/pfam_fl/pfam_fl" OUT_DB="./data/dbs/pfam_fl_sample1/pfam_fl"
python scripts/convert_tsv2cal.py  --input_basename ./data/dbs/pfam_fs_cut_clust/pfam
reseek -convert ./data/dbs/pfam_fs_cut_clust/pfam.cal -bca ./data/dbs/pfam_fs_cut_clust/pfam.bca
python scripts/convert_fsdb2cal.py --input_basename ./data/dbs/pfam_fl_clust/pfam_fl


python scripts/convert_fsdb2cal.py --input_basename ./data/dbs/pfam_fl_sample1/pfam_fl
reseek -convert ./data/dbs/pfam_fl_sample1/pfam_fl.cal -bca ./data/dbs/pfam_fl_sample1/pfam_fl.bca
foldseek easy-search --exhaustive-search 1 -e inf ./data/dbs/pfam_fl_sample1/pfam_fl ./data/dbs/pfam_fs_cut_clust/pfam ./data/sample1_pfam_foldseek_fscut.tsv tmp/pfam_cifcut
reseek -search ./data/dbs/pfam_fl_sample1/pfam_fl.bca -db ./data/dbs/pfam_fs_cut_clust/pfam.bca \
-output ./data/sample1_pfam_reseek.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq \
-verysensitive -evalue 1e99
```


## Investigating exhaustive pairwise alignments between members of a domain
```
python ./scripts/group_pfamclust_by_family.py #This will make a directory for each family in "data/dbs/pfam_fs_cut_clust/grp_by_family" directory inside each directory, tsv files of the database will be stored
mkdir -p tmp/tsv2fsdb_pfam_logs

ls data/dbs/pfam_fs_cut_clust/grp_by_family/PF*/PF*_ca.tsv | xargs -P 20 -I {} sh -c '
    base_name=$(basename {} _ca.tsv);
    ./scripts/convert_tsv2fsdb.sh ${1%_ca.tsv} tmp/tsv2fsdb_pfam_logs/${base_name}
' -- {}
# Log files inside tmp/tsv2fsdb_pfam_logs/ can be investigated to make sure that there is no empty file over there. Use the following command in this regard. The files with
# non-empty logs will go into the filepath specified by output.
python scripts/check_all_files_are_empty.py --dir "./tmp/tsv2fsdb_pfam_logs/" --output "./tmp/pfam_db_cpy_log.txt"

# Now, it is the time to cut the full length database:
./scripts/convert_fsdb2fasta.sh ./data/dbs/pfam_fl_clust/pfam_fl
python ./scripts/convert_fasta2tsv.py --input_basename "./data/dbs/pfam_fl_clust/pfam_fl"
python ./scripts/group_pfam_fl_clust_by_family.py #This will break the FL database so that FL proteins of each Pfam would be in a different folder
mkdir -p tmp/tsv2fsdb_pfamfl_logs
ls data/dbs/pfam_fl_clust/grp_by_family/PF*/PF*_ca.tsv | xargs -P 20 -I {} sh -c '
    base_name=$(basename {} _ca.tsv);
    ./scripts/convert_tsv2fsdb.sh ${1%_ca.tsv} tmp/tsv2fsdb_pfamfl_logs/${base_name}
' -- {}
python scripts/check_all_files_are_empty.py --dir "./tmp/tsv2fsdb_pfamfl_logs/" --output "./tmp/pfam_fl_db_cpy_log.txt"

# Now, we start the exhaustive search
all_pfs=$(ls data/dbs/pfam_fs_cut_clust/grp_by_family/PF*/PF*_ca)
mkdir -p tmp/pf_pf_exh_alis/
mkdir -p tmp/fs_exh_search_logs/

for pf_path in ${all_pfs}
do
    pf_basename=${pf_path/_ca/}
    pf_name=$(basename $pf_basename)
    foldseek easy-search --exhaustive-search 1 -e inf --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddtfull,qaln,taln $pf_basename $pf_basename tmp/pf_pf_exh_alis/${pf_name}.tsv tmp/tmp_pf_pf_exhaustive -v 1 > tmp/fs_exh_search_logs/${pf_name} 2>&1
done

# I stopped at this point. If I decided to continue, I must first make the databases from tsv files and then conduct the comprehensive search
# At this point, we have made one database per Pfam family. Next, we will conduct an exhaustive members against members alignment

```






















The following snippet is for cutting the clustered version of Pfam. It was needed only for troubleshooting TM-align failure
TM-align failed with intact target database for an unknown reason.

```
T_CHUNK_NUM=4
total_lines=$(wc -l < ./data/pfam_clust_reps.tsv)
lines_per_batch=$(( (total_lines + T_CHUNK_NUM - 1) / T_CHUNK_NUM ))                                                                                        split -l $lines_per_batch ./data/pfam_clust_reps.tsv ./tmp/clust_reps_
n=1
for file in ./tmp/clust_reps_*; do
    mv "$file" ./tmp/clust_reps_b${n}
    n=$((n + 1))
done

for i in $(seq 1 $T_CHUNK_NUM); do                                                                                                                              make -f ./scripts/Makefile.make_subdb REP_INFO="./tmp/clust_reps_b${i}" SRC_DB="./data/dbs/pfam_cif_cut_clust/pfam" OUT_DB="./data/dbs/pfam_cif_cut_clu$done                                                                                                                                                        ```
