# Introduction

This repository contains the scripts used for comparing Foldseek and Reseek 
for the structure-based domain annotation task.
You require the following libraries for this benchmarking:
* Pandas (version 2.2.2)
* Polars (version 1.0.0)


To start, make sure that both reseek, foldseek, and mmseqs have been added to 
your $PATH variable. 


The tar file of Pfam structures has to be in data/raw/pfam_cifs.tar file. 
All the databases will be stored in the data/raw/dbs directory.
The following db should be present in data/raw/dbs when the process starts:
* data/raw/dbs/pfam_fs_cut/pfam :This database can be made by InterProSDB repository. Please note that for making 
this database, tsv files are converted to the Foldseek database files. The tsv files used for this operation must be inside the copied database.


Using the following commands, make the database from cif files
```
mkdir -p data/raw/dbs/pfam_cif_cut
foldseek createdb data/raw/pfam_cifs.tar data/raw/dbs/pfam_cif_cut/pfam
```

Cluster the Pfam database and select the cluster rep with the highest average plddt

```
foldseek convert2fasta ./data/raw/dbs/pfam_cif_cut/pfam ./data/raw/dbs/pfam_cif_cut/pfam.fasta
mkdir -p tmp/fstmp/ tmp/alidb tmp/alis/sample_pf tmp/alis/pf_pf tmp/alis/pfam_clust tmp/logs/misc/
mmseqs easy-cluster ./data/raw/dbs/pfam_cif_cut/pfam.fasta tmp/alis/pfam_clust/pfam_clust tmp/fstmp/clust --min-seq-id 0.4 -c 0.8 --cov-mode 0 -s 9
mkdir -p ./data/processed/
python scripts/extract_plddt.py --input ./data/raw/pfam_cifs.tar --output ./data/processed/pfam_plddt.json # Extracts plddt from structures
python scripts/calc_plddt_average_from_json.py # It converts the json of plddts to a tsv file containing the average plddts
python ./scripts/select_clusterrep_by_plddt_and_size.py #It selects the one with highest pLDDT as cluster representative. It also removes domains that are shorter than 10 amino acids
# First, we select the representatives from cif_cut database
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/processed/pfam_clust_reps.tsv" SRC_DB="./data/raw/dbs/pfam_cif_cut/pfam" OUT_DB="./data/raw/dbs/pfam_cif_cut_clust/pfam"
# Then, we select the representatives from fs_cut database
python ./scripts/select_subtsv.py --input_basename "./data/raw/dbs/pfam_fs_cut/pfam" --list2sel ./data/processed/pfam_clust_reps.tsv --output_basename "./data/raw/dbs/pfam_fs_cut_clust/pfam" --remove_cif_extension True
mkdir -p tmp/logs/makedb/tsv2fsdb
./scripts/convert_tsv2fsdb.sh ./data/raw/dbs/pfam_fs_cut_clust/pfam tmp/logs/makedb/tsv2fsdb/make_fs_cut_clust_log
```

Select a sample of Pfam to query against PfamSDB

```
# This will select a sample of cif_cut database
python scripts/select_random_sample_from_pfam.py --number=10000  # Selects 10000 random seeds to be queried against clustered Pfam
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/processed/pfam_sample_reps.tsv" SRC_DB="./data/raw/dbs/pfam_cif_cut_clust/pfam" \
OUT_DB="./data/raw/dbs/pfam_cif_cut_sample/pfam"

# This will select a sample of the fs_cut database
python ./scripts/select_subtsv.py --input_basename "./data/raw/dbs/pfam_fs_cut_clust/pfam" \
--list2sel ./data/processed/pfam_sample_reps.tsv --output_basename "./data/raw/dbs/pfam_fs_cut_sample/pfam" --remove_cif_extension True
./scripts/convert_tsv2fsdb.sh ./data/raw/dbs/pfam_fs_cut_sample/pfam tmp/logs/makedb/tsv2fsdb/make_fs_cut_sample_log
```

All-against-all searching takes a long time. We had some limitations regarding the number of hours
that the job could be run on our computational server. Therefore, we cut each database into smaller
pieces and search each chunk separately.
```
CHUNK_NUM=16
total_lines=$(wc -l < ./data/processed/pfam_sample_reps.tsv)
lines_per_batch=$(( (total_lines + CHUNK_NUM - 1) / CHUNK_NUM ))
mkdir -p tmp/sample_list
split -l $lines_per_batch ./data/processed/pfam_sample_reps.tsv ./tmp/sample_list/sample_reps_

n=1  
for file in ./tmp/sample_list/sample_reps_*; do
    mv "$file" ./tmp/sample_list/sample_reps_b${n}
    n=$((n + 1))
done

for i in $(seq 1 $CHUNK_NUM); do
    make -f ./scripts/Makefile.make_subdb REP_INFO="./tmp/sample_list/sample_reps_b${i}" \
SRC_DB="./data/raw/dbs/pfam_cif_cut_sample/pfam" OUT_DB="./data/raw/dbs/pfam_cif_cut_sample/B${i}/pfam"
done


for i in $(seq 1 $CHUNK_NUM); do
    python ./scripts/select_subtsv.py --input_basename "./data/raw/dbs/pfam_fs_cut_sample/pfam" \
--list2sel ./tmp/sample_list/sample_reps_b${i} --output_basename "./data/raw/dbs/pfam_fs_cut_sample/B${i}/pfam" --remove_cif_extension True
    ./scripts/convert_tsv2fsdb.sh ./data/raw/dbs/pfam_fs_cut_sample/B${i}/pfam ./tmp/logs/makedb/tsv2fsdb/pfam_fs_cut_sample_B${i}
done
```


The following snippet, converts each fsdb to cal format and then converts the cal to bca format:
```
for file in $(find "./data/raw/dbs/" -iname "*.lookup"); do
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
    echo "reseek -search ./data/raw/dbs/pfam_cif_cut_sample/B${i}/pfam.bca -db ./data/raw/dbs/pfam_cif_cut_clust/pfam.bca -output data/alis/sample_vs_pfam/reseek_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -verysensitive -evalue 1e99" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/fscut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/raw/dbs/pfam_fs_cut_sample/B${i}/pfam ./data/raw/dbs/pfam_fs_cut_clust/pfam data/alis/sample_vs_pfam/fs_cut_B${i}.tsv tmp/pfam_fs_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


sh_path=./tmp/cifcut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/raw/dbs/pfam_cif_cut_sample/B${i}/pfam ./data/raw/dbs/pfam_cif_cut_clust/pfam data/alis/sample_vs_pfam/cif_cut_B${i}.tsv tmp/pfam_cif_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/mm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "mmseqs easy-search --prefilter-mode 2 -e inf ./data/raw/dbs/pfam_cif_cut_sample/B${i}/pfam.fasta ./data/raw/dbs/pfam_cif_cut_clust/pfam.fasta data/alis/sample_vs_pfam/mm_B${i}.tsv tmp/pfam_mm_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/tm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do  
    echo "foldseek easy-search --exhaustive-search 1 -e inf ./data/raw/dbs/pfam_cif_cut_sample/B${i}/pfam ./data/raw/dbs/pfam_cif_cut_clust/pfam data/alis/sample_vs_pfam/tm_B${i}.tsv tmp/pfam_tm_B${i} --alignment-type 1 --tmscore-threshold 0.0 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore" >> ${sh_path}
done  
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "12:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
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
Run time on Niagara Node of Compute Canada: 1 hour. 
```
mkdir -p ./data/first_label_occ

file_paths=$(find ./data/alis/sample_vs_pfam/ -type f -name "*.tsv")
echo "$file_paths" | parallel "python scripts/find_nonred_labels.py --input {}"
```


```
# The next lines will make the clustered version of pfam_fl, which I skip for now
# python scripts/find_pfam_fl_clust_reps.py #This will make a file containing the name of proteins that were used for making PfamSDB
# make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_fl_clust_reps.tsv" SRC_DB="./data/raw/dbs/pfam_fl/pfam_fl" OUT_DB="./data/raw/dbs/pfam_fl_clust/pfam_fl"
```

## Searching a subset of Pfam_fls to search against clustered PfamSDB
```
shuf --random-source=./data/pfam_fl_clust_reps.tsv -n 1000 ./data/pfam_fl_clust_reps.tsv -o ./data/pfam_fl_sample1.tsv
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/pfam_fl_sample1.tsv" SRC_DB="./data/raw/dbs/pfam_fl/pfam_fl" OUT_DB="./data/raw/dbs/pfam_fl_sample1/pfam_fl"
python scripts/convert_tsv2cal.py  --input_basename ./data/raw/dbs/pfam_fs_cut_clust/pfam
reseek -convert ./data/raw/dbs/pfam_fs_cut_clust/pfam.cal -bca ./data/raw/dbs/pfam_fs_cut_clust/pfam.bca
python scripts/convert_fsdb2cal.py --input_basename ./data/raw/dbs/pfam_fl_clust/pfam_fl


python scripts/convert_fsdb2cal.py --input_basename ./data/raw/dbs/pfam_fl_sample1/pfam_fl
reseek -convert ./data/raw/dbs/pfam_fl_sample1/pfam_fl.cal -bca ./data/raw/dbs/pfam_fl_sample1/pfam_fl.bca
foldseek easy-search --exhaustive-search 1 -e inf ./data/raw/dbs/pfam_fl_sample1/pfam_fl ./data/raw/dbs/pfam_fs_cut_clust/pfam ./data/sample1_pfam_foldseek_fscut.tsv tmp/pfam_cifcut
reseek -search ./data/raw/dbs/pfam_fl_sample1/pfam_fl.bca -db ./data/raw/dbs/pfam_fs_cut_clust/pfam.bca \
-output ./data/sample1_pfam_reseek.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq \
-verysensitive -evalue 1e99
```


## Investigating exhaustive pairwise alignments between members of a domain
```
python ./scripts/group_pfamclust_by_family.py #This will make a directory for each family in "data/raw/dbs/pfam_fs_cut_clust/grp_by_family" directory 
# inside each directory, tsv files of the database will be stored
mkdir -p tmp/logs/mkdb/fam_dbs/rs tmp/logs/mkdb/fam_dbs/fs

ls data/raw/dbs/pfam_fs_cut_clust/grp_by_family/PF*/PF*_ca.tsv | xargs -P 20 -I {} sh -c '
    name_wo_ext=$(basename {} _ca.tsv);
    path_wo_ext=${1%_ca.tsv};
    ./scripts/convert_tsv2fsdb.sh $path_wo_ext tmp/logs/mkdb/fam_dbs/fs/${name_wo_ext};
    python ./scripts/convert_tsv2cal.py --input_basename $path_wo_ext;
    reseek -convert ${path_wo_ext}.cal -bca ${path_wo_ext}.bca -threads 1 > tmp/logs/mkdb/fam_dbs/rs/${name_wo_ext} 2>&1;
    rm ${path_wo_ext}.cal;
    python scripts/convert_tsv2fasta.py --input_basename $path_wo_ext;
' -- {}

# Check log files in tmp/logs/mkdb/fam_dbs/ to make sure that everything has proceeded as expected. You can only focus on
# a few ones with the most/least size


# Now, we start the exhaustive search

# First, we create the directories for storing the alignment, temporary files, and log files

for search_type in fs fs3di tm mm rs
do
    mkdir -p data/alis/fam_alis/${search_type}/ tmp/fs_tmp/fam_alis/${search_type} tmp/logs/fam_alis/${search_type}
done
rm -r tmp/fs_tmp/fam_alis/rs #reseek doesn't need a tmp directory

# Now, we search each family against itself.

all_pfs=$(ls data/raw/dbs/pfam_fs_cut_clust/grp_by_family/PF*/PF*_ca)
for pf_path in ${all_pfs}
do
    pf_basename=${pf_path/_ca/}
    pf_name=$(basename $pf_basename)

    search_type=fs
    foldseek easy-search $pf_basename $pf_basename data/alis/fam_alis/${search_type}/${pf_name}.tsv tmp/fs_tmp/fam_alis/${search_type}/${pf_name} --exhaustive-search 1 -e inf --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddtfull,qaln,taln -v 1 > tmp/logs/fam_alis/${search_type}/${pf_name} 2>&1

    search_type=fs3di
    foldseek easy-search $pf_basename $pf_basename data/alis/fam_alis/${search_type}/${pf_name}.tsv tmp/fs_tmp/fam_alis/${search_type}/${pf_name} --alignment-type 0 --exhaustive-search 1 -e inf --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddtfull,qaln,taln -v 1 > tmp/logs/fam_alis/${search_type}/${pf_name} 2>&1

    search_type=tm
    foldseek easy-search $pf_basename $pf_basename data/alis/fam_alis/${search_type}/${pf_name}.tsv tmp/fs_tmp/fam_alis/${search_type}/${pf_name} --alignment-type 1 --tmscore-threshold 0.0 --exhaustive-search 1 -e inf --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddtfull,qaln,taln -v 1 > tmp/logs/fam_alis/${search_type}/${pf_name} 2>&1

    search_type=mm
    mmseqs easy-search ${pf_basename}.fasta ${pf_basename}.fasta data/alis/fam_alis/${search_type}/${pf_name}.tsv tmp/fs_tmp/fam_alis/${search_type}/${pf_name} --prefilter-mode 2 -e inf --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln -v 1 > tmp/logs/fam_alis/${search_type}/${pf_name} 2>&1

    search_type=rs
    reseek -search ${pf_basename}.bca -db ${pf_basename}.bca -output data/alis/fam_alis/${search_type}/${pf_name}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq+qrow+trow -verysensitive -evalue 1e99 > tmp/logs/fam_alis/${search_type}/${pf_name} 2>&1

done
```
