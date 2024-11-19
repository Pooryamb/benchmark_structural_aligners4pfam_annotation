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
python ./scripts/select_clusterrep_by_plddt.py
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
python ./scripts/select_subtsv.py --input_basename "./data/dbs/pfam_fs_cut_clust/pfam" --list2sel ./data/pfam_clust_reps.tsv --output_basename "./data/dbs/pfam_fs_cut_sample/pfam" --remove_cif_extension True
./scripts/convert_tsv2fsdb.sh ./data/dbs/pfam_fs_cut_sample/pfam ./tmp/make_fscut_sample_db_log
```

The following code snippet breaks the sampled data into batches. This is useful for those who want to 
run the search in different batches. Please note that this step is optional. If you have lots of resources, 
you can simply go to the search step. There is a python script that can run the search for both chunked and intact data.

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

```
mkdir -p data/run_times/sample_vs_pfam
mkdir -p data/alis/sample_vs_pfam


rs_search_command="touch data/run_times/sample_vs_pfam/reseek_started.txt?;
reseek -search \
./data/dbs/pfam_cif_cut_sample/pfam.bca# \
-db ./data/dbs/pfam_cif_cut_clust/pfam.bca \
-output data/alis/sample_vs_pfam/reseek.tsv? \
-columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq \
-verysensitive -evalue 1e99;
touch data/run_times/sample_vs_pfam/reseek_ended.txt?"

fs_cut_search_command="touch data/run_times/sample_vs_pfam/fs_cut_started.txt?;
foldseek easy-search --exhaustive-search 1 -e inf \
./data/dbs/pfam_fs_cut_sample/pfam# ./data/dbs/pfam_fs_cut_clust/pfam \
data/alis/sample_vs_pfam/fs_cut.tsv? tmp/pfam_fs_cut?;
touch data/run_times/sample_vs_pfam/fs_cut_ended.txt?"

cif_cut_search_command="touch data/run_times/sample_vs_pfam/cif_cut_started.txt?;
foldseek easy-search --exhaustive-search 1 -e inf \
./data/dbs/pfam_cif_cut_sample/pfam# ./data/dbs/pfam_cif_cut_clust/pfam \
data/alis/sample_vs_pfam/cif_cut.tsv? tmp/pfam_cif_cut?; 
touch data/run_times/sample_vs_pfam/cif_cut_ended.txt?"

tm_search_command="touch data/run_times/sample_vs_pfam/tm_started.txt?;
foldseek easy-search --exhaustive-search 1 -e inf \
./data/dbs/pfam_cif_cut_sample/pfam# ./data/dbs/pfam_cif_cut_clust/pfam \
data/alis/sample_vs_pfam/cif_cut.tsv? tmp/pfam_cif_cut? \
--alignment-type 1 --tmscore-threshold 0.0 \
--format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore;
touch data/run_times/sample_vs_pfam/tm_ended.txt?"

mm_search_command="touch data/run_times/sample_vs_pfam/mm_started.txt?;
mmseqs easy-search --prefilter-mode 2 -e inf \
./data/dbs/pfam_cif_cut_sample/pfam.fasta# ./data/dbs/pfam_cif_cut_clust/pfam.fasta \
data/alis/sample_vs_pfam/mm.tsv? tmp/pfam_mm?;
touch data/run_times/sample_vs_pfam/mm_ended.txt?"


python scripts/run_or_submit_command_job.py --command "$mm_search_command" --job_scheduler \
--time "00:30:00" --add_default_header --batches $CHUNK_NUM --job_name mm_search

python scripts/run_or_submit_command_job.py --command "$fs_cut_search_command" --job_scheduler \
--time "1:00:00" --add_default_header --batches $CHUNK_NUM --job_name fs_cut_search

python scripts/run_or_submit_command_job.py --command "$cif_cut_search_command" --job_scheduler \
--time "1:00:00" --add_default_header --batches $CHUNK_NUM --job_name cif_cut_search

python scripts/run_or_submit_command_job.py --command "$rs_search_command" --job_scheduler \
--time "6:00:00" --add_default_header --batches $CHUNK_NUM --job_name rs_search

python scripts/run_or_submit_command_job.py --command "$tm_search_command" --job_scheduler \
--time "12:00:00" --add_default_header --batches $CHUNK_NUM --job_name tm_search
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
