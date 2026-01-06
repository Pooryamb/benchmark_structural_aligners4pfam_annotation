#!/bin/sh

#SBATCH --time=00:15:00
#SBATCH --output=tmp/logs/search/timing_experiment.txt


module load StdEnv hmmer/3.4

dbs_path=./data/raw/dbs/
alis_path="tmp/alis/split_pf_timing/"
tmp_path="tmp/fstmp/split_pf_timing/"
timestamp_path="tmp/timestamps/split_split_single_batch/"
CPU_NUM=192


mkdir -p $alis_path $tmp_path $timestamp_path

#touch $timestamp_path/started_foldseek
#foldseek easy-search -e inf ${dbs_path}/pfam_split_query/pfam ${dbs_path}/pfam_split_target/pfam ${alis_path}/foldseek_timing.tsv ${tmp_path}/pfam_foldseek
#touch $timestamp_path/ended_foldseek


touch $timestamp_path/started_mmseqs
mmseqs easy-search -e inf ${dbs_path}/pfam_split_query/pfam.fasta ${dbs_path}/pfam_split_target/pfam ${alis_path}/mmseqs_timing.tsv ${tmp_path}/pfam_mmseqs
touch $timestamp_path/ended_mmseqs


#db_size=128502
#touch $timestamp_path/started_reseek
#reseek -search ${dbs_path}/pfam_split_query/pfam.bca -db ${dbs_path}/pfam_split_target/pfam.bca -output ${alis_path}/reseek_timing.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -sensitive -evalue 1e99 -dbsize ${db_size}
#touch $timestamp_path/ended_reseek


#touch $timestamp_path/started_hmmer
#hmmscan --cpu ${CPU_NUM} --tblout ${alis_path}/hmmscan_timing.tsv --noali ${dbs_path}/pfam_split_target/pfam.hmm ${dbs_path}/pfam_split_query/pfam.fasta > ./tmp/alis/split_pf_timing/hmmer_log.txt
#touch $timestamp_path/ended_hmmer
