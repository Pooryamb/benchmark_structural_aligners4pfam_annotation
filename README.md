# Introduction

This repository contains the scripts used for comparing Foldseek and Reseek 
for the structure-based domain annotation task.
You require the following libraries for this benchmarking:
* Pandas (version 2.2.2)
* Polars (version 1.0.0)
* p2rank (version 2.5.1-dev.3)



To start, make sure that both reseek, foldseek, and mmseqs have been added to 
your $PATH variable. Add the rest of the paths accordingly:
```
PRANK_PATH="/home/pooryamb/p2rank/distro/prank"
```

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
foldseek convert2fasta ./data/raw/dbs/pfam_cif_cut/pfam_ss ./data/raw/dbs/pfam_cif_cut/pfam_ss.fasta
mkdir -p tmp/fstmp/ tmp/alidb tmp/alis/sample_pf tmp/alis/pf_pf tmp/alis/pfam_clust tmp/logs/misc/
mmseqs easy-cluster ./data/raw/dbs/pfam_cif_cut/pfam.fasta tmp/alis/pfam_clust/pfam_clust tmp/fstmp/clust --min-seq-id 0.4 -c 0.8 --cov-mode 0 -s 9
mkdir -p ./data/processed/
python scripts/extract_plddt.py --input ./data/raw/pfam_cifs.tar --output ./data/processed/pfam_plddt.json # Extracts plddt from structures
python scripts/calc_avg_plddt_from_json.py # It converts the json of plddts to a tsv file containing the average plddts
python ./scripts/select_clusterrep_by_plddt_and_size.py #It selects the one with highest pLDDT as cluster representative. It also removes domains that are shorter than 10 amino acids
# First, we select the representatives from cif_cut database
make -f ./scripts/Makefile.make_subdb REP_INFO="./data/processed/pfam_clust_reps.tsv" SRC_DB="./data/raw/dbs/pfam_cif_cut/pfam" OUT_DB="./data/raw/dbs/pfam_cif_cut_clust/pfam"
# Then, we select the representatives from fs_cut database
python ./scripts/select_subtsv.py --input_basename "./data/raw/dbs/pfam_fs_cut/pfam" --list2sel ./data/processed/pfam_clust_reps.tsv --output_basename "./data/raw/dbs/pfam_fs_cut_clust/pfam" --remove_cif_extension True
mkdir -p tmp/logs/makedb/tsv2fsdb
./scripts/convert_tsv2fsdb.sh ./data/raw/dbs/pfam_fs_cut_clust/pfam tmp/logs/makedb/tsv2fsdb/make_fs_cut_clust_log
```

## Select a sample of Pfam to query against PfamSDB

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
mkdir -p tmp/timestamps/sample_pf
mkdir -p tmp/alis/sample_pf
mkdir -p tmp/alidbs/
mkdir -p tmp/logs/search/sample_pf
mkdir -p tmp/jobs
mkdir -p tmp/fstmp/sample_pf

dbs_path="./data/raw/dbs/"
alis_path="tmp/alis/sample_pf"
tmp_path="tmp/fstmp/sample_pf"

mkdir -p tmp/alis/sample_pf_sorted/
sh_path=./tmp/jobs/rs_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "reseek -search ${dbs_path}/pfam_cif_cut_sample/B${i}/pfam.bca -db ${dbs_path}/pfam_cif_cut_clust/pfam.bca -output ${alis_path}/reseek_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -verysensitive -evalue 1e99" >> ${sh_path}
done

python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


sh_path=./tmp/jobs/fscut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ${dbs_path}/pfam_fs_cut_sample/B${i}/pfam ${dbs_path}/pfam_fs_cut_clust/pfam ${alis_path}/fs_cut_B${i}.tsv ${tmp_path}/pfam_fs_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


sh_path=./tmp/jobs/cifcut_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ${dbs_path}/pfam_cif_cut_sample/B${i}/pfam ${dbs_path}/pfam_cif_cut_clust/pfam ${alis_path}/cif_cut_B${i}.tsv ${tmp_path}/pfam_cif_cut_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/jobs/mm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "mmseqs easy-search --prefilter-mode 2 -e inf ${dbs_path}/pfam_cif_cut_sample/B${i}/pfam.fasta ${dbs_path}/pfam_cif_cut_clust/pfam.fasta ${alis_path}/mm_B${i}.tsv ${tmp_path}/pfam_mm_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

sh_path=./tmp/jobs/tm_samp_ag_clust_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do  
    echo "foldseek easy-search --exhaustive-search 1 -e inf ${dbs_path}/pfam_cif_cut_sample/B${i}/pfam ${dbs_path}/pfam_cif_cut_clust/pfam ${alis_path}/tm_B${i}.tsv ${tmp_path}/pfam_tm_B${i} --alignment-type 1 --tmscore-threshold 0.0 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "12:00:00"; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

#To run on a single machine:
#bash $sh_path
#To run on a machine operated by JOB scheduler:
#python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "3:00:00"; sbatch ./tmp/${sh_path/.sh/_slurm_job.sh}
#

```
## Sort Reseek output
In the current version of Reseek, the hits are not necessarily sorted by e-value. So, we need to sort the files based on e-value.
```
for i in $(seq 1 $CHUNK_NUM); do
    sort --parallel=80 --buffer-size=50% -t$'\t' -k1,1 -k10,10g tmp/alis/sample_pf/reseek_B${i}.tsv -o tmp/alis/sample_pf_sorted/reseek_B${i}.tsv
done
rm -f tmp/alis/sample_pf/reseek_B*.tsv
mv tmp/alis/sample_pf_sorted/reseek_B*.tsv tmp/alis/sample_pf/
rm -rf tmp/alis/sample_pf_sorted/
```

## Seed characteristics
The following commands will be used to extract sees characteristics such as their secondary structure or contact number.
Calculation of the secondary structure requires downloading all full-lenth structures and might take a long time.
For the calculation of average sequence identity of each seed with other seeds of the same family, copy Pfam-A.seed.gz file (from Pfam website)
inside the ./data/raw directory

```
bash ./scripts/extract_ss_pfam.sh                  # This extracts the secondary structure from AlphaFold full length structures and writes the output in a fasta file.
python ./scripts/calc_ss_perc_and_trans_ratio.py   # This calculates the percantage of each secondary structure state (h/e/c) in the structure and also the length normalized number of transitions
python ./scripts/calc_avg_contact_numbers.py       # This calculates the contact number for the cif files
python ./scripts/extract_pfam_clust_alignments.py  # This will extract the pfam seeds alignments from the Pfam files
python ./scripts/find_avg_intrafam_pident.py       # This finds the average sequence identity of each seed with other members of its family
mkdir -p tmp/cifs                                  # This is for the scripts that require the cif files for feature extraction.
tar -xf ./data/pfam_cifs.tar -C ./tmp/cifs
python scripts/remove_non_clust_cifs.py            # This removes seeds that have not been selected as cluster representatives
find ./tmp/cifs/ -iname "*" -type f | parallel gunzip
python scripts/create_ds_file4prank.py             # Makes a file needed for running p2rank
mkdir -p ./tmp/prank_preds ./tmp/logs/prank/
mv ./tmp/extra/filepaths4prank.txt .; $PRANK_PATH predict -c alphafold filepaths4prank.txt -threads 20 -o ./tmp/prank_preds > ./tmp/logs/prank/prank_log.txt 2> ./tmp/logs/prank/prank_err.txt; mv ./filepaths4prank.txt tmp/extra/
mkdir -p ./tmp/residue_features/
python scripts/store_pocket_coordinates_json.py    # For each seed, writes the positions of highly confident pocket amino acids
python scripts/find_conserved_residues.py          # For each seed, writes the location of conserved residues
### python scripts/make_all_residues_json.py           # Creates a dumb json containing all locations of residues of seeds (we don't need it anymore)
```

To download the information about binding sites and active sites, the information has to be manually downloaded from UniProt website.
Although it is possible to use its api, due to the low number of batches and possibility of having a non-responsive response from UniProt,
it is better to download the data manually. In this regard, first, the UniProt accession of the seeds is written into multiple batches with
less than 100_000 entries in each, because UniProt website cannot handle more than 100_000 entries at a time.
First, use the following script to break the accessions into 3 files:
```
python scripts/write_unipids_in_batches2download_from_uniprot.py
```
Then, upload the ./tmp/extra/unip_ids_b*.txt files to https://www.uniprot.org/id-mapping website. Make sure that
the mapping is done from UniProtKB to UniProtKB. Customize the columns to show "Active site", "Binding site" , and "Sequence".
Then, download the tsv file for each batch, and each as unip_sites_b{i}.tsv.gz and copy them inside the ./data/raw directory.
Find the conserved residues in each sequence using the following commands:
```
gunzip ./data/raw/unip_sites_b*.tsv.gz
python ./scripts/parse_sites_from_unip.py      # Extracts the binding and active sites from UniProt tables and stores them in the json format for the ease of access
python ./scripts/map_sites_on_seeds.py         # Maps the active and binding sites to the seeds. This script simply corrects the offset effect. It does not expand positions from a seed to other family members
python ./scripts/expand_sites_to_all_family_members.py  # Only some seeds had residue annotation, this script will expand residue annotations to all other members of the family

```


## Preprocessing for sensitivity up to the first FP

Now, we find sensitivity by finding the number of TPs before the first FP.
We want to find labels both at family and clan level. So, we need to 
download the file containing the clan data, and make a fake clan id for 
those that don't belong to a clan. For such cases, we use the family id
as the clan id.

```
wget https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam37.0/Pfam-A.clans.tsv.gz -P ./data/raw/ # At the time of working on this project, version 37 was downloaded. Feel free to change the version if needed
gunzip ./data/raw/Pfam-A.clans.tsv.gz
python scripts/process_pf_clans_info.py # This will make a tsv file whose first column is Pfam id and second column is clan id
```

The next script preprocesses the files for finding the sensitivity based on the number of FPs before the first TP.
It finds the row number of the first occurence of each family/clan level label. The labels could be one of the followings:
** family: TP, clan: TP
** family: FP, clan: TP
** family: FP, clan: FP
Next, the jupyter notebook script uses the output of this step for the plots
Run time on Niagara Node of Compute Canada: 1 hour. 
```
mkdir -p data/processed/first_label_occ

file_paths=$(find ./tmp/alis/sample_pf/ -type f -name "*.tsv")
echo "$file_paths" | parallel "python scripts/find_nonred_labels.py --input {}"
```



## Preprocessing for CVE plot
For CVE plot, we plot sensitivity vs error for different e-value thresholds. To select e-value threshold,
we first find the number of rows with each e-value threshold. The following snippet can do this task:
```
mkdir -p data/processed/evalue_bins
file_paths=$(find ./tmp/alis/sample_pf/ -type f -name "*_B*.tsv")
echo "$file_paths" | parallel "python scripts/find_evalue_stratified_labels.py --input {}"

```



## Pfam_fl

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
foldseek lndb ./data/raw/dbs/pfam_cif_cut/pfam_h ./data/raw/dbs/pfam_cif_cut/pfam_ss_h
foldseek convert2fasta ./data/raw/dbs/pfam_cif_cut_clust/pfam_ss ./data/raw/dbs/pfam_cif_cut_clust/pfam_ss.fasta
python ./scripts/group_pfamclust_by_family.py  #This will make a directory for each family in "data/raw/dbs/pfam_cif_cut_clust/grp_by_family" directory
# inside each directory, tsv files of the database will be stored
mkdir -p tmp/logs/makedb/fam_dbs/rs tmp/logs/makedb/fam_dbs/fs

ls data/raw/dbs/pfam_cif_cut_clust/grp_by_family/PF*/PF*_ca.tsv | xargs -P 20 -I {} sh -c '
    name_wo_ext=$(basename {} _ca.tsv);
    path_wo_ext=${1%_ca.tsv};
    ./scripts/convert_tsv2fsdb.sh $path_wo_ext tmp/logs/makedb/fam_dbs/fs/${name_wo_ext};
    python ./scripts/convert_tsv2cal.py --input_basename $path_wo_ext;
    reseek -convert ${path_wo_ext}.cal -bca ${path_wo_ext}.bca -threads 1 > tmp/logs/makedb/fam_dbs/rs/${name_wo_ext} 2>&1;
    rm ${path_wo_ext}.cal;
    python scripts/convert_tsv2fasta.py --input_basename $path_wo_ext;
' -- {}

# Check log files in tmp/logs/makedb/fam_dbs/ to make sure that everything has proceeded as expected. You can only focus on
# a few ones with the most/least size


# Now, we start the exhaustive search

# First, we create the directories for storing the alignment, temporary files, and log files
for search_type in fs fs3di tm mm rs
do
    mkdir -p tmp/alis/fam_alis/${search_type}/ tmp/fstmp/fam_alis/${search_type} tmp/logs/search/fam_alis/${search_type}
done
rm -rf tmp/fstmp/fam_alis/rs #reseek doesn't need a tmp directory

# Now, we search each family against itself.

all_pfs=$(find data/raw/dbs/pfam_cif_cut_clust/grp_by_family/ -type f -name '*_ca')
for pf_path in ${all_pfs}  
do  
    pf_basename=${pf_path%"_ca"}  
    pf_name=$(basename $pf_basename)  

    # Search types  
    for search_type in fs fs3di tm mm rs  
    do  
        tsv_path=tmp/alis/fam_alis/${search_type}/${pf_name}.tsv  
        tmp_path=tmp/fstmp/fam_alis/${search_type}/${pf_name}  
        log_path=tmp/logs/search/fam_alis/${search_type}/${pf_name}_${search_type}.log  

        # Perform searches based on type  
        case $search_type in  
            fs)  
                foldseek easy-search $pf_basename $pf_basename $tsv_path $tmp_path \
                    --exhaustive-search 1 -e inf \
                    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddtfull,qaln,taln \
                    -v 1 > $log_path 2>&1  
                ;;  
            fs3di)  
                foldseek easy-search $pf_basename $pf_basename $tsv_path $tmp_path \
                    --alignment-type 0 --exhaustive-search 1 -e inf \
                    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,lddtfull,qaln,taln \
                    -v 1 > $log_path 2>&1  
                ;;  
            tm)  
                foldseek easy-search $pf_basename $pf_basename $tsv_path $tmp_path \
                    --alignment-type 1 --tmscore-threshold 0.0 --exhaustive-search 1 -e inf \
                    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore,lddtfull,qaln,taln \
                    -v 1 > $log_path 2>&1  
                ;;  
            mm)  
                mmseqs easy-search ${pf_basename}.fasta ${pf_basename}.fasta $tsv_path $tmp_path \
                    --prefilter-mode 2 -e inf \
                    --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln \
                    -v 1 > $log_path 2>&1  
                ;;  
            rs)  
                reseek -search ${pf_basename}.bca -db ${pf_basename}.bca -output $tsv_path \
                    -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq+qrow+trow \
                    -verysensitive -evalue 1e99 > $log_path 2>&1  
                ;;  
        esac  
    done  
done
```

## Investigating residue level alignment
Here, we aim to see how well residue mapping in the pairwise alignment of members of a family aligns with the curated alignments of Pfam.
In this regard, first, we run a code to identify the number of similar mappings between our pairwise aligners (MMseqs, Foldseek, or Reseek)
with Pfam MSAs. Besides investigating the overall alignment, we also check the mapping for specific residues, such as those located inside pockets or the
conserved residues.
```
for res_type in all pocket conserved active_site binding_site; do
    for tool in fs fs3di tm mm rs; do
        mkdir -p ./tmp/intrafam_residue_alignment_counts/${tool}/${res_type}
    done
done

python ./scripts/calc_aligned_residues_count.py      # Calculates the residue alignment counts for each tool, took ~9 hours on a 20 core machine
mkdir -p ./data/processed/residue_ali_frac_per_seed/
python ./scripts/calc_residue_alignment_per_seed.py  # Calculates the average fraction of aligned residues when each seed is used as the query

```


## Benchmarking by searching Pfam Split against Pfam Split

In this part, the Pfam database is splitted into two halves and a random sample of the members of the first half will be searched against the
second half using different tools. This new approach also allows to conduct searches using HMMER.

```
mkdir -p ./data/raw/dbs/pfam_split_query ./data/raw/dbs/pfam_split_target tmp/logs/makedb/splitted_pfam/
for i in {1..16}; do     mkdir -p ./data/raw/dbs/pfam_split_query/B${i}; done  # Query is split up into 16 batches to speed up the searches
python ./scripts/split_pfam_clust_into_tsvs.py   # Creates the tsv files for making the database

# Convert tsv files into files needed by other tools
find ./data/raw/dbs/pfam_split_*/ -iname "*_h.tsv" -type f | while read -r file; do
    wo_postfix=${file/_h.tsv/}
    base_name=$(basename $(dirname $wo_postfix))
    python ./scripts/convert_tsv2fasta.py --input_basename ${wo_postfix}
    bash ./scripts/convert_tsv2fsdb.sh ${wo_postfix} tmp/logs/makedb/splitted_pfam/${base_name}.txt
    python scripts/convert_tsv2cal.py  --input_basename ${wo_postfix}
    reseek -convert ${wo_postfix}.cal -bca ${wo_postfix}.bca
done

python scripts/make_target_split_pf_sto.py  #This will make the stockholm file of the target database
hmmbuild -o ./tmp/logs/makedb/splitted_pfam/pfam_hmm.txt ./data/raw/dbs/pfam_split_target/pfam.hmm ./data/raw/dbs/pfam_split_target/pfam.sto
hmmpress ./data/raw/dbs/pfam_split_target/pfam.hmm

```
The following code script can be used for searching a split of the database against its other split:

```
CHUNK_NUM=16
mkdir -p tmp/timestamps/split_pf
mkdir -p tmp/alis/split_pf
mkdir -p tmp/alidbs/
mkdir -p tmp/logs/search/split_pf
mkdir -p tmp/jobs
mkdir -p tmp/fstmp/split_pf

dbs_path="./data/raw/dbs/"
alis_path="tmp/alis/split_pf"
tmp_path="tmp/fstmp/split_pf"

mkdir -p tmp/alis/split_pf_sorted/

################################################################
###################Reseek-exhaustive############################

search_params=_exh
sh_path=./tmp/jobs/rs_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "reseek -search ${dbs_path}/pfam_split_query/B${i}/pfam.bca -db ${dbs_path}/pfam_split_target/pfam.bca -output ${alis_path}/reseek${search_params}_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -verysensitive -evalue 1e99" >> ${sh_path}
done

python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "1:00:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run

################################################################
###################Reseek-e_3-fast##############################

search_params=_e3_fast
sh_path=./tmp/jobs/rs_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "reseek -search ${dbs_path}/pfam_split_query/B${i}/pfam.bca -db ${dbs_path}/pfam_split_target/pfam.bca -output ${alis_path}/reseek${search_params}_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -fast -evalue 1e-3" >> ${sh_path}
done

python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}

################################################################
###################Reseek-e_3-sens##############################

search_params=_e3_sens
sh_path=./tmp/jobs/rs_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "reseek -search ${dbs_path}/pfam_split_query/B${i}/pfam.bca -db ${dbs_path}/pfam_split_target/pfam.bca -output ${alis_path}/reseek${search_params}_B${i}.tsv -columns query+target+qlo+qhi+ql+tlo+thi+tl+pctid+evalue+aq -sensitive -evalue 1e-3" >> ${sh_path}
done

python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}


################################################################
###################Foldseek-exhaustive##########################

search_params=_exh
sh_path=./tmp/jobs/fs_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ${dbs_path}/pfam_split_query//B${i}/pfam ${dbs_path}/pfam_split_target/pfam ${alis_path}/fs${search_params}_B${i}.tsv ${tmp_path}/pfam_fs${search_params}_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


################################################################
###################Foldseek-e3##################################

search_params=_e3
sh_path=./tmp/jobs/fs_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search -e 1e-3 ${dbs_path}/pfam_split_query//B${i}/pfam ${dbs_path}/pfam_split_target/pfam ${alis_path}/fs${search_params}_B${i}.tsv ${tmp_path}/pfam_fs${search_params}_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


#################################################################
########################MMseqs-exh###############################

search_params=_exh
sh_path=./tmp/jobs/mm_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "mmseqs easy-search --prefilter-mode 2 -e inf ${dbs_path}/pfam_split_query//B${i}/pfam.fasta ${dbs_path}/pfam_split_target/pfam.fasta ${alis_path}/mm${search_params}_B${i}.tsv ${tmp_path}/pfam_mm${search_params}_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


#################################################################
########################MMseqs-e3################################

search_params=_e3
sh_path=./tmp/jobs/mm_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "mmseqs easy-search -e 1e-3 ${dbs_path}/pfam_split_query//B${i}/pfam.fasta ${dbs_path}/pfam_split_target/pfam.fasta ${alis_path}/mm${search_params}_B${i}.tsv ${tmp_path}/pfam_mm${search_params}_B${i}" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "00:15:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


################################################################
######################TM-align##################################
search_params=_exh
sh_path=./tmp/jobs/tm_split_ag_split_exh_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "foldseek easy-search --exhaustive-search 1 -e inf ${dbs_path}/pfam_split_query//B${i}/pfam ${dbs_path}/pfam_split_target/pfam ${alis_path}/tm${search_params}_B${i}.tsv ${tmp_path}/pfam_tm_B${i} --alignment-type 1 --tmscore-threshold 0.0 --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,alntmscore,qtmscore,ttmscore" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "12:00:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}
#bash $sh_path   #This is for running on a single machine. The upper line should be commented if this one is going to be run


################################################################
######################hmmscan_exh############################### #This might be deleted

CPU_NUM=80
search_params=_exh
sh_path=./tmp/jobs/hmm_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "hmmscan --cpu ${CPU_NUM} --tblout ${alis_path}/hmmscan${search_params}_B${i}.tsv -o ${alis_path}/hmmscan${search_params}_large_file_B${i}.tsv --max ./data/raw/dbs/pfam_split_target/pfam.hmm ./data/raw/dbs/pfam_split_query/B${i}/pfam.fasta" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "3:00:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}


################################################################
######################hmmscan_e3################################
CPU_NUM=20
search_params=_e3
sh_path=./tmp/jobs/hmm_split_ag_split${search_params}_commands.sh
rm -f ${sh_path}
for i in $(seq 1 $CHUNK_NUM); do
    echo "hmmscan --cpu ${CPU_NUM} --tblout ${alis_path}/hmmscan${search_params}_B${i}.tsv -o ${alis_path}/hmmscan${search_params}_large_file_B${i}.tsv ./data/raw/dbs/pfam_split_target/pfam.hmm ./data/raw/dbs/pfam_split_query/B${i}/pfam.fasta" >> ${sh_path}
done
python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "3:00:00" --search_category split_pf; sbatch ${sh_path/.sh/_slurm_job.sh}


#To run on a single machine:
#bash $sh_path
#To run on a machine operated by JOB scheduler:
#python ./scripts/make_array_job_file.py --input_sh_path $sh_path --time "3:00:00" --search_category split_pf; sbatch ./tmp/${sh_path/.sh/_slurm_job.sh}
#
```
### Sort Reseek output

It takes ~ 1.5 hours on a system with 80 threads.

```
mkdir -p tmp/alis/split_pf_sorted/
files=$(ls tmp/alis/split_pf/reseek_*.tsv)
CPU_NUM=80

for file in $files; do
    file_base_name=$(basename $file)
    sort --parallel=${CPU_NUM} --buffer-size=50% -t$'\t' -k1,1 -k10,10g $file -o tmp/alis/split_pf_sorted/${file_base_name}
done

rm -f tmp/alis/split_pf/reseek_*.tsv
mv tmp/alis/split_pf_sorted/reseek*.tsv tmp/alis/split_pf/
rm -rf tmp/alis/split_pf_sorted/
```



Next, the first hit for each tool is identified. 
```
mkdir -p ./tmp/first_hits
python ./scripts/select_top_hits.py
```

### Residue level alignment for HMMER

```
python scripts/split_pfam_query_split.py # Makes a fasta file for members of each family in the query split
python scripts/split_pfam_target_split.py # Makes a sto file for members of each family in the target split
mkdir -p ./tmp/logs/makedb/splitted_pfam/hmmer_fam/
# Hmm files are built from sto files
parallel 'file={}; fam=$(basename "$file" .sto); hmmbuild -o ./tmp/logs/makedb/splitted_pfam/hmmer_fam/${fam}.txt ./data/raw/dbs/pfam_split_target/grp_by_family/${fam}.hmm "$file"; hmmpress ./data/raw/dbs/pfam_split_target/grp_by_family/${fam}.hmm' ::: ./data/raw/dbs/pfam_split_target/grp_by_family/*.sto
mkdir -p ./tmp/alis/fam_alis/hmmscan/
parallel 'fam=$(basename {} .fasta); hmmscan --cpu 1 --max --tblout ./tmp/alis/fam_alis/hmmscan/${fam}.tsv -o ./tmp/alis/fam_alis/hmmscan/${fam}_large_file.txt ./data/raw/dbs/pfam_split_target/grp_by_family/${fam}.hmm data/raw/dbs/pfam_split_query/grp_by_family/${fam}.fasta' ::: ./data/raw/dbs/pfam_split_query/grp_by_family/*.fasta


python scripts/find_msa_hmm_mapping.py   # Finds the correspondence between the msa and hmm columns
python scripts/find_q_hmm_mapping.py     # Finds the residue level alignment between random query sample and the hmm profiles
python scripts/find_seed_msa_mapping.py  # Finds the mapping between seed coordinates and curated MSA coordinates

mkdir -p data/processed/residue_ali_frac_per_seed_split_vs_split
python ./scripts/summarize_residue_alignment4hmmscan.py
python ./scripts/adjust_pwa_res_ali4split_vs_split.py

# Performance plots can be generated using the scripts/plot_residue_level_alignments-split_vs_split.ipynb notebook
For what percentage, the first hit was correct
Using different thresholds, what percentage becomes correct, what percentage is incorrect

```
