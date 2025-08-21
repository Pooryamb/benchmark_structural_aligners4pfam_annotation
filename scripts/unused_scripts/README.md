This file includes bried description of the files copied in this directory:



# Classification based on important residues
This part has scripts for finding structurally similar pairs of proteins:

```
mkdir -p tmp/structurally_similar_pairs
python scripts/select_hits_with_equal_num_of_tp_fp.py # selects an equal number of TPs and FPs from reseek and foldseek
# module load StdEnv python/3.11.5
# source ~/pooryamb/str_align_benchmark/bin/activate
python scripts/subselect_ali_seq.py # Selects the alignments for the tables made in the previous step
```
