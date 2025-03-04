mkdir -p ./tmp/pfam_fl
mkdir -p ./tmp/extra
python ./scripts/convert_pfam_fl_acc2links.py
cat ./tmp/extra/urls.txt | gsutil -m -q cp -I ./tmp/pfam_fl 2> ./tmp/logs/misc/pfam_fl_download_log
python extract_ss_from_cif.py
find ./tmp/pfam_fl -iname "*" | parallel rm 
rm -rf ./tmp/pfam_fl
