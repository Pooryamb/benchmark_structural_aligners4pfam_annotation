TSV_BASENAME=$1
LOG_PATH=$2

FOLDSEEK="foldseek"


${FOLDSEEK} tsv2db ${TSV_BASENAME}.tsv ${TSV_BASENAME} --output-dbtype 0 -v 2 >> ${LOG_PATH}
${FOLDSEEK} tsv2db ${TSV_BASENAME}_ss.tsv ${TSV_BASENAME}_ss --output-dbtype 0 -v 2 >> ${LOG_PATH}
${FOLDSEEK} tsv2db ${TSV_BASENAME}_h.tsv ${TSV_BASENAME}_h --output-dbtype 12 -v 2 >> ${LOG_PATH}
${FOLDSEEK} tsv2db ${TSV_BASENAME}_ca.tsv ${TSV_BASENAME}_ca --output-dbtype 12 -v 2 >> ${LOG_PATH}
${FOLDSEEK} compressca ${TSV_BASENAME} ${TSV_BASENAME}_ca2 --coord-store-mode 2 -v 2 >> ${LOG_PATH}
${FOLDSEEK} rmdb ${TSV_BASENAME}_ca -v 2 >> ${LOG_PATH}
${FOLDSEEK} mvdb ${TSV_BASENAME}_ca2 ${TSV_BASENAME}_ca -v 2 >> ${LOG_PATH}

