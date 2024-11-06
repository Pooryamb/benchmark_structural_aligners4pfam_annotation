input_basename=$1

FOLDSEEK="foldseek"

${FOLDSEEK} lndb ${input_basename}_h ${input_basename}_ca_h
${FOLDSEEK} convert2fasta ${input_basename} ${input_basename}.fasta
${FOLDSEEK} lndb ${input_basename}_h ${input_basename}_ss_h
${FOLDSEEK} convert2fasta ${input_basename}_ss ${input_basename}_ss.fasta
${FOLDSEEK} compressca ${input_basename} ${input_basename}_ca_f64 --coord-store-mode 3
${FOLDSEEK} lndb ${input_basename}_h ${input_basename}_ca_f64_h
${FOLDSEEK} convert2fasta ${input_basename}_ca_f64 ${input_basename}_ca.fasta
${FOLDSEEK} rmdb ${input_basename}_ca_f64
${FOLDSEEK} rmdb ${input_basename}_ca_f64_h
