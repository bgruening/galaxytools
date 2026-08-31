fastafile="$1";
tooldir="$2";

# BLAST
blastn -query $fastafile -db $tooldir/data/stx -task blastn -evalue 0.001 -out stx_blastn -outfmt '6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles' -num_threads 8 -strand both -dust yes -max_target_seqs 1 -perc_identity 95.0;

# SHIGATOXINTYPER: FILTER, CUT AND CONCATENATE BLAST OUTPUT
echo 'sseqid\tpident\tlength\tpositive' > blastn_shigatoxin_fct;
awk -F '\t' '($3>95 && $4>1200) { print $2 FS $3 FS $4 FS $16 }' stx_blastn | sort -nrk 4 -nrk 2 > blastn_shigatoxin_fc;
cat blastn_shigatoxin_fc >> blastn_shigatoxin_fct;
stxfilesize=$(wc -c "blastn_shigatoxin_fc" | awk '{print $1}');
if [ $stxfilesize -eq 0 ]
then
  echo '-\t-\t-\t-' >> blastn_shigatoxin_fct;
fi
