#!/usr/bin/env bash

known_te=($1)

echo 'get classification info and convert #unknown to #DNA/Helitron'
for j in *mod.EDTA.TElib.novel.fa; do 
    for i in `cat $j.real`; do 
        grep $i $j; 
    done| \
    perl -nle 's/#unknown/#DNA\/Helitron/; print $_' > $j.real.ori & 
done

echo 'aggregate novel TE libraries'
i=0
for j in *real.ori; do
  i=$(($i+5000));
  perl /EDTA/util/rename_TE.pl $j $i;
done > NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw
perl /EDTA/util/rename_TE.pl NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw > NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw2
mv NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw2 NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw

ls

echo 'remove redundant'
nohup perl /EDTA/util/cleanup_nested.pl \
    -in NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw \
    -cov 0.95 \
    -minlen 80 \
    -miniden 80 &

ls

echo 'remove a number of false TEs and rename IDs'
RepeatMasker -pa 36 -q -no_is -norna -nolow -div 40 -lib rm.fa -cutoff 225 NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw.cln
perl /EDTA/util/output_by_list.pl 1 NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw.cln 1 \
    <(awk '{print $5}' NAM.EDTA1.8.0.EDTA.TElib.novel.fa.raw.cln.out|grep TE) -ex -FA | \
    perl /EDTA/util/rename_TE.pl - > NAM.EDTA1.8.0.EDTA.TElib.novel.fa

echo 'make comprehensive TE library'
cat $known_te NAM.EDTA1.8.0.EDTA.TElib.novel.fa > NAM.EDTA1.8.0.TE11122019.TElib.fa

echo 'finished make_pan_library'

