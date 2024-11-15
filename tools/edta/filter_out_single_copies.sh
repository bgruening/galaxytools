#!/usr/bin/env bash

data=($1)

echo 'Filter out single-copy annotations'

for i in `echo $data|awk '{print $1}'`; do
    perl /EDTA/util/output_by_list.pl \
    1 \
    <(perl -nle 's/#.*//; print $_' $i.mod.EDTA.TElib.novel.fa) \
    1 \
    <(perl /EDTA/util/find_flTE.pl $i.mod.out | \
    awk '{print $10}'| \
    sort| \
    uniq -c |\
    perl -nle 'my ($count, $id) = (split); if ($id=~/LTR/){next if $count<=2} else {next if $count ==1} print $_' |\
    awk '{print $2}') -FA > $i.mod.EDTA.TElib.novel.fa.real &
done

