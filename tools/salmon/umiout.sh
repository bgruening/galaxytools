#!/usr/bin/bash

mkdir fixed;
for file in ./umiout/*;
do prefix="${file%.dot.gz}";
prefix=${prefix/.\/umiout\//};
gunzip $file;
sed "s/umiout\/$prefix.dot.gz/$prefix/" umiout/$prefix.dot > fixed/$prefix.dot;
dot -Tpdf fixed/$prefix.dot -o fixed/$prefix.pdf;
done
ls fixed
