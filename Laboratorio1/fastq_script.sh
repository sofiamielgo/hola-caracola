#!/bin/bash
for laboratorio1 in *.fastq
do echo $laboratorio1
wc -l $laboratorio1
echo "termiando" 
done

