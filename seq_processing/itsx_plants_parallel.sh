#!/bin/bash

#SETTING VARIABLES
#VSEARCH=$(which vsearch232)
SPLITS=40
mkdir itsx_cut

TMP_FASTA1=$(mktemp)


for f in `ls *.fas` ; do

    # Discard erroneous sequences and add expected error rates

       vsearch  --fastx_filter $f \
        --sizein \
        --sizeout \
        --fasta_width 0 \
        --relabel_sha1 \
        --fastaout itsx_cut/$f.relabel
done

for f in `ls itsx_cut/*.fas.relabel` ; do

	ITSx -i $f --complement F -t T --preserve T -o $f.out &
done

wait

for f in `ls itsx_cut/*out.ITS2.fasta` ; do

	# Dereplicate (vsearch)
	vsearch --derep_fulllength $f \
             --sizein \
             --sizeout \
             --fasta_width 0 \
             --minuniquesize 1 \
             --relabel_sha1 \
             --output "${f/fas.relabel.out.ITS2.fasta/fas}" > /dev/null
done

rm itsx_cut/*relabel*
