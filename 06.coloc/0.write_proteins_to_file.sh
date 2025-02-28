#!/bin/bash

# Write paths ----

echo Writing paths...

mkdir -p batch_files  # create output directory if it doesn't exist

ls /home/csu/scratch/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis/* > batch_files/all_strict_v2g_cis.path

count=1  # initialize counter variable
while read -r line; do
    if (( count % 2 == 1 )); then
        # create new output file for every 10 lines starting at count 1
        outfile="batch_$((count/2 + 1))"
        echo "$line" > "batch_files/$outfile"
    else
        # append line to current output file
        echo "$line" >> "batch_files/$outfile"
    fi
    ((count++))
done < batch_files/all_strict_v2g_cis.path

ls /home/csu/scratch/01.pQTL/01.aric/06.coloc/batch_files/batch* | sort -V > /home/csu/scratch/01.pQTL/01.aric/06.coloc/batch_files/listbatch.txt

