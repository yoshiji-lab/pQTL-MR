#!/bin/bash

mkdir -p all_strict_cis_batch  # create output directory if it doesn't exist

ls /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis/cis* > all_strict_cis_batch/all_strict_cis.path

count=1  # initialize counter variable
while read -r line; do
    if (( count % 10 == 1 )); then
        # create new output file for every 10 lines starting at count 1
        outfile="batch_$((count/10 + 1))"
        echo "$line" > "all_strict_cis_batch/$outfile"
    else
        # append line to current output file
        echo "$line" >> "all_strict_cis_batch/$outfile"
    fi
    ((count++))
done < all_strict_cis_batch/all_strict_cis.path

ls /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis_batch/batch* | sort -V > /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis_batch/listbatch.txt

mkdir -p log

