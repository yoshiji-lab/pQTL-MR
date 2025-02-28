#!/bin/bash


# Write deCODE paths ----

# echo Writing decode paths...
# 
# mkdir -p batch_files_decode  # create output directory if it doesn't exist
# 
# ls /scratch/richards/satoshi.yoshiji/11.pQTL/decode_proteomics_2021_modified_pval/*.txt.gz > batch_files_decode/decode.path
# 
# count=1  # initialize counter variable
# while read -r line; do
#     if (( count % 10 == 1 )); then
#         # create new output file for every 10 lines starting at count 1
#         outfile="batch_$((count/10 + 1))"
#         echo "$line" > "batch_files_decode/$outfile"
#     else
#         # append line to current output file
#         echo "$line" >> "batch_files_decode/$outfile"
#     fi
#     ((count++))
# done < batch_files_decode/decode.path
# 
# ls /scratch/richards/chen-yang.su/01.pQTL/01.decode/batch_files_decode/batch* | sort -V > /scratch/richards/chen-yang.su/01.pQTL/01.decode/batch_files_decode/listbatch.txt







# Write Fenland paths ----

# echo Writing Fenland paths...
# 
# mkdir -p batch_files_fenland  # create output directory if it doesn't exist
# 
# ls /scratch/richards/restricted/somalogic_mrc-epid.cam.ac.uk/*.txt.gz > batch_files_fenland/fenland.path
# 
# count=1  # initialize counter variable
# while read -r line; do
#     if (( count % 10 == 1 )); then
#         # create new output file for every 10 lines starting at count 1
#         outfile="batch_$((count/10 + 1))"
#         echo "$line" > "batch_files_fenland/$outfile"
#     else
#         # append line to current output file
#         echo "$line" >> "batch_files_fenland/$outfile"
#     fi
#     ((count++))
# done < batch_files_fenland/fenland.path
# 
# ls /scratch/richards/chen-yang.su/01.pQTL/01.fenland/batch_files_fenland/batch* | sort -V > /scratch/richards/chen-yang.su/01.pQTL/01.fenland/batch_files_fenland/listbatch.txt




# Write ARIC paths ----

echo Writing aric paths...

mkdir -p batch_files_aric  # create output directory if it doesn't exist

ls /scratch/richards/satoshi.yoshiji/10.Proteome/ARIC/EA/*.gz > batch_files_aric/aric.path

count=1  # initialize counter variable
while read -r line; do
    if (( count % 10 == 1 )); then
        # create new output file for every 10 lines starting at count 1
        outfile="batch_$((count/10 + 1))"
        echo "$line" > "batch_files_aric/$outfile"
    else
        # append line to current output file
        echo "$line" >> "batch_files_aric/$outfile"
    fi
    ((count++))
done < batch_files_aric/aric.path

ls /scratch/richards/chen-yang.su/01.pQTL/01.aric/batch_files_aric/batch* | sort -V > /scratch/richards/chen-yang.su/01.pQTL/01.aric/batch_files_aric/listbatch.txt


