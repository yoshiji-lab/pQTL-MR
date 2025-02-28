#!/bin/bash

### merge proxies across all outcomes and get a unique list of variants we need to find proxies for

# only parameter to change
#######################
cohort=01.aric
# cohort=01.decode
# cohort=01.fenland
# cohort=01.ukbppp
#######################


# arg3 - where proxy list for each outcome is stor
outdir=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/05.MR_strict_v2g/0.find_proxies
num_files=$(ls $outdir/raw_proxies | wc -l)

merged_file=${outdir}/proxy_list.txt 
touch $merged_file

# write all proxy list for each outcome to a merged file
for ((i=1; i<=$num_files; i++))
do
    echo ${outdir}/raw_proxies/${i}_*.txt
    cat ${outdir}/raw_proxies/${i}_*.txt >> $merged_file  
done

echo "Merged all outcomes to file: ${merged_file}!"
num_rsid=$(wc -l < "$merged_file")
echo "Merged file has $num_rsid variants!"


# get only unique rsids
unique_file=${outdir}/unique_proxy_list.txt
sort ${merged_file} | uniq > $unique_file
num_rsid=$(wc -l < "$unique_file")
echo Unique proxy list file has $num_rsid variants!

# rm -rf $outdir/raw_proxies
# rm $outdir/proxy_list.txt



### create batch files
batch_dir=$outdir/batch_files
mkdir -p $batch_dir
split -l 10 $unique_file $batch_dir/"batch_"  # this will create batch_XX files.
readlink -f $batch_dir/batch_* | sort -V > $batch_dir/listbatch.txt



