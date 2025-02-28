#!/bin/bash
#SBATCH --job-name=aric_cispqtl
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24GB
#SBATCH --output=/scratch/richards/chen-yang.su/01.pQTL/01.aric/log/%j.out
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chen-yang.su@mail.mcgill.ca
#SBATCH --array=1-466

# module load r/4.1.2

mkdir -p ld_clumped_cis_trans/  # directory for storing results
    
# read in a single batch of 10 proteins
list=$(head -n ${SLURM_ARRAY_TASK_ID} ./batch_files_aric/listbatch.txt | tail -n 1)

# IFS=$'\n' # Set IFS to only use newline as the delimiter
# # for the case where proteins have genes separated by whitespace: /scratch/richards/satoshi.yoshiji/24.meta_pQTL/01.2.gmetal_format/output/META-pQTL.2744_57.IGHG1 IGHG2 IGHG3 IGHG4 IGK@ IGL@.txt.gz

for prot_path in `cat $list`; do
    echo $prot_path
    /cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla ./01.format_ld_clump_cis_trans.R $prot_path

done


