#!/bin/bash
#SBATCH --job-name=v2g_aric
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=86GB
#SBATCH -t 8:00:00
#SBATCH --array=1-159
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err

cd $SLURM_SUBMIT_DIR


list=$(head -n ${SLURM_ARRAY_TASK_ID} /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_cis_batch/listbatch.txt | tail -n 1)
for protein_path in `cat $list` # protein_path now corresponds to only one protein's path (contains paths to decode's, fenland's, and aric's pQTL)
do
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla 02.v2g_filter.R $protein_path
done
