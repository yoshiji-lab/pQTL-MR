#!/bin/bash
#SBATCH --job-name=00proxy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16GB
#SBATCH -t 30:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --array=1-181

cd $SLURM_SUBMIT_DIR

# only parameter to change
#######################
cohort=01.aric
# cohort=01.decode
# cohort=01.fenland
# cohort=01.ukbppp
#######################

# arg1 = outcomes (from array id)

# arg2 - where exposure files are
path_to_exposures=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/03.cis_pQTL/all_strict_v2g_cis_res.tsv

# arg3 - where to store proxy outputs
outdir=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/05.MR_strict_v2g/0.find_proxies/raw_proxies/  # remember the end slash!

/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla /scratch/richards/chen-yang.su/01.pQTL/scripts/00.get_proxy_list_for_single_outcome.R ${SLURM_ARRAY_TASK_ID} ${path_to_exposures} ${outdir}


echo Finished Outcome ${SLURM_ARRAY_TASK_ID}!

