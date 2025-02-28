#!/bin/bash
#SBATCH --job-name=fasterMR
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48GB
#SBATCH -t 6:00:00
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

##### Longest job #18 takes 46 min#############################

# arg1 = outcomes (from array id)

# arg2 - where exposure files are
path_to_exposures=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/03.cis_pQTL/all_strict_v2g_cis/

# arg3 - where to store MR outputs
outdir=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/05.MR_strict_v2g/output/  # remember the end slash!

# arg4 - directory storing the proxies for each exposure
proxy_dir=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/05.MR_strict_v2g/0.find_proxies/found_proxies/


/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla /scratch/richards/chen-yang.su/01.pQTL/scripts/01.faster_MR_proxy_search.R ${SLURM_ARRAY_TASK_ID} ${path_to_exposures} ${outdir} ${proxy_dir}



