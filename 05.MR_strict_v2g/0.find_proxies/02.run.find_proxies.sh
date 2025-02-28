#!/bin/bash
#SBATCH --job-name=find_proxy
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48GB
#SBATCH -t 3:00:00
#SBATCH -o ./log/proxy.%j.out
#SBATCH -e ./log/proxy.%j.err
#SBATCH --array=1-146

cd $SLURM_SUBMIT_DIR

# only parameter to change
#######################
cohort=01.aric
# cohort=01.decode
# cohort=01.fenland
# cohort=01.ukbppp
#######################

#list=$(head -n ${SLURM_ARRAY_TASK_ID} /scratch/richards/chen-yang.su/01.pQTL/01.fenland/05.MR_strict_v2g/batch_files/listbatch.txt | tail -n 1)
list=$(head -n ${SLURM_ARRAY_TASK_ID} ./batch_files/listbatch.txt | tail -n 1)

# arg2 - where exposure files are
path_to_exposures=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/03.cis_pQTL/all_strict_v2g_cis_res.tsv

# arg3 - where to store proxy outputs
outdir=/scratch/richards/chen-yang.su/01.pQTL/${cohort}/05.MR_strict_v2g/0.find_proxies/found_proxies/  # remember the end slash!

# arg4 - path to ref panel to use for proxy search
ref_panel=/scratch/richards/chen-yang.su/data/02.ukb50k_maf0.01/

for rsid in `cat $list`;
do
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla /scratch/richards/chen-yang.su/01.pQTL/scripts/02.find_proxies.R $rsid $path_to_exposures $outdir $ref_panel
echo "done: $rsid"
done

echo Finished job $SLURM_SUBMIT_DIR



