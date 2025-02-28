#!/bin/bash
#SBATCH --job-name=summarize
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH -t 3:00:00
#SBATCH -o ./log/sum.%j.out
#SBATCH -e ./log/sum.%j.err

cd $SLURM_SUBMIT_DIR

# only parameter to change
#######################
cohort=01.aric
# cohort=01.decode
# cohort=01.fenland
# cohort=01.ukbppp
#######################


/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla /scratch/richards/chen-yang.su/01.pQTL/scripts/02.summarize.mr.R $cohort

echo Finished summarization for $cohort !

