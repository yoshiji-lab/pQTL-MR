#!/bin/bash
#SBATCH --job-name=pwcoco_aric
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48GB
#SBATCH -t 8:00:00
#SBATCH -o ./log/output.%j.out
#SBATCH -e ./log/output.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=chen-yang.su@mail.mcgill.ca
#SBATCH --array=1-551
#SBATCH --account=rrg-bengioy-ad

cd $SLURM_SUBMIT_DIR

list=$(head -n ${SLURM_ARRAY_TASK_ID} /home/csu/scratch/01.pQTL/01.aric/06.coloc/batch_files/listbatch.txt | tail -n 1)
for protein in `cat $list`;
do
protein_name_tmp=`basename $protein ".tsv"`  #in the format of "cis.A2ML1.15482_12"
protein_name=$(echo "$protein_name_tmp" | sed 's/^cis.//') #in the format of "A2ML1.15482_12"
echo "start: $protein_name"
/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/bin/Rscript --vanilla 01.colocalization_pQTL_trait.R $protein $protein_name
echo "done: $protein_name"
done

echo "Done!"