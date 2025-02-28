####
######### Find all proxies needed for a single outcome based on file with all the exposure IVs
####

start_time <- Sys.time()

.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(tidyverse)
#library(vroom)
library(TwoSampleMR)
library(stringr)
library(gwasvcf)
library(data.table)

# 1. Read in single outcome GWAS
args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
i <- as.integer(args[1])  # outcome

outcome_path <- fread('/scratch/richards/chen-yang.su/data/04.prep_MR_outcome/outcome.all.tsv', fill = T, sep = '\t')
ith_outcome_path <- outcome_path[i,] 
ith_outcome_name <- paste0(str_split(basename(ith_outcome_path$gwas_path),'[.]')[[1]][1], '.', str_split(basename(ith_outcome_path$gwas_path),'[.]')[[1]][2]) # e.g., "Age_Abdomen.35418184"
ith_outcome_df <- fread(ith_outcome_path$gwas_path, fill = T)
cat("Job Array Number:", i, "- Outcome:", ith_outcome_name, "\n")


# 2. Read in file with all variants to use as exposure.
all_exposure_file <- args[2]  
#all_exposure_file <- "/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis_res.tsv"
Ins_df <- fread(all_exposure_file)
cat("Path to all exposures:", all_exposure_file, "\n")

# 3. Give directory to store output of variants that require proxy search
outdir <- args[3]
#outdir <- '/scratch/richards/chen-yang.su/01.pQTL/01.aric/05.MR_strict_v2g/0.find_proxies/raw_proxies/'
system(paste0("mkdir -p ", outdir))
cat("Output directory:", outdir, '\n')

cat("Check which variants require proxy search... \n\n")






Ins<-format_data(Ins_df, type = "exposure", header = TRUE,
                 #phenotype_col = "", 
                 snp_col = "rs_id", beta_col = "beta",
                 se_col = "standard_error", eaf_col = "effect_allele_frequency", effect_allele_col = "effect_allele",
                 other_allele_col = "other_allele", pval_col = "p_value", min_pval = 1e-300, samplesize_col = 'n')

# retain variants in exposure only
ith_outcome <- ith_outcome_df %>% filter(rs_id %in% Ins$SNP)


################### Perform proxy search if input GWAS cis pQTLs are not in outcome GWAS ##################################
# MR proxy search
# Step 1. Check if variant is in outcome GWAS. If yes, then good and don't do anything, if not then do 2.
# Step 2. Perform proxy search with gwasvcf get_ld_proxies() function and filter out rows where the variant OR proxy of variant has MAF > 0.42
# Step 3. If proxy is in the outcome GWAS, relabel effect allele frequency (EAF) such that the minor allele is the minor allele



# Step 1. ----
proxy_list <- Ins$SNP[!Ins$SNP %in% ith_outcome$rs_id]  # list of cis-pQTLs that we need to find proxies for
cat("List of variants that we need to find proxies for: ", proxy_list, "\n")
cat("Total: ", length(unique(proxy_list)), "\n")

outfile <- paste0(outdir, "/", i, "_", ith_outcome_name, "_proxy_list.txt")
writeLines(proxy_list, outfile)

cat("Saved all variants to file :", outfile, " ---- \n\n")
