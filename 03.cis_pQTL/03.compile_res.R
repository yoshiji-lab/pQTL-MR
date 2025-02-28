.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1", 
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14", 
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(data.table)
library(dplyr)
library(tidyverse)


# TODO:
######### Check strict cis, V2G definition results




labeled_files <- list.files("/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis/", full.names = TRUE)

all_res <- NULL
i <- 1
for (file in labeled_files) {
  df <- fread(file, sep="\t", fill=TRUE)
  # df <- df[, c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency",
  #              "p_value", "minus_log10p", "z_score", #"direction", 
  #              "variant_id", "rs_id", "ci_lower", "ci_upper", "n", "cis_trans")]
  if (nrow(df) > 0) {
    # protein_file <- basename(file)
    # prot_name <- gsub("cis\\.|\\.tsv", "", protein_file)
    # 
    # df$protein <- prot_name
    all_res <- rbind(all_res, df)
    
  }
  print(i)
  i <- i + 1
}


nrow(all_res)  # Number of total strict v2g variants after meta-analysis: 1485
length(unique(all_res$protein))  # Number of proteins with at least one variant (cis or trans): 1102

nrow(all_res[all_res$cis_trans == "cis", ])  # Number of cis variants: 1485
length(unique(all_res[all_res$cis_trans == "cis", ]$protein)) # Number of proteins with at least one cis variant: 1102

write_tsv(all_res, "/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis_res.tsv")

num_cis_per_prot <- all_res[all_res$cis_trans == "cis", ] %>% group_by(protein) %>% dplyr::count() %>% arrange(desc(n))

write_tsv(num_cis_per_prot, "/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/num_strict_v2g_cis_per_prot.tsv") 

# Check no pqtls in MHC region
# df <- fread("/scratch/richards/chen-yang.su/01.pQTL/01.ukbppp/03.cis_pQTL/all_strict_v2g_cis_res.tsv")
# mhc_chr <- 6
# mhc_start <- 28510020
# mhc_end <- 33480577 
# mhc_prot <- df[df$chromosome == mhc_chr & df$base_pair_location >= mhc_start & df$base_pair_location <= mhc_end, ]



