.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1", 
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14", 
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))
library(data.table)
library(dplyr)
library(tidyverse)

# Compile cis trans results
system(paste0("mkdir -p /scratch/richards/chen-yang.su/01.pQTL/01.aric/cis_trans_pQTL_identification"))

protein_files <- list.files("/scratch/richards/chen-yang.su/01.pQTL/01.aric/ld_clumped_cis_trans", full.names = TRUE)

i <- 1

all_res <- NULL
seq_ids <- c()
gene_names <- c()
gene_seqids <- c()


for (protein_file in protein_files) {
  print(i)
  prot_id <- basename(protein_file)  # 10000_28.CRYBB2.txt.gz
  
  seq_id <- str_split(prot_id, "\\.")[[1]][1]  # 10000_28
  gene_name <- str_split(prot_id, "\\.")[[1]][2]  # CRYBB2
  gene_seqid <- paste0(gene_name, ".", seq_id)
  print(gene_seqid)
  
  # intialize empty dataframe with column headers
  all_variants <- fread(protein_file)
  #all_variants$chromosome <- sapply(str_split(all_variants$SNP, ":"), "[", 2)
  #all_variants$base_pair_location <- sapply(str_split(all_variants$SNP, ":"), "[", 3)
  
  
  if (nrow(all_variants) > 0) {

    all_res <- rbind(all_res, all_variants)
    seq_ids <- c(seq_ids, rep(seq_id, nrow(all_variants)))
    gene_names <- c(gene_names, rep(gene_name, nrow(all_variants)))
    gene_seqids <- c(gene_seqids, rep(gene_seqid, nrow(all_variants)))
  }
  i <- i + 1
}

all_res$seqid <- seq_ids
all_res$gene <- gene_names
all_res$protein <- gene_seqids

path_to_save <- paste0("/scratch/richards/chen-yang.su/01.pQTL/01.aric/cis_trans_pQTL_identification/all_cis_trans_res.tsv")
write_tsv(all_res, path_to_save)


#### Get stats ----
nrow(all_res)  # Number of total variants after meta-analysis: 2602
length(unique(all_res$protein))  # Number of proteins with at least one variant (cis or trans): 1601

# cis
all_res_cis <- all_res[all_res$cis_trans == "cis",]  
nrow(all_res_cis)  # Number of cis variants: 2586
length(unique(all_res_cis$protein)) # Number of proteins with at least one cis variant: 1596
num_cis_per_prot <- all_res_cis %>% group_by(protein) %>% dplyr::count() %>% arrange(desc(n))
write_tsv(all_res_cis, "/scratch/richards/chen-yang.su/01.pQTL/01.aric/cis_trans_pQTL_identification/all_cis_res.tsv")  
write_tsv(num_cis_per_prot, "/scratch/richards/chen-yang.su/01.pQTL/01.aric/cis_trans_pQTL_identification/num_cis_per_prot.tsv") 

all_res_trans <- all_res[all_res$cis_trans == "trans",]  
nrow(all_res_trans)  # 16 trans

nrow(all_res_cis) + nrow(all_res_trans) == nrow(all_res)


