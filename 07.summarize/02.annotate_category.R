.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(tidyverse)
library(data.table)
library(magrittr)
library(tidyr)

cohort_dir <- "01.aric"

setwd(paste0("/scratch/richards/chen-yang.su/01.pQTL/", cohort_dir, "/07.summarize"))

all_mr_pwcoco_pass <- fread('all_mr_pwcoco_pass.tsv', fill = T)

# outcomes
outcomes <- fread('/scratch/richards/chen-yang.su/data/04.prep_MR_outcome/outcome.all.tsv')
curated_outcomes <- data.frame()
for(i in 1:nrow(outcomes)){
  #i<-1
  ith_outcome <- outcomes[i,]
  #ith_outcome_path <- ith_outcome$gwas_path
  ith_outcome_name <- paste0(str_split(basename(ith_outcome$gwas_path),'[.]')[[1]][1], '.', str_split(basename(ith_outcome$gwas_path),'[.]')[[1]][2]) # e.g., "Age_Abdomen.35418184"
  #ith_outcome_type <- ith_outcome %>% select(Type) %>% pull() %>% unique()
  ith_outcome <- ith_outcome %>% mutate(outcome = ith_outcome_name, .after = 'outcome_name')
  curated_outcomes <- rbind(curated_outcomes, ith_outcome)
}

curated_outcomes %<>% select(outcome_name, outcome, sample_size, case_size, Type)

# combine
curated_all <- left_join(all_mr_pwcoco_pass, curated_outcomes, by = 'outcome') %>% relocate(outcome_name, .before = 'outcome')
curated_all %<>% separate(exposure, into = c('protein', 'seqid'), sep = '[.]', remove = F)

# add relevant columns
curated_all$olinkid <- NA
curated_all$protid <- curated_all$seqid

sl <- fread("/scratch/richards/chen-yang.su/data/01.SL_proteins/SL_proteins.txt")
sl %<>% separate(SeqId, into = c('seqid1', 'seqid2'), sep = '-') %>% unite(seqid, c('seqid1', 'seqid2'), sep = '_')
curated_all <- merge(curated_all, sl[, c("seqid", "UniProt")], all.x=T)
curated_all %<>% rename("uniprot" = "UniProt")

curated_all %<>% relocate(c(uniprot, seqid, olinkid, protid), .after = protein)

write_tsv(curated_all, file = 'all_mr_pwcoco_pass.category.tsv')

# select
curated_select <-  curated_all %>% select(uniprot, protein, protid, outcome_name,  b, lo_ci, up_ci, pval, Type)
write_tsv(curated_select, file = 'all_mr_pwcoco_pass.fornetworkplot.tsv')
