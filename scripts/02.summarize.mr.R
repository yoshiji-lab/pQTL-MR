.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(tidyverse)
library(data.table)
library(magrittr)
library(dplyr)

args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
cohort_dir <- args[1]  # cohort_dir
#cohort_dir <- "01.aric"

path <- paste0("/scratch/richards/chen-yang.su/01.pQTL/", cohort_dir, "/05.MR_strict_v2g/")
setwd(path)
system('mkdir -p output_summary/')
# system('mkdir -p output_summary/or/')
# system('mkdir -p output_summary/hetero/')
# system('mkdir -p output_summary/pleio/')
# system('mkdir -p output_summary/steiger/')

protein.list <- basename(list.dirs('output', recursive = F))

all_combine_mr_res <- data.frame()
for(protein_name in protein.list){
  print(protein_name)
  # e.g., protein_name<- 'PCSK9.5231_79'
  # e.g., protein_name <- 'HBZ.6919_3'
  protein_outcome_combinations <-list.files(paste0('output/', protein_name, '/or/'), recursive = F)  # collect Wald ratio or IVW results
  protein_outcome_combinations <- gsub("\\.or\\.txt", "", protein_outcome_combinations)
  
  combine_mr_res <- data.frame()
  for(ith_protein_outcome in protein_outcome_combinations){
      # e.g., ith_protein_outcome  <- 'PCSK9.5231_79--25OHD.32242144'
      # e.g., ith_protein_outcome  <- 'HBZ.6919_3--T2D.35551307'
      # name of ith_protein_outcome
      protein_outcome_name <- data.frame(protein_outcome = ith_protein_outcome)
      # collect wald ratio or IVW result
      ith_or <- fread(paste0('output/', protein_name, '/or/', ith_protein_outcome, '.or.txt'), fill = T) %>% filter(method == 'Wald ratio' | method == 'Inverse variance weighted')
      ith_or %<>% mutate(pval = ifelse(pval==0, 1e-300, pval))
      # collect weighted median
      ith_median <- fread(paste0('output/', protein_name, '/or/', ith_protein_outcome, '.or.txt'), fill = T) %>% filter(method == 'Weighted median')
      if(nrow(ith_median) > 0){
      ith_median %<>% transmute(b.median = b, se.median =se, pval.median = pval, lo_ci.median = lo_ci, up_ci.median = up_ci, or.median = or, or_lci95.median = or_lci95, or_uci95.median = or_uci95)
      ith_median %<>% mutate(pval.median = ifelse(pval.median ==0, 1e-300, pval.median))} else {
      ith_median <- data.frame(b.median = NA, se.median = NA, pval.median = NA, lo_ci.median = NA, up_ci.median =NA, or.median = NA, or_lci95.median= NA, or_uci95.median = NA)
      }
      
      # collect weighted mode
      ith_mode <- fread(paste0('output/', protein_name, '/or/', ith_protein_outcome, '.or.txt'), fill = T) %>% filter(method == 'Weighted mode')
      if(nrow(ith_mode) > 0){
      ith_mode %<>% transmute(b.mode = b, se.mode =se, pval.mode = pval, lo_ci.mode = lo_ci, up_ci.mode = up_ci, or.mode = or, or_lci95.mode = or_lci95, or_uci95.mode = or_uci95)
      ith_mode %<>% mutate(pval.mode = ifelse(pval.mode ==0, 1e-300, pval.mode))} else {
      ith_mode <- data.frame(b.mode = NA, se.mode = NA, pval.mode = NA, lo_ci.mode = NA, up_ci.mode =NA, or.mode = NA, or_lci95.mode= NA, or_uci95.mode = NA)
      }
      
      # collect Egger (beta)
      ith_egger <- fread(paste0('output/', protein_name, '/or/', ith_protein_outcome, '.or.txt'), fill = T) %>% filter(method == 'MR Egger')
      if(nrow(ith_egger) > 0){
      ith_egger %<>% transmute(b.egger = b, se.egger =se, pval.egger = pval, lo_ci.egger = lo_ci, up_ci.egger = up_ci, or.egger = or, or_lci95.egger = or_lci95, or_uci95.egger = or_uci95)
      ith_egger %<>% mutate(pval.egger = ifelse(pval.egger ==0, 1e-300, pval.egger)) } else {
      ith_egger <- data.frame(b.egger = NA, se.egger =NA, pval.egger = NA, lo_ci.egger = NA, up_ci.egger = NA, or.egger =NA, or_lci95.egger = NA, or_uci95.egger = NA)
      }
      
      # collect hetero
      hetero_path <- paste0('output/', protein_name, '/hetero/', ith_protein_outcome, '.hetero.txt')
      if(file.exists(hetero_path) & file.info(hetero_path)$size > 0){
        ith_hetero <- fread(paste0('output/', protein_name, '/hetero/', ith_protein_outcome, '.hetero.txt'), fill = T) %>% filter(method == 'Inverse variance weighted')
        ith_hetero %<>% dplyr::select(Q, Q_df, Q_pval, i_squared)
        ith_hetero %>% mutate(Q_pval = ifelse(Q_pval ==0, 1e-300, Q_pval))
      } else {
        ith_hetero <- data.frame(Q = NA, Q_df = NA, Q_pval = NA, i_squared = NA)
      }
      # collect pleio
      pleio_path <- paste0('output/', protein_name, '/pleio/', ith_protein_outcome, '.pleio.txt')
      if(file.exists(pleio_path) & file.info(pleio_path)$size > 0) {
      ith_pleio <- fread(paste0('output/', protein_name, '/pleio/', ith_protein_outcome, '.pleio.txt'), fill = T) 
      ith_pleio %<>% transmute(egger_intercept = egger_intercept, se.egger_intercept = se, pval.egger_intercept = pval)
      ith_pleio %<>% mutate(pval.egger_intercept = ifelse(pval.egger_intercept ==0, 1e-300, pval.egger_intercept))
      } else {
        ith_pleio <- data.frame(egger_intercept = NA, se.egger_intercept = NA, pval.egger_intercept = NA)
      }
      # collect steiger
      steiger_path <- paste0('output/', protein_name, '/steiger/', ith_protein_outcome, '.steiger.txt')
      if(file.exists(steiger_path) & file.info(steiger_path)$size > 0) {
        ith_steiger <- fread(paste0('output/', protein_name, '/steiger/', ith_protein_outcome, '.steiger.txt'), fill = T) 
        ith_steiger %<>% dplyr::select(snp_r2.exposure, snp_r2.outcome, correct_causal_direction, steiger_pval)
        ith_steiger %<>% mutate(steiger_pval = ifelse(steiger_pval ==0, 1e-300,steiger_pval))
        } else {
        ith_steiger <- data.frame(snp_r2.exposure = NA, snp_r2.outcome = NA, correct_causal_direction = NA, steiger_pval = NA)
        }
        
      # aggregate these results
        ith_aggregate <- cbind(protein_outcome_name, ith_or,  ith_hetero, ith_median, ith_mode, ith_egger,  ith_pleio, ith_steiger)
      # combine
        combine_mr_res <- rbind(combine_mr_res, ith_aggregate)
        
      } # closes for(ith_protein_outcome in protein_outcome_combinations){
      
  # add combine_mr_res (results for one particular protein) to all_combine_mr_res 
    all_combine_mr_res <- rbind(all_combine_mr_res, combine_mr_res)
} # closes for(protein_name in protein.list){

# save results
write_tsv(all_combine_mr_res, file = 'output_summary/all_combine_mr_res.tsv')

