.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(tidyverse)
library(data.table)
library(magrittr)
library(ggrepel)
library(dplyr)
library(tidyr)
	
###############
# parameter
###############

cohort_dir <- "01.aric"

setwd(paste0("/scratch/richards/chen-yang.su/01.pQTL/", cohort_dir, "/07.summarize"))

mr_res0 <- fread(paste0("/scratch/richards/chen-yang.su/01.pQTL/", cohort_dir, "/05.MR_strict_v2g/output_summary/all_combine_mr_res.tsv"), fill = T)
mr_res <- mr_res0 %>% filter(!is.na(pval))

# change outcome and exposure columns
mr_res %<>% separate(protein_outcome, into = c('exposure', 'outcome'), sep = '--', remove = F) 

# filter MR-significant and sensitivity analyses-passing proteins
# MR
mr_res$qval <- p.adjust(mr_res$pval, method = 'BH')
mr_res %<>% relocate(qval, .after = pval)
mr_res %<>% mutate(Significance = ifelse(qval < 0.05, 'Pass', 'Fail'))

# heterogeneity
mr_res %<>% mutate(Heterogeneity_test = case_when(
  i_squared >= 0.5 & Q_pval < 0.05 ~ 'Fail',
  i_squared < 0.5 | Q_pval >0.05 ~ 'Pass',
  is.na(i_squared) ~ 'NA'))

# consistency
mr_res %<>% mutate(Consistency_test = case_when(
  is.na(b.median) | is.na(b.mode) | is.na(b.egger) ~ 'NA',
  b*b.median > 0 & b*b.mode > 0 & b*b.egger> 0 ~ 'Pass',
  TRUE ~ 'Fail'
))

# pleiotropy
mr_res %<>%
  mutate(Pleiotropy_test = case_when(
    (pval.egger_intercept >= 0.05) ~ 'Pass',
    is.na(pval.egger_intercept) ~ 'NA',
    TRUE ~ 'Fail'
  ))


# combine causal estimate, Heterogeneity_test, and Pleiotropy_test
mr_res %>% filter(Significance == "Pass") %>% filter(Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') %>% nrow() 
mr_res %>% filter(Significance == "Pass") %>% filter(Heterogeneity_test == 'Pass'  | Heterogeneity_test == 'NA') %>% filter(Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA') %>% nrow() 
mr_res %>% filter(Significance == "Pass") %>% filter(Heterogeneity_test == 'Pass'  | Heterogeneity_test == 'NA') %>% filter(Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA') %>% filter(Consistency_test == 'Pass' | Consistency_test == 'NA') %>% nrow()
mr_res %>% filter(Significance == "Pass") %>% filter(Heterogeneity_test == 'Pass'  | Heterogeneity_test == 'NA') %>% filter(Pleiotropy_test == 'Pass' | Pleiotropy_test == 'NA') %>% filter(Consistency_test == 'Pass' | Consistency_test == 'NA') %>% filter(correct_causal_direction == TRUE) %>% nrow()

mr_res %<>% mutate(Sensitivity_test = case_when(
  (Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') & 
    (Pleiotropy_test == 'Pass'| Pleiotropy_test == 'NA') &
    (Consistency_test == 'Pass' | Consistency_test == 'NA') &
    correct_causal_direction == TRUE ~ 'Pass',
  TRUE ~ 'Fail'))

# count
table(mr_res[c('Significance', 'Sensitivity_test')])

mr_res %<>% mutate(All_pass = case_when(Significance == 'Pass' & 
                                         (Heterogeneity_test == 'Pass' | Heterogeneity_test == 'NA') & 
                                         (Pleiotropy_test == 'Pass'| Pleiotropy_test == 'NA') &
                                         (Consistency_test == 'Pass' | Consistency_test == 'NA') &
                                         correct_causal_direction == TRUE ~ 'Pass',
                                         TRUE ~ 'No'))

table(mr_res[c('Significance', 'All_pass')])


# create a new table only with proteins passing Sinificance and all sensitivity tests
mr_pass <- mr_res %>% filter(All_pass == 'Pass')

# master data frame to store all MR and coloc results
all_combine <- data.frame()

# loop through all MR-passing protein-outcome associations
for(i in 1:nrow(mr_pass)){
  # i <- 442
  ith_prot_out_res <- mr_pass[i,]
  # e.g., /scratch/richards/satoshi.yoshiji/24.meta_pQTL/06.coloc/01.coloc_res/PCSK9.5231_79/PCSK9.5231_79--25OHD.32242144/PCSK9.5231_79--25OHD.32242144.coloc
  protein_outcome_name <- ith_prot_out_res$protein_outcome
  protein_name <- ith_prot_out_res$exposure
  outcome_name <- ith_prot_out_res$outcome
  pwcoco_res_path <- paste0("/scratch/richards/chen-yang.su/01.pQTL/", cohort_dir, "/06.coloc/01.coloc_res/", protein_name, '/', protein_outcome_name, '/', protein_outcome_name, '.coloc')
  if(file.exists(pwcoco_res_path)){
      ith_pwcoco_res <- fread(pwcoco_res_path, fill = T)
      # ith_pwcoco_res <- fread('/scratch/richards/satoshi.yoshiji/24.meta_pQTL/06.coloc/01.coloc_res/PCSK9.5231_79/PCSK9.5231_79--25OHD.32242144/PCSK9.5231_79--25OHD.32242144.coloc')
      # which is the pair with highest PPH4?
      if(ith_pwcoco_res[1,]$H4 >= 0.8){ # unconditional colocalization's PPH4 is >= 0.8 (no pwcoco)
        ith_combine <- ith_prot_out_res %>% mutate(pwcoco_highest_PPH4 = ith_pwcoco_res[1,]$H4, pwcoco_SNP1 = ith_pwcoco_res[1,]$SNP1, pwcoco_SNP2 = ith_pwcoco_res[1,]$SNP2)
      } else{
        pair_maxH4 <- ith_pwcoco_res %>% arrange(desc(H4)) %>% filter(row_number() == 1)
        ith_combine <- ith_prot_out_res %>% mutate(pwcoco_highest_PPH4 = pair_maxH4$H4, pwcoco_SNP1 = pair_maxH4$SNP1, pwcoco_SNP2 = pair_maxH4$SNP2)
      } 
  } else{ 
       ith_combine <- ith_prot_out_res %>% mutate(pwcoco_highest_PPH4 = NA, pwcoco_SNP1 = NA, pwcoco_SNP2 = NA)
  } # close if-else
  # add ith_combine to master data frames 
  all_combine  <- rbind(all_combine, ith_combine) 
} # close for(ith_prot_out_res in mr_pass){

all_combine %<>% mutate(Colocalization_test = ifelse(pwcoco_highest_PPH4 >= 0.8, 'Pass', 'No'))
all_mr_pwcoco_pass <- all_combine %>% filter(Colocalization_test == 'Pass')

# count
table(all_combine[['Colocalization_test']])

# save
write_tsv(all_mr_pwcoco_pass, file = 'all_mr_pwcoco_pass.tsv') # save 

#########################
## volcano plot
#########################
#
## add name label
#sum5 <- mr_pwcoco_res %>% filter(Pleiotropy_test == 'Pass')
#
## top 10/3 positive-beta proteins
#top10_positivebeta <- sum5[sum5$b > 0,][1:10,1] %>%
#  unlist()
## top 10/3 negative-beta proteins
#top10_negativebeta <- sum5[sum5$b < 0,][1:5,1] %>%
#  unlist()
## proteins to be annotated
#top10_bothsides <- c(top10_negativebeta, top10_positivebeta)
#
#sum5$namelabel <- NA
#sum5$namelabel[sum5$protein %in% top10_bothsides] <- sum5$protein[sum5$protein %in% top10_bothsides]
##sum5$namelabel[(sum5$Significance == 'Bonferroni') & (abs(sum5$b) > 0.35)] <- sum5$protein[(sum5$Significance == 'Bonferroni') & (abs(sum5$b) > 0.35)]
##sum5$namelabel[sum5$Significance == 'Bonferroni' | sum5$Significance == 'Nominal'] <- sum5$protein[sum5$Significance == 'Bonferroni' | sum5$Significance == 'Nominal']
##sum5
#
#(
#  volcano <-
#    ggplot(data = sum5,
#           aes(x = b, y = -log10(pval),
#               col = Significance,
#               shape = Sensitivity_test,
#               label = namelabel
#           )) +
#    geom_point(size = 2) +
#    theme_classic() +
#    #geom_label_repel(aes(label = namelabel), box.padding = 0.6, fill = "white")+
#    geom_text_repel(max.overlaps = Inf) +  # needs this to show namelabel +
#    geom_hline(yintercept = c(-log10(0.05),-log10(1e-5)), col = c( "#3C5488FF", "#DC0000FF")) +
#    scale_color_manual(values=c( 'grey60', "#3C5488FF", "#DC0000FF")) + # grey blue red
#    theme(text = element_text(size = 20)) +
#    scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
#    scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
#    xlab('Beta') +
#    ylab(expression(-log[10](italic("P"))))
#)
#
#ggsave(filename = 'volcano_beta_nolabel.pdf',
#       volcano,
#       width = 10,
#       height = 8)
