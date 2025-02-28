library(tidyverse)
library(vroom)
library(magrittr)
library(coloc)
library(data.table)

cohort_dir <- "01.aric"

wd <- paste0('/home/csu/scratch/01.pQTL/', cohort_dir, '/06.coloc/01.coloc_res/')
system(paste0('mkdir -p ', wd))
setwd(wd)

##############
# The code has two layers to run colocalization for protein-outcome associations:
# Outer layer: for-loop for all outcomes (e.g., arg1 = '/scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.2190_55.F11.tsv', arg2 = 'F11.2190_55')
# Inner layer: do colocalization between one protein (designated by a shell script) and the outcome determined by the outer layer
##############

outcome_path <- fread('/home/csu/scratch/data/04.prep_MR_outcome/outcome.all.NARVAL.tsv')

# for-loop for all outcomes 
for(i in 1:nrow(outcome_path)){
  tryCatch({ # try catch outcome error
    #for(i in 1:5){  # for test: i <- 14 #i <-39 (CAD) # i <- 133
	ith_outcome_data <- outcome_path[i,] 
    # outcome name
    ith_outcome_name <- paste0(str_split(basename(ith_outcome_data$gwas_path),'[.]')[[1]][1], '.', str_split(basename(ith_outcome_data$gwas_path),'[.]')[[1]][2]) # e.g., "Age_Abdomen.35418184"
    print(paste0('working on ', ith_outcome_name))
    # outcome path and data
    ith_outcome_path <- ith_outcome_data %>% select(gwas_path) %>% as.character() # read outcome later (as this step is time-consuming)

    # sample size and case size
    ith_outcome_samplesize <- ith_outcome_data %>% select(sample_size) %>% as.numeric() # e.g., ith_outcome_samplesize <- 261984
    ith_outcome_casesize <- ith_outcome_data %>% select(case_size) %>% as.numeric() # e.g., ith_outcome_casesize <- 34541
    
    # pQTL data
    args <- commandArgs(trailingOnly = TRUE) # e.g., 1st column = /scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.PCSK9.5231_79.tsv; 2nd column = PCSK9.5231_79
    # cis-pQTL path
    cis_pqtl_path <- args[1] # cis_pqtl_path <- '/scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.PCSK9.5231_79.tsv'
    #cis_pqtl_path <- '/scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.COL6A3.11196_31.tsv'
    #cis_pqtl_path <- '/scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.F11.2190_55.tsv'
    # proteins (putatively causal proteins)
    protname <- args[2] # e.g. protname <- 'PCSK9.5231_79' protname <- 'F11.2190_55'
    #protname <- 'COL6A3.11196_31'
    seqid <- strsplit(protname, split = '[.]')[[1]][2]

    ################
    # To reduce computational time, proceed to pwcoco only if protein-disease association is nominally significant (pval < 0.05) in MR
    ################
    # mr_res_path <- paste0('/home/csu/scratch/01.pQTL/', cohort_dir, '/05.MR_strict_v2g/output/', protname, '/or/', protname, '--', ith_outcome_name, '.or.txt')
    # mr_res <- vroom(mr_res_path)
    
    all_res_path <- paste0('/home/csu/scratch/01.pQTL/', cohort_dir, '/05.MR_strict_v2g/output_summary/all_combine_mr_res.tsv')
    all_mr_res <- fread(all_res_path)
    prot_outcome <- paste0(protname, '--', ith_outcome_name)
    mr_res <- all_mr_res[all_mr_res$protein_outcome == prot_outcome,] 
    
    mr_res_waldivw_pval <- mr_res %>% filter(method == 'Inverse variance weighted' | method == 'Wald ratio') %>% select(pval) %>% as.numeric()
            if(mr_res_waldivw_pval < 0.05){
            print(paste0('proceed to pwcoco for ', protname, '--', ith_outcome_name))
            # create each protein's dir
            output_dir <- paste0('/home/csu/scratch/01.pQTL/', cohort_dir, '/06.coloc/01.coloc_res/', protname, '/')
            filename <- paste0(protname, '--', ith_outcome_name) # to make naming files easy in the downstream analyses
            system(paste0('mkdir -p ', output_dir, filename))
            setwd(paste0(output_dir, filename))
            
            # lead cis-pQTL
                  if(file.exists(cis_pqtl_path)){ # do the rest only if cis-pQTL exists
                  # start curation
                  ########### (i) lead cis-pQTL ###############
                  cispqtl <- vroom(cis_pqtl_path)
                  cispqtl %<>% arrange(p_value)
                  cispqtl_lowestpval <- cispqtl[1,] # row with minimum p val (ensure that only one row is selected; otherwise, we'll bump into error afterwards)
                  
                  # rsid
                  cispqtl_lowestpval_rsid <- cispqtl_lowestpval %>% dplyr::select(rs_id) %>% as.character() # e.g., rs145144275
                  # chr
                  cispqtl_lowestpval_chr <- cispqtl_lowestpval %>% dplyr::select(chromosome) %>% as.numeric() # e.g.,19
                  # pos
                  cispqtl_lowestpval_pos <- cispqtl_lowestpval %>% dplyr::select(base_pair_location) %>% as.numeric()
                  
                  ########### (ii) corresponding pGWAS ###############
                  pgwas_path <- system(paste0('ls /home/csu/scratch/data/03.cis_pgwas/', cohort_dir, '/cis_pgwas_tss_1mb/', protname, '.txt.gz'), intern = T) # 1 mb only
                  #pgwas_path <- system(paste0('ls /scratch/richards/satoshi.yoshiji/24.meta_pQTL/01.2.gmetal_format/output/*', seqid, '*.txt.gz'), intern = T)
                  
                  #check all
                  print(protname)
                  print(cis_pqtl_path)
                  print(pgwas_path)
                  
                  pgwas <- vroom(pgwas_path)
                  pgwas_sel <- pgwas %>%
                    filter(chromosome == cispqtl_lowestpval_chr) %>%
                    filter(base_pair_location >= cispqtl_lowestpval_pos-500000) %>%
                    filter(base_pair_location <= cispqtl_lowestpval_pos+500000)
                  # give varbeta
                  pgwas_sel <- pgwas_sel %>% mutate(varbeta = standard_error^2)
                  
                  # remove variants w/o rsid
                  pgwas_rsid <- pgwas_sel %>%
                    filter(!is.na(rs_id))
                  # remove duplicated snp
                  pgwas_nodup <- pgwas_rsid %>%
                    dplyr::arrange(p_value) %>%
                    dplyr::filter(!duplicated(rs_id))
                  
                  # rename
                  pgwas_nodup %<>%
                    transmute(
                      chr = chromosome,
                      position = base_pair_location,
                      beta = beta,
                      varbeta = varbeta,
                      SE = standard_error,
                      ALT = effect_allele,
                      REF = other_allele,
                      MAF = effect_allele_frequency,
                      N = n,
                      pvalues = p_value,
                      snp = rs_id
                    ) 
                  
                  # rename to pgwas_nodup and add id column
                  pgwas_nodup$id <- paste('chr', pgwas_nodup$chr, pgwas_nodup$position, pgwas_nodup$REF, pgwas_nodup$ALT, sep = ":")
                  
                  # drop_na
                  #pgwas_nodup <- drop_na(pgwas_nodup)
                  pgwas_nodup <- pgwas_nodup[complete.cases(pgwas_nodup), ]
                  
                  # deal w/ extreme MAF and pval
                  pgwas_nodup$MAF[pgwas_nodup$MAF == 0] <- 1e-4
                  pgwas_nodup$MAF[pgwas_nodup$MAF == 1] <- 0.9999
                  
                  pgwas_nodup$pvalues[pgwas_nodup$pvalues == 0] <- 1e-300
                  pgwas_nodup$pvalues[pgwas_nodup$pvalues == 1] <- 0.9999
                  
                  # # make it a list
                  # npnt_list <- as.list(npnt_nodup)
                  # npnt_list$type <- 'quant'
                  # str(npnt_list)
                  #
                  # # check dataset
                  # check_dataset(npnt_list)
                  
                  #####################
                  # outcome GWAS
                  #####################
                  outcome <- vroom(ith_outcome_path)
                  outcome_rename0 <- outcome %>%
                    dplyr::rename(
                      chr = chromosome,
                      position = base_pair_location,
                      beta = beta,
                      SE = standard_error,
                      ALT = effect_allele,
                      REF = other_allele,
                      pvalues = p_value,
                      MAF = effect_allele_frequency,
                      #  samplesize = SS,
                      snp = rs_id)
                  
                  # keep SNPs that are present in pGWAS
                  outcome_rename <- outcome_rename0 %>% filter(snp %in% unlist(pgwas_nodup$snp))  # to reduce computational time, restrict SNPs here
                  
                  outcome_rename %<>% mutate(id = paste0('chr', chr, ':', position, ':', REF, ':', ALT)) #ALT = effect_allele, REF = other_allele,
                  
                  outcome_sel <- outcome_rename %>% dplyr::select(id, chr, position, beta, SE, ALT, REF, pvalues, MAF, snp)
                  
                  # give varbeta
                  outcome_sel$varbeta <- (outcome_sel$SE)^2
                  #outcome_sel %<>% dplyr::select(-SE)
                  
                  # remove variants w/o rsid
                  outcome_rsid <- outcome_sel %>%
                    filter(!is.na(snp))
                  
                  # remove duplicated id
                  outcome_nodup <- outcome_rsid %>%
                    dplyr::arrange(pvalues) %>%
                    dplyr::filter(!duplicated(snp))
                  
                  # remove NA
                  #outcome_nodup <- drop_na(outcome_nodup)
                  
                  # deal w/ extreme MAF and pval to avoid errors
                  outcome_nodup$MAF[outcome_nodup$MAF == 0] <- 1e-4
                  outcome_nodup$MAF[outcome_nodup$MAF == 1] <- 0.9999
                  
                  outcome_nodup$pvalues[outcome_nodup$pvalues == 0] <- 1e-300
                  outcome_nodup$pvalues[outcome_nodup$pvalues == 1] <- 0.9999
                  
                  ####################
                  # keep only common SNPs
                  # inner_join will not only keep rows with common SNPs but make their order same in both datasets
                  ####################
                  print('working on common')
                  
                  common <- inner_join(pgwas_nodup, outcome_nodup, by = 'snp', suffix = c('.pgwas', '.outcome'))
                  common %<>% rename(ALT = ALT.pgwas, REF = REF.pgwas) %>% select(-ALT.outcome) %>% select(-REF.outcome)
                  
                  #common %<>% dplyr::rename(chr.pgwas = chr.pgwas, position.pgwas = position)
                  
                  # If there is no common SNPs, write error message out and skip the process. otherwise the entire process will stop.
                          if(nrow(common) == 0){
                            #resname <- paste0(protname, '.', outcome_name, '.coloc.tsv')
                            #write.table('no common snp', file = resname, quote = F, row.names = F, sep ='\t')
                            print('no common snp')
                          } else {
                            print('common snps exist')
                            # pQTL
							              pgwas_common <- common %>% dplyr::select(chr.pgwas, position.pgwas, ALT, REF, beta.pgwas,  pvalues.pgwas, MAF.pgwas, varbeta.pgwas, SE.pgwas, snp) %>%
                              dplyr::rename(chr = chr.pgwas, position = position.pgwas, ALT = ALT, REF = REF, beta = beta.pgwas, pvalues = pvalues.pgwas, MAF = MAF.pgwas, varbeta = varbeta.pgwas, SE = SE.pgwas)
                            
                            # make sure to obtain chr and pos derives from pgwas (not outcome, which can be based on GRCh37/38)
                            outcome_common <- common %>% dplyr::select(chr.pgwas, position.pgwas, ALT, REF, beta.outcome, pvalues.outcome, MAF.outcome, varbeta.outcome, SE.outcome, snp) %>%
                              dplyr::rename(chr = chr.pgwas, position = position.pgwas, ALT = ALT, REF = REF, beta = beta.outcome, pvalues = pvalues.outcome, MAF = MAF.outcome, varbeta = varbeta.outcome, SE = SE.outcome)
                            
                            # # reorder them as requested by pwcoco
                            # pgwas_common_pwcoco <- pgwas_common %>% dplyr::select(snp, ALT, REF, MAF, beta, SE, pvalues)
                            # outcome_common_pwcoco <- outcome_common %>% dplyr::select(snp, ALT, REF, MAF, beta, SE, pvalues)
                            # 
                            # # save results (with chr, pos)
                            # pgwas_common_filename <- paste0('04.pwcoco/pwcoco_gwas/', protname, '.', outcome_name, '.common.500kb.tsv')
                            # write.table(pgwas_common, file = pgwas_common_filename, sep = '\t', quote = F, row.names = F)
                            # outcome_common_filename <- paste0('04.pwcoco/pwcoco_gwas/', outcome_name, '.', protname, '.common.500kb.tsv')
                            # write.table(outcome_common, file = outcome_common_filename, sep = '\t', quote = F, row.names = F)
                            
                            ######
                            # run pwcoco
                            ######
                            sumstats1 <- pgwas_common
                            sumstats2 <- outcome_common
                            
                            ##################
                            # pQTL
                            ##################
                            # sumstats1_list <- as.list(sumstats1)
                            # sumstats1_list$type <- 'quant'
                            # sumstats1_list$N <- cispqtl_lowestpval$n
                            # sumstats1_list$MAF <- sumstats1$MAF
                            # #sumstats1_list$sdY <- 1
                            # str(sumstats1_list)
                            # # check dataset
                            # check_dataset(sumstats1_list)
                            # plot_dataset(sumstats1_list)
                            # 
                            #pwcoco
                            sumstats1_pwcoco <- sumstats1 %>% select(snp, ALT, REF, MAF, beta, SE, pvalues)
                            sumstats1_pwcoco %<>% mutate(n = cispqtl_lowestpval$n)
                            sumstats1_name <- paste0(wd, protname, '/', protname, '--', ith_outcome_name, '/', protname, '--', ith_outcome_name, '--sumstats1.tsv')
                            write_tsv(sumstats1_pwcoco, file = sumstats1_name)
                            
                                    ####################
                                    # outcome has two pattern: quantitative or binary outcomes
                                    ####################
                                    # outcome has two patterns: quantitative or binary outcomes
                                    if (is.na(ith_outcome_casesize)) {
                                      # ## If ith_outcome_casesize is NA or an empty vector, outcome is quantitative.
                                      # ## Make it a list
                                      # sumstats2_list <- as.list(sumstats2)
                                      # sumstats2_list$type <- 'quant'
                                      # sumstats2_list$N <- ith_outcome_samplesize
                                      # sumstats2_list$MAF <- sumstats2$MAF
                                      # # sumstats2_list$s <- casesize/samplesize
                                      # str(sumstats2_list)
                                      # 
                                      # # check dataset
                                      # check_dataset(sumstats2_list)
                                      # plot_dataset(sumstats2_list)
                                      
                                      #pwcoco
                                      sumstats2_pwcoco <- sumstats2 %>% select(snp, ALT, REF, MAF, beta, SE, pvalues)
                                      sumstats2_pwcoco %<>% mutate(n = cispqtl_lowestpval$n) 
                                      sumstats2_name <- paste0(wd, protname, '/', protname, '--', ith_outcome_name, '/', protname, '-', ith_outcome_name, '--sumstats2.tsv')
                                      write_tsv(sumstats2_pwcoco, file = sumstats2_name)
                                      
                                    } else {
                                      # # When ith_outcome_casesize is not NA and not an empty vector, the outcome is binary.
                                      # # # Make it a list
                                      # sumstats2_list <- as.list(sumstats2)
                                      # sumstats2_list$type <- 'cc'
                                      # sumstats2_list$N <- ith_outcome_samplesize
                                      # sumstats2_list$MAF <- sumstats2$MAF
                                      # sumstats2_list$s <- ith_outcome_casesize / ith_outcome_samplesize
                                      # str(sumstats2_list)
                                      # 
                                      # # check dataset
                                      # check_dataset(sumstats2_list)
                                      # plot_dataset(sumstats2_list)
                                      # 
                                      #pwcoco
                                      sumstats2_pwcoco <- sumstats2 %>% select(snp, ALT, REF, MAF, beta, SE, pvalues)
                                      sumstats2_pwcoco %<>% mutate(n = ith_outcome_samplesize) %>% mutate(case = ith_outcome_casesize)
                                      sumstats2_name <- paste0(wd, protname, '/', protname, '--', ith_outcome_name, '/', protname, '--', ith_outcome_name, '--sumstats2.tsv')
                                      write_tsv(sumstats2_pwcoco, file = sumstats2_name)
                                    }
                            
                            # Now both dataset 1 and 2 are ready, so let's perform pwcoco (instead of coloc)
                            # pwcoco
                            #system(paste0('/home/csu/scratch/tools/pwcoco/build/pwcoco --bfile /scratch/richards/restricted/ukb-general/scratch/genetic-data/genome/imputed.v3/plink.bycroftRand50k/',  cispqtl_lowestpval_chr , ' --sum_stats1 ',  sumstats1_name,  ' --sum_stats2 ', sumstats2_name, ' --top_snp 5 --maf 0.01 --out ', protname, '--', ith_outcome_name))
                            system(paste0('/home/csu/projects/rrg-vmooser/csu/software/pwcoco/build/pwcoco --bfile /home/csu/scratch/data/02.ukb50k_maf0.01/',  cispqtl_lowestpval_chr , ' --sum_stats1 ',  sumstats1_name,  ' --sum_stats2 ', sumstats2_name, ' --top_snp 5 --threads 4 --maf 0.01 --out ', protname, '--', ith_outcome_name))
                            print('pwcoco completed')
                            
                            # # # colocalizatiaon
                            # coloc_res <- coloc.abf(
                            #   dataset1 = sumstats1_list,
                            #   dataset2 = sumstats2_list
                            # )

                            # # write results
                            # system(paste0('mkdir -p ', protname))
                            # resname <- paste0(protname, '/', protname, '-', ith_outcome_name, '-coloc.tsv')
                            # write_tsv(as.data.frame(coloc_res$summary) %>% rownames_to_column(), file = resname)
                            
                            } #  close if else for if(nrow(common) == 0){
                  
                  } # closes if(file.exists(cis_pqtl_path)){ 
            } else { print(paste0('skip pwcoco for ', protname, '--', ith_outcome_name))}# closes if(mr_res_waldivw_pval < 0.05){                
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}) # tryCatch({ # try catch outcome error from for-loop of the outcome_path 
}  # close if
 
