.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))
library(dplyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(readxl)
#library(stringr)

#############
# prepare files
#############
# format seqid 
# dictionary file link seqid and EntrezGeneSymbol
sl <- fread('/scratch/richards/chen-yang.su/data/01.SL_proteins/SL_proteins.txt')
sl %<>% separate(SeqId, into = c('seqid1', 'seqid2'), sep = '-') %>% unite(seqid, c('seqid1', 'seqid2'), sep = '_')
sl %<>% dplyr::select(seqid, EntrezGeneSymbol) %>% filter(!is.na(EntrezGeneSymbol)) %>% filter(EntrezGeneSymbol != '')

# mapping file for ENSG and gene
ensg_dic <- read_excel('/scratch/richards/chen-yang.su/01.meta_pQTL/02.3.cis_pQTL/dynamicpqtl.supptable.v2.01Aug2023.xlsx', sheet = 'ST6', skip = 1)
ensg_dic %<>% dplyr::select(SeqId, 'Ensembl.Gene.ID') %>% dplyr::rename(gene_id = Ensembl.Gene.ID, seqid = SeqId)

# combine sl and ensg_dict
sl_ensg <- left_join(sl, ensg_dic, by = 'seqid')
# collapse them to unique proteins (e.g., two seqids for COL6A3 to one row)
sl_ensg %<>% dplyr::select(-seqid) %>% filter(!duplicated(EntrezGeneSymbol))

###############
# read cis-pQTL
###############
# path to the strict cis-pQTL file
args <- commandArgs(trailingOnly=TRUE)
cis_path <- args[1]
# cis_path <- '/scratch/richards/satoshi.yoshiji/24.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.COL6A3.11196_31.tsv'
# cis_path <- '/scratch/richards/chen-yang.su/01.meta_pQTL/02.3.cis_pQTL/all_strict_cis/cis.A1BG.16561_9.tsv'
cis_protname <- str_split(basename(cis_path), '[.]')[[1]][2]
cis_seqid <- str_split(basename(cis_path), '[.]')[[1]][3]
cis <- fread(cis_path)
# add ENSG
cis_ensg <- sl_ensg %>% filter(EntrezGeneSymbol == cis_protname) %>% dplyr::select(gene_id) %>% as.character()
cis %<>% mutate(gene_id = cis_ensg)

############
# prepare for further filtering with V2G
############
index <- fread(paste0('/scratch/richards/satoshi.yoshiji/database/opentarget/V2G/v2g_index_perchr/chr',  cis$chromosome[1], '.v2g_index.tsv.gz'))
index_select <- index %>% filter(rs_id %in% cis$rs_id) # filtering beforehand can speed up the process

scored <- fread(paste0('/scratch/richards/satoshi.yoshiji/database/opentarget/V2G/v2g_scored_perchr/chr', cis$chromosome[1], '.v2g_scored.tsv.gz'))
scored_select <- scored %>% dplyr::select(-source_score_list) %>% 
  filter(chr_id  == cis$chromosome[1]) %>% filter(position %in% index_select$position) # filtering beforehand can speed up the process

v2g <- left_join(index_select, scored_select, by = c('chr_id', 'position', 'ref_allele', 'alt_allele'))
v2g %<>% group_by(rs_id) %>% arrange(desc(overall_score)) %>% filter(!duplicated(gene_id)) # only retain genes with top V2G for each variant

# combine with sl_ensg
v2g_sl_ensg <- left_join(v2g, sl_ensg, by = 'gene_id')

###############
# loop through all variants for the protein of interest
###############
cis_v2g <- data.frame()
cis_v2g_pass <- data.frame()

for(i in 1:nrow(cis)){
  #i <- 1
  ith_variant <- cis[i,]
  # annotate the gene/seqid with top V2G score and pass/fail (pass = v2g highest for the target seqid)
  ith_v2g <- v2g_sl_ensg  %>% filter(rs_id == ith_variant$rs_id) %>% arrange(desc(overall_score)) %>% dplyr::select(-chr_id, -position, -ref_allele, -alt_allele) 
  ith_v2g_top <- ith_v2g[1,] %>% dplyr::rename(topv2g_gene_id = gene_id, topv2g_overall_score = overall_score, topv2g_EntrezGeneSymbol = EntrezGeneSymbol)
  
  ith_variant_v2g <- left_join(ith_variant, ith_v2g_top, by  ='rs_id')
  ith_variant_v2g %<>% mutate(v2g_filter = case_when(gene_id == topv2g_gene_id ~ 'pass', # is gene_id the same as topv2g_gene_id?
                                                    TRUE ~ 'fail'))
  cis_v2g <- rbind(cis_v2g, ith_variant_v2g)
  
  # filter
  ith_variant_v2g_pass <- ith_variant_v2g %>% filter(v2g_filter == 'pass')
  if(nrow(ith_variant_v2g_pass) != 0){
    cis_v2g_pass <- rbind(cis_v2g_pass, ith_variant_v2g_pass)
  }
} # close for loop: for(i in 1:nrow(cis)){

############
# save results include both V2G-passing and failing proteins
filename <- paste0('cis.', cis_protname, '.', cis_seqid, '.tsv')
system('mkdir -p /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_annotated_cis/')
write_tsv(cis_v2g, file = paste0('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_annotated_cis/',   filename))

# save results for V2G-passing variants only
filename <- paste0('cis.', cis_protname, '.', cis_seqid, '.tsv')
system('mkdir -p /scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis')
if(nrow(cis_v2g_pass) != 0){
write_tsv(cis_v2g_pass, file = paste0('/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis/',   filename))
}

