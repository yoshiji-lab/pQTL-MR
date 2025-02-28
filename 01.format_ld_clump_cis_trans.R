# Format to GWAS catalog columns, Label cis/trans, LD clumping

start_time <- Sys.time()

.libPaths(c("/home/richards/chen-yang.su/R/x86_64-pc-linux-gnu-library/4.1",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Compiler/gcc9/r-bundle-bioconductor/3.14",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx512/Core/r/4.1.2/lib64/R/library"))

library(magrittr)
library(dplyr)
library(tidyverse)
library(readxl)
library(data.table)
library(tidyr)
library(stringr)
library(ieugwasr)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

path <- args[1]  # "/scratch/richards/satoshi.yoshiji/10.Proteome/ARIC/EA/SeqId_10000_28.PHENO1.glm.linear.gz"
# path <- "/scratch/richards/satoshi.yoshiji/10.Proteome/ARIC/EA/SeqId_10000_28.PHENO1.glm.linear.gz"

protein <- basename(path)  # SeqId_10000_28.PHENO1.glm.linear.gz
prot_suffix <- str_replace(protein, "SeqId_", "")
seqid <- paste0(strsplit(prot_suffix, "\\.")[[1]][1])  # 10000_28
print(seqid)


# 1. Retain 4728 analyzed proteins only ----

file <- "/scratch/richards/chen-yang.su/data/01.1.SL_prot_TSS_4687/prot_tss_somalogic.v1.25Dec2023.xlsx"
dat <- read_excel(file, sheet = "Sheet1", skip = 1)
dat <- dat[dat$`Included in analysis` == "Yes", ]
# dat <- dat %>% select(SeqId, Chromosome, TSS)
dat['TSS_lb'] <- dat$TSS - 500000
dat['TSS_ub'] <- dat$TSS + 500000
dat <- dat %>% mutate(TSS_lb = if_else(TSS_lb < 0, 0, TSS_lb))  # set negative position to 0

if (seqid %in% dat$SeqId) {  # only process proteins we analyzed (4685)
  
  # 2. Retain variants on autosomes only ----
  
  
  df <- fread(path, sep = '\t', fill=TRUE)
  head(df)
  
  # >   head(df)
  # #CHROM      POS        ID REF ALT A1  A1_FREQ TEST OBS_CT        BETA        SE     T_STAT        P ERRCODE
  # 1:     22 24712800 rs1074487   A   G  G 0.204561  ADD   7213  0.02601580 0.0204240  1.2737800 0.202781       .
  # 2:     22 24713060 rs2007382   G   A  A 0.204631  ADD   7213  0.02513390 0.0204234  1.2306400 0.218498       .
  # 3:     22 24713102 rs2007386   T   C  C 0.204631  ADD   7213  0.02513390 0.0204234  1.2306400 0.218498       .
  # 4:     22 24713953  rs131483   C   T  T 0.306530  ADD   7213 -0.01792500 0.0181226 -0.9890960 0.322649       .
  # 5:     22 24714428  rs140347   T   G  T 0.496603  ADD   7213 -0.00204134 0.0167407 -0.1219390 0.902951       .
  # 6:     22 24714689  rs131482   T   G  T 0.494177  ADD   7213  0.00118302 0.0167352  0.0706908 0.943646       .
  
  ##################### speed up everything ######################
  str(df)
  df <- df[df$P < 5e-8, ]
  ################################################################
  
  unique(df$`#CHROM`)  
  
  # >   unique(df$`#CHROM`)  
  # [1] 22
  
  autosomes <- seq(1,22)
  
  # > autosomes
  # [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22
  
  autosomes_df <- df[df$`#CHROM` %in% autosomes, ]
  unique(autosomes_df$`#CHROM`)  
  
  
  # 3. Filter MAF > 0.01  ----
  autosomes_df <- autosomes_df[autosomes_df$A1_FREQ > 0.01 & autosomes_df$A1_FREQ < 0.99, ]

  
  # 4. Label cis/trans ----
  if (nrow(autosomes_df) > 0) {
    
    # get gene reference chr, TSS lower and upper bound
    gene_chr <- dat[dat$SeqId == seqid,]$Chromosome
    gene_tss_lb <- as.integer(dat[dat$SeqId == seqid,]$TSS_lb)
    gene_tss_up <- as.integer(dat[dat$SeqId == seqid,]$TSS_ub)

    autosomes_df$cis_trans <- ifelse( (autosomes_df$POS >= gene_tss_lb) & (autosomes_df$POS <= gene_tss_up) & (autosomes_df$`#CHROM` == gene_chr),
                                      "cis",
                                      "trans")
    # print(nrow(autosomes_df[autosomes_df$cis_trans == "cis", ]))
    # print(nrow(autosomes_df[autosomes_df$cis_trans == "trans", ]))
    # print(nrow(autosomes_df) == nrow(autosomes_df[autosomes_df$cis_trans == "cis", ]) + nrow(autosomes_df[autosomes_df$cis_trans == "trans", ]))




    # format columns to gwas catalog format ----

    # Convert 'chromosome' and 'base_pair_location' columns to numeric
    autosomes_df$chromosome <- as.numeric(autosomes_df$`#CHROM`, na.rm = TRUE)
    autosomes_df$base_pair_location <- as.numeric(autosomes_df$POS, na.rm = TRUE)

    # select effect allele and other allele
    autosomes_df$effect_allele <- autosomes_df$A1
    autosomes_df$other_allele <- ifelse(autosomes_df$A1 == autosomes_df$ALT, autosomes_df$REF, autosomes_df$ALT)


    # Capitalize 'effect_allele' and 'other_allele' columns
    autosomes_df$effect_allele <- toupper(autosomes_df$effect_allele)
    autosomes_df$other_allele <- toupper(autosomes_df$other_allele)

    autosomes_df$beta <- as.numeric(autosomes_df$BETA)
    autosomes_df$standard_error <- as.numeric(autosomes_df$SE)
    autosomes_df$effect_allele_frequency <- as.numeric(autosomes_df$A1_FREQ)

    autosomes_df$p_value <- as.numeric(autosomes_df$P)
    autosomes_df$minus_log10p <- as.numeric(-log10(autosomes_df$p_value))
    # z_score
    autosomes_df$z_score <- as.numeric(autosomes_df$beta / autosomes_df$standard_error)

    # Create the variant_id column
    autosomes_df$variant_id <- paste0(
      as.character(autosomes_df$chromosome), "_",
      as.character(autosomes_df$base_pair_location), "_",
      as.character(autosomes_df$other_allele), "_",
      as.character(autosomes_df$effect_allele)
    )

    autosomes_df$rs_id <- as.character(autosomes_df$ID)

    # Calculate ci_lower and ci_upper columns
    autosomes_df$ci_lower <- as.numeric(autosomes_df$beta - 1.96 * autosomes_df$standard_error)
    autosomes_df$ci_upper <- as.numeric(autosomes_df$beta + 1.96 * autosomes_df$standard_error)

    autosomes_df$n <- as.numeric(autosomes_df$OBS_CT)

    # Reorder the columns
    columns_order <- c(
      "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
      "standard_error", "effect_allele_frequency", "p_value", "minus_log10p",
      "z_score", #"direction",
      "variant_id", "rs_id", "ci_lower", "ci_upper", "n",
      "cis_trans"
    )
    autosomes_df %<>% dplyr::select(columns_order)

    # Sort the merged file by chromosome and base_pair_location
    autosomes_df <- autosomes_df[order(autosomes_df$chromosome, autosomes_df$base_pair_location), ]







    # LD clumping ----


    # Only clump variants with p < 5e-8
    autosomes_df <- autosomes_df[autosomes_df$p_value < 5e-8, ]

    # add new pval column (needs to be called pval for )
    autosomes_df$pval <- as.numeric(autosomes_df$p_value)
    typeof(autosomes_df$pval)

    max(autosomes_df$minus_log10p)  # here you see minus log 10 p actually shows min p value is not 0
    min(autosomes_df$pval)  # here you see min pval is 0

    # Rename for ld_clump_local function
    autosomes_df %<>% dplyr::rename("rsid" = "rs_id")


    # Initialize empty dataframe to store clumping results for all chromosomes
    header <- colnames(autosomes_df)
    new_df <- data.frame(matrix(ncol = length(header), nrow = 0))
    colnames(new_df) <- header
    head(new_df)

    chr <- as.integer(autosomes_df$chromosome[1])
    for (i in chr:chr) {
      tryCatch({
        df_clumped <- ld_clump_local(
          dat = autosomes_df,
          bfile = paste('/scratch/richards/chen-yang.su/data/02.ukb50k_maf0.01/', i, sep = ""),  # UKB 50k with MAF > 0.01
          plink_bin = '/scratch/richards/chen-yang.su/software/plink',
          clump_kb = 1000,
          clump_p = 5e-8,
          clump_r2 = 0.001
        )
        new_df <- rbind(new_df, df_clumped)
        print(new_df)

      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


    }

    # drop pval column since don't need it anymore
    new_df$pval <- NULL

    # rename back to original column names
    new_df <- new_df %>% dplyr::rename("rs_id" = "rsid")



    #################
    # prepare an output file name
    #################

    # protein symbol
    sldict <- fread('/scratch/richards/chen-yang.su/data/01.SL_proteins/SL_proteins.txt', sep = '\t')
    sldict <- sldict[sldict$Organism == 'Human', ]  # Filter the data for 'Organism' column equal to 'Human'
    sldict2 <- sldict[, c('SeqId', 'EntrezGeneSymbol')]  # Select columns 'SeqId' and 'EntrezGeneSymbol'
    sldict3 <- sldict2  # Create a copy of sldict2
    sldict3$SeqId <- gsub('-', '_', sldict2$SeqId)  # Modify the 'SeqId' column by replacing '-' with '_'
    proteinsymbol <- sldict3[sldict3$SeqId == seqid, 'EntrezGeneSymbol']  # Filter sldict3 to get 'EntrezGeneSymbol' for a specific 'SeqId'

    # output_path
    directory_path <- "/scratch/richards/chen-yang.su/01.pQTL/01.aric/ld_clumped_cis_trans/"
    output_path <- paste0(directory_path, seqid, '.', proteinsymbol, '.txt.gz')

    # Write the merged file to a compressed TSV file
    write_tsv(new_df, output_path)

    
  } else {
    print("no variants with MAF < 0.01")
    system('touch ./aric_unanalyzed_prots.txt')
    system(paste0('echo ', path, ' >> ./aric_unanalyzed_prots.txt'))
  }
  
} else {
  print("not in list of 4687 proteins")
  system('touch ./aric_unanalyzed_prots.txt')
  system(paste0('echo ', path, ' >> ./aric_unanalyzed_prots.txt'))
}

end_time <- Sys.time()

print(paste0("Time taken: ", end_time - start_time))
