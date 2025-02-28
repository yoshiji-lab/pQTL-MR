####
######### Given an rsid, find proxies for it 
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


args <- commandArgs(trailingOnly=TRUE) #input UNZIPPED file
rsid <- args[1]  # rsid to find proxy for
#rsid <- "rs1558526"

# 2. Read in file with all variants to use as exposure.
all_exposure_file <- args[2]  
#all_exposure_file <- "/scratch/richards/chen-yang.su/01.pQTL/01.aric/03.cis_pQTL/all_strict_v2g_cis_res.tsv"
Ins_df <- fread(all_exposure_file)
cat("Path to all exposures:", all_exposure_file, "\n")

# 3. Output dir to store raw proxies found
outdir <- args[3]
#outdir <- '/scratch/richards/chen-yang.su/01.pQTL/01.aric/05.MR_strict_v2g/0.find_proxies/found_proxies/'
system(paste0("mkdir -p ", outdir))
cat("Output directory:", outdir, '\n')

ref_panel <- args[4]  # path to plink files for reference panel
#ref_panel <- "/scratch/richards/chen-yang.su/data/02.ukb50k_maf0.01/"

    
# Get LD proxy
# https://mrcieu.github.io/gwasvcf/reference/get_ld_proxies.html
set_plink(path = "/scratch/richards/chen-yang.su/software/plink")  # else get_ld_proxy() may give error


chr <- Ins_df[Ins_df$rs_id == rsid, ]$chromosome  # extract chromosome
bfile <- paste0(ref_panel, chr)

raw_proxies <- get_ld_proxies(
  rsid,  # single rsid
  bfile,  # plink file corresponding to chromosome of the variant
  tag_kb = 5000,
  tag_nsnp = 5000,
  tag_r2 = 0.8,
  threads = 4,
  out = tempfile()
)
raw_proxies$R_sq <- raw_proxies$R^2  # compute rsq ourselves since function only returns R
  
# Example of first 2 rows of output (found 30 proxies) when querying rs564382789 NOTE: *_A are columns of original variant, *_B are proxy variant columns

# CHR_A     BP_A SNP_A       MAF_A CHR_B     BP_B SNP_B      PHASE MAF_B     R A1    B1    A2    B2     R_sq
# <int>    <int> <chr>       <dbl> <int>    <int> <chr>      <chr> <dbl> <dbl> <chr> <chr> <chr> <chr> <dbl>
#     19 39419059 rs564382789 0.401    19 39417652 rs10425217 CTTC  0.390 0.985 C     T     T     C     0.970
#     19 39419059 rs564382789 0.401    19 39417440 rs7508539  CGTA  0.390 0.985 C     G     T     A     0.970
  
  
# Filter out proxies where MAF > 0.42 (equivalently, retain proxies with MAF <= 0.42)
raw_proxies <- raw_proxies[raw_proxies$MAF_A <= 0.42 & raw_proxies$MAF_B <= 0.42, ]



outname <- paste0(outdir,"/", rsid, "_raw_proxies.txt")
write_tsv(raw_proxies, file=outname)


end_time <- Sys.time()

print(paste0("Time taken: ", end_time - start_time))