# pQTL-MR

## Overview

This repository contains scripts and resources for conducting Mendelian Randomization (MR) analyses using protein Quantitative Trait Loci (pQTL) data. The primary objective is to investigate the causal relationships between protein biomarkers and various complex traits or diseases.

## Repository Structure

```
pQTL-MR/
│-- README.md                # Project documentation
│-- scripts/                 # Utility scripts
│-- 0.write_proteins_to_file.sh  # Script to list proteins
│-- 01.format_ld_clump_cis_trans.R # Formatting and LD clumping for cis/trans pQTLs
│-- 01.launch.sh             # Launch script for analyses
│-- 02.compile_cis_trans.R   # Compilation of cis/trans pQTL results
│-- 03.cis_pQTL/             # Analysis of cis-pQTLs
│-- 05.MR_strict_v2g/        # MR analyses with strict variant-to-gene mapping
│-- 06.coloc/                # Colocalization analyses
│-- 07.summarize/            # Summarization of results
```

## Getting Started

### Prerequisites

Ensure you have the following software and packages installed:

- **R** (version 4.0 or higher) with the following packages:
  - `TwoSampleMR`
  - `data.table`
  - `dplyr`
  - `coloc`
- **PLINK** (for linkage disequilibrium clumping)
- **Shell environment** (bash) for executing shell scripts

### Installation

1. **Clone the repository:**

   ```bash
   git clone https://github.com/yoshiji-lab/pQTL-MR.git
   cd pQTL-MR
   ```

2. **Install R packages:**

   Open an R session and run:

   ```R
   install.packages(c("devtools", "data.table", "dplyr", "coloc"))
   devtools::install_github("MRCIEU/TwoSampleMR")
   ```

## Workflow

1. **Listing Proteins:**

   Use the `0.write_proteins_to_file.sh` script to generate a list of proteins for analysis.

   ```bash
   bash 0.write_proteins_to_file.sh
   ```

2. **Formatting and LD Clumping:**

   The `01.format_ld_clump_cis_trans.R` script formats pQTL data and performs linkage disequilibrium (LD) clumping to distinguish between cis and trans pQTLs.

   ```bash
   Rscript 01.format_ld_clump_cis_trans.R
   ```

3. **Compiling cis/trans pQTL Results:**

   Compile the results of cis and trans pQTL analyses using the `02.compile_cis_trans.R` script.

   ```bash
   Rscript 02.compile_cis_trans.R
   ```
4. **Compiling strict V2G results:**

   Create strict V2G cis-pQTLs using scripts in the `03.cis_pQTL` folder

5. **Mendelian Randomization Analysis:**

   Conduct MR analyses using the scripts in the `05.MR_strict_v2g/` directory. The scripts here use the utility scripts in the `scripts` directory

6. **Colocalization Analysis:**

   Perform colocalization analyses to assess whether the same genetic variant influences both the protein levels and the trait of interest. Use the scripts in the `06.coloc/` directory.

7. **Summarizing Results:**

   Summarize all findings using the scripts in the `07.summarize/` directory to generate comprehensive reports and visualizations.

## Results

The results from each analysis step are stored in their respective directories. Key findings, plots, and summary statistics are organized for easy interpretation and further exploration.

## Contributing

Contributions to this project are welcome. To contribute:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Commit your changes with descriptive messages.
4. Push your branch and open a Pull Request.

Please ensure that your code adheres to the project's coding standards and includes appropriate tests.

## License

This project is licensed under the MIT License. See the `LICENSE` file for more details.

## Contact

For questions, suggestions, or collaborations, please contact the project maintainer:

- **Yoshiji Lab** (Email: chen-yang.su@mail.mcgill.ca)

We appreciate your interest in our work and look forward to potential collaborations.
