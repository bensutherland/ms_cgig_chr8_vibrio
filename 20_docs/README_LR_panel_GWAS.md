# ms_cgig_chr8 amplicon-based low-resolution GWAS for OsHV-1 survivorship

### 01. Getting started ###
Data inputs:    
- BCF file from `amplitools` mhap workflow for genotyping, and through `impute_workflow` to isolate offspring only and rename;    
- phenotype file    

Put BCF file in `simple_pop_stats/02_input_data/`.      
Put pheno file in `simple_pop_stats/00_archive/`.    

Use bcftools view to convert the BCF file from BCF format to VCF format.    


### 02. Inspect and filter VCF, then prepare inputs for GEMMA analysis ###
Use script `01_scripts/chr8_oshv1_amp_02_vcf_to_gemma.R`       
Note: will have option to conduct filtering for genotyping rate of individuals, and to conduct mean imputation.   


### 03. Gemma analysis ###
Put GEMMA outputs from above, including geno, pheno(s), covariate into `03_results` of `ms_cgig_chr8`, then:     
```
# Change into the 03_results dir
cd 03_results

# Calculate kinship based on prepared GEMMA files
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam

# Calculate GWAS association values
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar
```

### 04. Plot results ###
Produce Manhattan plots from the GEMMA output by setting the filename and interactively running `01_scripts/chr8_oshv1_amp_03_gemma_results.R`.    

[Back to main README](https://github.com/bensutherland/ms_cgig_chr8)

