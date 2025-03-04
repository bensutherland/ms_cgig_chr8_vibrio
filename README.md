# ms_cgig_chr8_vibrio
Code repository to accompany response of CHR8 oyster families to _Vibrio aestuarianus_, which currently includes the following:     
- [1. ddRADseq of Vibrio challenge survivors (CHR8 families)](#01-ocv23-analysis)      
- [2. rhAmp assay for CHR8 genotypes in Vibrio challenge survivors (CHR8 families)](#02-OCV23-rhAmp-analysis)     

#### Requirements ####
[simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)    
[stacks_workflow](https://github.com/enormandeau/stacks_workflow)       
[GEMMA](https://github.com/genetics-statistics/GEMMA/tree/master)         

### 01. OCV23 analysis ###
#### a. Set up ####
Use `stacks_workflow` to analyze the data. Clone the repo in the same parent folder as the present repo, change directory into `stacks_workflow` and run all commands from `stacks_workflow`.       

```
git clone https://github.com/enormandeau/stacks_workflow.git
cd stacks_workflow

# Copy all raw data into 02-raw
cp -l ../00_raw_data/*.fastq.gz ./02-raw/
```

#### b. Check raw data ####
View raw data with fastqc and multiqc:      
```
mkdir 02-raw/fastqc_raw    
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 10   
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw   
```

#### c. Clean the raw data ####
Use cutadapt to clean the raw reads:     
```
# Prepare a utility text file with all filenames
./00-scripts/00_prepare_lane_info.sh

# Run cutadapt to trim reads based on quality, and remove adapters
./00-scripts/01_cutadapt.sh 8 
```

#### d. Check the trimmed data ####
View trimmed data with fastqc and multiqc:     
```
mkdir 02-raw/trimmed/fastqc_trimmed/    
fastqc -t 5 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc_trimmed/
multiqc -o 02-raw/trimmed/fastqc_trimmed/ 02-raw/trimmed/fastqc_trimmed       
```

#### e. Prepare sample metadata ####
Prepare a sample metadata file using necessary inputs as per `stacks_workflow` README.       


#### f. Demultiplex ####
Detect cut sites and barcodes to de-multiplex and truncate reads to 80 bp with `process_radtags` in parallel:     
`00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 nsiI mspI 8`    


#### g. Rename samples ####
`./00-scripts/03_rename_samples.sh`      


View trimmed data with fastqc and multiqc:     
```
mkdir 04-all_samples/fastqc/    
fastqc -t 48 04-all_samples/*.fq.gz -o 04-all_samples/fastqc/
multiqc -o 04-all_samples/fastqc/ 04-all_samples/fastqc/       
```

#### h. Alignment ####
Index the genome:      
`bwa index GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.gz`               

Align samples against the genome:     
`./00-scripts/bwa_mem_align_reads.sh 32`      


#### i. Remove control samples and low quality samples ####
We will do an analysis with all samples, except control samples, to see what we should use as a low read depth cutoff per sample.    

Move any control or unused samples to `04-all_samples/removed_samples` (newly created directory)

#### j. Genotype the samples ####
Prepare the population map file:
`./00-scripts/04_prepare_population_map.sh`    
Note: manually remove the blank population, since these files are no longer in the dataset.    

Genotype the samples:     
`00_scripts/stacks2_gstacks_reference.sh`     

Set flags and run the populations module to filter the data.     

```
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 6 -r 0.7 \
    --min-maf 0.01 \
    --ordered-export --plink --hwe --write-single-snp --fasta-loci --vcf

# Move the output to a new folder to store single-snp data
mkdir 05-stacks/popn_out_single_snp/
mv 05-stacks/populations.* 05-stacks/popn_out_single_snp/ 
```

#### k. Convert output plink data ####
Use the single-SNP per locus data.      
Convert plink files to a useable format for adegenet:        
`plink --ped 05-stacks/popn_out_single_snp/populations.plink.ped --map 05-stacks/popn_out_single_snp/populations.plink.map --maf 0.01 --recode A --allow-extra-chr --out 05-stacks/popn_out_single_snp/populations_single_snp`      

Note: these plink files are available from [FigShare](https://doi.org/10.6084/m9.figshare.26524321.v1)       

#### l. Population genetic analysis ####
Read the data into R using the script `01_scripts/01_import_plink_to_genind.R`           

#### m. Genome-wide association study (GWAS) ####
From `stacks_workflow` (or FigShare link above), copy the output plink map file, the VCF file, and the sample interpretaion file into `ms_cgig_chr8/02_input_data/`.    

Open `01_scripts/GWAS.R` in Rstudio interactively.    
Inputs:    
- `sample_interp_2024-07-18.csv`   
- non-imputed VCF `populations.snps.vcf`   
- imputed VCF `populations.snps.imputed.vcf`     

The following steps are taken:   
i. Changes linkage group (LG) RefSeq genome annotations to corresponding chromosome annotations and removes all other contigs in the RefSeq genome    
ii. Imputes missing genotypes within each family independently using mean imputation    

This Rscript will export the following files into `03_results`:    
- gwasanno.txt (annotation of the SNPs)
- gwascovar.txt (family covariate)
- gwasgeno.txt (genotype matrix)
- gwaspheno.txt (binary dead/alive vector)
- gwaspheno2.txt (numeric days to death vector)

Use GEMMA to run a GWAS using these files:     

```
# Run in command-line
# change directory into 03_results

# Calculate kinship matrix
gemma -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_allfam_pheno_dead_alive

# Run the GWAS analysis with the family covariate and binary phenotype
gemma -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam_pheno_dead_alive.cXX.txt  -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_covar_dead_alive

# Run the GWAS analysis with the family covariate and numeric day-to-death phenotype
# Calculate kinship matrix
gemma -g gwasgeno.txt -p gwaspheno2.txt -gk -maf 0.05 -o gwas_allfam_pheno_day_to_death

# Run the GWAS analysis with the family covariate and day-to-death phenotype
gemma -g gwasgeno.txt -p gwaspheno2.txt -k output/gwas_allfam_pheno_day_to_death.cXX.txt -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_covar_pheno_day_to_death

```


iii. Runs the following GWAS models using GEMMA (MAF threshold of 0.05):       

With all families (4 MBP families + VIU family) together:    
Phenotype (survival; 1-alive, 0-dead) = SNP (fixed) + Population structure (random; covariance matrix = G [genomic relationship] matrix) + error (random; covariance matrix = identity matrix)    
Phenotype (survival) = SNP (fixed) + Family effect (fixed) + Population structure (random; G matrix) + error (random; identity matrix)    

With all 4 MBP families together:    
Phenotype (survival) = SNP (fixed) + Population structure (random; G matrix) + error (random; identity matrix)    
Phenotype (survival) = SNP (fixed) + Family effect (fixed) + Population structure (random; G matrix) + error (random; identity matrix)    

With all families independently:    
Phenotype (survival) = SNP (fixed) + error (random; identity matrix)  

The SNPs are coded 0, 1, 2 to correspond to reference homozygote, heterozygote, and alternative homozygote genotypes, respectively. Therefore, the models assume the allele are additive (rather than exhibiting dominance, for example). 

iv. Plots GWAS results as Manhattan plots


### 02. OCV23 rhAmp analysis ###
Input files are csv files that are raw output from the genotyping platform (i.e., CFX96 instrument).     
Copy all csv files with rhAmp results into `02_input_data`. The following column names are required, and should be standard in qPCR machine output:    
`Well`, `Fluor`, `Content`, `Cq`    

Copy the file `OSU_CHR8_VC_Mapping_Family114-117_2024-01-21_MFrhAmpDNAid.txt` (available from FigShare (#TODO)) into `00_archive`. This will be used to connect the well and plate ID to the sample ID.    

Run the Rscript `01_scripts/rhamp_assay_analysis.R` interactively to do the following:    
- format as needed.    
- #todo: update this section

#### Plotting mortality by genotype ####
Put `mort_prop_by_fam.csv` and `sample_day_of_death_and_DNA_ID.csv` in `00_archive`.     
(#TODO: need to include these in FigShare repo)

After running the above script, you should have a file called `03_results/rhAmp_per_sample_genotype_summary.csv` that will be used by the plotting script.     

Use the following script interactively to plot the results:     
`01_scripts/rhamp_06_plot_results.R`     

