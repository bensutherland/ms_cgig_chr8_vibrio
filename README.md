# ms_cgig_chr8_vibrio
Code repository to accompany response of CHR8 oyster families to _Vibrio aestuarianus_, which currently includes the following:     
- [1. ddRADseq of Vibrio challenge survivors (CHR8 families)](#01-ocv23-analysis)      
- [2. rhAmp assay for CHR8 genotypes in Vibrio challenge survivors (CHR8 families)](#02-OCV23-rhAmp-analysis)     

#### Requirements ####
[simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)    
[stacks_workflow](https://github.com/enormandeau/stacks_workflow)       
[GEMMA](https://github.com/genetics-statistics/GEMMA/tree/master)         

#### Input materials ####
Metadata and intermediate files are available on FigShare for the following:    
- ddRADseq analysis [FigShare](https://doi.org/10.6084/m9.figshare.26524321.v1)        
- rhAmp assay analysis [FigShare](https://doi.org/10.6084/m9.figshare.26515438.v1)     

#### Citation ####
Please see the following publication:     
Surry LB, Sutherland BJG, et al. The presence of the Pacific oyster OsHV-1 resistance marker on chromosome 8 does not impact susceptibility to infection by *Vibrio aestuarianus*. bioRxiv 2024.11.25.625178; doi: https://doi.org/10.1101/2024.11.25.625178      


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
Use the single-SNP per locus data in VCF format.      

#### l. Population genetic analysis ####
Read the data into R using the script `01_scripts/01_sps_char_and_filt.R`.     
This will use the repository `simple_pop_stats` functions to read in the VCF file, prepare populations colour file, calculate and plot per-sample missing data, calculate allele frequencies after removing high missing data individuals, and save output.      


Analyze with `01_scripts/02_sps_analysis.R`     
This will run a principal components analysis (PCA).     


#### m. Genome-wide association study (GWAS) ####
Using the saved output of `01_scripts/02_sps_analysis.R` above, the script `01_scripts/03_GWAS.R` will use the loaded VCF file, filter to only the retained individuals, bring in the phenotype file, and create GEMMA input files.     
- `sample_interp_2024-07-18.csv`   
- non-imputed VCF `populations.snps.vcf`   

The script will do the following:   
i. Change linkage group (LG) RefSeq genome annotations to corresponding chromosome annotations and removes all other contigs in the RefSeq genome    
ii. Optional: impute missing genotypes within each family independently using mean imputation    
iii. Prepare GEMMA inputs      
iv. Plot GEMMA results in Manhattan plots.    

This Rscript will export the following files into `03_results`:    
- gwasanno.txt (annotation of the SNPs)
- gwasgeno.txt (genotype matrix)
- gwaspheno.txt (binary dead/alive vector)
- gwaspheno2.txt (numeric days to death vector)

Midway through the R script, use GEMMA to run a GWAS using these files:     
```
# Run in command-line
# change directory into 03_results

# Calculate kinship matrix
gemma -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_allfam

# Run the GWAS analysis with the family covariate and binary phenotype
gemma -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam.cXX.txt -n 1 -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam

# Go back to R script 01_scripts/04_GWAS.R to generate a Manhattan plot

# Then save the output folder and Manhattan plot into a subfolder, e.g., 03_results/gemma_dead_or_alive_no_impute/ 

# note: to analyze the second phenotype, simply replace gwaspheno.txt with gwaspheno2.txt, but do this after Manhattan plot was generated.    
```

Note: SNPs are coded 0, 1, 2 to correspond to reference homozygote, heterozygote, and alternative homozygote genotypes, respectively. Therefore, the models assume the allele are additive (rather than exhibiting dominance, for example). 


### 02. OCV23 rhAmp analysis ###
Input files are csv files that are raw output from the genotyping platform (i.e., CFX96 instrument).     
Copy all csv files with rhAmp results into `02_input_data`. The following column names are required, and should be standard in qPCR machine output:    
`Well`, `Fluor`, `Content`, `Cq`    

Copy the file `OSU_CHR8_VC_Mapping_Family114-117_2024-01-21_MFrhAmpDNAid.txt` (available from [FigShare](https://doi.org/10.6084/m9.figshare.26515438.v1)) into `00_archive`. This will be used to connect the well and plate ID to the sample ID.    

Run the Rscript `01_scripts/rhamp_assay_analysis.R` interactively to do the following:    
- format as needed.    
- #todo: update this section

#### Plotting mortality by genotype ####
Put `mort_prop_by_fam.csv` and `sample_day_of_death_and_DNA_ID.csv` in `00_archive`.     
(#TODO: need to include these in FigShare repo)

After running the above script, you should have a file called `03_results/rhAmp_per_sample_genotype_summary.csv` that will be used by the plotting script.     

Use the following script interactively to plot the results:     
`01_scripts/rhamp_06_plot_results.R`     

