# ms_cgig_chr8
Code repository to accompany all CHR8 analyses, which currently includes the following:     
- [1. wgrs of OsHV-1 exposure survivors (CHR8 families)](#01-wgrs-of-oshv-1-exposure)    
- [2. ddRADseq of Vibrio challenge survivors (CHR8 families)](#02-ocv23-analysis)      
- [3. rhAmp assay for CHR8 genotypes in Vibrio challenge survivors (CHR8 families)](#03-OCV23-rhAmp-analysis)     
- [4. OA exposure of families from crosses at VIU (VIU families)](#04-OA-exposed-families)      
- [5. Amplicon panel work for OsHV-1 trial](#05-Amplicon-based-association-of-OsHV-1-trial)       


#### Requirements ####
[amplitools](https://github.com/bensutherland/amplitools)      
[simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)    
[stacks_workflow](https://github.com/enormandeau/stacks_workflow)       
[wgrs_workflow](https://github.com/bensutherland/wgrs_workflow)        
[GEMMA](https://github.com/genetics-statistics/GEMMA/tree/master)         
[amplitargets](https://github.com/bensutherland/amplitargets)         

### 01. WGRS of OsHV-1 exposure ###
Follow `wgrs_workflow`.       


### 02. OCV23 analysis ###
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


#### l. Population genetic analysis ####
Read the data into R using the script `01_scripts/01_import_plink_to_genind.R`           

#### m. Genome-wide association study (GWAS) ####
The `01_scripts/GWAS.R` script does the following using `populations.snps_single-SNP_per_tag_2023-10-23.vcf` and `populations.plink_2023-10-23.map` as input:    
i. Changes linkage group (LG) RefSeq genome annotations to corresponding chromosome annotations and removes all other contigs in the RefSeq genome    
ii. Imputes missing genotypes within each family independently using mean imputation    
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


### 03. OCV23 rhAmp analysis ###
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




### 04. OA-exposed families ###
Genotype samples using `wgrs_workflow`.     

Put the LD-filtered VCF from `wgrs_workflow` in `ms_cgig_chr8/02_input_data`, then change into the `ms_cgig_chr8` main directory.       

Note: will need to convert the bcf to a vcf file, via:    
`bcftools view 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.001_AF_0.05_LD0.5w50kb.bcf  -Ov -o 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.001_AF_0.05_LD0.5w50kb.vcf`       

Note: you may be using a subset version of a VCF file here to test out the pipeline (e.g., 5% of total lines).     

Use the following interactive Rscripts (depend on `simple_pop_stats`):      
`01_scripts/COARL_01.R`      
...this script will:      
- allow user to set variables (VCF file, parent ID file, cross and phenotype file)
- load the genotype data and update identifiers 
- prepare a PCA and calculate relatedness of individuals
- plot allele frequencies of all loci
Output will be plots of the above, and `03_results/prepared_data.RData`.       

Second, use:       
`01_scripts/COARL_02.R`     
...this script will:     
- bring in cross information (i.e., which dam and sire were used for each family)
- infer expected offspring allele frequencies based on parental genotypes
Output will be `03_results/per_family_inferred_allele_frequency_data.RData`      

Third, use:     
`01_scripts/COARL_03.R`      
...this script will:     
- allow user to set variables (which phenotype to focus on)      
- prepare genotype, phenotype, and marker files in preparation for gemma (BIMBAM format)
Output will be gemma inputs, saved into `simple_pop_stats/03_results/<your_date-stamped_subdirectory>`.      

After the above scripts are run, change into the new directory that was created within `simple_pop_stats`, as given as the example below:     
```
cd ../simple_pop_stats/03_results/gemma_run_dw_size_mean_2024-05-14_14h53 

# compute Kinship matrix
gemma -g ./gemma_geno.txt -p gemma_pheno_dw_size_mean.txt -gk -maf 0.05 -o dw_size
# note: this will output into an output folder

# association test
gemma -g gemma_geno.txt -p gemma_pheno_dw_size_mean.txt -k output/dw_size.cXX.txt -n 1 -a gemma_geno_annot.txt -maf 0.05 -lmm 4 -o gwas_dw_size 

```

Next, use:     
`01_scripts/COARL_04.R`      
...this script will:     
- allow user to choose the gemma result folder and file
- load gemma results and view distribution of p-values
- generate Manhattan plot based on gemma output


### 05. Amplicon-based association of OsHV-1 trial ###
Requires the following inputs, put in `02_input_data`:      
- plink ped and map files from offspring genotyped by amplicon panel        
- VCF files from parents genotyped by amplicon panel (best replicate)
- genotypes from `wgrs_workflow` run of 20X parents 
- panel contig and SNP position info: `additional_file_S1_amp_panel_design_info.txt` (put in `00_archive`)     
- survival phenotype data from OsHV-1 trial: `qcat992_sample_mort_pheno_2024-06-17.txt` (put in `00_archive`)     

note: the .ped and .map files are used instead of the supplied VCF files for offspring due to a file format issue specific to the offspring VCF.     


#### 05.a. Prepare input data #### 
##### Offspring data #####
An issue arose where the .ped and .map files do not indicate the correct reference allele. To solve this, put the two files into `02_input_data`, clone amplitargets at the same level as this repo, and use the following script to prepare a ref allele file to inform the plink conversion:      
`01_scripts/correct_orientation_of_plink_map.R`     
...this will produce `02_input_data/G0923-21-VIUN_set_allele.txt`    

Change directory into `02_input_data` and run the following to correct the allele to the REF/ALT from the hotspot file:    
`plink2 --ped 02_input_data/G0923-21-VIUN.ped --map 02_input_data/G0923-21-VIUN.map --ref-allele 02_input_data/G0923-21-VIUN_set_allele.txt --recode vcf --out 02_input_data/G0923-21-VIUN_corr_alleles`     

It is also necessary to add the chromosome and positional information to the VCF file:     
To do this, download Additional File S1 from Sutherland et al. 2024, which provides the coordinates of each marker, save it as a .txt file in the present repo, `00_archive`. Then use the following script to update the contig and positional info in the provided VCF:    
`01_scripts/chr8_oshv1_trial_amp_01_prep_vcf.R`     

Output: `02_input_data/G0923-21-VIUN_corr_alleles_annot.vcf.gz`    
Decompress: `gunzip 02_input_data/G0923-21-VIUN_corr_alleles_annot.vcf.gz`


##### Parent data #####
Put the parent data in `02_input_data`. This data came with the correct ref allele and positional info.         
Use the following to merge parent data into a single file:    
```
# Decompress the compressed VCF files:    
gunzip 02_input_data/*.gz

# Compress the parent files using bgzip
ls 02_input_data/TSVC_variants_IonCode_0*.vcf | xargs -n 1 bgzip

# Index the parent files using bcftools
ls 02_input_data/TSVC_variants_IonCode_0*.vcf.gz | xargs -n 1 bcftools index

# Create filelist for merging the parent files
ls -1 02_input_data/TSVC_variants_IonCode_0*.vcf.gz > 02_input_data/parent_VCF_filelist.txt

# merge VCF files
bcftools merge --file-list ./02_input_data/parent_VCF_filelist.txt -Ov -o ./02_input_data/amp_panel_all_parents.vcf

# note: these have novel variants also, which will need to be removed eventually (below)    

```

#### 05.b. Convert to chromosome assembly coordinates ####
Clone a snplift repo for each of offspring and parent data.    
```
cd ..
git clone https://github.com/enormandeau/snplift.git snplift_offspring    
git clone https://github.com/enormandeau/snplift.git snplift_parents

cp ./ms_cgig_chr8/02_input_data/amp_panel_all_parents.vcf ./snplift_parents/04_input_vcf/     
cp ./ms_cgig_chr8/02_input_data/G0923-21-VIUN_corr_alleles_annot.vcf ./snplift_offspring/04_input_vcf/

```

Within each snplift repository, edit the `02_infos/snplift_config.sh` file with the following edits:    
- full path to the original genome (bwa indexed)
- full path to the target genome (bwa indexed)
- CORRECT_ALLELES=1 to convert the ref all allele when alignments are reverse complemented
- original VCF filename
- target (new) VCF filename
Note: setting the `CORRECT_ID` to 0 above prevents the ID column from being recalculated, so that your original IDs are carried through to the new VCF.       

Run snplift for each:    
```
cd snplift_offspring
time ./snplift 02_infos/snplift_config.sh      
cp ./G0923-21-VIUN_corr_alleles_annot_roslin.vcf ../ms_cgig_chr8/03_results/

cd ..

cd snplift_parents
time ./snplift 02_infos/snplift_config.sh      
cp ./amp_panel_all_parents_roslin.vcf ../ms_cgig_chr8/03_results/

```


#### 05.c. Combine offspring and parent data #####
```
# Change directory 
cd 03_results

# Fix header after snplift for both files
bcftools reheader G0923-21-VIUN_corr_alleles_annot_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output G0923-21-VIUN_corr_alleles_annot_roslin_rehead.vcf

bcftools reheader amp_panel_all_parents_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output amp_panel_all_parents_roslin_rehead.vcf

# Compress the snplift'ed, rehead'ed VCF files with bgzip to prepare for bcftools conversion
ls *_roslin_rehead.vcf | xargs -n 1 bgzip

# Convert each to .bcf file
bcftools view G0923-21-VIUN_corr_alleles_annot_roslin_rehead.vcf.gz -Ob -o G0923-21-VIUN_corr_alleles_annot_roslin_rehead.bcf
bcftools view amp_panel_all_parents_roslin_rehead.vcf.gz -Ob -o amp_panel_all_parents_roslin_rehead.bcf

# Index each .bcf file
bcftools index amp_panel_all_parents_roslin_rehead.bcf
bcftools index G0923-21-VIUN_corr_alleles_annot_roslin_rehead.bcf

# Run isec to compare between the files (note, see folder structure)
mkdir isec_output
bcftools isec ./G0923-21-VIUN_roslin_rehead.bcf ./amp_panel_all_parents_roslin_rehead.bcf -p isec_output/

```

Logically, it only makes sense to merge the parent and offspring data that are common between the two, which would mean files 0002.vcf (offspring) and 0003.vcf (parents) common to both.    
```
# Compress the VCF files
ls *.vcf | xargs -n 1 bgzip

# Index the VCF files
ls *.vcf.gz | xargs -n 1 bcftools index

```
(#TODO: add details about comparing genotypes between the technologies)      
 

##### Next phase #####
Now that the genotypes by panel have been confirmed to be finding the same results as the wgrs data, the next step will be as follows:      

a) exclude loci from the wgrs parent file that are present in the panel datafile:       
Assumes you have put the 'source' files into the folders as specified below.    
```
# copy the filtered parent wgrs file into the repo, and index
cp -l ../00_source_materials/parent_wgrs_genotypes/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./03_results/
bcftools index 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf

# the SNPlift'd parent panel file should already be in the repo
# 03_results/amp_panel_all_parents_roslin_rehead.bcf

# Create a folder for isec output
mkdir 03_results/isec_output_wgrs_filt_parents_and_panel_offspr

# Run isec with flag to collapse all loci regardless of matching alleles
bcftools isec --collapse all ./03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./03_results/amp_panel_all_parents_roslin_rehead.bcf -p 03_results/isec_output_wgrs_filt_parents_and_panel_offspr 

## Interpretation:    
# 0000.vcf = private to wgrs
# 0001.vcf = private to panel
# 0002.vcf = records from wgrs shared in both
# 0003.vcf = records from panel shared in both

# Therefore the desired file is 0000.vcf
cp -l 03_results/isec_output_wgrs_filt_parents_and_panel_offspr/0000.vcf ./03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci.vcf

```

b) Now that the wgrs has all overlapping loci with panel file removed, prepare the wgrs file to be combined:       
```
# Prepare to rename wgrs samples so they match between the files
bcftools query -l 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci.vcf > 03_results/wgrs_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 03_results/wgrs_parents_orgn_names.txt -o ./03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed.vcf ./03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci.vcf

# Collect the new wgrs parent names and sort the names into a file
bcftools query -l 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed.vcf | sort > ./03_results/sorted_wgrs_samples.txt

# Sort in the VCF file
bcftools view -S 03_results/sorted_wgrs_samples.txt ./03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed.vcf -o 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed_sorted_samples.vcf

# Compress and index 
bgzip 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed_sorted_samples.vcf
bcftools index 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed_sorted_samples.vcf.gz

```

c) Prepare the panel file to be combined:       
```
# Prepare to rename panel samples so they match between the files
bcftools query -l 03_results/amp_panel_all_parents_roslin_rehead.bcf > 03_results/amp_panel_parents_original_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 03_results/amp_panel_parents_original_names.txt 03_results/amp_panel_all_parents_roslin_rehead.bcf -o 03_results/amp_panel_all_parents_roslin_rehead_renamed.bcf

# Collect the new amp panel parent names and sort the names into a file
bcftools query -l amp_panel_all_parents_roslin_rehead_rename.bcf | sort > sorted_amp_parents_samples.txt

# Sort in the BCF file
bcftools view -S 03_results/sorted_amp_parents_samples.txt ./03_results/amp_panel_all_parents_roslin_rehead_renamed.bcf -o 03_results/amp_panel_all_parents_roslin_rehead_renamed_sorted_samples.vcf

# Compress and index
bgzip 03_results/amp_panel_all_parents_roslin_rehead_renamed_sorted_samples.vcf
bcftools index 03_results/amp_panel_all_parents_roslin_rehead_renamed_sorted_samples.vcf.gz

```

d) Combine the parent data with bcftools concat
```
bcftools concat 03_results/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_no_panel_loci_renamed_sorted_samples.vcf.gz 03_results/amp_panel_all_parents_roslin_rehead_renamed_sorted_samples.vcf.gz -Ob -o 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf

```


e) Next, want to combine the amp panel offspring into the wgrs+panel parent data:     
```
# Create output for isec
mkdir 03_results/isec_output_wgrs_panel_parents_and_panel_offspr

# Index the merged wgrs and panel parent datafile
bcftools index 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf

# Compare the amp panel offspring and amp panel + wgrs parent datafiles
bcftools isec ./03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead.bcf -p 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/

# Combine the amp panel offspring loci that are present in the amp panel + wgrs parent datafile with the amp panel + wgrs parent datafile
bgzip 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/0003.vcf
bcftools index 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/0003.vcf.gz
bcftools merge 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/0003.vcf.gz 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf -Ob -o 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci_and_offspr_shared_panel_loci.bcf

bcftools merge 
```





#### 05.c. Filter VCF and prepare for gemma analysis #### 
Use script `01_scripts/chr8_oshv1_amp_02_vcf_to_gemma.R`       


#### 05.d. Gemma analysis ####
Put output of the above script into `03_results`, then change directly into this folder to run the following commands.     
```
cd 03_results
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar
```

Then go to `chr8_oshv1_amp_03_gemma_results.R`.    


#### 05.e. Compare wgrs to amp panel output ####
Use the above instructions to create a corrected amp panel, snplifted VCF file, then run:      
```
# Compress the VCF then index with tabix (will get error if try to read the VCF file directly)      
bgzip G0923-21-VIUN_annot_snplift_to_roslin_corr.vcf && tabix -p vcf G0923-21-VIUN_annot_snplift_to_roslin_corr.vcf.gz      

# Convert to BCF file
bcftools view G0923-21-VIUN_annot_snplift_to_roslin_corr.vcf.gz -o G0923-21-VIUN_annot_snplift_to_roslin_corr.bcf

# Index
bcftools index G0923-21-VIUN_annot_snplift_to_roslin_corr.bcf
bcftools index mpileup_calls.bcf    

# Run isec to compare between the files (note, see folder structure)
bcftools isec ./00_source_materials/amp_panel_genotypes/G0923-21-VIUN_annot_snplift_to_roslin_corr.bcf ./00_source_materials/parent_wgrs_genotypes/mpileup_calls.bcf -p compare_amp_panel_and_wgrs_parents_all_loci/


```


