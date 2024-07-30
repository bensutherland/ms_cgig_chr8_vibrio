# ms_cgig_chr8 amplicon-based imputation and GWAS for OsHV-1 survivorship

### 00. Getting started ###
Clone the present repo, all commands will occur in the repo unless indicated.    

Requires the following inputs put in `10_impute_input`, each in its own subfolder as labeled:      
- parent wgrs (20X) filtered genotypes BCF file from `wgrs_workflow` in `parent_wgrs`   
- parent amplicon panel individual vcf.gz files in `parent_panel`   
- offspring amplicon panel individual vcf.gz files `offspring_panel`   
- survival phenotype data from OsHV-1 trial: `qcat992_sample_mort_pheno_2024-06-17.txt` (put in `00_archive`)     

Note: only include a single replicate per individual by the panel (pick the best replicate)   
Note: each type of amp panel data needed its own subfolder because they are named by the barcode and therefore will overlap with other datasets.     
Note: do not copy links of the VCF files, but rather full files.    

### 01. Prepare input data ### 
##### Offspring panel data #####
Merge the data into a single multi-sample VCF file:      
```
# decompress the files, then compress with bgzip
gunzip 10_impute_input/offspring_panel/*.gz
ls 10_impute_input/offspring_panel/*.vcf | xargs -n 1 bgzip

# index the files with bcftools
ls 10_impute_input/offspring_panel/*.vcf.gz | xargs -n 1 bcftools index

# create filelist for merging all VCF files
ls -1 10_impute_input/offspring_panel/*.vcf.gz > 10_impute_input/offspring_panel/sample_list.txt

# merge all VCF files
bcftools merge --file-list 10_impute_input/offspring_panel/sample_list.txt -Ov -o 10_impute_input/offspring_panel.vcf

```

##### Parent panel data #####
Merge the data into a single multi-sample VCF file:      
```
# decompress the files, then compress with bgzip
gunzip 10_impute_input/parent_panel/*.gz
ls 10_impute_input/parent_panel/*.vcf | xargs -n 1 bgzip

# index the files with bcftools
ls 10_impute_input/parent_panel/*.vcf.gz | xargs -n 1 bcftools index

# create filelist for merging all VCF files
ls -1 10_impute_input/parent_panel/*.vcf.gz > 10_impute_input/parent_panel/sample_list.txt

# merge all VCF files
bcftools merge --file-list 10_impute_input/parent_panel/sample_list.txt -Ov -o 10_impute_input/parent_panel.vcf

```
note: these have novel variants also, which will need to be removed eventually (below)    

##### Parent wgrs data #####
```
bcftools index 10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf
```

### 02. Convert to chromosome assembly coordinates ###
Move above the repo, then clone a snplift repo for each of offspring and parent data:      
```
cd ..
git clone https://github.com/enormandeau/snplift.git snplift_offspring    
git clone https://github.com/enormandeau/snplift.git snplift_parents

cp ./ms_cgig_chr8/10_impute_input/offspring_panel.vcf ./snplift_offspring/04_input_vcf/
cp ./ms_cgig_chr8/10_impute_input/parent_panel.vcf ./snplift_parents/04_input_vcf/
```

Within each snplift repository, edit `02_infos/snplift_config.sh`:    
```
- full path to the original genome (bwa indexed)
- full path to the target genome (bwa indexed)
- relative path to the original VCF filename
- relative path to the new VCF filename
- CORRECT_ALLELES=1 to convert the ref all allele when alignments are reverse complemented
Note: setting the `CORRECT_ID` to 0 above prevents the ID column from being recalculated, so that your original IDs are carried through to the new VCF.       
```

Run snplift for each:    
```
cd snplift_offspring
time ./snplift 02_infos/snplift_config.sh      
cp ./offspring_panel_roslin.vcf ../ms_cgig_chr8/10_impute_input/

cd ..

cd snplift_parents
time ./snplift 02_infos/snplift_config.sh      
cp ./parent_panel_roslin.vcf ../ms_cgig_chr8/10_impute_input/

```

Return to `ms_cgig_chr8` repo, then further prepare panel data post-SNPlift:    
```
# add headers
bcftools reheader 10_impute_input/offspring_panel_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output 10_impute_input/offspring_panel_roslin_rehead.vcf

bcftools reheader 10_impute_input/parent_panel_roslin.vcf --fai ~/genomes/GCF_902806645.1_cgigas_uk_roslin_v1_genomic.fna.fai --output 10_impute_input/parent_panel_roslin_rehead.vcf

# compress
ls 10_impute_input/*_rehead.vcf | xargs -n 1 bgzip

# convert to BCF file
bcftools view 10_impute_input/offspring_panel_roslin_rehead.vcf.gz -Ob -o 10_impute_input/offspring_panel_roslin_rehead.bcf

bcftools view 10_impute_input/parent_panel_roslin_rehead.vcf.gz -Ob -o 10_impute_input/parent_panel_roslin_rehead.bcf

# index
ls 10_impute_input/*_rehead.bcf | xargs -n 1 bcftools index

```

Remove novel variants from both panel files    
```
# Remove novel variants from offspring file 
# Identify hotspots (i.e., field 3 of VCF file has a marker name, and is not '.')
bcftools view 10_impute_input/offspring_panel_roslin_rehead.bcf | grep -vE '^#' - | awk '$3 != "." { print $1 "\t" $2 }' - > 10_impute_input/offspring_include_snps.txt

# Use bcftools to only keep these loci from VCF file
bcftools view --targets-file 10_impute_input/offspring_include_snps.txt ./10_impute_input/offspring_panel_roslin_rehead.bcf -Ob -o 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf

# Index
bcftools index 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf

## As above, but remove from parents file ##
bcftools view 10_impute_input/parent_panel_roslin_rehead.bcf | grep -vE '^#' - | awk '$3 != "." { print $1 "\t" $2 }' - > 10_impute_input/parent_include_snps.txt

bcftools view --targets-file 10_impute_input/parent_include_snps.txt ./10_impute_input/parent_panel_roslin_rehead.bcf -Ob -o 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf

bcftools index 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf
```


### 03. Exclude panel loci from parent wgrs data ###
Before we merge the panel and wgrs loci from the parents, we need to remove all panel loci (hotspot + novel) from the wgrs datafile to not have conflicting loci present.    

```
# prepare an output folder for bcftools isec
mkdir 11_impute_combine/isec_rem_panel_from_wgrs/

# run isec to identify loci private to wgrs data, incl. --collapse all flag to collapse regardless of alleles.    
bcftools isec --collapse all ./10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf -p 11_impute_combine/isec_rem_panel_from_wgrs/

## Interpretation:    
# 0000.vcf = private to parent wgrs
# 0001.vcf = private to parent panel
# 0002.vcf = records from parent wgrs shared in both
# 0003.vcf = records from parent panel shared in both

# Save the private records to parent wgrs
cp -l 11_impute_combine/isec_rem_panel_from_wgrs/0000.vcf 11_impute_combine/parent_wgrs_only.vcf

# TODO: add here how to see how many matching alleles there are
```


### 04. Concatenate parent panel loci into parent wgrs only data ###
Prepare the parent wgrs-only file to be combined:       
```
# Prepare to rename wgrs samples so they match between the files
bcftools query -l 11_impute_combine/parent_wgrs_only.vcf > 11_impute_combine/wgrs_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 11_impute_combine/wgrs_parents_orgn_names.txt -o 11_impute_combine/parent_wgrs_only_renamed.vcf 11_impute_combine/parent_wgrs_only.vcf

# Create sorted text file of new names
bcftools query -l 11_impute_combine/parent_wgrs_only_renamed.vcf | sort > 11_impute_combine/wgrs_parents_new_names_sorted.txt

# Sort in the VCF file
bcftools view -S 11_impute_combine/wgrs_parents_new_names_sorted.txt 11_impute_combine/parent_wgrs_only_renamed.vcf -o 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf

# Compress and index 
bgzip 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf

bcftools index 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz
```

Prepare the parent panel data to be combined:       
```
# Copy the parent panel data into the combined folder
cp -l 10_impute_input/parent_panel_roslin_rehead_hotspot_only.bcf 11_impute_combine/

# Prepare to rename panel samples so they match between the files
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead_hotspot_only.bcf > 11_impute_combine/panel_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 11_impute_combine/panel_parents_orgn_names.txt 11_impute_combine/parent_panel_roslin_rehead_hotspot_only.bcf -o 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf

# Collect the new amp panel parent names and sort the names into a file
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf | sort > 11_impute_combine/panel_parents_new_names_sorted.txt

# Sort in the BCF file
bcftools view -S 11_impute_combine/panel_parents_new_names_sorted.txt 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed.bcf -o 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf

# index
bcftools index 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf

```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz 11_impute_combine/parent_panel_roslin_rehead_hotspot_only_renamed_sorted.bcf -Ob -o 11_impute_combine/parent_wgrs_and_panel.bcf

# Index
bcftools index 11_impute_combine/parent_wgrs_and_panel.bcf

```

### 05. Merge parent wgrs and panel data with offspring panel data ###
a) Identify loci in the offspring panel file that overlap with the wgrs+panel parent file     
```
# Bring offspring data to folder
cp -l 10_impute_input/offspring_panel_roslin_rehead_hotspot_only.bcf* ./11_impute_combine/

# Create isec folder to capture output
mkdir 11_impute_combine/isec_combine_parents_and_offspring/

# Use isec to compare the files (no need for --collapse here, want only common shared REF alleles)
bcftools isec 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only.bcf -p 11_impute_combine/isec_combine_parents_and_offspring/

## Interpretation:    
# 0000.vcf = private to parents (wgrs+panel)
# 0001.vcf = private to offspring (panel)
# 0002.vcf = records from parents (wgrs+panel) shared in both
# 0003.vcf = records from offspring (panel) shared in both

# Save and rename 0003.vcf
cp 11_impute_combine/isec_combine_parents_and_offspring/0003.vcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf

# Compress and index
bgzip 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf
bcftools index 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz
```

b) Combine the wgrs+panel parent data with the panel offspring data
```
bcftools merge 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_hotspot_only_common_w_parents.vcf.gz -Ob -o 11_impute_combine/all_inds_wgrs_and_panel.bcf

bcftools index 11_impute_combine/all_inds_wgrs_and_panel.bcf

# Copy the all-data file into the imputation folder, and index
cp -l 11_impute_combine/all_inds_wgrs_and_panel.bcf* 12_impute_impute/

# Remove multiallelic sites
bcftools view --max-alleles 2 ./12_impute_impute/all_inds_wgrs_and_panel.bcf -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf    

```

### 06. Imputation ###
The data is now all in a single BCF file and is ready for the imputation process.     

Prepare a pedigree file:      
```
bcftools query -l 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic.bcf > 12_impute_impute/pedigree.txt 

# Annotate the above file as follows:    
# <indiv> <sire> <dam>     
# where if there is no sire or dam, put 0
# save as space-delimited with the suffix `_annot.txt`.    
```

Format from BCF to AlphaImpute2:       
```
# Prepare ai2 matrix by building a header, and extracting info from the BCF, and converting to ai2 format
./01_scripts/bcf_to_ai2.sh
# produces: 12_impute_impute/all_inds_wgrs_and_panel_no_multiallelic_ai2.txt

# the output will be in the following format:    
# mname \t ind1 \t ind2 \t (...)
# NC_047559.1 2945 \t 0 \t 0 \t 1 \t 9 (...)
# (...)
# where 0, 1, 2 are the number of alt alleles, and 9 is missing data  

# Split ai2 matrix into individual chr. This requires setting the input filename, and an identifiable chromosome string.
01_scripts/prep_geno_matrix_for_ai2.R   
# output will be in 12_impute_impute/ai2_input_<NC_047559.1>.txt, one file per chr 

```

Run imputation:     
```
# initialize the conda environment
conda activate ai2

# Run AlphaImpute2 on chromosome-separated datafiles
01_scripts/run_ai2.sh
# produces: 12_impute_impute/ai2_input_<NC_047559.1>.genotypes and *.haplotypes

# Transpose chromosome-separated imputed .genotypes files, and drop marker names on all matrices but the first to prepare for recombining the files back together
01_scripts/impute_rebuild_chr_lightweight.R
# produces: 13_impute_compare/*.genotypes_transposed_to_combine.txt

# Combine imputed, transposed, chr-separated files back together, then add marker names back in, based on the input ai2 file (before chromosome separation) 
01_scripts/combine_transposed_ai2_output_and_mnames.sh
# produces: 13_impute_compare/all_chr_combined.txt

```

### 07. Evaluate imputation ###
Evaluate results by comparing the imputed data with the 10X 'empirical' data:     
```
# Obtain 10X bcf file 
cp -l ~/Documents/cgig/CHR8_wgrs/wgrs_workflow_offspring/05_genotyping/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./13_impute_compare/
cp ~/Documents/cgig/CHR8_wgrs/wgrs_workflow/05_genotyping/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./13_impute_compare/ 

# Use bash script to pull out genotypes into text file in ai2 format
# Edit the following script to point to the above bcf file, and run
01_scripts/bcf_to_ai2.sh

# Compare the imputed and empirical ai2 file genotypes by chromosome using: 
01_scripts/eval_impute_lightweight.R
# produces plots of average concordance between methods per individual by chromosome, and other outputs to screen
```


### 08. Run GWAS on imputed data and plot results ###
Prepare GEMMA inputs:     
`01_scripts/imputed_ai2_to_gemma.R`      

Run GEMMA:    
```
cd 12_impute_impute
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar

```


Plot GEMMA outputs:    
`01_scripts/imputed_plot_gemma_results.R`    



#### 06. Filter VCF and prepare for gemma analysis #### 
Use script `01_scripts/chr8_oshv1_amp_02_vcf_to_gemma.R`       


#### 07. Gemma analysis ####
Put output of the above script into `03_results`, then change directly into this folder to run the following commands.     
```
cd 03_results
gemma -g gwas_geno.txt -p gwas_pheno.txt -gk -maf 0.05 -o gwas_all_fam
gemma -g gwas_geno.txt -p gwas_pheno.txt -k output/gwas_all_fam.cXX.txt -n 1 -c gwas_covar.txt  -maf 0.05 -lmm 4 -o gwas_all_fam_covar
```

Then go to `chr8_oshv1_amp_03_gemma_results.R`.    

Optional: create a reduced dataset for testing:     
```
# Subset only a single chr for testing
bcftools view 11_impute_combine/all_inds_wgrs_and_panel.bcf --regions NC_047567.1 -Ob -o 12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.bcf

# Show number of panel loci
bcftools view 12_impute_impute/all_inds_wgrs_and_panel_NC_047567_1.bcf | grep -vE '^#' - | awk '{ print $3 }' - | sort | uniq -c | sort -nk1 | less

```        
