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
##### Offspring data #####
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

##### Parent data #####
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


#### 02. Convert to chromosome assembly coordinates ####
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


#### 03. Combine offspring and parent data #####
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
 

#### 04. Combine offspring and parent data #####
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

# Index
bcftools index 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf

```

##### Combine panel offspring file with wgrs+panel parent file #####
a) Identify loci in the offspring panel file that overlap with the wgrs+panel parent file     
```
# Use isec to compare the files (no need for --collapse here, want only common shared REF alleles)
mkdir 03_results/isec_output_wgrs_panel_parents_and_panel_offspr
bcftools isec ./03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead.bcf -p 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/

## Interpretation:    
# 0000.vcf = private to parents (wgrs+panel)
# 0001.vcf = private to offspring (panel)
# 0002.vcf = records from parents (wgrs+panel) shared in both
# 0003.vcf = records from offspring (panel) shared in both

# Save and rename 0003.vcf
cp 03_results/isec_output_wgrs_panel_parents_and_panel_offspr/0003.vcf 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead_common_loci.vcf

# Compress and index
bgzip 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead_common_loci.vcf
bcftools index 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead_common_loci.vcf.gz
```

b) Combine the wgrs+panel parent data with the panel offspring data
```
bcftools merge 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci.bcf 03_results/G0923-21-VIUN_corr_alleles_annot_roslin_rehead_common_loci.vcf.gz -Ob -o 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci.bcf

```

#### 05. Imputation ####
Create a reduced dataset for testing:     
```
# Copy the file into a new folder, and index
cp -l 03_results/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci.bcf ./04_impute_panel/

# Subset only a single chr for testing
bcftools view 04_impute_panel/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci.bcf --regions NC_047559.1 -Ob -o 04_impute_panel/wgrs_filtered_parent_loci_amp_panel_parent_loci_amp_panel_offspring_loci_NC_047559.1.bcf   
 
```        

Prepare file for AlphaImpute2 using Rscript:    


`AlphaImpute2 -genotypes 04_impute_panel/genos.txt -pedigree 04_impute_panel/pedigree_annot.csv -out ai2test -maxthreads 1`



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


#### 08. Compare wgrs to amp panel output ####
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


