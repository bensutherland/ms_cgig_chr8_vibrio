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


### 03. Exclude panel loci from parent wgrs data ###
Before we merge the panel and wgrs loci from the parents, we need to remove all panel loci (hotspot + novel) from the wgrs datafile to not have conflicting loci present.    

```
# prepare an output folder for bcftools isec
mkdir 11_impute_combine/isec_rem_panel_from_wgrs/

# run isec to identify loci private to wgrs data, incl. --collapse all flag to collapse regardless of alleles.    
bcftools isec --collapse all ./10_impute_input/parent_wgrs/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1.bcf ./10_impute_input/parent_panel_roslin_rehead.bcf -p 11_impute_combine/isec_rem_panel_from_wgrs/

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
cp -l 10_impute_input/parent_panel_roslin_rehead.bcf 11_impute_combine/

# Prepare to rename panel samples so they match between the files
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead.bcf > 11_impute_combine/panel_parents_orgn_names.txt
# ...then manually annotate the file to separate the old_name and new_name with whitespace, single line per sample

# Rename in the VCF file
bcftools reheader --samples 11_impute_combine/panel_parents_orgn_names.txt 11_impute_combine/parent_panel_roslin_rehead.bcf -o 11_impute_combine/parent_panel_roslin_rehead_renamed.bcf

# Collect the new amp panel parent names and sort the names into a file
bcftools query -l 11_impute_combine/parent_panel_roslin_rehead_renamed.bcf | sort > 11_impute_combine/panel_parents_new_names_sorted.txt

# Sort in the BCF file
bcftools view -S 11_impute_combine/panel_parents_new_names_sorted.txt 11_impute_combine/parent_panel_roslin_rehead_renamed.bcf -o 11_impute_combine/parent_panel_roslin_rehead_renamed_sorted.bcf

# index
bcftools index 11_impute_combine/parent_panel_roslin_rehead_renamed_sorted.bcf

```

Combine the parent data with bcftools concat
```
bcftools concat --allow-overlaps 11_impute_combine/parent_wgrs_only_renamed_sorted.vcf.gz 11_impute_combine/parent_panel_roslin_rehead_renamed_sorted.bcf -Ob -o 11_impute_combine/parent_wgrs_and_panel.bcf

# Index
bcftools index 11_impute_combine/parent_wgrs_and_panel.bcf

```

### Merge parent wgrs and panel data with offspring panel data ###
a) Identify loci in the offspring panel file that overlap with the wgrs+panel parent file     
```
# Bring offspring data to folder
cp -l 10_impute_input/offspring_panel_roslin_rehead.bcf* ./11_impute_combine/

# Create isec folder to capture output
mkdir 11_impute_combine/isec_combine_parents_and_offspring/

# Use isec to compare the files (no need for --collapse here, want only common shared REF alleles)
bcftools isec 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead.bcf -p 11_impute_combine/isec_combine_parents_and_offspring/

## Interpretation:    
# 0000.vcf = private to parents (wgrs+panel)
# 0001.vcf = private to offspring (panel)
# 0002.vcf = records from parents (wgrs+panel) shared in both
# 0003.vcf = records from offspring (panel) shared in both

# Save and rename 0003.vcf
cp 11_impute_combine/isec_combine_parents_and_offspring/0003.vcf 11_impute_combine/offspring_panel_roslin_rehead_common_w_parents.vcf

# Compress and index
bgzip 11_impute_combine/offspring_panel_roslin_rehead_common_w_parents.vcf
bcftools index 11_impute_combine/offspring_panel_roslin_rehead_common_w_parents.vcf.gz
```

b) Combine the wgrs+panel parent data with the panel offspring data
```
bcftools merge 11_impute_combine/parent_wgrs_and_panel.bcf 11_impute_combine/offspring_panel_roslin_rehead_common_w_parents.vcf.gz -Ob -o 11_impute_combine/all_inds_wgrs_and_panel.bcf

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


