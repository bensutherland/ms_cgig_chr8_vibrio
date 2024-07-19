#!/bin/bash
# Convert BCF file to ai2 genotype file 

# Set user variables
INPUT_FOLDER="12_impute_impute"
INPUT_BCF="all_inds_wgrs_and_panel_no_multiallelic.bcf"

# Obtain individual names
bcftools query -l $INPUT_FOLDER/$INPUT_BCF > $INPUT_FOLDER/temp_ind_names.txt

echo "mname" | cat - $INPUT_FOLDER/temp_ind_names.txt |
        tr "\n" "\t" | 
        sed 's/\t$//g' > $INPUT_FOLDER/header.txt

# Add newline character
echo "" >> $INPUT_FOLDER/header.txt

# BCFtools extract chr, position, and genotypes, then convert to ai2 format
# bcftools query -f '%CHROM %POS[\t%GT]\n' $INPUT_FOLDER/$INPUT_BCF | 
#         sed 's/\b0\/0\b/0/g' | 
#         sed 's/\b0\/1\b/1/g' | 
#         sed 's/\b1\/1\b/2/g' | 
#         sed 's/\.\/\./9/g' > $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2.txt


cat $INPUT_FOLDER/header.txt $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2.txt > $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2_full.txt

# THEN REMOVE

