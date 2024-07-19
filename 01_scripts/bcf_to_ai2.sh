#!/bin/bash
# Convert BCF file to ai2 genotype file 

# Set user variables
INPUT_FOLDER="12_impute_impute"
INPUT_BCF="all_inds_wgrs_and_panel_no_multiallelic.bcf"

# BCFtools extract chr, position, and genotypes, then convert to ai2 format
bcftools query -f '%CHROM %POS[\t%GT]\n' $INPUT_FOLDER/$INPUT_BCF | 
        sed 's/\b0\/0\b/0/g' | 
        sed 's/\b0\/1\b/1/g' | 
        sed 's/\b1\/1\b/2/g' | 
        sed 's/\.\/\./9/g' > $INPUT_FOLDER/${INPUT_BCF%.bcf}_ai2.txt



