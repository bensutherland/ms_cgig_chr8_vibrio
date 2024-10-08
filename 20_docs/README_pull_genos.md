### Evaluate specific genotypes to understand GEMMA results ###
Determine what samples are present:    
`bcftools query -l 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.vcf  > 03_results/samples.txt`      

Inspect the GEMMA output in a spreadsheet to find a top outlier, such as 
`NC_047562.1__11997020`.        

Then search for this line in the input BCF file:    
`bcftools view 02_input_data/mpileup_calls_noindel5_miss0.1_SNP_q20_avgDP10_biallele_minDP4_maxDP100_miss0.1_offspring_only_rename.vcf | grep -vE '^#' - | awk '$1=="NC_047562.1" { print $0 } ' - | awk '$2 == "11997020" { print $0 }' - > ./03_results/NC_047562.1__11997020_genotypes.txt`

Note: will need to transpose this in a spreadsheet.    

Then, use the sample list to join to the genotypes (assumed correct order), then bring in the phenotype file to observe. 

