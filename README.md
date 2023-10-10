# ms_cgig_chr8
Code repository to accompany all CHR8 analyses. 

#### Requirements ####
amplitools      
simple_pop_stats


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


#### d. Prepare sample metadata ####



