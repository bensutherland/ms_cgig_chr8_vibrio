library(vcfR)
library(fastman)
library(missMethods)

vcf = read.vcfR("populations.snps_single-SNP_per_tag_2023-10-23.vcf")

map = read.table("populations.plink_2023-10-23.map")

#add linkage group (LG) info to map file, based on https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902806645.1/
map$LG = NA
map$LG[map$V1=="NC_047559.1"] = 1
map$LG[map$V1=="NC_047560.1"] = 2
map$LG[map$V1=="NC_047561.1"] = 3
map$LG[map$V1=="NC_047562.1"] = 4
map$LG[map$V1=="NC_047563.1"] = 5
map$LG[map$V1=="NC_047564.1"] = 6
map$LG[map$V1=="NC_047565.1"] = 7
map$LG[map$V1=="NC_047566.1"] = 8
map$LG[map$V1=="NC_047567.1"] = 9
map$LG[map$V1=="NC_047568.1"] = 10

#LG to Chr info can be found in Sup File in https://doi.org/10.1093/gigascience/giab020
map$Chr = NA
map$Chr[map$LG==1] = "Chr7"
map$Chr[map$LG==2] = "Chr1"
map$Chr[map$LG==3] = "Chr9"
map$Chr[map$LG==4] = "Chr6"
map$Chr[map$LG==5] = "Chr3"
map$Chr[map$LG==6] = "Chr2"
map$Chr[map$LG==7] = "Chr4"
map$Chr[map$LG==8] = "Chr5"
map$Chr[map$LG==9] = "Chr10"
map$Chr[map$LG==10] = "Chr8"

map$SNPname = paste(map$Chr,map$V4,sep="_")

#extract genotypes from vcf
geno = extract.gt(vcf, element = "GT")

#change SNP names to "Chr[1-10]_[location in bp]"
rownames(geno) = map$SNPname
#remove SNPs on contigs, i.e., not on Chr1-10
geno = geno[-grep("NA_",rownames(geno)),]

#missingness sanity check
plot(rowSums(is.na(geno))/ncol(geno))
plot(colSums(is.na(geno))/nrow(geno))


#change genotypes to numeric values
geno[geno=="0/0"] = 0
geno[geno=="0/1"] = 1
geno[geno=="1/1"] = 2
geno = t(geno)
mode(geno) = "numeric"

#remove parental samples
geno = geno[-grep("F0",rownames(geno)),]

#create family-specific genotype matrices
geno_F114 = geno[grep("F114",rownames(geno)),]
geno_F115 = geno[grep("F115",rownames(geno)),]
geno_F116 = geno[grep("F116",rownames(geno)),]
geno_F117 = geno[grep("F117",rownames(geno)),]
geno_OFR6.10 = geno[grep("OFR6.10",rownames(geno)),]

#run family-specific mean imputation on genotypes
geno_F114_impute = impute_mean(geno_F114)
geno_F115_impute = impute_mean(geno_F115)
geno_F116_impute = impute_mean(geno_F116)
geno_F117_impute = impute_mean(geno_F117)
geno_OFR6.10_impute = impute_mean(geno_OFR6.10)

#reconstruct full genotype matrix
geno = rbind(geno_F114_impute,
             geno_F115_impute,
             geno_F116_impute,
             geno_F117_impute,
             geno_OFR6.10_impute)


#create variable holding family information
var_family = rep(NA,nrow(geno))
var_family[grep("F114",rownames(geno))] = "F114"
var_family[grep("F115",rownames(geno))] = "F115"
var_family[grep("F116",rownames(geno))] = "F116"
var_family[grep("F117",rownames(geno))] = "F117"
var_family[grep("OFR6.10",rownames(geno))] = "OFR6.10"

#create variable holding live/dead information
var_status = rep(NA,nrow(geno))
var_status[grep("Alive",rownames(geno))] = "1"
var_status[grep("Dead",rownames(geno))] = "0"
var_status = as.numeric(var_status)


#GWAS with all families
gwaspheno = var_status
gwascovar = model.matrix(~as.factor(var_family))
gwasgeno = t(geno)
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwascovar,"gwascovar.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)


#Run in command-line
#gemma-0.98.5 available at https://github.com/genetics-statistics/GEMMA
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_allfam
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam.cXX.txt -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_covar
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_allfam.cXX.txt -n 1 -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_allfam_nocovar



#GWAS with only MBP families
gwaspheno = var_status[-which(var_family=="OFR6.10")]
gwascovar = model.matrix(~as.factor(var_family[-which(var_family=="OFR6.10")]))
gwasgeno = t(geno[-which(var_family=="OFR6.10"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwascovar,"gwascovar.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)


#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -gk -maf 0.05 -o gwas_mbpfam
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_mbpfam.cXX.txt -n 1 -c gwascovar.txt -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_mbpfam_covar
#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -k output/gwas_mbpfam.cXX.txt -n 1 -a gwasanno.txt -maf 0.05 -lmm 4 -o gwas_mbpfam_nocovar


#GWAS with only F114
gwaspheno = var_status[which(var_family=="F114")]
gwasgeno = t(geno[which(var_family=="F114"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)

#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F114

#GWAS with only F115
gwaspheno = var_status[which(var_family=="F115")]
gwasgeno = t(geno[which(var_family=="F115"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)

#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F115

#GWAS with only F116
gwaspheno = var_status[which(var_family=="F116")]
gwasgeno = t(geno[which(var_family=="F116"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)

#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F116

#GWAS with only F117
gwaspheno = var_status[which(var_family=="F117")]
gwasgeno = t(geno[which(var_family=="F117"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)

#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_F117

#GWAS with only OFR6.10
gwaspheno = var_status[which(var_family=="OFR6.10")]
gwasgeno = t(geno[which(var_family=="OFR6.10"),])
gwasgeno = cbind(rownames(gwasgeno),"X","Y",gwasgeno)
gwasanno = cbind(rownames(gwasgeno),
                 sapply(strsplit(rownames(gwasgeno),"_"), `[`, 2),
                 sapply(strsplit(sapply(strsplit(rownames(gwasgeno),"_"), `[`, 1),"Chr"), `[`, 2),
                 0)

write.table(gwaspheno,"gwaspheno.txt",row.names = F,col.names = F)
write.table(gwasanno,"gwasanno.txt",row.names = F,col.names = F,quote = F)
write.table(gwasgeno,"gwasgeno.txt",row.names = F,col.names = F,quote = F)

#./gemma-0.98.5 -g gwasgeno.txt -p gwaspheno.txt -n 1 -a gwasanno.txt -maf 0.05 -lm 2 -o gwas_OFR6.10




gemma_gwas2 = read.table("output/gwas_allfam_covar.assoc.txt")
gemma_gwas = read.table("output/gwas_allfam_covar.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="All familes + fixed covariable (n=165)")



gemma_gwas2 = read.table("output/gwas_allfam_nocovar.assoc.txt")
gemma_gwas = read.table("output/gwas_allfam_nocovar.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="All familes + no covariable (n=165)")


gemma_gwas2 = read.table("output/gwas_mbpfam_covar.assoc.txt")
gemma_gwas = read.table("output/gwas_mbpfam_covar.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="MBP familes + fixed covariable (n=139)")



gemma_gwas2 = read.table("output/gwas_mbpfam_nocovar.assoc.txt")
gemma_gwas = read.table("output/gwas_mbpfam_nocovar.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="MBP familes + no covariable (n=139)")



gemma_gwas2 = read.table("output/gwas_F114.assoc.txt")
gemma_gwas = read.table("output/gwas_F114.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="F114 (n=32)")



gemma_gwas2 = read.table("output/gwas_F115.assoc.txt")
gemma_gwas = read.table("output/gwas_F115.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="F115 (n=37)")


gemma_gwas2 = read.table("output/gwas_F116.assoc.txt")
gemma_gwas = read.table("output/gwas_F116.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="F116 (n=36)")



gemma_gwas2 = read.table("output/gwas_F117.assoc.txt")
gemma_gwas = read.table("output/gwas_F117.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="F117 (n=34)")



gemma_gwas2 = read.table("output/gwas_OFR6.10.assoc.txt")
gemma_gwas = read.table("output/gwas_OFR6.10.assoc.txt",skip = 1)
names(gemma_gwas) = gemma_gwas2[1,]

fastman(gemma_gwas,
        chr = "chr",
        bp = "ps",
        p="p_lrt",
        genomewideline = -log10(0.05/nrow(gemma_gwas)),
        suggestiveline = NULL,
        cex=1.5,cex.lab=1.5,cex.axis=1,
        ylim=c(0,6),
        main="OFR6.10 (n=26)")



