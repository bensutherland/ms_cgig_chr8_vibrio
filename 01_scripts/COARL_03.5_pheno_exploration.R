# Prepare data for GWAS
# B. Sutherland (initialized 2024-02-14)

# Requires that 01_scripts/COARL_02.R was already run
# Clear workspace, launch simple_pop_stats

# Install/ load additional libraries

load(file = "03_results/per_family_inferred_allele_frequency_data.RData")

head(cross_and_pheno.df)

par(mfrow=c(1,1), mgp = c(3,1,0), mar = c(5,4,4,2) + 0.1)

# Consider percent D-larvae total
hist(cross_and_pheno.df$sw_per_d_tot) # 'control'
hist(cross_and_pheno.df$dw_per_d_tot) # 'experimental'

hist(cross_and_pheno.df$dw_per_d_tot - cross_and_pheno.df$sw_per_d_tot)

plot(x = cross_and_pheno.df$sw_per_d_tot, y = cross_and_pheno.df$dw_per_d_tot)

hist(cross_and_pheno.df$dw_size_mean - cross_and_pheno.df$sw_size_mean, breaks = 20)

par(mfrow=c(2,2))
hist(cross_and_pheno.df$dw_size_mean)
hist(cross_and_pheno.df$sw_size_mean)
plot(x = cross_and_pheno.df$dw_size_mean, y = cross_and_pheno.df$sw_size_mean)
hist(cross_and_pheno.df$dw_size_mean - cross_and_pheno.df$sw_size_mean)

hist(cross_and_pheno.df$dw_size_mean - cross_and_pheno.df$sw_size_mean)

