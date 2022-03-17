library(ggplot2) # Version 3.3.5
library(ggpubr) # Version 0.4.0    
library(ggridges) # Version 0.5.2
library(scales) # Version 1.1.1 
library(R.matlab) # Version 3.6.2
library(MatchIt) # Version 4.2.0 
library(xtable) # Version 1.8-4
library(dplyr)

##############################################################################################################
#                                            DATA PREPARATION                                                #
##############################################################################################################
#
# ADULT CELL TYPES
#
# 1. Download the Lake et al. (2018) Supplementary Tables 1-13 from: https://doi.org/10.1038/nbt.4038
# 2. In the spreadsheet, navigate to Table S3
# 3. Copy the values in 'Gene' (Column B) and 'Cluster' (Column H) into a .txt file
# 4. Ensure that the columns in your txt file are named 'Gene' and 'Cluster'
# 5. Save the file as 'Lake2018.txt' into Data/gene_decoding/enrichment_analysis
lake <- read.table('Data/gene_decoding/enrichment_analysis/Lake2018.txt', header = T)
lake = lake[-which(lake$Cluster=='Ast_Cer'),] 
lake = lake[-which(lake$Cluster=='OPC_Cer'),]
lake = lake[-which(lake$Cluster=='Purk1'),]
lake = lake[-which(lake$Cluster=='Purk2'),]
lake = lake[-which(lake$Cluster=='Per'),]
lake = lake[-which(lake$Cluster=='Gran'),]
lake$Cluster <- as.factor(lake$Cluster)
#
# PRENATAL CELL TYPES
#
# 1. Download the Polioudakis et al. (2018) Supplementary Table S4 from: https://doi.org/10.1016/j.neuron.2019.06.011
# 2. In the spreadsheet, navigate to the tab 'Cluster enriched genes'
# 3. Copy the values in 'Gene' (Column B) and 'Cluster' (Column C) into a .txt file
# 4. Ensure that the columns in your .txt file are named 'Gene' and 'Cluster'
# 5. Save the file as 'Polioudakis2018.txt' into Data/gene_decoding/enrichment_analysis
polioudakis <- read.table('Data/gene_decoding/enrichment_analysis/Polioudakis2018.txt',header = T)
polioudakis$Cluster <- as.factor(polioudakis$Cluster)
#
# MAJOR DEPRESSION
#
# 1. Download the Li et al (2018) Supplementary Tables S1 to S16 from: https://doi.org/10.1126/science.aat7615
# 2. In the spreadsheet, navigate to Table S13
# 3. Copy the genes in column T ('MDD'-'Gene') into a .txt file
# 4. Ensure that the column in your .txt file is named 'Gene'
# 5. Save the file as 'Li2018.txt' into Data/gene_decoding/enrichment_analysis
li <- read.table('Data/gene_decoding/enrichment_analysis/Li2018.txt',header = T)
#
# SCHIZOPHRENIA
#
# 1. Download the Schizophrenia Working Group et al. (2020) Supplementary Tables from https://doi.org/10.1101/2020.09.12.20192922
# 2. 
PCG = read.table('Data/gene_decoding/enrichment_analysis/PGC2021.txt', header=T)

###
data= 'NSPN' # nspn.gsr, regr_FD_by_sex, glob.fc.regr, motionmatched, CV_stratified
if(data=='NSPN'){
  in.path = out.path = 'Results/gene_decoding/'
  out.path = 'Results/'
  load('Data/connectivity_matrices/nspn.main.RData')
}else{
  in.path = out.path = paste0('Results/replication/',data,'/gene_decoding/')
  out.path = paste0('Results/replication/',data,'/')
  load('Data/connectivity_matrices/',data,'.RData')
}

###
pls1 <- read.csv(paste0(in.path,"plsout/pls1.csv")) # Read in PLS outputs
colnames(pls1)[1]<- 'hgnc_symbol'

# # Load in available gene lengths and chromosome assignments
# genes = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/pls1_biomart.csv') 
# genes = genes[-which(is.na(genes$length)),]

pls1.genes = genes %>% left_join(pls1, by='hgnc_symbol') # Match chromosome assignment and gene length to PLS1
pls1.genes = pls1.genes[order(pls1.genes$z_uncorr)]

genes= pls1.genes[,c(1:7)]
pls1 = pls1.genes[,c(1,7,8,9)]

# create null models
source('Scripts/gene_decoding/gene_nulls_nearest_neighbour.R') # Gene enrichment function
source('Scripts/gene_decoding/run_gene_enrichment.R') # Wrapper script to submit gene lists to gene enrichment function
source('Scripts/gene_decoding/plot_enrichment_ggridges.R') # Plot enrichment function

n_perm = 5000
################ Chromosomal Enrichment ################
dir.create(paste0(out.path, 'tmp/')) # This is where gene length match plots will be stored

chromosome = data.frame(Gene  = genes$hgnc_symbol, Cluster = genes$chromosome_name)
chromosome = chromosome[-which(is.na(match(chromosome$Cluster,levels(chromosome$Cluster)[c(1:22,27:28)]))),]
chromosome$Cluster = droplevels(chromosome$Cluster)

chr_res = run_gene_enrichment(chromosome, pls1, TRUE, genes, n_perm, '~/Desktop') # Run gene enrichment

# Calculate chromosome median rank and SE
chr_median <- sapply(as.character(unique(chromosome$Cluster)),function(x) median(rank(pls1$z_uncorr)[which(genes$chromosome_name==x)]))-median(rank(pls1$z_uncorr))
chr_se <- sapply(as.character(unique(chromosome$Cluster)),function(x) sd(rank(pls1$z_uncorr)[which(genes$chromosome_name==x)])/sqrt(length(which(genes$chromosome_name==x))))
chr_n_genes <- sapply(as.character(unique(chromosome$Cluster)),function(x) length(which(genes$chromosome_name==x)))

chr_plot <- data.frame(chr=as.factor(names(chr_median)),
                       median=chr_median, # medians to plot
                       se_up=chr_median+chr_se,
                       se_down=chr_median-chr_se,
                       sign=ifelse(p.adjust(chr_res$gene_stats$p,method='fdr')<0.05,'*','') # FDR-corrected p-values
)

# Plot chromosomal results as bar plot
chr_plot$chr = factor(chr_plot$chr,levels=chr_plot$chr[order(-chr_plot$median)])
pdf(paste0(out.path,'chromosome.enrichment.pdf'),height=5,width=4)
ggplot(chr_plot,aes(y=median,x=chr))+
  geom_bar(data=chr_plot, aes(y=median,x=chr,fill=(chr_plot$median+1+abs(min(chr_plot$median)))),stat = 'identity') +
  geom_text(aes(label=sign),hjust=ifelse(chr_plot$median<0,2,-2),vjust=0.75,cex=10) +
  scale_fill_gradientn(colours = c("blue","white","red"),limits=c(0,2*abs(min(chr_plot$median))))+
  geom_pointrange(aes(ymin=se_down,ymax=se_up),colour="black")+
  coord_flip()+
  geom_hline(yintercept=0)+
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())+
  xlab('Chromosome') +
  ylab('Median Rank')
dev.off()

# Length stats table for Latex
print(xtable(chr_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)

################ Post Natal Cell Type Enrichment ################

postatal_res = run_gene_enrichment(lake, pls1, TRUE, genes, n_perm, paste0(out.path,'tmp/'))
p.postnatal<-plot_enrichment_ggridges(postatal_res$gene_stats, postatal_res$gene_null, TRUE, -3000, 3000) 
pdf(paste0(out.path, 'postnatal.enrichment.pdf'),4,4); p.postnatal; dev.off()
print(xtable(postnatal_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)

################ Pre Natal Cell Type Enrichment ################

prenatal_res = run_gene_enrichment(polioudakis, pls1, TRUE, genes, n_perm, paste0(out.path,'tmp/'))
p.prenatal<-plot_enrichment_ggridges(prenatal_res$gene_stats, prenatal_res$gene_null, TRUE, -3000, 3000) 
pdf(paste0(out.path,'prenatal.enrichment.pdf'),4,4); p.prenatal; dev.off()
print(xtable(prenatal_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)


################ MDD Enrichment ################
mdd_res = run_gene_enrichment(li, pls1, FALSE, genes, n_perm, paste0(out.path,'tmp/'),'MDD')
p.mdd <- plot_enrichment_ggridges(mdd_res$gene_stats, mdd_res$gene_null, FALSE, -4000, 4000)
pdf(paste0(out.path,'mdd.enrichment.pdf'),height=1.75,width=4); plot(p.mdd); dev.off()
print(xtable(mdd_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)


################ SCZ Enrichment ################
scz_res = run_gene_enrichment(PCG, pls1, FALSE, genes, n_perm, 'Results/gene_decoding/tmp/','SCZ')
p.scz <- plot_enrichment_ggridges(scz_res$gene_stats, scz_res$gene_null, FALSE, -4000, 4000)
pdf('Results/gene_decoding/scz.enrichment.pdf',height=1.75,width=4); plot(p.scz); dev.off()
print(xtable(scz_res$gene_length_stats, booktabs=TRUE), include.rownames=FALSE)

