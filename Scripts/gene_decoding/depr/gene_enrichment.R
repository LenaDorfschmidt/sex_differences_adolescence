#library(biomaRt)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(scales)
library(xtable)
library(R.matlab)
library(MatchIt)
library(viridis)
library(ggpointdensity)

script.path =  paste0('Scripts/')
data.path = paste0('Data/connectivity_matrices/')

source(paste0(script.path, '/helper/brainplot.R'))
load('Data/connectivity_matrices/nspn.main.RData')

sample = 'main' #sample = 'main' #sample = 'CV_stratified' # sample = 'glob.fc.regr' # sample = 'regr_FC_by_sex'

# load in data: NSPN, NSPN low motion or GSR
if (sample == 'main') {
  out_sample = paste0('Results/')
}else{
  out_sample = paste0('Results/replication/',sample,'/')
}

nm = nm.all[17:376]
pls_x = read.csv(paste0(out_sample, '/gene_decoding/plsout/pls_x.csv'), header = F)
nm.plot = nm[pls_x[,2]]
brainplot(c(rep(NA, 16),pls_x[,1]),c(rep('A',16),nm.plot),paste0(out_sample, '/gene_decoding/pls1'),-0.2,0.2,'RdBu','stacked')

if (sample == 'main') {
  diff.mi = read.table(paste0('Results/maturational_index/diff.mi.all.txt'))[17:376,1]
  
  df = data.frame(orig = pls_x$V1, mi = diff.mi[pls_x$V2])
  cor.orig = cor.test(df$orig, df$mi)
  pdf(paste0(out_sample, 'gene_decoding/corr.pls.delta.mi.pdf'),2.5,3)
  ggplot(df, aes(x=orig, y=mi)) + 
    geom_point() +
    geom_smooth(method='lm', color='black') +
    xlab('PLS1 Weights') + 
    ylab(expression(paste(Delta,"MI"))) + 
    theme_bw()
  dev.off()
}else{
  nm.plot = nm[pls_x$V2]
  pls_x_orig = read.csv(paste0('Results/gene_decoding/plsout/pls_x.csv'),header = F)
  pls_xs = data.frame(orig = pls_x_orig$V1, repl = pls_x$V1[match(pls_x_orig$V2, pls_x$V2)]) 
  cor.orig = cor.test(pls_xs$orig, pls_xs$repl)
  if(cor.orig$estimate <0) {
    pls_xs$repl = - pls_xs$repl
  }
  
  brainplot(c(rep(NA, 16),pls_xs$repl),c(rep('A',16),nm.plot),paste0(out_sample, '/gene_decoding/pls1'),-0.2,0.2,'RdBu','stacked')
  
 
  p<-ggplot(pls_xs, aes(x=orig, y=repl)) + 
    geom_point() +
    geom_smooth(method='lm', color='black') +
    xlab('Original PLS1 Weights') + 
    ylab('Replication PLS1 Weights') + 
    ggtitle(substitute(paste(rho," = ",m,",  p = ",p ), list(m=round(abs(cor.orig$estimate),2),p = round(cor.orig$p.value,5))))+
    theme_bw()
  
  pdf(paste0('Results/replication/',sample, '/gene_decoding/corr_pls1.pdf'),3,3)
  plot(p)
  dev.off()
  
}

source('Scripts/helper/FIQT.R')

# # Load PLS Results
# pls1_lena <- read.csv("~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/pls1.csv")
# 
# 
# pls1 <- read.csv("~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/pls1.csv")
# pls1_geneweights<- read.csv("~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/plsout/pls_geneWeights_1.csv")
# matchidx = match(pls1$gene_name, pls1_geneweights$gene_name)
# pls1_geneweights<-pls1_geneweights[matchidx,]
# 
# plot(pls1_geneweights$z_uncorr, pls1$z_uncorr)
# 
# 
# 
# pls1_lena<- read.csv("~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/plsout/pls_geneWeights_1.csv")
# #...
# pls1_lena <- pls1_lena[match(genes$hgnc_symbol,pls1_lena$gene_name),]
# pls1_lena<-pls1_lena[match(pls1$gene_name,pls1_lena$gene_name),]

pls1_lena_orig  <- read.csv("~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/plsout/pls_geneWeights_1.csv")
pls1_lena_orig$z_fdr = FIQT(pls1_lena_orig$z_uncorr)

pls1_lena <- read.csv(paste0('Results/replication/',sample,'/gene_decoding/plsout/pls_geneWeights_1.csv'))
pls1_lena$z_fdr = FIQT(pls1_lena$z_uncorr)

genes = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/pls1_biomart.csv')



in.path = paste0(out_sample, 'gene_decoding/')
out.path = paste0(out_sample, 'gene_decoding/')

# # get chromosome/ensembl info on genes from PLS1
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# genes <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'), 
#                filters = 'hgnc_symbol', 
#                values = pls1_lena$gene_name, 
#                mart = mart, 
#                uniqueRows = T)

genes <- genes[unique(match(pls1_lena$gene_name,genes$hgnc_symbol)),] #reorder genes based on pls1_lena
genes$length = genes$end_position-genes$start_position # calculate gene lengths
genes = genes[-which(rowSums(is.na(genes))>0),]

pls1_lena = cbind(pls1_lena[match(genes$hgnc_symbol, pls1_lena$gene_name),], genes)
pls1_lena_orig = cbind(pls1_lena_orig[match(genes$hgnc_symbol, pls1_lena_orig$gene_name),], genes)

#pls1_lena = pls1_lena[order(pls1_lena$z_fdr),]

# cort.rank = match('CORT',pls1_lena$gene_name)
# sst.rank = match('SST',pls1_lena$gene_name)
# scn1b.rank = match('SCN1B',pls1_lena$gene_name)
# 
# AHBA.gene.name = unlist(AHBA$probeInformation[,,1]$GeneSymbol)
# 
# diff.mi = read.table(paste0('Results/replication/',sample,'/maturational_index/diff.mi.all.txt'))[17:376,1]
# 
# ggplot(data.frame(diff.mi = diff.mi[1:180],AHBA = AHBA$parcelExpression[,match('SST',AHBA.gene.name)]), 
#        aes(y=AHBA,x=diff.mi))+
#   geom_point()+
#   scale_color_manual(values=colormap::colormap('RdBu')[c(1,0.45*72,0.6*72,72)]) +
#   geom_smooth(method = 'lm',formula = y ~ x,color='black') + 
#   ggtitle(paste(pls.z$name[n],'; ', 'r =', rho.val)) +
#   geom_smooth(method='lm')+
#   xlab(expression(paste(Delta,'MI'))) +
#   ylab('AHBA Gene Expression') + 
#   theme_bw()
# 
# 
# gene_selection = c(match('SCN1B',pls1_lena$gene_name), match('CORT',pls1_lena$gene_name), match('SST',pls1_lena$gene_name))
# 
# gene_selection = c('SCN1B','CORT','SST')
# 
# for (i in 1:length(gene_selection)) {
#   df = data.frame(diff.mi = diff.mi[1:180],AHBA = AHBA$parcelExpression[,match(pls1_lena$gene_name[i],AHBA.gene.name)])
#   rho.val = round(cor.test(df$AHBA, df$diff.mi)$estimate,2)
#   # pdf(paste0(plot.out,'plots/examplary_',pls.z$name[n],'.pdf'),height=2.5,width=3)
#   ggplot(df, aes(x=diff.mi, y= AHBA, color=diff.mi)) + 
#     geom_point(size=1) + 
#     #scale_color_manual(values=colormap::colormap('RdBu')[c(1,0.45*72,0.6*72,72)]) +
#     geom_smooth(method = 'lm',formula = y ~ x,color='black') + 
#     ggtitle(paste(pls1_lena$gene_name[i],'; ', 'r =', rho.val)) +
#     theme_bw(base_size = 12) +
#     xlab(expression(paste(Delta,'MI'))) +
#     ylab('AHBA Gene Expression') + 
#     theme(legend.position = "right",
#           legend.title = element_blank(),
#           text = element_text(family = 'Helvetica', size = 12))
#   #dev.off()
# }


if(cor.orig$estimate <0) {
  pls1_lena$z_fdr = -pls1_lena$z_fdr
  pls1_lena$z_uncorr = -pls1_lena$z_uncorr
}

pdf(paste0('Results/replication/',sample,'/gene_decoding/corr_plsweights.pdf'),4.5,3.5)
ggplot(data.frame(orig= pls1_lena_orig$z_uncorr,repl = pls1_lena$z_uncorr[match(pls1_lena_orig$gene_name, pls1_lena$gene_name)]), aes(x = orig, y=repl)) + 
  geom_pointdensity() + 
  scale_color_viridis() + 
  xlab('Replication PLS1 Weights') + ylab('Original PLS1 Weights')+
  ggtitle('Correlation with Original PLS1 Weights')+
  theme_bw() 
dev.off()  



# Calculate chromosome median rank and SE
# 105 vs 26?!
chr_median <- sapply(levels(as.factor(genes$chromosome_name))[-c(23:26)], 
                     function(x) median(rank(pls1_lena$z_uncorr)[which(genes$chromosome_name==x)]))-median(rank(pls1_lena$z_uncorr))
chr_se <- sapply(levels(as.factor(genes$chromosome_name))[-c(23:26)], 
                 function(x) sd(rank(pls1_lena$z_uncorr)[which(genes$chromosome_name==x)])/sqrt(length(which(genes$chromosome_name==x))))
chr_n_genes <- sapply(levels(as.factor(genes$chromosome_name))[-c(23:26)], 
                      function(x) length(which(genes$chromosome_name==x)))

# create null models
source('~/Documents/Education/Cambridge/fMRI_NSPN/Scripts/gene_decoding/gene_nulls_nearest_neighbour.R')

n_perm = 5000

################ Chromosomal Enrichment ################
chr_null = data.frame(gene_set_median_rand=NA,name=NA,real=NA,z=NA,p=NA)
chr_length_stats = data.frame(name=NA,p=NA,s=NA,z=NA, n=NA,approach=NA)
for (chr in levels(as.factor(genes$chromosome_name))[-c(23:26)]) {
  try({
    tmp = gene_nulls_nearest_neighbour(pls1_lena$gene_name, pls1_lena$z_uncorr,genes$hgnc_symbol[which(genes$chromosome_name==chr)],genes,chr,n_perm, rank=T,'~/Desktop/')  
    #tmp = gene_nulls_nearest_neighbour(pls1_lena$gene_name, pls1_lena$z_fdr,genes$hgnc_symbol[which(genes$chromosome_name==chr)],genes,chr,n_perm, rank=T,'~/Desktop/')  
    tmp2 = as.data.frame(tmp[[2]])
    tmp2$name <- chr
    chr_null = rbind(chr_null,tmp[[1]])
    chr_length_stats = rbind(chr_length_stats,tmp2)
  })
} 

chr_null = chr_null[-1,] # First entry was just dummy
chr_null$name = as.factor(chr_null$name) # format data

chr_stats = chr_null[match(unique(chr_null$name),chr_null$name),2:5] # to view stats/needed for fdr-correction

chr_plot <- data.frame(chr=as.factor(names(chr_median)),
                       median=chr_median, # medians to plot
                       se_up=chr_median+chr_se,
                       se_down=chr_median-chr_se,
                       sign=ifelse(p.adjust(chr_stats$p,method='fdr')<0.05,'*','') # FDR-corrected p-values
)
# Plot chromosomal results as bar plot
chr_plot$chr = factor(chr_plot$chr,levels=chr_plot$chr[order(-chr_plot$median)])
pdf(paste0(out.path, 'chromosome.enrichment.pdf'),height=5,width=4)
ggplot(chr_plot,aes(y=median,x=chr))+
  #geom_bar(data=chr_plot, aes(y=median,x=chr,fill=(chr_plot$median+1+abs(min(chr_plot$median)))),stat = 'identity') +
  geom_bar(data=chr_plot, aes(y=median,x=chr,fill=median),stat = 'identity') +
  geom_text(aes(label=sign),hjust=ifelse(chr_plot$median<0,2,-2),vjust=0.75,cex=10) +
  #scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-2*abs(min(chr_plot$median)),2*abs(min(chr_plot$median))))+
  scale_fill_gradientn(colours = c("blue","white","red"),limits=c(-abs(max(chr_plot$median)),abs(max(chr_plot$median))))+
  geom_pointrange(aes(ymin=se_down,ymax=se_up),colour="black")+
  ylim(-2000,2000)+
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

# Plot chromosomal results as joy plots
chr_joy = chr_null[with(chr_null, order(-real)),]
chr_joy$name = factor(chr_joy$name, levels=unique(chr_joy$name))
chr_joy = chr_joy[-which(chr_joy$name=='Y'),]
pdf(paste0(out.path, 'chromosomal.enrichment.null.pdf'),height=5,width=4)
ggplot(chr_joy,aes(x =  gene_set_median_rand, y = name)) +
  geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3,rel_min_height = 0.001) +
  geom_point(aes(y=name,x=real, fill=real),size=4,shape=21) +
  geom_vline(xintercept=0,linetype="dotted",cex=1.1)+
  scale_fill_gradientn(colours = c("blue", "white", "red"), limits=c(-2000,2000))+
  xlab('Median Rank')+
  ylab('') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
dev.off()

# Length stats table
print(xtable(chr_length_stats, booktabs=TRUE), include.rownames=FALSE)

################ Post Natal Cell Type Enrichment ################

lake <- read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Data/gene_decoding/enrichment_analysis/Lake18_celltypes.csv')
lake = lake[-which(lake$Cluster=='Ast_Cer'),]
lake = lake[-which(lake$Cluster=='OPC_Cer'),]
lake = lake[-which(lake$Cluster=='Purk1'),]
lake = lake[-which(lake$Cluster=='Purk2'),]
lake = lake[-which(lake$Cluster=='Per'),]
lake = lake[-which(lake$Cluster=='Gran'),]


lake$Cluster <- as.factor(lake$Cluster)

lake_null = data.frame(gene_set_median_rand=NA,name=NA,real=NA,z=NA,p=NA)
lake_length_stats = data.frame(cluster = NA, p=NA, s=NA,z=NA,n=NA,approach=NA)
for (cluster in unique(lake$Cluster)) {
  try({
    gene_set <- subset(lake, Cluster == cluster)$Gene
    out <- gene_nulls_nearest_neighbour(pls1_lena$gene_name,pls1_lena$z_uncorr,gene_set,genes,cluster,n_perm,rank=T,'~/Desktop')
    lake_null = rbind(lake_null,out[[1]])
    lake_length_stats = rbind(lake_length_stats,c(data.frame(cluster=cluster),out[[2]]))
  })
}   
lake_null = lake_null[-1,] # First entry was just dummy
lake_length_stats = lake_length_stats[-1,]

lake_stats = lake_null[match(unique(lake_null$name),lake_null$name),2:5] # to view stats/needed for fdr-correction

colnames(lake_null) <- c("null","list","real","z","p")
lake_null$list <- as.factor(lake_null$list)

sign.cells = lake_stats$name[p.adjust(lake_stats$p, method='fdr')<0.05]

lake.thresh = lake_null[!is.na(match(lake_null$list,sign.cells)),]
lake.thresh$class <- 'Other'
lake.thresh$class[grep('Ex',lake.thresh$list)] <- 'Exitatory'
lake.thresh$class[grep('In',lake.thresh$list)] <- 'Inhibitory'

lake.plot = lake.thresh[with(lake.thresh, order(class, real)),]
lake.plot$list = factor(lake.plot$list, levels=unique(lake.plot$list))

pdf(paste0(out.path, 'postnatal.celltype.enrichment.pdf'),height=3,width=4)
ggplot(lake.plot,aes(x = null, y = list)) +
  geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3,rel_min_height = 0.001,bandwidth=150) +
  geom_point(aes(y=list,x=real, fill=real),size=4,shape=21) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), lim = c(-3000,3000),oob=squish)+
  xlab('Median Rank')+
  ylab('') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
dev.off()

print(xtable(lake_length_stats, booktabs=TRUE), include.rownames=FALSE)


################ Pre Natal Cell Type Enrichment ################
polioudakis <- read.table('~/Documents/Education/Cambridge/fMRI_NSPN/Data/gene_decoding/enrichment_analysis/Polioudakis_celltypes.txt',header = T)
polioudakis$Cluster <- as.factor(polioudakis$Cluster)

polioudakis_null = data.frame(gene_set_median_rand=NA,name=NA,real=NA,z=NA,p=NA)
polioudakis_length_stats = data.frame(cluster = NA, p=NA, s=NA,z=NA,n=NA,approach=NA)
for (cluster in unique(polioudakis$Cluster)) {
  try({
    gene_set <- subset(polioudakis, Cluster == cluster)$Gene
    out <- gene_nulls_nearest_neighbour(pls1_lena$gene_name,pls1_lena$z_uncorr,gene_set,genes,cluster,n_perm,rank=T,'~/Desktop')
    polioudakis_null = rbind(polioudakis_null,out[[1]])
    polioudakis_length_stats = rbind(polioudakis_length_stats,c(data.frame(cluster=cluster),out[[2]]))
  })
}   
polioudakis_null = polioudakis_null[-1,] # First entry was just dummy
polioudakis_length_stats = polioudakis_length_stats[-1,]

polioudakis_stats = polioudakis_null[match(unique(polioudakis_null$name),polioudakis_null$name),2:5] # to view stats/needed for fdr-correction

colnames(polioudakis_null) <- c("null","list","real","z","p")
polioudakis_null$list <- as.factor(polioudakis_null$list)

sign.cells = polioudakis_stats$name[p.adjust(polioudakis_stats$p, method='fdr')<0.05]

polioudakis_thresh = polioudakis_null[!is.na(match(polioudakis_null$list,sign.cells)),]
polioudakis_thresh$class <- 'Other'
polioudakis_thresh$class[grep('Ex',polioudakis_thresh$list)] <- 'Exitatory'
polioudakis_thresh$class[grep('In',polioudakis_thresh$list)] <- 'Inhibitory'

polioudakis_plot = polioudakis_thresh[with(polioudakis_thresh, order(class, real)),]
polioudakis_plot$list = factor(polioudakis_plot$list, levels=unique(polioudakis_plot$list))

p<-ggplot(polioudakis_plot,aes(x = null, y = list)) +
  geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3,rel_min_height = 0.001,bandwidth=100) +
  geom_point(aes(y=list,x=real, fill=real),size=4,shape=21) +
  scale_fill_gradientn(colours = c("blue", "white", "red"), lim = c(-3000,3000),oob=squish)+
  xlab('Median Rank')+
  ylab('') +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())

pdf(paste0(out.path, 'prenatal.celltype.enrichment.pdf'),height=5,width=4)
print(p)
dev.off()

print(xtable(polioudakis_length_stats, booktabs=TRUE), include.rownames=FALSE)

################ MDD Enrichment ################
howard = read.table("~/Documents/Education/Cambridge/fMRI_NSPN/Data/gene_decoding/enrichment_analysis/Howard_MDD_GWAScatalog.txt")
median(rank(pls1_lena$z_uncorr)[match(howard$V1,pls1_lena$gene_name,nomatch=0)])-median(rank(pls1_lena$z_uncorr))

out <- gene_nulls_nearest_neighbour(pls1_lena$gene_name,pls1_lena$z_uncorr,howard$V1,genes,'Howard',n_perm,rank=T,'~/Desktop')

howard_null <- out[[1]]
howard_stat <- out[[2]]

howard_density <- density(howard_null$gene_set_median_rand)
howard_plot <- data.frame(x = howard_density$x, y = howard_density$y, fill=howard_density$x, real=unique(howard_null$real))

pdf(paste0(out.path, 'mdd.enrichment.pdf'),height=1.75,width=4)
ggplot(howard_plot) +
  geom_segment(aes(x = x, xend = x, y = 0, yend = y, color = fill)) +
  geom_line(aes(x=x, y=y)) +
  geom_point(aes(x=real, y=0,fill=real),size=4,shape=21) +
  scale_color_gradientn(colours = c("blue", "white", "red"),
                        limits=c(-3000,3000))+
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       limits=c(-3000,3000))+
  theme_minimal() +
  xlab('Median Rank')+   
  theme(legend.position="none",
        axis.text.x = element_text(size=12),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank())
dev.off()






#### X Chromosome Analysis ####
howard = read.table("~/Documents/Education/Cambridge/fMRI_NSPN/Data/gene_decoding/enrichment_analysis/Howard_MDD_GWAScatalog.txt")
MDD_genes = howard$V1
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
MDD_genes <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position'), 
                   filters = 'hgnc_symbol', values = MDD_genes, mart = mart, uniqueRows = T)

ggplot(MDD_genes, aes(x=forcats::fct_infreq(MDD_genes$chromosome_name)))+
  geom_bar()+
  ggtitle('Chromosomal Locations of Howard 2017 Genes')+
  theme_bw()+
  theme(axis.text.x =element_text(angle = 90),
        axis.title = element_blank())

XChr = genes[genes$chromosome_name=='X',]

AHBA=readMat('~/Documents/Education/Cambridge/Datasets/AllenBrainAtlas/Petra/AHBAdata.mat')
gene_names<-unlist(AHBA$probeInformation[[2]])

#mu_overall = mu_overall[which(!is.na(mu_overall))]
MDD_expression = AHBA$parcelExpression[,match(howard$V1, gene_names)]
#MDD_expression = MDD_expression[,which(!is.na(MDD_expression[1,]))]
XChr_expression = AHBA$parcelExpression[,match(XChr$hgnc_symbol,gene_names)]

#delX = which(rowSums(is.na(XChr_expression))>0)
#delY = which(rowSums(is.na(MDD_expression))>0)
#del = unique(c(delX, delY))

#XChr_expression = XChr_expression[-del,]
#MDD_expression = MDD_expression[-del,]

mu_exp_MDD = apply(MDD_expression, 1, function(x) mean(x, na.rm=T))
mu_exp_XChr = apply(XChr_expression, 1, function(x) mean(x, na.rm=T))
mu_overall = apply(AHBA$parcelExpression,1,function(x) mean(x,na.rm=T))
mu_overall[is.na(mu_exp_MDD)] = NA

load('~/Documents/Education/Cambridge/fMRI_NSPN/Data/connectivity_matrices/nspn_fmri_longitudinal_hcp_fdreg.RData')
diff.mi = read.table('~/Documents/Education/Cambridge/fMRI_NSPN/Results/maturational_index/diff.mi.txt')[,1]
hcp.all = rep(NA,1,376)
hcp.all[hcp.keep.id] = diff.mi
diff.mi = hcp.all[17:(17+179)]

#label = nm.all[seq(1,180)[-del]+16]
MDDcor = vector(length=dim(MDD_expression)[2])
MDD_lm = matrix(nrow = dim(MDD_expression)[2], ncol = 2)
for (i in 1:dim(MDD_expression)[2]) {
  try({
    df = data.frame(mi=diff.mi, gene_expr = (MDD_expression[,i]-mu_overall))
    l = lm(mi~gene_expr, data=df)
    MDD_lm[i,] = l$coefficients
  })
  tmp = try(print(cor.test(diff.mi, (MDD_expression[,i]-mu_overall))$estimate))
  if (is.numeric(tmp)) {
    MDDcor[i] = tmp
  }
}

MDD_lm = data.frame(MDD_lm[!is.na(MDD_lm[,1]),])
#MDD_lm = as.data.frame(MDD_lm)
colnames(MDD_lm) <- c('intercept','slope')
hist(MDDcor[!MDDcor==0])

XChrcor = vector(length=dim(XChr_expression)[2])
XChr_lm = matrix(nrow = dim(XChr_expression)[2], ncol = 2)
for (i in 1:dim(XChr_expression)[2]) {
  try({
    df = data.frame(mi=diff.mi, gene_expr = (XChr_expression[,i]-mu_overall))
    l = lm(mi~gene_expr, data=df)
    XChr_lm[i,] = l$coefficients
  })
  tmp = try(print(cor.test(diff.mi, (XChr_expression[,i]-mu_overall))$estimate))
  if (is.numeric(tmp)) {
    XChrcor[i] = tmp
  }
}
hist(XChrcor)

df = data.frame(x=diff.mi,y=mu_overall)
ggplot(df, aes(x=x, y=y)) +
  geom_point()+
  geom_abline(slope = MDD_lm[1,2] , intercept=MDD_lm[1,1])+
  geom_abline(slope = MDD_lm[2,2] , intercept=MDD_lm[2,1])+
  geom_abline(slope = MDD_lm[3,2] , intercept=MDD_lm[3,1])+
  geom_abline(slope = MDD_lm[4,2] , intercept=MDD_lm[4,1])

#cor.test(diff.mi[match(label,nm)],(mu_exp_MDD[match(label,nm)]-mu_overall), na.rm=T)
#cor.test(diff.mi[match(label,nm)],(mu_exp_XChr[match(label,nm)]-mu_overall), na.rm=T)

df = data.frame(MI = diff.mi, MDD = (mu_exp_MDD-mu_overall), XChr = (mu_exp_XChr-mu_overall) )
#df.MDD = data.frame(mi = diff.mi[match(label,nm)], mdd = (mu_exp_MDD[match(label,nm)]-mu_overall))
#df.xchr = data.frame(mi = diff.mi[match(label,nm)], xchr = (mu_exp_XChr[match(label,nm)]-mu_overall))

cor.MDD = cor.test(df$MI, df$MDD)
cor.XChr = cor.test(df$MI, df$XChr)

ggplot(df, aes(x=MI, y=MDD)) +geom_point() +geom_smooth(method='lm')+ 
  ggtitle(paste0('r = ',round(cor.MDD$estimate,2),'; p = ',round(cor.MDD$p.value,2)))+ 
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size=15,hjust=0.5))

ggplot(df, aes(x=MI, y=XChr)) +geom_point() +geom_smooth(method='lm')+ 
  ggtitle(paste0('r = ',round(cor.XChr$estimate,2),'; p = ',round(cor.XChr$p.value,2)))+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size=15,hjust=0.5))

label = rep(NA,1,376)
label[hcp.keep.id] = nm
label[c(1:16,(17+180):376)] = NA

plot.MDD = rep(NA,1,376)
plot.MDD[17:(17+179)] = df$MDD

plot.X = rep(NA,1,376)
plot.X[17:(17+179)] = df$XChr

plot.X[plot.X=='NaN']=NA
plot.MDD[plot.MDD=='NaN']=NA

source('~/Desktop/ggseg_usage/brainplot.R')
load('~/Desktop/ggseg_usage/glassersubcortex.bin')

brainplot(plot.X,label,'~/Desktop/tmp',-0.06, 0.06, 'seismic','stacked',TRUE)
brainplot(plot.MDD,label,'~/Desktop/tmp',-0.06, 0.06, 'seismic','stacked',FALSE)

diff.mi = read.table('~/Documents/Education/Cambridge/fMRI_NSPN/Results/maturational_index/diff.mi.txt')[,1]
brainplot(diff.mi,nm,'~/Desktop/tmp',-0.8, 0.8, 'seismic','stacked',FALSE)


#ggplot(df.xchr, aes(x=mi, y=xchr)) +geom_point() +geom_smooth(method='lm') + theme_bw()
#ggplot(df.MDD, aes(x=mi, y=mdd)) +geom_point() +geom_smooth(method='lm')+ theme_bw()





# library(plsdepot)
# source('~/Documents/Education/Cambridge/Katja/PLS_Paper/plot_perm_results.R')
# source('~/Documents/Education/Cambridge/Katja/PLS_Paper/plot_variance_explained.R')
# 
# ncomp = 2
# pls.model = plsreg2(X,Y, crosval = T, comps = ncomp)
# plot_variance_explained(pls.model$expvar[,2],pls.model$expvar[,4],'~/Desktop/')
# 
# nperm = 100 # Define number of permutation iterations
# perm.var.expl = matrix(NA, ncol=ncomp, nrow=nperm)
# #pb = txtProgressBar(min = 0, max = length(nperm), style=3)
# for (perm in 1:nperm) {
#   print(perm)
#   # Permute Y randomly
#   order = sample(1:dim(Y)[1]) 
#   # Recalculate PLS with permuted Y
#   perm.model = plsreg2(X[order,], Y, comps = ncomp, crosval = TRUE)
#   # Save explained variance in Y
#   perm.var.expl[perm,] = perm.model$expvar[,3]
# }
# # Calculate PLS component p-values
# component.pvalues = vector(length = ncomp)
# for (c in 1:ncomp) {
#   # Is the real variance explained higher than in the permutation distribution?
#   component.pvalues[c] = sum(pls.model$expvar[c,3]<perm.var.expl[,c])/nperm
# }
# 
# source('~/Desktop/ggseg_usage/brainplot.R')
# load('~/Desktop/ggseg_usage/glassersubcortex.bin')
# source('~/Documents/Education/Cambridge/general_scripts/brainplot.R')
# load('~/Documents/Education/Cambridge/general_scripts/glassersubcortex_R.bin')
# 
# 
# nm = nm.all[seq(1,180)[-del]+16]
# label = vector(length=length(glassersubcortex$label))
# label[match(nm,unique(glassersubcortex$label))] = nm
#   
# comp1 = vector(length=length(label))
# comp1[match(nm,unique(glassersubcortex$label))] = pls.model$x.scores[,1]
# 
# brainplot(comp1,label,'~/Desktop/tmp',-30, 30, 'seismic','stacked',TRUE)


