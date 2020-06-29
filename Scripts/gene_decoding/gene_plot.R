library(ggplot2)
library(scales)
# BrainSpan
neg.brainspan = read.csv('~/Desktop/negative.weights.CSEA.development.csv',header=T)
neg.brainspan$Period = gsub('_',' ',neg.brainspan$Period)
pos.brainspan = read.csv('~/Desktop/positive.weights.CSEA.development.csv',header=T)
pos.brainspan$Period = gsub('_',' ',pos.brainspan$Period)


neg.brainspan = neg.brainspan[order(neg.brainspan$Period_Label),]
neg.brainspan$Period = factor(neg.brainspan$Period, levels = unique(neg.brainspan$Period))
pos.brainspan = pos.brainspan[order(pos.brainspan$Period_Label),]
pos.brainspan$Period = factor(pos.brainspan$Period, levels = unique(pos.brainspan$Period))

# brainspan.joint = cbind(pos.brainspan[,c(1,2,3,5)],neg.brainspan$BH_05)
# colnames(brainspan.joint)[4:5] = c('positive','negative')
# 
# ggplot(brainspan.joint, aes(x=Period, y=ROI)) + 
#   geom_tile(alpha=-log10(brainspan.joint$positive),fill="red") +
#   geom_tile(alpha=-log10(brainspan.joint$negative),fill="blue") +
#   theme_bw() +
#   ylab('Brain Region')+
#   xlab('Developmental Period') +
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1,size=12),
#         axis.text.y = element_text(size=12),
#         axis.title.x = element_text(size=15),
#         axis.title.y = element_text(size=15))


brainspan.joint = rbind(pos.brainspan[,c(1,2,3,5)],neg.brainspan[,c(1,2,3,5)])
brainspan.joint$color = c(log10(pos.brainspan$BH_05),-log10(neg.brainspan$BH_05))
levels(brainspan.joint$Period)[c(2,3,5,8)] = c('Early Fetal','Late Fetal','Early Infancy','Late Childhood')

pdf('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/CSEA.development.pdf',height = 4, width = 4)
ggplot(brainspan.joint, aes(x=Period, y=ROI,fill=color,alpha=0.99)) + 
    geom_tile(stat="identity") +
    scale_fill_gradientn(colours = c("red",'white',"blue"), limits=c(-20,20),oob=squish)+ 
    theme_bw() +
   # ylab('Brain Region')+
   # xlab('Developmental Period') +
    theme(legend.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1,size=12),
          axis.text.y = element_text(size=12),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
dev.off()
  
# neg.brainspan = neg.brainspan[order(neg.brainspan$Period_Label),]
# neg.brainspan$Period = factor(neg.brainspan$Period, levels = unique(neg.brainspan$Period))
# ggplot(neg.brainspan, aes(x=Period, y=ROI,fill=-log10(neg.brainspan$BH_05),alpha=0.3)) + 
#   geom_tile(stat="identity") +
#   scale_fill_gradientn(colours = c("white","blue"))+
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1))
# 
# 
# pos.brainspan = pos.brainspan[order(pos.brainspan$Period_Label),]
# pos.brainspan$Period = factor(pos.brainspan$Period, levels = unique(pos.brainspan$Period))
# ggplot(pos.brainspan, aes(x=Period, y=ROI,fill=-log10(pos.brainspan$BH_05),alpha=0.3)) + 
#   geom_tile(stat="identity") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   scale_fill_gradientn(colours = c("white","red")) +
#   theme_bw() +
#   theme(legend.title = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1))


celltype.neg = read.csv('~/Desktop/CSEA_celltypes_pls_neg.csv')
celltype.pos = read.csv('~/Desktop/CSEA_celltypes_pls_pos.csv')

celltype.joint = rbind(celltype.neg[,c(1,2,3,5)],celltype.pos[,c(1,2,3,5)])
celltype.joint$color = c(-log10(celltype.neg$BH_0.05),log10(celltype.pos$BH_0.05))
celltype.joint = celltype.joint[celltype.joint$BH_0.05<0.05,]
celltype.joint = celltype.joint[!(celltype.joint$celltype==''),]
celltype.joint = celltype.joint[order(celltype.joint$structure),]

pdf('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/CSEA.celltype.pdf',height = 4, width = 4)
ggplot() + 
  geom_bar(data=celltype.joint, 
           aes(x=factor(celltype.joint$original,levels=unique(celltype.joint$original)),
               y=color, fill=color),stat = 'identity') +
  scale_fill_gradientn(colours = c("red","white","blue")) +
  coord_flip() +
  scale_x_discrete(labels=celltype.joint$celltype) +
  ylab(expression("Significance Level log"[10]*"(p-value)"))+
  xlab('Cell Type')+
  theme_bw()+
  theme(legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_blank())
dev.off()


###
source('~/Documents/Education/Cambridge/general_scripts/brainplot.R')
source('~/Documents/Education/Cambridge/fMRI_NSPN/Scripts/helper/forSurfacePlot.R')
load('~/Documents/Education/Cambridge/general_scripts/glassersubcortex_R.bin')
load('~/Documents/Education/Cambridge/fMRI_NSPN/Data/connectivity_matrices/nspn_fmri_longitudinal_hcp_fdreg.RData')

CORT = read.table('~/Desktop/AHBA_expression_CORT.txt')
NPY = read.table('~/Desktop/AHBA_expression_NPY.txt')
SST = read.table('~/Desktop/AHBA_expression_SST.txt')
  
brainplot(CORT[hcp.keep.id[17:346]-16,1],nm[17:346],'~/Desktop/CORT',0,1,'viridis','stacked')
brainplot(NPY[hcp.keep.id[17:346]-16,1],nm[17:346],'~/Desktop/NPY',0,1,'viridis','stacked')
brainplot(SST[hcp.keep.id[17:346]-16,1],nm[17:346],'~/Desktop/SST',0,1,'viridis','stacked')

#### Examplary Plots

AHBA = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/AHBA_PLS_X_data.csv', header = F)
gene_names = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/AHBA_gene_names.csv', header = F)[[1]]
gene_names = as.character(gene_names[2:length(gene_names)])
pls.z = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Results/gene_decoding/plsout/PLS_geneWeights_1_fdr.corrected.csv')
pls.z$name = as.character(pls.z$name)


plot.out = '~/Documents/Education/Cambridge/fMRI_NSPN/Results/maturational_index/'
mi.m = read.csv(paste0(plot.out,'/male/cc_mat_idx.txt'))[,2]
mi.f = read.csv(paste0(plot.out,'/female/cc_mat_idx.txt'))[,2]

diff.mi = mi.f-mi.m
mu.diff.mi = rep(NA,1,376)
mu.diff.mi[hcp.keep.id] = diff.mi
mu.diff.mi = mu.diff.mi[17:(17+179)]#rowMeans(cbind(mu.diff.mi[17:(17+179)],mu.diff.mi[(17+180):376]))

roi.trends = read.csv(paste0(plot.out,'trends.txt'),header=F)[,1]
trend = rep(NA, 1,376)
trend[hcp.keep.id] = roi.trends
trend = trend[17:(17+179)]

gene_selection = c(which(pls.z$name=='CORT'),which(pls.z$name=='NPY'),which(pls.z$name=='SST'),1,2)

for (i in 1:length(gene_selection)) {
  n = gene_selection[i]
  idx = match(pls.z$name[n],gene_names)
  df = data.frame(AHBA=AHBA[,idx],diff.mi=mu.diff.mi,trend=as.factor(trend))
  rho.val = round(cor.test(df$AHBA, df$diff.mi)$estimate,2)
  p <- ggplot(df, aes(x=diff.mi, y= AHBA, color=trend)) + 
    geom_point(size=1) + 
    scale_color_manual(values=colormap::colormap('RdBu')[c(1,0.45*72,0.6*72,72)]) +
    geom_smooth(method = 'lm',formula = y ~ x,color='black') + 
    ggtitle(paste(pls.z$name[n],'; ', 'r =', rho.val)) +
    theme_bw(base_size = 12) +
    xlab(expression(paste(Delta,'MI'))) +
    ylab('AHBA Gene Expression') + 
    theme(legend.position = "right",
          legend.title = element_blank(),
          text = element_text(family = 'Helvetica', size = 12))
  pdf(paste0(plot.out,'plots/examplary_',pls.z$name[n],'.pdf'),height=2.5,width=3)
  plot(p)
  dev.off()
}


pls_comp = read.table('~/Desktop/b4p.pls.xs.txt')
brainplot(pls_comp[hcp.keep.id[17:346]-16,1],nm[17:346],'~/Desktop/pls_comp',-0.18,0.18,'RdBu','stacked')
brainplot(roi.trends[hcp.keep.id[17:346]-16,1],nm[17:346],'~/Desktop/roi.trends',1,4,'RdBu','stacked')
brainplot(c(AHBA[,idx],AHBA[,idx])[hcp.keep.id[17:346]-16],nm[17:346],'~/Desktop/SNC1B',-3,3,'viridis','stacked')


