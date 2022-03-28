library(scales) # Version 1.1.1 
library(ggplot2) # Version 3.3.5
library(ggseg) # Version 1.5.4

data = 'NSPN' # data = 'CV_stratified' # data = 'glob.fc.regr' # data = 'regr_FD_by_sex' # data = 'motionmatched' # data = # data = 
script.path = paste0('Scripts/')

data.path = paste0('Data/connectivity_matrices/')
results.path = paste0('Results/')
# load in data: NSPN, NSPN low motion or GSR
if (data == 'NSPN') {
  load(paste0(data.path,'nspn.main.RData'))
  plot.out = paste0(results.path,'mdd_decoding/')
  data.path = paste0(results.path,'maturational_index/')
}else{
  load(paste0(data.path, 'nspn.motion.matched.RData'))
  plot.out = paste0(results.path,'replication/', sample,'/mdd_decoding/')
  data.path = paste0(results.path,'replication/', sample,'/maturational_index/')
}

ifelse(!dir.exists(file.path(plot.out)), dir.create(file.path(plot.out)), FALSE)

# Load in MDD case control map, generated from BioDep Data
p.case.control = read.table('Data/mdd/mdd.case.control.p.txt')[,1]
p.fdr.case.control = p.adjust(p.case.control, method='fdr')
t.case.control = read.table('Data/mdd/mdd.case.control.t.txt')[,1]

# Load in Delta MI
diff.mi = read.table(paste0(data.path, 'diff.mi.txt'))[[1]]
p.diff.MI = read.table(paste0(data.path, 'z.p.fdr.MI.txt'))[,1]

# Case control map subcortex
df.aseg$value = t.case.control[aseg.nm.idx]
pdf(paste0(plot.out,'/mdd.case.control.subc.pdf'))
ggseg(.data=df.aseg, mapping=aes(fill=value), atlas = "aseg", view = 'axial')+ 
  scale_fill_gradientn(colours= viridis_pal()(1000), limits = c(-4,1), oob=squish)   
dev.off()

# Case control map cortex
pdf(paste0(plot.out,'/mdd.case.control.cort.pdf'))
ggseg(.data=data.frame(value = t.case.control, label = nm.ggseg), mapping=aes(fill=value), atlas = 'glassersub', position='stacked')+ 
  scale_fill_gradientn(colours = viridis_pal()(1000), limits=c(-4,1), breaks = c(-4,1))+theme_void()+labs(fill='')
dev.off()

mi.colors = seismic.colorscale[round(rescale(diff.mi, to = c(1,1000)))]

# Significance of correlation between MDD case control difference and delta MI
t.test = cor.test(diff.mi,t.case.control)

df=data.frame(mi=diff.mi, mdd=t.case.control, color = ifelse(p.diff.MI<0.05,1,0),shape = ifelse(p.case.control<0.05,21,19))
pdf(paste0(plot.out, 'mdd.case.control.diff.mi.pdf'),2.3,2.7)
ggplot(df, aes(x=mi, y=mdd)) + 
  geom_point(fill = mi.colors,alpha=1, shape=df$shape) + 
  scale_color_discrete(c('red','black')) +
  geom_smooth(method='lm', color='black') + 
  ggtitle(substitute(paste(rho,'=',corr,', P=',p),list(corr=round(t.test$estimate,2),p=round(t.test$p.value,3))))+
  theme_bw() +
  xlab(expression(paste(Delta,'MI')))+
  ylab("MDD Case Control") +
  theme(text = element_text(size=12, family ="Helvetica"),
        plot.title = element_text(size=14,hjust=0.5,face='italic'))
dev.off()


# Spin test p-value for correlation
# Download spin test function from: https://github.com/frantisekvasa/rotate_parcellation
source('Scripts/external/rotate_parcellation-master/R/perm.sphere.p.R') 
p.spin.mdd = perm.sphere.p(df$mdd[17:346],df$mdd[17:346],perm.id.330)

print(paste0('MDD-Delta MI Correlation: ',round(cor.test(df$mdd,df$mi)$estimate,4), ' and p-value: ',round(cor.test(df$mdd,df$mi)$p.value,4)))
print(paste0('MDD-Delta MI Correlation spin p-value: ',round(p.spin.mdd,4)))


