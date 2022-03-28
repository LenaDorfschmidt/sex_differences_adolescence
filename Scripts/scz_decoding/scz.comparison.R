library(scales) # Version 1.1.1 
library(ggplot2) # Version 3.3.5
library(ggseg) # Version 1.5.4

data = 'NSPN' # data = 'CV_stratified' # data = 'glob.fc.regr' # data = 'regr_FD_by_sex' # data = 'motionmatched' # data = # data = 
script.path = paste0('Scripts/')

data.path = paste0('Data/connectivity_matrices/')
results.path = paste0('Results/')

# load in data
load(paste0(data.path,'nspn.main.RData'))
plot.out = paste0(results.path,'scz_decoding/')
data.path = paste0(results.path,'maturational_index/')

diff.mi = read.table(here('Results/gene_decoding/diff.mi.all.txt'), header=T)$x
p.diff.MI = read.table('Results/maturational_index/z.p.fdr.MI.txt')[,1]

ifelse(!dir.exists(file.path(plot.out)), dir.create(file.path(plot.out)), FALSE)

# Load in SCZ case control map, generated from Cobre Data
p.case.control = read.table('Data/scz/scz.case.control.p.txt')[,1]
p.fdr.case.control = p.adjust(p.case.control, method='fdr')
t.case.control = read.table('Data/scz/scz.case.control.t.txt')[,1]

# Load in Delta MI
diff.mi = read.table(paste0(data.path, 'diff.mi.all.txt'))[[1]]
p.fdr.diff.MI = read.table(paste0(data.path, 'z.p.fdr.MI.all.txt'))[,1]

mi.colors = seismic.colorscale[round(rescale(diff.mi, to = c(1,1000)))]

# Correlation between SCZ case control map and delta MI
t.test = cor.test(diff.mi,t.case.control)

df=data.frame(mi=diff.mi, scz=t.case.control, color = ifelse(p.fdr.diff.MI<0.05,1,0),shape = ifelse(p.fdr.diff.MI<0.05,21,19))
pdf(paste0(plot.out, 'scz.case.control.diff.mi.pdf'),2.3,2.7)
ggplot(df, aes(x=mi, y=scz)) + 
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
df.no.na = df[-which(is.na(df$mi)),]
p.spin.scz = perm.sphere.p(df.no.na$scz[17:346],df.no.na$mi[17:346],perm.id.330)

print(paste0('SCZ-Delta MI Correlation: ',round(cor.test(df$scz,df$mi)$estimate,4), ' and p-value: ',round(cor.test(df$scz,df$mi)$p.value,4)))
print(paste0('SCZ-Delta MI Correlation spin p-value: ',round(p.spin.scz,4)))

