library(here) # Version 1.0.0 
library(corrplot) # Version 0.84

############################################################
#            All Motion Sensitivity Analyses              #
############################################################
orig.path = paste0('Results/maturational_index/')
load(paste0('Data/connectivity_matrices/nspn.main.RData'))


# Supplementary Figure S42A
repl = c('orig','motionmatched','regr_FD_by_sex','nspn.gsr','CV_stratified','glob.fc.regr')
stats = c('all_bl14','all_delta_age','cc_mat_idx')
sex = c('female','male')
mult.fact = 6
for (i in 1:length(stats)) {
  repl.stats = matrix(ncol=nroi, nrow = mult.fact*2)
  for (s in 1:2) {
    idx = mult.fact * (s-1)
    repl.stats[idx+1,] = read.csv(paste0('Results/maturational_index/',sex[s],'/fc14-fc1426/',stats[i],'.txt'))[,2]
    for (j in 2:length(repl)) {
      repl.stats[j+idx,] = read.csv(paste0('Results/replication/',repl[j],'/maturational_index/',sex[s],'/fc14-fc1426/',stats[i],'.txt'))[,2]
    }
  }
  cor.mat = cor(t(repl.stats))
  colnames(cor.mat) = c(paste0(sex[1],paste0('_',repl)),paste0(sex[2],paste0('_',repl)))
  rownames(cor.mat) = c(paste0(sex[1],paste0('_',repl)),paste0(sex[2],paste0('_',repl)))
  
  pdf(paste0('Results/replication/corr.mat.all.',stats[i],'.pdf'),7,7)
  corrplot(cor.mat,type="lower",diag = FALSE,col=rev(brewer.pal(11,'RdBu')),number.digits = 1,addCoef.col = 'white',tl.pos=0)
  dev.off()
}

# Supplementary Figure S42B
repl.stats.mi = matrix(ncol=nroi, nrow = mult.fact)
repl.stats.mi[1,] = read.table(paste0('Results/maturational_index/diff.mi.txt'))[,1]
for (r in 2:length(repl)) {
  repl.stats.mi[r,] = read.table(paste0('Results/replication/',repl[r],'/maturational_index/diff.mi.txt'))[,1]
}

cor.mat.mi = cor(t(repl.stats.mi))
row.names(cor.mat.mi)=  repl
pdf(paste0('Results/replication/corr.mat.diff.mi.all.pdf'),4,4)
corrplot(cor.mat.mi,type="lower",diag = FALSE,number.digits = 2,addCoef.col = 'white',tl.pos=0,col = brewer_pal(palette='RdBu',direction=-1)(10))
dev.off()


############################################################
#            Individual Sensitivity Analyses              #
############################################################

sample = 'glob.fc.regr' #sample = 'GSR' 

if(sample == 'GSR'){
  repl.path = paste0('Results/replication/gsr/')
  load(paste0(base.path, '/Data/connectivity_matrices/gsr.RData'))
  l1 = -0.1
  h1 = 0.1
  l2 = -0.005
  h2=0.005
}else if(sample == 'regr_FD_by_sex'){
  repl.path = paste0('Results/replication/regr_FD_by_sex/')
  load(paste0('Data/connectivity_matrices/regr_FD_by_sex.RData'))
  l1=0.2
  h1=0.8
  l2=-0.015
  h2=0.015
}else if(sample == 'motionmatched'){
  repl.path = paste0(base.path, 'Results/replication/motionmatched/')
  load(paste0('Data/connectivity_matrices/motionmatched.RData'))
  l1=0.2
  h1=0.8
  l2=-0.015
  h2=0.015
}else if(sample == 'glob.fc.regr'){
  repl.path = paste0('Results/replication/glob.fc.regr/maturational_index/')
  load(paste0(here('Data/connectivity_matrices/nspn.main.RData')))
  l1=-0.01
  h1=0.8
  l2=-0.02
  h2=0.02
}else if(sample == 'CV_stratified'){
  repl.path = paste0('Results/replication/CV_stratified/maturational_index/')
  load(paste0(here('Data/connectivity_matrices/nspn.main.RData')))
  l1=0.1
  h1=0.8
  l2=-0.02
  h2=0.02
}
if(sample == 'interaction_1'){
  repl.path = paste0('Results/replication/interaction_1/')
  load(paste0('Data/connectivity_matrices/nspn.main.RData'))
  l1=0.2
  h1=0.8
  l2=-0.015
  h2=0.015
}


sex = c('female','male')
for (s in sex) {
  orig.bl14 = read.csv(paste0(orig.path,s,'/all_bl14.txt'))[,2]
  repl.bl14 = read.csv(paste0(repl.path,s,'/all_bl14.txt'))[,2]
  orig.sl = read.csv(paste0(orig.path,s,'/all_delta_age.txt'))[,2]
  repl.sl = read.csv(paste0(repl.path,s,'/all_delta_age.txt'))[,2]
  
  brainplot(repl.bl14,nm,paste0(repl.path,s,'/all_bl14'),l1,h1,'plasma','stacked')
  brainplot(repl.sl,nm,paste0(repl.path,s,'/all_delta_age'),l2,h2,'seismic','stacked')
  
  
  cor=cor.test(orig.bl14, repl.bl14)
  g<-ggplot(data.frame(orig=orig.bl14,repl=repl.bl14),aes(x=orig,y=repl))+
    geom_point()+
    geom_smooth(method = 'lm',color='black')+
    ggtitle(substitute(paste(rho,'=',corr,', P=',p),list(corr=round(cor$estimate,2),p=round(cor$p.value,5))))+
    xlab(expression('Original ' * BL[14]))+
    ylab(expression('Replication ' *  BL[14])) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          text= element_text(size=12, family='Helvetica'))
  pdf(paste0(repl.path, s,'/corr.plot.bl14.pdf'),height=3.5, width = 3.5)
  plot(g)
  dev.off()
  
  cor=cor.test(orig.sl, repl.sl)
  g <-ggplot(data.frame(orig=orig.sl,repl=repl.sl),aes(x=orig,y=repl))+
    geom_point()+
    geom_smooth(method = 'lm',color='black')+
    ggtitle(substitute(paste(rho,'=',corr,', P=',p),list(corr=round(cor$estimate,2),p=round(cor$p.value,5))))+
    xlab(expression('Original ' * FC[14-26]))+
    ylab(expression('Replication ' *  FC[14-26])) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          text= element_text(size=12, family='Helvetica'))
  pdf(paste0(repl.path, s,'/corr.plot.sl.pdf'),height=3.5, width = 3.5)
  plot(g)
  dev.off()
}



# Correlation MI
orig.diff.mi = read.table(paste0(orig.path, 'diff.mi.txt'))[,1]
repl.diff.mi = read.table(paste0(repl.path,'diff.mi.txt'))[,1]

cor=cor.test(repl.diff.mi, orig.diff.mi)
pdf(paste0(repl.path, 'corr.plot.orig.pdf'),height=3.5, width = 3.5)
ggplot(data.frame(orig=orig.diff.mi,repl=repl.diff.mi),aes(x=orig,y=repl))+
  geom_point()+
  geom_smooth(method = 'lm',color='black')+
  ggtitle(substitute(paste(rho,'=',corr,', P=',p),list(corr=round(cor$estimate,2),p=round(cor$p.value,5))))+
  xlab(expression('Original ' * Delta * 'MI'))+
  ylab(expression('Replication ' * Delta * 'MI')) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
dev.off()

orig.z = read.table(paste0(orig.path,'/z.MI.txt'))[,1]
repl.z = read.table(paste0(repl.path,'/z.MI.txt'))[,1]

cor=cor.test(orig.z, repl.z)
pdf(paste0(repl.path, 'corr.plot.z.pdf'),height=3.5, width = 3.5)
ggplot(data.frame(orig=orig.z,repl=repl.z),aes(x=orig,y=repl))+
  geom_point()+
  geom_smooth(method = 'lm',color='black')+
  ggtitle(substitute(paste(rho,'=',corr,', P=',p),list(corr=round(cor$estimate,2),p=round(cor$p.value,5))))+
  xlab(expression('Original Signifiance ' * Delta * 'MI'))+
  ylab(expression('Replication Signifiance' * Delta * 'MI')) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))
dev.off()




