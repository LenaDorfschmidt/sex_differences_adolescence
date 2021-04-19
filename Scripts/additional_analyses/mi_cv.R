library(vows)       
library(MuMIn)      
library(AICcmodavg) 
library(ggplot2) 
library(tidyr)
library(lattice)
library(RColorBrewer)
library(gplots)
library(Rfast)

source(paste0(script.path,'helper/ztest_unpooled.R'))
source(paste0(script.path,'helper/predict.LME.SE.R'))
source(paste0(script.path,'helper/brainplot.R'))


data = 'CV' 

base.path = '~/Documents/Education/Cambridge/fMRI_NSPN/'
script.path = paste0(base.path,'Scripts/')

source(paste0(script.path,'helper/brainplot.R'))

data.path = paste0(base.path,'Data/connectivity_matrices/')
results.path = paste0(base.path,'Results/')
plot.out = paste0(results.path,'replication/CV/')

load(paste0(data.path,'nspn.main.RData'))

ifelse(!dir.exists(file.path(paste0(plot.out))), dir.create(file.path(plot.out)), FALSE)


MRI_table$sex = as.factor(MRI_table$sex)
MRI_table$mri_center = as.factor(MRI_table$mri_center)

str = apply(fc,c(3,1),function(x) mean(x,na.rm=TRUE))
str.cort = apply(fc[,17:346,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc = apply(fc[,1:16,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc.cort = lapply(1:8, function(i) apply(fc[c(i,i+8),,],c(2,3),function(x) mean(x,na.rm=TRUE)))
fc.m <- rowMeans(str)
nroi=346

# Global
fc.glob=cbind(MRI_table,rowMeans(str))
colnames(fc.glob)[ncol(fc.glob)] = 'FC'
glob.m = lme(FC~mri_pseudo_age+sex+mri_center+vol+vol*sex, random=~1|id_nspn,data=fc.glob)

newdat = expand.grid(mri_center=as.factor(3),sex=unique(MRI_table$sex),mri_pseudo_age=range(MRI_table$mri_pseudo_age), vol = range(MRI_table$vol))  
pred= predictSE.lme(glob.m, newdata=newdat, se.fit=T, level=0)   

pdf(paste0(plot.out,'spaghetti.pdf'),height = 3,width = 4)
ggplot(fc.glob,aes(x=mri_pseudo_age,y=FC,color=as.factor(sex),fill=as.factor(sex))) + 
  geom_point(aes(alpha=0.5)) + 
  geom_line(aes(group=id_nspn,alpha=0.5)) +
  geom_ribbon(data=newdat[newdat$sex==0,], aes(x=mri_pseudo_age, ymin=pred$fit[newdat$sex==0]-2*pred$se.fit[newdat$sex==0], ymax=pred$fit[newdat$sex==0]+2*pred$se.fit[newdat$sex==0]), alpha=0.5,fill="grey", inherit.aes=F)  + 
  geom_ribbon(data=newdat[newdat$sex==1,], aes(x=mri_pseudo_age, ymin=pred$fit[newdat$sex==1]-2*pred$se.fit[newdat$sex==1], ymax=pred$fit[newdat$sex==1]+2*pred$se.fit[newdat$sex==1]), alpha=0.5,fill="grey", inherit.aes=F)  + 
  geom_line(data=newdat, aes(y=pred$fit),size=1)+
  xlab('Age') +
  ylab('Global FC')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12))
dev.off()


# ROI-wise
roi.sex.p = roi.sex.t = roi.age.p = roi.age.t = roi.vol.p = roi.vol.t = roi.vol.age.p = 
  roi.vol.age.t = roi.vol.sex.p = roi.vol.sex.t = vector(length=nroi)
for (roi in 1:nroi) {
  df = data.frame(fc=str[,roi], sex = MRI_table$sex, age=MRI_table$mri_pseudo_age, 
                  vol = MRI_table$vol, center = MRI_table$mri_center, id = MRI_table$id_nspn)
  l1 = lme(fc~age+sex+center+vol+vol*sex+age*vol, random=~1|id,data=df)
  
  roi.age.p[roi] = summary(l1)$tTable[2,5]
  roi.age.t[roi] = summary(l1)$tTable[2,4]
  roi.sex.p[roi] = summary(l1)$tTable[3,5]
  roi.sex.t[roi] = summary(l1)$tTable[3,4]
  roi.vol.p[roi] = summary(l1)$tTable[6,5]
  roi.vol.t[roi] = summary(l1)$tTable[6,4]
  roi.vol.sex.p[roi] = summary(l1)$tTable[7,5]
  roi.vol.sex.t[roi] = summary(l1)$tTable[7,4]
  roi.vol.age.p[roi] = summary(l1)$tTable[8,5]
  roi.vol.age.t[roi] = summary(l1)$tTable[8,4]
}

brainplot(roi.vol.age.t,nm,paste0(plot.out,'/roi.vol.age.t'),-max(abs(range(roi.vol.age.t))),max(abs(range(roi.vol.age.t))),'RdBu','stacked', FALSE)
brainplot(roi.vol.sex.t,nm,paste0(plot.out,'/roi.vol.sex.t'),-max(abs(range(roi.vol.sex.t))),max(abs(range(roi.vol.sex.t))),'RdBu','stacked', FALSE)
brainplot(roi.vol.t,nm,paste0(plot.out,'/roi.vol.t'),-max(abs(range(roi.vol.t))),max(abs(range(roi.vol.t))),'RdBu','stacked', FALSE)
brainplot(roi.sex.t,nm,paste0(plot.out,'/roi.sex.t'),-max(abs(range(roi.sex.t))),max(abs(range(roi.sex.t))),'RdBu','stacked', FALSE)
brainplot(roi.age.t,nm,paste0(plot.out,'/roi.age.t'),-max(abs(range(roi.age.t))),max(abs(range(roi.age.t))),'RdBu','stacked', FALSE)


# Edge-wise
edge.sex.p = edge.sex.t = edge.age.p = edge.age.t = edge.vol.p = edge.vol.t = 
  edge.vol.age.p = edge.vol.age.t = edge.vol.sex.p = edge.vol.sex.t = matrix(0, ncol=nroi,nrow=nroi)
for (r in 1:(nroi-1)) {
  for (c in (r+1):nroi) {
    df = data.frame(fc=fc[r,c,], sex = MRI_table$sex, age=MRI_table$mri_pseudo_age, 
                    vol = MRI_table$vol, center = MRI_table$mri_center, id = MRI_table$id_nspn)
    l1 = lme(fc~age+sex+center+vol+vol*sex+age*vol, random=~1|id,data=df)
    
    edge.age.p[r,c] = summary(l1)$tTable[2,4]
    edge.age.t[r,c] = summary(l1)$tTable[2,4]
    edge.sex.p[r,c] = summary(l1)$tTable[3,4]
    edge.sex.t[r,c] = summary(l1)$tTable[3,4]
    edge.vol.p[r,c] = summary(l1)$tTable[6,5]
    edge.vol.t[r,c] = summary(l1)$tTable[6,4]
    edge.vol.sex.p[r,c] = summary(l1)$tTable[7,5]
    edge.vol.sex.t[r,c] = summary(l1)$tTable[7,4]
    edge.vol.age.p[r,c] = summary(l1)$tTable[8,5]
    edge.vol.age.t[r,c] = summary(l1)$tTable[8,4]
  }
}

edge.sex.t = edge.sex.t + t(edge.sex.t)
edge.sex.p = edge.sex.p + t(edge.sex.p)
edge.age.p = edge.age.p + t(edge.age.p)
edge.age.t = edge.age.t + t(edge.age.t)
edge.vol.p = edge.vol.p + t(edge.vol.p)
edge.vol.t = edge.vol.t + t(edge.vol.t)
edge.vol.sex.p = edge.vol.sex.p + t(edge.vol.sex.p)
edge.vol.sex.t = edge.vol.sex.t + t(edge.vol.sex.t)
edge.vol.age.p = edge.vol.age.p + t(edge.vol.age.p)
edge.vol.age.t = edge.vol.age.t + t(edge.vol.age.t)


write.table(edge.t,paste0(plot.out,'edge.t.txt'), row.names = F, col.names=F)
write.table(edge.p,paste0(plot.out,'edge.p.txt'), row.names = F, col.names=F)
write.table(edge.age.p,paste0(plot.out,'edge.age.p.txt'), row.names = F, col.names=F)
write.table(edge.age.t,paste0(plot.out,'edge.age.t.txt'), row.names = F, col.names=F)
write.table(edge.vol.p,paste0(plot.out,'edge.vol.p.txt'), row.names = F, col.names=F)
write.table(edge.vol.t,paste0(plot.out,'edge.vol.t.txt'), row.names = F, col.names=F)
write.table(edge.vol.int.p,paste0(plot.out,'edge.vol.int.p.txt'), row.names = F, col.names=F)
write.table(edge.vol.int.t,paste0(plot.out,'edge.vol.int.t.txt'), row.names = F, col.names=F)



###############################################################
#####################  Maturational Index #####################
###############################################################

if(data=='GSR'){
  low = -0.15
  high = 0.3
} else{
  low=0.2
  high=0.8
}

# Cortical Volume at age 14 
lm.CV1 = lme(vol~mri_pseudo_age+sex + mri_pseudo_age*sex, random=~1|id_nspn,data=MRI_table)
lm.CV2 = lme(vol~mri_pseudo_age+sex + mri_pseudo_age, random=~1|id_nspn,data=MRI_table)
AIC(lm.CV1) < AIC(lm.CV2)

# Choosing interaction effect model based on AIC
lm.CV = lme(vol~mri_pseudo_age+sex+  mri_pseudo_age+sex + mri_pseudo_age*sex, random=~1|id_nspn,data=MRI_table)
CV.m = lm.CV$coefficients$fixed[1]+14*lm.CV$coefficients$fixed[2]+lm.CV$coefficients$fixed[3] + 14*lm.CV$coefficients$fixed[4]
CV.f = lm.CV$coefficients$fixed[1]+14*lm.CV$coefficients$fixed[2]
CV = list(CV.f, CV.m)
 
g_str = c('female','male')
for (s in (1:2)) {
  # set up folder structure
  gender_str = g_str[s]
  gender_idx = (MRI_table$sex == (s-1))
  ifelse(!dir.exists(file.path(paste0(plot.out, gender_str))), dir.create(file.path(plot.out, gender_str)), FALSE)
  
  # GLOBAL
  df.glob = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], x3 = MRI_table$vol[gender_idx], y = fc.m[gender_idx], id = MRI_table$id_nspn[gender_idx])
  l_glob = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df.glob) # fit lme
  glob.str.l.sl.p = summary(l_glob)$tTable[2,5]    # p-value for effect of age
  glob.str.l.sl.t = summary(l_glob)$tTable[2,4] 
  glob.str.l.sl = l_glob$coefficients$fixed[2]
  glob.str.l.14 = l_glob$coefficients$fixed[1] + l_glob$coefficients$fixed[2]*14 + l_glob$coefficients$fixed[3]*(1/3) + l_glob$coefficients$fixed[4]*(1/3) + CV[[s]] * l_glob$coefficients$fixed[5]
  
  write.csv(glob.str.l.14, file = paste0(plot.out,gender_str,'/glob_bl14.txt')) # BL 14 
  write.csv(glob.str.l.sl, file = paste0(plot.out,gender_str,'/glob_delta_age.txt')) # Delta Age
  
  # LOCAL
  # overall, cortex-all, subcortex-all
  nsub=8
  
  str.l.sl = str.l.14 = str.l.sl = str.l.sl.p = str.l.rsq = str.l.sl.t = vector(length = nroi) # overall node strength
  str.cort.l.sl = str.cort.l.14 = str.cort.l.sl.p = str.cort.l.rsq = str.cort.l.sl.t = vector(length = nroi) # cortical node strength
  str.subc.l.sl = str.subc.l.14 = str.subc.l.sl.p = str.subc.l.rsq = str.subc.l.sl.t = vector(length = nroi) # cortical node strength
  str.subc.cort.l.sl = str.subc.cort.l.14 = str.subc.cort.l.sl.p = str.subc.cort.l.sl.t = matrix(nrow=nsub, ncol=nroi)
  
  # LOCAL
  for (n in 1:nroi) {
    # str
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], x3 = MRI_table$vol[gender_idx], y = str[gender_idx,n], id = MRI_table$id_nspn[gender_idx]) # data frame for lme
      l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df) # fit lme
      str.l.sl.p[n] = summary(l)$tTable[2,5]    # p-value for effect of age
      str.l.sl.t[n] = summary(l)$tTable[2,4]
      str.l.sl[n] = l$coefficients$fixed[2]
      str.l.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3) + CV[[s]] * l_glob$coefficients$fixed[5]
    })
    
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], x3 = MRI_table$vol[gender_idx], y = str.cort[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
      str.cort.l.sl.p[n] = summary(l)$tTable[2,5]
      str.cort.l.sl.t[n] = summary(l)$tTable[2,4]
      str.cort.l.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3) +  CV[[s]] * l_glob$coefficients$fixed[5]
      str.cort.l.sl[n] = l$coefficients$fixed[2]
    })
    
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx],  x3 = MRI_table$vol[gender_idx], y = str.subc[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
      str.subc.l.sl.p[n] = summary(l)$tTable[2,5]
      str.subc.l.sl.t[n] = summary(l)$tTable[2,4]
      str.subc.l.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3) +  CV[[s]] * l_glob$coefficients$fixed[5]
      str.subc.l.sl[n] = l$coefficients$fixed[2]
    })
  }
  for (i in 1:nsub) {
    for (n in 1:nroi) {
      try({
        df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], x3 = MRI_table$vol[gender_idx], y = str.subc.cort[[i]][n,gender_idx], id = MRI_table$id_nspn[gender_idx])
        l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)
        str.subc.cort.l.sl.p[i,n] = summary(l)$tTable[2,5]
        str.subc.cort.l.sl.t[i,n] = summary(l)$tTable[2,4]
        str.subc.cort.l.14[i,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3) +  CV[[s]] * l_glob$coefficients$fixed[5]
        str.subc.cort.l.sl[i,n] = l$coefficients$fixed[2]
      })
    }
  }
  
  # correct p-values for multiple comparisons
  str.l.sl.p.fdr = p.adjust(str.l.sl.p, method = 'fdr')           # overall
  str.cort.l.sl.p.fdr = p.adjust(str.cort.l.sl.p, method = 'fdr') # cortex
  str.subc.l.sl.p.fdr = p.adjust(str.subc.l.sl.p, method = 'fdr') # subcortex
  
  # BL 14
  write.csv(str.cort.l.14, file = paste0(plot.out,gender_str,'/cort_bl14.txt')) # Cortex
  write.csv(str.subc.l.14, file = paste0(plot.out,gender_str,'/subc_bl14.txt')) # Subcortex
  write.csv(str.l.14, file = paste0(plot.out,gender_str,'/all_bl14.txt')) # Overall
  write.csv(str.subc.cort.l.14,paste0(plot.out,gender_str,'/sub_cort_bl14.csv'))
  
  # Delta Age
  write.csv(str.cort.l.sl, file = paste0(plot.out,gender_str,'/cort_delta_age.txt')) # Cortex
  write.csv(str.subc.l.sl, file = paste0(plot.out,gender_str,'/subc_delta_age.txt')) # Subcortex
  write.csv(str.l.sl, file = paste0(plot.out,gender_str,'/all_delta_age.txt')) # Overall
  write.csv(str.subc.cort.l.sl,paste0(plot.out,gender_str,'/sub_cort_delta_age.csv'))
  
  brainplot(str.cort.l.14[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_bl14'),low,high,'plasma','stacked', FALSE)
  brainplot(str.subc.l.14[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_bl14'),low,high,'plasma','stacked', FALSE)
  
  brainplot(str.cort.l.sl[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_delta_age'),-0.015,0.015,'seismic','stacked')
  brainplot(str.subc.l.sl[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/subc_delta_age'),-0.015,0.015,'seismic','stacked')
  
  for (i in 1:8){
    brainplot((as.vector(t(str.subc.cort.l.14[i,]))[17:346]),nm[17:346],paste0(plot.out,gender_str,'/subc.cort.bl14',nm[i]),low,high,'plasma','stacked',FALSE)
    brainplot((as.vector(t(str.subc.cort.l.sl[i,]))[17:346]),nm[17:346],paste0(plot.out,gender_str,'/subc.cort.sl',nm[i]),-0.015,0.015,'seismic','stacked',FALSE)
  }
  
  ## MATURATIONAL INDEX
  ##### Fit lme at level of individual FC edges
  err.idx = array(FALSE,dim=c(nroi,nroi))
  fc.l.sl.t = fc.l.sl.p = fc.l.14 = fc.l.sl = fc.l.interaction.p = fc.l.sl.t = fc.l.rsq = array(0,dim=c(nroi,nroi))
  for (i in 1:(nroi-1)) { # only loop over upper triangular (as matrices are symmetric), and then fill lower triangular with the transpose
    print(i) # track progress
    for (j in (i+1):nroi) {
      possibleError <- tryCatch({
        df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], x3 = MRI_table$vol[gender_idx], y = fc[i,j,gender_idx], id = MRI_table$id_nspn[gender_idx])  # data frame for model fitting
        l = lme(y ~ x1 + x2 + x3, random = ~ 1|id, data = df)   # fit model
        fc.l.sl.t[i,j] = summary(l)$tTable[2,4]            # t-statistic of effect of age
        fc.l.sl.p[i,j] = summary(l)$tTable[2,5]
        # predicted value at 14 (includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
        fc.l.14[i,j] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3) +  CV[[s]] * l_glob$coefficients$fixed[5] 
        fc.l.sl[i,j] = l$coefficients$fixed[2]   # beta coefficient of effect of age
        
        # fc.l.rsq[i,j] = r.squaredGLMM(l)[1,1]     
      }
      ,
      error=function(e) {
        print(paste("Oops! --> Error in Row ",i," Column ",j,sep = ""))
        err.idx[i,j] = TRUE
      }
      )
      if(inherits(possibleError, "error")) next
    }
  }
  # fill empty lower triangulars of matrices of lme parameters
  fc.l.sl.t = fc.l.sl.t + t(fc.l.sl.t)
  fc.l.sl.p = fc.l.sl.p + t(fc.l.sl.p)
  fc.l.14 = fc.l.14 + t(fc.l.14)
  fc.l.sl = fc.l.sl + t(fc.l.sl)
  fc.l.rsq = fc.l.rsq + t(fc.l.rsq)
  
  write.csv(fc.l.14, file = paste0(plot.out,gender_str,'/fc14.txt'))
  write.csv(fc.l.sl, file = paste0(plot.out,gender_str,'/delta_age_beta.txt'))
  write.csv(fc.l.sl.t, file = paste0(plot.out,gender_str,'/delta_age_t.txt'))
  
  ##### Relationship between edge-wise FC at 14 and change as a function of age at each node
  # regional relationships (region itself is removed)
  rr.slt_14.rho = rr.slt_14.p = vector(length=nroi)
  cc.slt_14.rho = cc.slt_14.p = vector(length=nroi)
  
  for (i in 1:nroi) {
    rr.slt_14.rho[i] = cor.test(fc.l.14[i,-i],fc.l.sl.t[i,-i],method='spearman')$estimate   # spearman rho
    rr.slt_14.p[i] = cor.test(fc.l.14[i,-i],fc.l.sl.t[i,-i],method='spearman')$p.value      # p-value
    
    cc.slt_14.rho[i] = cor.test(fc.l.14[i,-i],fc.l.sl[i,-i],method='spearman')$estimate   # spearman rho
    cc.slt_14.p[i] = cor.test(fc.l.14[i,-i],fc.l.sl[i,-i],method='spearman')$p.value      # p-value
  }
  
  # correct p-values for multiple comparisons
  rr.slt_14.p.fdr = p.adjust(rr.slt_14.p, method = 'fdr')
  cc.slt_14.p.fdr = p.adjust(cc.slt_14.p, method = 'fdr')
  
  write.csv(rr.slt_14.rho, file = paste0(plot.out,gender_str,'/rr_mat_idx.txt'))
  write.csv(cc.slt_14.rho, file = paste0(plot.out,gender_str,'/cc_mat_idx.txt'))
  
}
mi.m = read.csv(paste0(plot.out,'/male/cc_mat_idx.txt'))[,2]
mi.f = read.csv(paste0(plot.out,'/female/cc_mat_idx.txt'))[,2]
write.table(mi.f-mi.m,paste0(plot.out,'/diff.mi.txt'), row.names = F, col.names = F)

m.str.subc.cort.l.14 = read.csv(paste0(plot.out,'male/sub_cort_bl14.csv'))
f.str.subc.cort.l.14 = read.csv(paste0(plot.out,'female/sub_cort_bl14.csv'))
m.str.subc.cort.l.sl = read.csv(paste0(plot.out,'male/sub_cort_delta_age.csv'))
f.str.subc.cort.l.sl = read.csv(paste0(plot.out,'female/sub_cort_delta_age.csv'))

diff.str.subc.cort.l.14 = f.str.subc.cort.l.14-m.str.subc.cort.l.14
diff.str.subc.cort.l.sl = f.str.subc.cort.l.sl-m.str.subc.cort.l.sl

for (i in 1:8){
  brainplot((as.vector(t(diff.str.subc.cort.l.14[i,2:347]))[17:346]),nm[17:346],paste0(plot.out,'/diff.subc.cort.bl14',nm[i]),-0.15,0.15,'blue2red_harsh','stacked',FALSE)
  brainplot((as.vector(t(diff.str.subc.cort.l.sl[i,2:347]))[17:346]),nm[17:346],paste0(plot.out,'/diff.subc.cort.sl',nm[i]),-0.025,0.025,'blue2red_harsh','stacked',FALSE)
}

bl14.roi.m = read.csv(paste0(plot.out,'/male/all_bl14.txt'))[,2]
delta.roi.m = read.csv(paste0(plot.out,'/male/all_delta_age.txt'))[,2]
bl14.roi.f = read.csv(paste0(plot.out,'/female/all_bl14.txt'))[,2]
delta.roi.f = read.csv(paste0(plot.out,'/female/all_delta_age.txt'))[,2]

subc.bl14.roi.m = read.csv(paste0(plot.out,'/male/subc_bl14.txt'))[,2]
subc.bl14.roi.f = read.csv(paste0(plot.out,'/female/subc_bl14.txt'))[,2]
cort.bl14.roi.m = read.csv(paste0(plot.out,'/male/cort_bl14.txt'))[,2]
cort.bl14.roi.f = read.csv(paste0(plot.out,'/female/cort_bl14.txt'))[,2]

brainplot(cort.bl14.roi.f[17:nroi]-cort.bl14.roi.m[17:nroi],nm[17:nroi],paste0(plot.out,'/diff_cort_bl14'),-0.15,0.15,'blue2red_harsh','stacked')
brainplot(subc.bl14.roi.f[17:nroi]-subc.bl14.roi.m[17:nroi],nm[17:nroi],paste0(plot.out,'/diff_subc_bl14'),-0.15,0.15,'blue2red_harsh','stacked')

subc.sl.roi.m = read.csv(paste0(plot.out,'/male/subc_delta_age.txt'))[,2]
subc.sl.roi.f = read.csv(paste0(plot.out,'/female/subc_delta_age.txt'))[,2]
cort.sl.roi.m = read.csv(paste0(plot.out,'/male/cort_delta_age.txt'))[,2]
cort.sl.roi.f = read.csv(paste0(plot.out,'/female/cort_delta_age.txt'))[,2]

brainplot(cort.sl.roi.f[17:nroi]-cort.sl.roi.m[17:nroi],nm[17:nroi],paste0(plot.out,'/diff_cort_delta_age'),-0.025,0.025,'blue2red_harsh','stacked')
brainplot(subc.sl.roi.f[17:nroi]-subc.sl.roi.m[17:nroi],nm[17:nroi],paste0(plot.out,'/diff_subc_delta_age'),-0.025,0.025,'blue2red_harsh','stacked')

#======================================================================================================
#                                Parametric  Sex Difference Significance                                      
#======================================================================================================
m.14 = as.matrix(read.csv(paste0(plot.out,'/male/fc14.txt'))[,2:347])
colnames(m.14) <- NULL
f.14 = as.matrix(read.csv(paste0(plot.out,'/female/fc14.txt'))[,2:347])
colnames(f.14) <- NULL
m.delta = as.matrix(read.csv(paste0(plot.out,'/male/delta_age_beta.txt'))[,2:347])
colnames(m.delta) <- NULL
f.delta = as.matrix(read.csv(paste0(plot.out,'/female/delta_age_beta.txt'))[,2:347])
colnames(f.delta) <- NULL
nroi = 346

m.mi.beta = f.mi.beta = m.mi.se = f.mi.se = vector(length = nroi)
# Comparing Coefficients
for (i in 1:nroi) {
  m.lm = lm(scale(m.delta[i,-i],center=TRUE, scale=TRUE)~scale(m.14[i,-i],center=TRUE, scale=TRUE))
  m.mi.beta[i] = summary(m.lm)$coefficients[2,1]
  m.mi.se[i] = summary(m.lm)$coefficients[2,2]
  f.lm = lm(scale(f.delta[i,-i],center=TRUE, scale=TRUE)~scale(f.14[i,-i],center=TRUE, scale=TRUE))
  f.mi.beta[i] = summary(f.lm)$coefficients[2,1]
  f.mi.se[i] = summary(f.lm)$coefficients[2,2]
}

z = vector(length=nroi)
for (i in (1:346)) {
  z[i] = (m.mi.beta[i]-f.mi.beta[i])/sqrt(m.mi.se[i]^2+f.mi.se[i]^2)
}

brainplot(z,nm,paste0(plot.out, 'z'),-10,10,'viridis','stacked')

pvalue2sided=2*pnorm(-abs(z))
write.table(z, file = paste0(plot.out,'/z.MI.txt'),row.names = FALSE, col.names = FALSE)
write.table(pvalue2sided, file = paste0(plot.out,'/z.p.MI.txt'),row.names = FALSE, col.names = FALSE)

fdr.pvalue2sided = p.adjust(pvalue2sided,method = 'fdr')
write.table(fdr.pvalue2sided, file = paste0(plot.out,'/z.p.fdr.MI.txt'),row.names = FALSE, col.names = FALSE)

print(paste('Number of significantly different ROIs: ',sum(fdr.pvalue2sided<(0.05))))

z[fdr.pvalue2sided>=0.05]=NA
brainplot(z,nm,paste0(plot.out, 'z.thresh.'),-10,10,'viridis','stacked')
write.table(z, file = paste0(plot.out,'/z.MI.thresholded.txt'),row.names = FALSE, col.names = FALSE)

brainplot(diff.mi,nm,paste0(plot.out, 'diff.mi'),-0.8,0.8,'seismic','stacked')

diff.mi.thresh = diff.mi
diff.mi.thresh[fdr.pvalue2sided>=0.05] = NA
brainplot(diff.mi.thresh,nm,paste0(plot.out, 'diff.mi.thresh'),-0.8,0.8,'seismic','stacked')



