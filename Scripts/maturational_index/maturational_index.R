library(vows)       
library(MuMIn)      
library(AICcmodavg) 
library(ggplot2) 
library(tidyr)
library(lattice)
library(RColorBrewer)
library(gplots)
library(Rfast)

data = 'NSPN' 
script.path = paste0('Scripts/')

source(paste0(script.path,'helper/ztest_unpooled.R'))
source(paste0(script.path,'helper/predict.LME.SE.R'))
source(paste0(script.path,'helper/brainplot.R'))

data.path = paste0('Data/connectivity_matrices/')
results.path = paste0('Results/')
# load in data: NSPN, NSPN low motion or GSR
if (data == 'NSPN') {
  load(paste0(data.path,'nspn.main.RData'))
  plot.out = paste0(results.path,'maturational_index/')
}else if(data == 'motion.matched'){
  load(paste0(data.path, 'nspn.motion.matched.RData'))
  plot.out = paste0(results.path,'replication/motionmatched/')
}else if(data == 'nspn.repl.FD.by.sex'){
  load(paste0(data.path, 'nspn.repl.FD.by.sex.RData'))
  plot.out = paste0(results.path,'replication/regr_FD_by_sex/')
}else if(data == 'nspn.gsr'){
  load(paste0(data.path, 'nspn.gsr.RData'))
  plot.out = paste0(results.path,'replication/nspn.gsr/')
}

ifelse(!dir.exists(file.path(paste0(results.path,'replication'))), dir.create(file.path(results.path,'replication')), FALSE)
ifelse(!dir.exists(file.path(paste0(results.path))), dir.create(file.path(results.path)), FALSE)
ifelse(!dir.exists(file.path(paste0(plot.out))), dir.create(file.path(plot.out)), FALSE)

MRI_table$sex = as.factor(MRI_table$sex)
MRI_table$mri_center = as.factor(MRI_table$mri_center)

str = apply(fc,c(3,1),function(x) mean(x,na.rm=TRUE))
str.cort = apply(fc[,17:346,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc = apply(fc[,1:16,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc.cort = lapply(1:8, function(i) apply(fc[c(i,i+8),,],c(2,3),function(x) mean(x,na.rm=TRUE)))
fc.m <- rowMeans(str)
nroi=346

# Global FC development
# Figure 1A
fc.glob=cbind(MRI_table,rowMeans(str))
colnames(fc.glob)[ncol(fc.glob)] = 'fc.glob'
glob.m = lme(fc.glob~mri_pseudo_age+sex+mri_center, random=~1|id_nspn,data=fc.glob)

fc.glob=cbind(MRI_table,rowMeans(str))
colnames(fc.glob)[dim(fc.glob)[2]] = 'FC'
pdf(paste0(plot.out,'spaghetti.pdf'),height = 3,width = 4)
ggplot(fc.glob,aes(x=mri_pseudo_age,y=FC,color=as.factor(sex),fill=as.factor(sex))) + 
  geom_point(aes(alpha=0.5)) + 
  geom_line(aes(group=id_nspn,alpha=0.5)) +
  geom_smooth(method='lm', alpha = 0.7) + 
  xlab('Age') +
  ylab('Global FC')+
  theme_bw()+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12))
dev.off()

###############################################################
#####################  Maturational Index #####################
###############################################################

if(data=='nspn.gsr'){
  low = -0.1
  high = 0.1
} else{
  low=0.2
  high=0.8
}


# Allocate space for arrays
str.l.sl = str.l.14 = str.l.sl = str.l.sl.p = str.l.sl.t = str.l.sl.p.fdr = array(dim = c(2,nroi)) # overall 
str.cort.l.sl = str.cort.l.14 = str.cort.l.sl.p = str.cort.l.sl.t = str.cort.l.sl.p.fdr = array(dim = c(2,nroi)) # cortico-cortical 
str.subc.l.sl = str.subc.l.14 = str.subc.l.sl.p = str.subc.l.sl.t = str.subc.l.sl.p.fdr = array(dim = c(2,nroi)) # subcortico-cortcial 
str.subc.cort.l.sl = str.subc.cort.l.14 = str.subc.cort.l.sl.p = str.subc.cort.l.sl.t = array(dim = c(2,nsub,nroi))
str.l.sl.sign = str.l.14.sign = str.cort.l.sl.sign = str.cort.l.14.sign = str.subc.l.sl.sign = str.subc.l.14.sign = array(dim = c(2,nroi,2))
str.subc.cort.l.sl.sign = str.subc.cort.l.14.sign = array(NA, dim= c(2,nsub,nroi,2))

coefvec = c(1,14,1/3,1/3)
g_str = c('female','male')
for (s in (1:2)) {
  # set up folder structure
  gender_str = g_str[s]
  gender_idx = (MRI_table$sex == (s-1))
  ifelse(!dir.exists(file.path(paste0(plot.out, gender_str))), dir.create(file.path(plot.out, gender_str)), FALSE)
  
  # REGION-WISE
  # overall, cortex-all, subcortex-all
  nsub=8
  for (n in 1:nroi) {
    # str
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str[gender_idx,n], id = MRI_table$id_nspn[gender_idx]) # data frame for lme
      l = lme(y ~ x1 + x2, random = ~ 1|id, data = df) # fit lme
      str.l.sl.p[s,n] = summary(l)$tTable[2,5]    # p-value for effect of age
      str.l.sl.t[s,n] = summary(l)$tTable[2,4]
      str.l.sl[s,n] = l$coefficients$fixed[2]
      str.l.14[s,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      str.l.14.sign[s,n,] = predict.LME.SE(l,coefvec)
      str.l.sl.sign[s,n,] = summary(l)$tTable[2,1:2]
    })
    
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str.cort[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
      str.cort.l.sl.p[s,n] = summary(l)$tTable[2,5]
      str.cort.l.sl.t[s,n] = summary(l)$tTable[2,4]
      str.cort.l.14[s,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      str.cort.l.sl[s,n] = l$coefficients$fixed[2]
      str.cort.l.14.sign[s,n,] = predict.LME.SE(l,coefvec)
      str.cort.l.sl.sign[s,n,] = summary(l)$tTable[2,1:2]
    })
    
    try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str.subc[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
      str.subc.l.sl.p[s,n] = summary(l)$tTable[2,5]
      str.subc.l.sl.t[s,n] = summary(l)$tTable[2,4]
      str.subc.l.14[s,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      str.subc.l.sl[s,n] = l$coefficients$fixed[2]
      str.subc.l.14.sign[s,n,] = predict.LME.SE(l,coefvec)
      str.subc.l.sl.sign[s,n,] = summary(l)$tTable[2,1:2]
    })
  }
  for (i in 1:nsub) {
    for (n in 1:nroi) {
      try({
        df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str.subc.cort[[i]][n,gender_idx], id = MRI_table$id_nspn[gender_idx])
        l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
        str.subc.cort.l.sl.p[s,i,n] = summary(l)$tTable[2,5]
        str.subc.cort.l.sl.t[s,i,n] = summary(l)$tTable[2,4]
        str.subc.cort.l.14[s,i,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
        str.subc.cort.l.sl[s,i,n] = l$coefficients$fixed[2]
        str.subc.cort.l.14.sign[s,i,n,] = predict.LME.SE(l,coefvec)
        str.subc.cort.l.sl.sign[s,i,n,] = summary(l)$tTable[2,1:2]
      })
    }
  }
 
  # BL 14
  write.csv(str.cort.l.14[s,], file = paste0(plot.out,gender_str,'/cort_bl14.txt')) # Cortex
  write.csv(str.subc.l.14[s,], file = paste0(plot.out,gender_str,'/subc_bl14.txt')) # Subcortex
  write.csv(str.l.14[s,], file = paste0(plot.out,gender_str,'/all_bl14.txt')) # Overall
  write.csv(str.subc.cort.l.14[s,,],paste0(plot.out,gender_str,'/sub_cort_bl14.csv'))
  
  # Delta Age
  write.csv(str.cort.l.sl[s,], file = paste0(plot.out,gender_str,'/cort_delta_age.txt')) # Cortex
  write.csv(str.subc.l.sl[s,], file = paste0(plot.out,gender_str,'/subc_delta_age.txt')) # Subcortex
  write.csv(str.l.sl[s,], file = paste0(plot.out,gender_str,'/all_delta_age.txt')) # Overall
  write.csv(str.subc.cort.l.sl[s,,],paste0(plot.out,gender_str,'/sub_cort_delta_age.csv'))
  
  ## Figure 1B
  brainplot(str.l.14,nm,paste0(plot.out,gender_str,'/all_bl14'),low,high,'plasma','stacked', FALSE)
  brainplot(str.l.sl,nm,paste0(plot.out,gender_str,'/all_delta_age'),-0.015,0.015,'seismic','stacked')
  
  ## SI Figure S8 and S10
  brainplot(str.cort.l.14[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_bl14'),low,high,'plasma','stacked', FALSE)
  brainplot(str.subc.l.14[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_bl14'),low,high,'plasma','stacked', FALSE)
  brainplot(str.cort.l.sl[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/cort_delta_age'),-0.015,0.015,'seismic','stacked')
  brainplot(str.subc.l.sl[17:nroi],nm[17:nroi],paste0(plot.out,gender_str,'/subc_delta_age'),-0.015,0.015,'seismic','stacked')
  for (i in 1:8){
    brainplot((as.vector(t(str.subc.cort.l.14[s,i,]))[17:346]),nm[17:346],paste0(plot.out,gender_str,'/subc.cort.bl14',nm[i]),low,high,'plasma','stacked',FALSE)
    brainplot((as.vector(t(str.subc.cort.l.sl[s,i,]))[17:346]),nm[17:346],paste0(plot.out,gender_str,'/subc.cort.sl',nm[i]),-0.015,0.015,'seismic','stacked',FALSE)
  }
  
  ## MATURATIONAL INDEX
  ##### Fit lme at level of individual FC edges
  err.idx = array(FALSE,dim=c(nroi,nroi))
  fc.l.sl.t = fc.l.sl.p = fc.l.14 = fc.l.sl = fc.l.interaction.p = fc.l.sl.t = fc.l.rsq = array(0,dim=c(nroi,nroi))
  for (i in 1:(nroi-1)) { # only loop over upper triangular (as matrices are symmetric), and then fill lower triangular with the transpose
    print(i) # track progress
    for (j in (i+1):nroi) {
      possibleError <- tryCatch({
        df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = fc[i,j,gender_idx], id = MRI_table$id_nspn[gender_idx])  # data frame for model fitting
        l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)   # fit model
        fc.l.sl.t[i,j] = summary(l)$tTable[2,4]            # t-statistic of effect of age
        fc.l.sl.p[i,j] = summary(l)$tTable[2,5]
        # predicted value at 14 (includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
        fc.l.14[i,j] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)  # y = mx+b (+ average effect of gender)
        fc.l.sl[i,j] = l$coefficients$fixed[2]   # beta coefficient of effect of age
        fc.l.rsq[i,j] = r.squaredGLMM(l)[1,1]     
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
  
  # Figure 2B
  brainplot(cc.slt_14.rho,nm,paste0(plot.out,gender_str,'/mi'),-0.8,0.8,'RdBu','stacked',sig=NULL, disp_nan =FALSE)
  
  # correct p-values for multiple comparisons
  rr.slt_14.p.fdr = p.adjust(rr.slt_14.p, method = 'fdr')
  cc.slt_14.p.fdr = p.adjust(cc.slt_14.p, method = 'fdr')
  write.csv(rr.slt_14.rho, file = paste0(plot.out,gender_str,'/rr_mat_idx.txt'))
  write.csv(cc.slt_14.rho, file = paste0(plot.out,gender_str,'/cc_mat_idx.txt'))
}

N.females = sum(MRI_table$sex==0); N.males = sum(MRI_table$sex==1)

diff.bl14 = ztest_unpooled(str.l.14.sign[1,,1],str.l.14.sign[2,,1],str.l.14.sign[1,,2],str.l.14.sign[2,,2],N.females,N.males)
diff.cort.bl14 = ztest_unpooled(str.cort.l.14.sign[1,,1],str.cort.l.14.sign[2,,1],str.cort.l.14.sign[1,,2],str.cort.l.14.sign[2,,2],N.females,N.males)
diff.subc.bl14 = ztest_unpooled(str.subc.l.14.sign[1,,1],str.subc.l.14.sign[2,,1],str.subc.l.14.sign[1,,2],str.subc.l.14.sign[2,,2],N.females,N.males)
diff.sl = ztest_unpooled(str.l.sl.sign[1,,1],str.l.sl.sign[2,,1],str.l.sl.sign[1,,2],str.l.sl.sign[2,,2],N.females,N.males)
diff.cort.sl = ztest_unpooled(str.cort.l.sl.sign[1,,1],str.cort.l.sl.sign[2,,1],str.cort.l.sl.sign[1,,2],str.cort.l.sl.sign[2,,2],N.females,N.males)
diff.subc.sl = ztest_unpooled(str.subc.l.sl.sign[1,,1],str.subc.l.sl.sign[2,,1],str.subc.l.sl.sign[1,,2],str.subc.l.sl.sign[2,,2],N.females,N.males)

## Figure 1C, 1D and Supplemantary Figures S8, S9, S10, S11
try(brainplot(ifelse(diff.bl14$p.fdr.unpooled[17:nroi] < 0.05,(str.l.14[1,]-str.l.14[2,])[17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.cort.bl14_parametric'),-0.15,0.15,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
try(brainplot(ifelse(diff.sl$p.fdr.unpooled[17:nroi] < 0.05,(str.cort.l.sl[1,]-str.cort.l.sl[2,])[17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.subc.sl_parametric'),-0.025,0.025,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
try(brainplot(ifelse(diff.cort.bl14$p.fdr.unpooled[17:346] < 0.05,(str.cort.l.14[1,]-str.cort.l.14[2,])[17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.cort.bl14_parametric'),-0.15,0.15,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
try(brainplot(ifelse(diff.cort.sl$p.fdr.unpooled[17:346] < 0.05,(str.cort.l.sl[1,]-str.cort.l.sl[2,])[17:nroi],NA),nm[17:346],paste0(plot.out,'/diff.cort.sl_parametric'),-0.025,0.025,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
try(brainplot(ifelse(diff.subc.bl14$p.fdr.unpooled[17:nroi] < 0.05,(str.subc.l.14[1,]-str.subc.l.14[2,])[17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.subc.bl14_parametric'),-0.15,0.15,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
try(brainplot(ifelse(diff.subc.sl$p.fdr.unpooled[17:nroi] < 0.05,(str.subc.l.sl[1,]-str.subc.l.sl[2,])[17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.subc.sl_parametric'),-0.025,0.025,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))

diff.subc.cort.bl14 = diff.subc.cort.sl = list()
for (i in 1:nsub) {
  diff.subc.cort.bl14[[i]] = ztest_unpooled(str.subc.cort.l.14.sign[1,i,,1],str.subc.cort.l.14.sign[2,i,,1],str.subc.cort.l.14.sign[1,i,,2],str.subc.cort.l.14.sign[2,i,,2],N.females,N.males)
  diff.subc.cort.sl[[i]] = ztest_unpooled(str.subc.cort.l.sl.sign[1,i,,1],str.subc.cort.l.sl.sign[2,i,,1],str.subc.cort.l.sl.sign[1,i,,2],str.subc.cort.l.sl.sign[2,i,,2],N.females,N.males)
  try(brainplot(ifelse(diff.subc.cort.bl14[[i]]$p.fdr.unpooled[17:nroi] < 0.05,(str.subc.cort.l.14[1,,]-str.subc.cort.l.14[2,,])[i,17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.subc.cort.bl14_parametric',nm[i]),-0.15,0.15,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
  try(brainplot(ifelse(diff.subc.cort.sl[[i]]$p.fdr.unpooled[17:nroi] < 0.05,(str.subc.cort.l.sl[1,,]-str.subc.cort.l.sl[2,,])[i,17:nroi],NA),nm[17:nroi],paste0(plot.out,'/diff.subc.cort.sl_parametric',nm[i]),-0.025,0.025,'blue2red_harsh','stacked',sig=NULL, disp_nan =FALSE))
}

#======================================================================================================
#                                  Sex Difference in Maturational Index                                      
#======================================================================================================
mi.m = read.csv(paste0(plot.out,'/male/cc_mat_idx.txt'))[,2]
mi.f = read.csv(paste0(plot.out,'/female/cc_mat_idx.txt'))[,2]
write.table(mi.f-mi.m,paste0(plot.out,'/diff.mi.txt'), row.names = F, col.names = F)
diff.mi = mi.f - mi.m

# Figure 2C
brainplot(diff.mi,nm,paste0(plot.out,'/diff.mi'),-0.8,0.8,'seismic','stacked',sig=NULL, disp_nan =FALSE)

m.14 = as.matrix(read.csv(paste0(plot.out,'/male/fc14.txt'))[,2:347])
f.14 = as.matrix(read.csv(paste0(plot.out,'/female/fc14.txt'))[,2:347])
m.delta = as.matrix(read.csv(paste0(plot.out,'/male/delta_age_beta.txt'))[,2:347])
f.delta = as.matrix(read.csv(paste0(plot.out,'/female/delta_age_beta.txt'))[,2:347])
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

pvalue2sided=2*pnorm(-abs(z))
write.table(z, file = paste0(plot.out,'/z.MI.txt'),row.names = FALSE, col.names = FALSE)
fdr.pvalue2sided = p.adjust(pvalue2sided,method = 'fdr')
write.table(fdr.pvalue2sided, file = paste0(plot.out,'/z.p.fdr.MI.txt'),row.names = FALSE, col.names = FALSE)
write.table(ifelse(fdr.pvalue2sided<0.05,z,NA), file = paste0(plot.out,'/z.MI.thresholded.txt'),row.names = FALSE, col.names = FALSE)

print(paste('Number of significantly different ROIs: ',sum(fdr.pvalue2sided<(0.05))))

diff.mi.thresh = ifelse(fdr.pvalue2sided<0.05,diff.mi,NA)

# Supplementary Figure S13A
brainplot(ifelse(fdr.pvalue2sided<0.05,z,NA),nm,paste0(plot.out, 'z.diff.mi'),-10,10,'viridis','stacked')
brainplot(diff.mi.thresh,nm,paste0(plot.out,'/diff.mi.thresholded'),-0.8,0.8,'seismic','stacked',sig=NULL, disp_nan =FALSE)

#======================================================================================================
#                                Disecting Directionality in MI                                     
#======================================================================================================
bl14 = list(m.14,f.14)
delta.age = list(m.delta,f.delta)
growth = matrix(nrow=2, ncol=nroi)
# Sex: 1 -> male, 2 -> female
for (sex in 1:2) {
  tmp.bl14 = bl14[[sex]]
  tmp.delta.age = delta.age[[sex]]
  gender_idx = (MRI_table$sex == (sex-1))
  
  for (roi in 1:nroi) {
    growth[sex,roi] = sum(tmp.delta.age[roi,-roi] > 0)
  }
}

# Supplementary Figure S14
brainplot(ifelse(mi.f<0,(growth[1,]/345),NA),nm,paste0(plot.out,'f.growth.disr'),0,1,'PuOr','stacked')
brainplot(ifelse(mi.m<0,(growth[1,]/345),NA),nm,paste0(plot.out,'m.growth.disr'),0,1,'PuOr','stacked')
brainplot(ifelse(mi.f>0,(growth[1,]/345),NA),nm,paste0(plot.out,'f.growth.cons'),0,1,'PuOr','stacked')
brainplot(ifelse(mi.m>0,(growth[1,]/345),NA),nm,paste0(plot.out,'m.growth.cons'),0,1,'PuOr','stacked')

write.table(m.growth,paste0(plot.out,'growth.disruptive.male.txt'), row.names = FALSE, col.names = FALSE)
write.table(f.growth,paste0(plot.out,'growth.disruptive.female.txt'), row.names = FALSE, col.names = FALSE)
write.table(diff.growth,paste0(plot.out,'growth.disruptive.diff.txt'), row.names = FALSE, col.names = FALSE)

df = data.frame(diff = diff.mi, mi.f=mi.f, mi.m=mi.m)
df$more.disruptive = with(df, (diff<0)&(mi.f<0)); df$less.conservative = with(df, (diff<0)&(mi.f>0))
df$more.conservative = with(df, (diff>0)&(mi.f>0)); df$less.disruptive = with(df, (diff>0)&(mi.f<0))

df$all[df$more.disruptive]=1; df$all[df$less.conservative]=2; df$all[df$less.disruptive]=3; df$all[df$more.conservative]=4

# Supplementary Figure S15
brainplot(df$all,nm,paste0(plot.out,'/trends'),1,4,'RdBu','stacked')
write.table(df$all,paste0(plot.out,'trends.txt'), row.names = FALSE, col.names = FALSE)

df$diff.mi.more.disr[df$more.disruptive] = diff.mi[df$more.disruptive]; df$diff.mi.more.disr[!(fdr.pvalue2sided<0.05)] = NA
df$diff.mi.less.cons[df$less.conservative] = diff.mi[df$less.conservative]; df$diff.mi.less.cons[!(fdr.pvalue2sided<0.05)] = NA
df$diff.mi.less.disr[df$less.disruptive] = diff.mi[df$less.disruptive]; df$diff.mi.less.disr[!(fdr.pvalue2sided<0.05)] = NA
df$diff.mi.more.cons[df$more.conservative] = diff.mi[df$more.conservative]; df$diff.mi.more.cons[!(fdr.pvalue2sided<0.05)] = NA

# Figure 2D
brainplot(df$diff.mi.less.disr,nm,paste0(plot.out,'less.disruptive'),-0.8,0.8,'seismic','stacked')
brainplot(df$diff.mi.more.disr,nm,paste0(plot.out,'more.disruptive'),-0.8,0.8,'seismic','stacked')
brainplot(df$diff.mi.less.cons,nm,paste0(plot.out,'less.conservative'),-0.8,0.8,'seismic','stacked')
brainplot(df$diff.mi.more.cons,nm,paste0(plot.out,'more.conservative'),-0.8,0.8,'seismic','stacked')

print(paste('Frequency of more disruptive changes: ', 
            toString(sum(!is.na(df$diff.mi.more.disr))/(sum(!is.na(df$diff.mi.less.cons))+sum(!is.na(df$diff.mi.more.disr))))))
print(paste('Frequency of more conservative changes: ',
            toString(sum(!is.na(df$diff.mi.more.cons))/(sum(!is.na(df$diff.mi.more.cons))+sum(!is.na(df$diff.mi.less.disr))))))

# examplary more disruptive
roi=which(!is.na(df$diff.mi.more.disr))[6]
roi.name=nm[which(!is.na(df$diff.mi.more.disr))[6]]
df.plot = data.frame(bl14 = c(m.14[roi,-roi],f.14[roi,-roi]), sl=c(m.delta[roi,-roi],f.delta[roi,-roi]),sex=as.factor(c(rep('male',345),rep('female',345))))

# Figure 2D
pdf(paste0(plot.out, 'plot.example.more.disruptive.pdf'),3.6,2)
ggplot(df.plot,aes(x=bl14,y=sl,fill=sex,group=sex,color=sex)) + 
  geom_point(alpha=0.4) +
  theme_bw() +
  xlab(expression('BL'[14]))+
  ylab(expression(paste(Delta,'FC'[14-25])))+
  geom_smooth(method='lm') +
  ggtitle(roi.name) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size=14,hjust=0.5))
dev.off()


# examplary less conservative
roi=which(!is.na(df$diff.mi.less.cons))[19]
df.plot = data.frame(bl14 = c(m.14[roi,-roi],f.14[roi,-roi]), sl=c(m.delta[roi,-roi],f.delta[roi,-roi]),sex=as.factor(c(rep('male',345),rep('female',345))))
roi.name=nm[which(!is.na(df$diff.mi.less.cons))[19]]
# Figure 2D
pdf(paste0(plot.out,'plot.example.less.conservative.pdf'),3.6,2)
ggplot(df.plot,aes(x=bl14,y=sl,fill=sex,group=sex,color=sex)) + 
  geom_point(alpha=0.4) +
  theme_bw() +
  xlab(expression('BL'[14]))+
  ylab(expression(paste(Delta,'FC'[14-25])))+
  geom_smooth(method='lm') +
  ggtitle(roi.name) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        legend.text=element_text(size=12),
        legend.title = element_text(size=12),
        plot.title = element_text(size=14,hjust=0.5))
dev.off()


#======================================================================================================
#                                          Network Decoding                                     
#======================================================================================================

# von Economo network decoding
df.vonEconomo = data.frame(mi = diff.mi,class = vonEconomo.classes)
df.vonEconomo$class = as.factor(df.vonEconomo$class)

classmu = sapply(1:length(unique(vonEconomo)),function(class) mean(diff.mi[vonEconomo==class],na.rm=TRUE))
for (class in 1:length(unique(vonEconomo))) {
  df.vonEconomo$classmu[df.vonEconomo$class == class] = classmu[class]
}

coloroption = seismic

# Figure 3A
pdf(paste0(plot.out, "vonEconomoMI.pdf"), height=3, width=3.5)
ggplot(df.vonEconomo, aes(x=class, y=mi, fill = classmu)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.1,fill="white") +
  scale_fill_gradientn(colours=coloroption, limits = c(-0.8, 0.8) ) +
  ylab('Maturational Index') +
  xlab('Von Economo Classes') +
  scale_x_discrete(rep('.',8)) +
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold", color = vonEconomo.colors, size = 12),
        axis.title=element_text(size=12))
dev.off()


# Yeo Network Decoding
df.yeo = data.frame(mi = diff.mi, class = yeo.classes)
colnames(df.yeo)[2] = 'class' 
df.yeo$class = as.factor(df.yeo$class)

classmu = sapply(1:length(unique(yeo)),function(class) mean(diff.mi[yeo==class],na.rm=TRUE))
for (class in 1:length(unique(yeo))) {
  df.yeo$classmu[df.yeo$class == class] = classmu[class]
}

# Figure 3B
pdf(paste0(plot.out, "yeoMI.pdf"), height=3, width=3.5)
ggplot(df.yeo, aes(x=class, y=mi, fill = classmu)) + 
  geom_violin(scale='width',color="white") +
  geom_boxplot(width=0.1,fill="white") +
  scale_fill_gradientn(colours=coloroption, limits = c(-0.8, 0.8) ) +
  ylab('Maturational Index') +
  scale_x_discrete(rep('.',8)) +
  theme_bw()+
  theme(axis.text.x = element_text(face = "bold", color = yeo.colors, size = 12),
        axis.title=element_text(size=12))
dev.off()

