library(scales)
library(ggplot2)
library(abind)
library(nlme)
library(extrafont)

data = 'NSPN' # data = 'CV_stratified' # data = 'glob.fc.regr' # data = 'regr_FD_by_sex' # data = 'motionmatched' # data = # data = 
script.path = paste0('Scripts/')

data.path = paste0('Data/connectivity_matrices/')
results.path = paste0('Results/')
# load in data: NSPN, NSPN low motion or GSR
load(paste0(data.path,'nspn.main.RData'))
plot.out = paste0(results.path,'mdd_decoding/')
data.path = paste0(results.path,'maturational_index/')


ifelse(!dir.exists(file.path(plot.out)), dir.create(file.path(plot.out)), FALSE)

load('Data/biodep/biodep.rds')
fc.all = readRDS('Data/biodep_noregression/postregSubWavCorMat_LD.rds')
demographic = readRDS('Data/biodep/Sociodemo129.rds')


nroi=346
nsub = dim(fc.all)[3]
str = matrix(ncol=dim(fc.all)[3],nrow=376)
for (s in 1:nsub) {
  FC = fc.all[,,s]
  str[,s] = rowMeans(FC, na.rm=TRUE)
}
str = str[hcp.keep.id,!(demographic$Categorisation.by.CRP.and.HAM.D == 'Depressed High CRP')]
demographic = demographic[!(demographic$Categorisation.by.CRP.and.HAM.D == 'Depressed High CRP'),]
chisq.test(table(demographic$Sex, demographic$MainArm))


t.case.control = p.case.control = AIC.m= vector(length=nroi)

for (roi in 1:nroi) {
  df = data.frame(fc = str[roi,], group = demographic$Categorisation.by.CRP.and.HAM.D, sex = demographic$Sex)
  m = lm(fc~group+sex, data = df)
  t.case.control[roi] = summary(m)$coefficients[2,3]
  AIC.m[roi] = AIC(m)
}

brainplot(t.case.control,nm,paste0('/Results/mdd_decoding/mdd.case.control'),-4,0,'plasma','stacked', FALSE)

p.case.control = pt(t.case.control,df=96)
p.fdr.case.control = p.adjust(pt(t.case.control,df=96),method = 'fdr')
print(sum(p.fdr.case.control<0.05))

write.table(p.fdr.case.control,paste0('Data/mdd//p.fdr.case.control.txt'),row.names = FALSE,col.names = FALSE)
write.table(p.case.control,paste0('Data/mdd/p.case.control.txt'),row.names = FALSE,col.names = FALSE)
write.table(t.case.control,paste0('Data/mdd/case.control.biodep.346.txt'),row.names = FALSE,col.names = FALSE)

