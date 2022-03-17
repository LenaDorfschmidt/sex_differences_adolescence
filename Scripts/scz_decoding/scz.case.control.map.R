library(here)
library(RNifti)
library(ggseg)
library(ggplot2)
library(viridis)
library('RColorBrewer')

load(here('Data/connectivity_matrices/nspn.main.RData'))
demographics = read.csv(here('Data/scz/COBRE_demographics.csv'))

region.IDs = read.table(here('Data/SCZ/fMRI/HPC+subcort.region.IDs.txt'), header = T)
col.idx = region.IDs$roi_value[region.IDs$include==1]

nsub=length(demographics$SubjectID)
nts = 150
nroi=376
ts = array(NA, dim=c(150,376,nsub))
for (s in 1:nsub) {
  ts.tmp = read.table(here(paste0('data/SCZ/fMRI/preproc/',demographics$SubjectID[s],'/ppp_HPC__ts.txt')))
  ts[,,s] = as.matrix(ts.tmp[,col.idx])
}

sd.signal=apply(ts, c(2,3), sd)
subj.id.excl = c(144,38)  # Which subjects are leading to dropout? Subject 28 and 144 are causing all 1 and 2 pp dropout ROIs
# which(rowSums(sd.signal==0)==1) # 
# which(sd.signal[80,]==0)

demographics = demographics[-subj.id.excl,]
ts = ts[,,-subj.id.excl]

sd.signal=apply(ts, c(2,3), sd)
dropout.roi = which(rowSums(sd.signal==0)>0)
ts = ts[,-dropout.roi,]


mean.signal=apply(ts, c(2,3), mean)
Z = scale(mean.signal) #rowSums(Z>2.58) Is this concerning?


nsub = dim(ts)[3]
nroi=dim(ts)[2]

nedges = nroi*(nroi-1)/2
n.10.perc = round(nedges*1/10)

fc = fc.thresh = fc.thresh.bin = array(NA, dim=c(nroi,nroi,nsub))
for (s in 1:(nsub)) {
  mat = cor(scale(ts[,,s]))
  diag(mat) <- NA
  fc[,,s] = mat
  upper.tri.vec = mat[upper.tri(mat)]
  mat[mat<sort(upper.tri.vec)[nedges-n.10.perc]] = 0 
  fc.thresh[,,s] = mat
  fc.thresh.bin[,,s] = as.numeric(mat!=0)
}

str.thresh = apply(fc.thresh,c(2,3),function(x) mean(x,na.rm=T))
str = apply(fc,c(2,3),function(x) mean(x,na.rm=T))

tmp = cor(str.thresh, str)
pdf(here('Results/scz/effect_of_thresholding_ROI.pdf'),4,4)
ggplot(data.frame(cor.val = diag(tmp),group = demographics$group=='No_Known_Disorder'),aes(x=cor.val,fill=group))+
  geom_histogram()+theme_bw()+xlab('Correlation between Subject ROI Strength') + ylab('Count')+labs(fill='SCZ')
dev.off()

p.roi = t.roi = array(NA, dim=c(nroi,5))

for (r in 1:nroi) {
  df = data.frame(fc=str[r,], group=!(demographics$group=='No_Known_Disorder'),fd = demographics$mean_fd,age=demographics$age,sex=demographics$sex)
  lm.roi = lm(fc~group+fd+age+sex,data=df)
  p.roi[r,] = summary(lm.roi)$coefficients[,4]
  t.roi[r,] = summary(lm.roi)$coefficients[,3]
}

p.scz.fdr = p.adjust(p.roi[,2], method='fdr')
t.scz = t.roi[,2]

t.scz.376 = p.scz.fdr.376 = p.scz.376 = array(NA, dim=c(376,1))
p.scz.fdr.376[which(is.na(match(seq(1:376),dropout.roi))),1] = p.scz.fdr
p.scz.376[which(is.na(match(seq(1:376),dropout.roi))),1] = p.roi[,2]
t.scz.376[which(is.na(match(seq(1:376),dropout.roi))),1] = t.scz

write.table(p.scz.fdr.376,paste0('Data/scz/p.fdr.case.control.txt'),row.names = FALSE,col.names = FALSE)
write.table(p.scz.376,paste0('Data/scz/p.case.control.txt'),row.names = FALSE,col.names = FALSE)
write.table(t.scz.376,paste0('Data/scz/t.case.control.txt'),row.names = FALSE,col.names = FALSE)

pdf('Results/SCZ/scz.cort.match.Sarah.pdf',5,2)
ggseg(.data=data.frame(value = t.scz,label=nm.ggseg),mapping=aes(fill=value), atlas=glassersub, position = 'stacked') + scale_fill_gradientn(colours=brewer.pal(100,'RdBu'))
dev.off()

pdf('Results/SCZ/scz.cort.viridis.pdf',5,2)
ggseg(.data=data.frame(value = t.scz,label=nm.ggseg),mapping=aes(fill=value), atlas=glassersub, position='stacked') + scale_fill_gradientn(colours=viridis(100))
dev.off()




