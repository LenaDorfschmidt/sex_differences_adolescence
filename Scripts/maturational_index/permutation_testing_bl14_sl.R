library("lme4")
library(vows)       


data = 'NSPN' 
base.path = '~/Documents/Education/Cambridge/fMRI_NSPN/'
script.path = paste0(base.path,'Scripts/')

source(paste0(script.path,'helper/brainplot.R'))
load(paste0(script.path,'helper/glassersubcortex_R.bin'))

data.path = paste0(base.path,'Data/connectivity_matrices/')
results.path = paste0(base.path,'Results/')
plot.out = paste0(results.path,'maturational_index/perm.bl14.sl/')

load(paste0(data.path,'main.RData'))

## Setup permutation matrix
nPerm <- 1000
MAT <- matrix( NA_integer_, nrow=NROW(MRI_table), ncol=nPerm )
for(IDX in 1:NCOL(MAT) ) {
  for( occ in 1:3 ) {
    for( grp in 1:5 ) {
      WHICH <- which( MRI_table$event==occ & MRI_table$agebin==grp) 
      MAT[ WHICH, IDX ] <- sample( MRI_table$sex[ WHICH ] )
    }
  }
}

# Make sure all went well
FTAB1 <- ftable(xtabs( ~ sex + event + agebin, MRI_table ),col.vars=c(2,1))
FTAB2 <- ftable(xtabs( ~ MAT[,1] + MRI_table$event + MRI_table$agebin, MRI_table ),col.vars=c(2,1))
xtabs( ~ MAT[,1] + MRI_table$event )

write.table(MAT,paste0(results.path, 'perm.bl14.sl/perm.matrix.loc.bl14.sl.txt'),row.names = FALSE,col.names = FALSE)

str = apply(fc,c(3,1),function(x) mean(x,na.rm=TRUE))
str.cort = apply(fc[,17:346,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc = apply(fc[,1:16,],c(3,1),function(x) mean(x,na.rm=TRUE))
str.subc.cort = lapply(1:8, function(i) apply(fc[c(i,i+8),,],c(2,3),function(x) mean(x,na.rm=TRUE)))

fc.m <- colMeans(str)
nroi=346
nsub=8
  
source(paste0(script.path, 'maturational_index/perm_bl14_sl.R'))
perm.mat = read.table(paste0(results.path, 'maturational_index/perm_bl14_sl/perm.matrix.loc.bl14.sl.txt'))

nPerm = 1000
perm.bl14 = perm.sl = cort.bl14 = cort.sl = sub.bl14 = sub.sl = array(NA,dim=c(nPerm,2,346))
sub.cort.sl = sub.cort.bl14 = array(NA, dim=c(nPerm,2,8,346))
for (perm in 901:nPerm) {
  print(perm)
  rand.sex = perm.mat[,perm]
  output <- try({
    output <- perm_bl14_sl(rand.sex)
    perm.sl[perm,,] <- output[[1]]
    perm.bl14[perm,,] <- output[[2]]
    cort.bl14[perm,,] <- output[[3]]
    cort.sl[perm,,] <- output[[4]]
    sub.bl14[perm,,] <- output[[5]]
    sub.sl[perm,,] <- output[[6]]
    sub.cort.sl[perm,,,] <- output[[7]]
    sub.cort.bl14[perm,,,] <- output[[8]]
  })
}

all.cort.bl14 = all.cort.sl = all.sub.bl14 = all.sub.sl = array(NA,dim=c(nPerm,2,346))
all.sub.cort.sl = all.sub.cort.bl14 = array(NA, dim=c(nPerm,2,8,346))

for (perm in seq(100,1000,100)) {
  load(paste0('~/Desktop/perm_',perm,'.RData'))
  all.cort.bl14[1:perm,,] = cort.bl14[((perm-99):perm),,]
  all.cort.sl[1:perm,,] = cort.sl[((perm-99):perm),,]
  all.sub.bl14[1:perm,,] = sub.bl14[((perm-99):perm),,]
  all.sub.sl[1:perm,,] = sub.sl[((perm-99):perm),,]
  all.sub.cort.bl14[1:perm,,,] = sub.cort.bl14[((perm-99):perm),,,]
  all.sub.cort.sl[1:perm,,,] = sub.cort.sl[((perm-99):perm),,,]
}

save(all.cort.bl14,all.cort.sl,all.sub.bl14,all.sub.sl,
     all.sub.cort.bl14,all.sub.cort.sl,
     file=paste0(results.path, 'maturational_index/perm_bl14_sl/perm_bl14_sl.RData'))

#nPerm=sum(apply(is.na(perm.sl[,1,]),1,sum)==0)
rand.diff.bl14 = perm.bl14.1-perm.bl14.2
rand.diff.cort.bl14 = cort.bl14[,1,]-cort.bl14[,2,]
rand.diff.sub.bl14 = sub.bl14[,1,]-sub.bl14[,2,]

rand.diff.sl = perm.sl[,1,]-perm.sl[,2,]
rand.diff.cort.sl = cort.sl[,1,]-cort.sl[,2,]
rand.diff.sub.sl = sub.sl[,1,]-sub.sl[,2,]

rand.diff.sub.cort.bl14 = rand.diff.sub.cort.sl = array(NA, dim=c(nPerm,8,346))
for (sub in 1:nsub) {
  rand.diff.sub.cort.bl14[,sub,] = all.sub.cort.bl14[,1,sub,]-all.sub.cort.bl14[,2,sub,]
  rand.diff.sub.cort.sl[,sub,] = all.sub.cort.sl[,1,sub,]-all.sub.cort.sl[,2,sub,]
}


plot.out = '~/Documents/Education/Cambridge/fMRI_NSPN/Results/maturational_index/perm_bl14_sl/'
plot.in = '~/Documents/Education/Cambridge/fMRI_NSPN/Results/maturational_index/'
m.14 = as.matrix(read.csv(paste0(plot.in,'/male/all_bl14.txt')))[,2]
f.14 = as.matrix(read.csv(paste0(plot.in,'/female/all_bl14.txt')))[,2]
m.sl = as.matrix(read.csv(paste0(plot.in,'/male/all_delta_age.txt')))[,2]
f.sl = as.matrix(read.csv(paste0(plot.in,'/female/all_delta_age.txt')))[,2]

sub.m.14 = as.matrix(read.csv(paste0(plot.in,'/male/subc_bl14.txt')))[,2]
sub.f.14 = as.matrix(read.csv(paste0(plot.in,'/female/subc_bl14.txt')))[,2]
sub.m.sl = as.matrix(read.csv(paste0(plot.in,'/male/subc_delta_age.txt')))[,2]
sub.f.sl = as.matrix(read.csv(paste0(plot.in,'/female/subc_delta_age.txt')))[,2]

cort.m.14 = as.matrix(read.csv(paste0(plot.in,'/male/cort_bl14.txt')))[,2]
cort.f.14 = as.matrix(read.csv(paste0(plot.in,'/female/cort_bl14.txt')))[,2]
cort.m.sl = as.matrix(read.csv(paste0(plot.in,'/male/cort_delta_age.txt')))[,2]
cort.f.sl = as.matrix(read.csv(paste0(plot.in,'/female/cort_delta_age.txt')))[,2]

sub.cort.m.14 = as.matrix(read.csv(paste0(plot.in,'/male/sub_cort_bl14.csv')))[,2:347]
sub.cort.f.14 = as.matrix(read.csv(paste0(plot.in,'/female/sub_cort_bl14.csv')))[,2:347]
sub.cort.m.sl = as.matrix(read.csv(paste0(plot.in,'/male/sub_cort_delta_age.csv')))[,2:347]
sub.cort.f.sl = as.matrix(read.csv(paste0(plot.in,'/female/sub_cort_delta_age.csv')))[,2:347]

diff.14 = f.14-m.14
diff.sl = f.sl-m.sl
sub.diff.14 = sub.f.14-sub.m.14
sub.diff.sl = sub.f.sl-sub.m.sl
cort.diff.14 = cort.f.14-cort.m.14
cort.diff.sl = cort.f.sl-cort.m.sl
sub.cort.diff.14 = sub.cort.f.14-sub.cort.m.14
sub.cort.diff.sl = sub.cort.f.sl-sub.cort.m.sl

p.bl.14 = p.sl = cort.p.bl.14 = cort.p.sl = sub.p.bl.14 = sub.p.sl = vector(length=nroi)
sub.cort.p.bl.14 = sub.cort.p.sl = matrix(ncol=nroi, nrow=nsub)
for (roi in 1:nroi) {
  #p.bl.14[roi] = sum(abs(rand.diff.bl14[,roi])>abs(diff.14[roi]),na.rm=TRUE)/nPerm  
  #p.sl[roi] = sum(rand.diff.sl[,roi]>diff.sl[roi],na.rm=TRUE)/nPerm 
  cort.p.bl.14[roi] = sum(abs(rand.diff.cort.bl14[,roi])>abs(cort.diff.14[roi]),na.rm=TRUE)/nPerm  
  cort.p.sl[roi] = sum(abs(rand.diff.cort.sl[,roi])>abs(cort.diff.sl[roi]),na.rm=TRUE)/nPerm 
  sub.p.bl.14[roi] = sum(abs(rand.diff.sub.bl14[,roi])>abs(sub.diff.14[roi]),na.rm=TRUE)/nPerm  
  sub.p.sl[roi] = sum(abs(rand.diff.sub.sl[,roi])>abs(sub.diff.sl[roi]),na.rm=TRUE)/nPerm 
  for (sub in 1:nsub) {
    sub.cort.p.bl.14[sub,roi] = sum(abs(rand.diff.sub.cort.bl14[,sub,roi])>abs(sub.cort.diff.14[sub,roi]),na.rm=TRUE)/nPerm  
    sub.cort.p.sl[sub,roi] = sum(abs(rand.diff.sub.cort.sl[,sub,roi])>abs(sub.cort.diff.sl[sub,roi]),na.rm=TRUE)/nPerm 
  }
}


tmp.cort.sl = ifelse((p.adjust(cort.p.sl[17:346],method='fdr')<0.05), cort.diff.sl[17:346], NA)
tmp.cort.bl14 = ifelse((p.adjust(cort.p.bl.14[17:346],method='fdr')<0.05), cort.diff.14[17:346], NA)
tmp.sub.sl = ifelse((p.adjust(sub.p.sl[17:346],method='fdr')<0.05), sub.diff.sl[17:346], NA)
tmp.sub.bl14 = ifelse((p.adjust(sub.p.bl.14[17:346],method='fdr')<0.05), sub.diff.14[17:346], NA)



### For Publ

brainplot(ifelse(p.adjust(cort.p.bl.14[17:346],method='fdr')<0.05,cort.diff.14[17:346],NA),nm[17:346],paste0(plot.out,'plots/1000_diff_cort_bl14_sign_only'),-0.15,0.15,
          'blue2red_harsh','stacked',FALSE)
brainplot(ifelse(p.adjust(sub.p.bl.14[17:346],method='fdr')<0.05,sub.diff.14[17:346],NA),nm[17:346],paste0(plot.out,'plots/1000_diff_sub_bl14_sign_only'),-0.15,0.15,
          'blue2red_harsh','stacked',FALSE)


for (sub in 1:nsub) {
  tmp.sl = ifelse((p.adjust(sub.cort.p.sl[sub,17:346],method='fdr')<0.05), sub.cort.diff.sl[sub,17:346], NA)
  tmp.bl14 = ifelse((p.adjust(sub.cort.p.bl.14[sub,17:346],method='fdr')<0.05), sub.cort.diff.14[sub,17:346], NA)
  if (any(!(is.na(tmp.bl14)))) {
    brainplot(sub.cort.diff.14[sub,17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_subc_cort_bl14_',nm[sub]),-0.15,0.15,
              'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(sub.cort.p.bl.14[sub,17:346],method='fdr')<0.05,1,NA))
    
    #   brainplot(tmp.bl14,nm[17:346],paste0(plot.out,'plots/diff_subc_cort_bl14_',nm[sub]),-0.25,0.25,'blue2red_harsh','stacked')
  }
  if (any(!(is.na(tmp.sl)))) {
    brainplot(ifelse(p.adjust(sub.cort.p.sl[sub,17:346],method='fdr')<0.05,sub.cort.diff.sl[sub,17:346],NA),nm[17:346],paste0(plot.out,'plots/1000_diff_subc_cort_sl_sign_only',nm[sub]),-0.025,0.025,
              'blue2red_harsh','stacked',FALSE)
  }
}




if (any(!(is.na(tmp.cort.bl14)))) {
  brainplot(cort.diff.14[17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_cort_bl14_'),-0.15,0.15,
            'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(cort.p.bl.14[17:346],method='fdr')<0.05,1,NA))
}
if(any(!(is.na(tmp.sub.bl14)))){
  brainplot(sub.diff.14[17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_sub_bl14_'),-0.15,0.15,
            'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(sub.p.bl.14[17:346],method='fdr')<0.05,1,NA))
}
if (any(!(is.na(tmp.cort.sl)))) {
  brainplot(cort.diff.sl[17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_cort_sl_'),-0.025,0.025,
            'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(cort.p.sl[17:346],method='fdr')<0.05,1,NA))
}
if (any(!(is.na(tmp.sub.sl)))) {
  brainplot(sub.diff.sl[17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_sub_sl_'),-0.025,0.025,
            'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(sub.p.sl[17:346],method='fdr')<0.05,1,NA))
}


for (sub in 1:nsub) {
  tmp.sl = ifelse((p.adjust(sub.cort.p.sl[sub,17:346],method='fdr')<0.05), sub.cort.diff.sl[sub,17:346], NA)
  tmp.bl14 = ifelse((p.adjust(sub.cort.p.bl.14[sub,17:346],method='fdr')<0.05), sub.cort.diff.14[sub,17:346], NA)
  if (any(!(is.na(tmp.bl14)))) {
    brainplot(sub.cort.diff.14[sub,17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_subc_cort_bl14_',nm[sub]),-0.15,0.15,
              'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(sub.cort.p.bl.14[sub,17:346],method='fdr')<0.05,1,NA))
    
    #   brainplot(tmp.bl14,nm[17:346],paste0(plot.out,'plots/diff_subc_cort_bl14_',nm[sub]),-0.25,0.25,'blue2red_harsh','stacked')
  }
  if (any(!(is.na(tmp.sl)))) {
    brainplot(sub.cort.diff.sl[sub,17:346],nm[17:346],paste0(plot.out,'plots/1000_diff_subc_cort_sl_',nm[sub]),-0.025,0.025,
              'blue2red_harsh','stacked',FALSE,ifelse(p.adjust(sub.cort.p.sl[sub,17:346],method='fdr')<0.05,1,NA))
    
    #   brainplot(tmp.sl,nm[17:346],paste0(plot.out,'plots/diff_subc_cort_sl_',nm[sub]),-0.025,0.025,'blue2red_harsh','stacked')
  }
}


binary.bl.14 = ifelse(fdr.p.bl.14<0.05,1,2)
binary.sl = ifelse(fdr.p.sl<0.05,1,2)

write.table(binary.bl.14,paste0(plot.out, 'sign.bl.14.txt'),row.names = FALSE,col.names = FALSE)
write.table(binary.sl,paste0(plot.out, 'sign.sl.txt'),row.names = FALSE,col.names = FALSE)
write.table(p.bl.14,paste0(plot.out, 'p.bl.14.14.txt'),row.names = FALSE,col.names = FALSE)
write.table(p.sl,paste0(plot.out, 'p.sl.txt'),row.names = FALSE,col.names = FALSE)

write.table(forSurfacePlot(binary.bl.14),paste0(plot.out, 'b4p.sign.bl.14.txt'),row.names = FALSE,col.names = FALSE)
write.table(forSurfacePlot(binary.sl),paste0(plot.out, 'b4p.sign.sl.txt'),row.names = FALSE,col.names = FALSE)

p.bl.14 = read.csv(paste0(plot.out, 'p.bl.14.14.txt'))[,1]
thresh.bl.14 = ifelse(fdr.p.bl.14<0.05,fdr.p.bl.14,NA)


brainplot(m.14,nm,paste0(plot.out,'plots/m.14'),0.2,0.6,'plasma','stacked')
brainplot(f.14,nm,paste0(plot.out,'plots/f.14'),0.2,0.6,'plasma','stacked')
brainplot(m.sl,nm,paste0(plot.out,'plots/m.sl'),-0.01,0.01,coloroption = NULL,'stacked')
brainplot(f.sl,nm,paste0(plot.out,'plots/f.sl'),-0.01,0.01,coloroption = NULL,'stacked')

brainplot(diff.14,nm,paste0(plot.out,'plots/diff.14'),-0.1,0.1,'viridis','stacked')
bl14.thresh = diff.14    
bl14.thresh[!(fdr.p.bl.14<0.05)] = NA
brainplot(bl14.thresh,nm,paste0(plot.out,'plots/diff.14.thresholded'),-0.1,0.11,'viridis','stacked')

brainplot(diff.sl,nm,paste0(plot.out,'plots/diff.sl'),-0.01,0.01,'viridis','stacked')
sl.thresh = diff.sl  
sl.thresh[!(fdr.p.sl<0.05)] = NA
brainplot(sl.thresh,nm,paste0(plot.out,'plots/diff.sl.thresholded'),-0.1,0.1,'viridis','stacked')
