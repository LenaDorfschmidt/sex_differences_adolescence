load('Data/downloaded/nspn.subj.info.RData')
load('Data/downloaded/nspn.fmri.general.vars.RData')

nm.all = nm$V1
nm = nm$V1[hcp.keep.id]

load('Data/downloaded/nspn.fmri.main.RData')
MRI_table$sex = as.factor(as.numeric(sex.main)-1)
MRI_table$mri_center = as.factor(as.numeric(site.main))

#### Main Data #####
fc = fc.main
if(!all(id.main==MRI_table$id_nspn)) errorCondition('Wait, something went wrong. Subject order does not match!')
save(MRI_table,fc,nroi,pos,nm,nm.all,hcp.keep.id,hcp.drop.id,yeo.colors, yeo.classes, vonEconomo.classes, vonEconomo.colors, glassersub,file='Data/connectivity_matrices/nspn.main.RData')

#### GSR #####
load('Data/downloaded/nspn.fmri.gsr.RData')
if(!all(id.gsr==MRI_table$id_nspn)) errorCondition('Wait, something went wrong. Subject order does not match!')
fc = fc.gsr
save(MRI_table,fc,nroi,pos,nm,nm.all,hcp.keep.id,hcp.drop.id,yeo.colors, yeo.classes, vonEconomo.classes, vonEconomo.colors, glassersub,file='Data/connectivity_matrices/nspn.gsr.RData')

#### Regr by sex #####
load('Data/downloaded/nspn.repl.FD.by.sex.RData')
fc = fc.regr.by.sex
save(MRI_table,fc,nroi,pos,nm,nm.all,hcp.keep.id,hcp.drop.id,yeo.colors, yeo.classes, vonEconomo.classes, vonEconomo.colors, glassersub,file='Data/connectivity_matrices/nspn.repl.FD.by.sex.RData')

#### Motion Matched ####
load('Data/downloaded/nspn.motion.matched.RData')
fc = fc.main[,,motion.matched.ids]
MRI_table = MRI_table[motion.matched.ids,]
save(MRI_table,fc,nroi,pos,nm,nm.all,hcp.keep.id,hcp.drop.id,yeo.colors, yeo.classes, vonEconomo.classes, vonEconomo.colors, glassersub,file='Data/connectivity_matrices/nspn.motion.matched.RData')




