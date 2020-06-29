


perm_bl14_sl <- function(rand.sex){
  all.MRI_table = MRI_table
  all.MRI_table$sex = as.factor(rand.sex)
  nsub=8
  # overall, cortex-all, subcortex-all
  str.l.sl = str.l.14 = str.l.sl = str.l.sl.p = str.l.sl.t = matrix(NA,nrow=2,ncol=nroi) # overall node strength
  str.cort.l.sl = str.cort.l.14 = str.cort.l.sl.p = str.cort.l.sl.t = matrix(NA,nrow=2,ncol=nroi) # cortical node strength
  str.subc.l.sl = str.subc.l.14 = str.subc.l.sl.p = str.subc.l.sl.t = matrix(NA,nrow=2,ncol=nroi) # cortical node strength
  str.subc.cort.l.sl = str.subc.cort.l.14 = str.subc.cort.l.sl.p = str.subc.cort.l.sl.t = array(NA,dim=c(2,8,nroi)) # cortical node strength
  
  g_str = c('female','male')
  for (s in (1:2)) {
    gender_idx = (all.MRI_table$sex == (s-1))

    # LOCAL
    for (n in 1:nroi) {
      # str
      # df = data.frame(x1 = age[gender_idx], x2 = center[gender_idx], y = str[gender_idx,n], id = id[gender_idx]) # data frame for lme
      # l = lme(y ~ x1 + x2, random = ~ 1|id, data = df) # fit lme
      # str.l.sl.p[s,n] = summary(l)$tTable[2,5]    # p-value for effect of age
      # str.l.sl.t[s,n] = summary(l)$tTable[2,4]
      # str.l.sl[s,n] = l$coefficients$fixed[2]
      
      # # predicted value at 14 (includes "average" effect of sex (0.5), and average effect of each of three scanner sites (1/3))
      # str.l.14[n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      #str.cort.cort (variable names analogous to above)
      
      try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str.cort[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
      str.cort.l.sl.p[s,n] = summary(l)$tTable[2,5]
      str.cort.l.sl.t[s,n] = summary(l)$tTable[2,4]
      str.cort.l.14[s,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      str.cort.l.sl[s,n] = l$coefficients$fixed[2]
      })

      
      # str.cort.subc (variable names analogous to above)
      # str.subc.cort (variable names analogous to above)
      try({
      df = data.frame(x1 = MRI_table$mri_pseudo_age[gender_idx], x2 = MRI_table$mri_center[gender_idx], y = str.subc[gender_idx,n], id = MRI_table$id_nspn[gender_idx])
      l = lme(y ~ x1 + x2, random = ~ 1|id, data = df)
      str.subc.l.sl.p[s,n] = summary(l)$tTable[2,5]
      str.subc.l.sl.t[s,n] = summary(l)$tTable[2,4]
      str.subc.l.14[s,n] = l$coefficients$fixed[1] + l$coefficients$fixed[2]*14 + l$coefficients$fixed[3]*(1/3) + l$coefficients$fixed[4]*(1/3)
      str.subc.l.sl[s,n] = l$coefficients$fixed[2]
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
        })
      }
    }
  }
  
  return(list(str.l.sl,str.l.14,str.cort.l.sl,str.cort.l.14,str.subc.l.sl,str.subc.l.14,str.subc.cort.l.sl,str.subc.cort.l.14))
}




