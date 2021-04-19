predict.LME.SE = function(l,coefvec){
  coefs = l$coefficients$fixed
  cc = as.matrix(expand.grid(data.frame(t(coefvec))))
  pred = coefs %*% t(cc)
  vc <- vcov(l)
  SE<-sqrt(diag(cc %*% vc %*% t(cc)))
  return(c(pred, SE))
}