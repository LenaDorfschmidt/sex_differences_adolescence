ztest_unpooled = function(beta1, beta2, SE1, SE2,N1,N2){
  unpooled_variance = (SE1/N1) + (SE2/N2) ## unpooled variance
  p.unpooled=(1-pnorm(q=abs(unname( beta1 - beta2 ))/sqrt(unpooled_variance)))/2
  p.fdr.unpooled = p.adjust(p.unpooled,method = 'fdr')
  return(data.frame(diff=(beta1-beta2), SE.unpooled = unpooled_variance, 
                    p.unpooled=p.unpooled,p.fdr.unpooled=p.fdr.unpooled)) 
}