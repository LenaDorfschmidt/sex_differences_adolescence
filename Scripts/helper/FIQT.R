FIQT <- function(z=z, min.p=10^-300){
  pvals<-2*pnorm(abs(z),low=F)
  pvals[pvals<min.p]<- min.p
  adj.pvals<-p.adjust(pvals,method="fdr")
  mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
  mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
  mu.z
}
