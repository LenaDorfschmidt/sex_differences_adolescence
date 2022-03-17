run_gene_enrichment <- function(geneslist, pls1, clustered, genes, n_perm, out_path, listname=NULL){
  if (clustered==TRUE) {
    gene_null = data.frame(gene_set_median_rand=NA,name=NA,real=NA,z=NA,p=NA)
    gene_length_stats = data.frame(cluster = NA, p=NA, s=NA,z=NA,n=NA,approach=NA)
    for (cluster in unique(geneslist$Cluster)) {
      try({
        gene_set <- subset(geneslist, Cluster == cluster)$Gene
        out <- gene_nulls_nearest_neighbour(pls1$hgnc_symbol,pls1$z_uncorr,gene_set,genes,cluster,n_perm,rank=T,out_path)
        gene_null = rbind(gene_null,out[[1]])
        gene_length_stats = rbind(gene_length_stats,c(data.frame(cluster=cluster),out[[2]]))
      })
    }   
    gene_null = gene_null[-1,] # First entry was just dummy
    gene_length_stats = gene_length_stats[-1,] # First entry was just dummy
    
    gene_stats = gene_null[match(unique(gene_null$name),gene_null$name),2:5] # to view stats/needed for fdr-correction
    
    colnames(gene_null) <- c("null","list","real","z","p")
    gene_null$list <- as.factor(gene_null$list)
    gene_null = gene_null[-1,] # First entry was just dummy
  }else{
    out <- gene_nulls_nearest_neighbour(pls1$hgnc_symbol,pls1$z_uncorr,geneslist$Gene,genes,listname,n_perm,rank=T,out_path)
    gene_null <- out[[1]]
    gene_length_stats <- out[[2]]
    gene_stats=gene_null[1,c(2:5)]
  }
  return(list(gene_null = gene_null, gene_stats = gene_stats, gene_length_stats = gene_length_stats))
}
