gene_nulls_nearest_neighbour <- function(gene_list, gene_stat, gene_set,genes,name, n_perm, rank=T, outdir){
  if (rank==T){
    gene_stat <- rank(gene_stat)
  }
  print(paste0("Current geneset: ", name))
  # Calculate real median ranl
  gene_set_median <- median(gene_stat[match(gene_set,gene_list,nomatch=0)],na.rm=T)-median(gene_stat)
  print(paste0("Empirical Median Rank = ", gene_set_median))
  # Create null distribution matched for gene length
  gene_set_length <- median(genes$length[match(gene_set,gene_list,nomatch=0)],na.rm=T) # Real median gene length
  gene_set_length_sd <- sd(genes$length[match(gene_set,gene_list,nomatch=0)],na.rm=T) # Real sd gene length
  # get a matched sample of genes not including the ones in the set with a ratio of 3:!
  genes2 <- genes
  genes2$group <- (genes$hgnc_symbol %in% gene_set)
  zz <- match.data(matchit(group ~ length, data=genes2, method="nearest", 
                           distance="mahalanobis", replace=TRUE,full = T, ratio = 5))
  genes2 <- subset(zz, group == 'FALSE')
  
  # set a default p-value before entering the loop
  p.value = 1
  count = 1
  while ((!((p.value > 0.025)&(p.value < 0.975)))&(count<11)) {
    gene_set_rand <- sapply(1:n_perm, function(x) list(sample(genes2$hgnc_symbol, size = length(gene_set),replace = T)))
    gene_set_lengths_rand <- sapply(1:n_perm, function(x) median(genes$length[match(gene_set_rand[[x]],gene_list)],na.rm=T))
    #hist(gene_set_lengths_rand,50)
    p.value = sum(gene_set_lengths_rand < gene_set_length)/n_perm
    Norm.svalue = shapiro.test(gene_set_lengths_rand)$p.value
    Norm.wvalue = shapiro.test(gene_set_lengths_rand)$statistic
    z.value = (gene_set_length - mean(gene_set_lengths_rand))/sd(gene_set_lengths_rand)
    count = count+1
    approach = 'standard'
  }
  
  # Restrict to only 2 nearest neighbours
  if(!((p.value > 0.025)&(p.value < 0.975))) {
    genes2 <- genes
    genes2$group <- (genes$hgnc_symbol %in% gene_set)
    zz <- match.data(matchit(group ~ length, data=genes2, method="nearest", 
                             distance="mahalanobis", replace=TRUE,full = T, ratio = 2))
    genes2 <- subset(zz, group == 'FALSE')
    
    # set a default p-value before entering the loop
    p.value = 1
    count = 1
    #while(!((p.value > 0.025)&(p.value < 0.975))){
    while ((!((p.value > 0.025)&(p.value < 0.975)))&(count<11)) {
      gene_set_rand <- sapply(1:n_perm, function(x) list(sample(genes2$hgnc_symbol, size = length(gene_set),replace = T)))
      gene_set_lengths_rand <- sapply(1:n_perm, function(x) median(genes$length[match(gene_set_rand[[x]],gene_list)],na.rm=T))
      #hist(gene_set_lengths_rand,50)
      p.value = sum(gene_set_lengths_rand < gene_set_length)/n_perm
      Norm.svalue = shapiro.test(gene_set_lengths_rand)$p.value
      Norm.wvalue = shapiro.test(gene_set_lengths_rand)$statistic
      z.value = (gene_set_length - mean(gene_set_lengths_rand))/sd(gene_set_lengths_rand)
      count = count+1
      approach = '2 nearest neighbours'
    }
  }
  
  # Rerun without replacement
  if(!((p.value > 0.025)&(p.value < 0.975))) {
    genes2 <- genes
    genes2$group <- (genes$hgnc_symbol %in% gene_set)
    zz <- match.data(matchit(group ~ length, data=genes2, method="nearest", 
                             distance="mahalanobis", replace=FALSE,full = T, ratio = 3))
    genes2 <- subset(zz, group == 'FALSE')
    
    # set a default p-value before entering the loop
    p.value = 1
    count = 1
    #while(!((p.value > 0.025)&(p.value < 0.975))){
    while ((!((p.value > 0.025)&(p.value < 0.975)))&(count<11)) {
      gene_set_rand <- sapply(1:n_perm, function(x) list(sample(genes2$hgnc_symbol, size = length(gene_set),replace = T)))
      gene_set_lengths_rand <- sapply(1:n_perm, function(x) median(genes$length[match(gene_set_rand[[x]],gene_list)],na.rm=T))
      #hist(gene_set_lengths_rand,50)
      p.value = sum(gene_set_lengths_rand < gene_set_length)/n_perm
      Norm.svalue = shapiro.test(gene_set_lengths_rand)$p.value
      Norm.wvalue = shapiro.test(gene_set_lengths_rand)$statistic
      z.value = (gene_set_length - mean(gene_set_lengths_rand))/sd(gene_set_lengths_rand)
      count = count+1
      approach = 'without replacement'
    }
  }
  
  # print some stats to the console
  print(paste0("Number of Genes to pull from for null sets = ",nrow(genes2)))
  print(paste0("Number of Genes to pull = ",sum(!is.na(match(gene_set,gene_list)))))
  print(paste0("Gene Length P = ",p.value))
  # Output gene length histogram
  pdf(paste0(outdir,'/hist_gene_length_',name,'.pdf'))
  hist(gene_set_lengths_rand, breaks = 15)
  abline(v=gene_set_length,col="red") 
  dev.off()
  # Save stats to object
  length_stats = data.frame(p = p.value, s = Norm.svalue, z = z.value, n = nrow(genes2), approach = approach)
  # Calculate random median ranks difference
  gene_set_median_rand <- sapply(1:n_perm, function(x) median(gene_stat[match(gene_set_rand[[x]],gene_list,nomatch=0)],na.rm=T)-median(gene_stat))
  # Create ouput object
  enrich <- as.data.frame(gene_set_median_rand)
  enrich$name <- name
  enrich$real <- gene_set_median
  # Calculate one-sided p-value and z=scpre
  if (gene_set_median < mean(gene_set_median_rand)){ 
    enrich$z = qnorm(sum(gene_set_median_rand < gene_set_median)/n_perm)
    enrich$p = sum(gene_set_median_rand < gene_set_median)/n_perm
  }
  else {
    enrich$z = 1-qnorm(sum(gene_set_median_rand > gene_set_median)/n_perm)
    enrich$p = sum(gene_set_median_rand > gene_set_median)/n_perm
  }
  print(paste0("Enrichment p: ", enrich$p[1]))
  print(paste0("Enrichment z: ", enrich$z[1]))
  # Output enrichment histogram
  pdf(paste0(outdir,'/hist_enrichment_',name,'.pdf'))
  hist(gene_set_median_rand)
  abline(v=gene_set_median,col="red") 
  dev.off()
  return(list(enrich, length_stats))
}