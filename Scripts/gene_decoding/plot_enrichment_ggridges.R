plot_enrichment_ggridges <- function(gene_stats, gene_null, clustered, min, max){
  if(clustered==T){
    sign.cells = gene_stats$name[p.adjust(gene_stats$p, method='fdr')<0.05]
    gene.thresh = gene_null[!is.na(match(gene_null$list,sign.cells)),]
    
    gene.thresh$class <- 'Other'
    gene.thresh$class[grep('Ex',gene.thresh$list)] <- 'Exitatory'
    gene.thresh$class[grep('In',gene.thresh$list)] <- 'Inhibitory'
    
    gene.plot = gene.thresh[with(gene.thresh, order(class, real)),]
    #gene.plot = gene.thresh[plot.order,]
    gene.plot$list = factor(gene.plot$list, levels=unique(gene.plot$list))
    
    #'Results/gene_decoding/postnatal.celltype.enrichment.pdf' # min=-3000 max = 3000
    p<-ggplot(gene.plot,aes(x = null, y = list)) +
      geom_density_ridges_gradient(aes(fill = ..x..), scale = 3, size = 0.3,rel_min_height = 0.001,bandwidth=150) +
      geom_point(aes(y=list,x=real, fill=real),size=4,shape=21) +
      scale_fill_gradientn(colours = c("blue", "white", "red"), lim = c(min,max),oob=squish)+
      xlab('Median Rank')+
      ylab('') +
      theme_minimal() +
      theme(legend.title = element_blank(),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(size=12),
            axis.title.y = element_blank())
  }else{
    gene_density <- density(gene_null$gene_set_median_rand)
    gene_plot <- data.frame(x = gene_density$x, y = gene_density$y, fill=gene_density$x, real=unique(gene_null$real))
    
    p <- ggplot(gene_plot)+ geom_segment(aes(x = x, xend = x, y = 0, yend = y, color = fill)) +
      geom_line(aes(x=x, y=y)) +
      geom_point(aes(x=real, y=0,fill=real),size=4,shape=21) +
      scale_color_gradientn(colours = c("blue", "white", "red"),
                            limits=c(min,max))+
      scale_fill_gradientn(colours = c("blue", "white", "red"),
                           limits=c(min,max))+
      theme_minimal() +
      xlab('Median Rank')+
      theme(legend.position="none",
            axis.text.x = element_text(size=12),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size=12),
            axis.title.y = element_blank())
  }
  return(p)
}


