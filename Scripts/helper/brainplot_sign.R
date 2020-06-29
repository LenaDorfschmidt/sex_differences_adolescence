library(ggseg)
library(ggsegExtra)
library(RColorBrewer)
library(scales)
library(viridis)
library(colormap)

brainplot = function(values,nm.labels,outpath=NULL,lower=NULL,upper=NULL,coloroption=NULL,position=NULL,disp_nan=FALSE){
  #brainplot = function(values,nm.labels,outpath=NULL,lower=NULL,upper=NULL,coloroption=NULL,position=NULL,disp_nan=FALSE,sig=NULL){
  if(is.null(min)){
    lower=min(values,na.rm=T)
  }
  if(is.null(position)){
    position="stacked"
  }
  if(is.null(upper)){
    upper=max(values,na.rm=T)
  }
  if(is.null(coloroption)){
    coloroption = c("blue","white","red")
  }
  if(all(coloroption=='seismic')){
    coloroption = read.csv('~/Documents/Education/Cambridge/fMRI_NSPN/Scripts/helper/seismic.txt',header = F)[,1]
  }
  if(all(coloroption=='PuOr')){
    coloroption = brewer.pal(10,'PuOr')
  }
  if(all(coloroption=='viridis')){
    coloroption = viridis_pal()(1000)
  }
  if(all(coloroption=='blue2green')){
    coloroption=colorRamps::blue2green(1000) 
  }
  if(all(coloroption=='plasma')){
    coloroption=viridis_pal(option = "plasma")(1000)
  }
  if(all(coloroption=='RdBu')){
    coloroption=colormap(colormap='RdBu')
  }
  if(all(coloroption=='yeo')){
    coloroption = c('#660066','#0066CC','#009900','#FF00FF','#FFFF99','#FF8000','#FF6666','#404040')
  }
  if(all(coloroption=='vonEconomo')){
    coloroption =c('#6f057a','#220ced','#067a08','#f59505','#fffb00','#00fffb','#f200ff','#404040')
  }
  if(all(coloroption=='blue2red_harsh')){
    coloroption = colorRamps::blue2red(1000)
  }
  # if (is.null(sig)) {
  #   sig = vector(length=length(nm))
  # }
  load('~/Documents/Education/Cambridge/general_scripts/glassersubcortex_R.bin')
  
  cortex.values = values[which(!is.na(values))]
  cortex.lables = nm.labels[which(!is.na(values))]
  #cortex.lables = glassersubcortex$label[match(toupper(nm.labels[!is.na(values)]), toupper(glassersubcortex$label))]
  
  cortex = data.frame(
    label = cortex.lables,
    val = cortex.values,
    stringsAsFactors = FALSE)#,
  # sig=sig)
  if (disp_nan) {
    c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val), position=position) +
      scale_fill_gradientn(colours=coloroption,limits=c(lower,upper))
  } else{
    c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val), position=position) +
      scale_fill_gradientn(colours=coloroption,limits=c(lower,upper),oob=squish)
  }
  
  # if (disp_nan) {
  #   c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val,color=sig), position=position) +
  #     scale_fill_gradientn(colours=coloroption,limits=c(lower,upper))+
  #     scale_color_continuous(na.value="transparent")
  # } else{
  #   c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val,color=sig), position=position) +
  #     scale_fill_gradientn(colours=coloroption,limits=c(lower,upper),oob=squish)+
  #     scale_color_continuous(low='red',high='red',na.value="transparent")
  # }
  # 
  
  
  if(!is.null(outpath)){
    pdf(paste0(outpath,'_cortex.pdf'))
    plot(c)
    dev.off()
  }
  return(c)
}
