library(ggseg) # We are using ggseg 1.5.5
library(ggsegExtra) # We are using ggseg 1.5.5
library(RColorBrewer)
library(scales)
library(viridis)
library(colormap)

#brainplot = function(values,nm.labels,outpath=NULL,lower=NULL,upper=NULL,coloroption=NULL,position=NULL,disp_nan=FALSE){
brainplot = function(values,nm.labels,outpath=NULL,lower=NULL,upper=NULL,coloroption=NULL,position=NULL,disp_nan=FALSE,sig=NULL){
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

  # Convert Labels 
  #label.idx = match(nm.labels,nm)
  nm.labels = gsub('  ','_',nm)
  nm.labels = paste(nm.labels, "_ROI", sep="") # uncomment for orig
  nm.labels[1:16] = c("Thalamus","Caudate","Putamen","Pallidum","Hippocampus","Amygdala","Accumbens","Diencephalon",
                   "Thalamus","Caudate","Putamen","Pallidum","Hippocampus","Amygdala","Accumbens","Diencephalon")
  #nm.labels = labels[label.idx]
  
  ###
  
  sig = if(!is.null(sig)){sig[which(!is.na(values))]}
  cortex.values = values[which(!is.na(values))]
  cortex.labels = nm.labels[which(!is.na(values))]
  #cortex.labels = glassersubcortex$label[match(toupper(nm.labels[!is.na(values)]), toupper(glassersubcortex$label))]
  
  if(!is.null(sig)){
    cortex = data.frame(
      label = cortex.labels,
      val = cortex.values,
      stringsAsFactors = FALSE,
      sig=sig)
    
    #c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val,color=ifelse(sig == 1, 1, NA)), position=position) +
    c <- ggseg(.data=cortex, atlas="glassersub", mapping=aes(fill=val,color=ifelse(sig == 1, 1, NA)), position=position) +
      scale_color_continuous(low='red',high='red',na.value="transparent") +
      if(disp_nan){
        scale_fill_gradientn(colours=coloroption,limits=c(lower,upper),oob=squish)
      } else{
        scale_fill_gradientn(colours=coloroption,limits=c(lower,upper))
      }
    
  } else{
    cortex = data.frame(
      label = cortex.labels,
      val = cortex.values,
      stringsAsFactors = FALSE)
    
    c <- ggseg(.data=cortex, atlas="glassersub", mapping=aes(fill=val), position=position) +
      if(disp_nan){
        scale_fill_gradientn(colours=coloroption,limits=c(lower,upper),oob=squish)
      } else{
        scale_fill_gradientn(colours=coloroption,limits=c(lower,upper))
      }
  }
  # if (disp_nan) {
  #   c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val), position=position) +
  #     scale_fill_gradientn(colours=coloroption,limits=c(lower,upper))
  # } else{
  #   c <- ggseg(.data=cortex, atlas="glassersubcortex", mapping=aes(fill=val), position=position) +
  #     scale_fill_gradientn(colours=coloroption,limits=c(lower,upper),oob=squish)
  # }
  
  
  
  
  
  if(!is.null(outpath)){
    pdf(paste0(outpath,'_cortex.pdf'))
    plot(c)
    dev.off()
  }
  return(c)
}
