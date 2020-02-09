## Function 1 - Heatmap for all the variables


#' Basic correlation function for measures with each cancer type
#' @name corr_spearman_type
#' @param xx exposure file merged with measures of interest
#' @param n_sig number of evolutionary dynamics signatures
#' @param varcol input dataframe containing all samples with corresponding types
#' @param output specify ouput correlation value, p value or sample size
#' @return depends on output choice
#' @importFrom psych corr.test
corr_spearman_type <- function(xx,n_sig,varcol,output){

  col <- function(data,n_sig,varcol){
    return(corr.test(data[,1:n_sig],data[,varcol],method = "spearman",use = "pairwise",adjust="holm"))
  }
  
  r_commands <- paste0("r <- data.frame(",paste0("cor=col(xx,",n_sig,",",varcol,")[['r']]",collapse=","),")")
  p_commands <- paste0("p <- data.frame(",paste0("p=col(xx,",n_sig,",",varcol,")[['p']]",collapse=","),")")
  n_commands <- paste0("n <- data.frame(",paste0("n=col(xx,",n_sig,",",varcol,")[['n']]",collapse=","),")")
  eval(parse(text=r_commands))
  eval(parse(text=p_commands))
  eval(parse(text=n_commands))
  if (output=="r") return(r)
  if (output=="p") return(p)
  if (output=="n") return(n)
}

#' Pre-processing before computing correlation and set the output formate
#' @name cor_evoSig
#' @param file exposure file merged with measures of interest
#' @param varcol number of evolutionary dynamics signatures
#' @param n_sig input dataframe containing all samples with corresponding types
#' @param output specify ouput correlation value, p value or sample size
#' @return depends on output choice
cor_evoSig <- function(file,varcol,n_sig,output="r")   {
  
  # delete NAs
  data <- file[which(!is.na(file[,varcol])),]
  
  if (nrow(data)<=2) {
    data <- NULL
  }  else {
    data <- as.data.frame(t(corr_spearman_type(data,n_sig=n_sig,varcol=varcol,output=output)))
    if (ncol(data)==n_sig) colnames(data) <- paste0("Signature ",1:n_sig)
    if (ncol(data)==1) colnames(data) <- "n_sample"
  }
  return(data)
}

#' Correlation of a single measure with evolutionary signatures across cancer types 
#' @name heatmap_plot_for_single_measure_across_types
#' @param file Ccube folder
#' @param varcol Output folder
#' @param n_sig number of signatures
#' @param keep_pl0.05 whether to keep correaltion with p>0.05
#' @param minsample the minimun samples
#' @param Type Specify whether to compute correlation for each cancer type
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
#' @importFrom plyr ddply "."
#' @import dplyr
#' @importFrom reshape2 melt
heatmap_plot_for_single_measure_across_types <- function(file,varcol,n_sig,keep_pl0.05 = FALSE,minsample=minsample,Type){     
 
  if (Type) {
    xx_r <- ddply(file, .(cancertype), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig) %>% melt(id="cancertype")
    xx_p <- ddply(file, .(cancertype), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig,output="p") %>% melt(id="cancertype")
    cor_n_type <- ddply(file, .(cancertype), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig,output="n")
    xx_cor <- left_join(xx_r,xx_p,by=c("cancertype","variable")) %>% left_join(cor_n_type,by="cancertype")
  } else {
    xx_r <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig) %>% mutate(cancertype="All") %>% melt(id="cancertype")
    xx_p <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig,output="p") %>% mutate(cancertype="All") %>% melt(id="cancertype")
    cor_n_type <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig,output="n") %>% mutate(cancertype="All")
    xx_cor <- left_join(xx_r,xx_p,by=c("cancertype","variable")) %>% left_join(cor_n_type,by="cancertype")
  }
  
  colnames(xx_cor)[3:5] <- c("r","p","n_sample")
  
  # choose whether to keep value with p>0.05
  if (keep_pl0.05 == FALSE & any(xx_cor$p>0.05)) 
    xx_cor[which(xx_cor$p>0.05),]$r <- 0
  
  # set min tumour samples, delete types with NA and sample size low than the threshold
  if (mean(is.na(xx_cor$r))>0) 
    xx_cor <- xx_cor[-which(is.na(xx_cor$r)),] 
  
  if (length(which(xx_cor$n_sample < minsample))>0) 
    xx_cor <-xx_cor[-which(xx_cor$n_sample < minsample),] 
  
  cor_n_type <- unique(xx_cor[,c("cancertype","n_sample")]) 
 
  varlabel <- colnames(file)[varcol]
  
  ###Set the background color to white / makes the background white
  par(bg = 'white')
  
  coul = colorRampPalette(brewer.pal(8, "RdBu"))(256)
  
  p1 <- ggplot(data=xx_cor,mapping=aes(x=cancertype,y=variable,fill=r)) + geom_tile() + #ggplot need each value be a rows in the new datasets
    scale_fill_gradient2(low=coul[256],mid="#ffffff",high=coul[1],midpoint = 0,name="Spearman Correlation")+
    geom_text(aes(x=cancertype,y=variable,label=round(r,2)),size=3,col="#ffffff")+
    labs(subtitle = paste0(varlabel))+
    annotate('text',x=1:nrow(cor_n_type),y=rep(n_sig+0.7,each=nrow(cor_n_type)),label=cor_n_type[,2],size=3)+
    ylab("Evolutionary Signatutres")+xlab("TCGA cancertype")+
    theme(axis.text.x = element_text(angle = 30, hjust = 1,size=10),
          axis.text.y=element_text(size=10),
          legend.position = "right",legend.text = element_text(size=6),
          legend.title = element_text(size=10,face="bold"),
          axis.title.y = element_text(size=12,face="bold"),
          plot.margin = unit(c(3,1,1,1),"lines"))+ 
    coord_cartesian(ylim=c(1,n_sig+0.4))
  
  return(list(p1,xx_cor))
}


#' Plot heatmap for correlation between measures and evolutionary signature
#' @name heatmap_plot_multi_variable
#' @param file exposure file with measures of interest
#' @param varcol index of columns of measures
#' @param n_sig number of signatures
#' @param Type Specify compute correlation for each cancer type or for all
#' @param minsample minimun samples required for computing correlation
#' @return return heatmap for correlation between measures and evolutionary signature
#' @export
heatmap_plot_multi_variable <- function(file,varcol,n_sig,Type=TRUE,minsample=10,output=NA){
  
  for (i in 1:length(varcol)){
    
    varlabel <- colnames(file)[varcol[i]]
    commands1 <- paste0("p_",i,"<- heatmap_plot_for_single_measure_across_types(file,varcol=",varcol[i],",n_sig=",n_sig,",Type=",Type,",minsample=",minsample,")[[1]]")
    commands2 <- paste0("table_",i,"<- heatmap_plot_for_single_measure_across_types(file,varcol=",varcol[i],",n_sig=",n_sig,",Type=",Type,",minsample=",minsample,")[[2]]")
    commands3 <- paste0("table_",i,"$varlable <- '",varlabel,"'")
    eval(parse(text=commands1))
    eval(parse(text=commands2))
    eval(parse(text=commands3))
    
    if (i==1) {
      total_table <- table_1} 
    else {
      commands <- paste0("total_table <- rbind(total_table,table_",i,")")
      eval(parse(text=commands))
    }
    
    if (!is.na(output)){
      multi_dir_create(output)
      pdf(paste0(output,"correlation_",varlabel,".pdf"), height = 10, width = 30)
      eval(parse(text=paste0("print(p_",i,")")))
      dev.off()
    }
  }
  return(total_table)
}

#' Correlation heatmap plotting
#' @name Cor_heatmap
#' @param file input exposure file with measures
#' @param output output folder path
#' @param n_sig number of signatures
#' @param width width of output plot
#' @param height height of output plot
#' @param facet specify facet 
#' @param minsample minimun samples required for computing correlation
#' @return return heatmap for correlation between measures and evolutionary signature
#' @export
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grid grid.draw
#' @import ggplot2
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @import dplyr
Cor_heatmap <- function(file,var_index,output="NA",n_sig,minsample,width = 12,height = 30,facet="signature") {
    
    file[,var_index] <- apply(file[,var_index],2,as.numeric)
    file[,(ncol(file)-1):ncol(file)] <- apply(file[,(ncol(file)-1):ncol(file)],2,as.character)
    
    p_type <- heatmap_plot_multi_variable(file,varcol=var_index,n_sig=n_sig,minsample=minsample,output=output)   
    
    p_all <- heatmap_plot_multi_variable(file,varcol=var_index,n_sig=n_sig,Type=F)    
    
    p_ave <- p_type %>%
      group_by(variable,varlable) %>%
      dplyr::summarize(value=mean(r))

    coul = colorRampPalette(brewer.pal(8, "RdBu"))(256)
    par(bg = 'white')
    
    if (facet=="cancertype") p1 <- ggplot(p_type,aes(x=variable,y=varlable,fill=r))+geom_tile()+facet_grid(cols = vars(cancertype))+
                                    geom_text(aes(x=variable,y=varlable,label=round(r,2)),size=3,col="#ffffff")
    if (facet=="signature") p1 <- ggplot(p_type,aes(x=cancertype,y=varlable,fill=r))+geom_tile()+facet_grid(rows = vars(variable))+
                                    geom_text(aes(x=cancertype,y=varlable,label=round(r,2)),size=3,col="#ffffff")
    if (facet=="measure") p1 <- ggplot(p_type,aes(x=variable,y=cancertype,fill=r))+geom_tile()+facet_grid(cols = vars(varlable))+
                                    geom_text(aes(x=variable,y=cancertype,label=round(r,2)),size=3,col="#ffffff")
    
    p1 <- p1 + scale_fill_gradient2(low=coul[256],mid="#ffffff",high=coul[1],midpoint = 0,name="Spearman Correlation")+
      #labs(subtitle = paste0(varlabel))+
      #annotate('text',x=1:nrow(cor_n_type),y=rep(n_sig+0.7,each=nrow(cor_n_type)),label=cor_n_type[,2],size=3)+
      ylab("")+xlab("")+ theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10),
            axis.text.y=element_text(size=10),
            legend.position = "right",legend.text = element_text(size=6),
            legend.title = element_text(size=10,face="bold"),
            axis.title.y = element_text(size=12,face="bold"),
            plot.margin = unit(c(3,1,1,1),"lines"),panel.background = element_rect(fill = "white",colour = "grey",size = 0.5, linetype = "solid"),
            panel.grid.minor.y = element_blank(),
            strip.text.y = element_text(color="white",size=12))
    #coord_cartesian(ylim=c(1,n_sig+0.4))
     
    g1 <- ggplotGrob(p1)
    my_color <- brewer.pal(8, "Set3")[c(1,3:8,2)]
    ### change the color of facet strips
    
    strip_both <- which(grepl('strip-', g1$layout$name))
    k <- 1
    
    for (i in strip_both) {
      j1 <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
      g1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- my_color[k]
      k <- k+1
    }
    
    if(!is.na(output)) {
      
      multi_dir_create(output)
      
      save(p_type,file=paste0(output,"cor_table.RData"))
      grid.draw(g1)
      ggsave(filename = paste0(output,"Integrated_correlation_facet(",facet,").pdf"),width=width,height=height)
    }
}  



 

       