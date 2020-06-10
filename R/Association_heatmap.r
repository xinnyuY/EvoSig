## Function 1 - Heatmap for all the variables

#' Basic correlation function for each cell
#' @name corr_spearman_type
#' @param xx exposure file merged with measures of interest
#' @param n_sig number of evolutionary dynamics signatures
#' @param varcol input dataframe containing all samples with corresponding types
#' @param output specify ouput correlation value, p value or sample size
#' @return depends on output choice
corr_spearman_type <- function(xx,n_sig,varcol,output){
  col <- function(data,n_sig,varcol){
    return(corr.test(data[,1:n_sig],data[,varcol],method = "spearman",use = "pairwise",adjust="holm"))
  }
  
  if (output=="r") {
    commands <- paste0("data.frame(",paste0("cor=col(xx,",n_sig,",",varcol,")[['r']]",collapse=","),")")
  }
  
  if (output=="p") {
    commands <- paste0("data.frame(",paste0("p=col(xx,",n_sig,",",varcol,")[['p']]",collapse=","),")")
  }
  
  if (output=="n") {
    commands <- paste0("data.frame(",paste0("n=col(xx,",n_sig,",",varcol,")[['n']]",collapse=","),")")
  }
  
  eval(parse(text=commands))
}

#' Pre-processing before computing correlation and set the output formate
#' @name cor_evoSig
#' @param file exposure file merged with measures of interest
#' @param varcol number of evolutionary dynamics signatures
#' @param n_sig input dataframe containing all samples with corresponding types
#' @param output specify ouput correlation value, p value or sample size
#' @return depends on output choice
#' #@import dplyr
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
#' @name heatmap_plot_for_single_measure_across_groups
#' @param file input file
#' @param varcol columns index of variable of interest 
#' @param n_sig number of signatures
#' @param keep_pl0.05 whether to keep correaltion with p>0.05
#' @param minsample the minimun samples
#' @param Type Specify whether to compute correlation for each cancer type
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
#' @importFrom plyr ddply "."
#' @import dplyr
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette
heatmap_plot_for_single_measure_across_groups <- function(file,varcol,n_sig,group=NA,minsample=10){   
  
  if (!is.na(group)) {
    command_r <- paste0("xx_r <- ddply(file, .(",group,"), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig) %>% melt(id='",group,"')")
    command_p <- paste0("xx_p <- ddply(file, .(",group,"), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig,output='p') %>% melt(id='",group,"')")
    command_n_type <- paste0("cor_n_type <- ddply(file, .(",group,"), .fun =cor_evoSig,varcol=varcol,n_sig=n_sig,output='n') ")
  } else {
    command_r <- paste0("xx_r <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig) %>% mutate(",group,"='All') %>% melt(id='",group,"')")
    command_p <- paste0("xx_p <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig,output='p') %>% mutate(",group,"='All') %>% melt(id='",group,"')")
    command_n_type <- paste0("cor_n_type <- file %>% cor_evoSig(varcol=varcol,n_sig=n_sig,output='n') %>% mutate(",group,"='All') ")
  }

  command <- paste0("xx_cor <- left_join(xx_r,xx_p,by=c('",group,"','variable')) %>% left_join(cor_n_type,by='",group,"')")
  
  eval(parse(text=command_r));eval(parse(text=command_p));eval(parse(text=command_n_type));eval(parse(text=command))
  
  colnames(xx_cor)[3:5] <- c("r","p","n_sample")
  
  xx_cor[which(xx_cor$p>0.05),]$r <- 0
  
  if (length(which(xx_cor$n_sample < minsample))>0) 
    xx_cor <-xx_cor[-which(xx_cor$n_sample < minsample),] 
  
  return(xx_cor)
}

#' Obtain correlation table with r,p,n between signature exposures and measures of interest (across group)
#' @name heatmap_plot_multi_groups
#' @param file exposure file with measures of interest
#' @param varcol index of columns of measures
#' @param n_sig number of signatures
#' @param group specify group variables
#' @return return correlation table with r,p,n between signature exposures and measures of interest (across group)
#' @export
#heatmap_plot_multi_variable_groups(file,var_idx=50:52,n_sig=4,group="cms_label")
cor_table_multi_variable_groups <- function(file,var_idx,n_sig,group=NA){
  
  cor_var_group <- function(file,var_idx,n_sig,group=NA) {
    varlabel <- colnames(file)[var_idx]
    commands2 <- paste0("table <- heatmap_plot_for_single_measure_across_groups(file,varcol=",var_idx,",n_sig=",n_sig,",group='",group,"')")
    commands3 <- paste0("table$varlable <- '",varlabel,"'")
    eval(parse(text=commands2));eval(parse(text=commands3))
    return(table)
  }
  command <- paste0("lapply(var_idx, function(x) cor_var_group(file,var_idx=x,n_sig=",n_sig,",group='",group,"')) %>% bind_rows()") 
  eval(parse(text=command))
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
#' @importFrom gridExtra grid.arrange
Cor_heatmap <- function(file,var_index,immune_measures_mapping=immune_measures_mapping,output="NA",n_sig,minsample,width = 30,height = 69,immune=TRUE,type=TRUE,var=TRUE,r_filter=0.2) {
    
    file[,var_index] <- apply(file[,var_index],2,as.numeric)
    file$cancertype <- as.character(file$cancertype)
    colnames(file)[1:n_sig] = paste0("Signature ",1:n_sig)
    
    colnames(immune_measures_mapping)[1] <- "varlable"
    
    p_type <- heatmap_plot_multi_variable(file,varcol=var_index,n_sig=n_sig,minsample=minsample) %>%
      left_join(immune_measures_mapping,by="varlable")
    p_type <- p_type[!is.na(p_type$Cluster),]
    
    my_xlab=paste0(names(table(file$cancertype)),"\n(N=",table(file$cancertype),")",sep="")
    
    #p_all <- heatmap_plot_multi_variable(file,varcol=var_index,n_sig=n_sig,Type=F)    
    
    # p_ave <- p_type %>%
    #   group_by(variable,varlable) %>%
    #   dplyr::summarize(value=mean(r))
  
    coul = colorRampPalette(brewer.pal(8, "RdBu"))(256)
    par(bg = 'white')
    sig_color <- brewer.pal(8, "Set3")[c(1,3:8,2)]
    facet_color <- c('#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')
    
   
    multi_dir_create(output)
    
    if (immune==TRUE) {
      for (i in 1:n_sig){
        p_type_filter <- subset(p_type,variable==paste0("Signature ",i))
        
        
        p <- ggplot(p_type_filter,aes(x=cancertype,y=varlable,fill=r))+geom_tile()+ facet_grid(rows = vars(Cluster),switch="y",scale="free_y",space="free_y") +
          scale_fill_gradient2(low=coul[256],mid="#ffffff",high=coul[1],midpoint = 0,name="Spearman Correlation")+
          geom_text(aes(x=cancertype,y=varlable,label=round(r,2)),size=6,col="#ffffff") +
          ylab("")+xlab("")+ theme(axis.text.x = element_text(size=17),
                                   axis.text.y=element_text(size=17),
                                   legend.position = "bottom",legend.text = element_text(size=10),
                                   legend.title = element_text(size=19,face="bold"),
                                   axis.title.y = element_blank(),
                                   axis.ticks.x = element_blank(),
                                   plot.margin = unit(c(0,2.5,0,2.5),"cm"),panel.background = element_rect(fill = "white",colour = "white",size=3,linetype = "solid"),
                                   panel.grid.minor.y = element_blank(),
                                   strip.text.y = element_text(color="black",size=14,margin=margin(0,0.4,0,0.4,"cm")))+
          scale_x_discrete(labels=my_xlab)
        
        
        p1 <- ggplotGrob(p)
        
        ### change the color of facet strips
        strip_both <- which(grepl('strip-', p1$layout$name))
        k <- 1
        
        for (j in strip_both) {
          j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
          p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- facet_color[k]
          k <- k+1
        }
        
        
        rect <- data.frame(y=1,color=sig_color[i],text=paste0("Signature ",i))
        
        g <- ggplot(rect,aes(x=0,y,fill=color,label=text),)+
          geom_tile(width=0.1,height=1)+
          geom_text(color="white",size=8)+ 
          scale_fill_identity((guide="none"))+
          theme_void()+theme(plot.margin=unit(c(1,0,0.3,0), "cm"))
        
        plot_grid(g)
        
        pdf(file=paste0(output,"immune_correlation_signature_",i,".pdf"),height=25,width=28)
        
        grid.arrange(g,p1,ncol =1,heights=c(1,20))
        dev.off()
      }}
      
      # plot scatterplot for correlation higher than `r_filter``
      cor_g <- subset(p_type,abs(r)>=r_filter) %>%
        mutate(p=round(p,2))
      
      if (type == TRUE) {
      typelist <- unique(cor_g$cancertype)
      
      for (i in 1:length(typelist)){
        type = typelist[i]
        cor <- subset(cor_g,cancertype==type)
        cor$variable <- as.character(cor$variable)
        cor <- cor[order(cor$varlable),]
        
        immune_expo_type <- subset(file,cancertype==type)
      
       
        for (j in 1:nrow(cor)){
          x = cor$variable[j]
          y = cor$varlable[j]
          r = cor$r[j]
          p = cor$p[j]
          command1 = paste0("n = length(which(!is.na(immune_expo_type$`",y,"`)==TRUE))")
          eval(parse(text=command1))
          
          # add scale_y_sqrt() for non-negative and large values
          command2 = paste0("sqrt <- (min(immune_expo_type$`",y,"`,na.rm = T)>=0 & max(immune_expo_type$`",y,"`,na.rm = T)>=10)")
          eval(parse(text=command2))
          
          command3 <- paste0("p <-ggplot(immune_expo_type)+geom_point(aes(x=`",x,"`,y=`",y,"`))+scale_x_sqrt()+xlab('Exposure (count/MB)')+labs(title= '",x,"',subtitle='spearman correaltion = ",round(r,3),", p = ",p,"\nn=",n,"')+theme(axis.text = element_text(size=6),plot.title=element_text(size=12,face='bold'),plot.caption=element_text(size=4,face = 'italic'))")
          eval(parse(text=command3))
          
          if (sqrt) p <- p + scale_y_sqrt()
          
          command4 <- paste0("p",j,"<- ggplotGrob(p)")
          eval(parse(text=command4))
        }
        
        command5 <- paste0("g",i,"<- grid.arrange(",paste0('p',1:nrow(cor),collapse=","),",nrow=",ceiling(nrow(cor)/4),",top = textGrob('- - - - - - - - ",type," - - - - - - - -',gp=gpar(fontsize=20,font=2)))")
        eval(parse(text=command5))
        
        pdf(file=paste0(output,"immune_scatter_",type,"_",Sys.Date(),".pdf"),width=12,height=(3*ceiling(nrow(cor)/4)+1))
        command5 <- paste0("grid.draw(g",i,")")
        eval(parse(text=command5))
        dev.off()
        
      }}
    
    # # plot scatterplot for correlation (signature~cancertype)
    if (var==TRUE) {
    var_list <- unique(cor_g$varlable)
   
    for (i in 1:length(var_list)){
      var = var_list[i]
      cor <- subset(cor_g,varlable==var) %>%
        mutate(variable=as.character(variable),cancertype=as.character(cancertype))%>%
        arrange(cancertype)
     
      for (j in 1:nrow(cor)){
        x = cor$cancertype[j]
        y = cor$variable[j]
        r = cor$r[j]
        p_val = cor$p[j]
        
        immune_expo_type <- subset(file,cancertype==x)
        
        command1 = paste0("n = length(which(!is.na(immune_expo_type$`",var,"`)==TRUE))")
        eval(parse(text=command1))
        
        # add scale_y_sqrt() for non-negative and large values
        command2 = paste0("sqrt_x <-  max(immune_expo_type$`",y,"`,na.rm = T)>=10")
        command3 = paste0("sqrt_y <- (min(immune_expo_type$`",var,"`,na.rm = T)>=0 & max(immune_expo_type$`",var,"`,na.rm = T)>=10)")
        eval(parse(text=command2))
        eval(parse(text=command3))
        
        command4 <- paste0("p <-ggplot(immune_expo_type)+geom_point(aes(x=`",y,"`,y=`",var,"`))+ylab('",var,"')+labs(title= '",x,"',subtitle='spearman correaltion = ",round(r,3),"\np = ",p_val,"\nn=",n,"')+theme(axis.text = element_text(size=6),plot.title=element_text(size=12,face='bold'),plot.caption=element_text(size=4,face = 'italic'))")
        eval(parse(text=command4 ))
        
        if (sqrt_x) p <- p + scale_x_sqrt()
        if (sqrt_y) p <- p + scale_y_sqrt()
        
        command5 <- paste0("p",j,"<- ggplotGrob(p)")
        eval(parse(text=command5))
      }
     
      command5 <- paste0("g",i,"<- grid.arrange(",paste0('p',1:nrow(cor),collapse=","),",nrow=",ceiling(nrow(cor)/4),",top = textGrob('- - - - - - - - ",var," - - - - - - - -',gp=gpar(fontsize=20,font=2)))")
      eval(parse(text=command5))
      if (nrow(cor)>=4){
      pdf(file=paste0(output,"immune_scatter_",var,"_",Sys.Date(),".pdf"),width=12,height=(3*ceiling(nrow(cor)/4)+1))
      } else {
        pdf(file=paste0(output,"immune_scatter_",var,"_",Sys.Date(),".pdf"),width=3*nrow(cor),height=(3*ceiling(nrow(cor)/4)+1))  
      }
      command5 <- paste0("grid.draw(g",i,")")
      eval(parse(text=command5))
      dev.off()
      
      
    }
  }

}  

#' Correlation heatmap plotting
#' @name Cor_heatmap_colon
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
#' @importFrom gridExtra grid.arrange
Cor_heatmap_colon <- function(file,var_index,immune_measures_mapping=immune_measures_mapping,output="NA",n_sig,immune=TRUE,var=TRUE,group="cancertype") {
  
  file[,var_index] <- apply(file[,var_index],2,as.numeric)
  file[,group] <- as.character(file[,group])
  colnames(file)[1:n_sig] = paste0("Signature ",1:n_sig)
  
  colnames(immune_measures_mapping)[1] <- "varlable"
  
  cor_table <- cor_table_multi_variable_groups(file,var_idx=var_index,n_sig=n_sig,group=group) %>%
    left_join(immune_measures_mapping,by="varlable") %>%
    drop_na(Cluster) 
  
  cor_table[is.na(cor_table$r),]$r = 0
  
  group_size_xlab = paste0(names(table(file[,group])),"\n(N=",table(file[,group]),")",sep="")
  
  coul = colorRampPalette(brewer.pal(8, "RdBu"))(256)
  par(bg = 'white')
  sig_color <- brewer.pal(8, "Set3")[c(1,3:8,2)]
  facet_color <- c('#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2')
  

  multi_dir_create(output)
  
  eval(parse(text=paste0("p <- ggplot(cor_table,aes(x=",group,",y=varlable,fill=r))")))    
  p <- p +geom_tile()+ facet_grid(rows = vars(Cluster),cols=vars(variable),switch="y",scale="free_y",space="free_y") +
        scale_fill_gradient2(low=coul[256],mid="#ffffff",high=coul[1],midpoint = 0,name="Spearman Correlation")
  eval(parse(text=paste0("p <- p + geom_text(aes(x=",group,",y=varlable,label=round(r,2)),size=3,col='#ffffff')")))   
  p <- p+labs(x=" ",y=" ",caption = "*Correlation value only shown with adjusted p value <0.05 and available samples>10")+ 
    theme(axis.ticks.x = element_blank(),
          panel.background = element_rect(fill = "white",colour = "white",size=3,linetype = "solid"),
          panel.grid.minor.y = element_blank())+ scale_x_discrete(labels=my_xlab)
      
      
  p1 <- ggplotGrob(p)
      
  ### change the color of facet strips
  strip_t <- which(grepl('strip-t', p1$layout$name))
  strip_l <- which(grepl('strip-l', p1$layout$name))
      
      k <- 1
      for (j in strip_t) {
        j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
        p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- sig_color[k]
        k <- k+1
      }
      
      
      k <- 1
      for (j in strip_l) {
        j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
        p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- facet_color[k]
        k <- k+1
      }
      
      grid.draw(p1)
    
}
  


       