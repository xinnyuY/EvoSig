

#' Plot signature matrix
#' @name sig_plot
#' @param sig signature matrix
#' @return mydata mutation data frame
#' @export
#' @importFrom reshape2 melt
#' @importFrom grid grid.draw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>% set_colnames
#' @import dplyr
sig_plot <- function(sig,tag=NULL,col="NA"){
  
  xx <- as.data.frame(t(apply(sig,2,function(x) x/sum(x)))) %>%
    set_colnames(1:ncol(.)) %>%
    mutate(signature=paste0("Signature ",1:nrow(.))) %>%
    melt(.,id=c("signature"))
  
  if (col=="NA"){
    if (ncol(sig)<=8) fills <- brewer.pal(8, "Set3")[c(1,3:8,2)] else
      fills <- brewer.pal(ncol(sig), "Set3")[c(1,3:8,2,9:ncol(sig))]
  } else {fills <- col}
  
  p1 <- ggplot(xx,aes(y=value,x=variable)) + geom_bar(aes(fill=signature),stat='identity') + scale_fill_manual(values = fills)+ theme_grey()+
    #ggtitle(paste0("rank = ",rank,", cancertype = ",type, ", MatrixType = ",MatType )) + 
    theme(legend.title = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(color= "white",size=10),
          panel.grid.minor.y = element_blank(),
          axis.title.x = element_text(color = "grey20", size = 8),
          axis.text.x = element_text(color = "grey20", size = 6),
          axis.text.y = element_text(color = "grey20", size = 6),
          plot.title = element_text(size= 10)) +
    facet_grid(cols = vars(signature))+ scale_x_discrete(breaks=c("1","50","100") ,labels=c("0", "0.5", "1"))+xlab("Cancer Cell Fraction") + ylab("") +
    labs(title="Consensus Signature of evolutionary dynamics",tag=tag)
  
  g1 <- ggplot_gtable(ggplot_build(p1))
  
  strip_both <- which(grepl('strip-', g1$layout$name))
  
  k <- 1
  for (i in strip_both) {
    j1 <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
    g1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
    #g1$layout$clip[j3] <- "off"
    k <- k+1
  }
  grid.draw(g1)
  return(g1)
}

#' Perform nmf for a cancer type
#' @name nmf_sig_plot_type
#' @param type cancer type
#' @param MatType matrix type
#' @param input_folder ccf file path
#' @param output output file path
#' @param rank rank
#' @return save nmf results in output folder and plot signature for this type
#' @import NMF
#' @importFrom magrittr %>% set_colnames
#' @import dplyr
nmf_sig_plot_type <- function(type,MatType="fraction",input_folder,output,rank,nrun){
  
  type_path <- paste0(input_folder,type,"/")
 
  if (MatType=="fraction") {
    file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfFractionMatrix_",dir(type_path))])
  } 
  
  if (MatType=="count"){
    file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfCountMatrix_",dir(type_path))])
  }

  if (!dir.exists(paste0(output,type))) {
    dir.create(paste0(output,type)) 
  }
  
  load(file=file_path)

  #format rank summary file
  if (exists("ccfFractionMatrix")) ccfMat <- ccfFractionMatrix
  if (exists("ccfCountMatrix")) ccfMat <- ccfCountMatrix
  
  n_sample <- ncol(ccfMat)
  
  ccf <- t(apply(ccfMat[1:100,],1,as.numeric))
  
  #preprocess for rows with all 0
  index_p <- which(rowSums(ccf)>0)
  index_n <- which(!rowSums(ccf)>0)
  ccf<- ccf[which(rowSums(ccf)>0),]
  
  #run NMF
  res <- nmf(ccf,rank,nrun=nrun,.opt='vp4')
  
  sig <- as.data.frame(matrix(0,nrow=length(index_p)+length(index_n),ncol=ncol(res@fit@W)))
  
  sig[c(index_p),] <- as.data.frame(res@fit@W) %>% set_colnames(paste0("sig_",1:ncol(.)))
  
  expo <- as.matrix(res@fit@H)
  
  #output sig and expo
  expo <- as.data.frame(t(expo)) 
  colnames(expo)[1:rank] <- paste0("sig_",1:rank)
  
  save(expo,file=paste0(output,type,'/',type,'_',n_sample,"_expo_",Sys.Date(),".RData"))
  save(sig,file=paste0(output,type,'/',type,'_',n_sample,"_sig_",Sys.Date(),".RData"))
  save(res,file=paste0(output,type,'/',type,'_',n_sample,"_res_",Sys.Date(),".RData"))
}

#' Perform nmf for multiple cancer types
#' @name nmf_sig_all_plot
#' @param input_folder ccf file path
#' @param output output file path
#' @param rank_summary rank summary file
#' @return save nmf results in output folder and plot signature for all cancer types
#' @export
#' @import dplyr
#' @import NMF
nmf_sig_all_plot <- function(input_folder,output,rank_summary,MatType="fraction",nrun=100) {
  
  if (!dir.exists(output)) dir.create(output)
  
  lapply(1:nrow(rank_suammry),function(x) nmf_sig_plot_type(rank_summary[x,1],input_folder=input_folder,output=output,rank=rank_summary[x,2],MatType=MatType,nrun=nrun))
}
