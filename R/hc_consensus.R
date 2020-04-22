
#' Plotting for HC results 
#' @name hc_consensus ccfMat ccf matrix for all samples
#' @param combine_sig combined signatures for all cancer types
#' @param cluster clustering number
#' @param output output folder path
#' @param ccfMatrix ccf matrix for all samples
#' @param distance distance function
#' @return consensus_sig
#' @import pheatmap
#' @importFrom rlang .data
#' @importFrom cowplot save_plot
#' @import dplyr
#' @import crayon
#' @importFrom magrittr %>% set_colnames
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @export
hc_consensus_sig <- function(input_folder,output_folder,cancertype,min.nc,max.nc,method="ward.D2",index="hubert",distance="euclidean"){
  
  cat(red("Start run hierarchical clustering to extract consensus signature. \n"))
  
  combine_sig <- combine_sig_nmf(input_folder=input_folder,output_folder=output_folder,cancertype=cancertype) 
  
  if (index=="hubert") {
    pdf(paste0(output_folder,"hc_hubert_plot.pdf")) 
    res.nbclust  <- NbClust(t(combine_sig),min.nc = min.nc, max.nc = max.nc,method=method,index = index,distance=distance)
    dev.off()
    
    second_diff = res.nbclust$All.index[2:10]-res.nbclust$All.index[1:9]
    cluster = which(second_diff==max(second_diff))+min.nc
    print(paste0("The suggested number of clusters from the Hubert index was - ",cluster))
    cat(paste0("-> The Hubert index was saved as: ",paste0(output_folder,"hc_hubert_plot.pdf \n")))
    
    upper = quantile(as.matrix(combine_sig),0.95);breaksList = seq(0, upper, by = 0.01)
    col <- colorRampPalette(rev(brewer.pal(n = cluster, name = "RdYlBu")))(length(breaksList))
  }
  
  out <- pheatmap(combine_sig, cutree_cols = cluster, fontsize_col = 5,fontsize_row = 0.4,color = col, breaks = breaksList,clustering_distance_cols=distance,cluster_rows=F,filename=paste0(output_folder,"/",distance,"_",cluster,"_hc_heatmap_",Sys.Date(),".pdf"),clustering_method = method)
  
  sig_label <- as.data.frame(cutree(out$tree_col,k=cluster)) %>%
    set_colnames("cluster") %>%
    mutate(sig=rownames(.))
  
  ## Compute consensu signatures for each cluster
  consensus_sig <- as.data.frame(t(combine_sig)) %>%
    mutate(sig=rownames(.)) %>%
    left_join(sig_label,by="sig") %>%
    group_by(cluster) %>%
    mutate(sig=NULL) %>%
    summarise_all(mean) 
  
  save(consensus_sig,file=paste0(output_folder,"/",distance,"_",cluster,"_consensus_sig_",Sys.Date(),".RData"))
  
  consensus_sig <- apply(t(consensus_sig[,2:101]),2,as.numeric) 
  p1 <- plot_grid(sig_plot(consensus_sig))
  save_plot(paste0(output_folder,"/",distance,"_",cluster,"_consensus_sig_",Sys.Date(),".pdf"),p1,base_asp = cluster)
  
  cat(paste0("-> The consensus signature result was saved at: ",output_folder,"\n"))
  
  return(consensus_sig) 
}

#' sig_assignment 
#' @name hc_consensus ccfMat ccf matrix for all samples
#' @param signature signature
#' @param ccfMatrix ccfmatrix
#' @param output output folder path
#' @return exposure
#' @import dplyr
#' @import crayon
#' @importFrom rlang .data
#' @importFrom cowplot save_plot
#' @importFrom YAPSA LCD
#' @export
sig_assignment <- function(signature,ccfMatrix,output=NA){
  
  ccfMatrix[is.na(ccfMatrix)] = 0
  
  exposure <-LCD(ccfMatrix,signature) 
  exposure <- as.data.frame(t(exposure)) %>% set_colnames(paste0("sig_",1:ncol(.))) 
  
  if (!is.na(output)) save(exposure,file=paste0(output,"/lcd_exposure",Sys.Data,".RData"))
  
  return(exposure)
}
 
