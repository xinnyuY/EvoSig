hc_consensus <-
function(combine_sig,cluster,output,ccfMatrix,distance="euclidean")  {
  library("RColorBrewer")
  library("pheatmap")
  library("cowplot")
  #library(factoextra)
  
  upper = quantile(combine_sig,0.95)
  breaksList = seq(0, upper, by = 0.01)
  col <- colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(length(breaksList))
    
  ## set `cutree_cols` based on suggested cluster number 
  out <- pheatmap(combine_sig, cutree_cols = cluster, fontsize_col = 5,fontsize_row = 0.4,color = col, breaks = breaksList,clustering_distance_cols=distance, cluster_rows=F,filename=paste0(output,"/",distance,"_",cluster,"_hc_heatmap.pdf"),clustering_method = "ward.D2")
 
  sig_label <- as.data.frame(cutree(out$tree_col,k=cluster)) %>%
      set_colnames("cluster") %>%
      mutate(sig=rownames(.))
    
  ## Compute consensu signatures for each cluster
  consensus_sig <- as.data.frame(t(combine_sig)) %>%
      mutate(sig=rownames(.)) %>%
      left_join(.,sig_label,by="sig") %>%
      group_by(cluster) %>%
      mutate(sig=NULL) %>%
      summarise_all(mean) 
  
  save(consensus_sig,file=paste0(output,"/",distance,"_",cluster,"_consensus_sig.RData"))
  consensus_sig <- apply(t(consensus_sig[,2:101]),2,as.numeric)
  
  p1 <- plot_grid(sig_plot(consensus_sig))
  save_plot(paste0(output,"/",distance,"_",cluster,"_consensus_sig.pdf"),p1,base_asp = cluster)
  
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("YAPSA",quietly = TRUE)) BiocManager::install("YAPSA")
  library(YAPSA)
 
  load(paste0(HC_folder,dir(HC_folder)[grep("samples",dir(HC_folder))]))
  
  ccfMatrix[is.na(ccfMatrix)] = 0
  idx_random <- sample(1:ncol(ccfMatrix))
  ccfMatrix <- ccfMatrix[,idx_random]
  samplelist <- post_summary$samplename[idx_random]
  
  write.csv(samplelist,file=paste0(HC_folder,distance,"_",cluster,"_sample_list_random.csv"))
  
  exposure <- YAPSA::LCD(ccfMatrix,consensus_sig) 
  table(rowSums(ccfMatrix)>0)
  exposure <- as.data.frame(t(exposure)) %>%
    set_colnames(paste0("sig_",1:ncol(.)))
  save(exposure,file=paste0(output,"/",distance,"_",cluster,"_lcd_exposure.RData"))
  
  return(consensus_sig)
}
