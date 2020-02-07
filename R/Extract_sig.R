
#' Perform NMF with consensus signature
#' @param ccfMat ccf matrix for all samples
#' @param consensus_sig consensus signature matrix
#' @param samplelist input sample list
#' @return exposure
Extract_sig <- function(ccfMat,consensus_sig,samplelist){
    
    n_sig <- nrow(consensus_sig)
    
    ccfMat[is.na(ccfMat)] <- 0
    
    EvoDynamics_exposure <- as.data.frame(t(YAPSA::LCD(ccfMat,t(consensus_sig[,2:101])))) %>%
      set_colnames(paste0("Evo_sig_",1:n_sig)) %>%
      mutate(sample=samplelist) %>%
      file_format(n_sig+1)
      
    save(EvoDynamics_exposure,file=paste0(Exposure_folder,"evoDynamic_exposure.RData"))
    
    EvoDynamics_exposure
}

  