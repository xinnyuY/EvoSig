
#' Perform NMF with consensus signature
#' @param ccfMat ccf matrix for all samples
#' @param consensus_sig consensus signature matrix(columns as signaguture,row as bins)
#' @param samplelist input sample list
#' @param output output folder path
#' @return exposure
#' @export
#' @importFrom YAPSA LCD
Extract_sig <- function(ccfMat,consensus_sig,output=NA){
    
    n_sig <- ncol(consensus_sig)
    
    Mat <- t(ccfMat[,1:100])
    Mat[is.na(Mat)] <- 0
    
    EvoDynamics_exposure <- as.data.frame(t(LCD(Mat,consensus_sig))) %>%
      set_colnames(paste0("Evo_sig_",1:n_sig)) %>%
      mutate(samplename=ccfMat[,101]) %>%
      file_format(n_sig+1)
    
    if (!is.na(output)) {
      
      if (!dir.exists(output)) {
        dir.create(output)
      } 
      save(EvoDynamics_exposure,file=paste0(output,"evoDynamic_exposure.RData"))
    }
    
    EvoDynamics_exposure
}

