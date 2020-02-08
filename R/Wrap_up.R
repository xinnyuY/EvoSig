EvoDynamicsSigExtract <- function(ccfMat=NA,consensus_sig,input,output,ccfMatBuild=F,post_summary=NA,postSummaryBuild=F) {
  
  if (is.na(post_summary) & postSummaryBuild) {
    post_summary <- Build_post_summary(input,output=NA, typefile="TCGA_type_sample_mapping.rda",minsample=30)
  }
  
  if (!is.na(post_summary)) load(post_summary)
  
  if (!is.na(ccfMat)) load(ccfMat)
  
  if (is.na(ccfMat) & ccfMatBuild) {
    EvoDynamics_exposure <- Extract_sig(ccfMat=ccfMat,consensus_sig = consensus_sig,samplelist=post_summary$samplename)
  }
  
  
 
}

