
#' Produce input ccf count matrix for rank estimate
#' @name EvoSig_input
#' @param Ccube_folder Ccube folder
#' @param ccfupper CCF upper bound
#' @param ccfMatBuild Specify whether construct ccf count matrix or not
#' @param post_summary specify post summary file path if choose to load from local file
#' @param postSummaryBuild Specify whether construct post summary from scratch or not
#' @param output output folder path
#' @param minsample set minimum samples needed for a single cancer type
#' @return CCF matrix for each cancer type
#' @export
#' @import dplyr
EvoSig_input <- function(Ccube_folder,output,ccfMatBuild=F,post_summary=NA,postSummaryBuild=F,TCGA=F,ICGC=F,ccfupper=1,minsample=30,minmutation=30,mindepth=100,typefile=NA) {
  
  # Step 1.1: build post_summary or load from file
  if (is.na(post_summary) | postSummaryBuild) {
    if (TCGA) post_summary <- Build_post_summary(Ccube_folder,output=output,minsample=30,multicore=FALSE,TCGA=T)
    if (ICGC) post_summary <- Build_post_summary(Ccube_folder,output=output,minsample=30,multicore=FALSE,ICGC=T)
  } else {
    post_summary <- Build_post_summary(Ccube_folder,output=output,minsample=30,multicore=FALSE,typefile=typefile)
  }
  
  # Step 1.2: filter post_summary
  post_summary <- subset(post_summary,ave_depth >=mindepth)
  post_summary <- subset(post_summary,n_mutations >=minmutation)
  
  types <- post_summary %>% group_by(cancertype) %>% summarize(n=n())  %>% filter(n>=minsample) 
  
  post_summary <- subset(post_summary,cancertype %in% types$cancertype)
  
  save(post_summary,file=paste0(output,nrow(post_summary),"_post_summary_x100_n",minsample,".RData"))
  write.csv(types,file=paste0("type_summary_x100_n",minsample,".csv"))
  
  # Step 2: Construct ccf matrix for samples in the post_summary by types and for all and save the result to Matrix_folder
  if (ccfMatBuild){
    ccfMatBuild_output(post_summary = post_summary, input_folder=Ccube_folder,output=output,ccfupper = ccfupper,
                       RankEstimateType="count",add_samplename=TRUE)  
  }
  
  return(post_summary)

}


