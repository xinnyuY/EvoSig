
ccube_input
library(EvoSig)

Ccube_folder = "D:/Project/Data/TCGA/Ccube/TCGA_0106/"
typefile= "D:/Project/Data/TCGA/Ccube/TCGA_cancertype.csv"
output = "D:/Project/Xinyu/EvoSig.test/3.10_x100/"


EvoSig_input <- function(ccfMat=NA,consensus_sig,Ccube_folder,output,ccfMatBuild=F,post_summary=NA,postSummaryBuild=F,TCGA=F,ICGC=F,ccfupper=1) {
  
  # Step1: build post_summary
  if (is.na(post_summary) | postSummaryBuild) {
    if (TCGA) post_summary <- Build_post_summary(Ccube_folder,output=output, typefile="TCGA_type_sample_mapping.rda",minsample=30,multicore=FALSE)
  }
  data()
  data("TCGA_type_sample_mapping.rda")
  if (!is.na(post_summary)) load(post_summary)
  
  data("TCGA_type")
  # Step 2: Construct ccf matrix for samples in the post_summary by types and for all and save the result to Matrix_folder
  Matrix_folder <- paste0(output,"ccfMat/")
  ccfMatBuild_output(post_summary = post_summary, input_folder=Ccube_folder,output=Matrix_folder,ccfupper = ccfupper,
                     RankEstimateType="count",add_samplename=TRUE)  
  if (!is.na(ccfMat)) load(ccfMat)
  
  # if (is.na(ccfMat) & ccfMatBuild) {
  #   EvoDynamics_exposure <- Extract_sig(ccfMat=ccfMat,consensus_sig = consensus_sig,samplelist=post_summary$samplename)
  # }
  
  
 
}

