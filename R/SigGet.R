# Unify TCGA sample barcode
file_format <- function(filename=filename,samplenamecol){

  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  
  return(filename)
  
}


# Construct ccf distribution count matrix for given samples
CountMatBuild <- function(samplelist){
  
  ccf_folder <- paste0(input_folder,type,"/")
  
  n_sample <- length(samplelist)
  ccfBand <- seq(0,1,length.out = 101)
  ccfBandCountsMat <- matrix(nrow=length(ccfBand),ncol=n_sample)
  
  if (n_sample !=0) stop("The number of choosen sample is 0")
  
  for (j in 1:n_sample){
    sample_name <- as.character(samplelist[j])
    
    if (file.exists(paste0(tcga_folder,sample_name,"/ccube_res_v0.3_final_new.RData"))) load(paste0(ccf_folder,sample_name,"/ccube_res_v0.3_final_new.RData"))
    if (file.exists(paste0(tcga_folder,sample_name,"/ccube_result.RData"))) load(paste0(ccf_folder,sample_name,"/ccube_result.RData"))
   
    ssm <- ccubeRes[['ssm']]
    
    matchBandCountDf <- data.frame(Var1=as.character(1:101))
    matchBandCountDf <- left_join(matchBandCountDf,as.data.frame(table(findInterval(ssm$ccube_ccf,ccfBand))),stringAsFactors = F,by="Var1")
    ccfBandCountsMat[,j] <- matchBandCountDf$Freq
    
  }
  
  return(ccfBandCountsMat)
  
}
