
#' Build count matrix for input samples
#' @name Build_post_summary
#' @param input ccf file foder path
#' @param output output path 
#' @param typefile type with sample-type mapping
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
#' @import dplyr
#' @import ggplot2
Build_post_summary <- function(input,output=NA, typefile="NA",minsample=30){
  
  post_summary_analyse <- function(samplename){
    
    suppressWarnings(rm(ssm))
    ssm <- load_ccf(samplename,input=input)
    ccf <- unique(ssm$ccube_ccf_mean)
    ccf_mean_order <- sort(ccf,decreasing = T)
    Ncluster <- length(ccf)
    
    post_summary <- data.frame(samplename <- samplename) %>%
      mutate(
        Tumor_Sample_Barcode = unique(ssm$Tumor_Sample_Barcode),
        n_mutations = nrow(ssm),
        ccf_01_percentage = mean(ssm$ccube_ccf<=1),
        ccf_02_percentage = mean(ssm$ccube_ccf<=2),
        Ncluster = Ncluster,
        purity = unique(ssm$purity),
        ccube_purity = ifelse(exists("ssm$ccube_purity"),unique(ssm$ccube_purity),NA),
        ccf_mean_cluster1 = ifelse(Ncluster>=1,ccf_mean_order[1],0),
        ccf_mean_cluster2 = ifelse(Ncluster>=2,ccf_mean_order[2],0),
        ccf_mean_cluster3 = ifelse(Ncluster>=3,ccf_mean_order[3],0),
        ccf_mean_cluster4 = ifelse(Ncluster>=4,ccf_mean_order[4],0),
        ccf_mean_cluster5 = ifelse(Ncluster>=5,ccf_mean_order[5],0)
      ) %>% set_colnames(colnames)
    
    return(post_summary)
    
  }
  sample_list <- dir(input)
  
  colnames <- c("samplename","Tumor_Sample_Barcode","n_mutations","ccf_0-1_percentage","ccf_0-2_percentage","Ncluster","purity","ccube_purity","ccf_mean_cluster1","ccf_mean_cluster2","ccf_mean_cluster3","ccf_mean_cluster4","ccf_mean_cluster5")
  
  post_summary <- do.call(rbind,lapply(sample_list,post_summary_analyse))
  
  if (!is.na(typefile)) {
    cancertype <- read.csv(file=typefile)[,-1] 
    colnames(cancertype)[1] <- "samplename"
    post_summary <- left_join(post_summary, cancertype,by="samplename") 
    if ("Types" %in% colnames(post_summary)) 
      colnames(post_summary)[which(colnames(post_summary)=="Types")] <- "cancertype"
  }
  
  if (!is.na(output)) {
    
    if (!dir.exists(output)) {
      dir.create(output)
    } 
    
    write.csv(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".csv"))
    save(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".RData"))
  }
  
  xx <- post_summary %>% group_by(cancertype) %>% summarize(n=n())
  
  ggplot(data=xx,aes(x=cancertype,y=n)) +
    geom_bar(position = 'dodge', stat='identity') +
    geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=-0.25)
  
  ggsave(paste0(output,"type_summary.pdf"),width = 50, height = 15, units = "cm")
  
  ## filter types with minimum samples
  types <- xx %>% filter(n>minsample) %>% select(cancertype) 
  
  post_summary <- subset(post_summary,cancertype %in% types$cancertype)
  
  save(post_summary,file=paste0(output,nrow(post_summary),"_filtered_post_summary(n>",minsample,").RData"))
  
  return(post_summary)
}


#' Build CCF Matrix for all samples within region 0-upper
#' @name CountMatBuild
#' @param samplelist samples' name
#' @param upper set CCF upper bound
#' @param input_folder folder path stores ccf files 
#' @param genelist specify whether to only count for certain genes
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
CountMatBuild <- function(samplelist,upper,input_folder,genelist=NA){
    
    n_sample <- length(samplelist)
    rows <- upper/0.01 + 1
    ccfBand <- seq(0,upper,length.out = rows)
    ccfBandCountsMat <- matrix(nrow = rows,ncol = n_sample)
    
    if ( n_sample == 0 ) stop("The number of choosen sample is 0")
    
    for (j in 1:n_sample){
      sample_name <- as.character(samplelist[j])
      ssm <- load_ccf(sample_name,input = input_folder)
      
      if (!is.na(genelist)) ssm <- subset(ssm,.data$SYMBOL %in% genelist)
      
      matchBandCountDf <- data.frame(Var1=as.character(1:rows)) %>%
        suppressWarnings(left_join(as.data.frame(table(findInterval(ssm$ccube_ccf,ccfBand))),stringAsFactors = F,by="Var1"))
      
      ccfBandCountsMat[,j] <- matchBandCountDf$Freq
    }
    
    # delete the final row and set NAs to 0
    ccfBandCountsMat <- ccfBandCountsMat[-rows,]
    ccfBandCountsMat[is.na(ccfBandCountsMat)] <- 0
    
    if (n_sample == 1) {
      ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
      ccfBandFractionMat.random <- ccfBandFractionMat[sample(1:(rows-1))]
      ccfBandCountsMat.random <- ccfBandCountsMat[sample(1:(rows-1))]
    } else {
      ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
      ccfBandCountsMat.random <- NMF::randomize(ccfBandCountsMat)
      ccfBandFractionMat.random <- apply(ccfBandCountsMat.random,2,function(x) x/sum(x))
    }
    
    return(list(ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat,ccfBandFractionMat.random))
    
  }

#' Build count matrix for input samples per cancer type
#' @name TypeMatrixBuild 
#' @param post_summary samples' name
#' @param ccfupper set CCF upper bound
#' @param input_folder folder path stores ccf files 
#' @param genelist specify whether to only count for certain genes
#' @param output output folder path
#' @param RankEstimateType output rank estimate matirx format
#' @return CCF matrix for each cancer type
#' @export
TypeMatrixBuild <- function(post_summary,input_folder,output,ccfupper=1,RankEstimateType="fraction"){
  
  type <- post_summary %>% group_by(.data$cancertype)%>% summarize(n=n())
  
  typelist <- as.character(subset(type,n>=30)$cancertype)
  
  ntype <- length(typelist)
  
  print(paste0("Construct ccf count matrix for ",ntype," types (> 30 samples)"))

  if (!dir.exists(output)) dir.create(output)
  
  for (i in 1:ntype){
    
    # get sample list for each type
    type <-  typelist[i]
    
    print(paste0("finish load ",i,'th type - ',type))
    
    samplelist <- subset(post_summary,cancertype==type)$samplename
    n_sample <- length(samplelist)
    
    ccfBandCountsMat <- suppressWarnings(CountMatBuild(samplelist,input_folder = input_folder,upper=ccfupper))
 
    if (!is.na(output)){
        
        ccfCountMatrix <-ccfBandCountsMat[[1]]
        ccfCountsMatrix.random <- ccfBandCountsMat[[2]]
        ccfFractionMatrix <- ccfBandCountsMat[[3]]
        ccfFractionMatrix.random <- ccfBandCountsMat[[4]]
        
        # create direction
        if (!dir.exists(paste0(output,type,"/"))) dir.create(paste0(output,type,"/"))  
        if (!dir.exists(paste0(output,"rank_estimate/"))) dir.create(paste0(output,"rank_estimate/"))  
                                                                    
        # Output Matrix 
        output_format <- paste0(output,type,"/",type,"_", n_sample,"_0-",ccfupper)
        output_rank <- paste0(output,"rank_estimate/",type,"_", n_sample,"_0-",ccfupper)
      
        write.csv(samplelist,file=paste0(output,type,"/",type,"_",n_sample,"_samplelist_",Sys.Date(),".csv"))
        
        save(ccfCountMatrix,file=paste0(output_format,"_ccfCountMatrix_",Sys.Date(),".RData"))
        save(ccfCountsMatrix.random,file=paste0(output_format,"_ccfCountsMatrix.random_",Sys.Date(),".RData"))
        save(ccfFractionMatrix,file=paste0(output_format,"_ccfFractionMatrix_",Sys.Date(),".RData"))
        save(ccfFractionMatrix.random,file=paste0(output_format,"_ccfFractionMatrix.random_",Sys.Date(),".RData"))
        
        # output to rank estimate folder
        if (RankEstimateType=="fraction") {
          write.csv(ccfFractionMatrix.random,file=paste0(output_rank,"_ccfFractionMatrix.random_",Sys.Date(),".csv"))
          write.csv(ccfFractionMatrix,file=paste0(output_rank,"_ccfFractionMatrix_",Sys.Date(),".csv"))
        }
        
        if (RankEstimateType=="count") {
          write.csv(ccfCountsMatrix.random,file=paste0(output_rank,"_ccfCountMatrix.random_",Sys.Date(),".csv"))
          write.csv(ccfCountMatrix,file=paste0(output_format,"_ccfCountMatrix_",Sys.Date(),".csv"))
        }
    }
  }
}
