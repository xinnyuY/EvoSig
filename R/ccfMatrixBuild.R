
#' Build count matrix for input samples
#' @name Build_post_summary
#' @param input ccf file foder path
#' @param output output path 
#' @param typefile type with sample-type mapping
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
#' @import dplyr
#' @import ggplot2
Build_post_summary <- function(input,output=NA, typefile=NA,minsample=30,multicore=FALSE,TCGA=F,ICGC=F){
  
  cat("Start building post summary for",length(dir(input)),"input CCF files \n")
  if (multicore) {
    if (!exists("cl")) cl <- makeCluster(detectCores() - 1, type="SOCK") 
    registerDoParallel(cl)  
  }
  

  post_summary_analyse_tcga <- function(samplename,TCGA=F,ICGC=F){
    
    colnames <- c("samplename","Tumor_Sample_Barcode","n_mutations","ave_depth","ccf_0-1_percentage","ccf_0-2_percentage","Ncluster","purity","ccube_purity","ccf_mean_cluster1","ccf_mean_cluster2","ccf_mean_cluster3","ccf_mean_cluster4","ccf_mean_cluster5")
    
    ssm <- load_ccf(samplename,input=input)
    ccf <- unique(ssm$ccube_ccf_mean)
    ccf_mean_order <- sort(ccf,decreasing = T)
    Ncluster <- length(ccf)
    ave_deapth <- mean(ssm$total_counts)
    
    post_summary <- data.frame(samplename <- samplename) %>%
      mutate(
        Tumor_Sample_Barcode = unique(ssm$Tumor_Sample_Barcode),
        n_mutations = nrow(ssm),
        ave_deapth <-ave_deapth,
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
    
     cat(">")
    return(post_summary)
  }
  
  post_summary_analyse_icgc <- function(samplename,TCGA=F,ICGC=F){
    
    colnames <- c("samplename","n_mutations","ave_depth","ccf_0-1_percentage","ccf_0-2_percentage","Ncluster","purity","ccube_purity","ccf_mean_cluster1","ccf_mean_cluster2","ccf_mean_cluster3","ccf_mean_cluster4","ccf_mean_cluster5")
    
    ssm <- load_ccf(samplename,input=input)
    ccf <- unique(ssm$ccube_ccf_mean)
    ccf_mean_order <- sort(ccf,decreasing = T)
    Ncluster <- length(ccf)
    ave_deapth <- mean(ssm$total_counts)

    post_summary <- data.frame(samplename <- samplename) %>%
        mutate(
          n_mutations = nrow(ssm),
          ave_deapth <- ave_deapth,
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
    
    cat(">")
    return(post_summary)
  }
  
  sample_list <- dir(input)
  
 
  n_sample <- nrow(sample_list)
  
  #if (multicore) post_summary <- do.call(rbind,foreach(1:n_sample) %dopar% post_summary_analyse(sample_list[i]))
  if (TCGA) post_summary <- do.call(rbind,lapply(sample_list,post_summary_analyse_tcga))
  if (ICGC) post_summary <- do.call(rbind,lapply(sample_list,post_summary_analyse_icgc))
  
  # add type variable
  if ( (!is.na(typefile) & !TCGA) | (!is.na(typefile) & !ICGC) ) 
    cancertype <- read.csv(file=typefile)[,-1] 
  
  if (TCGA) {
    data("TCGAtypefile") 
    cancertype <- TCGAtypefile
  }
  
  if (ICGC) {
    data("ICGCtypefile") 
    cancertype <- ICGCtypefile
  }
  
  if (colnames(cancertype)[1]=="Tumor_sample_barcode") 
    colnames(cancertype)[1] <- "samplename"
  
  post_summary <- suppressWarnings(left_join(post_summary, cancertype,by="samplename"))
  
  if ("Types" %in% colnames(post_summary)) 
      colnames(post_summary)[which(colnames(post_summary)=="Types")] <- "cancertype"


  # Output
  if (!is.na(output)) {
    
    if (!dir.exists(output)) dir.create(output,recursive = T)
    
    write.csv(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".csv"))
    save(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".RData"))

  }
  
  if (multicore) stopCluster(cl)
  
  cat("\n Done \n")
  post_summary
  
}


#' Build CCF Matrix for all samples within region 0-upper
#' @name ccfMatBuild
#' @param samplelist samples' nam+e
#' @param upper set CCF upper bound
#' @param input_folder folder path stores ccf files 
#' @param genelist specify whether to only count for certain genes
#' @param add_samplename specify whether to add samplename variables
#' @param binwidth ccf bin width
#' @param plot whether to plot mutation distribution for samples
#' @param random whether to output randomized ccf count matrix
#' @return A list containing ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat and ccfBandFractionMat.random
#' @export
#' @importFrom NMF randomize
#' @import dplyr
#' @import reshape2
#' @import ggplot2
#' @import ggpubr
ccfMatBuild <- function(samplelist,upper,input_folder,binwidth=0.01,add_samplename=FALSE,plot=FALSE,random=FALSE){
  
  n_sample <- length(samplelist)
  rows <- upper/binwidth + 1
  ccfBand <- seq(0,upper,length.out = rows)
  ccfBandCountsMat <- matrix(nrow = rows,ncol = n_sample)
  
  if (n_sample == 0) stop("The number of choosen sample is 0")
  
  for (j in 1:n_sample){
    sample_name <- as.character(samplelist[j])
    ssm <- load_ccf(sample_name,input = input_folder)
    
    #if (!is.na(genelist)) ssm <- subset(ssm,.data$SYMBOL %in% genelist)
    matchBandCountDf <- data.frame(Var1=as.character(1:rows))
    matchBandCountDf <- suppressWarnings(left_join(matchBandCountDf,as.data.frame(table(findInterval(ssm$ccube_ccf,ccfBand))),stringAsFactors = F,by="Var1"))
    
    ccfBandCountsMat[,j] <- matchBandCountDf$Freq
  }
  
  if (random==TRUE){
    
    if (n_sample == 1) {
      ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
      ccfBandFractionMat.random <- ccfBandFractionMat[sample(1:(rows-1))]
      ccfBandCountsMat.random <- ccfBandCountsMat[sample(1:(rows-1))]
    } else {
      ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
      ccfBandCountsMat.random <- randomize(ccfBandCountsMat)
      ccfBandFractionMat.random <- apply(ccfBandCountsMat.random,2,function(x) x/sum(x))
    }
    
    outputlist <- list(ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat,ccfBandFractionMat.random)
    
  } else {
    
    if (n_sample == 1) {
      ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
    } else {
      ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
    }
    outputlist <- list(ccfBandCountsMat,ccfBandFractionMat)
  }
  
  if (add_samplename==TRUE) {
    outputlist <- lapply(outputlist,function(x){x <- as.data.frame(t(x)) %>% mutate(samplename=samplelist)})
  }
  
  
  if (plot==TRUE) {
    
    ccf_plot <- outputlist[1] %>% as.data.frame()
    
    null_row = which(cumsum(rev(colSums(ccf_plot[,1:rows],na.rm=TRUE)))==0)
    
    
    if (length(null_row)>0) {
      min_rows = rows-max(null_row)
      ccf_plot[,(rows-max(null_row)+1):rows] = NULL } else {
        min_rows = rows
      }
    ccf_plot[is.na(ccf_plot)] <- 0
    
    p_ccf <- ccf_plot %>% melt(id='samplename') %>% mutate(variable=as.numeric(variable))%>%
      ggplot(aes(x=variable,y=value))+geom_bar(stat='identity')+ 
      theme_pubr()+
      labs(x="Cancer Cell Fraction",y="Counts")+
      facet_wrap(~samplename,ncol=11,nrow=70,scale="free")  +   
      theme(axis.text.x = element_text(size=10))
  }
  return(list(outputlist,p_ccf))
}



#' Build count matrix for input samples per cancer type
#' @name ccfMatrixBuild
#' @param post_summary samples' name
#' @param ccfupper set CCF upper bound
#' @param input_folder folder path stores ccf files 
#' @param genelist specify whether to only count for certain genes
#' @param output output folder path
#' @param RankEstimateType output rank estimate matirx format
#' @return CCF matrix for each cancer type
#' @export
#' @import dplyr
ccfMatBuild_output <- function(post_summary,input_folder,output=NA,ccfupper=1,RankEstimateType="fraction",add_samplename=TRUE,minsample=30){
  
      samplelist_all <- post_summary$samplename
      type <- post_summary %>% group_by(.data$cancertype)%>% summarize(n=n())
      typelist <- as.character(subset(type,n>=30)$cancertype)
      ntype <- length(typelist)
  
      cat(paste0("\n Start constructing ccf count matrix for ",ntype," types (>",minsample," samples) \n"))

      if (!is.na(output)) {
        ccfOutput_path <- paste0(output,"ccfMat/")
        ccfOutput_all_path <- paste0(output,"ccfMat/All/")
        ccfOutput_rank_path <- paste0(output,"ccfMat/rank_estimate/")
        
        multi_dir_create(c(ccfOutput_all_path,ccfOutput_rank_path))

        for (i in 1:ntype){
          
          type <-  typelist[i]
          cat(paste0("\n >>>> loading ",i,'th type - ',type," <<<< \n"))
          
          samplelist_type <- subset(post_summary,cancertype==type)$samplename
          n_sample <- length(samplelist_type)
          
          ccfBandCountsMat <- suppressWarnings(ccfMatBuild(samplelist_type,input_folder = input_folder,upper=ccfupper,add_samplename = add_samplename))
          
          ccfCountMatrix <-ccfBandCountsMat[[1]]
          ccfCountsMatrix.random <- ccfBandCountsMat[[2]]
          ccfFractionMatrix <- ccfBandCountsMat[[3]]
          ccfFractionMatrix.random <- ccfBandCountsMat[[4]]
          
          # create directio
          multi_dir_create(paste0(ccfOutput_path,type,"/"))  
         
          # Output Matrix 
          output_format <- paste0(ccfOutput_path,type,"/",type,"_", n_sample,"_0-",ccfupper)
          output_format_rank <- paste0(ccfOutput_rank_path,type,"_", n_sample,"_0-",ccfupper)
          
          save(ccfCountMatrix,file=paste0(output_format,"_ccfCountMatrix_",Sys.Date(),".RData"))
          save(ccfCountsMatrix.random,file=paste0(output_format,"_ccfCountMatrix.random_",Sys.Date(),".RData"))
          save(ccfFractionMatrix,file=paste0(output_format,"_ccfFractionMatrix_",Sys.Date(),".RData"))
          save(ccfFractionMatrix.random,file=paste0(output_format,"_ccfFractionMatrix.random_",Sys.Date(),".RData"))
          
          # output to rank estimate folder
          if (RankEstimateType=="fraction") {
            write.csv(ccfFractionMatrix.random,file=paste0(output_format_rank,"_ccfFractionMatrix.random_",Sys.Date(),".csv"))
            write.csv(ccfFractionMatrix,file=paste0(output_format_rank,"_ccfFractionMatrix_",Sys.Date(),".csv"))
          }
          
          if (RankEstimateType=="count") {
            write.csv(ccfCountsMatrix.random,file=paste0(output_format_rank,"_ccfCountMatrix.random_",Sys.Date(),".csv"))
            write.csv(ccfCountMatrix,file=paste0(output_format_rank,"_ccfCountMatrix_",Sys.Date(),".csv"))
          }
        }
      }
     
        # ouput ccf matrix for all samples
        ccfBandCountsMat_all <- suppressWarnings(ccfMatBuild(samplelist_all,input_folder = input_folder,upper=ccfupper,add_samplename = add_samplename))
    
        ccfCountMatrix_all <- ccfBandCountsMat_all[[1]]
        ccfFractionMatrix_all <- ccfBandCountsMat_all[[3]]
          
        output_all_format <- paste0(ccfOutput_all_path,"All_",length(samplelist_all),"_0-",ccfupper)
          
        # Output Matrix 
        save(ccfCountMatrix_all,file=paste0(output_all_format,"_ccfCountMatrix_",Sys.Date(),".RData"))
        save( ccfFractionMatrix_all,file=paste0(output_all_format,"_ccfFractionMatrix_",Sys.Date(),".RData"))
    
        cat(paste0("\n Done, saved results to ",output))
}
