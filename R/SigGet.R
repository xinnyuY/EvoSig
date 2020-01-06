## Basic function

file_format <- function(filename=filename,samplenamecol){
  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  return(filename)
}

file_format <- function(filename=filename,samplenamecol){
  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  return(filename)
}

load_ccf <- function(samplename,input){
  Check <- ArgumentCheck::newArgCheck()
  suppressWarnings(rm(ssm,res,ccubeRes))
  
  format1 <- paste0(input,samplename,"/ccube_result.RData")
  format2 <- paste0(input,samplename,"/ccube_res_v0.3_final_new.RData")
  format3 <- paste0(input,samplename,"/ccube_res_run1.RData")
  format4 <- paste0(input,samplename,"/ccubeRes.RData")
  
  if (file.exists(format4 )) load(format4)  
    else if (file.exists(format1 )) load(format1)  
     else if(file.exists(format2)) load(format2) 
      else if(file.exists(format3)) load(format3) 
        else{
          ArgumentCheck::addError(
          msg = "No file has been loaded",
          argcheck = Check)
          }
  
  if (exists("ccubeRes")) return(ccubeRes$ssm)  else
     if (exists("res")) return(res$ssm) else
       if (exists("ssm")) return(ssm)
}

## Step1: Ccube input format
  ParseSnvCnaPcawgFormat <- function (ssm, cna) {
      library(dplyr)
      
      ssm <- ssm %>%
        mutate(chr= substr(chr,4,length(chr)),
               cn_frac = NA,
               major_cn = NA,
               minor_cn = NA,
               mutation_id = NA)
      
      for (jj in seq_len(nrow(cna)) ) {
        cc = cna[jj,]
        
        idx = which(ssm$chr == cc$chromosome &  (ssm$Start_Position >= cc$start & ssm$End_Position <= cc$end) )
        
        if (length(idx) > 0) {
          ssm[idx,] <- ssm[idx,] %>% 
            mutate( major_cn=cc$major_cn,minor_cn =cc$minor_cn, cn_frac = 1)
        }
      }
      
      ssm$mutation_id = paste0(ssm$chr, "_", ssm$Start_Position )
      
      ssm <- ssm %>%
        select(-chr,-Start_Position,-End_Position,-df.n_alt_count,-n_ref_count) %>%
        rename(var_counts='t_alt_count',ref_counts='t_ref_count') %>%
        mutate(total_counts=var_counts+ref_counts,normal_cn=2) %>%
        filter(!is.na(major_cn) ,!is.na(minor_cn),!is.na(cn_frac),major_cn > 0)
      
      return(ssm)
    }
  
 ## Step2: Build post summary for ccube samples 
Build_post_summary <- function(input,output=NA,mergeType= F, typefile="NA"){
  library(magrittr)
  library(dplyr)
  
  sample_list <- dir(input)
  
  colnames <- c("samplename","Tumor_Sample_Barcode","n_mutations","ccf_0-1_percentage","ccf_0-2_percentage","Ncluster","purity","ccube_purity","ccf_mean_cluster1","ccf_mean_cluster2","ccf_mean_cluster3","ccf_mean_cluster4","ccf_mean_cluster5")

  # if (i==1) cat("Start building post summary for",length(sample_list),"files","\n")
  # if (i%%1000==0) print(paste0("----- Finish loading ",i," files -----"))

  post_summary_analyse <- function(samplename){
      library(dplyr)
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
     
  post_summary <- do.call(rbind,lapply(sample_list,post_summary_analyse))
  
  if (mergeType==T & !is.na(typefile)) {
    cancertype <- read.csv(file=typefile)[,-1] 
    colnames(cancertype)[1] <- "samplename"
    post_summary <- left_join(post_summary, cancertype,by="samplename")
  }
  
  if (!is.na(output)) {
    write.csv(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".csv"))
    save(post_summary,file=paste0(output,"post_summary_",length(sample_list),"_",Sys.Date(),".RData"))
  }
  return(post_summary)
}

 ## Step2: Construct ccf matrix
 CountMatBuild <- function(samplelist,upper,input_folder,genelist=NA){
  library(dplyr)
  
  n_sample <- length(samplelist)
  rows <- upper/0.01 + 1
  
  ccfBand <- seq(0,upper,length.out = rows)
  ccfBandCountsMat <- matrix(nrow=rows,ncol=n_sample)
  
  if (n_sample ==0) stop("The number of choosen sample is 0")
  
  for (j in 1:n_sample){
    sample_name <- as.character(samplelist[j])
    ssm <- load_ccf(sample_name,input=input_folder)
    if (!is.na(genelist)) ssm <- subset(ssm,SYMBOL %in% genelist)
    matchBandCountDf <- data.frame(Var1=as.character(1:rows))
    matchBandCountDf <- suppressWarnings(left_join(matchBandCountDf,as.data.frame(table(findInterval(ssm$ccube_ccf,ccfBand))),stringAsFactors = F,by="Var1"))
    ccfBandCountsMat[,j] <- matchBandCountDf$Freq
  }
  # delete the final row
  ccfBandCountsMat <- ccfBandCountsMat[-rows,]
  ccfBandCountsMat[is.na(ccfBandCountsMat)] <- 0
  
  if (n_sample==1) {ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
                    ccfBandFractionMat.random <- ccfBandFractionMat[sample(1:(rows-1))]
                    ccfBandCountsMat.random <- ccfBandCountsMat[sample(1:(rows-1))]} else 
    {ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
     ccfBandCountsMat.random <- NMF::randomize(ccfBandCountsMat)
     ccfBandFractionMat.random <- apply(ccfBandCountsMat.random,2,function(x) x/sum(x))
     }
    
  return(list(ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat,ccfBandFractionMat.random))
}
                                        
TypeMatrixBuild <- function(post_summary,input_folder,output,ccfupper=1,RankEstimateType="fraction"){
  
  type <- post_summary %>%
    group_by(cancertype)%>%
    summarize(n=n())
  
  typelist <- as.character(subset(type,n>=30)$cancertype)
  
  ntype <- length(typelist)
  
  print(paste0("Construct ccf count matrix for ",ntype," types (> 30 samples)"))

  if (!dir.exists(output)) dir.create(output)
  
  for (i in 1:ntype){
    
    library(NMF)
    library(dplyr)
    
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
                                        
## Step 5: NMF for Each type
sig_plot <- function(sig){
  library(gridExtra)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(dplyr)
  library(reshape2)
  library(cowplot)
  library(magrittr)
  
  xx <- as.data.frame(t(apply(sig,2,function(x) x/sum(x)))) %>%
    set_colnames(1:ncol(.)) %>%
    mutate(signature=paste0("Signature ",1:nrow(.))) %>%
    melt(.,id=c("signature"))
  
  fills <- brewer.pal(8, "Set3")[c(1,3:8,2)]
  
  p1 <- ggplot(xx,aes(y=value,x=variable)) + geom_bar(aes(fill=signature),stat='identity') + scale_fill_manual(values = fills)+ theme_grey()+
    #ggtitle(paste0("rank = ",rank,", cancertype = ",type, ", MatrixType = ",MatType )) + 
    theme(legend.title = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(color= "white",size=10),
          panel.grid.minor.y = element_blank(),
          axis.title.x = element_text(color = "grey20", size = 8),
          axis.text.x = element_text(color = "grey20", size = 6),
          axis.text.y = element_text(color = "grey20", size = 6),
          plot.title = element_text(size= 10)) +
    facet_grid(cols = vars(signature))+ scale_x_discrete(breaks=c("1","50","100") ,labels=c("0", "0.5", "1"))+xlab("Cancer Cell Fraction") + ylab("") 
  
  g1 <- ggplot_gtable(ggplot_build(p1))
  
  strip_both <- which(grepl('strip-', g1$layout$name))
  
  k <- 1
  for (i in strip_both) {
    j1 <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
    g1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
    #g1$layout$clip[j3] <- "off"
    k <- k+1
  }
  grid.draw(g1)
  #output sig and expo
  return(g1)
}
   
nmf_sig_plot_type <- function(type,MatType="fraction",input_folder,output,rank_summary){
     library(magrittr)
  
     type_path <- paste0(input_folder,type,"/")
     if (MatType=="fraction"){
     file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfFractionMatrix_",dir( type_path))])
     } 
  
     if (MatType=="count"){
     file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfCountsMatrix_",dir( type_path))])
     }
     
     samplename_path <- paste0(input_folder,type,"/",dir(type_path)[grep("samplelist",dir(type_path))]) 
     
     if (!dir.exists(paste0(output,type))) dir.create(paste0(output,type))
       
     load(file=file_path)
     samplename <- as.character(read.csv(file=samplename_path)[,-1])
     n_sample <- length(samplename)
     
     #format rank summary file
     rank <- as.numeric(subset(rank_summary,Type == type)$Rank)
     
     if (exists("ccfFractionMatrix")) ccfMat <- ccfFractionMatrix
     if (exists("ccfCountsMatrix")) ccfMat <- ccfCountsMatrix
     
     index <- sample(1:ncol(ccfMat))
     ccfMat <- ccfMat[,index]
     samplename_random <- samplename[index]
     
     #preprocess for rows with all 0
     index_p <- which(rowSums(ccfMat)>0)
     index_n <- which(!rowSums(ccfMat)>0)
     ccfMat <- ccfMat[which(rowSums(ccfMat)>0),]
     
     #run NMF
     res <- nmf(ccfMat,rank,.opt='vp4')
     sig <- as.data.frame(matrix(0,nrow=length(index_p)+length(index_n),ncol=ncol(res@fit@W)))
     sig[c(index_p),] <- as.data.frame(res@fit@W) %>%
       set_colnames(paste0("sig_",1:ncol(.)))
     
     expo <- as.matrix(res@fit@H)
     
     #sig plot
     
     p_sig <- sig_plot(sig)
   
     #output sig and expo
     expo <- as.data.frame(t(expo)) %>%
        cbind(.,samplename_random) 
     colnames(expo)[1:rank] <- paste0("sig_",1:rank)
     
     save(expo,file=paste0(output,type,'/',type,'_',n_sample,"_expowithsample_",Sys.Date(),".RData"))
     save(sig,file=paste0(output,type,'/',type,'_',n_sample,"_sig_",Sys.Date(),".RData"))
     save(res,file=paste0(output,type,'/',type,'_',n_sample,"_res_",Sys.Date(),".RData"))
   
     return(list(p_sig,n_sample))
}                         
                              
nmf_sig_all_plot <- function(input_folder,output,rank_summary) {
  library(gridExtra)
  library(ggplot2)
  library(grid)
  library(dplyr)
  library(reshape2)
  library(cowplot)
  library(NMF)
  
  if (!dir.exists(output)) dir.create(output)
  
  colnames(rank_summary) <- c("cancertype","rank")
  cancertype <- as.character(rank_summary$cancertype)
 
  j = 0
  for (i in 1:length(cancertype)) {
    tryCatch({
    type <- cancertype[i]
   
    print(paste0("load ",i,"th type - ",type))
    
    nmf_result <- nmf_sig_plot_type(type,input_folder=input_folder,output=output,rank_summary=rank_summary)
    p_sig <- nmf_result[[1]]
    n_sample <- nmf_result[[2]]
  
    p <- plot_grid(p_sig)
    title <- ggdraw() + draw_label(paste0("Extracting Signatures for ",type," with ",n_sample," samples "),fontface='bold')
    j <- j+1
    commands <- paste0("p",j," <- plot_grid(title,p,align='V',ncol=1,rel_heights = c(0.1,0.15,0.65))")
    eval(parse(text=commands))
    },error = function(e) print(paste0("Failed run NMF on this type")))
  }
    #plot_grid(p1)
    
    ntype <- length(cancertype)
    
    n_file <-  j %% 8
    
    if (n_file==0) n_pdf <- (j %/% 8) else
      n_pdf <- (j %/% 8)+1
    
    for (i in 1:n_pdf){
     
      if (i != n_pdf) {
        commands1 <- paste0("pdf(file=paste0(output,'Sig_summary_p",(8*i-7),"-",8*i,"_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("g",i,"<- plot_grid(",paste0("p",(8*i-7):(8*i),collapse=","),",align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))") 
      } else {
        
        commands1 <- paste0("pdf(file=paste0(output,'Sig_summary_p",(8*i-8),"-",8*i-8+n_file,"_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("g",i," <- plot_grid(",paste0("p",(8*i-8):(8*i-8+n_file),collapse=","),",align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))")
      }
      
      commands3 <- paste0("print(g",i,")")
      commands4 <- "dev.off()"
      eval(parse(text=commands1))
      eval(parse(text=commands2))
      eval(parse(text=commands3))
      eval(parse(text=commands4))
    }
  
    return( eval(parse(text=paste0("list(n_pdf,",paste0("g",1:n_pdf,collapse=","),")"))))
}

## Step 6: Combine sig and run Hierarchical cluster
combine_sig_nmf <- function(input_folder,cancertype){
  
  for (i in 1:length(cancertype)){
  tryCatch({
    type <- cancertype[i]
    print(paste0("load ",i," th type : ",type))
    
    sig_file <- dir(paste0(input_folder,type,"/"))[grep("sig",dir(paste0(input_folder,type,"/")))]
    load(paste0(input_folder,type,"/",sig_file))
  
    colnames(sig) <- paste0(type,"_sig",1:ncol(sig))
    if (i==1) combine_sig <- sig else
       combine_sig <- cbind(combine_sig,sig)
  },error=function(e) print("Fail load this type"))
  }  
   save(combine_sig,file=paste0(input_folder,"combine_sig_",Sys.Date(),".RData"))
   return(combine_sig)
}
                              
## Step 7: test the number of cluster
 hc_cluster_test <- function(data,methods,distance,min.nc = 2,max.nc = 10){
  library(NbClust)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  tabla = as.data.frame(matrix(ncol = length(distance), nrow = length(methods)))
    names(tabla) = distance
    
    for (j in 2:length(distance)){
      for(i in 1:length(methods)){
        nb = NbClust(data,distance = distance[j],
                     min.nc = min.nc, max.nc = max.nc, 
                     method = "complete", index =methods[i])
        tabla[i,j] = nb$Best.nc[1]
        tabla[i,1] = methods[i]
      }}
  
  tabla <- rbind(tabla,c("Most_common",apply(tabla[,2:5],2,getmode)))
  
  return(tabla)
} 

## Step 8: run Hierarchical clustering and plot (save to 8.HC_consensus)
hc_consensus <- function(combine_sig,cluster,output,distance="euclidean")  {
  library("RColorBrewer")
  library("pheatmap")
  #library(factoextra)
  
  upper = quantile(combine_sig,0.95)
  breaksList = seq(0, upper, by = 0.01)
  col <- colorRampPalette(rev(brewer.pal(n = 6, name = "RdYlBu")))(length(breaksList))
    
  ## set `cutree_cols` based on suggested cluster number 
  out <- pheatmap(combine_sig, cutree_cols = cluster, fontsize_col = 5,fontsize_row = 0.4,color = col, breaks = breaksList,clustering_distance_cols=distance, cluster_rows=F,filename=paste0(output,"/hc_heatmap.pdf"))
 
  sig_label <- as.data.frame(cutree(out$tree_col,k=6)) %>%
      set_colnames("cluster") %>%
      mutate(sig=rownames(.))
    
  ## Compute consensu signatures for each cluster
  consensus_sig <- as.data.frame(t(combine_sig)) %>%
      mutate(sig=rownames(.)) %>%
      left_join(.,sig_label,by="sig") %>%
      group_by(cluster) %>%
      mutate(sig=NULL) %>%
      summarise_all(mean) 
  
  save(consensus_sig,file=paste0(output,"/consensus_sig.RData"))
  consensus_sig <- apply(t(consensus_sig[,2:101]),2,as.numeric)
  
  p1 <- plot_grid(sig_plot(consensus_sig))
  save_plot(paste0(output,"/consensus_sig.pdf"),p1,base_asp = cluster)
 
  return(consensus_sig)
}
  
