# Unify TCGA sample barcode
file_format <- function(filename=filename,samplenamecol){

  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  
  return(filename)
  
}


##################################################################
### --- Step1: Construct Count Matrix by type ----- ####

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

# Construct ccf distribution count matrix for given cancer types 
TypeCountBuildTCGA <- function(cancertype,input_folder,output){
    
    x <- load(paste0(input_folder,"TCGA_summary.RData"))
    post_summary <- get(x)
    
    if (!dir.exists(output)) dir.create(output)
    
    for (i in 1:length(cancertype)){
      
      library(NMF)
      library(dplyr)
      
      # get sample list for each type
      type <- cancertype[i]
      
      print(paste0("finish load ",i,'th type - ',type))
    
      samplelist <- subset(post_summary,Types==type)$samplename
      n_sample <- length(samplelist)
  
      ccfBandCountsMat <- CountMatBuild(samplelist)
    
      if (!dir.exists(paste0(output,type,"/"))) dir.create(paste0(output,type,"/"))  
      if (!dir.exists(paste0(output,"Working/"))) dir.create(paste0(output,"Working/")) 
      
      write.csv(samplelist,file=paste0(output,type,"/",type,"_TCGA_", n_sample,"_samplelist.csv"))
      
      ccfCountsMat <- ccfBandCountsMat[[1]]
      ccfFractionMatrix <- ccfBandCountsMat[[2]]
      ccfFractionRandomMatrix <-ccfBandCountsMat[[3]]
      
      save(ccfCountsMat,file=paste0(output,type,"/",type,"_TCGA_",n_sample,"_0-1_countMatrix.RData"))
      write.csv( ccfCountsMat,file=paste0(output,type,"/",type,"_TCGA_",n_sample,"_0-1_countMatrix.csv"))
      write.csv(ccfFractionMatrix ,file=paste0(output,type,"/",type,"_TCGA_",n_sample,"_0-1_fractionMatrix.csv"))
      write.csv(ccfFractionRandomMatrix,file=paste0(output,type,"/",type,"_TCGA_",n_sample,"_0-1_fractionMatrix_random.csv"))
      write.csv(ccfFractionMatrix ,file=paste0(output,"Working/",type,"_TCGA_",n_sample,"_0-1_fractionMatrix.csv"))
      write.csv(ccfFractionRandomMatrix,file=paste0(output,"Working/",type,"_TCGA_",n_sample,"_0-1_fractionMatrix_random.csv"))
    }
}
TypeCountBuildICGC <- function(cancertype,input_folder,output){
  
  x <- load(paste0(input_folder,"ICGC_summary.RData"))
  post_summary <- get(x)
  
  if (!dir.exists(output)) dir.create(output)
  

  for (i in 1:length(cancertype)){
    
    library(NMF)
    library(dplyr)
    
    # get sample list for each type
    type <- cancertype[i]
    
    print(paste0("finish load ",i,'th type - ',type))
    
    samplelist <- subset(post_summary,Types==type)$samplename
    n_sample <- length(samplelist)
    
    ccfBandCountsMat <- CountMatBuild(samplelist)
    
    if (!dir.exists(paste0(output,type,"/"))) dir.create(paste0(output,type,"/"))  
    
    write.csv(samplelist,file=paste0(output,type,"/",type,"_ICGC_", n_sample,"_samplelist.csv"))
    
    ccfCountsMat <- ccfBandCountsMat[[1]]
    ccfFractionMatrix <- ccfBandCountsMat[[2]]
    ccfFractionRandomMatrix <-ccfBandCountsMat[[3]]
    
    save( ccfCountsMat,file=paste0(output,type,"/",type,"_ICGC_", n_sample,"_0-1_countMatrix.RData"))
    write.csv( ccfCountsMat,file=paste0(output,type,"/",type,"_ICGC_",n_sample,"_0-1_countMatrix.csv"))
    write.csv(ccfFractionMatrix ,file=paste0(output,type,"/",type,"_ICGC_",n_sample,"_0-1_fractionMatrix.csv"))
    write.csv(ccfFractionRandomMatrix,file=paste0(output,type,"/",type,"_ICGC_",n_sample,"_0-1_fractionMatrix_random.csv"))
  }
}

### --- Step3: Input rank estimate csv and plot ----- ####
rank_estimate_plot <- function(rank_folder,output,cancertype){
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  
  for (i in 1:length(cancertype)){
    
    type <- cancertype[i]

    estimate <- read.csv(file=paste0(rank_folder,type,".csv"))
    estimate_random <- read.csv(file=paste0(rank_folder,type,"_random.csv"))
    
    estimate$type <- 'normal'
    estimate_random$type <- 'random'
    estimate_rank <- rbind(estimate,estimate_random)
    rss_decrease <-  (estimate[order(estimate$rank),]$rss[1:10]-estimate[order(estimate$rank),]$rss[2:11]) - (estimate_random[order(estimate_random$rank),]$rss[1:10]-estimate_random[order(estimate_random$rank),]$rss[2:11])
    write.csv(estimate_rank,paste0(output,type,'_estimate_rank_',Sys.Date(),'.csv'))
    write.csv( rss_decrease,paste0(output,type,'_rss_decrease_',Sys.Date(),'.csv'))
    
    pdf(paste0(output,type,'_estimate_rank_',Sys.Date(),".pdf"),width=11,height=8)
    p1 <- ggplot(data=estimate_rank,aes(x=rank,y=cophenetic,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
    p2 <- ggplot(data=estimate_rank,aes(x=rank,y=dispersion,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
    p3 <- ggplot(data=estimate_rank,aes(x=rank,y=evar,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
    p4 <- ggplot(data=estimate_rank,aes(x=rank,y=rss,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
    p5 <- ggplot(data=estimate_rank,aes(x=rank,y=euclidean,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
    p6 <- ggplot(data=estimate_rank,aes(x=rank,y=kl,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "bottom",axis.title.x=element_blank())
    
    g_legend<-function(a.gplot){
      tmp <- ggplot_gtable(ggplot_build(a.gplot))
      leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
      legend <- tmp$grobs[[leg]]
      return(legend)}
    
    mylegend<-g_legend(p6)
    
    print(grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6+theme(legend.position = "none"),nrow=2),mylegend,nrow=2,top=paste0("NMF rank estimate for ",type," (rss_decrease: ",rss_decrease,")"),heights=c(10, 1)))
    dev.off()
    print(paste0("finish ",i," th file"))
  }
}

### --- Step4: Run NMF based on chosen rank and output sig/expo ----- 

# only plot for choosen types
                          nmf_sig_plot <- function(type,MatType="count",input_folder,output,rank_summary){
  
  
     file_path <- dir(paste0(input_folder,type,"/"))[grep(paste0(MatType,"Matrix.csv"),dir(paste0(input_folder,type,"/")))]
     samplename_path <- dir(paste0(input_folder,type,"/"))[grep("samplelist",dir(paste0(input_folder,type,"/")))]
     if (!dir.exists(paste0(output,type))) dir.create(paste0(output,type))
       
     CountMat <- read.csv(file=paste0(input_folder,type,'/',file_path))[,-1]
     samplename <- as.character(read.csv(file=paste0(input_folder,type,'/',samplename_path))[,-1])
     n_sample <- length(samplename)
     rank <- subset(rank_summary,TCGA.Study == type)$rank
     
     index <- sample(1:ncol(CountMat))
     CountMat <- CountMat[,index]
     samplename_random <- samplename[index]
     
     #preprocess for rows with all 0
     index_p <- which(rowSums(CountMat)>0)
     index_n <- which(!rowSums(CountMat)>0)
     CountMat <- CountMat[which(rowSums(CountMat)>0),]
     
     #run NMF
     res <- nmf(CountMat,rank,.opt='vp4')
     sig <- as.data.frame(res@fit@W)
     sig <- rbind(sig,matrix(0,nrow=length(index_n),ncol=ncol(sig)))
     sig <- sig[c(index_p,index_n),]
     
     expo <- res@fit@H
     expo <- as.matrix(expo)
     
     #sig plot
     sig1 <- as.data.frame(t(apply(sig,2,function(x) x/sum(x))))
     colnames(sig1) <- 1:ncol(sig1)
     sig1$sig <- paste0("Signature ",1:rank)
     xx <- melt(sig1,id=c("sig"))
     
     fills <- brewer.pal(8, "Set3")[c(1,3:8,2)]
     
     pdf(file=paste0(output,type,"/",type,"_",n_sample,"_sig_plot_",MatType,"_",Sys.Date(),".pdf"),height=3,width=12)
     
     p1 <- ggplot(xx,aes(y=value,x=variable)) + geom_bar(aes(fill=sig),stat='identity') + scale_fill_manual(values = fills)+ theme_grey()+
       ggtitle(paste0("rank = ",rank,", cancertype = ",type, ", MatrixType = ",MatType )) + 
       theme(legend.title = element_blank(),
                legend.position = "none",
                strip.text.x = element_text(color= "white",size=10),
                panel.grid.minor.y = element_blank(),
                axis.title.x = element_text(color = "grey20", size = 8),
                axis.text.x = element_text(color = "grey20", size = 6),
                axis.text.y = element_text(color = "grey20", size = 6),
                plot.title = element_text(size= 10)) +
       facet_grid(cols = vars(sig))+ scale_x_discrete(breaks=c("1","50","100") ,labels=c("0", "0.5", "1"))+xlab("Cancer Cell Fraction") + ylab("") 
       
      
     # theme(strip.background = element_blank(), strip.text  = element_text(color = 'white',size=8),
     #       legend.title=element_blank(),legend.text=element_blank(),legend.position = "none",plot.title = element_text(size=10),
     #       axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
     #       axis.text.x = element_text(color = "grey20", size = 6),
     #       axis.text.y = element_text(color = "grey20", size = 6))
     
     
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
     
     dev.off()
     #output sig and expo
     expo <- as.data.frame(t(expo))
     expo <- cbind(expo,samplename_random)
     colnames(expo)[1:rank] <- paste0("sig_",1:rank)
     
     save(expo,file=paste0(output,type,'/',type,'_',n_sample,"_expowithsample_",Sys.Date(),".RData"))
     save(sig,file=paste0(output,type,'/',type,'_',n_sample,"_sig_",Sys.Date(),".RData"))
     save(res,file=paste0(output,type,'/',type,'_',n_sample,"_res_",Sys.Date(),".RData"))
     
     return(list(g1,CountMat))
  }
                                   
# overall plot for all given type
nmf_sig_all_plot <- function(cancertype,input_folder,output,rank_summary) {
  library(gridExtra)
  library(ggplot2)
  library(RColorBrewer)
  library(grid)
  library(dplyr)
  library(NMF)
  library(reshape2)
  library(cowplot)
  
  i <- 1
  for (i in 1:length(cancertype)) {
    type <- cancertype[i]
    
    nmf_count <- nmf_sig_plot(type,input_folder=input_folder,output=output,rank_summary=rank_summary)
    p_count <- nmf_count[[1]]
    CountMat <- nmf_count[[2]]
    n_sample <- ncol(CountMat)
    p_fraction <- nmf_sig_plot(type,MatType="fraction",input_folder=input_folder,output=output,rank_summary=rank_summary)[[1]]
    
    p <- plot_grid(p_count,p_fraction, align="V",ncol=2)
    table <- t(as.data.frame(round(quantile(colSums(CountMat),c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)),2)))
    table1 <- tableGrob(table,rows=NULL,theme=ttheme_default(colhead=list(fg_params = list(parse=TRUE)),base_size = 8))
    title <- ggdraw() + draw_label(paste0("Extracting Signatures for ",type," with ",ncol(CountMat)," samples "),fontface='bold')
  
    commands <- paste0("p",i," <- plot_grid(title,table1,p,align='V',ncol=1,rel_heights = c(0.1,0.15,0.65))")
    eval(parse(text=commands))
    
    }
    
    n_pdf <- (length(cancertype) %/% 8)+1
    
    for (i in 1:n_pdf){
      commands1 <- paste0("pdf(file=paste0(output,'Sig_summary_p",(8*i-7),"-",8*i,"_',Sys.Date(),'.pdf'),height=25,width=12)")
      if (i != n_pdf) commands2 <- paste0("g",i,"<- plot_grid(",paste0("p",(8*i-7):(8*i),collapse=","),",align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))") else
      {
        n_file <-  length(cancertype) %% 8
        commands2 <- paste0("g",i," <- plot_grid(",paste0("p",(8*i-7):(8*i-8+n_file),collapse=","),",align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))")
        paste0("p",(8*i-7):(8*i-8+n_file),seq=",")
      }
      
      commands3 <- paste0("print(g",i,")")
      commands4 <- "dev.off()"
      eval(parse(text=commands1))
      eval(parse(text=commands2))
      eval(parse(text=commands3))
      eval(parse(text=commands4))
    }
    
    return(list(g1,g2,g3,g4))
}
                                   
### --- Step5: Association Analysis ----- 
