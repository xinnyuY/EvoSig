#' Plot signature
#' @param sig signature matrix
#' @return mydata mutation data frame
#' @export


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

nmf_sig_plot_type <- function(type,MatType="fraction",input_folder,output,rank){
  if (!requireNamespace("NMF", quietly = TRUE)) install.packages("NMF")
  library(NMF)
  
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

nmf_sig_all_plot <-
function(input_folder,output,rank_summary) {
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
    
    rank = as.numeric(subset(rank_summary,cancertype == type)$rank)
    
    nmf_result <- nmf_sig_plot_type(type,input_folder=input_folder,output=output,rank=rank)
    
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