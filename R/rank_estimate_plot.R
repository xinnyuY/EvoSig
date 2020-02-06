rank_estimate_plot <-
function(outputFolder,rankfilepath) {
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
  library(dplyr)
  file <- unique(unlist(lapply(dir(outputFolder),function(x) strsplit(x,"_")[[1]][[1]])))
  rank_blank <- data.frame(cancertype=file,rank=NA,rss_suggested_rank=NA)
  

  for (i in 1:length(file)) {
    tryCatch({
        filename <- file[i]
        estimate <- read.csv(paste0(outputFolder ,filename,"_ccfFractionMatrix.csv")) %>% mutate(type='normal')
        estimate_random <- read.csv(paste0(outputFolder ,filename,"_ccfFractionMatrix.random.csv")) %>% mutate(type='random')
        
        print(paste0("plot for ",filename)) 
        estimate_rank <- rbind(estimate,estimate_random)
        
        rss_decrease <-  which((estimate[order(estimate$rank),]$rss[1:10]-estimate[order(estimate$rank),]$rss[2:11]) - (estimate_random[order(estimate_random$rank),]$rss[1:10]-estimate_random[order(estimate_random$rank),]$rss[2:11])<0)[1]
        rank_blank[which(rank_blank$cancertype == filename),]$rss_suggested_rank <- rss_decrease + 1
         
        write.csv(estimate_rank,paste0(outputFolder,filename,'_rank_summary.csv'))
        
        pdf(paste0(outputFolder,filename,'_rank_summary.pdf'),width = 12,height = 6)
        p1 <- ggplot(data=estimate_rank,aes(x=rank,y=cophenetic,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        p2 <- ggplot(data=estimate_rank,aes(x=rank,y=dispersion,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        p3 <- ggplot(data=estimate_rank,aes(x=rank,y=evar,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        p4 <- ggplot(data=estimate_rank,aes(x=rank,y=rss,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        p5 <- ggplot(data=estimate_rank,aes(x=rank,y=euclidean,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        p6 <- ggplot(data=estimate_rank,aes(x=rank,y=kl,group=type)) + geom_line(aes(color=type))+geom_point(aes(color=type))+theme(legend.position = "bottom",axis.title.x=element_blank())
        p7 <- ggplot(data=estimate_rank) + geom_line(aes(x=rank,y=sparseness1,group=type,color=type))+geom_line(aes(x=rank,y=sparseness2,group=type,color=type))+geom_point(aes(x=rank,y=sparseness2,group=type,color=type))+theme(legend.position = "none",axis.title.x=element_blank())
        
        g_legend<-function(a.gplot){
          tmp <- ggplot_gtable(ggplot_build(a.gplot))
          leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
          legend <- tmp$grobs[[leg]]
          return(legend)}
        
        mylegend <-g_legend(p6)
        
        print(grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6+theme(legend.position = "none"),nrow=2),mylegend,nrow=2,top=paste0("NMF rank estimate for ",filename),heights=c(10, 1)))
        dev.off()
    },error=function(e) print("error!"))
  }
  
  write.csv(rank_blank,file=rankfilepath)
}
