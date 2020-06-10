#' Plotting rank estimate
#' @name rank_estimate_plot
#' @param outputFolder folder stores rank estimate files
#' @param rankfilepath rank file path to output
#' @return save nmf results in output folder and plot signature for all cancer types
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @import dplyr
#' @import ggplot2
#' @import reshape2
rank_estimate_plot <- function(outputFolder,rankfilepath,format) {
  
  typelist <- unique(unlist(lapply(dir(outputFolder),function(x) strsplit(x,"_")[[1]][[1]])))
  rank_blank <- data.frame(cancertype=typelist,rank=NA,rss_suggested_rank=NA)
  
  i = 1
  for (i in 1:length(typelist)) {
    tryCatch({
        type <- typelist[i]
        if (format=="fraction"){
          estimate <- read.csv(paste0(outputFolder,type,"_ccfFractionMatrix.csv")) %>% mutate(Data='normal')
          estimate_random <- read.csv(paste0(outputFolder,type,"_ccfFractionMatrix.random.csv")) %>% mutate(Data='random')
        }
        
        if (format=="count"){
          estimate <- read.csv(paste0(outputFolder,type,"_ccfCountMatrix.csv")) %>% mutate(Data='normal')
          estimate_random <- read.csv(paste0(outputFolder,type,"_ccfCountMatrix.random.csv")) %>% mutate(Data='random')
        }
        
        estimate_rank <- rbind(estimate,estimate_random) %>% filter(rank>2)
        xx <- reshape2::melt(estimate_rank,id=c("rank","Data")) %>% mutate(Measure=NA)
        xx$variable <- as.character(xx$variable)
       
        xx[which(xx$variable=="dispersion"),]$Measure <- "Best fit"
        xx[which(xx$variable=="cophenetic"),]$Measure <- "Consensus"
        xx[which(xx$variable=="sparseness1"),]$Measure <- "Basis"
        xx[which(xx$variable=="sparseness2"),]$Measure <- "Coefficients"
        xx[c(which(xx$variable=="sparseness2"),which(xx$variable=="sparseness1")),]$variable <- "sparseness"
        xx[which(xx$variable=="rss"),]$Measure <- "Consensus"
        xx <- subset(xx,variable %in% c("cophenetic","dispersion","rss","sparseness"))
        
        min_rank = min(xx$rank)
        
        idx = 1:(nrow(estimate)-1)
        rss_decrease <-  min_rank + which((estimate[order(estimate$rank),]$rss[idx]-estimate[order(estimate$rank),]$rss[idx+1]) - (estimate_random[order(estimate_random$rank),]$rss[idx]-estimate_random[order(estimate_random$rank),]$rss[idx+1])<0)[1]
        rank_blank[which(rank_blank$cancertype == type),]$rss_suggested_rank <- rss_decrease 
         
        write.csv(estimate_rank,paste0(outputFolder,type,'_rank_summary.csv'))
        
        command1 <- paste0('g',i,'<-ggplot(data=xx,aes(x=rank,y=value))+ geom_line(aes(lty=Data,color=Measure),cex=1)+geom_point(aes(shape=Data,color=Measure),cex=2)+facet_wrap(~variable,scales="free",nrow=1)+
          labs(title=paste0("Rank estimate for ",type),subtitle=paste0("- rss suggested rank = ",rss_decrease))+
          scale_x_continuous(breaks = min(xx$rank):max(xx$rank))+theme_grey()+ theme(strip.background = element_rect(fill="orange"),strip.text = element_text(colour ="white",size=14))')
        
        command2 <- paste0("print(g",i,")")
        eval(parse(text=command1))
        eval(parse(text=command2))
      
    },error=function(e) print("error!"))
  }
    if (exists("g1")) {
    # output signature plot for all cancer types
    j <- length(typelist) 
    n_file <-  j %% 8
    
    if (n_file==0) {
      n_pdf <- (j %/% 8) 
    } else {
      n_pdf <- (j %/% 8)+1
    }
    
    for (i in 1:n_pdf) {
      if (i != n_pdf) {
        commands1 <- paste0("pdf(file=paste0(outputFolder",",'rank_estimate_g',(8*i-7),'-',8*i,'_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("p",i,"<- plot_grid(paste0('g',(8*i-7):(8*i),collapse=','),align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))") 
      } else {
        commands1 <- paste0("pdf(file=paste0(outputFolder",",'rank_estimate_g',(8*i-8),'-',8*i-8+n_file,'_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("p",i," <- plot_grid(paste0('g',(8*i-8):(8*i-8+n_file),collapse=','),align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))")
      }
      commands3 <- paste0("print(p",i,")")
      commands4 <- "dev.off()"
      
      eval(parse(text=commands1))
      eval(parse(text=commands2))
      eval(parse(text=commands3))
      eval(parse(text=commands4))
    }
  write.csv(rank_blank,file=rankfilepath)
  }
}
