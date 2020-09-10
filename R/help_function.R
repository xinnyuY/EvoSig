#' Correlation calculation function 
#' @name cor
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @return A table including correlation r,p,n between variables A and variables B
#' @import dplyr
#' @import reshape2
cor = function(df,va,vb){
  res = psych::corr.test(df[,va],df[,vb],method = "spearman",use = "pairwise",adjust="holm",ci=FALSE) 
  if (!is.null(res)) {
    r = melt(res['r']) %>% dplyr::rename(r = value) %>% dplyr::select(-L1)
    p = melt(res['p']) %>% dplyr::rename(adj.p = value) %>% dplyr::select(-L1)
    n = melt(res['n']) %>% dplyr::rename(n = value) %>% dplyr::select(-L1)
    
    if (ncol(n)>1) {
      left_join(left_join(r,p,by=c("Var1","Var2")),n,by=c("Var1","Var2")) %>% mutate(significant=ifelse(adj.p<0.05,1,0))
    } else {
      left_join(r,p,by=c("Var1","Var2")) %>% mutate(n=as.numeric(n),significant=ifelse(adj.p<0.05,1,0))
    }
  } else {rep(NA,5)}
}

#' Correlation calculation function by group and plot association heatmap
#' @name cor_facet
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @param facet name of group variable
#' @param heatmap specify whether plot association heatmap
#' @return A table including correlation r,p,n between variables A and variables B during group variable and association heatmap (optional)
#' @import dplyr
#' @import plyr
#' @import reshape2
#' @import ggpubr
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
cor_facet = function(df,va,vb,facet,heatmap=FALSE){
  coul = colorRampPalette(brewer.pal(8, "RdBu"))(256)
  facet_idx=which(colnames(df)==facet)

  command1 <- paste0("cor_table_all <- cor(df,va=va,vb=vb) %>% mutate(",facet,"='All')")
  eval(parse(text=command1))
  command2 <- paste0("cor_table <- ddply(df, .(",facet,") , .fun =cor,va=va,vb=vb) %>% mutate(",facet,"=paste0('",facet," ',",facet,")) %>% rbind(cor_table_all) %>% mutate(label=paste0(",facet,",' \n (n=',n,')'))") 
  eval(parse(text=command2))
  
  cor_table_plot <- cor_table %>%
    mutate(r=replace(r,(is.na(r)|adj.p>=0.05),0)) %>% # set correlation with na or p>0.05 invisible
    mutate(r=round(r,2),Var1=paste0("Signature ",as.numeric(Var1)),Var2=as.character(Var2))
  
  # delete empty rows
  non_empty_gene <- cor_table_plot %>% group_by(Var2) %>% dplyr::summarise(mean_r=mean(r)) %>% filter(mean_r!=0) %>% as.data.frame() %>% .[,1] %>% as.character()
  
  cor_table_plot <- subset(cor_table_plot,Var2 %in% non_empty_gene)
  
  if (heatmap==TRUE) {
    commands_p <- paste0("p <- cor_table_plot %>% ggplot(aes(y=Var2,x=",facet,",fill=r)) +
         geom_tile()+ geom_text(aes(y=Var2,x=",facet,",label=r),size=4,col='#ffffff')")
    eval(parse(text=commands_p))
    
    p_heatmap <- cor_table_plot %>% ggplot(aes(y=Var2,x=Var1,fill=r)) +
      geom_tile()+ geom_text(aes(y=Var2,x=Var1,label=r),size=4,col='#ffffff')+
      facet_grid(cols=vars(label),switch="y",scale="free",space="free")+ 
      scale_fill_gradient2(low=coul[256],mid="#ffffff",high=coul[1],midpoint = 0,name="Spearman Correlation") +
      labs(x="",y="")+ theme_pubclean()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5))
    p_heatmap
    return(list(cor_table,p_heatmap))
  } else return(cor_table)
}



#' Load multiple format ccf file
#' @name load_ccf
#' @param samplename cancer type
#' @param input 
#' @export
#' @return ssm
load_ccf <- function(samplename,input){
  Check <- ArgumentCheck::newArgCheck()
  suppressWarnings(rm(ssm,res,ccubeRes))
  
  format1 <- paste0(input,samplename,"/ccube_result.RData")
  format2 <- paste0(input,samplename,"/ccube_res_v0.3_final_new.RData")
  format3 <- paste0(input,samplename,"/ccube_res_run1.RData")
  format4 <- paste0(input,samplename,"/ccubeRes.RData")
  
  if (file.exists(format4 )) load(format4) else
    if (file.exists(format1 )) load(format1) else
     if(file.exists(format2)) load(format2) else  
      if(file.exists(format3)) load(format3) else{
          ArgumentCheck::addError(
          msg = "No file has been loaded",
          argcheck = Check)
          }
  if (exists("ssm")) return(ssm) else
    if (exists("ccubeRes")) return(ccubeRes$ssm) 
}

#' create multiple dir
#' @name multi.dir.create
#' @param list list of directory
#' @return create multiple directory
#' @export
multi_dir_create <- function(dirlist){
  for (i in dirlist) {if (!dir.exists(i)) dir.create(i,recursive = T)}
}

#' Unify format of data frame 
#' @name file_format
#' @param filename data frame
#' @param samplenamecol column index of samplename
#' @export
file_format <- function(filename=filename,samplenamecol){
    names <- colnames(filename)
    names[samplenamecol] <- "samplename"
    names -> colnames(filename)
    filename$samplename <- substr(filename$samplename,1,12)
    filename$samplename <- gsub("[.]","-",filename$samplename)
    return(filename)
  }

#' Merge ssm and cna files
#' @name ParseSnvCnaPcawgFormat
#' @param ssm ssm
#' @param cna cna
#' @export
#' @import dplyr
ParseSnvCnaPcawgFormat <- function (ssm, cna) {
    
    ssm <- ssm %>%
      mutate(chr = substr(chr,4,length(chr)),
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

#' Customize top strip color
#' @name strip_color
#' @param p plot
#' @param col customized color
#' @param draw whether to display in plots panal
#' @param direction strip location
#' @import ggplot2 
#' @import grid
#' @export
strip_color <- function(p,col=signature_col,draw=FALSE,direction='top'){
  p1 <- ggplot_gtable(ggplot_build(p))
  k <- 1
  
  if (direction=='top') strip_col <- which(grepl('strip-t', p1$layout$name))
  if (direction=='bottom') strip_col <- which(grepl('strip-b', p1$layout$name))
  if (direction=='left') strip_col <- which(grepl('strip-l', p1$layout$name))
  if (direction=='right') strip_col <- which(grepl('strip-r', p1$layout$name))
  
  for (j in strip_col) {
    j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
    p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col[k]
    k <- k+1
  }
  if (draw) grid.draw(p1)
  return(p1)
}


