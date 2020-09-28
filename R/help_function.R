#' Correlation calculation function 
#' @name cor.table
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @return A table including correlation r,p,n between variables A and variables B
#' @import dplyr
#' @import reshape2
#' @export
cor.table = function(df,va,vb){
  
  df[,c(vb)] <- apply(df[,c(vb)],2,function(x) x=replace(x,x==0,NA)) # avoid 0 and zero-deviation
  
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

#' Plot scatterplot for specific measure with evo sig and subtype - 20200915
#' @name p_scatter
#' @param df exposure file merged with measures of interest
#' @param sig Column indexs of evo signature
#' @param vb Column indexs of all measures in association heatmap in df
#' @param facet name of facet variable consistent with association heatmap
#' @return Scatterplots
#' @import dplyr
#' @import ggpubr
#' @import scales
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
p_scatter <- function(df,vb,var,sig,facet) {
  
  var_name = colnames(df)[var]
  n_facet = nchar(facet)+2
  facet_name = paste0(facet," ")
  
  cor_table = cor_facet(df=df,va=sig,vb=vb,facet=facet) %>%
    mutate(p_label=paste0("R =",r," ,p = ",round(adj.p,4)),significant=as.factor(significant),
           facet=as.character(facet),Var1=as.character(Var1)) %>%
    mutate(facet=gsub(facet_name,"",facet)) %>%
    filter(Var2==var_name) 
  
  if (length(sig)==1) cor_table$Var1 = colnames(df)[sig]
  
  colnames(df)[var]="Var"
  
  df_new <- df %>%
    melt(id=colnames(df)[-sig]) %>%
    dplyr::rename(Var1=variable) %>%
    drop_na(Var) %>%
    mutate(Var=as.numeric(Var),facet=as.character(facet),Var1=as.character(Var1)) %>%
    filter(value>0,Var>0) %>%
    left_join(cor_table,by=c("Var1","facet")) %>%
    left_join(subset(cor_table,facet=="All") %>% 
                dplyr::rename(all_p_label= p_label,all_significant=significant,all_label=label) %>% .[,c(2,7:9)],by="Var1")
  
  fun_median_y <- function(x){
    return(data.frame(y=median(x),label=round(median(x,na.rm=T),2)))}
  
  ymin = min(df_new$Var); ymax=max(df_new$Var);ymean <- mean(df_new$Var)
  xmax = max(df_new$value); xmin=min(df_new$value);xmean=mean(df_new$value)
  
  p1 <- ggplot(df_new)+
    scale_x_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
    scale_y_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
    geom_point(aes(x=value,y=Var,colour = significant))+
    geom_smooth(method = "lm",aes(value,Var,color=significant))+
    geom_boxplot(aes(x=0,y=Var),width=xmax/100,varwidth=TRUE)+  # y_axis boxplot
    geom_boxplot(aes(x=value,y=ymin),orientation="y",width=ifelse(which.min(c(ymax/120,ymin))==1,ymax/120,ymin/3))+ # x_axis boxplot
    facet_grid(cols=vars(label),rows=vars(Var1),scale="free")+
    labs(y="",x="")+
    theme_pubclean()+
    geom_text(aes(y=ymax*0.75,x=xmean*1.25,label=p_label,color=significant),size=4,check_overlap=TRUE)+
    scale_colour_manual(values=c("grey","#1A9993FF"))+
    theme(legend.position = "right",plot.margin = unit( c(0.2,0,0,0) , units = "lines" ) )
  
  
  p_all <- df_new %>% ggplot()+
    scale_x_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
    scale_y_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
    geom_point(aes(x=value,y=Var,color=all_significant))+
    geom_smooth(method = "lm",aes(value,Var,color=all_significant))+
    geom_boxplot(aes(x=0,y=Var),width=xmax/100,varwidth=TRUE)+   # y_axis boxplot
    geom_boxplot(aes(x=value,y=ymin),orientation="y",width=ifelse(which.min(c(ymax/120,ymin))==1,ymax/120,ymin/3))+ # x_axis boxplot
    facet_grid(cols=vars(all_label),rows=vars(Var1),scale="free")+
    labs(y=var_name,x=" ")+
    theme_pubclean()+
    geom_text(aes(y=ymax*0.75,x=xmean*1.25,label=all_p_label,color=all_significant),size=4,check_overlap=TRUE)+
    scale_colour_manual(values=c("grey","#1A9993FF"))+
    theme(strip.text.y = element_blank() , 
          strip.background.y = element_blank(),
          plot.margin = unit( c(0.2,0.2,0,0) , units = "lines" ),
          legend.position = "none")
  
  p <- grid.arrange(p_all,p1,nrow=1,widths=c(1,5),bottom="Evolutionary Dynamics Signature Exposure")
  return(p)
}
#' Correlation calculation function by group and plot association heatmap - 20200915
#' @name cor_facet
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @param facet name of group variable
#' @param heatmap specify whether plot association heatmap
#' @param title Title for association heatmap
#' @param empty_row_delete whether to delete rows without any significant correlation value
#' @param flip whether to flip facet and x
#' @param keep_all whether to show all facet
#' @param col_low color for correlation r=-1
#' @param col_high color for correlation r=1
#' @return A table including correlation r,p,n between variables A and variables B during group variable and association heatmap (optional)
#' @import dplyr
#' @importFrom plyr ddply
#' @import reshape2
#' @import ggpubr
#' @import ggplot2
#' @export
cor_facet = function(df,va,vb,facet,heatmap=FALSE,title="",empty_row_delete=FALSE,flip=FALSE,keep_all=TRUE,col_low="#4a7b94",col_high="#bb5a39"){
  facet_idx=which(colnames(df)==facet)
  
  vb = vb[which(apply(df[,vb],2,function(x) length(unique(x))==1)==FALSE)]
  
  if (length(vb)!=0) {
  
  eval(parse(text=paste0("cor_table_all <- cor.table(df,va=va,vb=vb) %>% mutate(",facet,"='All') %>% dplyr::rename(facet=",facet,")")))
  
  command2 <- paste0("cor_table <- ddply(df, .(",facet,") , .fun =cor.table,va=va,vb=vb) %>% mutate(",facet,"=paste0('",facet," ',",facet,")) %>% dplyr::rename(facet=",facet,") %>% rbind(cor_table_all)") 
  eval(parse(text=command2))
  
  cor_table <- cor_table %>%
    mutate(r=replace(r,is.na(r),0),adj.p=replace(adj.p,is.na(adj.p)|is.nan(adj.p),1)) %>% # set correlation with na or p>0.05 invisible
    mutate(r=round(r,2),Var2=as.character(Var2),label=paste0(facet," \n (n=",n,")"))
  
  if (heatmap==TRUE) {
    cor_table_plot <- cor_table %>% mutate(r=replace(r,adj.p>0.05|n<=10,0)) 
    
    # delete empty rows
    if (keep_all==FALSE)  cor_table_plot <- subset( cor_table_plot,facet!="All")
    
    if (empty_row_delete==TRUE) {
      non_empty_gene <- cor_table_plot %>% group_by(Var2) %>% 
        dplyr::summarise(empty_r=all(r==0)) %>% filter(empty_r==FALSE) %>% as.data.frame() %>% .[,1] %>% as.character()
      
      cor_table_plot <- subset(cor_table_plot,Var2 %in% non_empty_gene)
    }
    
    if (flip==TRUE) {
      
      n_facet <- length(non_empty_gene)
      
      p_heatmap <- cor_table_plot %>% 
        ggplot(aes(y=Var2,x=facet,fill=r)) +
        geom_tile()+ geom_text(aes(y=Var2,x=facet,label=r),size=4,col='#ffffff')+
        facet_grid(cols=vars(Var1),switch="y",scale="free",space="free")
      
    } else {
      
      p_heatmap <- cor_table_plot %>% ggplot(aes(y=Var2,x=Var1,fill=r)) +
        geom_tile()+ geom_text(aes(y=Var2,x=Var1,label=r),size=4,col='#ffffff')
      
      if (length(unique(cor_table_plot$n))<=5) {
        
        p_heatmap <- p_heatmap + facet_grid(cols=vars(label),switch="y",scale="free",space="free")
        
      } else {
        
        n_var1 <- length(unique(cor_table_plot$Var1))
        p_heatmap <- p_heatmap + 
          geom_text(aes(y=Var2,x=n_var1+1,label=n),size=3,col="grey50")+
          facet_grid(cols=vars(facet),switch="y",scale="free",space="free")+ 
          coord_cartesian( xlim=c(1,n_var1+0.5),clip = "off")
      }
    }
 
    p_heatmap <- p_heatmap +
      scale_fill_gradient2(low=col_low,mid="#ffffff",high=col_high,midpoint = 0,name="Spearman Correlation")+
      labs(x="",y="",title=title)+ theme_pubclean()+
      theme(axis.text.x = element_text(angle=90,vjust=0.5),plot.margin = unit(c(1, 7, 1, 1), "lines"),legend.position = "right")
    
    p_heatmap
    ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,non_empty_gene)),return(list(cor_table,p_heatmap)))
    
  } else return(cor_table)
  
  } else {print("the standard deviation is zero")}
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


