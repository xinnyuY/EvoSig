load_ccf <-
function(samplename,input){
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

multi.dir.create <- function(list){
  for (i in list) {if (!dir.exists(i)) dir.create(i,recursive = T)}
}

file_format <-
  function(filename=filename,samplenamecol){
    names <- colnames(filename)
    names[samplenamecol] <- "samplename"
    names -> colnames(filename)
    filename$samplename <- substr(filename$samplename,1,12)
    filename$samplename <- gsub("[.]","-",filename$samplename)
    return(filename)
  }

ParseSnvCnaPcawgFormat <-
  function (ssm, cna) {
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

