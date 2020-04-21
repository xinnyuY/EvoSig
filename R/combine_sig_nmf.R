
#' Combine signatures from different cancer type together
#' @name combine_sig_nmf
#' @param input_folder folder stores signatures for each type
#' @param cancertype cancer type list
#' @param output_folder output_folder
#' @export
#' @return matrix combining all signature matrix
combine_sig_nmf <- function(input_folder,output_folder=NA,cancertype){
  
  for (i in 1:length(cancertype)){
      tryCatch({
        type <- cancertype[i]
        
        cat(paste0("-> loading  CCF matirx for ",i,"th type : ",type),"\n")
        
        sig_file <- dir(paste0(input_folder,type,"/"))[grep("sig",dir(paste0(input_folder,type,"/")))]
        
        load(paste0(input_folder,type,"/",sig_file))
      
        colnames(sig) <- paste0(type,"_sig",1:ncol(sig))
        if (i==1) {
          combine_sig <- sig
        } else {
           combine_sig <- cbind(combine_sig,sig)
        }
      },error=function(e) print("Fail load this type"))
  }  
   
   if (!is.na(output_folder)) save(combine_sig,file=paste0(output_folder,"combine_sig_",Sys.Date(),".RData"))
   
   combine_sig
}
