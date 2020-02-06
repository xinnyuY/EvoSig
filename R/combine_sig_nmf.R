combine_sig_nmf <-
function(input_folder,cancertype){
  
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
