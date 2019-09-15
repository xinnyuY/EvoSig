file_format <- function(filename=filename,samplenamecol){

  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  return(filename)
  
}
