#' Perform Hierarchical clustering number estimate
#' @name hc_cluster_test ccfMat ccf matrix for all samples
#' @param data combined signatures for all cancer types
#' @param methods clustering method
#' @param distance distance funciton
#' @param min minimum clustering number
#' @param max maximum clustering number
#' @return exposure
#' @import NbClust
hc_cluster_test <- function(data,methods,distance,min = 2,max = 10){

  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  tabla = as.data.frame(matrix(ncol = length(distance), nrow = length(methods)))
    names(tabla) = distance
  
  for (j in 2:length(distance))
      for(i in 1:length(methods)){
        tryCatch({
        nb = NbClust(data,distance = distance[j],
                     min.nc = min, max.nc = max, 
                     method = "complete", index =methods[i])
        tabla[i,j] = nb$Best.nc[1]
        tabla[i,1] = methods[i]
      },error=function(e) print("error"))
    } 
  
  tabla <- rbind(tabla,c("Most_common",apply(tabla[,2:5],2,getmode)))
  
  return(tabla)
}
