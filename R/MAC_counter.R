MAC_counter <-
function(data, max_lag_prop=1/3, MAC_cutoff=0.2, file_name=F, lag_matrix=T, symmetric=F){
  
  results=MAC_lags(data, max_lag_prop=1/3, symmetric=F)
  
  results_MACs=results[,1:nrow(data)]
  results_lags=results[,-(1:nrow(data))]
  
  # write files
  
  if(file_name!=F){
    if(symmetric==F){
    
      if(lag_matrix==T){
        lag_name = paste(c("lag_", file_name, ".csv"), sep="", collapse="")
        write.csv(results_lags, file=lag_name, row.names=F)
       }

      MAC_name = paste(c("MAC_", file_name, ".csv"), sep="", collapse="")
      write.csv(results_MACs, file=MAC_name, row.names=F)

    }else{
    
    # write files
    
      if(lag_matrix==T){
          lag_name = paste(c("lag_symmetric_", file_name, ".csv"), sep="", collapse="")
          write.csv(results_lags, file=lag_name, row.names=F)
      }

      MAC_name = paste(c("MAC_symmetric_", file_name, ".csv"), sep="", collapse="")
      write.csv(results_MACs, file=MAC_name, row.names=F)
    
    }
  }
  
  final_results = c()
  
  cors = sort(as.matrix(abs(results_MACs)), decreasing=T, na.last=NA)
  max_cors = cors[which(cors>=MAC_cutoff)]
  unique_cors = unique(max_cors)
  
  ind_cor = which(abs(results_MACs)>=0.2)
  ind_rc = which(abs(results_MACs)>=0.2, arr.ind=T)
  
  data_inds=cbind(abs(results_MACs[ind_cor]),results_MACs[ind_cor], results_lags[ind_cor],ind_rc[,1], ind_rc[,2])
  
  final_results=data_inds[order(data_inds[,1], decreasing = T),-1]
  
  colnames(final_results) = c("Correlation", "Lag", "Row gene index", "Column gene index")
  
  return(final_results)
  
}
