MAC_lags <-
function(data, max_lag_prop=1/3, symmetric=F){
  
  ## setup lag info ##
  
  num_time_int = ncol(data)
  max_lag = floor(num_time_int*(max_lag_prop))
  window = num_time_int-max_lag
  num_genes = nrow(data)
  
  ## Calculate Mean Vector and Covariance Matrix for Lag-0 ##
  
  data_0 = data[,c(1:window)]
  means_0 = rowMeans(data_0)
  cent_0 = t(data_0-means_0)
  Cor_0 = as.matrix(cor(cent_0))
  
  group = rep(1, window)
  rowsumx_0 = rowsum(cent_0, group)
  rowsumx2_0 = rowsum(cent_0^2, group)
  
  ## Calculate Cross-Correlation for Lags ##
  
  max_lag_cor = Cor_0
  lag_max = matrix(0,num_genes,num_genes)
  
  for(i in (1:max_lag)){
    
    data_lag = data[,c((1+i):(window+i))]
    means_lag = rowMeans(data_lag)
    cent_lag = t(data_lag-means_lag)

    rowsumx_lag = rowsum(cent_lag, group)
    rowsumx2_lag = rowsum(cent_lag^2, group)
    
    Cor_lag = (t(cent_0)%*%cent_lag - (1/window)*t(rowsumx_0)%*%rowsumx_lag)/(t(sqrt(rowsumx2_0-rowsumx_0^2/window))%*%sqrt(rowsumx2_lag-rowsumx_lag^2/window))
    
    # Keep cor only if greater than previous lag cor
    
    ind=which(abs(max_lag_cor)<abs(Cor_lag))
          
    max_lag_cor[ind] = Cor_lag[ind]
          
    lag_max[ind] = i
        
  }
  
  MACs_greatest = max_lag_cor
  lag_greatest = lag_max
  
  if(symmetric==T){
    
    max_lag_cor_t = t(max_lag_cor)
    lag_max_t = t(lag_max)
    ind2 = which(abs(max_lag_cor)>=abs(max_lag_cor_t))
    MACs_greatest[ind2] = max_lag_cor[ind2]
    lag_greatest[ind2] = lag_max[ind2]
    MACs_greatest[-ind2] = max_lag_cor_t[-ind2]
    lag_greatest[-ind2] = lag_max_t[-ind2]
    
  }
  
  diag(MACs_greatest) = NA
  diag(lag_greatest) = NA
  
  results = cbind(MACs_greatest, lag_greatest)
  
  return(results)
  
}
