MAC_perm <-
function(data, MACs_observ, num_perms=100, max_lag_prop=1/3, FDR_cutoffs=101, perm_file_name=F){

  MACs_perm = c()

  # simplified MAC_counter function
  

  for(n in (1:num_perms)){

    set.seed(n)
    samp_size=min(100, nrow(data))
    data_perm = data[1:samp_size,]
    inds = sample(1:nrow(data),samp_size)
    
    for(z in (1:samp_size)){

      data_perm[z,] = sample(data[inds[z],])
    }
    
    results=MAC_lags(data=data_perm, max_lag_prop=1/3, symmetric=F)
    results_MACs=results[,1:nrow(data_perm)]
    MACs_perm = c(MACs_perm, results_MACs)

  }


### Calculate FDR ###

  cors = seq(0,1,length.out=FDR_cutoffs)

  num.cors.perm = rep(0, FDR_cutoffs)
  MACs_observed = rep(0, FDR_cutoffs)

  for(r in (1:FDR_cutoffs)){

    num.cors.perm[r] = num.cors.perm[r] + length(which(MACs_perm >= cors[r]))

    MACs_observed[r] = MACs_observed[r] + length(which(MACs_observ >= cors[r]))

  }

  perm.size=length(which(is.na(results_MACs)==F))
  obs.size=length(which(is.na(MACs_observ)==F))
  
  MACs_ave_perm = num.cors.perm/num_perms*(obs.size/perm.size)

  fdr = rep(NA,FDR_cutoffs)

  for(s in (1:FDR_cutoffs)){

    fdr[s] = MACs_ave_perm[s]/MACs_observed[s]

  }

  fdr[is.na(fdr)] = 0

  results = cbind(cors, MACs_observed, MACs_ave_perm, fdr)
  results_rev = results[FDR_cutoffs:1,]

  if(perm_file_name!=F){
    perm_name = paste(c("perm_", perm_file_name, ".csv"), sep="", collapse="")
    write.csv(results_rev, file = perm_name, row.names=F)
    }else{return(results_rev)}

}
