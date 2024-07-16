
########FUNCTIONS 

glmnet_shuffeled_function1000x_get_all_1000R2<- function(forx,fory){
  nrun <- 1000  # no. of shuffling
  pb <- txtProgressBar(style=3)
  r2n <- rep(0,nrun)
  r2_s_all_pathwas_shuffeled<-vector(mode = "list") # this is list because each pathway will have 1000 repeats
  for (k in 1:ncol(fory)) { #this is taking each of pathways that we will estimate r2 for
    for(i in seq(nrun)){
      xi <- forx[sample(nrow(forx)),]
      f <- cv.glmnet(x=xi, y=fory[,k], family='binomial')  # linear regression fit
      idx <- f$index[2]                       # index of lambda.min and lambda.1se
      r2 <- f$glmnet.fit$dev.ratio[idx]    # R2 for the fit with two choices of lambda
      r2n[i] <- r2
      #setTxtProgressBar(pb,i/nrun)
    }
    
    r2_s_all_pathwas_shuffeled[[k]]<- r2n  # collecting all 1000 observations from shuffling
    
  }
  close(pb)
  r2_s_all_pathwas_shuffeled_df<-data.frame(do.call(rbind,r2_s_all_pathwas_shuffeled))   #these are now mean r2 values for 10 runs for each pathway 
  rownames(r2_s_all_pathwas_shuffeled_df)<-colnames(fory)
  Thusend_perm_means<-data.frame(rowMeans(r2_s_all_pathwas_shuffeled_df))
  return(r2_s_all_pathwas_shuffeled_df)   #these are values for 1000 runs for each pathway
} 


########### this is normal function , without shufflking for 1000 times 

glmnet_binary_function1000x<- function(forx,fory){
  nrun <- 1000  # no. of replicates for how many times we run cv.glmnet
  r2s_cellularity <- rep(0,nrun)
  r2_s_all_pathwas_cellularity<-vector(mode = "list") 
  for (k in 1:ncol(fory)) { #this is takimng each of pathways that we will estimate r2 for
    for(i in seq(nrun)){
      
      f <- cv.glmnet(x=forx, y=fory[,k], family='binomial')  # linear regression fit, nfold is always 10 here, default
      idx <- f$index[2]                       # index of lambda.min and lambda.1se
      r2 <- f$glmnet.fit$dev.ratio[idx]    # R2 for the fit with two choices of lambda
      r2s_cellularity[i] <- r2
    }
    r2_s_all_pathwas_cellularity[[k]]<- r2s_cellularity  # estimated r2 for the original data
    
  }
  return(r2_s_all_pathwas_cellularity)   #these are now 1000 r2 values for 10 runs for each pathway
} 


#######################################   FUNCTIONS FOR REPEATED GLMNET RUNS AND GET r2 values


glmnet_shuffeled_function1000x<- function(forx,fory){
  nrun <- 1000  # no. of shuffling
  pb <- txtProgressBar(style=3)
  r2n <- rep(0,nrun)
  r2_s_all_pathwas_shuffeled<-vector(mode = "list") # this is list because each pathway will have 1000 repeats
  for (k in 1:ncol(fory)) { #this is taking each of pathways that we will estimate r2 for
    for(i in seq(nrun)){
      xi <- forx[sample(nrow(forx)),]
      f <- cv.glmnet(x=xi, y=fory[,k], family='binomial')  # linear regression fit
      idx <- f$index[2]                       # index of lambda.min and lambda.1se
      r2 <- f$glmnet.fit$dev.ratio[idx]    # R2 for the fit with two choices of lambda
      r2n[i] <- r2
    #  setTxtProgressBar(pb,i/nrun)
    }
    
    r2_s_all_pathwas_shuffeled[[k]]<- r2n  # collecting all 1000 observations from shuffling
    
  }
  close(pb)
  r2_s_all_pathwas_shuffeled_df<-data.frame(do.call(rbind,r2_s_all_pathwas_shuffeled))   #these are now mean r2 values for 10 runs for each pathway 
  rownames(r2_s_all_pathwas_shuffeled_df)<-colnames(fory)
  Thusend_perm_means<-data.frame(rowMeans(r2_s_all_pathwas_shuffeled_df))
  return(Thusend_perm_means)   #these are now mean r2 values for 10 runs for each pathway
} 

######################################   FUNCTON TO RUN And EXGTRACT GLMNET RESULTS######################################  
Function_for_glmnet <- function(predictor_matrix,to_be_predicted_matrix) {
  
  
  glmnet_results<- vector(mode = "list")
  glmnet_results_r2<- vector(mode = "list")
  
  for (i in 1:ncol(to_be_predicted_matrix)) {
    
    f <- cv.glmnet(as.matrix(predictor_matrix),as.matrix(to_be_predicted_matrix[,i]), family='binomial')
    #plot(f) # display cv outcome
    cf <- as.matrix(coef(f))
    cf <- cf[cf!=0,,drop=F]
    out_table<-data.frame(print(cf))
    out_table$pathway<-colnames(to_be_predicted_matrix)[i]
    
    glmnet_results[[i]]<-out_table 
    
    #ALSO get r2 of the model to see how good is the model 
    idx <- f$index       # index of lambda.min and lambda.1se
    r2 <- data.frame(f$glmnet.fit$dev.ratio[idx])
    r2$pathway<-colnames(to_be_predicted_matrix)[i]
    glmnet_results_r2[[i]]<-r2
  }
  
  
  glmnet_results_df<-data.frame(do.call(rbind,glmnet_results))
  glmnet_results_r2_df<-data.frame(do.call(rbind,glmnet_results_r2))
  
   
  #we need to return list of this 2 tables
  return(list(glmnet_results_df = glmnet_results_df, glmnet_results_r2_df =glmnet_results_r2_df))
}
##################################################################################################################  

#usage of the function
#forx: matrix to be predicted; rows individuals, colm pathways or whatever need to be predicted 
#fory matrix which we have faetures in that we ask which once can predict x; rows individuals, colm features (genes)
#glmnet_shuffeled_function1000x(as.matrix(forx) ,as.matrix(fory))





######################################   FUNCTON TO RUN And EXTRACT GLMNET RESULTS######################################  
r <- function(predictor_matrix,to_be_predicted_matrix) {
  
  
  glmnet_results<- vector(mode = "list")
  glmnet_results_r2<- vector(mode = "list")
  
  for (i in 1:ncol(to_be_predicted_matrix)) {
    
    f <- cv.glmnet(as.matrix(predictor_matrix),as.matrix(to_be_predicted_matrix[,i]), family='binomial')
    #plot(f) # display cv outcome
    cf <- as.matrix(coef(f))
    cf <- cf[cf!=0,,drop=F]
    out_table<-data.frame(print(cf))
    out_table$pathway<-colnames(to_be_predicted_matrix)[i]
    
    glmnet_results[[i]]<-out_table 
    
    #ALSO get r2 of the model to see how good is the model 
    idx <- f$index       # index of lambda.min and lambda.1se
    r2 <- data.frame(f$glmnet.fit$dev.ratio[idx])
    r2$pathway<-colnames(to_be_predicted_matrix)[i]
    glmnet_results_r2[[i]]<-r2
  }
  
  
  glmnet_results_df<-data.frame(do.call(rbind,glmnet_results))
  glmnet_results_r2_df<-data.frame(do.call(rbind,glmnet_results_r2))
  
  
  #we need to return list of this 2 tables 
  return(list(glmnet_results_df = glmnet_results_df, glmnet_results_r2_df =glmnet_results_r2_df))
}
