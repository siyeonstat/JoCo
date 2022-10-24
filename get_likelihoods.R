
get_likelihoods <- function(data, stage, outcome, death, treatment, weight = 1, robust_int = F, isvar = F, varlist){
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  # (Equal) Weights of each individual for the weighted cox model
  weight_vec <- rep(weight, nrow(data))
  
  ###########################################
  ##### Variables used for the analysis #####
  ###########################################
  
  # Extract variables for the current analysis(loop) by identifying the stage number
  stage_end <- paste(".", stage, sep ="")
  varname <- names(data)[which(substrRight(colnames(data), 2) == stage_end )]
  
  # Put the stage number to the outcome, death and treatment
  # Get the list of main terms and interaction terms for the analysis
  outcome_name <- paste(outcome, stage_end, sep = "")  
  death_name <- paste(death, stage_end, sep = "")  
  treatment_name <- paste(treatment, stage_end, sep = "")  
  
  
  # Obtain the list of covariates
  if(isvar == T){
    
    cova_names <- varlist
    
  } else{
    
    cova_names <- varname[-which(varname %in% c(outcome_name, death_name, treatment_name))] #main terms
    
  }
  
  int_names <-  paste(cova_names, "*", treatment_name)  #interaction terms
  
  
  ###########################################
  ######### Create lists for output #########
  ###########################################
  
  loglik_df_list <- list()
  min_var <- list()
  prev_loglik_rec <- list()
  omit_num_rec <- list()
  
  
  ###########################################
  ############# Fit a full model ############
  ###########################################
  
  # Fit the cox model
  formula_temp_0 <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"), paste(c(cova_names, int_names, treatment_name), collapse="+")))
  fit_temp_0 <- tryCatch(coxph(formula_temp_0, data = data, weights = weight_vec, robust = robust_int),
                         error = function(c) "error",
                         warning = function(c) "warning",
                         message = function(c) "message"
  )
  
  
  
  
  # Save the output
  loglik_df <- matrix(, nrow = 1, ncol = 4)
  colnames(loglik_df) <- c("Var", "Loglik" ,"omit_num", "warning")
  loglik_df[1,1] <- "None"
  
  if(length(fit_temp_0) == 1){
    
    loglik_df[1,4] <- "warn" 
    prev_loglik <- NA
    
  } else {
    
    loglik_df[1,2] <- fit_temp_0$loglik[2]
    prev_loglik <- fit_temp_0$loglik[2] #Save the maximum likelihood
    
  }
  
  loglik_df[1,3] <- 0
  loglik_df_list[[1]] <- loglik_df
  
  
  # Save the record of log-likelihood and the number of omitted variables 
  prev_loglik_rec[[1]] <- prev_loglik
  omit_num_rec[[1]] <- 0
  
  
  ##################################################################################
  ############# Start the analysis with a feature removed one at a time ############
  ##################################################################################
  
  all_cova <-  c(cova_names, int_names)
  cand_cova_vec <- all_cova
  
  j <- 1 # Indicates the current loop number
  na_ind <- 1 # If greater than 1, remove more than 1 variable in the current cox model
  
  while(length(cand_cova_vec)>1 | length(cand_cova_vec)==1){
    
    
    print(j)
    
    
    if(na_ind > 1){
      
      if(length(cand_cova_vec)>1){
        
        
        tempmat <- combn(x = cand_cova_vec, m = na_ind)
        
        tempflag <- rep(0,ncol(tempmat))
        
        for(k in 1:ncol(tempmat)){
          
          tempvec <- tempmat[,k]
          tempvec_int <- tempvec[which(substrRight(tempvec,5) == treatment_name)]
          tempvec_main <- tempvec[which(substrRight(tempvec,5) != treatment_name)]
          
          if(any(tempvec_main %in% substr(tempvec_int,1,nchar(tempvec_int)-8))){
            
            tempflag[k] <- 1
            
          }
          
          #  if(((substrRight(tempmat[1,k],5) == "trt.3") & (substrRight(tempmat[2,k],5) != "trt.3") & 
          #    (tempmat[2,k] == substr(tempmat[1,k],1,nchar(tempmat[1,k])-8))) |
          #     ((substrRight(tempmat[2,k],5) == "trt.3") & (substrRight(tempmat[1,k],5) != "trt.3") & 
          #      (tempmat[1,k] == substr(tempmat[2,k],1,nchar(tempmat[2,k])-8)))) {
          #    
          #    tempflag[k] <- 1
          #    
          #  }
        }
        
        flagind <- which(tempflag == 1)
        tempmat2 <- tempmat[,-flagind]
        
        
        loglik_df <- matrix(,nrow = dim(tempmat2)[2], ncol = 4)
        colnames(loglik_df) <- c("Var", "Loglik" ,"omit_num", "warning")
        
        loglik_temp_vec <- c()
        
        for(i in 1:dim(tempmat2)[2]){
          
          remove_cova <- tempmat2[,i]
          
          temp_remove_cova <- remove_cova[which(substrRight(remove_cova,5) != treatment_name)]
          temp_remove_cova2 <- temp_remove_cova[which(paste(temp_remove_cova, paste("*", treatment_name)) %in% cand_cova_vec)]
          
          if(length(temp_remove_cova2) > 0){
            
            remove_cova <- c(remove_cova, paste(temp_remove_cova2, paste("*", treatment_name)))
          }
          
          cova_names_temp <- cand_cova_vec[-which(cand_cova_vec %in% remove_cova)]
          
          formula_temp <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"),
                                           paste(c(cova_names_temp, treatment_name), collapse="+")))
          
          
          fit_temp <- tryCatch(coxph(formula_temp, data = data, weights = weight_vec, robust = robust_int),
                               error = function(c) "error",
                               warning = function(c) "warning",
                               message = function(c) "message"
          )
          
          
          if(length(fit_temp) == 1){
            
            logLik_temp <- fit_temp
            loglik_df[i,4] <- "warn" 
          } else{
            
            logLik_temp <- fit_temp$loglik[2]
            
          }
          
          
          loglik_df[i,1] <- paste(remove_cova, collapse = " , ")
          loglik_df[i,2] <- logLik_temp
          loglik_df[i,3] <- length(all_cova) - length(cand_cova_vec) + length(remove_cova)
          
          loglik_temp_vec <- c(loglik_temp_vec, logLik_temp)
          
          print(i)
          
        }
        
        loglik_df_list[[j+1]] <- loglik_df
        
        loglik_temp_vec <- as.numeric(loglik_temp_vec[is.na(loglik_df[,4])])
        loglik_df <- loglik_df[is.na(loglik_df[,4]),]
        
        if(is.vector(loglik_df)){
          loglik_df <- t(as.matrix(loglik_df))
          
        }
        
        if(length(loglik_temp_vec) == 0){
          
          min_var_temp <- NA
          min_var_temp2 <- NA
          
          prev_loglik <- -Inf
          prev_loglik_rec[[j+1]] <- -Inf
          omit_num_rec[[j+1]] <- NA
          
        } else {
          
          
          min_var_temp <- loglik_df[ which(loglik_temp_vec == max(loglik_temp_vec)), "Var"]
          min_var_temp2 <- unlist(strsplit(min_var_temp, " , "))
          
          prev_loglik <- max(loglik_temp_vec)
          prev_loglik_rec[[j+1]] <- max(loglik_temp_vec)
          omit_num_rec[[j+1]] <- as.numeric(loglik_df[which(loglik_df[,"Loglik"] == max(loglik_temp_vec)), "omit_num"])
          
        }
        
        
        #  prev_loglik <- max(loglik_temp_vec)
        # prev_loglik_rec[[j+1]] <- max(loglik_temp_vec)
        #omit_num_rec[[j+1]] <- loglik_df[which(loglik_df[,"Loglik"] == max(loglik_temp_vec)), "omit_num"]
        
        min_var[[j]] <- min_var_temp2
        min_ind <- which(cand_cova_vec %in% min_var[[j]])
        
        if(length(min_ind) > 0 ){
          
          cand_cova_vec <- cand_cova_vec[-min_ind]
          na_ind <- 1
        } 
        
        if(length(min_ind) == 0){
          
          na_ind <- na_ind + 1
          
        }
        # cand_cova_vec <- cand_cova_vec[-min_ind]
        
        j <- j+1
        
        
      }else{
        
        ###    
        
        
        
        loglik_df <- matrix(,nrow = 1, ncol = 4)
        colnames(loglik_df) <- c("Var", "Loglik" ,"omit_num", "warning")
        
        loglik_temp_vec <- c()
        
        
        
        cova_names_temp <- c()
        
        formula_temp <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"),
                                         paste(c(cova_names_temp, treatment_name), collapse="+")))
        
        
        fit_temp <- tryCatch(coxph(formula_temp, data = data, weights = weight_vec, robust = robust_int),
                             error = function(c) "error",
                             warning = function(c) "warning",
                             message = function(c) "message"
        )
        
        
        if(length(fit_temp) == 1){
          
          logLik_temp <- fit_temp
          loglik_df[1,4] <- "warn" 
        } else{
          
          logLik_temp <- fit_temp$loglik[2]
          
        }
        
        
        loglik_df[1,1] <- paste("everything left", collapse = " , ")
        loglik_df[1,2] <- logLik_temp
        loglik_df[1,3] <- length(all_cova) - length(cand_cova_vec) + length(remove_cova)
        
        loglik_temp_vec <- c(loglik_temp_vec, logLik_temp)
        
        print(i)
        
        
        
        loglik_df_list[[j+1]] <- loglik_df
        
        loglik_temp_vec <- as.numeric(loglik_temp_vec[is.na(loglik_df[,4])])
        loglik_df <- loglik_df[is.na(loglik_df[,4]),]
        
        if(is.vector(loglik_df)){
          loglik_df <- t(as.matrix(loglik_df))
          
        }
        
        if(length(loglik_temp_vec) == 0){
          
          min_var_temp <- NA
          min_var_temp2 <- NA
          
          prev_loglik <- -Inf
          prev_loglik_rec[[j+1]] <- -Inf
          omit_num_rec[[j+1]] <- NA
          
        } else {
          
          min_var_temp <- loglik_df[which(prev_loglik - loglik_temp_vec == min(prev_loglik - loglik_temp_vec)), "Var"]
          min_var_temp2 <- unlist(strsplit(min_var_temp, " , "))
          
          prev_loglik <- max(loglik_temp_vec)
          prev_loglik_rec[[j+1]] <- max(loglik_temp_vec)
          omit_num_rec[[j+1]] <- loglik_df[which(loglik_df[,"Loglik"] == max(loglik_temp_vec)), "omit_num"]
          
        }
        
        
        #  prev_loglik <- max(loglik_temp_vec)
        # prev_loglik_rec[[j+1]] <- max(loglik_temp_vec)
        #omit_num_rec[[j+1]] <- loglik_df[which(loglik_df[,"Loglik"] == max(loglik_temp_vec)), "omit_num"]
        
        min_var[[j]] <- min_var_temp2
        min_ind <- which(cand_cova_vec %in% min_var[[j]])
        
        if(length(min_ind) > 0 ){
          
          cand_cova_vec <- cand_cova_vec[-min_ind]
          na_ind <- 1
        } 
        
        if(length(min_ind) == 0){
          
          na_ind <- na_ind + 1
          
        }
        # cand_cova_vec <- cand_cova_vec[-min_ind]
        
        j <- j+1
        
      }
      
      
    } else {
      
      
      print(j) 
      
      loglik_df <- matrix(,nrow = length(cand_cova_vec), ncol = 4)
      colnames(loglik_df) <- c("Var", "Loglik" ,"omit_num", "warning")
      
      
      loglik_temp_vec <- c()
      
      for(i in 1:length(cand_cova_vec)){
        
        # Obtain the covariate to be removed in the current step. If the selected covariate is main feature, remove the corresponding interaction as well.
        remove_cova <- cand_cova_vec[i]
        if(substr(remove_cova, nchar(remove_cova)-5+1, nchar(remove_cova)) != paste("*", treatment_name) & (paste(remove_cova,"*", treatment_name) %in% cand_cova_vec)){
          
          remove_cova <- c(remove_cova, paste(remove_cova, paste("*", treatment_name)))
          
        }
        
        
        # Fit a cox model with one of the covariates removed
        cova_names_temp <- cand_cova_vec[-which(cand_cova_vec %in% remove_cova)]
        
        formula_temp <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"),
                                         paste(c(cova_names_temp, treatment_name), collapse="+")))
        

        fit_temp <- tryCatch(coxph(formula_temp, data = data, weights = weight_vec, robust = robust_int),
                             error = function(c) "error",
                             warning = function(c) "warning",
                             message = function(c) "message"
        )
        
        # Store the result to the dataframe of the current loop
        if(length(fit_temp) == 1 ){
          
          logLik_temp <- fit_temp
          loglik_df[i,4] <- "warn" 
          
        } else if((length(fit_temp) > 1 & all(is.na(coef(fit_temp))))){
          
          logLik_temp <- "warning"
          loglik_df[i,4] <- "warn" 
          
        } else {
          
          logLik_temp <- fit_temp$loglik[2]
          
        }
        
        
        loglik_df[i,1] <- paste(remove_cova, collapse = " , ")
        loglik_df[i,2] <- logLik_temp
        loglik_df[i,3] <- length(all_cova) - length(cand_cova_vec) + length(remove_cova)
        
        loglik_temp_vec <- c(loglik_temp_vec, logLik_temp)
        print(i)
        
      }
      
      loglik_df_list[[j+1]] <- loglik_df
      
      loglik_temp_vec <- as.numeric(loglik_temp_vec[is.na(loglik_df[,4])])
      loglik_df <- loglik_df[is.na(loglik_df[,4]),]
      
      if(is.vector(loglik_df)){
        loglik_df <- t(as.matrix(loglik_df))
      }
      
      if(length(loglik_temp_vec) == 0){
        
        min_var_temp <- NA
        min_var_temp2 <- NA
        
        prev_loglik <- -Inf
        prev_loglik_rec[[j+1]] <- -Inf
        omit_num_rec[[j+1]] <- NA
        
      } else {
        
        if(is.na(prev_loglik)){
          
          min_var_temp <- loglik_df[which(loglik_temp_vec == max( loglik_temp_vec)), "Var"]
          
        } else{
          
          min_var_temp <- loglik_df[which(prev_loglik - loglik_temp_vec == min(prev_loglik - loglik_temp_vec)), "Var"]
          
        }
        
        min_var_temp2 <- unique(unlist(strsplit(min_var_temp, " , ")))
        
        prev_loglik <- max(loglik_temp_vec)
        prev_loglik_rec[[j+1]] <- max(loglik_temp_vec)
     
        if(length(min_var_temp2) > 1){
          
          omit_num_rec[[j+1]] <- omit_num_rec[[j]] + length(min_var_temp2)
          
        }else if(((length(min_var_temp2) == 1) & (!is.na(min_var_temp2)))){
          
          omit_num_rec[[j+1]] <- omit_num_rec[[j]] + length(min_var_temp2)
          
        }else{
          
          omit_num_rec[[j+1]] <- omit_num_rec[[j]] 
        }
        
      }
      

      min_var[[j]] <- min_var_temp2
      min_ind <- which(cand_cova_vec %in% min_var[[j]])
      
      if(length(min_ind) > 0 ){
        
        cand_cova_vec <- cand_cova_vec[-min_ind]
        na_ind <- 1
      } 
      
      if(length(min_ind) == 0){
        
        na_ind <- na_ind + 1
        
      }
      
      j <- j+1
      
    }
    
  }
  

  output <- list()
  output$loglik_df_list <-  loglik_df_list
  output$min_var <- min_var
  output$prev_loglik_rec <- prev_loglik_rec
  output$omit_num_rec <- omit_num_rec
  
  return(output)
  
  
}

