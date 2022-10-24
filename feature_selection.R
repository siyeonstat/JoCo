
feature_selection <- function(y = unlist(out$prev_loglik_rec), min_var = out$min_var, omit_num = out$omit_num_rec, treatment, outcome,
                              death = "death", data, stage , weight = 1, robust_int = F, isvar = F, varlist, slopepoint){
  
  weight_vec = rep(weight, nrow(data))
  ###########################################
  ##### Variables used for the analysis #####
  ###########################################
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  # Extract variables for the current analysis(loop) by identifying the stage number
  stage_end <- paste(".", stage, sep ="")
  varname <- names(data)[which(substrRight(colnames(data), 2) == stage_end )]
  
  # Put the stage number to the outcome, death and treatment
  # Get the list of main terms and interaction terms for the analysis
  outcome_name <- paste(outcome, stage_end, sep = "")  
  death_name <- paste(death, stage_end, sep = "")  
  treatment_name <- paste(treatment, stage_end, sep = "")  
  
  
  if(isvar == T){
    
    cova_names <- varlist
    
  } else{
    
    outcome_name2 <- paste(outcome2, stage_end, sep = "")  
    
    cova_names <- varname[-which(varname %in% c(outcome_name, death_name, treatment_name, outcome_name2))] #main terms
    
  }
  
  int_names <-  paste(cova_names, "*", treatment_name)  #interaction terms
  
  all_cova <-  c(cova_names, int_names)
  
  omit_num_na <- which(is.na(omit_num))
  
  if(length(omit_num_na) > 0){
    
    y <- y[-omit_num_na]
    omit_num <- omit_num[-omit_num_na]
    
  }
  
  yna <- ifelse(is.na(y[1]), 1, 0)
  yinf <- length(y[y==-Inf & !is.na(y)])
  y_old <- y[!is.na(y)]
  y <- y[!is.na(y) & y>-Inf]
  
  
  # x <- 1:length(y)
  if(yna == 1){
    
    x <- omit_num[-1]
    
  } else {
    
    x <- omit_num
  }
  
  ######
  # x<- 1:length(y)
  
  z <- c() #z is index
  
  
  z <- c(1, z)
  
  slope1 <- vector(mode = 'numeric', length = length(y))
  slope1[1] <- 0
  
  
  l_temp <- vector(mode = 'numeric', length = length(y))
  l_temp[1] <- y[1]
  
  j = 1
  
  while(j < length(y)){
    #  while(j < 46){  
    for (i in (j+1):length(slope1)){
      
      slope1[i] <- (y[i]-y[j])/(x[i]-x[j])
      l_temp[i] <- y[i]
      
    }
    
    j_prev <- j
    
    #order1 <- order(abs(slope1))
    #j <- order1[!(order1 %in% (1:j_prev))][1]
    
    j <- which.min(abs(slope1[-(1:j_prev)])) + j_prev
      

    z <- c(z, j)
    
    if((j-j_prev) > 1){
      
      j_ind <- (j_prev+1):(j-1)
      
    } else {
      
      j_ind <- j
    }
    
    slope1[j_ind] <- slope1[j] 
    l_temp[j_ind] <- y[j]
    
  }
  
  l_temp_x <- c()
  l_temp_y <- c()
  
  for(i in 1:length(l_temp)){
    
    if(!is.na(l_temp[i+1]) & (l_temp[i+1]-l_temp[i] != 0)){
      
      l_temp_x <- c(l_temp_x, x[i])
      l_temp_y <- c(l_temp_y,  l_temp[i])
      
    }  else if(i == length(l_temp)){
      
      l_temp_x <- c(l_temp_x, x[i])
      l_temp_y <- c(l_temp_y, l_temp[i] )
      
      
    } else if(i == 1){
      
      l_temp_x <- c(l_temp_x, x[i])
      l_temp_y <- c(l_temp_y, l_temp[i] )
    }
  }
  
  
  dat_plot <- data.frame(x = x, y = y, slope= slope1, l=l_temp)  
  dat_plot2 <- data.frame(x = l_temp_x, y=l_temp_y)
  
  # fit_plot <- ggplot(dat_plot, aes(x, y)) + xlab("Number of covariates removed from the cox model") + ylab("Log-likelihood") +
  #  geom_point() +
  #  geom_line(aes(y = slope), color = "red") 
  
 fit_plot <- ggplot(dat_plot, aes(x, y)) + xlab("Number of covariates removed from the cox model") + ylab("Log-likelihood") +
      geom_point() +
      geom_line(data=dat_plot2, color = "red" )  
    
 
  
  
  z1 <- z
  
  # slope1 <- slope1[-1]
  #for(i in 1:(length(slope1)-1)){
  
  #  if(slope1[i] > slope1[i+1]){
  #    
  #    z1 <- c(z1, i)
  
  #  }
  
  #  }
  diff1 <- diff(slope1, 1)
  
  if(length(diff1)>3){
    
    diff1temp <- diff1
    diff1temp[diff1temp > slopepoint] <- 0 
    
    zeronum <- sum(diff1temp == 0)
    
    
    indexslope <- order(diff1temp)[1:(length(diff1temp)-zeronum)]
    # indexslope <- which.min(diff1temp)
    
  } else {
    
    indexslope <- 2
  }
  
  z1 <- indexslope
  
  
  #orderslope <- dat_plot[z,][order(dat_plot[z,"slope"]),]
  
  
  #addtoz1 <- as.numeric(rownames(orderslope[1:floor(nrow(orderslope)/2),]))
  
  #z1 <- c(z1, addtoz1)
  #z1 <- sort(unique(z1), decreasing = T)
  
  
  
  #  slope2 <- c()
  #  for(i in 1:(length(slope1)-1)){
  
  #   slope2 <- c(slope2, slope1[i+1]-slope1[i]) 
  
  #  }
  
  
  #z and z1 are index
  
  
  
  
  s <-1
  ind <- 0
  
  
  # z2 <- rev( z[!(z %in% z1)])
  s2 <- 1
  
  x1 <- x[z1]
  
  while(ind == 0){
    
    if(!is.null(z1) & length(z1) > 0 & !(s > length(z1))){
      
      rem_final <- c()
      if(yna == 1){
        
        if(z1[s] == length(unlist(min_var))){
          
          rem_final <- all_cova
          
        } else{
          
          for(j in 1:(z1[s])){
            
            rem_final <- c(rem_final, min_var[[j]])
            
          }
        }
      } else{
        
        if(z1[s] > 1){
          
          for(j in 1:(z1[s]-1)){
            
            rem_final <- c(rem_final, min_var[[j+yinf]])
          }
        }
        
        
        
        
      }
      
      cova_final <- all_cova[-which(all_cova %in% rem_final)]
        
      
      
      formula_final <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"),
                                        paste(c(cova_final, treatment_name), collapse="+")))
      
      fit_final <- tryCatch(coxph(formula_final, data = data, weights = weight_vec, robust = robust_int),
                            error = function(c) "error",
                            warning = function(c) "warning",
                            message = function(c) "message"
      )
      
      if(length(fit_final)>1){
        
        summary_fit <- summary(fit_final)$coefficient
        
        p_int <- summary_fit[substr(rownames(summary_fit),1,7) %in%  paste(treatment_name, "1:", sep="") 
                             | substrRight(rownames(summary_fit),7) %in% paste(":", treatment_name,1,  sep=""), "Pr(>|z|)"]
        
        
        
        
      } else{
        
        summary_fit <- NA
        
        p_int <- NA
        
      }
      
      
      # p_val_int <- min(summary_fit[substr(rownames(summary_fit),1,7) %in%  paste(treatment_name, 1, ":", sep="") 
      #            | substrRight(rownames(summary_fit),7) %in% paste(":", treatment_name, 1,  sep=""), "Pr(>|z|)"])
      
      
      # main_ind_max <- max(which(!(substrRight(cova_final, 5) %in% treatment_name)))
      
      #p_val_int <- summary(fit_final)$coefficient[-(1:(main_ind_max+1)) , "Pr(>|z|)"]
      
      
      fit_plot <- fit_plot + geom_point(data = dat_plot[z1[s],], pch=21, fill=NA, size=4, colour="blue", stroke=1) +
        annotate("text", x = x[z1[s]]+1, y = y[z1[s]]+5, label = paste("M", z1[s]-1, sep = ""))
      
      if((length(p_int) >0) & any(!is.na(p_int))){
        
        ind <- 1
        
      } else {
        
        s <- s+1
        
      }
      
    } else{
      
      rem_final <- c()
      if(yna == 1){
        
        if(z2[s2] == length(unlist(min_var))){
          
          rem_final <- all_cova
          
        } else{
          
          for(j in 1:(z2[s2])){
            
            rem_final <- c(rem_final, min_var[[j]])
            
          }
        }
      } else{
        
        if(z2[s2] > 1){
          
          for(j in 1:(z2[s2]-1)){
            
            rem_final <- c(rem_final, min_var[[j+yinf]])
          }
        }
        
        
        
        
      }
      
    
        
        if(is.null(rem_final)){
          
          cova_final <- all_cova
          
        } else{
          
          cova_final <- all_cova[-which(all_cova %in% rem_final)]
          
        }
        
       
      
      
      formula_final <- as.formula(paste(paste("Surv(", outcome_name, ",", death_name,") ~"),
                                        paste(c(cova_final, treatment_name), collapse="+")))
      
      fit_final <- tryCatch(coxph(formula_final, data = data, weights = weight_vec, robust = robust_int),
                            error = function(c) "error",
                            warning = function(c) "warning",
                            message = function(c) "message"
      )
      #fit_final <- coxph(formula_final, data = data, weights = weight_vec, robust = robust_int)
      
      if(length(fit_final)>1){
        
        summary_fit <- summary(fit_final)$coefficient
        
        p_int <-summary_fit[substr(rownames(summary_fit),1,7) %in%  paste(treatment_name, "1:", sep="") 
                            | substrRight(rownames(summary_fit),7) %in% paste(":", treatment_name,1,  sep=""), "Pr(>|z|)"]
        
        
      } else{
        
        summary_fit <- NA
        
        p_int <- NA
        
      }
      
      fit_plot <- fit_plot + geom_point(data = dat_plot[z2[s2],], pch=21, fill=NA, size=4, colour="blue", stroke=1) +
        annotate("text", x = x[z2[s2]]+1, y = y[z1[s]]+5, label = paste("M", z2[s2]-1, sep = ""))
      
      if((length(p_int) >0) & any(!is.na(p_int))){
        
        ind <- 1
        
      } else {
        
        s2 <- s2+1
        
      }
      
      if(length(z2) < s2){
        
        stop("ERROR!!!!!NO MODEL!!!!!!!!!!")
      }
      
    }
    
  }
  
  
  output <- list()
  output$fit <- fit_final
  output$covariates <- cova_final
  output$plot <- fit_plot
  output$cova_final <- cova_final
  
  return(output)
  
}


