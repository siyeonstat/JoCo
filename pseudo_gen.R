psuedo_gen <- function(data, treatment, outcome, id, idcreate = T,
                        stage, sample_time = 10, fit, truncate = F, tau, seed = NULL){
  
  set.seed(seed)
  
  dat <- data.frame(data, surv = 0)
  
  stage_end <- paste(".", stage, sep ="")
  surv_new_name <- paste("survnew", stage_end, sep = "")  
  colnames(dat)[length(colnames(dat))] <- surv_new_name
  death_name <- paste("death", stage_end, sep="")
  
  
  #Sample 10 times, repeat the same row 10 times.
  datn <- dat
  datn <- datn[rep(seq_len(nrow(datn)), each = sample_time),]
  
  treatment_name <- paste(treatment, stage_end, sep = "")  
  old_treatment_name <- paste("old", treatment, stage_end, sep = "")  
  surv_new_name <- paste("survnew", stage_end, sep = "")  
  outcome_name <- paste(outcome, stage_end, sep = "")  
  
  outcome_name1 <- paste(outcome, ".1", sep = "")  
  outcome_name2 <- paste(outcome, ".2", sep = "")  
  
  if(truncate == T){
    
    for(i in 1:nrow(dat)){
      print(i)
      if((!(is.na(dat[, treatment_name][i]))) & (dat[, treatment_name][i] != dat[, old_treatment_name][i])){
        
        surv <-survfit(fit, dat[i,])
        
        temp_sampletime <-  sample(summary(surv)$time, size = sample_time, prob = summary(surv)$n.event/sum(summary(surv)$n.event))
        
        if(stage == 3){
          
          temp_totaltime <- temp_sampletime + dat[i, outcome_name2] + dat[i, outcome_name1]
          temp_sampletime[temp_totaltime>tau]  <- tau -  dat[i, outcome_name2] - dat[i, outcome_name1] 
          
        } else if(stage == 2){
          
          temp_totaltime <- temp_sampletime + dat[i, outcome_name1]
          temp_sampletime[temp_totaltime>tau]  <- tau - dat[i, outcome_name1] 
          
        } else if(stage == 1){
          
          temp_totaltime <- temp_sampletime 
          temp_sampletime[temp_totaltime>tau]  <- tau 
          
        }
        
        datn[, surv_new_name][(sample_time*i-(sample_time-1)):(sample_time*i)] <- temp_sampletime
        
      } else {
        
        
        datn[, surv_new_name][(sample_time*i-(sample_time-1)):(sample_time*i)] <- dat[, outcome_name][i]
        
      }
      
      
    }
  } else if(truncate == F){
    
    for(i in 1:nrow(dat)){
      print(i)
      if( (!(is.na(dat[, old_treatment_name][i]))) & (!(is.na(dat[, treatment_name][i]))) & (dat[, treatment_name][i] != dat[, old_treatment_name][i])){
        #browser()
        surv <-survfit(fit, dat[i,])
        #datn[, surv_new_name][(sample_time*i-(sample_time-1)):(sample_time*i)] <- sample(summary(surv)$time, size = sample_time, prob = summary(surv)$n.event/sum(summary(surv)$n.event))
        
        
        survvec <- -diff(surv$surv, lag = 1)
        #survvec <- diff(rev(surv$surv), lag = 1)
        survprobdf <- data.frame(time = surv$time[-1], survvec)
        
        survprobdf<- survprobdf[survprobdf$survvec>0,]
        
        if(nrow(survprobdf) == 0){
          
          survprobdf[1,1] <- surv$time[1]
          survprobdf[1,2] <- 1
          
          colnames(survprobdf) <- c("time", "survvec")
          
        }
        
        datn[, surv_new_name][(sample_time*i-(sample_time-1)):(sample_time*i)] <- sample(survprobdf$time, sample_time, prob= survprobdf$survvec, replace = T)
        
      } else {
        
        datn[, surv_new_name][(sample_time*i-(sample_time-1)):(sample_time*i)] <- dat[, outcome_name][i]
        
      }
      
      
    }
    
  }
  
  
  #datn[, id]
  
  if(idcreate == T){
    for(i in 1:(nrow(datn)/sample_time)){
      print(i)
      datn[, id][sample_time*i-(sample_time-1)] <- paste(datn[, id][sample_time*i-(sample_time-1)], ".1", sep="")
      datn[, id][sample_time*i-(sample_time-2)] <- paste(datn[, id][sample_time*i-(sample_time-2)], ".2", sep="")
      datn[, id][sample_time*i-(sample_time-3)] <- paste(datn[, id][sample_time*i-(sample_time-3)], ".3", sep="")
      datn[, id][sample_time*i-(sample_time-4)] <- paste(datn[, id][sample_time*i-(sample_time-4)], ".4", sep="")
      datn[, id][sample_time*i-(sample_time-5)] <- paste(datn[, id][sample_time*i-(sample_time-5)], ".5", sep="")
      datn[, id][sample_time*i-(sample_time-6)] <- paste(datn[, id][sample_time*i-(sample_time-6)], ".6", sep="")
      datn[, id][sample_time*i-(sample_time-7)] <- paste(datn[, id][sample_time*i-(sample_time-7)], ".7", sep="")
      datn[, id][sample_time*i-(sample_time-8)] <- paste(datn[, id][sample_time*i-(sample_time-8)], ".8", sep="")
      datn[, id][sample_time*i-(sample_time-9)] <- paste(datn[, id][sample_time*i-(sample_time-9)], ".9", sep="")
      datn[, id][sample_time*i] <- paste(datn[, id][sample_time*i], ".10", sep="")
      
    }
  }
  #change id for new data
  
  
  return(datn)
  
}

