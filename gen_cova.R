





gen_cova <- function(outcome = "time_intR", outcome2="time_intR2", death = "death", treatment = "trt", stage = 1, data = dat0_6, only = T){
  
  stage_end = paste(".", stage, sep = "")
  
  outcome_name <- paste(outcome, stage_end, sep = "")  
  death_name <- paste(death, stage_end, sep = "")  
  treatment_name <- paste(treatment, stage_end, sep = "")  
  outcome_name2 <- paste(outcome2, stage_end, sep = "")  
  
  substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
  }
  
  varname <- names(data)[which(substrRight(colnames(data), 2) == stage_end )]
  
  if(only == T){
  
    cova_names1 <- varname[-which(varname %in% c(outcome_name, death_name, treatment_name, outcome_name2))] #main terms
    
  } else if(only ==F){
    
    cova_names1 <- varname
  }
  
  return(cova_names1)
  
}

