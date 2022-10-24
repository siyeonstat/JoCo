

censor_func <- function(dat = dat1_0, outcome_name = "total_time", censor_name = "censor", 
                        variable_name = cova_cen){
  
  formula_cen <- as.formula(paste(paste("Surv(", outcome_name, ",", censor_name,") ~"),
                                  paste(c(variable_name), collapse="+")))
  
  fit_cen <- coxph(formula_cen, data = dat)
  
  
  
  output <-  ifelse(dat[,censor_name] == "1", predict(fit_cen, data = dat, type = "survival"), 
                    ifelse(dat[,censor_name] == "0", 1-predict(fit_cen, data = dat, type = "survival"),NA))
  
  return(output)  
  
  
}
