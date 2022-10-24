
opt_treatment <- function(data, treatment, stage, fit){
  
  stage_end <- paste(".", stage, sep ="")
  treatment_name <- paste(treatment, stage_end, sep = "")  
  
  dat0 <- data
  dat0[,treatment_name] <- factor(0, levels = levels(data[,treatment_name]))
  
  dat1 <- data
  dat1[,treatment_name] <- factor(1, levels = levels(data[,treatment_name]))
  
  datn <- data.frame(data, old = data[,treatment_name])
  old_treatment_name <- paste("old", treatment, stage_end, sep = "")  
  colnames(datn)[length(colnames(datn))] <- old_treatment_name
  
  datn[, treatment_name] <- ifelse(predict(fit, dat0, type = "risk") > predict(fit, dat1, type = "risk"),
                                   1, 0)
  datn[, treatment_name] <- as.factor(datn[, treatment_name])
  
  
  return(datn)
  
}