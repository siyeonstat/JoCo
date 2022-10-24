
prop_func_rf <- function(data, variable_name, treatment, seed = NULL){
  
  set.seed(seed)
  formula_prop1 <- as.formula(paste( treatment, "~" , paste(variable_name, collapse="+")))
  
  #data[, treatment] <- as.numeric(data[, treatment])-1
  
  prop_model1 <- randomForest(formula_prop1, data = data)
  
  # prop_model1 <- glm(formula_prop1, data = data,
  #                   family = binomial(link = "logit"))
  
  #predict(prop_model1, data = data, type = "prob")
  output <-  ifelse(data[,treatment] == "1", predict(prop_model1, data = data, type = "prob")[,2], predict(prop_model1, data = data, type = "prob")[,1])
  
  return(output)
  
  
}


