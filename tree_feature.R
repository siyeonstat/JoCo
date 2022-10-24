tree_feature <- function(x, y, censor, number = 20, seed = NULL){
  
  set.seed(seed)
  temp_rlt <- RLT(x, y, model = "survival", censor = censor)
  imp <-temp_rlt$VarImp[,order(temp_rlt$VarImp[1,], decreasing = T)]
  impcan <- names(imp)[1:number]
  
  return(impcan)
}