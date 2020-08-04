makeCompVarEM = function(x, ...){
  UseMethod("makeCompVarEM")
}

makeCompVarEM.default = function(x, f, num.steps = 10){
  if(class(x) == "data.frame"){
    x = as.matrix(x)
  }
  
  if(class(x) != "matrix" || class(f) != factor){
    stop("x must be a matrix and f must be a factor")
  }
  
  if(length(f) != nrow(x)){
    stop("f must have as many entries as x has rows")
  }
  
  split.data = split(x, f)
  
  ## mu - grand mean
  
  mu = colMeans(x)
  
  theta0 = do.call("rbind", lapply(split.data, colMeans))
  
  for(k in 1:num.steps){
    calcSW = function(item, mx){
      d = as.matrix(sweep(item, 2, mx))
      return(nrow(item) * t(d) %*% d)
    }
    
    s.w = Reduce('+', mapply(split.data, calcSW))
  }
  
}