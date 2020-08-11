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
  
  n = nrow(x)
  m = length(split.data)
  
  ## mu - grand mean
  
  mu = colMeans(x)
  
  ## theta - group means
  theta0 = lapply(split.data, colMeans)
  
  ## define the within group sum of squares function
  calcSW = function(item, mx){
    d = as.matrix(sweep(item, 2, mx))
    return(nrow(item) * t(d) %*% d)
  }
  
  calcSB = function()
  
  
  
  for(k in 1:num.steps){
    U = Reduce('+', mapply(calcSW, split.data, theta0)) / 1 # Note 1 is a placeholder - clearly this is wrong
    
    
  }
  
}