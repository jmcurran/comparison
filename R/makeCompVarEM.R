makeCompVarEM = function(x, ...){
  UseMethod("makeCompVarEM")
}

#' @importFrom CVglasso CVglasso
makeCompVarEM.default = function(x, f, num.steps = 10, ...){
  if(!data.frame(x) || !is.factor(f)){
    stop("x must be a data.frame and f must be a factor")
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
  r = lapply(split.data, nrow)
  
  ## define the within group sum of squares function
  calcSW = function(item, mx){
    d = as.matrix(sweep(item, 2, mx))
    return(t(d) %*% d)
  }
  
  calcSB = function(item, mx){
    d = mx - mu
    return(nrow(item) * outer(d, d))
  }
  
  theta = function(mx, ni, Cinv, Uinv){
    as.vector(solve(Cinv + ni * Uinv ) %*% (Cinv %*% mu +  ni * Uinv %*% mx))
  }
  
  for(k in 1:num.steps){
    Uhat = Reduce('+', mapply(calcSW, split.data, theta0, SIMPLIFY = FALSE)) / (n - m)
    Chat = Reduce('+', mapply(calcSB, split.data, theta0, SIMPLIFY = FALSE)) / (m - 1)
    gl = CVglasso(S = Chat, path = TRUE)
    C = gl$Sigma
    Cinv = gl$Omega
    Uinv = solve(Uhat)
    theta0 = mapply(theta, theta0, r, MoreArgs = list(Cinv = Cinv, Uinv = Uinv), SIMPLIFY = FALSE)
  }
  
  obj = list(v.within = Uhat, 
             v.between = C)
             # n.observations = n.observations, 
             # n.items = n.items, 
             # item.n = item.n, 
             # item.means = item.means, 
             # n.vars = n.vars, 
             # overall.means = overall.means, 
             # multivariate = multivariate.flag, 
             # balanced = balanced.flag, 
             # s.within = s.w, 
             # s.between = s.b, 
             # warn.type = warn.type)
  class(obj) = "compvar"
  return(obj)
  
}