timing = function(){
  set.seed(123)

  Y = matrix(rnorm(600), ncol = 3)
  Sigma = matrix(rnorm(9), ncol = 3)
  
  res = rep(0, 200)
  
  system.time({
    for(r in 1:10000){
      for(i in 1:200){
        res[i] = Y[i,,drop = FALSE] %*% Sigma %*% t(Y[i,,drop=FALSE])
      }
      res = exp(-0.5 * res)
    }
  })
  
  
  
  system.time({
    for(r in 1:10000){
      res = apply(Y, 1, function(row){t(row) %*% Sigma %*% row})
    }
  })
}