#' Calculate the likelihood ratio
#' 
#' Takes a `compitem` object which represents some control item, and a
#' `compitem` object which represents a recovered item, then uses information
#' from a `compcovar` object, which represents the information from the
#' population, to calculate a likelihood ratio (LR) as a measure of the evidence
#' given by the observations for the same/different source propositions.

#'
#' @param control a `compitem` object with the control item information.
#' @param recovered a `compitem` object with the recovered item information.
#' @param background a `compcovar` object with the population information.
#' @param method a choice of the method used to calculate the LR. Presently there
#' are three methods, `"mvn"` - multivariate normal approximation, `"kde"` - (multivariate) kernel
#' density estimates and `"lindely"` which uses the method published by Lindley (1977).
#'
#' @return an estimate of the likelihood ratio
#' 
#' @references Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.
#' @export
#'
#' @examples
#' data(glass)
#' 
#' controlMeasurements = subset(glass, item == "s1")
#' control = makeCompItem(item ~ logKO + logCaO + logFeO, 
#'                        data = controlMeasurements[1:6,])
#' recovered.1 = makeCompItem(item ~ logKO + logCaO + logFeO, 
#'                        data = controlMeasurements[7:12,])
#' recoveredMeasurements = subset(glass, item == "s2")
#' recovered.2 = makeCompItem(item ~ logKO + logCaO + logFeO,
#'                            data = recoveredMeasurements[7:12,])
#'                            
#' background = makeCompVar(item ~ logKO + logCaO + logFeO, data = glass)
#'                            
#' ## Same source comparison using a multivariate normal (MVN) approximation
#' calcLR(control, recovered.1, background)
#' 
#' ## Same source comparison using a multivariate kernel density estimate (MVK) approximation
#' calcLR(control, recovered.1, background, "kde")
#' 
#' ## Different source comparison using a multivariate normal (MVN) approximation
#' calcLR(control, recovered.2, background)
#' 
#' ## Different source comparison using a multivariate kernel density estimate (MVK) approximation
#' calcLR(control, recovered.2, background, "kde")
#' 
calcLR = function(control, recovered, background, method = c("mvn","kde","lindley")){
  if(!is(control, "compitem") || !is(recovered, "compitem")){
    stop("control and recovered must be of class compitem")
  }
  
  if(!is(background, "compvar")){
    stop("background must be of class compvar")
  }
  
  method = match.arg(method)
  
  if(method == "mvn"){ ## I might be able to do this a smarter way with switch, but I can't see how just yet
    calcLR_MVN(control, recovered, background)
  }else if(method == "kde"){
    calcLR_KDE(control, recovered, background)
  }else{
    
  }
}

calcLR_MVN = function(control, recovered, background){
  ## Calculates the likelihood ratio for a multivariate random effect model with
  ## between items modelled normal rather than kernal REQUIRES control - a
  ## compitem object calculated from the observations from the item considered to
  ## be the control item - calculated from two.level.comparison.items() from the
  ## file items_two_level.r recovered - a compitem object calculated from the
  ## observations from the item considered to be the recovered item - calculated
  ## from two.level.comparison.items() from the file items_two_level.r
  ## 
  ## background - a \code{compcovar} object calculated from the observations of the population as a whole -
  ## calculated from the two.level.components() function from the file
  ## components_two_level.r RETURNS LR - an estimate of the likelihood ratio ' '
  
  # redefine some of the items sent first the object with the population information
  U = background$v.within
  C = background$v.between
  mu = background$overall.means
  
  # then the objects with the contro and recovered item information
  mean.cont = control$item.means
  n.cont = control$n.replicates
  mean.rec = recovered$item.means
  n.rec = recovered$n.replicates
  
  # weighted mean for the control and recovered item
  y.star = ((n.cont * mean.cont) + (n.rec * mean.rec))/(n.cont + n.rec)
  
  # calculate some of the repeated components
  diff.cont.rec = mean.cont - mean.rec
  diff.cont.mu = mean.cont - mu
  diff.rec.mu = mean.rec - mu
  diff.y.star.mu = y.star - mu
  
  
  # This code independently written by Agnieszka Rzepecka of the IFR, and tidied
  # up a bit by David Lucy. The output agrees with previous efforts - thus we develop
  # some confidence it's doing the right thing I keep on expecting it to fall
  # over in the univariate case but it seems not to so we'll use it until it
  # starts to be a problem 
  
  ## Numerator calculation Cleaned up and hopefully sped up by James Curran
  k1 = (n.cont * n.rec) / (n.cont + n.rec)
  logNumerator1 = k1 * t(diff.cont.rec) %*% solve(U) %*% 
                  (diff.cont.rec) + log(det(U / k1))
  logNumerator2 = t(diff.y.star.mu) %*% 
                  solve(U / (n.cont + n.rec) + C) %*%
                  (diff.y.star.mu) +  
                  log(det(U / (n.cont + n.rec) + C))
  logNumerator = -0.5 * (logNumerator1 + logNumerator2)
  
  # Denominator calculation
  # denom1 = exp(-1/2 * t(diff.cont.mu) %*% solve(U/n.cont + C) %*% (diff.cont.mu)) * (det(U/n.cont + C))^(-1/2)
  # denom2 = exp(-1/2 * t(diff.rec.mu) %*% solve(U/n.rec + C) %*% (diff.rec.mu)) * (det(U/n.rec + C))^(-1/2)
  # denom = denom1 * denom2
  # 
  logDenominator1 = -0.5 * (t(diff.cont.mu) %*% solve(U / n.cont + C) 
                            %*% diff.cont.mu + log(det(U / n.cont + C)))
  logDenominator2 = -0.5 * (t(diff.rec.mu) %*% solve(U / n.rec + C) 
                            %*% diff.rec.mu + log(det(U / n.rec + C)))
  logDenominator = logDenominator1 + logDenominator2
  
  LR = as.numeric(exp(logNumerator - logDenominator))
  #browser()
  ## original code from 2004 which follows the Aitken & Lucy paper (p.115) closely 
  ## this is here so people trying to follow the code with
  ## the paper can do so 
  ## D1 = U / n.cont 
  ## D2 = U / n.rec 
  ## H2 = t(y.star - mu) %*% solve(U/(n.rec + n.cont) + C) %*% (y.star - mu) 
  ## H3 = t(diff.cont.rec) %*% solve(D1 + D2) %*% (diff.cont.rec) 
  ## num1 = exp(- 0.5 * (H2 + H3)) 
  ## num2 = sqrt(det(2 * pi * solve(((n.cont + n.rec) * solve(U) ) + solve(C)))) 
  ## num = num1 * num2
  ## mu.star = solve(solve(D1 + C) + solve(D2 + C)) %*% ((solve(D1 + C) %*% cont.means) +
  ## (solve(D2 + C) %*% rec.means)) 
  ## H4 = t(mu - mu.star) %*% (solve(D1 + C) + solve(D2 + C)) %*% (mu - mu.star) 
  ## H5 = t(diff.cont.rec) %*% solve(D1 + D2 + (2 * C)) %*% diff.cont.rec 
  ## den1 = 1 / sqrt(det(2 * pi * C)) 
  ## den2 = sqrt(det(2 * pi * solve(n.cont * solve(U) + solve(C)))) 
  ## den3 = sqrt(det(2 * pi * solve(n.rec * solve(U) + solve(C)))) 
  ## den4 = exp(-0.5 * (H4 + H5)) lr = num / prod(den1, den2, den3, den4)
  
  return(LR)
}

calcLR_KDE = function(control, recovered, background) {
  ## Calculates the likelihood ratio for a multivariate random effects with
  ## between items modelled as kernel densities This function could still do with
  ## being made a bit quicker I have tried using apply() type formulations where
  ## appropriate but they are slightly slower than the counted iterations - I
  ## suspect the fact that they are dealing with matrix operations has something
  ## to do with it we're using the prod(eigen(A)$values) form rather than det(A)
  ## as it seems to be a little more reliable as some of the matricies tend to
  ## singularity REQUIRES control - a compitem object calculated from the
  ## observations from the item considered to be the control item - calculated
  ## from two.level.comparison.items() from the file items_two_level.r recovered -
  ## a compitem object calculated from the observations from the item considered
  ## to be the recovered item - calculated from two.level.comparison.items() from
  ## the file items_two_level.r
  
  
  
  ## calculates the difference between two numbers intended for use as an
  ## apply() functionette x is a matrix from whose rows we wish to subtract y
  ## y is a vector
  minu = function(x, y) {
    y - x
  }
  
  
  # convert to local definitions names changed to conform to three level code
  control.mean = control$item.means
  recovered.mean = recovered$item.means
  Nc = control$n.replicates
  Nr = recovered$n.replicates
  
  n.variables = background$n.vars
  
  
  
  # try to trap any errors here
  multivariate.flag = n.variables > 1
  
  if (!exists("multivariate.flag")) {
    stop("undefined number of variables")
  }
  
  n.groups = background$n.items
  group.means = background$item.means
  
  U = background$v.within
  C = background$v.between
  
  
  # window width calculation was formerly calculated by a separate function
  h.opt = (((4 / ((2 * n.variables) + 1)) ^ (1 / (n.variables + 4))) * (n.groups ^ (-(1 / (n.variables + 4)))))
  
  # print(h.opt)
  
  # stop()
  
  # debug code assign('bit', control.mean, env=.GlobalEnv)
  
  # do some essential type redefinition Nc = as.numeric(Nc) Nr = as.numeric(Nr) control.mean = as.matrix(control.mean) recovered.mean =
  # as.matrix(recovered.mean)
  
  
  # stop('systematic halt')
  
  
  if (multivariate.flag) {
    D.control = U / Nc
    D.recovered = U / Nr
    
    if (!is.numeric(try(solve(U))
    )) {
      value = "NA"
      return(value)
    }
    # if(! is.numeric(try(solve(C)))){value = 'NA'; return(value)}
    
    # component way of getting the inverses of the D matrices
    inv.U = solve(U)
    inv.D.control = inv.U * Nc
    inv.D.recovered = inv.U * Nr
    
    control.minus.recovered = control.mean - recovered.mean
    
    A = (Nc + Nr) * inv.U
    inv.A = U / (Nc + Nr)
    
    # assign('A', A, .GlobalEnv) # useful bit of debug code
    
    # print(control.mean) assign('A', control.mean, .GlobalEnv) assign('B', inv.D.control, .GlobalEnv)
    
    
    y.star = inv.A %*% ((inv.D.control %*% (control.mean)) + (inv.D.recovered %*% (recovered.mean)))
    ##
    
    # print(prod(eigen(C)$values)) print(det(C)) stop('HHHH')
    
    
    
    #top1 = sqrt(abs(prod(eigen(C)$values)))
    top1 = sqrt(det(C))
    
    ## I have removed this because it is just wrong.
    ## If the determinant of C is negative, then it means the estimate of the covariance matrix
    ## is negative definite - i.e. it is no longer a covariance matrix - and we should not be using it
    ## for computation. The function should fall over when top1 is computed, as it will be taking a sqrt of a negative number
    # trap dodgy covariance matrices
    # if (prod(eigen(C)$values) < 0) {
    #   warning(
    #     "negative determinant - taking absolute value",
    #     call. = TRUE,
    #     immediate. = TRUE
    #   )
    # }
    
    
    top2 = n.groups * (h.opt ^ n.variables)
    inv.h.opt.squared.times.C = solve((h.opt ^ 2) * C)
    
    # top3 = 1 / sqrt(prod(abs(
    #   eigen(A + inv.h.opt.squared.times.C)$values
    # )))
    top3 = det(A + inv.h.opt.squared.times.C)^-0.5
    
    
    top4 = exp(
      -0.5 * (control.minus.recovered) %*% solve(D.control + D.recovered) %*% (control.minus.recovered)
    )
    matt1 = rep(0, n.groups)
    
    # stop()
    
    # this invocation of minu makes a little difference to speed of execution
    # y.star.minus.mean = t(apply(group.means, 1, minu, y = y.star)) 
    ## Replaced by James - this is almost certainly faster, but it is unlikely to make much difference
    y.star.minus.mean = sweep(-group.means, 2, y.star, '+')
    numerator.constant = solve(inv.A + (h.opt ^ 2) * C)
    
    
    for (ctr in 1:n.groups) {
      matt1[ctr] = y.star.minus.mean[ctr, , drop = FALSE] %*% numerator.constant %*% t(y.star.minus.mean[ctr, , drop = FALSE])
    }
    
    top5 = sum(exp(-0.5 * matt1))
    
    numerator = prod(top1, top2, top3, top4, top5)
    ##
    
    
    
    
    ##
    # bot1 = 1 / sqrt(abs(prod(
    #   eigen(inv.D.control + inv.h.opt.squared.times.C)$values
    # )))
    bot1 = det(inv.D.control + inv.h.opt.squared.times.C)^-0.5
    
    matt2 = rep(0, n.groups)
    control.konstant = solve(D.control + ((h.opt ^ 2) * C))
    
    
    control.mean.minus.group.means = sweep(-group.means, 2, control.mean, '+')
    
    for (ctr in 1:n.groups) {
      matt2[ctr] = control.mean.minus.group.means[ctr,,drop = FALSE] %*% control.konstant %*% t(control.mean.minus.group.means[ctr,,drop = FALSE])
    }
    
    
    ## the matrix parts can equally well be done as the following apply() type of formulation - unfortunately despite its cleverness it is
    ## slower than the iterative approach ma2 = t(apply(group.means, 1, minu, y=control.mean)) konstant = solve(D.control + ((h.opt ^ 2) * C))
    ## ma3 = apply(ma2, 1, mulp, y=konstant) ma3 = exp(-0.5 * ma3) ma3 is now equal to matt2
    
    bot2 = sum(exp(-0.5 * matt2))
    # bot3 = 1 / sqrt(abs(prod(
    #   eigen(inv.D.recovered + inv.h.opt.squared.times.C)$values
    # )))
    bot3 = det(inv.D.recovered + inv.h.opt.squared.times.C)^-0.5
    
    matt3 = rep(0, n.groups)
    recovered.konstant = solve(D.recovered + ((h.opt ^ 2) * C))
    
    # assign('cc', inv.A, envir=.GlobalEnv)
    
    recovered.mean.minus.group.means = sweep(-group.means, 2, recovered.mean, '+')
    
    for (ctr in 1:n.groups) {
      matt3[ctr] = recovered.mean.minus.group.means[ctr, , drop = FALSE] %*% recovered.konstant %*% t(recovered.mean.minus.group.means[ctr, , drop = FALSE])
    }
    
    bot4 = sum(exp(-0.5 * matt3))
    
    denominator = prod(bot1, bot2, bot3, bot4)
    
    LR = numerator / denominator
    return(LR)
    ##
  } else{
    ## Univariate
    h = h.opt
    k = nrow(group.means)
    a.sq = (1 / Nc) + (1 / Nr)
    
    
    w = ((Nc * control.mean) + (Nr * recovered.mean)) / (Nc + Nr)
    
    K.num = k * sqrt(Nc + Nr) * sqrt(U + (Nc * C * (h ^ 2))) * sqrt(U + (Nr * C * (h ^ 2)))
    K.den = sqrt(a.sq) * sqrt(U) * sqrt(Nc * Nr) * sqrt(U + ((Nc + Nr) * C * (h ^ 2)))
    K = K.num / K.den
    
    
    bit1 = ((control.mean - recovered.mean) ^ 2) / (2 * a.sq * U)
    bit2 = 2 * (U + ((Nc + Nr) * C * (h ^ 2)))
    bit3 = 2 * (U + (Nc * C * (h ^ 2)))
    bit4 = 2 * (U + (Nr * C * (h ^ 2)))
    
    num1 = 0
    den1 = 0
    den2 = 0
    
    for (ctr in 1:k) {
      tmp = ((Nc + Nr) * ((w - group.means[ctr]) ^ 2)) / bit2
      num1 = num1 + exp(-tmp)
      
      tmp = (Nc * ((control.mean - group.means[ctr]) ^ 2)) / bit3
      den1 = den1 + exp(-tmp)
      
      tmp = (Nr * ((recovered.mean - group.means[ctr]) ^ 2)) / bit4
      den2 = den2 + exp(-tmp)
    }
    
    numerator = K * exp(-bit1) * num1
    denominator = den1 * den2
    
    LR = as.numeric(numerator / denominator)
    
  }  #
  
  return(LR)
}

calcLR_Lindley = function(control, recovered, background){
  ## Calculates the likelihood ratio for a univariate random effects with
  ## between items modelled normal this is taken from Lindley's 1977 work and
  ## really forms the precursor to the rest of this package. REQUIRES control -
  ## a compitem object calculated from the observations from the item considered
  ## to be the control item - calculated from makeCompItem from the file
  ## items_two_level.r recovered - a compitem object calculated from the
  ## observations from the item considered to be the recovered item - calculated
  ## from two.level.comparison.items() from the file items_two_level.r
  
  
  # first check to make sure that the observations are univariate
  if (background$multivariate) {
    stop("Data are multivariate - univariate only allowed for this function")
  }
  
  # U - is within variance C - is between variance mu - all means
  
  # redefine some of the items sent first the object with the population information
  U = background$v.within
  C = background$v.between
  mu = background$overall.means
  
  # then the objects with the control and recovered item information
  x = control$item.means
  m = control$n.replicates
  y = recovered$item.means
  n = recovered$n.replicates
  
  
  a = sqrt((1/m) + (1/n))
  w = ((m * x) + (n * y))/(m + n)
  
  delta.1 = C + (U/m)
  delta.2 = C + (U/n)
  delta.3 = C + (U/(n + m))
  z = ((delta.2 * x) + (delta.1 * y))/(delta.1 + delta.2)
  
  bit1 = (sqrt(delta.1) * sqrt(delta.2))/(a * sqrt(U) * sqrt(delta.3))
  bit2 = (((x - y)^2) * C)/(a^2 * U * (delta.1 + delta.2))
  bit3 = ((w - mu)^2)/(2 * delta.3)
  bit4 = (((z - mu)^2) * (delta.1 + delta.2))/(2 * delta.1 * delta.2)
  bit5 = bit4 - bit3
  
  LR = as.numeric(bit1 * exp(-bit2) * exp(bit5))
  
  return(LR)
} 
