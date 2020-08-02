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
calcLR = function(control, recovered, background, method = c("mvn","kde","lindley")){
  
}

#' @export
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
  cont.means = control$item.means
  n.cont = control$n.replicates
  rec.means = recovered$item.means
  n.rec = recovered$n.replicates
  
  # weighted mean for the control and recovered item
  y.star = ((n.cont * cont.means) + (n.rec * rec.means))/(n.cont + n.rec)
  
  # calculate some of the repeated components
  diff.cont.rec = cont.means - rec.means
  diff.cont.mu = cont.means - mu
  diff.rec.mu = rec.means - mu
  diff.y.star.mu = y.star - mu
  
  
  # This code independently written by Agnieszka Rzepecka of the IFR, and tidied
  # up a bit by David Lucy. The output agrees with previous efforts - thus we develop
  # some confidence it's doing the right thing I keep on expecting it to fall
  # over in the univariate case but it seems not to so we'll use it until it
  # starts to be a problem 
  
  ## Numerator calculation Cleaned up and hopefully sped up by James Curran
  k1 = (n.cont * n.rec) / (n.cont + n.rec)
  k2 = (1 / k1)^nrow(U)
  logNumerator1 = k1 * t(diff.cont.rec) %*% solve(U) %*% (diff.cont.rec) + log(k2 * det(U))
  logNumerator2 = t(diff.y.star.mu) %*% solve(U/(n.cont + n.rec) + C) %*% (diff.y.star.mu) + log(det(U/(n.cont + n.rec) + C))
  logNumerator = logNumerator1 + logNumerator2
  
  # Denominator calculation
  logDenominator1 = t(diff.cont.mu) %*% solve(U/n.cont + C) %*% (diff.cont.mu) + log(det(U/n.cont + C))
  logDenominator2 = t(diff.rec.mu) %*% solve(U/n.rec + C) %*% (diff.rec.mu) + log(det(U/n.rec + C))
  logDenominator = logDenominator1 + logDenominator2
  
  LR = as.numeric(exp(-0.5 *(logNumerator - logDenominator)))
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