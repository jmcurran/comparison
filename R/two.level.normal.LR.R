#' Likelihood ratio calculation - normal
#'
#' Takes a `compitem` object which represents some control item, and a
#' `compitem` object which represents a recovered item, then uses information
#' from a `compcovar` object, which represents the information from the
#' population, to calculate a likelihood ratio as a measure of the evidence
#' given by the observations for the same/different source propositions.
#'
#' @param control a `compitem` object with the control item information
#' @param recovered a `compitem` object with the recovered item information
#' @param background a `compcovar` object with the population information
#'
#' @details Does the likelihood ratio calculations for a two-level model
#'   assuming that the between item distribution is uni/multivariate normal.
#'
#' @return an estimate of the likelihood ratio
#' 
#' @author Agnieszka Martyna and David Lucy
#' 
#' @references Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.
#' @export
#' @examples
#' # load Greg Zadora's glass data
#' data(glass)
#' 
#' # calculate a compcovar object based upon glass
#' # using K, Ca and Fe - warning - could take time
#' # on slower machines
#' Z <- two.level.components(glass, c(7,8,9), 1)
#' 
#' # calculate a compitem object representing the control item
#' control <- two.level.comparison.items(glass[1:6,], c(7,8,9))
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from the same item (item 1)
#' recovered.1 <- two.level.comparison.items(glass[7:12,], c(7,8,9))
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from a different item (item 2)
#' recovered.2 <- two.level.comparison.items(glass[19:24,], c(7,8,9))
#' 
#' 
#' # calculate the likelihood ratio for a known
#' # same source comparison - should be 51.16539
#' # This value is 51.14243 in this version and the last version David wrote (1.0-4)
#' lr.1 <- two.level.normal.LR(control, recovered.1, Z)
#' lr.1
#' # calculate the likelihood ratio for a known
#' # different source comparison - should be 0.02901532
#' # This vsalue is 0.02899908 in this version and the last version David wrote (1.0-4)
#' lr.2 <- two.level.normal.LR(control, recovered.2, Z)
#' lr.2
two.level.normal.LR = function(control, recovered, background) {
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
    
    
    # This code independently written by Agnieszka Rzepecka of the IFR, and tidied up a bit by David Lucy output agrees with previous efforts -
    # thus we develop some confidence it's doing the right thing I keep on expecting it to fall over in the univariate case but it seems not to
    # so we'll use it until it starts to be a problem Numerator calculation
    nom1 = exp(-1/2 * t(diff.cont.rec) %*% solve(U/n.cont + U/n.rec) %*% (diff.cont.rec)) * (det(U/n.cont + U/n.rec))^(-1/2)
    nom2 = exp(-1/2 * t(diff.y.star.mu) %*% solve(U/(n.cont + n.rec) + C) %*% (diff.y.star.mu)) * (det(U/(n.cont + n.rec) + C))^(-1/2)
    nom = nom1 * nom2
    
    # Denominator calculation
    denom1 = exp(-1/2 * t(diff.cont.mu) %*% solve(U/n.cont + C) %*% (diff.cont.mu)) * (det(U/n.cont + C))^(-1/2)
    denom2 = exp(-1/2 * t(diff.rec.mu) %*% solve(U/n.rec + C) %*% (diff.rec.mu)) * (det(U/n.rec + C))^(-1/2)
    denom = denom1 * denom2
    
    LR = as.numeric(nom/denom)
    
    ## original code from 2004 which follows the Aitken & Lucy paper (p.115) closely ## this is here so people trying to follow the code with
    ## the paper can do so ## D1 = U / n.cont D2 = U / n.rec H2 = t(y.star - mu) %*% solve(U/(n.rec + n.cont) + C) %*% (y.star - mu) H3 =
    ## t(diff.cont.rec) %*% solve(D1 + D2) %*% (diff.cont.rec) num1 = exp(- 0.5 * (H2 + H3)) num2 = sqrt(det(2 * pi * solve(((n.cont + n.rec)
    ## * solve(U) ) + solve(C)))) num = num1 * num2 mu.star = solve(solve(D1 + C) + solve(D2 + C)) %*% ((solve(D1 + C) %*% cont.means) +
    ## (solve(D2 + C) %*% rec.means)) H4 = t(mu - mu.star) %*% (solve(D1 + C) + solve(D2 + C)) %*% (mu - mu.star) H5 = t(diff.cont.rec) %*%
    ## solve(D1 + D2 + (2 * C)) %*% diff.cont.rec den1 = 1 / sqrt(det(2 * pi * C)) den2 = sqrt(det(2 * pi * solve(n.cont * solve(U) +
    ## solve(C)))) den3 = sqrt(det(2 * pi * solve(n.rec * solve(U) + solve(C)))) den4 = exp(-0.5 * (H4 + H5)) lr = num / prod(den1, den2,
    ## den3, den4)
    
    return(LR)
} 

