#' Likelihood ratio calculation using Lindley's approach
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
#'   assuming that the between item distribution is univariate normal. This
#'   function is taken from the approach devised by Denis Lindley in his 1977
#'   paper (details below) and represents the progenitor of all the functions in
#'   this package.
#'
#' @return an estimate of the likelihood ratio
#' 
#' @references Lindley, D. (1977) A problem in forensic Science. \emph{Biometrika}: \bold{64}; 207-213.
#' @author David Lucy
#' 
#' @export
#'
#' @examples
#' # load Greg Zadora's glass data
#' data(glass)
#' 
#' # calculate a compcovar object based upon dat
#' # using K
#' Z = two.level.components(glass, 7, 1)
#' 
#' # calculate a compitem object representing the control item
#' control = two.level.comparison.items(glass[1:6,], 7)
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from the same item (item 1)
#' recovered.1 = two.level.comparison.items(glass[7:12,], 7)
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from a different item (item 2)
#' recovered.2 = two.level.comparison.items(glass[19:24,], 7)
#' 
#' 
#' # calculate the likelihood ratio for a known
#' # same source comparison - should be 6.323941
#' # This value is 6.323327 in this version and in the last version written by David (1.0-4)
#' lr.1 = two.level.lindley.LR(control, recovered.1, Z)
#' lr.1
#'
#' # calculate the likelihood ratio for a known
#' # different source comparison - should be 0.004422907
#' # This value is 0.004421978 in this version and the last version written by David (1.0-4)
#' lr.2 = two.level.lindley.LR(control, recovered.2, Z)
#' lr.2
two.level.lindley.LR = function(control, recovered, background) {
    ## Calculates the likelihood ratio for a univariate random effects with between
    ## items modelled normal this is taken from Lindley's 1977 work and really forms
    ## the precursor to the rest of this package REQUIRES control - a compitem
    ## object calculated from the observations from the item considered to be the
    ## control item - calculated from two.level.comparison.items() from the file
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
}  #


