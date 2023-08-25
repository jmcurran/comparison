#' Calculate the calibrated set of idea LRs
#'
#' Calculates and returns the calibrated set of `ideal' LRs from the
#' observed LRs using the penalised adjacent violators algorithm. This is very
#' much a rewrite of Nico Brummer's `optloglr()` function for Matlab.
#' 
#' This is an internal function, and is not meant to be called directly. However
#' it has been exported just in case.
#' 
#'
#' @param LR.ss a vector of likelihood ratios for the comparisons of items known to be from the same source
#' @param LR.ds a vector of likelihood ratios for the comparisons of items known to be from different sources
#' @param method the method used to perform the calculation, either `"raw"` or `"laplace"`
#' @param ties method to solve ties in the predictors list, either `"none"` (not solved) or `"primary"`, `"secondary"` or `"tertiary"` (passed to the isotone::gpava() function)
#' 
#' @author David Lucy
#'
#' @return a `list` with two items: \describe{
#'   \item{LR.cal.ss}{calibrated LRs for the comparison for same set}
#'   \item{LR.cal.ds}{calibrated LRs for the comparison for different set}
#' }
#' 
#' @references Ramos, D. & Gonzalez-Rodriguez, J. (2008) Cross-entropy analysis of the information in forensic speaker recognition; IEEE Odyssey.
#' @references de Leeuw, J. & Hornik, K. & Mair, P., (2009), Isotone Optimization in R: Pool-Adjacent-Violators Algorithm (PAVA) and Active Set Methods, https://www.jstatsoft.org/article/view/v032i05
#' 
#' @seealso [isotone::gpava()], [calc.ece()]
#' @export
calibrate.set = function(LR.ss, LR.ds, method = c("raw", "laplace"), ties = c("none", "primary", "secondary", "tertiary")) {
    
    method = match.arg(method)
    ties = match.arg(ties)

    # get the lengths of the LR vectors and prior vector
    n.ss = length(LR.ss)
    n.ds = length(LR.ds)
    n = n.ss + n.ds
    
    # produce vectors for the indicator function and LRs use convention that 1 signifies same set, and 0 different set
    indicator.array = c(rep(0, length = n.ds), rep(1, length = n.ss))
    LR.array = c(LR.ds, LR.ss)
    
    # order the vectors
    ordered.indicies = order(LR.array)
    ordered.indicator.array = indicator.array[ordered.indicies]
    
    # implement Laplace's rule of sucession in the way of Brummer
    if (method == "laplace") {
        ordered.indicator.array = c(1, 0, ordered.indicator.array, 1, 0)
    }
    
    # compute the calibrated.set values with the gpava function
    if (ties == "none") {
        # this function call dont take into account ties values in the input data.
        # this is the default behavior if the user dont specify a method to solve the ties values.
        calibrated.set = gpava(
            seq(length(ordered.indicator.array)),
            ordered.indicator.array
        )
        
    } else {
        # this function call the gpava function taking into account the parameter passed to solve the ties values
        # check the publication of Leeuw et al (2009)
        ordered.LR.array = LR.array[ordered.indicies]
        
        if (method == "laplace") {
            ordered.LR.array = c(0, 0, ordered.LR.array, Inf, Inf)
        }
        
        calibrated.set = gpava(
            ordered.LR.array,
            ordered.indicator.array,
            ties = ties
        )
    }
    
    calibrated.posterior.probabilities = calibrated.set$x
    
    # disentangle all the bits and pieces
    if (method == "laplace") {
        calibrated.posterior.probabilities = calibrated.posterior.probabilities[3:(n + 2)]
        ordered.indicator.array = ordered.indicator.array[3:(n + 2)]
    }
    
    
    # compute the calibrated LR based upon the prior odds and the observed data
    prior.odds = n.ss / n.ds
    
    log.calibrated.posterior.LRs = log(calibrated.posterior.probabilities/(1 - calibrated.posterior.probabilities)) - log(prior.odds)
    
    calibrated.posterior.LRs = exp(log.calibrated.posterior.LRs)
    
    # unpack all the calibrated LR values
    LR.cal.ss = calibrated.posterior.LRs[ordered.indicator.array == 1]
    LR.cal.ds = calibrated.posterior.LRs[ordered.indicator.array == 0]
    
    out = list(LR.cal.ss = LR.cal.ss, LR.cal.ds = LR.cal.ds)
    
    # debug code assign('kk', calibrated.posterior.LRs, .GlobalEnv)
    
    return(out)
}

