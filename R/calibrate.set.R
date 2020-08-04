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
#' 
#' @author David Lucy
#'
#' @return a `list` with two items: \describe{
#'   \item{LR.cal.ss}{calibrated LRs for the comparison for same set}
#'   \item{LR.cal.ds}{calibrated LRs for the comparison for different set}
#' }
#' 
#' @references D. Ramos and J. Gonzalez-Rodrigues, (2008) "Cross-entropy analysis of the information in forensic speaker recognition," in Proc. IEEE Odyssey, Speaker Lang. Recognit. Workshop.
#' 
#' @seealso [isotone::gpava()], [calc.ece()]
#' @export
calibrate.set = function(LR.ss, LR.ds, method = c("raw", "laplace")) {
    
    method = match.arg(method)

    # correction for sorting - possibly not needed for R's sorting leave it in anyway
    LR.ss = LR.ss - 1e-06
    
    # get the lengths of the LR vectors and prior vector
    n.ss = length(LR.ss)
    n.ds = length(LR.ds)
    n = n.ss + n.ds
    
    # set an arbitary cutoff to prevent division by zero errors but bodgy this and I'm not too keen on it - however so long as the value is
    # close to, but less than, 1 it makes little difference cutoff = 0.99999999999999999999999999
    
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
    
    
    # gpava doesn't really care what you send it as the first argument so long as it is ascending and of the appropriate length - here I just
    # sent an integer array
    if (method == "raw") {
        calibrated.set = gpava(1:n, ordered.indicator.array)
    }
    if (method == "laplace") {
        calibrated.set = gpava(1:(n + 4), ordered.indicator.array)
    }
    
    calibrated.posterior.probabilities = calibrated.set$x
    
    # disentangle all the bits and pieces
    if (method == "laplace") {
        calibrated.posterior.probabilities = calibrated.posterior.probabilities[3:(n + 2)]
        ordered.indicator.array = ordered.indicator.array[3:(n + 2)]
    }
    
    
    # prior odds
    prior.odds = n.ss / n.ds
    
    # arbitrarily set a cutoff to prevent division by zero errors when we
    # calculate the odds not very keen on this - don't think it is needed
    # calibrated.posterior.probabilities[calibrated.posterior.probabilities > cutoff] = cutoff
    
    log.calibrated.posterior.LRs = log(calibrated.posterior.probabilities/(1 - calibrated.posterior.probabilities)) - log(prior.odds)
    
    # Brummer adds this in the ensure the idempotent property of the logLRs - not too sure how important this is
    bit.to.add.on = 1:n * 1e-06/n
    log.calibrated.posterior.LRs = log.calibrated.posterior.LRs + bit.to.add.on
    calibrated.posterior.LRs = exp(log.calibrated.posterior.LRs)
    
    # unpack all the calibrated LR values
    LR.cal.ss = calibrated.posterior.LRs[ordered.indicator.array == 1]
    LR.cal.ds = calibrated.posterior.LRs[ordered.indicator.array == 0]
    
    out = list(LR.cal.ss = LR.cal.ss, LR.cal.ds = LR.cal.ds)
    
    # debug code assign('kk', calibrated.posterior.LRs, .GlobalEnv)
    
    return(out)
}

