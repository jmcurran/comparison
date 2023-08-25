#' Calculate the calibrated set of LRs with the logistic regression
#'
#' @param LR.ss a vector of likelihood ratios for the comparisons of items known to be from the same source
#' @param LR.ds a vector of likelihood ratios for the comparisons of items known to be from different sources
#'
#' @author Marco De Donno
#'
#' @seealso [logistic.apply.calibration()]
#'
#' @return a `list` with multiple items: \describe{
#'   \item{prior.odds}{prior odds for the input data}
#'   \item{coefficients}{coefficients of the fitted model}
#'   \item{data}{The input and calibrated data}
#'   \item{LR.cal.ss}{The calibrated data for the same source list}
#'   \item{LR.cal.ds}{The calibrated data for the different-sources list}
#' }
#'
#' @examples
#' LR.same = c(0.5, 2, 4, 6, 8, 10)              # the list of LRs for the same source proposition
#' LR.different = c(0.2, 0.4, 0.6, 0.8, 1.1)     # the list of LRs for the different source proposition
#' logistic.calibrate.set(LR.same, LR.different) # compute the logistic calibration on the data
#'
#' @export
logistic.calibrate.set = function(LR.ss, LR.ds) {
    fit = logistic.calibrate.get.model(LR.ss, LR.ds)
    
    # compute the calibrated value for all input LRs
    LR.cal.ss = logistic.apply.calibration(LR.ss, fit)
    LR.cal.ds = logistic.apply.calibration(LR.ds, fit)
    
    out = list(
        LR.cal.ss = LR.cal.ss,
        LR.cal.ds = LR.cal.ds
    )
    
    return(out)
}

