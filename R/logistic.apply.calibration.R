#' Calculate the calibrated LRs with the model precomputed
#'
#' This function perform the logistic calibration on the provided data.
#' In the context of likelihood ratios, the 'ideal' value for the LR is Infinity for the same source dataset, and 0 for the different-sources dataset.
#' The 'post' values are fixed to 1 for the same source and 0 for the same different-sources datasets (corresponding to the posterior probability P(H_ss|E)).
#'
#' @param LR a vector of likelihood ratios to be calibrated (raw values).
#' @param model a logistic.calibrate.set() fitted model to be applied. This variable can be the reture of the logistic.calibrate.set() or the logistic.calibrate.set()$fit variable.
#'
#' @author Marco De Donno
#'
#' @seealso [logistic.calibrate.set()]
#'
#' @return a `list` with the calibrated LR values
#'
#' @examples
#'  # the list of LRs for the same source proposition
#' LR.same = c(0.5, 2, 4, 6, 8, 10)
#' # the list of LRs for the different source proposition
#' LR.different = c(0.2, 0.4, 0.6, 0.8, 1.1)
#' # compute the logistic calibration on the data
#' model = logistic.calibrate.get.model(LR.same, LR.different) 
#'  # the list of news LRs (to be calibrated)
#' LR.unknown = c(0.6, 0.7, 1.2, 5)
#' # compute the calibrated LRs for the list with the model
#' logistic.apply.calibration(LR.unknown, model)
#'
#' @export
logistic.apply.calibration = function(LR, model) {
    # extract the coefficients of the model variable
    coef = coefficients(model)
    
    # the predictors are the log10(LR), to be coherent with the model fitting function
    predictors = log10(LR)
    
    # compute the calibrated values
    calibrated.posterior.ratio = exp(predictors * coef[2] + coef[1])
    
    # compute the correct LR depending upon the prior odds
    calibrated.posterior.probabilities = calibrated.posterior.ratio/(calibrated.posterior.ratio + 1)
    calibrated.posterior.LRs = (calibrated.posterior.probabilities/(1 - calibrated.posterior.probabilities)) / (model$prior.odds)
    
    return(calibrated.posterior.LRs)
}

