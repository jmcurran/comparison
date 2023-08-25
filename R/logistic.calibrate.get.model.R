#' Compute and returns the logistic regression for a dataset
#'
#' @param LR.ss a vector of likelihood ratios for the comparisons of items known 
#' to be from the same source
#' @param LR.ds a vector of likelihood ratios for the comparisons of items known
#'  to be from different sources
#'
#' @author Marco De Donno
#'
#' @seealso [logistic.apply.calibration()]
#'
#' @return a `list` with multiple items: \describe{
#'   \item{coefficients}{coefficients of the fitted model}
#'   \item{prior.odds}{prior odds for the input data}
#' }
#' 
#' @importFrom stats binomial
#'
#' @examples
#' # the list of LRs for the same source proposition
#' LR.same = c(0.5, 2, 4, 6, 8, 10)
#' # the list of LRs for the different source proposition
#' LR.different = c(0.2, 0.4, 0.6, 0.8, 1.1)
#' # compute the logistic calibration on the data
#' logistic.calibrate.get.model(LR.same, LR.different) 
#'
#' @export
logistic.calibrate.get.model = function(LR.ss, LR.ds) {
    # preparation of the data.frame to store all the mendatory data to make the calibration
    LR.ss.dataframe = data.frame(
        lr   = LR.ss,
        post = 1
    )
    LR.ds.dataframe = data.frame(
        lr   = LR.ds,
        post = 0
    )
    data = rbind(
        LR.ss.dataframe,
        LR.ds.dataframe
    )
    
    # get the lengths of the LR vectors and prior vector
    n.ss = length(LR.ss)
    n.ds = length(LR.ds)
    n = n.ss + n.ds
    prior.odds = n.ss / n.ds
    
    # compute the predictor list
    data$loglr = log10(data$lr)
    
    # fit the logistic function on the data
    fit = glm(data$post ~ data$loglr, family = binomial(link = "logit"))
    
    # warn the user if the logistic regression is not increasing with increasing LRs
    # this can happen if the same-source and difference-sources are inverted, or if the model is not sufficiently performant
    if(fit$coefficients[2] < 0)
        warning( "The logistic regression is decreasing. Check that the LR.ss and LR.ds variables are in the right order (the LR values should be bigger for the same source that for the different source proposition) or that your model is sufficiently performant" )
    
    # prepare the return value
    out = list(
        coefficients = coefficients(fit),
        prior.odds   = prior.odds,
        fit          = fit
    )
    
    return(out)
}

