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
    fit$prior.odds = prior.odds
    
    # warn the user if the logistic regression is not increasing with increasing LRs
    # this can happen if the same-source and difference-sources are inverted, or if the model is not sufficiently performant
    if(fit$coefficients[2] < 0)
        warning( "The logistic regression is decreasing. Check that the LR.ss and LR.ds variables are in the right order (the LR values should be bigger for the same source that for the different source proposition) or that your model is sufficiently performant" )
    
    # compute the calibrated value for all input LRs
    data$callr = logistic.apply.calibration(data$lr, fit)
    LR.cal.ss  = logistic.apply.calibration(LR.ss, fit)
    LR.cal.ds  = logistic.apply.calibration(LR.ds, fit)
    
    # prepare the return value
    outdata = subset(data, select = c(lr, post, callr))
    
    out = list(
        coefficients = coefficients(fit),
        prior.odds   = prior.odds,
        data         = outdata,
        LR.cal.ss    = LR.cal.ss,
        LR.cal.ds    = LR.cal.ds
    )
    
    return(out)
}

