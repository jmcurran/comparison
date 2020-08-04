#' Empirical cross-entropy (ECE) calculation
#'
#' Calculates the empirical cross-entropy (ECE) for likelihood ratios from a
#' sequence same and different item comparisons.
#'
#' ### Acknowledgements
#' The function to calculate the values of the likelihood ratio for the
#' `calibrated.set` draws heavily upon the `opt_loglr.m` function from
#' Niko Brummer's FoCal package for Matlab.
#' 
#' @seealso [isotone::gpava()], [calibrate.set()]
#'
#' @param LR.ss a vector of likelihood ratios (LRs) from same source
#'   calculations
#' @param LR.ds a vector of LRs from different source calculations
#' @param prior a vector of ordinates for the prior in ascending order, and
#'   between 0 and 1. Default is 99 divisions of 0.01 to 0.99.
#'
#' @author David Lucy
#'
#' @importFrom isotone gpava
#'
#'
#'
#' @references @references D. Ramos and J. Gonzalez-Rodrigues, (2008) "Cross-entropy analysis of the information in forensic speaker recognition," in Proc. IEEE Odyssey, Speaker Lang. Recognit. Workshop.
#'   Zadora, G. & Ramos, D. (2010) Evaluation of glass samples for forensic purposes -  an application of likelihood ratio model and information-theoretical
#'   approach. Chemometrics and Intelligent Laboratory: 102; 63-83.
#'
#' @return Returns an S3 object of class ```ece```
#'
#' @examples
#' LR.same = c(0.5, 2, 4, 6, 8, 10) 		# the same has 1 LR < 1
#' LR.different = c(0.2, 0.4, 0.6, 0.8, 1.1) 	# the different has 1 LR > 1
#' ece.1 = calc.ece(LR.same, LR.different)	# simplest invocation
#' plot(ece.1)					# use plot method
#' @export
calc.ece = function(LR.ss, LR.ds, prior = seq(from = 0.01, to = 0.99, length = 99)) {
    # get the lengths of the LR vectors and prior vector
    n.ss = length(LR.ss)
    n.ds = length(LR.ds)
    n.prior = length(prior)
    
    # convert the prior to prior-odds
    odds = prior/(1 - prior)
    
    # set up arrays of null LRs to calculate the null ECE
    LR.null.ss = rep(1, n.ss)
    LR.null.ds = rep(1, n.ds)
    
    # calculate a set of calibrated LRs
    cal.set = calibrate.set(LR.ss, LR.ds, method = "raw")
    LR.cal.ss = cal.set$LR.cal.ss
    LR.cal.ds = cal.set$LR.cal.ds
    
    # set up the empirical cross entropy vector and the null vector and calibrated array
    ECE = NULL
    ECE.null = NULL
    ECE.cal = NULL
    
    # for all values of the prior
    for (ctr in 1:n.prior) {
        bit.1 = prior[ctr]/n.ss
        # for all same source LRs - do as a vector - sum later
        bit.2a = log2(1 + (1/(LR.ss * odds[ctr])))
        bit.2b = log2(1 + (1/(LR.null.ss * odds[ctr])))
        bit.2c = log2(1 + (1/(LR.cal.ss * odds[ctr])))
        
        bit.3 = (1 - prior[ctr])/n.ds
        # for all different source LRs - do as a vector - sum later
        bit.4a = log2(1 + (LR.ds * odds[ctr]))
        bit.4b = log2(1 + (LR.null.ds * odds[ctr]))
        bit.4c = log2(1 + (LR.cal.ds * odds[ctr]))
        
        # the ECE for the ctr'th value of the prior and the CE for the null for each value in the prior Ramos IEEE 2008 Equation 15
        ECE[ctr] = (bit.1 * sum(bit.2a)) + (bit.3 * sum(bit.4a))
        ECE.null[ctr] = (bit.1 * sum(bit.2b)) + (bit.3 * sum(bit.4b))
        ECE.cal[ctr] = (bit.1 * sum(bit.2c)) + (bit.3 * sum(bit.4c))
    }
    
    # debug code assign('aa', prior.ratio[ctr], .GlobalEnv)
    
    # S3 classes - function now S4 put together a list of the prior and ECEs out = list(prior, ECE, ECE.null, ECE.cal) names(out) =
    # c('prior', 'ECE', 'ECE.null', 'ECE.cal') and return it return(out)
    
    #return(new("ece", prior = prior, ece.null = ECE.null, ece = ECE, ece.cal = ECE.cal))
    
    result = list(prior = prior,
                  ece.null = ECE.null,
                  ece = ECE,
                  ece.cal = ECE.cal)
    class(result) = "ece"
    return(result)
}

