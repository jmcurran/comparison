#' Create a `compitem` object.
#'
#' This function creates a `compitem` object from a `data.frame` or `matrix` of
#' observations from an item to be deemed a control, or a recovered, item.
#'
#' @param data a `matrix` or `data.frame` of observed properties from either
#' the control item, or the recovered item
#' @param data.columns vector of integers giving which columns in `data` are
#'   the observations of the properties
#'
#' @return an object of class `compitem`
#' 
#' @importFrom stats complete.cases
#' @export
#'
#' @examples
#' # load Greg Zadora's glass data
#' data(glass)
#'
#' # calculate a compitem object representing the control item
#' control = two.level.comparison.items(glass[1:6,], c(7,8,9))
two.level.comparison.items = function(data, data.columns) {
    ## A `compitem` object has the data from
    ## the item to be compared - typically for any specific comparison problem there
    ## will be two of these objects requires: dat a data frame of variables and
    ## factors data.columns integer vector indicating which columns in dat contain
    ## the measurements returns: item.means real p array where p is the number of
    ## variables item.n integer: number of replicate observations for each item
    ## n.vars integer: number of continuous variables multivariate logical `TRUE` or `FALSE`
    ## indicates whether multivariate or not observed the raw observations in the
    ## form of an i*p matrix where there are i replicated observations and p
    ## variables missing values: traps and removes `NA`s - issues warnings about `NA`s
    
    # set the warn.type as none - then add in as warnings accrue
    warn.type = "none"
    
    # clean the data a bit - get rid of NA rows - crude and may lead to cases with n<2 which is tested for later
    if (any(is.na(data))) {
        warning("data contains NAs - cases removed", immediate. = FALSE, call. = FALSE)
        warn.type = "NAs"
        data = data[complete.cases(data), ]
    }
    
    
    # requires a cast to matrix for some reason
    observed = as.matrix(data[, data.columns], rownames.force = TRUE)
    
    # simple counts of the data
    item.means = apply(observed, 2, mean)
    n.replicates = nrow(observed)
    n.vars = length(data.columns)
    
    # test for, and set the multivariate flag
    if (n.vars > 1) {
        multivariate.flag = TRUE
    } else {
        multivariate.flag = FALSE
    }
    
    # tidy up the rownames property of the observations
    rownames(observed) = as.character(1:n.replicates)
    
    
    #return(new("compitem", item.means = item.means, n.replicates = n.replicates, n.vars = n.vars, multivariate = multivariate.flag, observed = observed, 
    #    warn.type = warn.type))
    
    obj = list(item.means = item.means, 
               n.replicates = n.replicates, 
               n.vars = n.vars, 
               multivariate = multivariate.flag, 
               observed = observed, 
               warn.type = warn.type)
    class(obj) = "compitem"
    return(obj)
    
}

