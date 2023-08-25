#' Create a `compitem` object.
#'
#' This function creates a `compitem` from a set of observations
#' on items to be deemed control, or a recovered, items. For example,
#' a set of elemental concentration measurements on a sample of glass
#' fragments taken from a crime scene source such as a window.
#'
#' @param x a `matrix` or `data.frame` or a `formula`
#' @param data if \code{x} is a formula, then the user must supply a \code{data.frame}
#' containing the observations.
#' @param \dots other arguments that may be passed to the function.
#'
#' @return an object of class `compitem`
#' 
#' @author David Lucy and James Curran
#' 
#' @importFrom stats complete.cases
#' @importFrom methods is
#' @export
#'
#' @examples
#' # load Greg Zadora's glass data
#' data(glass)
#'
#' # calculate a compitem object representing the control item
#' controlMeasurements = subset(glass, item == "s1", select = c(logKO, logCaO, logFeO))
#' control = makeCompItem(controlMeasurements)
#' 
#' # example using the formula interface
#' controlMeasurements = subset(glass, item == "s1")
#' control = makeCompItem(item ~ logKO + logCaO + logFeO, data = controlMeasurements)
#'  
makeCompItem = function(x, ...){
  UseMethod("makeCompItem")
}

#' @export
makeCompItem.default = function(x, ...){
  if(!(is(x, "matrix") || is(x, "data.frame"))){
    stop("x must be a matrix or a data.frame")
  }
  
  warn.type = "none"
  
  x = as.matrix(x, rownames.force = TRUE)
  
  # clean the data a bit - get rid of NA rows - crude and may lead to cases
  # with n < 2 which is tested for later
  if (any(is.na(x))) {
    cc = complete.cases(x)
    numRemoved = sum(!cc)
    x = x[cc, ]
    
    msg = paste0("x contains missing values\n", numRemoved, " cases have been removed")
    warning(msg, immediate. = FALSE, call. = FALSE)
    warn.type = "NAs"
  }
  
  # simple counts of the data
  item.means = colMeans(x)
  n.replicates = nrow(x)
  n.vars = ncol(x)
  
  # test for, and set the multivariate flag
  multivariate.flag = n.vars > 1

  # tidy up the rownames property of the observations
  rownames(x) = 1:n.replicates
  
  
  #return(new("compitem", item.means = item.means, n.replicates = n.replicates, n.vars = n.vars, multivariate = multivariate.flag, observed = observed, 
  #    warn.type = warn.type))
  
  obj = list(item.means = item.means, 
             n.replicates = n.replicates, 
             n.vars = n.vars, 
             multivariate = multivariate.flag, 
             observed = x, 
             warn.type = warn.type)
  class(obj) = "compitem"
  return(obj)
}

#' @describeIn makeCompItem Create a `compitem` object using a formula.
#' @export
makeCompItem.formula = function(x, data =  NULL, ...){
    if(missing(x) || !is(x, "formula")){
      stop("missing or incorrect formula")
    }
    
    form = x

    mf = model.frame(form, data)
    
    group = model.response(mf)
    variables = mf[,-1, drop = FALSE]
    
    makeCompItem(variables, ...)
}