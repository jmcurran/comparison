#' Compute integrated means and covariances
#'
#' Takes a large sample from the background population and calculates the within
#' and between covariance matrices, a vector of means, a vector of the counts of
#' replicates for each item from the sample, and other bits needed to make up a
#' `compcovar` object.
#'
#' @param x a `matrix`, or `data.frame`, of observations, with cases in rows,
#'   and properties as columns, or a `formula`.
#' @param item.column item.column an integer scalar indicating which column contains 
#' the item label.
#' @param data if \code{x} is a \code{formula}, then the user must supply a \code{data.frame}
#' containing the observations.
#' @param \dots other arguments.
#' @param item.column an integer indicating which column gives the item.
#'
#' @details Uses ML estimation at the moment - this will almost certainly change in the future and
#' hopefully allow regularisation methods to get a more stable (and non-singular) estimate.
#' 
#' @author David Lucy and James Curran
#'
#' @return an object of class `compvar`
#' @importFrom stats coefficients glm model.frame model.response rnorm
#' @export
#'
#' @examples
#' # load Greg Zadora's glass data
#' data(glass)
#' 
#' # calculate a compcovar object based upon glass
#' # using K, Ca and Fe - warning - could take time
#' # on slower machines
#' background = subset(glass, select = c(item, logKO, logCaO, logFeO))
#' Z1 = makeCompVar(background, 1)
#' 
#' # Use the formula interface
#' Z2 = makeCompVar(item ~ logKO + logCaO + logFeO, data = glass)
makeCompVar = function(x, ...){
    UseMethod("makeCompVar")
}


#' @describeIn makeCompVar Create a `compvar` object using a formula.
#' @export
makeCompVar.default = function(x, item.column, ...) {
    ## Compute integrated means and covariances function does univariate and
    ## multivariate in one function. This function may be used for balanced and
    ## unbalanced data so long as the imbalance is not extreme. Requires: data, a data
    ## frame of variables and factors data.columns integer vector indicating which
    ## columns in data contain the measurements item.column integer scalar indicating
    ## which column contains the item labels: items are the top level of variability
    ## returns: v.within real p*p matrix where p is the number of variables,
    ## estimate of the within fragment covariance matrix v.between real p*p matrix
    ## where p is the number of variables - the between items covariance matrix
    ## n.observations total number of observations - 1*1 integer n.items number of
    ## items - 1*1 integer item.n integer: number of fragments for each item - of
    ## length n.fragments item.means real j*p matrix where j is the number of groups
    ## and p number of variables n.vars integer: number of contnuous variables
    ## overall.means p*1 real vector of means for all variables for all observations
    ## multivariate logical T/F indicates whether multivariate balanced logical T/F
    ## indicates whether the replication is balanced warn.type character - the type
    ## of warning(s) issued - crude for the moment notes: works by getting the
    ## nested means first, then calculating the covariance components from these
    ## using weighted estimates missing values: crudely handled by mean(x,
    ## na.rm=TRUE) to stop any NA's from crashing the function - not very well
    ## implemented in this version
    
    
    
    # set the warn.type as none - then add in as warnings accrue
    warn.type = "none"
    
    ## COMMON 
    ## clean the data a bit - get rid of NA rows - crude and may lead to
    ## cases with n<2 which is tested for later
    if (any(is.na(x))) {
        cc = complete.cases(x)
        numRemoved = sum(!cc)
        x = x[cc, ]
        
        msg = paste0("x contains missing values\n", numRemoved, " cases have been removed")
        warning(msg, immediate. = FALSE, call. = FALSE)
        warn.type = "NAs"
    }
    
    
    
    # definitions and declarations
    vars = names(x[,-item.column])
    n.vars = length(vars)
    
    multivariate.flag = n.vars > 1
    
    
    n.observations = nrow(x)
    
    items = unique(x[, item.column])
    n.items = length(items)
    
    item.n = matrix(0, nrow = n.items, ncol = 1)
    rownames(item.n) = items
    colnames(item.n) = "n"
    
    item.means = matrix(0, nrow = n.items, ncol = n.vars)
    rownames(item.means) = items
    colnames(item.means) = vars
    
    
    # means for all cases not affected by imbalance - names attribute awkward
    # overall mean calculation now even more awkward as colMeans now fussy about
    # getting a matrix rather than a vector so colMeans OK if multivariate -
    # otherwise mean
    if (multivariate.flag) {
        overall.means = colMeans(x[, -item.column], na.rm = TRUE)
    }else{
        overall.means = mean(x[, -item.column], na.rm = TRUE)
        names(overall.means) = vars
    }
    # overall.means = colMeans(dat[,data.columns], na.rm=TRUE); if(!multivariate.flag){names(overall.means) = vars}
    
    
    
    # debugging assign('kk', overall.means, env=.GlobalEnv)
    
    
    
    s.w = matrix(0, nrow = n.vars, ncol = n.vars)
    row.names(s.w) = vars
    colnames(s.w) = vars
    s.star = s.w
    ## 
    
    ## split the data up into groups by item
    f = factor(x[,item.column], levels = unique(x[,item.column]))
    split.data = split(x[, -item.column], f)
    
    
    ## get the number of obs per group
    item.n = sapply(split.data, nrow)
    if(any(item.n < 2)){
        stop("There must be at two observations in each group/item")
    }
    
    
    ## multivariate ##
    if (multivariate.flag) {
        
        ## compute the group means
        item.means = lapply(split.data, colMeans, na.rm = TRUE)
        
        calcSW = function(item, mx){
            d = as.matrix(sweep(item, 2, mx))
            return(t(d) %*% d)
        }

        s.w = Reduce('+', mapply(calcSW, split.data, item.means, SIMPLIFY = FALSE))
        
        calcSB = function(item, mx){
            d = mx - overall.means
            #browser()
            return(nrow(item) * outer(d, d))
        }
        
        s.b = Reduce('+', mapply(calcSB, split.data, item.means, SIMPLIFY = FALSE))
        
        ## compute the group means
        item.means = do.call("rbind", item.means)

        # # pick out the observations for each level of the grouping factor then
        # # get the mean and the sum of squared deviations (S.w)
        # for (i in 1:n.items) {
        #     current.item = x[x[, item.column] == items[i], ]
        #     item.n[i] = nrow(current.item)
        # 
        #     # trap cases with too few replicates to calculate a mean from - fatal if it occurs
        #     # if (item.n[i] < 2) {
        #     #     stop("too few replicates in an item")
        #     # }
        # 
        #     #item.means[i, ] = apply(current.item[, -item.column], 2, mean, na.rm = TRUE)
        # 
        # 
        #     # for each observation calculate the SS dev component
        #     for (j in 1:item.n[i]) {
        #         s.w = s.w + (item.n[i] * ((as.numeric(current.item[j, -item.column] - item.means[i, ]) %*%
        #                                         t(as.numeric(current.item[j, -item.column] - item.means[i, ])))))
        #     }
        # 
        #     # the between sum of squared deviations between the item means and
        #     # overall mean - OK for unbalanced at the moment
        #     s.star = s.star + (item.n[i] * (((item.means[i, ] - overall.means) %*% t(item.means[i, ] - overall.means))))
        # }
        # browser()
    }
    ## 
    
    
    
    
    
    
    
    
    ## univariate ##
    if (!multivariate.flag) 
        # slightly different for univariate case
    {
        # define the outputs as zero
        s.w = 0
        s.star = 0
        
        for (i in 1:n.items) {
            current.item = x[x[, item.column] == items[i], ]
            item.n[i] = nrow(current.item)
            
            
            
            # differs from multivariate case
            item.means[i] = mean(current.item[, -item.column], na.rm = TRUE)
            
            
            
            # for each observation calculate the SS dev component
            for (j in 1:item.n[i]) {
                s.w = s.w + (item.n[i] * ((current.item[j, -item.column] - item.means[i])^2))
            }
            # the between sum of squared deviations between the item means and overall mean - OK for unbalanced at the moment
            s.star = s.star + (item.n[i] * ((item.means[i, ] - overall.means)^2))
        }
        
        s.w = as.matrix(s.w)
        rownames(s.w) = vars
        colnames(s.w) = vars
        s.star = as.matrix(s.star)
        rownames(s.star) = vars
        colnames(s.star) = vars
        
        
    }
    ## 
    
    
    # both s.w and s.b are pretty well certain to be correct as I have gotten the
    # same values from different independently written functions - the final
    # evaluation of C and U is less certain - these need a good checking 
    # givethe sums of squared deviations as part of the output 
    ##########
    # NOTE: JMC - I have removed these lines because I do not think
    # they are necessary or correct
    # s.w = s.w * (n.items/n.observations)
    # s.b = s.star * (n.items/n.observations)
    ############
    
    # set the status of whether the data are balanced between items
    balanced.flag = length(unique(item.n)) == 1

    
    
    #U = s.w/(n.observations - n.items)#(n.items * s.w)/(n.observations * (n.observations - n.items))
    #U = #U * (n.observations/n.items)  # this bit may be wrong - in so U agrees with previous code
    # NOTE: I have rewritten this 
    ## The WGMS = U = WGSS / (n.observations - n.items)
    U = s.w / (n.observations - n.items)
    
    
    # this may be the correct one C = ((s.star) / (n.observations / n.items * (n.items - 1))) - (s.w / ((n.observations ^ 2 / n.items ^ 2) *
    # (n.observations - n.items)))
    
    ## this one may also be wrong - in so U agrees with previous code 
    ## C = (s.b / (n.items - 1)) - (s.w / ((n.observations^2 / n.items) - n.items)) 
    ## this (below) is correct by A&L2004 - thanks to Hanjing Zhang and Colin Aitken for this revision
    ## C = (s.b/(n.items - 1)) - (s.w/((n.observations^2/n.items) - n.observations))
    ## NOTE: I (JMC) have rewritten this as it is wrong
    ## The BG Covariance = (BGMS - WGMS) / r where r is the number of observations per group. T
    ## This is ONLY true in the balanced case
    C = (s.b / (n.items - 1) - U) / (n.observations / n.items)
    
    # return(new("compcovar", v.within = U, v.between = C, n.observations = n.observations, n.items = n.items, item.n = item.n, item.means = item.means, 
    #     n.vars = n.vars, overall.means = overall.means, multivariate = multivariate.flag, balanced = balanced.flag, s.within = s.w, s.between = s.b, 
    #     warn.type = warn.type))
    
    obj = list(v.within = U, 
               v.between = C, 
               n.observations = n.observations, 
               n.items = n.items, 
               item.n = item.n, 
               item.means = item.means, 
               n.vars = n.vars, 
               overall.means = overall.means, 
               multivariate = multivariate.flag, 
               balanced = balanced.flag, 
               s.within = s.w, 
               s.between = s.b, 
               warn.type = warn.type)
    class(obj) = "compvar"
    return(obj)
}

#' @describeIn makeCompVar Create a `compvar` object using a formula.
#' @export
makeCompVar.formula = function(x, data =  NULL, ...){
    if(missing(x) || !is(x, "formula"))
        stop("missing or incorrect formula")
    
    form = x
    
    mf = model.frame(form, data)

    makeCompVar.default(x = mf, item.column = 1)
}
