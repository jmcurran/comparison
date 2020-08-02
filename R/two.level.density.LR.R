#' Calculate the likelihood ratio using multivariate KDEs
#'
#'
#' Takes a `compitem` object which represents some control item, and a
#' `compitem` object which represents a recovered item, then uses information
#' from a `compcovar` object, which represents the
#' information from the population, to calculate a likelihood ratio as a measure
#' of the evidence given by the observations for the same/different source
#' propositions.
#' 
#' @param control a `compitem` object with the control item information.
#' @param recovered a `compitem` object with the recovered item information.
#' @param background a `compcovar` object with the population information.
#'
#' @return an estimate of the likelihood ratio
#' 
#' @references Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.
#' 
#' @export
#'
#' @examples
#' library(comparison)
#' # load Greg Zadora's glass data
#' data(glass)
#' 
#' # calculate a compcovar object based upon glass
#' # using K, Ca and Fe - warning - could take time
#' # on slower machines
#' Z = two.level.components(glass, c(7,8,9), 1)
#' 
#' # calculate a compitem object representing the control item
#' control = two.level.comparison.items(glass[1:6,], c(7,8,9))
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from the same item (item 1)
#' recovered.1 = two.level.comparison.items(glass[7:12,], c(7,8,9))
#' 
#' # calculate a compitem object representing the recovered item
#' # known to be from a different item (item 2)
#' recovered.2 = two.level.comparison.items(glass[19:24,], c(7,8,9))
#' 
#' 
#' # calculate the likelihood ratio for a known
#' # same source comparison - should be 20.59322
#' # 2020-08-01 Both this version and the previous version return 20.58967
#' lr.1 = two.level.density.LR(control, recovered.1, Z)
#' lr.1
#' 
#' # calculate the likelihood ratio for a known
#' # different source comparison - should be 0.02901532
#' # 2020-08-01 Both this version and the previous version return 0.01161392
#' lr.2 = two.level.density.LR(control, recovered.2, Z)
#' lr.2
two.level.density.LR = function(control, recovered, background) {
    ## Calculates the likelihood ratio for a multivariate random effects with
    ## between items modelled as kernel densities This function could still do with
    ## being made a bit quicker I have tried using apply() type formulations where
    ## appropriate but they are slightly slower than the counted iterations - I
    ## suspect the fact that they are dealing with matrix operations has something
    ## to do with it we're using the prod(eigen(A)$values) form rather than det(A)
    ## as it seems to be a little more reliable as some of the matricies tend to
    ## singularity REQUIRES control - a compitem object calculated from the
    ## observations from the item considered to be the control item - calculated
    ## from two.level.comparison.items() from the file items_two_level.r recovered -
    ## a compitem object calculated from the observations from the item considered
    ## to be the recovered item - calculated from two.level.comparison.items() from
    ## the file items_two_level.r
    
    
    
    ## calculates the difference between two numbers intended for use as an
    ## apply() functionette x is a matrix from whose rows we wish to subtract y
    ## y is a vector
    minu = function(x, y) {
        y - x
    }
    
    
    # convert to local definitions names changed to conform to three level code
    control.mean = control$item.means
    recovered.mean = recovered$item.means
    Nc = control$n.replicates
    Nr = recovered$n.replicates
    
    n.variables = background$n.vars
    
    
    
    # try to trap any errors here
    if (n.variables > 1) {
        multivariate.flag = TRUE
    } else {
        multivariate.flag = FALSE
    }
    if (!exists("multivariate.flag")) {
        stop("undefined number of variables")
    }
    
    n.groups = background$n.items
    group.means = background$item.means
    
    U = background$v.within
    C = background$v.between
    
    
    # window width calculation was formerly calculated by a separate function
    h.opt = (((4 / ((2 * n.variables) + 1)) ^ (1 / (n.variables + 4))) * (n.groups^(-(1 / (n.variables + 4)))))
    
    # print(h.opt)
    
    # stop()
    
    # debug code assign('bit', control.mean, env=.GlobalEnv)
    
    # do some essential type redefinition Nc = as.numeric(Nc) Nr = as.numeric(Nr) control.mean = as.matrix(control.mean) recovered.mean =
    # as.matrix(recovered.mean)
    
    
    # stop('systematic halt')
    
    
    if (multivariate.flag) {
        D.control = U / Nc
        D.recovered = U / Nr
        
        if (!is.numeric(try(solve(U)))) {
            value = "NA"
            return(value)
        }
        # if(! is.numeric(try(solve(C)))){value = 'NA'; return(value)}
        
        # component way of getting the inverses of the D matrices
        U.inv = solve(U)
        inv.D.control = U.inv * Nc
        inv.D.recovered = U.inv * Nr
        
        control.minus.recovered = control.mean - recovered.mean
        A = inv.D.control + inv.D.recovered
        
        # clever way of calculating the inverse of A
        inv.A = U / (Nc + Nr)
        
        # assign('A', A, .GlobalEnv) # useful bit of debug code
        
        # print(control.mean) assign('A', control.mean, .GlobalEnv) assign('B', inv.D.control, .GlobalEnv)
        
        
        y.star = inv.A %*% ((inv.D.control %*% (control.mean)) + (inv.D.recovered %*% (recovered.mean)))
        
        
        
        ##
        
        # print(prod(eigen(C)$values)) print(det(C)) stop('HHHH')
        
        
        
        top1 = sqrt(abs(prod(eigen(C)$values)))
        
        # trap dodgy covariance matricies
        if (prod(eigen(C)$values) < 0) {
            warning(
                "negative determinant - taking absolute value",
                call. = TRUE,
                immediate. = TRUE
            )
        }
        
        
        top2 = n.groups * (h.opt ^ n.variables)
        inv.h.opt.squared.times.C = solve((h.opt ^ 2) * C)
        
        top3 = 1 / sqrt(prod(abs(
            eigen(A + inv.h.opt.squared.times.C)$values
        )))
        
        
        top4 = exp(
            -0.5 * (control.minus.recovered) %*% solve(D.control + D.recovered) %*% (control.minus.recovered)
        )
        matt1 = rep(0, n.groups)
        
        # stop()
        
        # this invocation of minu makes a little difference to speed of execution
        y.star.minus.mean = t(apply(group.means, 1, minu, y = y.star))
        numerator.constant = solve(inv.A + (h.opt ^ 2) * C)
        
        
        for (ctr in 1:n.groups) {
            matt1[ctr] = exp(
                -0.5 * t(y.star.minus.mean[ctr,]) %*% numerator.constant %*% y.star.minus.mean[ctr,]
            )
        }
        
        
        
        
        top5 = sum(matt1)
        
        numerator = prod(top1, top2, top3, top4, top5)
        ##
        
        
        
        
        ##
        bot1 = 1 / sqrt(abs(prod(
            eigen(inv.D.control + inv.h.opt.squared.times.C)$values
        )))
        
        matt2 = rep(0, n.groups)
        control.konstant = solve(D.control + ((h.opt ^ 2) * C))
        
        
        
        for (ctr in 1:n.groups) {
            matt2[ctr] = exp(
                -0.5 * (control.mean - group.means[ctr,]) %*% control.konstant %*% (control.mean - group.means[ctr,])
            )
        }
        
        
        ## the matrix parts can equally well be done as the following apply() type of formulation - unfortunately despite its cleverness it is
        ## slower than the iterative approach ma2 = t(apply(group.means, 1, minu, y=control.mean)) konstant = solve(D.control + ((h.opt ^ 2) * C))
        ## ma3 = apply(ma2, 1, mulp, y=konstant) ma3 = exp(-0.5 * ma3) ma3 is now equal to matt2
        
        bot2 = sum(matt2)
        bot3 = 1 / sqrt(abs(prod(
            eigen(inv.D.recovered + inv.h.opt.squared.times.C)$values
        )))
        
        matt3 = rep(0, n.groups)
        recovered.konstant = solve(D.recovered + ((h.opt ^ 2) * C))
        
        # assign('cc', inv.A, envir=.GlobalEnv)
        
        
        
        for (ctr in 1:n.groups) {
            matt3[ctr] = exp(
                -0.5 * (recovered.mean - group.means[ctr,]) %*% recovered.konstant %*% (recovered.mean - group.means[ctr,])
            )
        }
        
        
        
        bot4 = sum(matt3)
        
        denominator = prod(bot1, bot2, bot3, bot4)
        
        LR = numerator / denominator
        return(LR)
        ##
    }  #
    
    
    
    
    
    
    
    
    
    
    
    
    if (!multivariate.flag)
    {
        h = h.opt
        k = nrow(group.means)
        a.sq = (1 / Nc) + (1 / Nr)
        
        
        w = ((Nc * control.mean) + (Nr * recovered.mean)) / (Nc + Nr)
        
        K.num = k * sqrt(Nc + Nr) * sqrt(U + (Nc * C * (h ^ 2))) * sqrt(U + (Nr * C * (h ^
                                                                                           2)))
        K.den = sqrt(a.sq) * sqrt(U) * sqrt(Nc * Nr) * sqrt(U + ((Nc + Nr) * C * (h ^
                                                                                      2)))
        K = K.num / K.den
        
        
        bit1 = ((control.mean - recovered.mean) ^ 2) / (2 * a.sq * U)
        bit2 = 2 * (U + ((Nc + Nr) * C * (h ^ 2)))
        bit3 = 2 * (U + (Nc * C * (h ^ 2)))
        bit4 = 2 * (U + (Nr * C * (h ^ 2)))
        
        num1 = 0
        den1 = 0
        den2 = 0
        
        for (ctr in 1:k) {
            tmp = ((Nc + Nr) * ((w - group.means[ctr]) ^ 2)) / bit2
            num1 = num1 + exp(-tmp)
            
            tmp = (Nc * ((control.mean - group.means[ctr]) ^ 2)) / bit3
            den1 = den1 + exp(-tmp)
            
            tmp = (Nr * ((
                recovered.mean - group.means[ctr]
            ) ^ 2)) / bit4
            den2 = den2 + exp(-tmp)
        }
        
        numerator = K * exp(-bit1) * num1
        denominator = den1 * den2
        
        LR = as.numeric(numerator / denominator)
        
    }  #
    
    
    
    
    return(LR)
}  #
