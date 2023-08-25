## The function UC provides the distributional parameters for likelihood ratio computation in
## comparison problems. For such parameters estimation it is necessary to use a database of
## similar objects described by the same variables gathered during analyses of the same type.
## The function calculates the within-object variability (denoted by U) and between-object
## variability (denoted by C) of the variables of interest. The database consists of m objects
## for which n measurements were obtained as a result of the analysis. In the case of univariate
## data, the number of variables (p) is 1, whereas in the case of multivariate data, the number
## of variables is p>1. The files containing data MUST be organised in such a way that rows
## correspond to n measurements performed on m objects (which means nm rows) and columns
## correspond to p variables.

## Terms used in the function: population - a matrix (dimensions:(m*n) x p) modified from a
## database according to the jackknife procedure; variables - a vector indicating which columns
## are taken into account for calculations (columns depict the variables measured in the
## analysis).

UC = function(population, variables, p, n) {
    items = unique(population$Item)
    m = length(items)  ## m corresponds to the number of objects creating a population

    ## Defining S.star and Sw matrices initially filled with 0 at the beginning of the loops
    S.star = matrix(0, nrow = p, ncol = p)
    variable.names = colnames(population[, variables])
    rownames(S.star) = variable.names
    colnames(S.star) = variable.names

    Sw = matrix(0, nrow = p, ncol = p)
    rownames(Sw) = variable.names
    colnames(Sw) = variable.names

    ## Dealing with multivariate data (p>1)
    if (p > 1) {
        ## mean.all corresponds to the vector of means of all variables calculated using all
        ## measurements for all objects
        mean.all = matrix(apply(population[, variables], 2, mean), nrow = 1)
        colnames(mean.all) = variable.names

        ## object.mean corresponds to a matrix of means of all variables calculated using all
        ## measurements for each object
        object.mean = matrix(0, nrow = m, ncol = p)
        rownames(object.mean) = as.character(items)
        colnames(object.mean) = variable.names

        ## i runs through all objects from the population
        for (i in 1:m) {
            Sw2 = matrix(0, nrow = p, ncol = p)

            ## creating a matrix of measurements for the ith object
            ith.object = as.matrix(population[which(population$Item == items[i]), variables])

            ## calculating the object.mean
            object.mean = matrix(apply(ith.object, 2, mean), nrow = 1)

            ## j runs through n measurements for the chosen ith object
            for (j in 1:nrow(ith.object)) {
                xij.minus.xi = ith.object[j, ] - object.mean
                Sw2 = Sw2 + t(xij.minus.xi) %*% xij.minus.xi
            }

            ## creating S.star matrix
            S.star = S.star + t(object.mean - mean.all) %*% (object.mean - mean.all)

            ## creating Sw matrix
            Sw = Sw + Sw2
        }

        ## creating U and C matrices
        U = Sw/(m * (n - 1))
        C = S.star/(m - 1) - U/n
    }

    ## Dealing with univariate data (p=1)
    if (p == 1) {
        mean.all = matrix(mean(population[, variables], nrow = 1))
        colnames(mean.all) = variable.names

        object.mean = matrix(0, nrow = m, ncol = p)
        rownames(object.mean) = as.character(items)
        colnames(object.mean) = variable.names

        for (i in 1:m) {
            Sw2 = matrix(0, nrow = p, ncol = p)

            ith.object = as.matrix(population[which(population$Item == items[i]), variables])
            object.mean = matrix(mean(ith.object), nrow = 1)

            for (j in 1:nrow(ith.object)) {
                xij.minus.xi = ith.object[j, ] - object.mean
                Sw2 = Sw2 + t(xij.minus.xi) %*% xij.minus.xi
            }

            S.star = S.star + t(object.mean - mean.all) %*% (object.mean - mean.all)
            Sw = Sw + Sw2
        }

        U = Sw/(m * (n - 1))
        C = S.star/(m - 1) - U/n
    }

    result = list(U = U, C = C, mean.all = mean.all, object.mean = object.mean, m = m)
    return(result)
}
