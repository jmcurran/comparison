LR.KDE.function = function(y.mean.1, y.mean.2, y.star, U, C, h, population, variables, p, m, n.1,
    n.2) {

    ## The function 'LR.KDE.function' calculates likelihood ratio results in comparison problems
    ## when KDE is used for probability density function estimation.

    ## Numerator calculation

    nom2 = (2 * pi)^(-p/2) * exp(-1/2 * (y.mean.1 - y.mean.2) %*% solve(U/n.1 + U/n.2) %*% t(y.mean.1 -
        y.mean.2)) * (det(U/n.1 + U/n.2))^(-1/2)

    nom1 = 0

    ## i runs through all objects from the population
    for (i in 1:m) {
        items = unique(population$Item)

        ## creating a matrix of measurements for the ith object
        ith.object = as.matrix(population[which(population$Item == items[i]), variables])

        ## calculating the 'object.mean'
        object.mean = matrix(apply(ith.object, 2, mean), nrow = 1)

        exp.1.1 = exp(-(y.star - object.mean) %*% (solve(U/(n.1 + n.2) + C * h^2)) %*% t(y.star -
            object.mean)/2)

        nom1 = nom1 + exp.1.1
    }

    nom2.1 = nom1/m

    exp.1.2 = (2 * pi)^(-p/2) * det(U/(n.1 + n.2) + C * h^2)^(-1/2)
    nom3 = nom2.1 * exp.1.2

    nom = nom2 * nom3

    ## Denominator calculation
    denom1 = 0

    for (i in 1:m) {
        items = unique(population$Item)
        ith.object = as.matrix(population[which(population$Item == items[i]), variables])
        object.mean = matrix(apply(ith.object, 2, mean), nrow = 1)

        exp.2.1 = exp(-(y.mean.1 - object.mean) %*% (solve(U/n.1 + C * h^2)) %*% t(y.mean.1 - object.mean)/2)

        denom1 = denom1 + exp.2.1
    }

    denom2 = denom1/m

    exp.2.2 = (2 * pi)^(-p/2) * det(U/n.1 + C * h^2)^(-1/2)
    denom3 = denom2 * exp.2.2

    denom4 = 0

    for (i in 1:m) {
        items = unique(population$Item)
        ith.object = as.matrix(population[which(population$Item == items[i]), variables])
        object.mean = matrix(apply(ith.object, 2, mean), nrow = 1)

        exp.2.3 = exp(-(y.mean.2 - object.mean) %*% (solve(U/n.2 + C * h^2)) %*% t(y.mean.2 - object.mean)/2)


        denom4 = denom4 + exp.2.3
    }

    denom5 = denom4/m

    exp.2.4 = (2 * pi)^(-p/2) * det(U/n.2 + C * h^2)^(-1/2)
    denom6 = denom5 * exp.2.4

    denom = denom3 * denom6

    LR.KDE = nom/denom

    result = list(LR.KDE = LR.KDE)
    return(result)
}























