##The function 'LR.Nor.function' calculates likelihood ratio results in comparison problems assuming a normal between-object distribution.

LR.Nor.function = function(y.mean.1, y.mean.2, y.star, mean.all, U, C, n.1, n.2, p)
{
	##Numerator calculation

	nom1 = (2*pi)^(-p)*exp(-1/2*(y.mean.1-y.mean.2) %*% solve(U/n.1+U/n.2) %*% t(y.mean.1-y.mean.2))*(det(U/n.1+U/n.2))^(-1/2)
	nom2 = exp(-1/2*(y.star-mean.all) %*% solve(U/(n.1+n.2)+C) %*% t(y.star-mean.all))*(det(U/(n.1+n.2)+C))^(-1/2)
	nom = nom1*nom2

	##Denominator calculation

	denom1 = (2*pi)^(-p)*exp(-1/2*(y.mean.1-mean.all) %*% solve(U/n.1+C) %*% t(y.mean.1-mean.all))*(det(U/n.1+C))^(-1/2)
	denom2 = exp(-1/2*(y.mean.2-mean.all) %*% solve(U/n.2+C) %*% t(y.mean.2-mean.all))*(det(U/n.2+C))^(-1/2)

	denom = denom1*denom2

	LR.Nor = nom/denom

	result = list(LR.Nor = LR.Nor)
	return (result)
}

