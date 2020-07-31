\name{compcovar-class}
\Rdversion{1.1}
\docType{class}
\alias{compcovar-class}
%\alias{compcovar,compcovar-method}

\title{Class of object \code{compcovar}}

\description{A set of estimates for the hiearchical covariance matricies, and means, calculated for a background population needed for the comparison of two items.}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("compcovar", ...)}. The function \code{\link{two.level.components}} is recommended for the creation of these objects.
}

\section{Slots}{
  \describe{
    \item{\code{v.within}:}{estimate of the within item covariance matrix - object of class \code{"matrix"}.}
    \item{\code{v.between}:}{estimate of the between item covariance matrix - object of class \code{"matrix"}}
    \item{\code{n.observations}:}{count of the total number of observations in the sample from the background population - object of class \code{"numeric"}}
    \item{\code{n.items}:}{count of the number of items in the sample from the background population of - object class \code{"numeric"}}
    \item{\code{item.n}:}{an array of the number of replicated observations for each item in the sample - object of class \code{"matrix"}}
    \item{\code{item.means}:}{an array of the means of the replicated observations for each item in the sample - object of class \code{"matrix"}}
    \item{\code{n.vars}:}{count of the number of properties measured - object of class \code{"numeric"} ~~ }
    \item{\code{overall.means}:}{means of all the properties measured for each item in the background population - object of class \code{"numeric"}}
    \item{\code{multivariate}:}{flag indicating whether the observations are multivariate - object of class \code{"logical"}}
    \item{\code{balanced}:}{flag indicating whether the observations are balanced - Object of class \code{"logical"}}
    \item{\code{s.within}:}{sum of squared deviations for within all items - object of class \code{"matrix"}}
    \item{\code{s.between}:}{sum of squared deviations between all items - object of class \code{"matrix"}}
    \item{\code{warn.type}:}{whether any warnings were issued \code{"character"}}
  }
}



\section{Methods}{
  \describe{
    \item{two.level.density.LR}{\code{signature(control = "compitem", recovered = "compitem", background = "compcovar")}}
	 }
}

\references{Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.}

\author{David Lucy \email{d.lucy@lancaster.ac.uk} - \url{http://www.maths.lancs.ac.uk/~lucy}}

\note{ML estimates for variance components - need something a bit more robust}

\seealso{
\code{\linkS4class{compcovar}}\cr 
\code{\linkS4class{compitem}}\cr
\code{\link{two.level.comparison.items}}\cr
\code{\link{two.level.components}}\cr
\code{\link{two.level.density.LR}}
}



\examples{showClass("compcovar")}


\keyword{classes}
\keyword{multivariate}

