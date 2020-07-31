\name{compitem-class}
\Rdversion{1.1}
\docType{class}
\alias{compitem-class}
%\alias{compitem,compitem-method}

\title{Class of object \code{compitem}}

\description{An object which represents the information required for comparison calculated from either the control, or recovered items.}


\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("compitem", ...)}. Using \code{\link{two.level.comparison.items}} is recommended.
}

\section{Slots}{
\describe{
\item{\code{item.means}:}{a vector of the means for the observations of the item - object of class \code{"numeric"}}
\item{\code{n.replicates}:}{count of the number of replicated observations of the item - object of class \code{"numeric"}}
\item{\code{n.vars}:}{count of the number of propoerties observed for the item - object of class \code{"numeric"}}
\item{\code{multivariate}:}{flag indicating whether the observations are multivariate - object of class \code{"logical"}}
\item{\code{observed}:}{the raw observations for the item - object of class \code{"matrix"}}
\item{\code{warn.type}:}{whether any warnings were issued \code{"character"}}
}
}




\section{Methods}{
\describe{
\item{two.level.density.LR}{\code{signature(control = "compitem", recovered = "compitem", background = "compcovar")}: ... }
}
}


\references{Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.}

\author{David Lucy \email{d.lucy@lancaster.ac.uk} - \url{http://www.maths.lancs.ac.uk/~lucy}}

\note{Best not to create these directly, use \code{\link{two.level.comparison.items}} instead.}

\seealso{
\code{\linkS4class{compcovar}}\cr 
\code{\linkS4class{compitem}}\cr
\code{\link{two.level.comparison.items}}\cr
\code{\link{two.level.components}}\cr
\code{\link{two.level.density.LR}}
}

\examples{
showClass("compitem")
}

\keyword{classes}
\keyword{multivariate}


