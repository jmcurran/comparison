\name{comparison-package}
\alias{comparison-package}
%\alias{comparison}
\docType{package}
\title{Calculates and compares likelihood ratios for uni/multivariate comparison problems}

\description{Takes uni/multivariate observations of a physical object, for example trace element measurements, and a similar set of uni/multivariate observations from another object, and calculates a likelihood ratio for the propositions that the first item came from the same source as the second given some populational data.}


\details{
  \tabular{ll}{
  Package: \tab comparison\cr
  Type: \tab Package\cr
  Version: \tab 1.0\cr
  Date: \tab 11-07-2012\cr
  License: \tab GPL2\cr
  LazyLoad: \tab yes\cr
  }
}

\author{David Lucy \email{d.lucy@lancaster.ac.uk} - \url{http://www.maths.lancs.ac.uk/~lucy}}

\section{Acknowledgements}{In a package of functions such as these which have undergone a long development over a number of years, it is inevitable that a number of people, besides those directly cited, have helped to correct and add to the code. These people are (in alphabetical order):Ivo Alberink (NFI), Anabel Bolck (NFI), Sonja Menges (BKA), Geoff Morrison (UNSW), Tereza Neocleous (Glasgow), Anders Nordgaard (SKL), Brad Patterson (George Mason), Phil Rose (ANU), Agnieszka Rzepecka (Jagiellonian), Marjan Sjerps (NFI) and Hanjing Zhang (Edinburgh).}


\references{Aitken, C.G.G. & Lucy, D. (2004) Evaluation of trace evidence in the form of multivariate data. \emph{Applied Statistics}: \bold{53}(1); 109-122.}

\keyword{package}

\seealso{
\code{\linkS4class{compcovar}}\cr 
\code{\linkS4class{compitem}}\cr
\code{\link{two.level.comparison.items}}\cr
\code{\link{two.level.components}}\cr
\code{\link{two.level.density.LR}}\cr
\code{\link{two.level.normal.LR}}\cr
\code{\link{two.level.lindley.LR}}\cr
\code{\link[isotone]{gpava}}\cr
\code{\linkS4class{ece}}\cr
\code{\link{calibrate.set}}\cr
\code{\link{calc.ece}}\cr
\code{\link{plot}}
}


\examples{
ls()
}

\keyword{multivariate}

