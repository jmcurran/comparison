\name{ece-class}

\docType{class}

\alias{ece-class}

\title{Class of object \code{ece}}

\description{Objects of class \code{"ece"} are empirical cross-entropies.}
 
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("ece", ...)} or, more commonly, via the \code{\link{calc.ece}} function which is recommended.}

\section{Slots}{
\describe{
\item{\code{prior}:}{an arbitary array of ascending values between zero and one representing each value for which the \code{ece} is to be evaluated}
\item{\code{ece.null}:}{\code{ece} values calculated assuming the observations give precisely no information about the membership of the }
\item{\code{ece}:}{\code{ece} values calculated from the likelihood ratios derived from the observations which are related to membership}
\item{\code{ece.cal}:}{values calculated from the likelihood ratios derived from a calibrated set using the function \code{calibrate_set}}
}
}

\section{Methods}{
\describe{
\item{plot}{variation on the generic \code{plot} function to plot up all three sets of ece values}
}
}

\references{Ramos, D. & Gonzalez-Rodriguez, J. (2008) Cross-entropy analysis of the information in forensic speaker recognition; IEEE Odyssey.

Zadora, G. & Ramos, D. (2010) Evaluation of glass samples for forensic purposes - an application of likelihood ratio model and information-theoretical approach. \emph{Chemometrics and Intelligent Laboratory}: \bold{102}; 63-83.
}

\seealso{
\code{\link[isotone]{gpava}}\cr
\code{\linkS4class{ece}}\cr
\code{\link{calibrate.set}}\cr
\code{\link{calc.ece}}\cr
\code{\link{plot}}
}

\keyword{classes}
