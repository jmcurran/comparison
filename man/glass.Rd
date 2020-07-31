\name{glass}
\alias{glass}
\docType{data}
\title{Data - elemental composition of glass}

\description{Greg Zadora's glass composition data for seven elements from 200 glass items.}

\usage{data(glass)}

\format{
\tabular{lll}{
\code{item} \tab factor - \tab 200 levels - which item the measurements came from\cr
\code{fragment} \tab factor - \tab 4 levels - which of the four fragments from each item the observations were made upon\cr
\code{logNaO} \tab numeric - \tab log of sodium concentration to oxygen concentration\cr
\code{logMgO} \tab numeric - \tab log of magnesium concentration to oxygen concentration\cr
\code{logAlO} \tab numeric - \tab log of aluminium concentration to oxygen concentration\cr
\code{logSiO} \tab numeric - \tab log of silicon concentration to oxygen concentration\cr
\code{logKO} \tab numeric - \tab log of potassium concentration to oxygen concentration\cr
\code{logCaO} \tab numeric - \tab log of calcium concentration to oxygen concentration\cr
\code{logFeO} \tab numeric - \tab log of iron concentration to oxygen concentration
}
}

\details{
These data are from Greg Zadora at the  \href{http://ies.krakow.pl/}{Institute of Forensic Research} in Cracow, Poland. They are the log of the ratios of each element to oxygen, so logNaO is the log(10) of the Sodium to Oxygen ratio, and logAlO is the log of the Aluminium to Oxygen ratio. The instrumental method was SEM-EDX.\cr

The item indicates the object the glass came from. The levels for each item are unique to that item. The fragment can be considered a sub-item. When collecting these observations Greg took a glass object, say a jam jar, he would then break it, and extract four fragments. Each fragment would be measured three times upon different parts of that fragment. The fragment labels are repeated, so, for example, fragment "f1" from item "s2" has nothing whatsoever to do with fragment "f1" from item "s101".\cr

These data are a \code{list}. If some object more matrix like is needed convert using \code{as.data.frame()}.\cr

For two level models use "item" as the lower level - three level models can use the additional information from the individual fragments.
}

\source{Grzegorz Zadora \href{http://ies.krakow.pl/}{Institute of Forensic Research} Cracow, Poland.}

\references{Aitken, C.G.G. Zadora, G. & Lucy, D. (2007) A Two-Level Model for Evidence Evaluation. \emph{Journal of Forensic Sciences}: \bold{52}(2); 412-419.}

\examples{
data(glass)
}

\keyword{datasets}
