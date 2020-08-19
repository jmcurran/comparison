#' Glass composition data for seven elements from 200 glass items.
#'
#' These data are from Grzegorz (Greg) Zadora at the [Institute of Forensic Research](http://ies.krakow.pl/) in 
#' Krakow, Poland. They are the log of the ratios of each
#' element to oxygen, so logNaO is the log(10) of the Sodium to Oxygen ratio,
#' and logAlO is the log of the Aluminium to Oxygen ratio. The instrumental
#' method was SEM-EDX.
#'
#' The `item` indicates the object the glass came from. The levels for each item
#' are unique to that item. The `fragment` can be considered a sub-item. When
#' collecting these observations Greg took a glass object, say a jam jar, he
#' would then break it, and extract four fragments. Each fragment would be
#' measured three times upon different parts of that fragment. The fragment
#' labels are repeated, so, for example, fragment "f1" from item "s2" has
#' nothing whatsoever to do with fragment "f1" from item "s101".
#'
#' For two level models use `item` as the lower level - three level models can
#' use the additional information from the individual fragments.
#'
#' @name glass
#' @docType data
#' 
#' @format a `data.frame` with 2400 rows and 9 columns.
#' \describe{
#'   \item{item}{factor}{200 levels - which item the measurements came from}
#'   \item{fragment}{factor}{4 levels - which of the four fragments from each item the observations were made
#'   upon}
#'   \item{logNaO}{numeric}{log of sodium concentration to oxygen concentration}
#'   \item{logMgO}{numeric}{log of magnesium concentration to oxygen concentration}
#'   \item{logAlO}{numeric}{log of aluminium concentration to oxygen concentration}
#'   \item{logSiO}{numeric}{log of silicon concentration to oxygen concentration}
#'   \item{logKO}{numeric}{log of potassium concentration to oxygen concentration}
#'   \item{logCaO}{numeric}{log of calcium concentration to oxygen concentration}
#'   \item{logFeO}{numeric}{log of iron concentration to oxygen concentration}
#' }
#' @usage 
#' data(glass)
#' @references Aitken, C.G.G. Zadora, G. & Lucy, D. (2007) A Two-Level Model for
#'   Evidence Evaluation. \emph{Journal of Forensic Sciences}: \bold{52}(2);
#'   412-419.
#' @source Grzegorz Zadora [Institute of Forensic Research](http://ies.krakow.pl/), Krakow, Poland.
#' @keywords datasets
"glass"
NULL
