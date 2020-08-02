#' S3 method for class `compitem`
#'
#' @param x an object of class `compitem` created by [makeCompItem()]
#' @param ... further arguments passed to or from other methods.
#'
#' @export
print.compitem = function(x, ...){
  cat(paste0("Comparison Item: ", deparse(substitute(x)), "\n"))
  cat(paste0("No. Obs. : ", x$n.replicates, "\n"))
  cat(paste0("No. Vars : ", x$n.vars, "\n"))
  cat(paste0("Warnings : ", paste0(x$warn.type, collapse = ", "), "\n\n"))
}