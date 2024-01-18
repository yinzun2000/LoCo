#' @export
print.htlm  <-  function(x, prefix="\t", ...) {
    cat('\n')
    cat(strwrap('Two-legged coupling strength between', prefix = prefix), sep = "\n")
    cat(strwrap('the land and the atmosphere (LoCo)', prefix = prefix), sep = "\n")
    cat('\n')
    cat("Data: ",x$data.name,"\n",sep="")
    cat("The land leg: ",x$land,"; with p-value = ",x$land.s,"\n",sep="")
    cat("The atmos leg: ",x$atmos,"; with p-value = ",x$atmos.s,"\n",sep="")
    cat("The total leg: ",x$tot,"\n",sep="")
    invisible(x)
}
