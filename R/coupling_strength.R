#' Evaluate the coupling strength between two variables. The strength is determined by their correlation (causality) and amplitude (influence).

#' @param x a numeric vector representing time series of the driver 
#' @param y a numeric vector representing time series of the passive variable
#' @return a list containing (1) the coupling strenghs and its significance
#' @export
coupling_strength <-  function(x,y)
{
    #remove NA from the time series
    ind.na  <-  is.na(x*y)

    #check if available observations are sufficient
    if (length(ind.) < 50)
        warning('coupling_strength: available observations are too short.')

    x.t <-  x[!ind.na]
    y.t <-  y[!ind.na]

    r.t <-  cor.test(x.t,y.t,
                     method=method)

    s.y <-  sd(y.t)

    couple.o    <-  list(coupling=s.y*
                         as.numeric(r.t$estimate),
                         p.value=as.numeric(r.t$p.value))
    return(couple.o)
}

