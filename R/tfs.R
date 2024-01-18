#' Triggering feedback strength.
#' K. Findell, P. Gentine, B. R. Linter, and C. Kerr, 2011: Probability of afternoon precipitation in eastern United States and Mexico enhanced by high  evaporation. Nat. Geosci., 4, 434-439, https://doi.org/10.1038 /ngeo1174.

#' @param ef.am time series of evaporative fraction in the early morning.
#' @param pre.am precipitation in the morning [mm]
#' @param pre.pm precipitation in the afternoon [mm]
#' @param pre.thr the threshold to identify preciptation events (=0.1 mm by default)
#' @param n.bin number of bins (=10 by default).
#' @return the value of TFS 

#' @export
TFS <-  function(ef.am,pre.am,pre.pm,pre.thr=.1,n.bin=10){
    #select available dates
    #select days without morning rainfall
    idx.nona  <-  seq(1,length(ef.am))[!is.na(ef.am) &
                                       !is.na(pre.am) &
                                       !is.na(pre.pm) &
                                       pre.am < pre.thr]
    ef.t  <-  ef.am[idx.nona]
    pre.t <-  pre.pm[idx.nona]
    pev.t <-  pre.t
    pev.t[pre.t >= pre.thr] <-  1
    pev.t[pre.t < pre.thr] <-  0
    len.tot <-  length(pre.t)

    q.ef  <-  as.numeric(quantile(ef.t,
                                  probs=seq(0,1,length.out=n.bin+1)))
    r.o <-  array(NA,dim=n.bin)
    e.o <-  array(NA,dim=n.bin)
    ef.o <-  array(NA,dim=n.bin)

    for (nz in 1:n.bin)
    {
      ind.t  <-  seq(1,len.tot)[ef.t > q.ef[nz] &
                                  ef.t <= q.ef[nz+1]]
      ef.o[nz] <-  mean(ef.t[ind.t])
      r.o[nz]  <-  sum(pev.t[ind.t])/length(ind.t)
      pre.tt <-  pre.t[ind.t]
      e.o[nz]  <-  mean(pre.tt[pre.tt >= pre.thr])
    }

    e.o[is.na(e.o)] <-  0
    q.ef.t  <-  (q.ef[1:n.bin]+q.ef[2:(n.bin+1)])/2

    sum.tfs <-  0
    for (nz in 2:n.bin)
        sum.tfs <-  sum.tfs + (r.o[nz] - r.o[nz-1])/
                    (n.bin*(ef.o[nz] - ef.o[nz-1]))

    return(sum.tfs*sd(ef.t))
}
