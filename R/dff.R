#' The Duan-Findell-Fueglistaler (DFF) framework.
#' Reference: Duan S. Q., Findell K. L., and Fueglistaler S. A.: Coherent mechanistic patterns of land hydroclimatic changes, GRL, 50, e2022GL102285. https://doi.org/10.1029/2022GL102285

#' @param sm a 3D array of soil moisture (lon,lat,time).
#' @param indx a 3D array of specific index wanted to show in the DFF framework (lon,lat,time).
#' @param ai a 2D array of arid index (lon,lat).
#' @param lsm a 2D array of land mask (lon,lat). 1=land, 0=ocean.
#' @param lat a vector of latitude.
#' @param n.bin number of bins (default=50).
#' @return a 2D matrix (n.bin,n.bin) showing the patterns of area-weighted indx with x-axis the percentile of soil moisture and y-axis the percentile of arid index.
#' @export

DFF <-  function(sm,indx,ai,lsm,lat,n.bin=50){
    #get the dimension of longitude, latitude, and time series
    n.lon <-  dim(sm)[1]
    n.lat <-  dim(sm)[2]
    n.t   <-  dim(sm)[3]

    v.o <-  array(NA,dim=c(n.bin,n.bin))
    sm.t <-  array(NA,dim=c(n.lon,n.lat,n.bin))
    vv.t <-  array(NA,dim=c(n.lon,n.lat,n.bin))

    for (nx in 1:n.lon)
    for (ny in 1:n.lat)
    {
        if (lsm[nx,ny] < .5)
            next

        v.t <-  CalBin(sm[nx,ny,],indx[nx,ny,],n.bin=n.bin)
        sm.t[nx,ny,] <-  v.t$sm
        vv.t[nx,ny,] <-  v.t$v
    }

    vv.t  <-  AreaWeighted(vv.t,lat)
    sm.tt <-  matrix(sm.t,prod(dim(sm.t)[1:2]),dim(sm.t)[3])
    vv.tt <-  matrix(vv.t,prod(dim(vv.t)[1:2]),dim(vv.t)[3])
    ai.tt <-  matrix(ai,prod(dim(ai)[1:2]))

    sm.tt <-  sm.tt[lsm > .5,]
    vv.tt <-  vv.tt[lsm > .5,]
    ai.tt <-  ai.tt[lsm > .5]
    ind.ai  <-  as.numeric(CutQuantile(ai.tt,n.bin,labels=T))

    for (i in 1:n.bin)
        v.o[,i] <-  apply(vv.tt[ind.ai == i,],2,mean,na.rm=T)
    return(v.o)
}

CalBin  <-  function(sm.i,v.i,n.bin=50){
    sm.o  <-  as.vector(tapply(sm.i,CutQuantile(sm.i,n.bin),
                               mean,na.rm=T))
    v.o   <-  as.vector(tapply(v.i,CutQuantile(sm.i,n.bin),
                               mean,na.rm=T))
    v.oo  <-  list(sm=sm.o,v=v.o)
    return(v.oo)
}

CutQuantile <-  function(v.i,n.bin=50,labels=F){
    v.max <-  max(abs(v.i),na.rm=T)
    v.t <-  v.i + rnorm(length(v.i),0,sd=v.max*1e-5)
    v.o <-  cut(v.t, breaks=c(quantile(v.t,
                                       probs = seq(0,1,by=1/n.bin),
                                       na.rm = T)),
                include.lowest=TRUE,labels=seq(1,n.bin))
    return(v.o)
}


AreaWeighted    <-  function(v.i,lat.i){
    lat.w   <-  WeightAlongLat(lat.i)
    v.o <-  sweep(v.i,2,lat.w,`*`)
    return(v.o)
}

WeightAlongLat  <-  function(lat.i){
    lat.r   <-  lat.i * pi /180
    d.t     <-  (lat.r[2] - lat.r[1])/2
    ww.t    <-  1+cos(lat.r - d.t)*cos(lat.r + d.t)
    return(ww.t/sum(ww.t))
}
