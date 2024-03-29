#' Select the proper time series of interest from the input data. It is a pre-processing process of the DFF framework.

#' @param var.i the input 3D variable (lon,lat,time).
#' @param lon.i longitude.
#' @param lat.i latitude.
#' @param ts.i  time series in the format of Date.
#' @param lon.range the range of longitude of interest (default=c(-180,180)).
#' @param lat.range the range of latitude of interest (default=c(-30,30)).
#' @param north.time.range the range of time of interest in the north  hemisphere (default is JJA).
#' @param south.time.range the range of time of interest in the south hemisphere (default is DJF).
#' @return a list: var: lon: masked longitude; lat: masked latitude; var: 3D array of variable.
#' @export

SelRegTime  <-  function(var.i,lon.i,lat.i,ts.i,
                         lon.range=c(-180,180),
                         lat.range=c(-30,30),
                         north.time.range=c('06-01','08-31'),
                         south.time.range=c('12-01','02-28')){

    #collect the region of interest first
    ind.lon.t   <-  which(lon.i >= lon.range[1] &
                          lon.i <= lon.range[2],
                          arr.ind=T)
    ind.lat.t   <-  which(lat.i >= lat.range[1] &
                          lat.i <= lat.range[2],
                          arr.ind=T)
    var.sub     <-  var.i[ind.lon.t,ind.lat.t,]
    
    lon.o   <-  lon.i[ind.lon.t]
    lat.o   <-  lat.i[ind.lat.t]
    
    #select the desired time series
    ind.lat.s   <-  which(lat.o < 0, arr.ind=T)
    ind.lat.n   <-  which(lat.o >= 0, arr.ind=T)

    #get the time index of north.time.range
    doy.n   <-  as.numeric(format(as.Date(paste0('2000-',
                                                 north.time.range)),
                                  '%j'))
    doy.s   <-  as.numeric(format(as.Date(paste0('2000-',
                                                 south.time.range)),
                                  '%j'))

    doy.ts  <-  as.numeric(format(as.Date(ts.i),'%j'))

    if (doy.s[2] < doy.s[1])
    {
        ind.ts.s    <- which(doy.ts >= doy.s[1] | doy.ts <= doy.s[2],
                             arr.ind=T)
    } else
        ind.ts.s    <- which(doy.ts >= doy.s[1] & doy.ts <= doy.s[2],
                             arr.ind=T)

    if (doy.n[2] < doy.n[1])
    {
        ind.ts.n    <- which(doy.ts >= doy.n[1] | doy.ts <= doy.n[2],
                             arr.ind=T)
    } else
        ind.ts.n    <- which(doy.ts >= doy.n[1] & doy.ts <= doy.n[2],
                             arr.ind=T)

    len.ts.o    <-  max(length(ind.ts.s),length(ind.ts.n))

    v.o <-  array(NA,dim=c(length(lon.o),length(lat.o),len.ts.o))

    v.o[,1:length(ind.lat.s),1:length(ind.ts.s)]    <-
        var.sub[,ind.lat.s,ind.ts.s]
    v.o[,(length(ind.lat.s)+1):length(lat.o),1:length(ind.ts.n)] <-
        var.sub[,ind.lat.n,ind.ts.n]

    return(list(lon=lon.o,lat=lat.o,var=v.o))
}
