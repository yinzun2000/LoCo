#' Calculate two-legged matrix for evaluating land-atmosphere coupling strength.

#' @param land a numeric vector or 3D array representing time series of land status (e.g., soil moisture)
#' @param flux a numeric vector or 3D array representing time series of energy fluxes (e.g., latent heat flux)
#' @param atmos a numeric vector or 3D array representing time series of atmospheric status (e.g., height of planet boundary height)
#' @param sign.level a number 0~1 determining the threshold of correlation significance
#' @param parallel TRUE or FALSE (default). if parallel computation will be triggered.
#' @param core.frac a number 0~1 controling the maximum cores utilized for parallel computation (fraction of total cores available). It works only when parallel=T. The number of cores are contrainted by both core.frac and cor.max.
#' @param core.max a integer (=5 by default) controling the maximum cores utilized for parallel computation. It works only when parallel=T. The number of cores are contrainted by both core.frac and cor.max.
#' @return a list containing (1) the coupling strenghs of land, atmosphere, and total leg; (2) significance of land and atmosphere coupling
#' @export
two_legged_coupling <-  function(land,flux,atmos,sign.level=0.05,parallel=F,core.frac=.5,core.max=10)
{
    #check inputs types
    t.land  <-  loco_check_type(land)
    t.flux  <-  loco_check_type(flux)
    t.atmos <-  loco_check_type(atmos)

    if (!(t.land$t == 'n') |
        !(t.flux$t == 'n') |
        !(t.atmos$t == 'n'))
        stop('two_legged_coupling: all inputs must be numeric')

    #check if the dimensions of these inputs are the same
    if (!all(t.land$d == t.flux$d,
             t.flux$d == t.atmos$d))
        stop('two_legged_coupling: dimensions of inputs are mismatch.')

    #1-D time series
    if (length(t.land$d) == 1)
    {
        tlm.o   <-  tlm.1d(land,flux,atmos)
        dname   <-  paste(deparse(substitute(land)), ", ", deparse(substitute(flux)),", and", deparse(substitute(atmos)))
        tlm.o$data.name <-  dname
        class(tlm.o)    <-  'htlm'
    } else if (length(t.land$d) == 3) #3D [lon,lat,time]
    {
        tlm.o   <-  tlm.3d(land,flux,atmos,sign.level,parallel,core.frac=core.frac,core.max=core.max)
    } else
        stop('two_legged_coupling: inputs dimension must be 1D [time] or 3D [lon,lat,time].')

    return(tlm.o)
}

# This function calculate 1-D two-legged metric
# and provide detailed diagnostic summary
tlm.1d  <-  function(land,flux,atmos)
{

    #remove NA from the time series
    ind.na  <-  is.na(land*flux*atmos)

    #check if available observations are sufficient
    if (length(ind.na) < 50)
        stop('two_legged_coupling: no enough available observations (must no less than 50)')

    land.t  <-  land[!ind.na]
    flux.t  <-  flux[!ind.na]
    atmos.t <-  atmos[!ind.na]

    c.land    <-  cor.test(land.t,flux.t)
    c.atmos   <-  cor.test(flux.t,atmos.t)

    cor.land    <-  as.numeric(c.land$estimate)
    cor.atmos   <-  as.numeric(c.atmos$estimate)
    p.land    <-  c.land$p.value
    p.atmos   <-  c.atmos$p.value

    s.flux  <-  sd(flux.t)
    s.atmos <-  sd(atmos.t)


    tlm.o   <-  list(land=s.flux*cor.land,
                     atmos=s.atmos*cor.atmos,
                     tot=s.atmos*cor.atmos*cor.land,
                     land.s=p.land,
                     atmos.s=p.atmos,
                     data.name=NULL)
    return(tlm.o)
}


tlm.3d  <-  function(land,flux,atmos,sign.level=0.05,parallel=F,core.frac,core.max)
{

    nlon    <-  dim(land)[1]
    nlat    <-  dim(land)[2]
    ntime   <-  dim(land)[3]

    if (ntime < 50)
        stop('two_legged_coupling: no enough available observations (must no less than 50)')

    #values of TLM [lat,lon,(land,atmos,total,p.land,p.atmos)]

    if (parallel)
    {
        n.cores <-  detectCores(logical=F)
        cores   <-  min(core.frac*n.cores,core.max)
        cl  <-  makeCluster(cores)
        registerDoParallel(cl, cores=cores)

        tlm.v   <-  foreach(ny=1:nlat,
                            .combine='acomb')%dopar%
        {
            tlm.t   <-  array(NA,dim=c(nlon,5))
            for (nx in 1:nlon)
                tlm.t[nx,]   <-  tlm.fast(land[nx,ny,],
                                          flux[nx,ny,],
                                          atmos[nx,ny,])
            tlm.t
        }
        gc()
        stopImplicitCluster()
        stopCluster(cl)
    } else
    {
        tlm.v   <-  array(NA,dim=c(nlon,nlat,5))
        #set the progress bar
        pb.t    <-  txtProgressBar(min=0,max=nlon,style=2,
                                   width=50,char="=")
        for (nx in 1:nlon)
        {
            for (ny in 1:nlat)
            {
                tlm.v[nx,ny,]   <-  tlm.fast(land[nx,ny,],
                                             flux[nx,ny,],
                                             atmos[nx,ny,])
                                             
            }
            setTxtProgressBar(pb.t,nx)
        }

    }

    #threshold for significance
    tt  <-  qt(sign.level/2,ntime-2)
    r.th    <-  abs(tt)/sqrt(ntime-2+tt**2)

    if (parallel)
    {
        tlm.o   <-  list(land=tlm.v[,1,],
                         atmos=tlm.v[,2,],
                         tot=tlm.v[,3,],
                         s.land=abs(tlm.v[,4,]) > r.th,
                         s.atmos=abs(tlm.v[,5,]) > r.th)
    } else
        tlm.o   <-  list(land=tlm.v[,,1],
                         atmos=tlm.v[,,2],
                         tot=tlm.v[,,3],
                         s.land=abs(tlm.v[,,4]) > r.th,
                         s.atmos=abs(tlm.v[,,5]) > r.th)

    return(tlm.o)
}

# This function calculate 1-D two-legged metric
# with fast. Note that this function will not 
# filt NA in the time series for calculation.
tlm.fast  <-  function(land,flux,atmos)
{

    c.land  <-  cor(land,flux)
    c.atmos <-  cor(flux,atmos)

    tlm.o   <-  c(c.land*sd(flux),
                  c.atmos*sd(atmos),
                  c.land*c.atmos*sd(atmos),
                  c.land,c.atmos)
    return(tlm.o)
}

acomb   <-  function(...) abind(...,along=3)
