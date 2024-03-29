#This function convert time from UTC to LST
UTC2LST <-  function(ts.i,lon.i)
{
        #ts.i is the input time based on UTC
        #lon.i is the input longitude
        ts.o    <-  ts.i + as.numeric(lon.i)*3600/15
    return(ts.o)
}

#This function convert time from LST to UTC
LST2UTC <-  function(ts.i,lon.i)
{
        #ts.i is the input time based on LST
        #lon.i is the input longitude
        ts.o    <-  ts.i - lon.i*3600/15
    return(ts.o)
}

#This function estimates the sun rise and sun set time (LST) for given location
#the date.i must be in Date format
GetSunRiseSet  <-  function(date.i,lon.i,lat.i)
{
    #date.utc    <-  LST2UTC(date.i)
    v.tmp   <-  getSunlightTimes(date.i,lat=lat.i,
                                 lon=lon.i,
                                 keep=c("sunrise","sunset"))
    if (is.na(v.tmp$sunrise))
    {
        #polar night or day
        mon.t   <-  as.numeric(format(date.i,'%m'))
        if (abs(mon.t-7) < min(abs(mon.t - 1),abs(mon.t - 13)))
        {
            return(ifelse(lat.i < 0, -1, 1))
        } else
            return(ifelse(lat.i < 0, 1, -1))
    } else
    {
        #sunrise
        srs  <-  c(as.PCICt(v.tmp$sunrise,cal="365_day"),
                   as.PCICt(v.tmp$sunset,cal="365_day"))
        return(UTC2LST(srs,lon.i))
    }
}

GetTimeNum  <-  function(time.i){
    h.o    <-  as.numeric(format(time.i,'%H'))
    m.o    <-  as.numeric(format(time.i,'%M'))
    s.o    <-  as.numeric(format(time.i,'%S'))

    return(h.o*3600+m.o*60+s.o)
}
