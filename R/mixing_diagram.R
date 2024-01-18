#' Calculate the mixing diagram to illustrate the co-evolution of potential temperature and humidity during a entire day process.

#' @param sp a numeric vector of surface pressure [Pa]
#' @param t2m a numeric vector of 2m temperature [K]
#' @param q2m a numeric vector of 2m specific humidity [kg.kg-1]
#' @param shf a numeric vector  of upward sensible heat flux [W.m-2]
#' @param lhf a numeric vector of upward latent heat flux [W.m-2]
#' @param hbpl a numeric vector of planetary boundary layer height [m]
#' @param time information. It can be either in PCICt format (e.g., "1980-01-01 10:50:33") or in numeric format (hours, e.g., 12, 12.5, etc). If the numeric format is applied, there will no daytime information generated.
#' @param time.type whether the time is LST or UTC (default).
#' @param lon if the time is UTC
#' @return a list containing the evolution of humidity and temperature at a entire day circle.
#' @export
mixing_diagram  <-  function(sp,t2m,q2m,shf,lhf,hpbl,times,lon=NULL,lat=NULL,time.type='UTC',day.lst=NULL){
    Lv  <-  2.5e6
    cp  <-  1005.7
    Rd  <-  287.04
    ep  <-  .622
    n.t <-  length(sp)
    dt  <-  as.numeric(times[2] - times[1])

    if (time.type == 'UTC')
    {
        time.utc    <-  times
        time.lst    <-  UTC2LST(time.utc,lon)
    } else if (time.type == 'LST')
    {
        time.lst    <-  times
        time.utc    <-  LST2UTC(time.utc,lon)
    } else
        stop('mixing_diagram: the time format cannot be recognized.')

    #key outputs
    cp.t    <-  cp*t2m
    lv.q    <-  Lv*q2m
    rho <-  sp/(Rd*t2m*(1+q2m/ep)/(1+q2m))
    f.sfc.ts <-  shf*dt/(rho*hpbl)
    m.sfc.ts <-  lhf*dt/(rho*hpbl)

    if (class(times[1]) == "PCICt")
    #calculate sunrise and sunset
    {

        dates.lst   <- format(time.lst,"%Y-%m-%d")
        days.lst   <- unique(dates.lst)
        days.utc   <- unique(format(time.utc,"%Y-%m-%d"))
        num.day   <-  length(days.utc)
        idx.day <-  array(NA,dim=n.t)
        f.sfc   <-  array(NA,dim=num.day)
        m.sfc   <-  array(NA,dim=num.day)
        f.atm   <-  array(NA,dim=num.day)
        m.atm   <-  array(NA,dim=num.day)

        for (i.d in 1:num.day)
        {
           idd.t   <-  seq(1,n.t)[dates.lst == days.utc[i.d]]
           n.d.t   <-  length(idd.t)
           if (n.d.t == 0) #not contain certain days
            next

           idx.day[idd.t]   <-  i.d
           sun.rs  <- GetSunRiseSet(as.Date(days.utc[i.d]),lon,lat)
           if (class(sun.rs[1]) == 'PCICt') #not polar day or night
           {
               id.night<-  idd.t[GetTimeNum(time.lst[idd.t]) <=
                                 GetTimeNum(sun.rs[1]) |
                                 GetTimeNum(time.lst[idd.t]) >=
                                 GetTimeNum(sun.rs[2])]
               id.day  <-  idd.t[GetTimeNum(time.lst[idd.t]) >
                                 GetTimeNum(sun.rs[1]) &
                                 GetTimeNum(time.lst[idd.t]) <
                                 GetTimeNum(sun.rs[2])]
               idx.day[id.day]     <-  i.d
               idx.day[id.night]   <-  -i.d
               f.sfc[i.d]  <-  sum(f.sfc.ts[id.day])
               m.sfc[i.d]  <-  sum(m.sfc.ts[id.day])
               f.atm[i.d]  <-  sum(cp*(t2m[id.day+1] - t2m[id.day])[1:(length(id.day)-1)]) - f.sfc[i.d]
               m.atm[i.d]  <-  sum(Lv*(q2m[id.day+1] - q2m[id.day])[1:(length(id.day)-1)]) - m.sfc[i.d]
           } else
           {
               idx.day[idd.t]  <-  sun.rs*i.d
           }
           idx.day[is.na(idx.day)]  <-  0
        }
        #determine the daytime
    } else if (class(mean(times)) == "numeric")
    #the input hours is in numeric format (e.g., 1,2,3,...)
    {
        stop('mixing_diagram: this format is not supported currently.')
        if (is.null(day.lst))
        warning("mixing_diagram: no specific daytime information. The sun rise and sun set time are set to 6am and 6pm LST.")
        sun.rs  <-  c(6,18)
    } else
        stop('mixing_diagram: Time format is unacceptable.')

    #get the index of daytime starts and ends
    mdiag.o <-  list(cpt=cp.t,lvq=lv.q,f.sfc=f.sfc,
                     m.sfc=m.sfc,f.atm=f.atm,sp=sp,
                     m.atm=m.atm,times=time.lst,idx.day=idx.day)
    class(mdiag.o)  <-  'hmdiag1d'
    return(mdiag.o)
}

#' @export
mixing_diagram_simple   <-  function(sp,t2m,q2m,shf,lhf,hpbl,dt){
    #reference pressure [Pa]
    Lv  <-  2.5e6
    cp  <-  1005.7
    Rd  <-  287.04
    ep  <-  .622

    n.t <-  length(t2m.i)
    rho <-  sp.i/(Rd*t2m.i*(1 + (q2m.i/ep))/(1 + q2m.i))

    cp.t    <-  cp*t2m.i
    lv.q    <-  Lv*q2m.i

    F.sfc   <-  (shf.i*dt)/(rho*pbl.i)
    M.sfc   <-  (lhf.i*dt)/(rho*pbl.i)

    F.atm   <-  c(cp*(t2m.i[2:n.t] - t2m.i[1:(n.t-1)]),0) - F.sfc
    M.atm   <-  c(Lv*(q2m.i[2:n.t] - q2m.i[1:(n.t-1)]),0) - M.sfc

    return(list(lq=lv.q,ct=cp.t,F.sfc=sum(F.sfc),
                M.sfc=sum(M.sfc),F.atm=sum(F.atm),
                M.atm=sum(M.atm)))
}

#' @export
plot_mdiag_1d   <-  function(mdiag,id.day,day.only=T,rh=F,FM=F,...){
    Lv  <-  2.5e6
    cp  <-  1005.7
    if (class(id.day) == 'PCICt' |
        class(id.day) == 'Date')
    {
        dates   = unique(format(mdiag$times,'%Y-%m-%d'))
        idd.day = unique(abs(mdiag$idx.day)[mdiag$idx.day != 0])
        id.date = format(id.day,'%Y-%m-%d')
        id.day  = idd.day[dates == id.date]
    } else if (!is.numeric(id.day))
        stop('plot_diag_1d: the id.day cannot be recognized.')

    n.t <-  length(mdiag$cpt)
    id.t    <-  seq(1,n.t)[abs(mdiag$idx.day) == id.day]
    id.dn   <-  mdiag$idx.day[id.t]
    cp.t    <-  mdiag$cpt[id.t]
    lv.t    <-  mdiag$lvq[id.t]
    sp.t    <-  mdiag$sp[id.t]
    f.sfc   <-  mdiag$f.sfc[id.day]
    m.sfc   <-  mdiag$m.sfc[id.day]
    f.atm   <-  mdiag$f.atm[id.day]
    m.atm   <-  mdiag$m.atm[id.day]
    hrs.t   <-  format(mdiag$times[id.t],'%H:%M')
    sp.m    <-  mean(sp.t)

    col.night   <-  1
    col.day <-  2

    if (day.only)
    {
        if (length(lv.t[id.dn > 0]) <= 1)
            stop('plot_mdiag_1d: the daytime of this day is too short (step <= 1).')
        if (FM)
        {
            x.ra    <-  range(lv.t[id.dn > 0],
                              (lv.t[id.dn > 0])[1] + m.sfc)
            y.ra    <-  range(cp.t[id.dn > 0],
                              (cp.t[id.dn > 0])[1] + f.sfc)
        } else
        {
            x.ra    <-  range(lv.t[id.dn > 0])
            y.ra    <-  range(cp.t[id.dn > 0])
        }
        plot(lv.t[id.dn >0],cp.t[id.dn > 0],type='l',
             ylab=expression('2-m Temperature ('*L[v]*q*' [J/kg])'),
             xlab=expression('2-m Humidity ('*C[p]*K*'[J/kg])'),
             xlim=x.ra,ylim=y.ra,col=col.night)
        text((lv.t[id.dn > 0])[1],(cp.t[id.dn > 0])[1],
             labels=(hrs.t[id.dn > 0])[1])
        len.dnd <-  length(cp.t[id.dn > 0])
        text((lv.t[id.dn > 0])[len.dnd],(cp.t[id.dn > 0])[len.dnd],
             labels=(hrs.t[id.dn > 0])[len.dnd])
    } else
    {
        if (FM)
        {
            x.ra    <-  range(lv.t,
                              (lv.t[id.dn > 0])[1] + m.sfc)
            y.ra    <-  range(cp.t,
                              (cp.t[id.dn > 0])[1] + f.sfc)
        } else
        {
            x.ra    <-  range(lv.t)
            y.ra    <-  range(cp.t)
        }
        plot(lv.t,cp.t,type='l',
             ylab=expression('2-m Temperature ('*L[v]*q*' [J/kg])'),
             xlab=expression('2-m Humidity ('*C[p]*K*'[J/kg])'),
             xlim=x.ra,ylim=y.ra,col=col.night)
        if (length(lv.t[id.dn > 1]))
        {
            lines(lv.t[id.dn > 0], cp.t[id.dn > 0],col=col.day)
            text((lv.t[id.dn > 0])[1],(cp.t[id.dn > 0])[1],
                 labels=(hrs.t[id.dn > 0])[1])
            len.dnd <-  length(cp.t[id.dn > 0])
            text((lv.t[id.dn > 0])[len.dnd],(cp.t[id.dn > 0])[len.dnd],
                 labels=(hrs.t[id.dn > 0])[len.dnd])
        } else
        {
            text(lv.t[1],cp.t[1],
                 labels=hrs.t[1])
            len.dnd <-  length(cp.t)
            text(lv.t[len.dnd],cp.t[len.dnd],
                 labels=hrs.t[len.dnd])
        }
    }

    if (FM)
    {
        #start point
        p.s.x <-  (lv.t[id.dn >0])[1]
        p.s.y <-  (cp.t[id.dn >0])[1]
        arrows(x0=p.s.x,y0=p.s.y,
               x1=p.s.x+m.sfc,
               y1=p.s.y+f.sfc,
               col=rgb(0,0,1,.5),lwd=3,
               angle=15,length=.2)
        arrows(x0=p.s.x+m.sfc,
               y0=p.s.y+f.sfc,
               x1=p.s.x+m.sfc+m.atm,
               y1=p.s.y+f.sfc+f.atm,
               col=rgb(0,1,0,.5),lwd=3,
               angle=15,length=.2)
    }

    if (rh)
    {
        xx  <-  seq(x.ra[1],x.ra[2],length.out=100)
        yy  <-  seq(y.ra[1],y.ra[2],length.out=100)
        mat <-  outer(xx,yy,
                      Vectorize(function(x,y) SH2RH(x/Lv,y/cp,sp.m)))
        contour(xx,yy,mat,add=T,levels=seq(0,100,10),col='grey',
                labels=paste0(seq(0,100,10),'%'))
    }
}

t_contour  <-  function(lv.i,rh.i,p.i){
    cp  <-  1005.7
    Lv  <-  2.5e6
    L   <-  2.5e6
    Rw  <-  461.52
    T0  <-  273.15
    Es.T0   <-  6.11
    q2m.t   <-  lv.i/Lv

    med.t   <-  q2m.t*(p.i)/(rh.i*(.622+.378*q2m.t))

    t.o <-  1/(1/T0 - log(med.t/Es.T0)*Rw/L)

    return(t.o*cp)
}
