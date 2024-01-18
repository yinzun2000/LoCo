#' Plot one-day mixing diagram.

#' @param mdiag the input object as an output of mixing_diagram
#' @param id.day It can be in integer format to define the number of the day for
#' plot. Or it can be in PCICt or Date format to define the specific date for
#' plot.
#' @param day.only If only plot the daytime period (default) or the entire
#' day. Note that day.only=T may lead to error if there is a polar night.
#' @param rh Whether relative humidity is plotted or not (default) as contour
#' lines in the background.
#' @param FM Whether the vectors of surface and atmosphere are plotted or not
#' (default).
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
