#' Calculate the low-level humidity index (Hilow) as a part of the CTP-Hilow framework.

#' @param p.prof 1D array of pressure profile [Pa].
#' @param t.prof 1D array of temperature profile [Pa].
#' @param q.prof 1D array of specific humidity profile [kg/kg].
#' @param ps value of surface pressure [Pa] (default=1e5 Pa).
#' @param nseg number of segments for integration (default=20).
#' @return a list containing ctp [J/kg] and hilow [C].

#' @export
CTP_Hilow <-  function(p.prof,t.prof,q.prof,ps=1e5,nseg=20){
  #parameters
  grav  <-  9.81
  Rd  <-  287.04
  cp  <-  1005.7
  Rcp <-  Rd/cp
  Lv  <-  2.5e6
  Lv_cp <-  Lv/cp
  Rv  <-  461.5
  ep  <-  0.622

  #check if profile is upside down
  if (p.prof[1] < p.prof[length(p.prof)])
  {
      p.prof  <-  rev(p.prof)
      t.prof  <-  rev(t.prof)
      q.prof  <-  rev(q.prof)
  }

  #1. Calculate CTP
  #get the increment of pressure
  d.p <-  (3e4 - 1e4)/nseg

  #initialize the variables
  ctp.o <-  0
  t.par.old <-  LocLin(ps-1e4,p.prof,t.prof)

  for (i in 1:nseg)
  {
    #the pressure of the bottom and top
    p.bot <-  ps -1e4 - d.p*(i-1)
    p.top <-  ps -1e4 - d.p*i

    #find the corresponding t and q [K]
    t.bot <-  LocLin(p.bot,p.prof,t.prof)
    t.top <-  LocLin(p.top,p.prof,t.prof)
    #[Pa]
    q.bot <-  LocLin(p.bot,p.prof,q.prof)
    q.top <-  LocLin(p.top,p.prof,q.prof)

    #----------------------------------------------------
    #--- Get moist adiabatic lapse rate [K/m] and
    #--- the depth of the layer from lower to upper level
    #----------------------------------------------------
    p.mid <-  PrsIntPol(p.bot,p.top,p.bot,p.top)
    t.mid <-  PrsIntPol(t.bot,t.top,p.bot,p.top)
    q.mid <-  PrsIntPol(q.bot,q.top,p.bot,p.top)
    dz    <-  (p.bot - p.top)/
              (grav * p.mid/(Rd*t.mid*(1+(q.mid/ep))/
                             (1+q.mid)))
    #calculate saturated specific humidity [kg.kg-1]
    q.sat  <-  SSH(p.mid,t.mid)
    moist.lapse <-  (grav/cp) *
                    (1 + (Lv * q.sat)/(Rd*t.mid)) /
                    (1 + (Lv**2 * q.sat)/(cp*Rv*t.mid**2))

    #----------------------------------------------------
    #--- Get the parcel temperature
    #----------------------------------------------------
    t.par <-  t.par.old - moist.lapse*dz
    
    #----------------------------------------------------
    #--- Get mid-point temps from environment and parcel
    #----------------------------------------------------
                                                                          t.par.mid <-  0.5 * (t.par + t.par.old)
                                                                          t.mid     <-  0.5 * (t.bot + t.top)
    
    #----------------------------------------------------
    #--- Integrate from old increment to increment level
    #----------------------------------------------------
    ctp.o <-  ctp.o + (Rd*(t.par.mid - t.mid))*log(p.bot/p.top)

    #----------------------------------------------------
    #--- Update the last increment values
    #----------------------------------------------------
    t.par.old <-  t.par
    t.bot     <-  t.top
    q.bot     <-  q.top
    p.bot     <-  p.top
  }

  #2. calculate Hilow
  #calculate t50, t150, q50, q150
  t50   <-  LocLin(ps-5e3,p.prof,t.prof)
  t150  <-  LocLin(ps-15e3,p.prof,t.prof)
  q50   <-  LocLin(ps-5e3,p.prof,q.prof)
  q150  <-  LocLin(ps-15e3,p.prof,q.prof)

  #calculate the dew point temperature (K)
  t.dew50   <-  DewT(q50,ps-5e3)
  t.dew150  <-  DewT(q150,ps-15e3)

  hilow.o <-  (t50 - t.dew50) + (t150 - t.dew150)
  return(list(ctp=ctp.o,hilow=hilow.o))

}

LocLin  <-  function(r.i,r.a,v.a){
  #get linear interpolate according to reference value and arrays
  #r.i: input reference
  #r.a: reference arrays
  #v.a: value arrays
  if ((r.i > max(r.a)) | (r.i < min(r.a)))
      stop("CalCTP: reference input is out of the reference array")

  #find the index besides the r.i
  ind.a <-  which.min(r.a[r.a >= r.i])
  ind.b <-  which.max(r.a[r.a < r.i])

  if (r.a[ind.a] == r.i)
  {
      v.o <-  v.a[ind.a]
  } else
      v.o <-  v.a[ind.a] + (v.a[ind.b] - v.a[ind.a])*
          (r.i - r.a[ind.a])/(r.a[ind.b] - r.a[ind.a])
  return(v.o)
}

DewT  <-  function(q.i,p.i){
  #dew point temperature [K]
  #q.i: specific humidity [kg.kg-1]
  #p.i: pressure [Pa]

  #some parameters
  a.t <-  610.8
  b.t <-  237.3
  c.t <-  17.2693882

  #--------------------------------------------
  #--- Vapor pressure and convert to Pa
  #--------------------------------------------
  e.t <-  q.i*p.i/(0.622+0.378*q.i)

  #--------------------------------------------
  #--- Vapor pressure and convert to Pa
  #--------------------------------------------
  v.o <-  (log(e.t/a.t)*b.t) / (c.t-log(e.t/a.t)) + 273.15

  return(v.o)
}


