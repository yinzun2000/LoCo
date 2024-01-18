#' Calculate the mixing diagram to illustrate the co-evolution of potential temperature and humidity during a selected time period.

#' @param sp a numeric vector of surface pressure [Pa]
#' @param t2m a numeric vector of 2m temperature [K]
#' @param q2m a numeric vector of 2m specific humidity [kg.kg-1]
#' @param shf a numeric vector  of upward sensible heat flux [W.m-2]
#' @param lhf a numeric vector of upward latent heat flux [W.m-2]
#' @param hbpl a numeric vector of planetary boundary layer height [m]
#' @param dt time interval [seconds]
#' @return a list containing the evolution of humidity and temperature at a entire day circle.
#' @export
mixing_diagram_simple   <-  function(sp,t2m,q2m,shf,lhf,hpbl,dt){
    #reference pressure [Pa]
    Lv  <-  2.5e6
    cp  <-  1005.7
    Rd  <-  287.04
    ep  <-  .622

    n.t <-  length(t2m)
    rho <-  sp/(Rd*t2m*(1 + (q2m/ep))/(1 + q2m))

    cp.t    <-  cp*t2m
    lv.q    <-  Lv*q2m

    F.sfc   <-  (shf*dt)/(rho*hpbl)
    M.sfc   <-  (lhf*dt)/(rho*hpbl)

    F.atm   <-  c(cp*(t2m[2:n.t] - t2m[1:(n.t-1)]),0) - F.sfc
    M.atm   <-  c(Lv*(q2m[2:n.t] - q2m[1:(n.t-1)]),0) - M.sfc

    return(list(lq=lv.q,ct=cp.t,F.sfc=sum(F.sfc),
                M.sfc=sum(M.sfc),F.atm=sum(F.atm),
                M.atm=sum(M.atm)))
}
