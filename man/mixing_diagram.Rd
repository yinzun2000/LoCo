% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixing_diagram.R
\name{mixing_diagram}
\alias{mixing_diagram}
\title{Calculate the mixing diagram to illustrate the co-evolution of potential temperature and humidity during a entire day process.}
\usage{
mixing_diagram(
  sp,
  t2m,
  q2m,
  shf,
  lhf,
  hpbl,
  times,
  lon = NULL,
  lat = NULL,
  time.type = "UTC",
  day.lst = NULL
)
}
\arguments{
\item{sp}{a numeric vector of surface pressure \link{Pa}}

\item{t2m}{a numeric vector of 2m temperature \link{K}}

\item{q2m}{a numeric vector of 2m specific humidity \link{kg.kg-1}}

\item{shf}{a numeric vector  of upward sensible heat flux \link{W.m-2}}

\item{lhf}{a numeric vector of upward latent heat flux \link{W.m-2}}

\item{lon}{if the time is UTC}

\item{time.type}{whether the time is LST or UTC (default).}

\item{hbpl}{a numeric vector of planetary boundary layer height \link{m}}

\item{time}{information. It can be either in PCICt format (e.g., "1980-01-01 10:50:33") or in numeric format (hours, e.g., 12, 12.5, etc). If the numeric format is applied, there will no daytime information generated.}
}
\value{
a list containing the evolution of humidity and temperature at a entire day circle.
}
\description{
Calculate the mixing diagram to illustrate the co-evolution of potential temperature and humidity during a entire day process.
}
