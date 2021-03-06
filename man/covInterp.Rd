% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.R
\name{covInterp}
\alias{covInterp}
\title{Covariance matrix interpolation}
\usage{
covInterp(tt, phi.fun, omega, nugg.fun, t.grid)
}
\arguments{
\item{tt}{time vector to interpolate on}

\item{phi.fun}{eigen function}

\item{omega}{eigen value}

\item{nugg.fun}{nugget function}

\item{t.grid}{time vector on which the covariance matrix was calculated}
}
\value{
covariance matrix interpretated on tt
}
\description{
Covariance matrix interpolation
}
