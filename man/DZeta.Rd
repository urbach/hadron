\name{DZeta}
\alias{DZeta}
\title{Computes Lueschers Zeta Function}
\description{
}
\usage{
  DZeta(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1., A=1.,
tol=0.000001, Lmax = 5)
}
\arguments{
  \item{qsq}{
    a vector containing the squared scattering momenta to evaluate the
    Zeta function on
  }
  \item{l}{
    the angular momentum quantum number l
  }
  \item{m}{
    the angular momentum quantum number m
  }
  \item{dvec}{
    integer vector with three elements indicating the total momentum
  }
  \item{gamma}{
    boost factor
  }
  \item{A}{
    mass difference parameter for two particles of different mass
  }
  \item{tol}{
    error tolerance in numerical integrations
  }
  \item{Lmax}{
    summation limit box size
  }
}
\value{
  a vector of complex values of the same length as \code{qsq} containing
  the value of the Zeta function
}
\details{
  The Zeta function is evaluated following 
  Phys. Rev., D85:114507, 2012. The
  evaluation is implemented in C using GSL routines.
}
\references{
  See Phys. Rev., D85:114507, 2012
}
\examples{
  to be added
}
\author{Marcus Petschlies, \email{marcus.petschlies@hiskp.uni-bonn.de}}
