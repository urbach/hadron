## interface function for the C implementation of Lueschers Zeta function, using quadmath library

DZeta <- function(qsq, l=0, m=0, dvec=c(0,0,0), gamma=1, A=1, tol=0.000001, Lmax=6) {
  if(!is.numeric(qsq) || !is.numeric(l) || !is.numeric(m) || !is.numeric(dvec) || !is.numeric(gamma) || !is.numeric(A) || !is.numeric(tol))
    stop("argument qsq to LuescherZeta must be numeric!\n")
  n <- length(qsq)
  if(length(gamma) == 1) {
    gamma <- rep(gamma, times=n)
  }
  if(n != length(gamma)) {
    stop("gamma must be either a scalar or a vector of length identical with the length of qsq\n")
  }
  if(length(A) == 1) {
    A <- rep(A, times=n)
  }
  if(n != length(A)) {
    stop("A must be either a scalar or a vector of length identical with the length of qsq\n")
  }

  return(.Call( "DZetaFunction",    as.double(qsq), as.integer(n), as.integer(l[1]), as.integer(m[1]), as.integer(dvec[1:3]),  as.double(gamma), as.double(A), as.double(tol), as.integer(Lmax) ) )
}
