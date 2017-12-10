#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "dzeta_function.h"

SEXP DZetaFunction(SEXP qsq_, SEXP n_, SEXP l_, SEXP m_, SEXP dvec_, SEXP gamma_, SEXP A_, SEXP tol_, SEXP Lmax_ ) {
/*  int dzeta_function (double z[2], double q2, int l, int m, int*dvec, double gamma, double A, double epsAbs, double epsRel, int Lmax); */
  

  PROTECT( qsq_   = AS_NUMERIC(qsq_)   );
  PROTECT( l_     = AS_INTEGER(l_)     );
  PROTECT( m_     = AS_INTEGER(m_)     );
  PROTECT( gamma_ = AS_NUMERIC(gamma_) );
  PROTECT( A_     = AS_NUMERIC(A_)     );
  PROTECT( dvec_  = AS_INTEGER(dvec_)  );
  PROTECT( n_     = AS_INTEGER(n_)     );
  PROTECT( tol_   = AS_NUMERIC(tol_)   );
  PROTECT( Lmax_  = AS_INTEGER(Lmax_)  );

  double       *qsq    = NUMERIC_POINTER(qsq_);
  const int    l       = INTEGER_POINTER(l_)[0];
  const int    m       = INTEGER_POINTER(m_)[0];
  double       *gamma  = NUMERIC_POINTER(gamma_);
  double       *A      = NUMERIC_POINTER(A_);
  const int    n       = INTEGER_POINTER(n_)[0];
  const double tol     = NUMERIC_POINTER(tol_)[0];
  int          *dvec   = INTEGER_POINTER(dvec_);
  const int    Lmax    = INTEGER_POINTER(Lmax_)[0];
  SEXP res;
  Rcomplex * resp;
  double ires[2];
  PROTECT(res = NEW_COMPLEX(n));
  resp = COMPLEX_POINTER(res);

  for(int i = 0; i < n; i ++) {
    int exitstatus = 0;

    /* check for NaN or NA in qsq */
    if(ISNAN(qsq[i]) || ISNA(qsq[i])) {
      resp[i].r = NA_REAL;
      resp[i].i = NA_REAL;
    }
    else {
      exitstatus = dzeta_function (ires, qsq[i], l, m, dvec, gamma[i], A[i], tol, 1.e12, Lmax);
      resp[i].r = ires[0];
      resp[i].i = ires[1];
    }

    /* check for anything went wrong in dzeta_function evaluation */
    if( exitstatus != 0 ) {
      warning("[QZetaFunction] Error from dzeta_function for momentum number %d, status was %d\n", i, exitstatus);
      resp[i].r = NA_REAL;
      resp[i].i = NA_REAL;
    }
  }

  UNPROTECT(10);
  return(res);

}
