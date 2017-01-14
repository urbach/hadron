#include <quadmath.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_integration.h>
#include "init_zeta_function.h"
#include "qzeta_function.h"
/*******************************************************************************
 * dvec = integer vector, P = 2 pi d / L
 *
 * P^m_l(cos(theta)) 
 *  double gsl_sf_legendre_sphPlm (int l, int m, double x) 
 *
 * check convergence with epsAbs and epsRel
 *
 * typedef for integration function double (*intfunction)    (double x, void*params)
 * __float128 to string: int quadmath_snprintf (char *s, size_t size, const char *format, ...) 
 * or try printf("%Qe", 1.2Q)
 *
 * We characterize the n-vectors into 8 clases; we use 0 < h < k < l
 *     repres.  0-type <=-type
 * (0) (0,0,0)  3      0
 * ======================
 * (1) (0,0,h)  2      1
 * ======================
 * (2) (0,h,h)  1      2
 * ----------------------
 * (3) (0,h,k)  1      3
 * ======================
 * (4) (h,h,h)  0      0
 * ----------------------
 * (5) (h,h,k)  0      1
 * ----------------------
 * (6) (h,k,k)  0      2
 * ----------------------
 * (7) (h,k,l)  0      3
 *
 *  normalized associated Legendre polynomicals with negative m
 *
 *  sqrt{(2l+1)/(4 pi)}  x sqrt{(l-m)/(l+m) } P_l^m =
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l+|m|)/(l-|m|) } P_l^{-|m|} =
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l+|m|)/(l-|m|) } x (-1)^m  (l-|m|)/(l+|m|)  P_l^{|m|}
 *    sqrt{(2l+1)/(4 pi)}  x sqrt{(l-|m|)/(l+|m|) } x (-1)^m  P_l^{|m|}
 *    (-1)^m x associated Legendre polynomial with (l, |m|)
 *******************************************************************************/

__float128 qplm( unsigned int l, int m, __float128 x);

#ifdef QZETA_USE_KAHAN

#define _QKAHAN_SUM_CO_PL_CO(_a, _b, _c) { \
  __float128 _ks_y, _ks_t;                 \
  _ks_y  = (_c)[0] - (_b)[0];              \
  _ks_t  = (_a)[0] + _ks_y;                \
  (_b)[0] = (_ks_t - (_a)[0] ) - _ks_y;    \
  (_a)[0] = _ks_t;                         \
  _ks_y  = (_c)[1] - (_b)[1];              \
  _ks_t  = (_a)[1] + _ks_y;                \
  (_b)[1] = (_ks_t - (_a)[1] ) - _ks_y;    \
  (_a)[1] = _ks_t;                         \
}
    
#define _QKAHAN_SUM_CO_PL_RE(_a, _b, _r) { \
  __float128 _ks_y, _ks_t;                 \
  _ks_y  = (_r) - (_b)[0];                 \
  _ks_t  = (_a)[0]  + _ks_y;               \
  (_b)[0] = (_ks_t - (_a)[0] ) - _ks_y;    \
  (_a)[0] = _ks_t;                         \
}

#else

#define _QKAHAN_SUM_CO_PL_CO(_a, _b, _c) { \
  (_a)[0] += (_c)[0];\
  (_a)[1] += (_c)[1];\
}

#define _QKAHAN_SUM_CO_PL_RE(_a, _b, _r) { \
  (_a)[0] += (_r);\
}

#endif  /* of if def QZETA_USE_KAHAN */


/* #define _V3_EQ_V3(_a, _b) ( ( (_a)[0]==(_b)[0] ) && ( (_a)[1]==(_b)[1] ) && ( (_a)[2]==(_b)[2] ) ) */

#define QZETA_EPSILON 1.e-14

#define _QSQR(_x) ((_x)[0]*(_x)[0]+ (_x)[1]*(_x)[1]+ (_x)[2]*(_x)[2])

#define _QSCP(_x, _y) ((_x)[0]*(_y)[0] + (_x)[1]*(_y)[1] + (_x)[2]*(_y)[2])

#define _QCO_EQ_CO_TI_CO(_c, _a, _b) { \
  (_c)[0] = (_a)[0] * (_b)[0] - (_a)[1] * (_b)[1];\
  (_c)[1] = (_a)[0] * (_b)[1] + (_a)[1] * (_b)[0];\
}

#define _QCO_EQ_CO_TI_RE(_c, _a, _r) {\
  (_c)[0] = (_a)[0] * (_r);           \
  (_c)[1] = (_a)[1] * (_r);           \
}

#define _QCO_PL_EQ_CO_TI_RE(_c, _a, _r) {\
  (_c)[0] += (_a)[0] * (_r);             \
  (_c)[1] += (_a)[1] * (_r);             \
}

#define _QCO_TI_EQ_RE(_c, _r) { \
  (_c)[0] *= (_r);              \
  (_c)[1] *= (_r);              \
}

#define _QCO_EQ_ZERO(_c)  {\
  (_c)[0] = 0.0Q;          \
  (_c)[1] = 0.0Q;          \
}

#define _QCO_EQ_CO(_a, _b) {\
  (_a)[0] = (_b)[0];        \
  (_a)[1] = (_b)[1];        \
}


/***********************************************************
 * angular momentum functions
 ***********************************************************/
#define _PLM_L_MAX 4
static __float128 plm_norm[(_PLM_L_MAX+1)*(_PLM_L_MAX+1)];

int qinit_plm_norm (void) {

  int k, l, m;
  int idxp, idxm, idx0;
  __float128 qone_over_sqrt4pi = 0.25Q * M_2_SQRTPIq;
  __float128 z;

  for(l=0; l <= _PLM_L_MAX; l++) {

    idx0 = l * (l + 1);
    plm_norm[idx0] = sqrtq(2.Q * l + 1.Q) * qone_over_sqrt4pi;

    for(m = 1; m<=l; m++) {

      idxm = idx0 - m;
      idxp = idx0 + m;

      plm_norm[idxp] = plm_norm[idx0];

      z = (__float128)(l+m);
      for(k=l+m-1; k>l-m; k--) {
        z *= (__float128)k;
      }
      plm_norm[idxp] = plm_norm[idx0] / sqrtq( z );

      plm_norm[idxm] =  (__float128)(1 - 2*(m%2)) * plm_norm[idxp];
    }
  }

  /* TEST */
/*
  for(l=0; l <= _PLM_L_MAX; l++) {
    for(m = -l; m <= l; m++) {
      idx0 = l*l + (l+m);
      fprintf(stdout, "# [qinit_plm_norm] l=%3d m=%3d idx=%3d norm=%25.16Qe\n", l, m, idx0, plm_norm[idx0]);
    }
    fprintf(stdout, "# [qinit_plm_norm] ----------------------------------------------------------------------------\n");
  }
*/
  /* END OF TEST */

  return(0);
}

/***************************************************************************
 * normalized associated Legendre polynomials
 * 
 * - input: l and m; l >= 0 and -l <= m <= l
 *          x with  -1 <= x <= 1; x = cos(theta)
 * - output: sqrt{(2l+1)/(4 pi)} x sqrt{(l-m)/(l+m) } P_l^m (x)
 * - explicit form of P_l^m from
 *     https://en.wikipedia.org/wiki/Associated_Legendre_polynomials
 ***************************************************************************/
__float128 qplm( unsigned int l, int m, __float128 x) {

  int idx = l * ( l + 1 ) + m;
  int mabs = m<0 ? -m : m;

  __float128 norm = plm_norm[idx];
  __float128 p = sqrtq(-1.Q);
  __float128 xx;

  if(l == 0) {
    p = 1.Q;
  } else if(l == 1) {
    if(mabs == 0) {
      p = x;
    } else if (mabs == 1) {
      p = -sqrtq( 1.Q - x*x );
    }
  } else if(l == 2) {
    if(mabs == 0) {
      p = 0.5Q * (3.Q*x*x - 1. );
    } else if (mabs == 1) {
      p = -3.Q * x * sqrtq( 1.Q - x*x );
    } else if (mabs == 2) {
      p = 3.Q * ( 1.Q - x*x );
    }
  } else if(l == 3) {
    if(mabs == 0) {
      p = 0.5Q * x * ( 5.Q*x*x - 3.Q );
    } else if(mabs == 1) {
      p = -1.5Q * (5.Q*x*x - 1.Q) *  sqrtq( 1.Q - x*x );
    } else if(mabs == 2) {
      p = 15.Q * x * (1.Q - x*x );
    } else if(mabs == 3) {
      xx = 1.Q - x*x;
      p = -15.Q * xx * sqrtq(xx);
    }
  } else if(l == 4) {
    if(mabs == 0) {
      xx = x*x;
      p = 0.125Q * ( ( 35.Q*xx -30.Q ) * xx +3.Q );
    } else if (mabs == 1) {
      xx = x*x;
      p = -2.5Q * x * (7.Q * xx - 3.Q) * sqrtq(1.Q - xx);
    }
    } else if (mabs == 2) {
      xx = x*x;
      p = 7.5Q * (7.Q * xx - 1.Q) * (1.Q - xx);
    } else if (mabs == 3) {
      xx = 1.Q - x*x;
      p = -105.Q * x*  xx * sqrtq(xx);
    } else if (mabs == 3) {
      xx = 1.Q - x*x;
      p = 105.Q * xx * xx;
  }
  return(p * norm);
}

/**************************************************************************
 * rotations and reflections
 **************************************************************************/

#if 0
int init_refrot_list (void) {
  int k, i, iclass, have_rotref, nrotref;
  int n0, n1, n2;
  int s0, s1, s2;
  int p[3], q[3];
  int etype, ztype;

  int vector_types[8][3] = {
    {0,0,0},
    {0,0,1},
    {0,1,1},
    {0,1,2},
    {1,1,1},
    {1,1,2},
    {1,2,2},
    {1,2,3} };

  k = 0;
  for(n0=0; n0<3; n0++) {
    /* cyclic permutation */
    n1 = ( n0 + 1) % 3;
    n2 = ( n0 + 2) % 3;
    /* all sign changes */
    for( s0=1; s0>=-1; s0-=2) {
    for( s1=1; s1>=-1; s1-=2) {
    for( s2=1; s2>=-1; s2-=2) {

      rotref_list[k][0][0] = n0;
      rotref_list[k][0][1] = n1;
      rotref_list[k][0][2] = n2;

      rotref_list[k][1][0] = s0;
      rotref_list[k][1][1] = s1;
      rotref_list[k][1][2] = s2;
      
      k++;
    }}}

    /* anti-cyclic permutation */
    n1 = ( n0 + 2) % 3;
    n2 = ( n0 + 1) % 3;
    /* all sign changes */
    for( s0=1; s0>=-1; s0-=2) {
    for( s1=1; s1>=-1; s1-=2) {
    for( s2=1; s2>=-1; s2-=2) {

      rotref_list[k][0][0] = n0;
      rotref_list[k][0][1] = n1;
      rotref_list[k][0][2] = n2;

      rotref_list[k][1][0] = s0;
      rotref_list[k][1][1] = s1;
      rotref_list[k][1][2] = s2;
      
      k++;
    }}}
  }  /* end of loop on n0 */

  for(iclass = 0; iclass<8; iclass++) {

    etype = 2*(int)(vector_types[iclass][0] < vector_types[iclass][1]) + (int)(vector_types[iclass][1] < vector_types[iclass][2]);
    ztype = (int)(vector_types[iclass][0] == 0) + (int)(vector_types[iclass][1] == 0) + (int)(vector_types[iclass][2] == 0);

    /* fprintf(stdout, "# [] vector %3d%3d%3d e-type %d z-type %d\n", vector_types[iclass][0], vector_types[iclass][1], vector_types[iclass][2], etype, ztype); */

    nrotref = 0;
    for(k=0; k<48; k++) {
      p[0] = vector_types[iclass][ rotref_list[k][0][0] ] * rotref_list[k][1][0];
      p[1] = vector_types[iclass][ rotref_list[k][0][1] ] * rotref_list[k][1][1];
      p[2] = vector_types[iclass][ rotref_list[k][0][2] ] * rotref_list[k][1][2];

      /* check against previous rot-refs */
      have_rotref = 0;
      for(i=nrotref-1; i>=0; i--) {
        q[0] = vector_types[iclass][ rotref_selection_permutation[etype][ztype][i][0] ] * rotref_selection_sign[etype][ztype][i][0];
        q[1] = vector_types[iclass][ rotref_selection_permutation[etype][ztype][i][1] ] * rotref_selection_sign[etype][ztype][i][1];
        q[2] = vector_types[iclass][ rotref_selection_permutation[etype][ztype][i][2] ] * rotref_selection_sign[etype][ztype][i][2];
        have_rotref = _V3_EQ_V3(p, q);
        /* fprintf(stdout, "\t%d\t%3d%3d\t%3d%3d%3d\t%3d%3d%3d\t%d\n", iclass, k, i, p[0], p[1], p[2], q[0], q[1], q[2], have_rotref); */
        if( have_rotref ) break;
      }
      /* if we do not have it yet, add the rotref */
      if( !have_rotref ) {
        rotref_selection_permutation[etype][ztype][nrotref][0] = rotref_list[k][0][0];
        rotref_selection_permutation[etype][ztype][nrotref][1] = rotref_list[k][0][1];
        rotref_selection_permutation[etype][ztype][nrotref][2] = rotref_list[k][0][2];

        rotref_selection_sign[etype][ztype][nrotref][0] = rotref_list[k][1][0];
        rotref_selection_sign[etype][ztype][nrotref][1] = rotref_list[k][1][1];
        rotref_selection_sign[etype][ztype][nrotref][2] = rotref_list[k][1][2];

        nrotref++;
        /* fprintf(stdout, "\t added\n"); */
      } /* else {
        fprintf(stdout, "\t discarded\n");
      }
      fprintf(stdout, "# --------------------------------------------------------------------------------------\n"); */

    }  /* end of loop on rot-refs */

    /* set the number of rotrefs */
    rotref_number[etype][ztype] = nrotref;

    /* fprintf(stdout, "# =======================================================================================\n"); */

  }  /* end of loop on vector classes */


  /* TEST */
/*
  for(iclass = 0; iclass<8; iclass++) {

    etype = 2*(int)(vector_types[iclass][0] < vector_types[iclass][1]) + (int)(vector_types[iclass][1] < vector_types[iclass][2]);
    ztype = (int)(vector_types[iclass][0] == 0) + (int)(vector_types[iclass][1] == 0) + (int)(vector_types[iclass][2] == 0);

    fprintf(stdout, "# [init_refrot_list] vector %3d%3d%3d e-type %d z-type %d rotref %d \n", vector_types[iclass][0], vector_types[iclass][1], vector_types[iclass][2], etype, ztype, rotref_number[etype][ztype]);

    for(k=0; k<rotref_number[etype][ztype]; k++) {
      fprintf(stdout, "\t%2d\t%3d%3d%3d\t%3d%3d%3d\n", k, 
          rotref_selection_permutation[etype][ztype][k][0],
          rotref_selection_permutation[etype][ztype][k][1],
          rotref_selection_permutation[etype][ztype][k][2],
          rotref_selection_sign[etype][ztype][k][0],
          rotref_selection_sign[etype][ztype][k][1],
          rotref_selection_sign[etype][ztype][k][2]);
    }
    fprintf(stdout, "# [init_refrot_list]\n# [init_refrot_list]\n");
  }  
*/
  /* END OF TEST */


  return(0);
}  /* end of init_refrot_list */
#endif

double qintegrand_12 (double t, void*params) {
  __float128 qt   = (__float128)t;
  __float128 q2   = (__float128)((double*)params)[0];
  __float128 a    = (__float128)((double*)params)[1];
  return( (double)(expq( q2 * qt - a / qt) *sqrtq(qt)) ); 
}  /* end of qintegrand_12 */

double qintegrand_32 (double t, void*params) {
  __float128 qt   = (__float128)t;
  __float128 q2   = (__float128)((double*)params)[0];
  __float128 a    = (__float128)((double*)params)[1];
  return( (double)(expq( q2 * qt - a / qt) *sqrtq(qt) * qt) ); 
}  /* end of qintegrand_32 */

double qintegrand_lp32 (double t, void*params) {
  __float128 qt   = (__float128)t;
  __float128 q2   = (__float128)((double*)params)[0];
  __float128 a    = (__float128)((double*)params)[1];
  __float128 l    = (__float128)((double*)params)[2];
  return( (double)(expq( q2 * qt - a / qt)  / ( sqrtq(qt) * powq(qt,l+1) )) ); 
}  /* end of qintegrand_lp32 */


double qint_l0m0_kernelFunction (double t, void*params) {
  __float128 qt   = (__float128)t;
  __float128 q2   = (__float128)((double*)params)[0];
  return( (double)(expq(qt*q2) * sqrtq( qt )) );
}  /* end of qint_l0m0_kernelFunction */


double qint_l0m0_kernelFunction2 (double t, void*params) {
  __float128 qt  = (__float128)t;
  __float128 q2  = (__float128)((double*)params)[0];
  return( (double)( ( expq(qt*q2) - 1. ) / ( qt * sqrtq(qt) )) );
}  /* end of qint_l0m0_kernelFunction */

/*****************************************************************************************************************
 * Zeta function
 * input:
 *   q2 - the value for q^2
 *   l, m - angular momentum quantum numbers
 *   dvec - vector d related to total linear momentum P = 2 pi/L x d
 *   gamma - gamma from Lorentz transformation, gamma = 1 / sqrt( 1 - v^2 )
 *   A - numerical factor for case of different masses
 *   epsAbs, epsRel - convergence criteria, not checked here at the moment, just set Lmax large enough from
 *   experience
 *****************************************************************************************************************/

int qzeta_function (double z[2], double q2, int l, int m, int*dvec, double gamma, double A, double epsAbs, double epsRel, int Lmax ) {

  const double int_epsabs = 1.e-08; 
  const double int_epsrel = 1.e-07;
  const size_t int_limit  = 500000;
  const size_t int_key    = 6;

  static int rotref_list[48][2][3];
  static int rotref_selection_permutation[4][4][48][3];
  static int rotref_selection_sign[4][4][48][3];
  static int rotref_number[4][4];

  static int refrot_list_is_initialized = 0;
  static int plm_norm_is_initialized = 0;

  int k1, k2, k3, status, i, irotref;
  int nvec[3];
  __float128 qq2       = (__float128)q2;
  __float128 qgamma    = (__float128)gamma;
  __float128 qgammaInv = 1.0Q / qgamma;
  __float128 qA        = (__float128)A;

  gsl_integration_workspace *int_workspace = NULL;
  gsl_function F;

  int parity_factor =  1 - (l % 2);
  int lessequals_type, zeros_type;
  double int_parameters[3], int_l0m0_parameters;
  double int_integral_l0m0_value, int_integral_l0m0_error;
  double int_integral_12_value, int_integral_12_error, int_integral_32_value, int_integral_32_error;

  __float128 qdvec[3], qnd, qddInv, qr_rrInv, qdd;
  __float128 qterm1[2], qterm2[2], qterm1c[2], qterm2c[2], qterm3[2];
  __float128 qdvecGammaInvMinusOne[3];
  __float128 qint_norm_const[2], qint_norm_var[2];
  __float128 qint_iterate_const=0.0Q, qint_integral_value=0.0Q;
  __float128 qint_l0m0_norm_const, qint_l0m0_add;
  __float128 qshift[3], qn[3];
  __float128 qtmp, qtmp1, qtmp2, qw[2], qw2[2];
  __float128 qr[3], qr_r, qr_rr, qr_costheta, qr_phi;


  if ( refrot_list_is_initialized == 0 ) {
    fprintf(stdout, "# [qzeta_function] initialize refrot_list\n");
    init_refrot_list (rotref_list, rotref_selection_permutation, rotref_selection_sign, rotref_number);
    refrot_list_is_initialized = 1;
  }
  if ( plm_norm_is_initialized == 0 ) {
    fprintf(stdout, "# [qzeta_function] initialize plm_norm\n");
    qinit_plm_norm();
    plm_norm_is_initialized = 1;
  }


  qdvec[0] = (__float128)dvec[0];
  qdvec[1] = (__float128)dvec[1];
  qdvec[2] = (__float128)dvec[2];
  qdd      = _QSQR(qdvec);
  qddInv   = (qdd == 0.0Q) ? 0.0Q : 1.0Q / qdd;

  qshift[0] = -0.5Q * (qgammaInv * qA * qdvec[0]);
  qshift[1] = -0.5Q * (qgammaInv * qA * qdvec[1]);
  qshift[2] = -0.5Q * (qgammaInv * qA * qdvec[2]);

  qdvecGammaInvMinusOne[0] = qdvec[0] * ( qgammaInv - 1.0Q );
  qdvecGammaInvMinusOne[1] = qdvec[1] * ( qgammaInv - 1.0Q );
  qdvecGammaInvMinusOne[2] = qdvec[2] * ( qgammaInv - 1.0Q );


  /* initialize first integration */
  if( (int_workspace = gsl_integration_workspace_alloc ( (size_t)(int_limit+2))) == NULL ) {
    fprintf(stderr, "# [qzeta_function] Error from gsl_integration_workspace_alloc\n");
    return(1);
  }

  qtmp = qgamma * M_PIq * sqrtq(M_PIq);
  i = l % 4;
  if(i == 0) {
    qint_norm_const[0] = qtmp;
    qint_norm_const[1] = 0.0Q;
  } else if (i == 1) {
    qint_norm_const[0] = 0.0Q;
    qint_norm_const[1] = -qtmp;
  } else if (i == 2) {
    qint_norm_const[0] = -qtmp;
    qint_norm_const[1] = 0.0Q;
  } else {
    qint_norm_const[0] = 0.0Q;
    qint_norm_const[1] = qtmp;
  }

  /* fprintf(stdout, "# [qzeta_function] int_norm_const = %Qe + I %Qe\n", qint_norm_const[0], qint_norm_const[1]); */


  /* q2 value */
  int_parameters[0] = q2;
  
  /* initialize second integration */
  qint_l0m0_norm_const = -qgamma * M_PIq * 2.0Q * qq2 * qq2;
  qint_l0m0_add        =  qgamma * M_PIq * ( 2.0Q * qq2 - 1.0Q ) * expq(qq2);
  int_l0m0_parameters  = q2;

  /* initialize real and imaginary part of term 1 and 2 */
  _QCO_EQ_ZERO( qterm1  );
  _QCO_EQ_ZERO( qterm1c );
  _QCO_EQ_ZERO( qterm2  );
  _QCO_EQ_ZERO( qterm2c );
  _QCO_EQ_ZERO( qterm3  );


  for(k3=0; k3 <= Lmax; k3++) {
    nvec[2] = k3;
  for(k2=0; k2 <= k3;    k2++) {
    nvec[1] = k2;
  for(k1=0; k1 <= k2;    k1++) {
    nvec[0] = k1;

    lessequals_type = 2*(int)(k1 < k2) + (int)(k2 < k3);
    zeros_type  = (int)(k1 ==  0) + (int)(k2 ==  0) + (int)(k3 ==  0);

    /* fprintf(stdout, "# [qzeta_function] %3d%3d%3d\t%3d%3d\t%3d\n", k1, k2, k3, lessequals_type, zeros_type, rotref_number[lessequals_type][zeros_type]); */

    for(irotref = 0; irotref < rotref_number[lessequals_type][zeros_type]; irotref++) {

      qn[0] = (__float128)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][0] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][0] );
      qn[1] = (__float128)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][1] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][1] );
      qn[2] = (__float128)( nvec[ rotref_selection_permutation[lessequals_type][zeros_type][irotref][2] ] * rotref_selection_sign[lessequals_type][zeros_type][irotref][2] );


      /* fprintf(stdout, "# [qzeta_function] %3d%3d%3d\t%3d\t%3.0Qf %3.0Qf %3.0Qf\n", k1, k2, k3, irotref, qn[0], qn[1], qn[2]); */

      qnd   = _QSCP(qn, qdvec);

      /***************
       * SECOND TERM *
       ***************/
      qtmp  = qnd * qddInv;

      qr[0] = qn[0] + qtmp * qdvecGammaInvMinusOne[0] + qshift[0];
      qr[1] = qn[1] + qtmp * qdvecGammaInvMinusOne[1] + qshift[1];
      qr[2] = qn[2] + qtmp * qdvecGammaInvMinusOne[2] + qshift[2];

      qr_rr  = _QSQR(qr);
      qr_r   = sqrtq(qr_rr);

      if(qr_r < QZETA_EPSILON) {
        qr_costheta = 1.0Q;
        qr_phi      = 0.0Q;
      } else {
        qr_costheta = qr[2] / qr_r;
        qtmp = qr[0] * qr[0] + qr[1] * qr[1];

        if(qtmp < QZETA_EPSILON) {
          qr_phi = 0.0Q;
        } else {
          qr_phi = atan2q(qr[1], qr[0]);
        }
      }

      /* fprintf(stdout, "# [qzeta_function] %Qe %Qe %Qe \t %Qe %Qe\t %Qe %Qe \n", qr[0], qr[1], qr[2], qr_r, qr_rr, qr_costheta, qr_phi); */

      qtmp  = (__float128)m * qr_phi;
      qtmp1 = powq(qr_r, (__float128)l ) * qplm(l, m, qr_costheta);
      qtmp2 = qr_rr - qq2;
      /* fprintf(stdout, "# [qzeta_function]  %Qe %Qe  %Qe\n", qtmp, qtmp1, qtmp2); */

      qw[0] = cosq( qtmp ) * qtmp1;
      qw[1] = sinq( qtmp ) * qtmp1;

      qtmp1 = expq( -qr_rr ) / qtmp2;
      _QCO_TI_EQ_RE( qw, qtmp1 );
      /* fprintf(stdout, "# [qzeta_function] qw = %Qe   %Qe\n", qw[0], qw[1]); */

      _QKAHAN_SUM_CO_PL_CO(qterm2, qterm2c, qw);
      /* fprintf(stdout, "# [qzeta_function] term2 %3.0Qf %3.0Qf %3.0Qf \t %Qe \t %Qe\n", qn[0], qn[1], qn[2], qterm2[0], qterm2[1]); */



      /***************
       * FIRST TERM  *
       ***************/
      if(zeros_type < 3) {
        qtmp = qnd * qddInv * ( qgamma - 1.0Q );

        qr[0] = -M_PIq * ( qn[0] + qtmp * qdvec[0] ); 
        qr[1] = -M_PIq * ( qn[1] + qtmp * qdvec[1] ); 
        qr[2] = -M_PIq * ( qn[2] + qtmp * qdvec[2] ); 

        qr_rr    = _QSQR(qr);
        qr_rrInv = 1.0Q / qr_rr;
        qr_r     = sqrtq(qr_rr);

        if(qr_r < QZETA_EPSILON) {
          qr_costheta = 1.0Q;
          qr_phi      = 0.0Q;
        } else {
          qr_costheta = qr[2] / qr_r;

          qtmp = qr[0] * qr[0] + qr[1] * qr[1];
          if(qtmp < QZETA_EPSILON) {
            qr_phi = 0;
          } else {
            qr_phi = atan2q( qr[1], qr[0]);
          }
        }

        /* fprintf(stdout, "# [qzeta_function] qr=(%Qe, %Qe, %Qe) qr_r = %Qe qr_costheta = %Qe qr_phi = %Qe qr_rr = %Qe\n", qr[0], qr[1], qr[2], qr_r, qr_costheta, qr_phi, qr_rr); */

        int_parameters[1] = (double)qr_rr;
        F.params          = (void*)int_parameters;
  
        F.function = qintegrand_12;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_12_value, &int_integral_12_error);
        if(status != 0) {
          fprintf(stderr, "[qzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        /* fprintf(stdout, "# [qzeta_function] integration 12 result   = %25.16e%25.16e\n", int_integral_12_value, int_integral_12_error); */

        F.function = qintegrand_32;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_32_value, &int_integral_32_error);
        if(status != 0) {
          fprintf(stderr, "[qzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        /* fprintf(stdout, "# [qzeta_function] integration 32 result   = %25.16e%25.16e\n", int_integral_32_value, int_integral_32_error); */

        qint_iterate_const = expq(-(qr_rr - qq2));

        qtmp1 = (__float128)int_integral_32_value;
        qtmp2 = (__float128)int_integral_12_value;

        for(i=0; i<l+2; i++) {
          qtmp = ( qint_iterate_const - qq2 * qtmp1 + ((__float128)i - 1.5Q)  * qtmp2 ) * qr_rrInv;
          qtmp1 = qtmp2;
          qtmp2 = qtmp;
        }
        qint_integral_value = qtmp2;
        /* fprintf(stdout, "# [qzeta_function] integral value          = %25.16Qe\n", qint_integral_value); */

        /* TEST */
/*
        F.function = qintegrand_lp32;
        int_parameters[2] = (double)l;
        status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_32_value, &int_integral_32_error);
        if(status != 0) {
          fprintf(stderr, "[qzeta_function] Error from qag, status was %d\n", status);
          return(1);
        }
        fprintf(stdout, "# [qzeta_function] integration lp32 result = %25.16e%25.16e\n", int_integral_32_value, int_integral_32_error);
*/

        qtmp  = (__float128)m * qr_phi;
        qtmp1 = powq(qr_r, l) * qplm (l, m, qr_costheta);

        qw[0] = qtmp1 * qint_integral_value;
        qw[1] = sinq( qtmp ) * qw[0];
        qw[0] *= cosq( qtmp );
        /* fprintf(stdout, "# [qzeta_function] qw = %25.16Qe + I %25.16Qe\n", qw[0], qw[1]); */

        qtmp = M_PIq * qA * qnd;
        qint_norm_var[0] =  cosq(qtmp) * (__float128)      parity_factor;
        qint_norm_var[1] =  sinq(qtmp) * (__float128)( 1 - parity_factor );
        /* fprintf(stdout, "# [qzeta_function] qint_norm_var = %25.16Qe + I %25.16Qe\n", qint_norm_var[0], qint_norm_var[1]); */

            
        _QCO_EQ_CO_TI_CO(qw2, qw, qint_norm_var);
        /* fprintf(stdout, "# [qzeta_function] qw2 = %25.16Qe + I %25.16Qe\n", qw2[0], qw2[1]); */

        _QKAHAN_SUM_CO_PL_CO(qterm1, qterm1c, qw2);


        /* fprintf(stdout, "# [qzeta_function] term1 %3.0Qf %3.0Qf %3.0Qf \t %25.16Qe \t %25.16Qe\n", qr[0], qr[1], qr[2], qterm1[0], qterm1[1]); */
      }  /* end of if zeros_type < 3 */





    }  /* end of loop on rotations and reflections */



  }  /* end of loop on k1 */
  }  /* end of loop on k2 */

    /* fprintf(stdout, "# [qzeta_function] qconvergence %3d \t %25.16Qe \t %25.16Qe\n", k3, qterm2[0], qterm2[1]); */

  }  /* end of loop on k3 */

  /* subtraction for l = 0 and m = 0 */
  if(l == 0 && m == 0) {

    F.function = qint_l0m0_kernelFunction;
    F.params   = (void*)(&int_l0m0_parameters);
    status = gsl_integration_qag (&F, 0., 1., int_epsabs, int_epsrel, int_limit, int_key, int_workspace, &int_integral_l0m0_value, &int_integral_l0m0_error);
    if(status != 0) {
      fprintf(stderr, "[qzeta_function] Error from qag, status was %d\n", status);
      return(1);
    }
    qtmp = qint_l0m0_norm_const * (__float128)int_integral_l0m0_value + qint_l0m0_add;
    /* TEST */
    /* fprintf(stdout, "# [qzeta_function] (l=0 m=0) modified integral value %16.7e  %25.16Qe\n", q2, qtmp); */
    qterm3[0] = qtmp;
    qterm3[1] = 0.;
  }


  if(int_workspace != NULL) {
    gsl_integration_workspace_free(int_workspace);
  }

  /* multiply qterm2 with exp( q2 ) */
  qtmp = expq(  qq2 );
  _QCO_TI_EQ_RE(qterm2, qtmp);

  /* multiply qterm1 with constant integral normalization */
  _QCO_EQ_CO(qw, qterm1);
  _QCO_EQ_CO_TI_CO(qterm1, qw, qint_norm_const);

  z[0] = (double)( qterm1[0] + qterm2[0] + qterm3[0] );
  z[1] = (double)( qterm1[1] + qterm2[1] + qterm3[1] );

  return(0);
}  /* end of qzeta_function */
