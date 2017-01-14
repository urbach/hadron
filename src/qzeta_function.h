#ifndef _QZETA_FUNCTION_H
#define _QZETA_FUNCTION_H


double qintegrand_12 (double t, void*params);

double qintegrand_32 (double t, void*params);

double qintegrand_lp32 (double t, void*params);

double qint_l0m0_kernelFunction (double t, void*params);

double qint_l0m0_kernelFunction2 (double t, void*params);

int qinit_plm_norm (void);

int qzeta_function (double z[2], double q2, int l, int m, int*dvec, double gamma, double A, double epsAbs, double epsRel, int Lmax);

#endif
