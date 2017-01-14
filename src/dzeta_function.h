#ifndef _DZETA_FUNCTION_H
#define _DZETA_FUNCTION_H


double dintegrand_12 (double t, void*params);

double dintegrand_32 (double t, void*params);

double dintegrand_lp32 (double t, void*params);

double dint_l0m0_kernelFunction (double t, void*params);

double dint_l0m0_kernelFunction2 (double t, void*params);

int dinit_plm_norm (void);

int dzeta_function (double z[2], double q2, int l, int m, int*dvec, double gamma, double A, double epsAbs, double epsRel, int Lmax);

#endif
