#include "init_zeta_function.h"
/**************************************************************************
 * rotations and reflections
 **************************************************************************/

#define _V3_EQ_V3(_a, _b) ( ( (_a)[0]==(_b)[0] ) && ( (_a)[1]==(_b)[1] ) && ( (_a)[2]==(_b)[2] ) )

int init_refrot_list (int rotref_list[48][2][3], int rotref_selection_permutation[4][4][48][3], int rotref_selection_sign[4][4][48][3], int rotref_number[4][4]) {
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
