/***************************************************************************
                          common.h  -  description
                             -------------------
    begin                : Thu Sep 5 2002
    copyright            : (C) 2002 by Mirela Andronescu
    email                : andrones@cs.ubc.ca
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef COMMONPK_H
#define COMMONPK_H

#include <stdio.h>
#include <math.h>
#include "structs.h"

/*
void giveup (char *string1, char *string2);
// to add: variable nb of parameters, as in scanf, printf

void empty_string (char * str);
    
int penalty_by_size (int size, char type);
// PRE:  size is the size of the loop
//       type is HAIRP or INTER or BULGE
// POST: return the penalty by size of the loop

int penalty_by_size_enthalpy (int size, char type);

void substr (char *source, int begin, int end, char *dest);
// PRE:  begin and end are smaller than strlen(source)
// POST: Put in dest what is in source between position begin and position end

*/

// from simfold/init.cpp
double ascii_to_doublePK (char *string);

int ascii_to_intPK(char *string);

double ascii_to_CCdouble_PK (char *string, double multiplier);

// from simfold/init.cpp
void giveupPK (char *string1, char *string2);

PARAMTYPE dangling_energy (int *sequence, char* structure, int i1, int i2, int i3, int i4);
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
// POST: return dangling energy between i2 and i3


PARAMTYPE dangling_energy_left (int *sequence, char* structure, int i1, int i2, int i3, int i4);
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i1 and i3


PARAMTYPE dangling_energy_right (int *sequence, char* structure, int i1, int i2, int i3, int i4);
//      (    (    )...)
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2

int dangling_energy_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34);
//      (   )...(   )
//      i1..i2..i3..i4
// PRE:  (i1, i2) and (i3, i4) are pairs, i2 and i3 are neighbours, i2 < i3
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy() if d_12 or d_13 are 1, they will be consider as restricted 
//			and not included in the dangling calculation
// POST: return dangling energy between i2 and i3

int dangling_energy_left_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34);
//      (....(    )   )
//      i1   i3  i4  i2
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy_left() if d_12 or d_13 are 1, they will be consider as restricted 
//		   and not included in the dangling calculation
// POST: return dangling energy between i1 and i3

int dangling_energy_right_res (int *sequence, int i1, int i2, int i3, int i4, int d_12, int d_34);
//      (    (    )...)
//      i1   i3  i4  i2
//		 d_12 = dangling end corresponding to base pair i1.i2
//		 d_34 = dangling end corresponding to base pair 13.i4
//       similar to dangling_energy_right() if d_12 or d_13 are 1, they will be consider as restricted 
//		   and not included in the dangling calculation
// PRE:  (i1, i2) and (i3, i4) are pairs, i1 and i3 are neighbours, i3 < i2
// POST: return dangling energy between i4 and i2

int LEcoax_stack_energy_flush_a (int i, int j, int ip, int jp, int *sequence);
// Pre: ip is i+1 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
// Post: Returns coaxial stacking energy from Walter et al 1999, or Turner parameters
//       if the given combination is not calculated in Walter et al. 1999, or 0 if
//       the returned value is INF (base pair is not canonical)

int LEcoax_stack_energy_flush_b (int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int other_j, int other_jp, int *sequence, int ignore_dangles);
// Pre: ip is i+1 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking energy from Turner lab / Mathews parameters (coaxial.dat).
//		 Currently, this is simply the Turner parameters for a regular stacked pair, or 0 if
//       the returned value is INF (base pair is not canonical)

int LEcoax_stack_energy_mismatch(int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int *sequence, char* structure, int other_j, int other_jp, int& stack_dangle, int ignore_dangles);
// Pre: ip is i+2 and jp = bp(ip) and j = bp(i); that is, i is always closest to 5' end compared to ip
//		flag = COAX_MULTI, COAX_PSEUDO, or COAX_OTHER
//		dangle_i corresponds to the base pair i.j
//		dangle_ip corresponds to the base pair ip.jp
// Post: Returns coaxial stacking mismatch energy from Turner lab / Mathews parameters.
//		 stack_dangle is set to 0 or 1 to indicate whether the i or ip dangle, respectively, is involved in the stacking


void count_LEcoax_stack_energy_flush_b (int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int other_j, int other_jp, int *sequence, double* counter, int ignore_dangles);

void count_LEcoax_stack_energy_mismatch(int i, int j, int ip, int jp, int flag, int dangle_i, int dangle_ip, int *sequence, char* structure, int other_j, int other_jp, int& stack_dangle, double* counter, int ignore_dangles);

PARAMTYPE LEstacked_pair_energy (int i, int j, int *sequence);

PARAMTYPE LEhairpin_loop_energy (int i, int j, int* sequence, char *csequence);
// PRE:  seq is the sequence of the loop; important for tetraloops
//       I assume i-j can pair
// POST: Help function to compute_hairpin_matrix
//       Return the energy of a hairpin

PARAMTYPE LEinternal_loop_energy (int i, int j, int ip, int jp, int *sequence);
// PRE:  The energy matrix was calculated
// POST: Read from the read pool, write to nodes;
//       Store the node and return the energy.

void detect_original_PKed_pairs_limited(char *structure, int *p_table);
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair

void detect_original_PKed_pairs_many(char *structure, short *p_table);
// PRE:  structure contains the desired structure
// POST: pairs will contain the index of each base pair
//               or -1 if it does not pair

void h_init (stack_ds *st);
void h_push (stack_ds *st, int el);
int h_pop (stack_ds *st);

void create_string_params_PK_CC();
int structure_type_index_PK(char type[]);
int structure_type_index_PK_CC(char type[]);

double compute_PK_sensitivity (char *ref_structure, char *pred_structure);
double compute_PK_ppv (char *ref_structure, char *pred_structure);

///////////////// FUNCTIONS FROM SIMFOLD ////////////////////
void count_LEdangling_energy (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter);
void count_LEdangling_energy_left (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter);
void count_LEdangling_energy_right (int *sequence, char *structure, int link, int i1, int i2, int i3, int i4, double *counter);

#endif
