/***************************************************************************
                          externs.h  -  description
                             -------------------
    begin                : Thu Apr 11 2002
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

// this file contains the "extern" version of the global variables. It's the mirror of globals.h.
// Include globals.h in the driver, and externs.h in the library files. If you include globals.h 
// twice, you get linking errors.
 
 
#ifndef EXTERNS_H
#define EXTERNS_H

#include "structs.h"
#include "constants.h"


extern int debug;
extern int ignore_internal;      // for debugging purposes
extern int ignore_multi;      // for debugging purposes
extern int fix_dangles;
extern int simple_internal_energy;
extern int simple_dangling_ends;
extern int no_dangling_ends;    // if 1, don't add dangling ends at all
extern int max_internal_loop;
extern int *constraints;

extern char similarity_rule[MAXNUMPARAMS][2000];

extern char string_params[MAXNUMPARAMS][MAXPNAME];    // for playing with the parameters
extern char string_params_human_readable[MAXNUMPARAMS][MAXPNAME]; 
extern int num_params;

extern char bbseq_left[MAXNUMPARAMS][10];   // the 5' sequence of the corresponding building block
extern char bbseq_right[MAXNUMPARAMS][10];   // the 3' sequence of the corresponding building block
extern char bbstr_left[MAXNUMPARAMS][10];   // the 5' structure of the corresponding building block
extern char bbstr_right[MAXNUMPARAMS][10];   // the 3' structure of the corresponding building block

// this is used for parameter learning only
extern int counter_min_dangle[NUM_DANG][NUM_DANG];                    // cell[0][1] says how many times min(dangle0,dangle1) appear, where dangle0 is x213 in the 318 nomenclature
// there are 48 dangling ends: 24 top and 24 bottom. 


extern PARAMTYPE stack [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE tstackh [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE tstacki [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE int11   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE int21   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE int22   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE dangle_top  [NUCL] [NUCL] [NUCL];
extern PARAMTYPE dangle_bot  [NUCL] [NUCL] [NUCL];
extern PARAMTYPE internal_penalty_by_size [MAXLOOP+1];
extern PARAMTYPE bulge_penalty_by_size [MAXLOOP+1];
extern PARAMTYPE hairpin_penalty_by_size [MAXLOOP+1];
extern miscinfo misc;

#if (MODEL == SIMPLE)
extern hairpin_tloop triloop[MAXTRILOOPNO];
extern hairpin_tloop tloop[MAXTLOOPNO];
extern int nb_triloops;
extern int nb_tloops;

#elif (MODEL == EXTENDED)
extern hairpin_tloop special_hl[MAX_SPECIAL_LOOP_NO];
extern int nb_special_hl;

// middle of asymmetric internal loops 2x2
extern PARAMTYPE int22mid[NUCL] [NUCL] [NUCL] [NUCL]; 
extern PARAMTYPE int22mid_group1;      // group 1 according to Christiansen_Znosko_2008
// That is: a U · U pair adjacent to an R · R pair, a G · A  or A · G pair adjacent to a Y · Y pair, or  any combination of A · C, U · C, C · U,  C · C, C · A, or A · A pairs
extern PARAMTYPE int22mid_group2;      // group 2 according to Christiansen_Znosko_2008
// That is: any combination of adjacent G · A and A · G pairs or two U · U pairs
extern PARAMTYPE int22mid_group3;      // group 3 according to Christiansen_Znosko_2008
// That is: a U · U pair adjacent to a Y · Y (not U · U),  C · A, or A · C pair
extern PARAMTYPE int22mid_group4;      // group 4 according to Christiansen_Znosko_2008
// That is: a G · G pair not adjacent to a U · U pair

extern PARAMTYPE int11_experimental_addition   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];       // values to be added to the simple 10-parameter model proposed by Davis_Znosko_2007, so that we use the experimental values for these parameters
extern PARAMTYPE int21_experimental_addition   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];       // values to be added to the simple 6-parameter model proposed by Badhwar_Znosko_2007, so that we use the experimental values for these parameters

//extern PARAMTYPE internal_asymmetry [MAXLOOP+1];
//extern PARAMTYPE internal_symmetry [MAXLOOP/2+1];
extern PARAMTYPE internal_penalty_by_size_2D[MAXLOOP_I][MAXLOOP_I];

extern PARAMTYPE bulgeA;
extern PARAMTYPE bulgeC;
extern PARAMTYPE bulgeG;
extern PARAMTYPE bulgeU;
extern PARAMTYPE bulge1[NUCL] [NUCL] [NUCL] [NUCL] [NUCL];     // bulge of size 1
#endif


// enthalpies information
extern PARAMTYPE enthalpy_stack [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_tstackh [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_tstacki [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_int11   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_int21   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_int22   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_dangle_top  [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_dangle_bot  [NUCL] [NUCL] [NUCL];
extern PARAMTYPE enthalpy_internal_penalty_by_size [MAXLOOP+1];
extern PARAMTYPE enthalpy_bulge_penalty_by_size [MAXLOOP+1];
extern PARAMTYPE enthalpy_hairpin_penalty_by_size [MAXLOOP+1];
extern miscinfo enthalpy_misc;
extern hairpin_tloop enthalpy_triloop[MAXTRILOOPNO];
extern hairpin_tloop enthalpy_tloop[MAXTLOOPNO];
extern int enthalpy_nb_triloops;
extern int enthalpy_nb_tloops;

// parameters from the configuration file
extern int nb_params;
extern char par_name [57] [100];
extern char par_value [57] [100];
extern char * std_dir_par;
extern char * rna_stack_energy37_v31_par;
extern char * rna_stack_energy37_v23_par;
extern char * rna_stack_enthalpy_par;
extern char * rna_loop_energy37_v31_par;
extern char * rna_loop_energy37_v23_par;
extern char * rna_loop_enthalpy_par;
extern char * rna_triloop_energy37_v31_par;
extern char * rna_triloop_energy37_v23_par;
extern char * rna_triloop_enthalpy_par;
extern char * rna_tloop_energy37_v31_par;
extern char * rna_tloop_energy37_v23_par;
extern char * rna_tloop_enthalpy_par;
extern char * rna_tstackh_energy37_v31_par;
extern char * rna_tstackh_energy37_v23_par;
extern char * rna_tstackh_enthalpy_par;
extern char * rna_tstacki_energy37_v31_par;
extern char * rna_tstacki_energy37_v23_par;
extern char * rna_tstacki_enthalpy_par;
extern char * rna_int11_energy37_v31_par;
extern char * rna_int11_energy37_v23_par;
extern char * rna_int11_enthalpy_par;
extern char * rna_int21_energy37_v31_par;
extern char * rna_int21_energy37_v23_par;
extern char * rna_int21_enthalpy_par;
extern char * rna_int22_energy37_v31_par;
extern char * rna_int22_energy37_v23_par;
extern char * rna_int22_enthalpy_par;
extern char * rna_miscloop_energy37_v31_par;
extern char * rna_miscloop_energy37_v23_par;
extern char * rna_miscloop_enthalpy_par;
extern char * rna_dangle_energy37_v31_par;
extern char * rna_dangle_energy37_v23_par;
extern char * rna_dangle_enthalpy_par;

extern char * dna_stack_energy37_par;
extern char * dna_stack_enthalpy_par;
extern char * dna_loop_energy37_par;
extern char * dna_loop_enthalpy_par;
extern char * dna_triloop_energy37_par;
extern char * dna_triloop_enthalpy_par;
extern char * dna_tloop_energy37_par;
extern char * dna_tloop_enthalpy_par;
extern char * dna_tstackh_energy37_par;
extern char * dna_tstackh_enthalpy_par;
extern char * dna_tstacki_energy37_par;
extern char * dna_tstacki_enthalpy_par;
extern char * dna_int11_energy37_par;
extern char * dna_int11_enthalpy_par;
extern char * dna_int21_energy37_par;
extern char * dna_int21_enthalpy_par;
extern char * dna_int22_energy37_par;
extern char * dna_int22_enthalpy_par;
extern char * dna_miscloop_energy37_par;
extern char * dna_miscloop_enthalpy_par;
extern char * dna_dangle_energy37_par;
extern char * dna_dangle_enthalpy_par;
extern char * rna_special_hl_energy_par;



extern double stack_double [NUCL] [NUCL] [NUCL] [NUCL];
extern double tstackh_double [NUCL] [NUCL] [NUCL] [NUCL];
extern double tstacki_double [NUCL] [NUCL] [NUCL] [NUCL];
extern double int11_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern double int21_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern double int22_double   [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL] [NUCL];
extern double dangle_top_double  [NUCL] [NUCL] [NUCL];
extern double dangle_bot_double  [NUCL] [NUCL] [NUCL];
extern double internal_penalty_by_size_double [MAXLOOP+1];
extern double bulge_penalty_by_size_double [MAXLOOP+1];
extern double hairpin_penalty_by_size_double [MAXLOOP+1];
extern miscinfo_double misc_double;
extern hairpin_tloop_double triloop_double[MAXTRILOOPNO];
extern hairpin_tloop_double tloop_double[MAXTLOOPNO];



#endif
