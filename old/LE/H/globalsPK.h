/***************************************************************************
                          globals.h  -  description
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

/***************************************************************************
 * If modifications are made to this file, change externs.h accordingly. 
 * 
 ***************************************************************************/

#ifndef GLOBALSPK_H
#define GLOBALSPK_H

#include "constants.h"
#include "structsPK.h"

//#include "structs.h"  // July 16 - removed - unnecessary


// *** Add new energy model code here
pkmodelinfoDP pkmodelDP;  // energy constants used for the pseudoknotted energy model (Dirks & Pierce) 

pkmodelinfoRE pkmodelRE;  // energy constants used for the pseudoknotted energy model (Rivas & Eddy)

pkmodelinfoCC2006 pkmodelCC2006;  // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
float cc2006_s2_l1 [CC2006_STEMSIZE] [CC2006_LOOPSIZE];  // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
float cc2006_s1_l2 [CC2006_STEMSIZE] [CC2006_LOOPSIZE];  // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
float cc2006_s2_formula [4] [CC2006_STEMSIZE_FORMULA];    // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)
float cc2006_s1_formula [4] [CC2006_STEMSIZE_FORMULA];    // energy constants used for the pseudoknotted energy model (Cao & Chen 2006)

int coaxstack_f_a [NUCL] [NUCL] [NUCL] [NUCL];  // coaxial stacking energy parameters (as described by Cao & Chen)
int coaxstack_f_b [NUCL] [NUCL] [NUCL] [NUCL];  // flush coaxial stacking energy parameters (as described by Mathews)
int coaxstack_m1 [NUCL] [NUCL] [NUCL] [NUCL];  // mismatch coaxial stacking for continuous backbone (as described by Mathews)
int coaxstack_m2 [NUCL] [NUCL] [NUCL] [NUCL];  // mismatch coaxial stacking for discontinuous backbone (as described by Mathews)

char string_params_PK_CC[5000][MAXPNAME];    // for playing with the parameters

// parameters from the configuration file

// *** Add new energy model code here to the end of this matrix (change the first index to +2)
char par_namePK [14] [100] = {
            "RNA_PKMODEL_DP",
            "DNA_PKMODEL_DP",
            "RNA_PKMODEL_RE",
            "DNA_PKMODEL_RE",
            "RNA_PKMODEL_CC2006",
            "DNA_PKMODEL_CC2006",
            "RNA_COAXSTACK_F_A_ENERGY37",  // flush, model a
            "DNA_COAXSTACK_F_A_ENERGY37",  // flush, model a
            "RNA_COAXSTACK_F_B_ENERGY37",  // flush, model b
            "DNA_COAXSTACK_F_B_ENERGY37",  // flush, model b
            "RNA_COAXSTACK_M1_ENERGY37",  // mismatch, stack spanning continuous backbone
            "DNA_COAXSTACK_M1_ENERGY37",  // mismatch, stack spanning continuous backbone
            "RNA_COAXSTACK_M2_ENERGY37",  // mismatch, stack spanning broken backbone
            "DNA_COAXSTACK_M2_ENERGY37"  // mismatch, stack spanning broken backbone
            };

// *** Add new energy model code here to the end of this matrix (change the first index to +2)
char par_valuePK [14] [100];
char * rna_pkmodelDP_par   			= par_valuePK[0];
char * dna_pkmodelDP_par   			= par_valuePK[1];
char * rna_pkmodelRE_par   			= par_valuePK[2];
char * dna_pkmodelRE_par   			= par_valuePK[3];
char * rna_pkmodelCC2006_par   		= par_valuePK[4];
char * dna_pkmodelCC2006_par   		= par_valuePK[5];
char * rna_coaxstack_f_a_par   		= par_valuePK[6];
char * dna_coaxstack_f_a_par   		= par_valuePK[7];
char * rna_coaxstack_f_b_par   		= par_valuePK[8];
char * dna_coaxstack_f_b_par   		= par_valuePK[9];
char * rna_coaxstack_m1_par	   		= par_valuePK[10];
char * dna_coaxstack_m1_par	   		= par_valuePK[11];
char * rna_coaxstack_m2_par	   		= par_valuePK[12];
char * dna_coaxstack_m2_par	   		= par_valuePK[13];

// PARAMETER TUNING

int g_count_Ps = 0;		// pseudoloop initiation energy
int g_count_Psm = 0;		// penalty for introducing pseudoknot inside a multiloop
int g_count_Psp = 0;		// penalty for introducting pseudoknot inside a pseudoloop
int g_count_Pb = 0;		// penalty for band
int g_count_Pup = 0;		// penalty for unpaired base in pseudoloop or band
int g_count_Pps = 0;		// penalty for nested closed region inside either a pseudoloop or a multiloop that spans a band
int g_count_stP = 0;		// multiplicative penalty for stacked pair that spans a band 
int g_count_intP = 0;	// multiplicative penalty for internal loop that spans a band
int g_count_a = 0;		// penalty for introducing a multiloop
int g_count_a_p = 0;		// penalty for introducing a multiloop that spans a band
int g_count_b = 0;		// penalty for multiloop base pair 
int g_count_b_p = 0;		// penalty for multiloop base pair when the multiloop spans a band
int g_count_c = 0;		// penalty for unpaired base in multiloop
int g_count_c_p = 0;	



#endif
